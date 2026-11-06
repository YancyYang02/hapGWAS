"""Real genotypes, two kinship regimes, same haplotype windows.

The scanned chromosome is chrom 40 in both the phased VCF and the PLINK bed.
Background SNPs are drawn from the other chromosomes of the 8,251,083-SNP
chip. Each window is then tested three ways:

  clean  kinship = 300 background SNPs only
  full   kinship = those 300 plus up to 40 chip SNPs inside the window +/- 20 kb
  lowo   the same full kinship, with those local SNPs removed by the
         Woodbury downdate

lowo should match clean. full should move only where those local SNPs tag
the haplotype. A second pass repeats the strongest windows with 2000
background SNPs, where the same 40 tags are a much smaller fraction of the
kinship. One further window is downdated inside the real 8,251,083-SNP GCTA
GRM; that change is expected to be negligible because 40/8251083 is tiny.
"""

from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np
import pandas as pd

from hapgwas.blocks import fixed_windows
from hapgwas.dosage import build_dosage, haplotypes_from_phased_snps
from hapgwas.exact_lrt import block_exact_lrts
from hapgwas.genotype import load_geno
from hapgwas.grm import load_gcta_grm, vanraden_factors
from hapgwas.io import align_ids
from hapgwas.lmm import fit_null_lmm
from hapgwas.omnibus import lambda_gc_from_p
from hapgwas.window_kinship import lowo_joint_lrt, prepare_window_basis, region_columns

BED = "/fs2/home/clwen/YancongYang/Genetic_Parameter_Estimation/gemma_463_gwas_updated"
IDMAP = "/fs2/home/clwen/YancongYang/Genetic_Parameter_Estimation/id_update_final.txt"
VCF = "/fs5/home/clwen/YancongYang/hapgwas_try/demo_data/phased_chr4_1p5k.vcf.gz"
PHENO = "/fs5/home/clwen/YancongYang/hapgwas_try/data/pheno_key.tsv"
GRM = "/fs2/home/clwen/YancongYang/Genetic_Parameter_Estimation/gcta_grm"
OUT = Path("/fs5/home/clwen/YancongYang/hapgwas_try/runs/lowo_tagged")
TRAITS = ["BW6", "AFWR"]
N_BG = 300
N_BG_WIDE = 2000
TAG_CAP = 40
FLANK_BP = 20_000
M_GCTA = 8_251_083


def _impute(dosage: np.ndarray) -> np.ndarray:
    out = np.array(dosage, dtype=np.float64, copy=True)
    for j in range(out.shape[1]):
        col = out[:, j]
        good = np.isfinite(col)
        if not np.any(good):
            col[:] = 0.0
        elif not np.all(good):
            col[~good] = float(np.mean(col[good]))
        out[:, j] = col
    return out


def _decode(raw: bytes, n: int) -> np.ndarray:
    byte = np.frombuffer(raw, dtype=np.uint8)
    code = np.empty(byte.size * 4, dtype=np.uint8)
    code[0::4] = byte & 3
    code[1::4] = (byte >> 2) & 3
    code[2::4] = (byte >> 4) & 3
    code[3::4] = (byte >> 6) & 3
    code = code[:n]
    dose = np.full(n, np.nan, dtype=np.float64)
    dose[code == 0] = 2.0
    dose[code == 2] = 1.0
    dose[code == 3] = 0.0
    return dose


def _read_bed(path: str, n: int, indices: np.ndarray) -> np.ndarray:
    indices = np.asarray(indices, dtype=np.int64)
    rec = (n + 3) // 4
    order = np.argsort(indices, kind="mergesort")
    sorted_idx = indices[order]
    out = np.empty((n, indices.size), dtype=np.float64)
    with open(path, "rb") as handle:
        magic = handle.read(3)
        if magic != b"\x6c\x1b\x01":
            raise SystemExit(f"unexpected PLINK bed magic: {magic!r}")
        for col, snp_i in enumerate(sorted_idx):
            handle.seek(3 + int(snp_i) * rec)
            out[:, order[col]] = _decode(handle.read(rec), n)
    return out


def _maf_keep(dosage: np.ndarray, lo: float = 0.05, hi: float = 0.95) -> np.ndarray:
    af = np.nanmean(dosage, axis=0) / 2.0
    return np.isfinite(af) & (af > lo) & (af < hi)


def _hap(geno, block):
    h1 = geno.hap1[:, block.snp_start : block.snp_end]
    h2 = geno.hap2[:, block.snp_start : block.snp_end]
    pairs = haplotypes_from_phased_snps(h1, h2)
    hap = build_dosage(pairs, min_freq=0.01, min_count=5, collapse_mode="freq", drop_ref=True)
    snp = h1.astype(float) + h2.astype(float)
    return hap.dosage, snp


def _pair(y: np.ndarray, background: np.ndarray, tags: np.ndarray, hap_dosage: np.ndarray, snp_dosage: np.ndarray):
    dosage = np.column_stack([background, tags]) if tags.size else background
    _z, _keep, kinship = vanraden_factors(dosage, maf_min=0.01)
    null = fit_null_lmm(y, K=kinship, method="lmm")
    stat, df, p_full, *_rest = block_exact_lrts(null, hap_dosage, snp_dosage)
    basis = prepare_window_basis(y, dosage, maf_min=0.01)
    n_bg = background.shape[1]
    cols = region_columns(basis, n_bg, n_bg + tags.shape[1], 0, n_bg + tags.shape[1])
    stat_l, df_l, p_lowo = lowo_joint_lrt(basis, hap_dosage, cols)
    h2 = float(null.vg / (null.vg + null.ve)) if (null.vg + null.ve) > 0 else None
    return {
        "p_full": p_full,
        "df_full": int(df),
        "p_lowo": p_lowo,
        "df_lowo": int(df_l),
        "n_removed": int(cols.size),
        "h2_full": h2,
        "stat_full": stat,
        "stat_lowo": stat_l,
    }


def _identity(background: np.ndarray, tags: np.ndarray) -> float:
    dosage = np.column_stack([background, tags])
    z_all, keep, k_all = vanraden_factors(dosage, maf_min=0.01)
    n_bg = background.shape[1]
    in_bg = keep < n_bg
    z_bg = z_all[:, in_bg]
    k_clean = (z_bg @ z_bg.T) / z_bg.shape[1]
    region = np.flatnonzero(~in_bg)
    rank = int(region.size)
    m = z_all.shape[1]
    scale = m / (m - rank)
    update = z_all[:, region] / np.sqrt(m - rank)
    k_lowo = scale * k_all - update @ update.T
    k_lowo = 0.5 * (k_lowo + k_lowo.T)
    return float(np.max(np.abs(k_lowo - k_clean)))


def _log10(p: float) -> float:
    if p is None or not np.isfinite(p) or p <= 0:
        return float("nan")
    return float(-np.log10(p))


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    print("reading bim", flush=True)
    n_snps = 8_251_083
    chrom = np.empty(n_snps, dtype=np.int16)
    pos = np.empty(n_snps, dtype=np.int64)
    a1 = np.empty(n_snps, dtype="U16")
    with open(BED + ".bim", encoding="utf-8") as handle:
        for i, line in enumerate(handle):
            parts = line.split()
            chrom[i] = int(parts[0])
            pos[i] = int(parts[3])
            a1[i] = parts[4]
            if i + 1 == n_snps:
                break
        else:
            raise SystemExit(f"bim ended at {i + 1}, expected {n_snps}")
    on40 = np.flatnonzero(chrom == 40)
    print(f"chrom40 snps={on40.size} other={n_snps - on40.size}", flush=True)

    fam = pd.read_csv(BED + ".fam", sep=r"\s+", header=None, names=["FID", "IID", "P", "M", "SEX", "PHE"], dtype=str)
    fam_index = {f"{fid} {iid}": i for i, (fid, iid) in enumerate(zip(fam["FID"], fam["IID"]))}
    a_to_row = {}
    with open(IDMAP, encoding="utf-8") as handle:
        for line in handle:
            parts = line.split()
            if len(parts) < 4:
                continue
            key = f"{parts[2]} {parts[3]}"
            if key not in fam_index:
                continue
            a_to_row[parts[0]] = fam_index[key]
    print(f"fam={len(fam)} mapped={len(a_to_row)}", flush=True)

    geno0 = load_geno(VCF)
    geno0.ids = [str(i) for i in geno0.ids]
    keep_ids = [iid for iid in geno0.ids if iid in a_to_row]
    id_pos = {iid: i for i, iid in enumerate(geno0.ids)}
    geno = geno0.subset_samples(np.asarray([id_pos[iid] for iid in keep_ids], dtype=int))
    bed_rows = np.asarray([a_to_row[iid] for iid in keep_ids], dtype=int)
    print(f"analysis samples={geno.n_samples}", flush=True)

    # Coding check at the first shared site.
    vcf_pos = int(geno.positions[0])
    hit = on40[pos[on40] == vcf_pos]
    if hit.size != 1:
        raise SystemExit(f"expected one bim hit at 40:{vcf_pos}, found {hit.size}")
    bed_dose = _read_bed(BED + ".bed", len(fam), hit)[:, 0][bed_rows]
    alt = geno.snp_dosage(0)
    ok = np.isfinite(bed_dose)
    corr = float(np.corrcoef(bed_dose[ok], alt[ok])[0, 1])
    ref = str(geno.alleles[0, 0]) if geno.alleles is not None else "?"
    alt_allele = str(geno.alleles[0, 1]) if geno.alleles is not None else "?"
    coding = {
        "vcf_pos": vcf_pos,
        "bim_a1": str(a1[int(hit[0])]),
        "vcf_ref": ref,
        "vcf_alt": alt_allele,
        "n_compared": int(ok.sum()),
        "corr_a1_vs_alt_dosage": corr,
    }
    print("coding", coding, flush=True)
    if abs(corr) < 0.98:
        raise SystemExit("PLINK dosage does not match the phased VCF; stopping")
    flip = corr < 0
    coding["flipped_to_match_alt"] = flip

    def load_aligned(indices: np.ndarray) -> np.ndarray:
        raw = _read_bed(BED + ".bed", len(fam), indices)
        if flip:
            raw = 2.0 - raw
        return _impute(raw[bed_rows])

    other = np.flatnonzero(chrom != 40)
    step = max(1, other.size // 8000)
    candidates = other[::step][:8000]
    print(f"reading {candidates.size} background candidates", flush=True)
    bg_all = load_aligned(candidates)
    passed = _maf_keep(bg_all)
    bg_idx = np.flatnonzero(passed)[:N_BG_WIDE]
    if bg_idx.size < N_BG:
        raise SystemExit(f"only {bg_idx.size} background SNPs passed MAF")
    bg_wide = bg_all[:, bg_idx]
    bg = bg_wide[:, :N_BG]
    print(f"background kept wide={bg_wide.shape[1]} narrow={bg.shape[1]}", flush=True)

    blocks = fixed_windows(geno.n_snps, window=4, step=2, chrom="40", positions=geno.positions)
    pos40 = pos[on40]
    tag_index = []
    for block in blocks:
        mid = 0.5 * (int(block.start_bp) + int(block.end_bp))
        left = int(np.searchsorted(pos40, mid - FLANK_BP, side="left"))
        right = int(np.searchsorted(pos40, mid + FLANK_BP, side="right"))
        if right - left < 5:
            left = int(np.searchsorted(pos40, mid - 100_000, side="left"))
            right = int(np.searchsorted(pos40, mid + 100_000, side="right"))
        choose = on40[left:right]
        if choose.size > TAG_CAP:
            dist = np.abs(pos[choose].astype(np.float64) - mid)
            choose = choose[np.argsort(dist)[:TAG_CAP]]
            choose.sort()
        tag_index.append(choose)
    unique = np.unique(np.concatenate(tag_index))
    print(f"reading {unique.size} tag SNPs for {len(blocks)} windows", flush=True)
    tag_dose = load_aligned(unique)
    tag_col = {int(snp_i): j for j, snp_i in enumerate(unique)}

    identity = _identity(bg, tag_dose[:, [tag_col[int(i)] for i in tag_index[0]]])
    print("identity_max_abs", identity, flush=True)

    pheno = pd.read_csv(PHENO, sep=None, engine="python")
    pheno["IID"] = pheno["IID"].astype(str)
    report = {
        "coding": coding,
        "identity_max_abs_K": identity,
        "n_background": N_BG,
        "n_background_wide": int(bg_wide.shape[1]),
        "tag_cap": TAG_CAP,
        "flank_bp": FLANK_BP,
        "n_samples_geno": int(geno.n_samples),
        "traits": {},
    }

    for trait in TRAITS:
        sub = pheno[["IID", trait]].dropna().rename(columns={trait: "PHENO"})
        y, keep, _x, ids = align_ids(sub, geno.ids, None)
        y_ids = ids
        # align_ids keep is into geno order, which is keep_ids order
        g = geno.subset_samples(keep)
        bg_t = bg[keep]
        bg_w = bg_wide[keep]
        tags_t = tag_dose[keep]
        _z, _k, k_clean = vanraden_factors(bg_t, maf_min=0.01)
        null_clean = fit_null_lmm(y, K=k_clean, method="lmm")
        h2_clean = float(null_clean.vg / (null_clean.vg + null_clean.ve))
        rows = []
        print(f"{trait} n={null_clean.n} h2_clean={h2_clean:.4f}", flush=True)
        for i, block in enumerate(blocks):
            hap_dosage, snp_dosage = _hap(g, block)
            these = tags_t[:, [tag_col[int(j)] for j in tag_index[i]]]
            compared = _pair(y, bg_t, these, hap_dosage, snp_dosage)
            _st, _df, p_clean, *_ = block_exact_lrts(null_clean, hap_dosage, snp_dosage)
            rows.append(
                {
                    "block_id": block.block_id,
                    "start_bp": block.start_bp,
                    "end_bp": block.end_bp,
                    "n_tags": int(tag_index[i].size),
                    "p_clean": p_clean,
                    "p_full": compared["p_full"],
                    "p_lowo": compared["p_lowo"],
                    "n_removed": compared["n_removed"],
                    "h2_full": compared["h2_full"],
                }
            )
            if i % 100 == 0:
                print(f"  {trait} {i}/{len(blocks)}", flush=True)
        frame = pd.DataFrame(rows)
        frame.to_csv(OUT / f"{trait}_windows.tsv", sep="\t", index=False, float_format="%.6g")
        d_lowo = np.abs(frame["p_lowo"].map(_log10) - frame["p_clean"].map(_log10))
        d_full = frame["p_full"].map(_log10) - frame["p_clean"].map(_log10)
        order = frame.sort_values("p_clean").head(15)
        _z, _k, k_wide = vanraden_factors(bg_w, maf_min=0.01)
        null_wide = fit_null_lmm(y, K=k_wide, method="lmm")
        h2_wide = float(null_wide.vg / (null_wide.vg + null_wide.ve))
        wide = []
        top_ids = set(order["block_id"])
        for i, block in enumerate(blocks):
            if block.block_id not in top_ids:
                continue
            hap_dosage, snp_dosage = _hap(g, block)
            these = tags_t[:, [tag_col[int(j)] for j in tag_index[i]]]
            compared = _pair(y, bg_w, these, hap_dosage, snp_dosage)
            _st, _df, p_clean_w, *_ = block_exact_lrts(null_wide, hap_dosage, snp_dosage)
            wide.append(
                {
                    "block_id": block.block_id,
                    "start_bp": block.start_bp,
                    "end_bp": block.end_bp,
                    "p_clean": p_clean_w,
                    "p_full": compared["p_full"],
                    "p_lowo": compared["p_lowo"],
                    "n_removed": compared["n_removed"],
                    "h2_full": compared["h2_full"],
                }
            )
        wide_df = pd.DataFrame(wide)
        wide_df.to_csv(OUT / f"{trait}_wide_top15.tsv", sep="\t", index=False, float_format="%.6g")
        report["traits"][trait] = {
            "n": int(null_clean.n),
            "n_ids": len(y_ids),
            "h2_clean_300": h2_clean,
            "h2_clean_wide": h2_wide,
            "n_blocks": int(len(frame)),
            "lambda_clean": float(lambda_gc_from_p(frame["p_clean"])),
            "lambda_full": float(lambda_gc_from_p(frame["p_full"])),
            "lambda_lowo": float(lambda_gc_from_p(frame["p_lowo"])),
            "median_abs_log10_lowo_vs_clean": float(np.nanmedian(d_lowo)),
            "p95_abs_log10_lowo_vs_clean": float(np.nanpercentile(d_lowo, 95)),
            "median_log10_full_minus_clean": float(np.nanmedian(d_full)),
            "top10_clean": frame.sort_values("p_clean").head(10).to_dict(orient="records"),
            "wide_background": N_BG_WIDE,
            "wide_top": wide,
        }
        print(
            trait,
            "median |log10 lowo-clean|",
            round(report["traits"][trait]["median_abs_log10_lowo_vs_clean"], 4),
            "median log10 full-clean",
            round(report["traits"][trait]["median_log10_full_minus_clean"], 4),
            flush=True,
        )

    # Dense genome-wide GRM: remove 40 tags from the real GCTA matrix at the
    # AFWR window with the smallest clean p. This is the negative control.
    trait = "AFWR"
    top = min(report["traits"][trait]["top10_clean"], key=lambda r: r["p_clean"])
    block = next(b for b in blocks if b.block_id == top["block_id"])
    bi = blocks.index(block)
    sub = pheno[["IID", trait]].dropna().rename(columns={trait: "PHENO"})
    y, keep, _x, ids = align_ids(sub, geno.ids, None)
    g = geno.subset_samples(keep)
    grm_ids, kinship = load_gcta_grm(GRM)
    index = {iid: i for i, iid in enumerate(grm_ids)}
    ix = np.asarray([index[iid] for iid in ids], dtype=int)
    k_use = kinship[np.ix_(ix, ix)]
    hap_dosage, snp_dosage = _hap(g, block)
    null = fit_null_lmm(y, K=k_use, method="lmm")
    _st, _df, p_gcta, *_ = block_exact_lrts(null, hap_dosage, snp_dosage)
    tags = tag_dose[keep][:, [tag_col[int(j)] for j in tag_index[bi]]]
    z_tag, tag_keep, _k_tag = vanraden_factors(tags, maf_min=0.01)
    rank = z_tag.shape[1]
    scale = M_GCTA / (M_GCTA - rank)
    update = z_tag / np.sqrt(M_GCTA - rank)
    k_down = scale * k_use - update @ update.T
    k_down = 0.5 * (k_down + k_down.T)
    null_d = fit_null_lmm(y, K=k_down, method="lmm")
    _st, _df, p_down, *_ = block_exact_lrts(null_d, hap_dosage, snp_dosage)
    report["dense_grm"] = {
        "trait": trait,
        "block_id": block.block_id,
        "start_bp": block.start_bp,
        "end_bp": block.end_bp,
        "n_tags_removed": rank,
        "m_gcta": M_GCTA,
        "max_abs_K_change": float(np.max(np.abs(k_down - k_use))),
        "p_gcta": p_gcta,
        "p_downdated": p_down,
        "h2_gcta": float(null.vg / (null.vg + null.ve)),
        "h2_downdated": float(null_d.vg / (null_d.vg + null_d.ve)),
    }
    print("dense", report["dense_grm"], flush=True)
    report["seconds"] = round(time.time() - t0, 1)
    (OUT / "lowo_tagged.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    print("WROTE", OUT / "lowo_tagged.json", flush=True)


if __name__ == "__main__":
    main()
