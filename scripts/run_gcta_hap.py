"""Haplotype scan of the broiler panels under the genome-wide GCTA GRM.

The kinship is the 8,251,083-SNP GCTA matrix (465 birds), not the VanRaden
matrix of the scanned window. Region exclusion is not applied: that file has
no SNP factors to downdate. Discovery columns are the joint score, the exact
joint likelihood ratio, and omni_p. Lambda is computed on raw p-values.
"""

from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np
import pandas as pd

from hapgwas.blocks import fixed_windows
from hapgwas.genotype import load_geno
from hapgwas.grm import load_gcta_grm
from hapgwas.io import align_ids, write_assoc
from hapgwas.joint_test import run_genome
from hapgwas.lmm import fit_null_lmm
from hapgwas.omnibus import clump_overlapping, lambda_gc_from_p

PANELS = [
    ("chr4_1500", Path("/fs5/home/clwen/YancongYang/hapgwas_try/demo_data/phased_chr4_1p5k.vcf.gz")),
    ("first5000", Path("/fs5/home/clwen/YancongYang/hapgwas_try/data/phased_first5000.vcf.gz")),
]
PHENO = Path("/fs5/home/clwen/YancongYang/hapgwas_try/data/pheno_key.tsv")
GRM = "/fs2/home/clwen/YancongYang/Genetic_Parameter_Estimation/gcta_grm"
OUT = Path("/fs5/home/clwen/YancongYang/hapgwas_try/runs/gcta_hap")
TRAITS = ["BW6", "AFWR", "BW3_6", "EBW6", "GBW6", "ERFI"]
WINDOW, STEP = 4, 2
PCOLS = ("omni_p", "joint_p", "joint_lrt_p", "beyond_p", "beyond_lrt_p", "snp_bonf_p", "allele_bonf_p")


def _lambda(p: pd.Series) -> float:
    arr = pd.to_numeric(p, errors="coerce").to_numpy(dtype=float)
    return float(lambda_gc_from_p(arr))


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    grm_ids, kinship = load_gcta_grm(GRM)
    diag = float(np.mean(np.diag(kinship)))
    print(f"GRM n={len(grm_ids)} mean_diag={diag:.4f}", flush=True)
    if not 0.5 < diag < 2.0:
        raise SystemExit(f"GCTA GRM diagonal mean {diag} is not near 1; refusing to scan")
    pheno = pd.read_csv(PHENO, sep=None, engine="python")
    pheno["IID"] = pheno["IID"].astype(str)
    rows = []
    for panel, path in PANELS:
        t_panel = time.time()
        geno0 = load_geno(path)
        geno0.ids = [str(i) for i in geno0.ids]
        chrom = str(geno0.chrom[0]) if geno0.chrom is not None else "40"
        blocks = fixed_windows(geno0.n_snps, window=WINDOW, step=STEP, chrom=chrom, positions=geno0.positions)
        bp = {b.block_id: (b.start_bp, b.end_bp) for b in blocks}
        print(f"PANEL {panel} n={geno0.n_samples} snps={geno0.n_snps} blocks={len(blocks)}", flush=True)
        for trait in TRAITS:
            if trait not in pheno.columns:
                print("skip", trait, flush=True)
                continue
            t0 = time.time()
            sub = pheno[["IID", trait]].dropna().rename(columns={trait: "PHENO"})
            y, keep, _X, ids = align_ids(sub, geno0.ids, None)
            index = {iid: i for i, iid in enumerate(grm_ids)}
            missing = [iid for iid in ids if iid not in index]
            if missing:
                raise SystemExit(f"{trait}: {len(missing)} IDs missing from the GRM, e.g. {missing[:3]}")
            ix = np.asarray([index[iid] for iid in ids], dtype=int)
            k_use = kinship[np.ix_(ix, ix)]
            geno = geno0.subset_samples(keep)
            null = fit_null_lmm(y, K=k_use, method="lmm")
            block_df, allele_df = run_genome(
                null,
                blocks,
                geno.hap1,
                geno.hap2,
                min_freq=0.01,
                min_count=5,
                apply_fdr=True,
            )
            block_df["start_bp"] = block_df["block_id"].map(lambda i: bp.get(i, (None, None))[0])
            block_df["end_bp"] = block_df["block_id"].map(lambda i: bp.get(i, (None, None))[1])
            out = OUT / panel / trait
            out.mkdir(parents=True, exist_ok=True)
            write_assoc(block_df, out / "assoc_blocks.tsv")
            write_assoc(allele_df, out / "assoc_alleles.tsv")
            loci = clump_overlapping(block_df, "joint_lrt_p") if len(block_df) else block_df
            write_assoc(loci, out / "assoc_loci_lrt.tsv")
            n_blocks = max(len(block_df), 1)
            thr = 0.05 / n_blocks
            h2 = float(null.vg / (null.vg + null.ve)) if (null.vg + null.ve) > 0 else None
            summary = {
                "panel": panel,
                "trait": trait,
                "kinship": "gcta_grm_8251083",
                "n": int(null.n),
                "n_snps": int(geno.n_snps),
                "n_blocks": int(len(block_df)),
                "n_loci_lrt": int(len(loci)),
                "h2": h2,
                "vg": float(null.vg),
                "ve": float(null.ve),
                "threshold": thr,
                "seconds": round(time.time() - t0, 1),
            }
            for col in PCOLS:
                if col not in block_df.columns:
                    continue
                ps = pd.to_numeric(block_df[col], errors="coerce")
                summary[col] = {
                    "hits_blocks": int((ps < thr).sum()),
                    "hits_loci": int((pd.to_numeric(loci[col], errors="coerce") < thr).sum()) if col in loci.columns else None,
                    "min_p": float(ps.min()) if ps.notna().any() else None,
                    "lambda_gc": _lambda(ps),
                }
            if len(loci):
                top = loci.sort_values("joint_lrt_p").head(8)
                keep_cols = [
                    c
                    for c in (
                        "block_id",
                        "chrom",
                        "start_bp",
                        "end_bp",
                        "snp_start",
                        "snp_end",
                        "joint_lrt_p",
                        "joint_p",
                        "omni_p",
                        "beyond_lrt_p",
                        "snp_bonf_p",
                        "allele_bonf_p",
                    )
                    if c in top.columns
                ]
                summary["top_loci"] = top[keep_cols].to_dict(orient="records")
            (out / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
            rows.append(summary)
            brief = {
                c: summary[c]["hits_loci"]
                for c in ("omni_p", "joint_p", "joint_lrt_p", "snp_bonf_p")
                if c in summary
            }
            lams = {
                c: round(summary[c]["lambda_gc"], 3)
                for c in ("joint_p", "joint_lrt_p", "omni_p")
                if c in summary
            }
            print(trait, panel, "n", summary["n"], "h2", None if h2 is None else round(h2, 4), "loci", brief, "lam", lams, flush=True)
        print(f"PANEL {panel} done in {time.time() - t_panel:.1f}s", flush=True)
    (OUT / "gcta_hap.json").write_text(json.dumps(rows, indent=2), encoding="utf-8")
    print("WROTE", OUT / "gcta_hap.json", flush=True)


if __name__ == "__main__":
    main()
