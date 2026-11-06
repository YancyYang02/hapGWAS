"""Command-line interface for hapGWAS."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from hapgwas import __version__
from hapgwas.blocks import (
    blocks_from_frame,
    fixed_windows,
    gabriel_blocks,
    read_plink_blocks,
    summarize_blocks,
)
from hapgwas.demo import run_demo
from hapgwas.genotype import PhasedGeno, load_geno, save_npz
from hapgwas.grm import load_gcta_grm, vanraden_grm
from hapgwas.window_kinship import prepare_window_basis
from hapgwas.io import align_ids, read_covariates, read_phenotype, write_assoc, write_phased_vcf
from hapgwas.joint_test import _add_bh, apply_block_fdr, run_genome, score_snp_matrix
from hapgwas.lmm import NullFit, fit_null_lmm
from hapgwas.omnibus import clump_overlapping, lambda_gc_from_p


def _resolve_blocks(args, geno: PhasedGeno):
    n_snps = geno.n_snps
    positions = geno.positions
    if getattr(args, "blocks", None):
        # autodect PLINK vs index table
        peek = pd.read_csv(args.blocks, sep=None, engine="python", nrows=2)
        cols_u = {c.upper() for c in peek.columns}
        if "SNPS" in cols_u or ({"CHR", "BP1", "BP2"} <= cols_u):
            return read_plink_blocks(args.blocks, snp_ids=geno.snp_ids, positions=positions, chrom=geno.chrom)
        bdf = pd.read_csv(args.blocks, sep=None, engine="python")
        return blocks_from_frame(bdf, n_snps)

    mode = getattr(args, "block_mode", "window")
    chrom0 = str(geno.chrom[0]) if geno.chrom is not None else args.chrom
    if mode == "gabriel":
        return gabriel_blocks(
            geno.hap1,
            geno.hap2,
            positions=positions,
            chrom=chrom0,
            max_block_snps=args.gabriel_max_snps,
            min_snps=2,
        )
    if mode == "window":
        return fixed_windows(
            n_snps,
            window=args.window,
            step=args.step or args.window,
            chrom=args.chrom if geno.chrom is None else chrom0,
            positions=positions,
        )
    raise ValueError(f"Unknown block mode {mode}")


def _cmd_demo(args: argparse.Namespace) -> int:
    out = Path(args.out)
    summary = run_demo(out, use_lmm=not args.ols)
    print("hapGWAS demo finished")
    print(f"  out_dir:     {summary['out_dir']}")
    print(f"  null model:  {summary['null_method']}  vg={summary['vg']:.4g}  ve={summary['ve']:.4g}")
    print(f"  blocks:      {summary['n_blocks']}")
    print(f"  causal p:    {summary['causal_block_joint_p']}")
    print("  top blocks:")
    for row in summary["top_blocks"]:
        print(f"    {row['block_id']}  snps[{row['snp_start']}:{row['snp_end']}]  joint_p={row['joint_p']:.3g}")
    print("  wrote: assoc_blocks.tsv, assoc_alleles.tsv, blocks.tsv, pheno.tsv, phased_demo.npz")
    return 0


def _cmd_convert(args: argparse.Namespace) -> int:
    geno = load_geno(args.geno, companion=args.map or args.sample)
    out = Path(args.out)
    if out.suffix.lower() == ".vcf":
        write_phased_vcf(geno, out)
    else:
        if out.suffix.lower() != ".npz":
            out = out.with_suffix(".npz")
        save_npz(geno, out)
    print(f"Converted {geno.source} -> {out}  (n={geno.n_samples}, p={geno.n_snps})")
    if geno.meta.get("phase_assumed"):
        print("WARNING: PED allele order treated as phase (phase_assumed=True).")
    return 0


def _snp_table(geno: PhasedGeno, null: NullFit, columns: np.ndarray) -> pd.DataFrame:
    columns = np.asarray(columns, dtype=int)
    if columns.size == 0:
        return pd.DataFrame(columns=["snp_index", "chrom", "pos", "snp_id", "snp_chisq", "snp_p"])
    dosage = geno.snp_dosage_matrix()[:, columns]
    chisq, p = score_snp_matrix(null, dosage)
    chrom = np.asarray(["NA"] * geno.n_snps, dtype=object) if geno.chrom is None else geno.chrom.astype(str)
    pos = np.full(geno.n_snps, -1) if geno.positions is None else geno.positions.astype(int)
    sid = np.asarray([f"snp{i}" for i in range(geno.n_snps)], dtype=object) if geno.snp_ids is None else geno.snp_ids.astype(str)
    return pd.DataFrame(
        {
            "snp_index": columns,
            "chrom": chrom[columns],
            "pos": pos[columns],
            "snp_id": sid[columns],
            "snp_chisq": chisq,
            "snp_p": p,
        }
    )


def _kinship(args, geno: PhasedGeno, kept_ids: list[str]) -> tuple[np.ndarray | None, str]:
    """Return (K, description). OLS is (None, 'ols')."""
    if args.ols or args.grm == "none":
        if args.ols and getattr(args, "gcta_grm", None):
            raise ValueError("--ols ignores kinship; do not combine it with --gcta-grm")
        return None, "ols"
    if getattr(args, "gcta_grm", None):
        if args.kin_geno:
            raise ValueError("use either --gcta-grm or --kin-geno, not both")
        ids, kinship = load_gcta_grm(args.gcta_grm)
        index = {iid: i for i, iid in enumerate(ids)}
        missing = [iid for iid in kept_ids if iid not in index]
        if missing:
            raise ValueError(f"--gcta-grm is missing {len(missing)} phenotype IDs, e.g. {missing[:3]}")
        ix = np.asarray([index[iid] for iid in kept_ids], dtype=int)
        return kinship[np.ix_(ix, ix)], "gcta-grm"
    if args.kin_geno:
        kin = load_geno(args.kin_geno, companion=args.kin_map or args.kin_sample)
        index = {iid: i for i, iid in enumerate(kin.ids)}
        missing = [iid for iid in kept_ids if iid not in index]
        if missing:
            raise ValueError(f"--kin-geno is missing {len(missing)} phenotype IDs, e.g. {missing[:3]}")
        kin = kin.subset_samples(np.asarray([index[iid] for iid in kept_ids], dtype=int))
        return vanraden_grm(kin.snp_dosage_matrix()), "vanraden-kin-geno"
    if args.grm == "provided":
        if geno.K is None:
            raise ValueError("--grm provided requires a K array in the genotype file")
        return geno.K, "provided"
    if args.grm == "auto" and geno.K is not None:
        return geno.K, "provided"
    return None, "vanraden"


def _chrom_labels(geno: PhasedGeno) -> np.ndarray:
    if geno.chrom is None:
        return np.asarray(["NA"] * geno.n_snps, dtype=object)
    return geno.chrom.astype(str)


def _cmd_run(args: argparse.Namespace) -> int:
    geno = load_geno(args.geno, companion=args.map or args.sample)
    pheno = read_phenotype(args.pheno)
    cov = read_covariates(args.cov) if args.cov else None
    y, keep, X, kept_ids = align_ids(pheno, geno.ids, cov)
    geno = geno.subset_samples(keep)

    blocks = _resolve_blocks(args, geno)
    if args.write_blocks:
        summarize_blocks(blocks).to_csv(args.write_blocks, sep="\t", index=False)
    cond_snp = geno.find_snp(args.condition_snp) if args.condition_snp is not None else None

    K_fixed, kin_name = _kinship(args, geno, kept_ids)
    labels = _chrom_labels(geno)
    uniq = list(dict.fromkeys(labels.tolist()))
    loco = (
        kin_name == "vanraden"
        and K_fixed is None
        and args.grm in ("auto", "vanraden", "loco")
        and len(uniq) > 1
        and "NA" not in uniq
    )
    scan_kwargs = dict(
        min_freq=args.min_freq,
        min_count=args.min_count,
        collapse_mode=args.collapse,
        condition_snp=cond_snp,
        condition_auto_topsnp=args.condition_auto_topsnp,
        snp_ids=geno.snp_ids,
        apply_fdr=False,
    )
    lowo_basis = None
    if getattr(args, "lowo", False) and not args.ols:
        if kin_name == "gcta-grm":
            raise ValueError(
                "--lowo removes SNP columns from a VanRaden factorisation. "
                "A GCTA .grm.bin has no SNP factors, so it cannot be downdated. "
                "Leave --lowo off with --gcta-grm, or build the kinship with --kin-geno."
            )
        lowo_basis = prepare_window_basis(y, geno.snp_dosage_matrix(), X=X, maf_min=args.min_freq)
        scan_kwargs["lowo_basis"] = lowo_basis
        scan_kwargs["lowo_flank"] = args.lowo_flank
        scan_kwargs["n_snps"] = int(geno.n_snps)

    block_parts = []
    allele_parts = []
    snp_parts = []
    nulls = []
    if loco:
        for chrom in uniq:
            in_chrom = np.where(labels == chrom)[0]
            background = np.where(labels != chrom)[0]
            if background.size < 30:
                K = vanraden_grm(geno.snp_dosage_matrix())
                how = f"vanraden-all(fallback:{chrom})"
            else:
                K = vanraden_grm(geno.snp_dosage_matrix()[:, background])
                how = f"vanraden-loco:{chrom}"
            null = fit_null_lmm(y, X=X, K=K, method="lmm")
            nulls.append((how, null))
            these = [b for b in blocks if str(b.chrom) == chrom]
            if these:
                bdf, adf = run_genome(null, these, geno.hap1, geno.hap2, **scan_kwargs)
                block_parts.append(bdf)
                allele_parts.append(adf)
            snp_parts.append(_snp_table(geno, null, in_chrom))
    else:
        if kin_name == "vanraden":
            K_fixed = vanraden_grm(geno.snp_dosage_matrix())
            kin_name = "vanraden-same-file"
        method = "ols" if K_fixed is None else "lmm"
        null = fit_null_lmm(y, X=X, K=K_fixed, method=method)
        nulls.append((kin_name if K_fixed is not None else "ols", null))
        bdf, adf = run_genome(null, blocks, geno.hap1, geno.hap2, **scan_kwargs)
        block_parts.append(bdf)
        allele_parts.append(adf)
        snp_parts.append(_snp_table(geno, null, np.arange(geno.n_snps)))

    block_df = pd.concat(block_parts, ignore_index=True) if block_parts else pd.DataFrame()
    allele_df = pd.concat(allele_parts, ignore_index=True) if allele_parts else pd.DataFrame()
    snp_df = pd.concat(snp_parts, ignore_index=True) if snp_parts else pd.DataFrame()
    block_df, allele_df = apply_block_fdr(block_df, allele_df)
    if len(snp_df) and "snp_p" in snp_df.columns:
        snp_df = _add_bh(snp_df, "snp_p", "snp_q")
    loci = clump_overlapping(block_df, "omni_p") if len(block_df) else block_df

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    write_assoc(block_df, out / "assoc_blocks.tsv")
    write_assoc(allele_df, out / "assoc_alleles.tsv")
    write_assoc(snp_df, out / "assoc_snps.tsv")
    write_assoc(loci, out / "assoc_loci.tsv")

    n_blocks = max(len(block_df), 1)
    thr = 0.05 / n_blocks
    lambdas = {
        "joint": lambda_gc_from_p(block_df["joint_p"]) if "joint_p" in block_df else None,
        "acj": lambda_gc_from_p(block_df["acj_p"]) if "acj_p" in block_df else None,
        "omni": lambda_gc_from_p(block_df["omni_p"]) if "omni_p" in block_df else None,
        "joint_lrt": lambda_gc_from_p(block_df["joint_lrt_p"]) if "joint_lrt_p" in block_df else None,
        "joint_lowo": lambda_gc_from_p(block_df["joint_lowo_p"]) if "joint_lowo_p" in block_df else None,
        "beyond_lrt": lambda_gc_from_p(block_df["beyond_lrt_p"]) if "beyond_lrt_p" in block_df else None,
        "snp": lambda_gc_from_p(snp_df["snp_p"]) if "snp_p" in snp_df else None,
    }

    def _hits(frame, col):
        if frame is None or len(frame) == 0 or col not in frame.columns:
            return 0
        return int(np.sum(frame[col].to_numpy(dtype=float) < thr))

    meta = {
        "version": __version__,
        "n": int(len(y)),
        "n_snps": int(geno.n_snps),
        "n_blocks": int(len(block_df)),
        "n_loci": int(len(loci)) if loci is not None else 0,
        "bonferroni": thr,
        "kinship": [name for name, _n in nulls],
        "null": [
            {"kinship": name, "method": n.method, "vg": n.vg, "ve": n.ve,
             "h2": (n.vg / (n.vg + n.ve)) if (n.vg + n.ve) > 0 else None}
            for name, n in nulls
        ],
        "lambda_gc_raw_p": lambdas,
        "hits_blocks": {col: _hits(block_df, col) for col in ("omni_p", "joint_p", "acj_p", "allele_bonf_p", "snp_bonf_p")},
        "hits_loci": {col: _hits(loci, col) for col in ("omni_p", "joint_p", "acj_p", "allele_bonf_p", "snp_bonf_p")},
        "snp_hits": _hits(snp_df.assign(snp_p=snp_df["snp_p"]) if len(snp_df) else snp_df, "snp_p") if len(snp_df) else 0,
    }
    # SNP genome-wide threshold is 0.05/n_snps, not the block threshold.
    if len(snp_df):
        snp_thr = 0.05 / max(len(snp_df), 1)
        meta["snp_bonferroni"] = snp_thr
        meta["snp_hits"] = int(np.sum(snp_df["snp_p"].to_numpy(dtype=float) < snp_thr))
    (out / "run_meta.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")

    print(f"hapGWAS {__version__}  n={meta['n']}  snps={meta['n_snps']}  blocks={meta['n_blocks']}  loci={meta['n_loci']}")
    print(f"kinship: {', '.join(meta['kinship'])}")
    for item in meta["null"]:
        h2 = item["h2"]
        h2s = "NA" if h2 is None else f"{h2:.3f}"
        print(f"  {item['kinship']}: {item['method']}  h2={h2s}  vg={item['vg']:.4g}  ve={item['ve']:.4g}")
    print("lambda_GC on raw p (not Bonferroni): " + ", ".join(
        f"{k}={v:.3f}" if isinstance(v, float) and v == v else f"{k}=NA" for k, v in lambdas.items()
    ))
    print(f"block Bonferroni {thr:.3g}  hits (blocks / loci):")
    for col in ("omni_p", "joint_p", "acj_p", "allele_bonf_p", "snp_bonf_p"):
        print(f"  {col:16s}  {meta['hits_blocks'][col]:4d} / {meta['hits_loci'][col]:4d}")
    if len(snp_df):
        print(f"SNP GWAS hits at 0.05/n_snps: {meta['snp_hits']}")
    print(f"Wrote {out / 'assoc_blocks.tsv'}")
    print(f"Wrote {out / 'assoc_loci.tsv'}")
    print(f"Wrote {out / 'assoc_snps.tsv'}")
    print(f"Wrote {out / 'assoc_alleles.tsv'}")
    show = [
        "block_id", "omni_p", "joint_p", "joint_lrt_p", "joint_lowo_p", "beyond_p", "beyond_lrt_p",
        "tree_p", "allele_bonf_p", "snp_bonf_p",
    ]
    if len(block_df) and "omni_p" in block_df.columns:
        top = block_df.sort_values("omni_p").head(5)
        print(top[[c for c in show if c in top.columns]].to_string(index=False))
    return 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="hapgwas",
        description="Unified haplotype GWAS: block joint multi-df + HapAllele dosage tests",
    )
    p.add_argument("--version", action="version", version=f"hapgwas {__version__}")
    sub = p.add_subparsers(dest="cmd", required=True)

    d = sub.add_parser("demo", help="Run synthetic end-to-end demo")
    d.add_argument("--out", default="hapgwas_demo_out", help="Output directory")
    d.add_argument("--ols", action="store_true", help="Use OLS null instead of spectral LMM")
    d.set_defaults(func=_cmd_demo)

    c = sub.add_parser("convert", help="Convert VCF/haps/PED+MAP -> npz (or npz->vcf)")
    c.add_argument("--geno", required=True, help="input .vcf/.vcf.gz/.haps/.ped/.npz")
    c.add_argument("--map", default=None, help="PLINK .map companion for .ped")
    c.add_argument("--sample", default=None, help="Oxford .sample companion for .haps")
    c.add_argument("--out", required=True, help="output .npz or .vcf")
    c.set_defaults(func=_cmd_convert)

    r = sub.add_parser("run", help="Run hapGWAS on phased genotypes + phenotype")
    r.add_argument("--geno", required=True, help=".npz / phased .vcf / .haps / .ped")
    r.add_argument("--map", default=None, help="PLINK .map if --geno is .ped")
    r.add_argument("--sample", default=None, help="Oxford .sample if --geno is .haps")
    r.add_argument("--pheno", required=True, help="phenotype table with IID, PHENO")
    r.add_argument("--cov", default=None, help="optional covariates table")
    r.add_argument("--blocks", default=None, help="block TSV or PLINK --blocks file")
    r.add_argument(
        "--block-mode",
        dest="block_mode",
        choices=["window", "gabriel"],
        default="window",
        help="block partition when --blocks not given (default: window)",
    )
    r.add_argument("--gabriel-max-snps", dest="gabriel_max_snps", type=int, default=30)
    r.add_argument("--write-blocks", dest="write_blocks", default=None, help="write resolved blocks TSV")
    r.add_argument("--window", type=int, default=3, help="fixed window size if block-mode=window")
    r.add_argument("--step", type=int, default=None, help="window step (default=window)")
    r.add_argument("--chrom", default="1")
    r.add_argument("--min-freq", dest="min_freq", type=float, default=0.01)
    r.add_argument("--min-count", dest="min_count", type=int, default=5)
    r.add_argument(
        "--collapse",
        choices=["freq", "keep_mass", "none"],
        default="freq",
        help="rare haplotype collapse mode",
    )
    r.add_argument(
        "--condition-snp",
        dest="condition_snp",
        default=None,
        help="condition hap joint test on SNP (index / snp_id / chrom:pos)",
    )
    r.add_argument(
        "--condition-auto-topsnp",
        dest="condition_auto_topsnp",
        action="store_true",
        help="per block, condition on the strongest 1-df SNP inside the block",
    )
    r.add_argument("--ols", action="store_true", help="force OLS null (ignore relatedness)")
    r.add_argument(
        "--lowo",
        action="store_true",
        help="also report the region-excluded haplotype likelihood ratio (one GRM decomposition, Woodbury downdate per window)",
    )
    r.add_argument(
        "--lowo-flank",
        dest="lowo_flank",
        type=int,
        default=20,
        help="SNPs either side of the window to remove from the GRM with the window",
    )
    r.add_argument(
        "--grm",
        choices=["auto", "vanraden", "loco", "provided", "none"],
        default="auto",
        help="auto: use K in the genotype file, otherwise a VanRaden GRM (LOCO when one file has 2+ chromosomes)",
    )
    r.add_argument(
        "--gcta-grm",
        dest="gcta_grm",
        default=None,
        help="GCTA --make-grm prefix (.grm.id + .grm.bin). IIDs must match the phenotype. Not compatible with --lowo.",
    )
    r.add_argument("--kin-geno", default=None, help="phased file used only to build the GRM (leave the test chromosome out)")
    r.add_argument("--kin-map", default=None, help="PLINK .map for --kin-geno if it is a .ped")
    r.add_argument("--kin-sample", default=None, help="Oxford .sample for --kin-geno if it is .haps")
    r.add_argument("--out", default="hapgwas_out")
    r.set_defaults(func=_cmd_run)
    return p


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return int(args.func(args))


if __name__ == "__main__":
    sys.exit(main())
