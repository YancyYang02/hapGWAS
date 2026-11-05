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
from hapgwas.io import align_ids, load_grm, read_covariates, read_phenotype, write_assoc, write_phased_vcf
from hapgwas.joint_test import run_genome, score_snp
from hapgwas.lmm import fit_null_lmm

# All block-level scores produced in one run (integrated tool)
ALL_METHODS = ("joint", "acj", "sacj", "allele_bonf", "snp_bonf", "pca", "ridge")
METHOD_HELP = (
    "Comma-separated methods to report in assoc_blocks.tsv. "
    "Always computed together under one null; this only controls which "
    "hit/summary columns are emphasized. Choices: "
    + ",".join(ALL_METHODS)
    + " (default: all)."
)


def _parse_methods(s: str | None) -> list[str]:
    if not s or s.strip().lower() in {"all", "*"}:
        return list(ALL_METHODS)
    out = []
    for part in s.split(","):
        m = part.strip().lower().replace("-", "_")
        aliases = {"ghap": "allele_bonf", "ghap_like": "allele_bonf", "allelebonf": "allele_bonf",
                   "snp": "snp_bonf", "snpbonf": "snp_bonf"}
        m = aliases.get(m, m)
        if m not in ALL_METHODS:
            raise SystemExit(f"Unknown method {part!r}; choose from {ALL_METHODS}")
        if m not in out:
            out.append(m)
    return out



def _attach_grm(geno, args):
    """Attach external GRM/kinship to geno if --grm given; else keep geno.K."""
    path = getattr(args, "grm", None)
    if not path:
        return geno
    if getattr(args, "ols", False):
        print("WARNING: --grm ignored because --ols was set")
        return geno
    K = load_grm(path, sample_ids=list(geno.ids))
    geno.K = K
    geno.meta = dict(geno.meta or {})
    geno.meta["grm_path"] = str(path)
    geno.meta["grm_n"] = int(K.shape[0])
    print(f"Attached GRM {path} shape={K.shape} diag_mean={float(K.diagonal().mean()):.4f}")
    return geno


def _resolve_blocks(args, geno: PhasedGeno):
    n_snps = geno.n_snps
    positions = geno.positions
    if getattr(args, "blocks", None):
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


def enrich_block_baselines(
    block_df: pd.DataFrame,
    allele_df: pd.DataFrame,
    null,
    geno: PhasedGeno,
    methods: list[str],
) -> pd.DataFrame:
    """Add allele-Bonferroni (GHap-like) and SNP-Bonferroni columns into block table."""
    out = block_df.copy()
    if "allele_bonf" in methods and len(allele_df):
        gh = []
        for bid, sub in allele_df.groupby("block_id"):
            ps = pd.to_numeric(sub["marg_p"], errors="coerce").dropna()
            if len(ps) == 0:
                continue
            m = float(ps.min())
            gh.append({"block_id": bid, "allele_bonf_p": min(1.0, m * len(ps)), "allele_min_p": m})
        if gh:
            out = out.merge(pd.DataFrame(gh), on="block_id", how="left")
            # alias used in papers / case notes
            out["ghap_like_bon_p"] = out["allele_bonf_p"]
    if "snp_bonf" in methods and len(out):
        snp_ps = []
        for j in range(geno.n_snps):
            dos = geno.hap1[:, j].astype(float) + geno.hap2[:, j].astype(float)
            _, p = score_snp(null, dos)
            snp_ps.append(p if np.isfinite(p) else np.nan)
        snp_ps = np.asarray(snp_ps, dtype=float)
        sb = []
        for _, r in out.iterrows():
            s0, s1 = int(r.snp_start), int(r.snp_end)
            sub = snp_ps[s0:s1]
            sub = sub[np.isfinite(sub)]
            if len(sub) == 0:
                sb.append(np.nan)
                continue
            m = float(np.min(sub))
            sb.append(min(1.0, m * len(sub)))
        out["snp_bonf_p"] = sb
    return out


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
    geno = _attach_grm(geno, args)
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


def _run_one_trait(
    geno: PhasedGeno,
    pheno: pd.DataFrame,
    cov,
    args,
    methods: list[str],
    blocks,
) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    y, keep, X, kept_ids = align_ids(pheno, geno.ids, cov)
    g = geno.subset_samples(keep)
    null = fit_null_lmm(
        y,
        X=X,
        K=None if args.ols else g.K,
        method="ols" if args.ols else "auto",
    )
    block_df, allele_df = run_genome(
        null,
        blocks,
        g.hap1,
        g.hap2,
        min_freq=args.min_freq,
        min_count=args.min_count,
        collapse_mode=args.collapse,
        condition_snp=None,
        condition_auto_topsnp=args.condition_auto_topsnp,
        snp_ids=g.snp_ids,
    )
    block_df = enrich_block_baselines(block_df, allele_df, null, g, methods)
    n_blocks = len(block_df)
    thr = 0.05 / max(n_blocks, 1)
    colmap = {
        "joint": "joint_p",
        "acj": "acj_p",
        "sacj": "sacj_p",
        "allele_bonf": "allele_bonf_p",
        "snp_bonf": "snp_bonf_p",
        "pca": "pca_p",
        "ridge": "ridge_p",
    }
    hits, mins = {}, {}
    for m in methods:
        col = colmap[m]
        if col not in block_df.columns:
            hits[m] = 0
            mins[m] = float("nan")
            continue
        ps = pd.to_numeric(block_df[col], errors="coerce")
        hits[m] = int((ps < thr).sum())
        mins[m] = float(ps.min()) if ps.notna().any() else float("nan")
    summary = {
        "n": int(len(y)),
        "n_snps": int(g.n_snps),
        "n_blocks": n_blocks,
        "thr": thr,
        "hits": hits,
        "min_p": mins,
        "methods": methods,
        "null_method": null.method,
        "vg": float(null.vg),
        "ve": float(null.ve),
        "has_K": g.K is not None and not bool(getattr(args, "ols", False)),
    }
    return block_df, allele_df, summary


def _cmd_run(args: argparse.Namespace) -> int:
    methods = _parse_methods(getattr(args, "methods", None))
    geno = load_geno(args.geno, companion=args.map or args.sample)
    geno = _attach_grm(geno, args)
    pheno = read_phenotype(args.pheno, pheno_col=args.pheno_col)
    cov = read_covariates(args.cov) if args.cov else None
    blocks = _resolve_blocks(args, geno)
    if args.write_blocks:
        summarize_blocks(blocks).to_csv(args.write_blocks, sep="\t", index=False)

    if args.condition_snp is not None:
        # keep legacy single-SNP conditioning path via re-call below
        pass

    block_df, allele_df, summary = _run_one_trait(geno, pheno, cov, args, methods, blocks)

    # optional legacy condition_snp (re-run with conditioning on joint only)
    if args.condition_snp is not None:
        cond_snp = geno.find_snp(args.condition_snp)
        y, keep, X, _ = align_ids(pheno, geno.ids, cov)
        g = geno.subset_samples(keep)
        null = fit_null_lmm(y, X=X, K=None if args.ols else g.K, method="ols" if args.ols else "auto")
        block_df, allele_df = run_genome(
            null, blocks, g.hap1, g.hap2,
            min_freq=args.min_freq, min_count=args.min_count, collapse_mode=args.collapse,
            condition_snp=cond_snp, condition_auto_topsnp=False, snp_ids=g.snp_ids,
        )
        block_df = enrich_block_baselines(block_df, allele_df, null, g, methods)

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    write_assoc(block_df, out / "assoc_blocks.tsv")
    write_assoc(allele_df, out / "assoc_alleles.tsv")
    (out / "run_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(f"Wrote {out / 'assoc_blocks.tsv'} ({len(block_df)} blocks)")
    print(f"Wrote {out / 'assoc_alleles.tsv'} ({len(allele_df)} allele rows)")
    print(f"methods={methods} thr={summary['thr']:.3g} hits={summary['hits']}")
    cols = ["block_id", "snp_start", "snp_end", "joint_df", "joint_p", "acj_p", "sacj_p",
            "allele_bonf_p", "snp_bonf_p", "pca_p", "ridge_p"]
    if len(block_df):
        top = block_df.sort_values("joint_p").head(5)
        print(top[[c for c in cols if c in top.columns]].to_string(index=False))
    return 0


def _cmd_run_multi(args: argparse.Namespace) -> int:
    """Multi-trait batch: one geno, wide phenotype table (IID + many trait columns)."""
    methods = _parse_methods(getattr(args, "methods", None))
    geno = load_geno(args.geno, companion=args.map or args.sample)
    geno = _attach_grm(geno, args)
    path = Path(args.pheno)
    sep = "," if path.suffix.lower() == ".csv" else r"\s+"
    wide = pd.read_csv(path, sep=sep, engine="python")
    cols_l = {c.lower(): c for c in wide.columns}
    id_col = args.id_col
    if id_col not in wide.columns and id_col.lower() in cols_l:
        id_col = cols_l[id_col.lower()]
    if id_col not in wide.columns:
        raise SystemExit(f"Need id column {args.id_col!r} in {path}")
    wide = wide.rename(columns={id_col: "IID"})
    wide["IID"] = wide["IID"].astype(str)
    skip = {"iid", "fid", "id", "family", "pid", "mid", "sex"}
    if args.traits:
        traits = [t.strip() for t in args.traits.split(",") if t.strip()]
    else:
        traits = [c for c in wide.columns if c != "IID" and c.lower() not in skip]
    if args.max_traits and args.max_traits > 0:
        traits = traits[: args.max_traits]
    cov = read_covariates(args.cov) if args.cov else None
    blocks = _resolve_blocks(args, geno)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    rows = []
    for i, trait in enumerate(traits, 1):
        sub = wide[["IID", trait]].dropna().rename(columns={trait: "PHENO"})
        # numeric only
        sub["PHENO"] = pd.to_numeric(sub["PHENO"], errors="coerce")
        sub = sub.dropna(subset=["PHENO"])
        if len(sub) < 50:
            print(f"[{i}/{len(traits)}] SKIP {trait} n={len(sub)}")
            continue
        print(f"[{i}/{len(traits)}] {trait} n={len(sub)}", flush=True)
        trait_out = out / "by_trait" / str(trait).replace("/", "_")
        trait_out.mkdir(parents=True, exist_ok=True)
        block_df, allele_df, summary = _run_one_trait(geno, sub, cov, args, methods, blocks)
        write_assoc(block_df, trait_out / "assoc_blocks.tsv")
        write_assoc(allele_df, trait_out / "assoc_alleles.tsv")
        (trait_out / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
        row = {"trait": trait, "n": summary["n"], "n_blocks": summary["n_blocks"], "thr": summary["thr"]}
        for m in methods:
            row[f"hits_{m}"] = summary["hits"].get(m, 0)
            row[f"min_{m}"] = summary["min_p"].get(m, float("nan"))
        # convenience flags
        row["sacj_beats_snp"] = int(row.get("hits_sacj", 0) > row.get("hits_snp_bonf", 0))
        row["acj_beats_snp"] = int(row.get("hits_acj", 0) > row.get("hits_snp_bonf", 0))
        rows.append(row)
    summary_df = pd.DataFrame(rows)
    summary_df.to_csv(out / "multi_trait_hits.tsv", sep="\t", index=False)
    meta = {
        "n_traits_attempted": len(traits),
        "n_traits_done": len(rows),
        "n_snps": int(geno.n_snps),
        "n_samples_geno": int(geno.n_samples),
        "methods": methods,
        "window": args.window,
        "step": args.step or args.window,
        "ols": bool(args.ols),
    }
    (out / "multi_meta.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")
    print(f"Wrote {out / 'multi_trait_hits.tsv'} ({len(rows)} traits)")
    if len(summary_df):
        print(summary_df.head(20).to_string(index=False))
    return 0


def _add_run_shared(r: argparse.ArgumentParser) -> None:
    r.add_argument("--geno", required=True, help=".npz / phased .vcf / .haps / .ped")
    r.add_argument("--map", default=None, help="PLINK .map if --geno is .ped")
    r.add_argument("--sample", default=None, help="Oxford .sample if --geno is .haps")
    r.add_argument("--pheno", required=True, help="phenotype table")
    r.add_argument("--cov", default=None, help="optional covariates table")
    r.add_argument("--blocks", default=None, help="block TSV or PLINK --blocks file")
    r.add_argument(
        "--block-mode", dest="block_mode", choices=["window", "gabriel"], default="window",
        help="block partition when --blocks not given (default: window)",
    )
    r.add_argument("--gabriel-max-snps", dest="gabriel_max_snps", type=int, default=30)
    r.add_argument("--write-blocks", dest="write_blocks", default=None, help="write resolved blocks TSV")
    r.add_argument("--window", type=int, default=3, help="fixed window size if block-mode=window")
    r.add_argument("--step", type=int, default=None, help="window step (default=window)")
    r.add_argument("--chrom", default="1")
    r.add_argument("--min-freq", dest="min_freq", type=float, default=0.01)
    r.add_argument("--min-count", dest="min_count", type=int, default=5)
    r.add_argument("--collapse", choices=["freq", "keep_mass", "none"], default="freq")
    r.add_argument("--grm", "--kinship", dest="grm", default=None,
                   help="kinship/GRM path (.npy/.npz/text/GCTA .grm.bin+.id); enables spectral LMM")
    r.add_argument("--ols", action="store_true", help="force OLS null (ignore K)")
    r.add_argument("--methods", default="all", help=METHOD_HELP)
    r.add_argument("--out", default="hapgwas_out")


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="hapgwas",
        description=(
            "Unified haplotype GWAS: joint + ACJ + SACJ + allele-Bonf (GHap-like) "
            "+ SNP-Bonf + optional PCA/ridge in one run"
        ),
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

    r = sub.add_parser("run", help="Single-trait hapGWAS (all methods → one table)")
    _add_run_shared(r)
    r.add_argument("--pheno-col", dest="pheno_col", default="PHENO", help="phenotype column name")
    r.add_argument("--condition-snp", dest="condition_snp", default=None)
    r.add_argument("--condition-auto-topsnp", dest="condition_auto_topsnp", action="store_true")
    r.set_defaults(func=_cmd_run)

    m = sub.add_parser("run-multi", help="Multi-trait batch: wide pheno + all methods per trait")
    _add_run_shared(m)
    m.add_argument("--id-col", dest="id_col", default="IID")
    m.add_argument("--traits", default=None, help="comma-separated trait columns (default: all numeric non-ID)")
    m.add_argument("--max-traits", dest="max_traits", type=int, default=0, help="cap traits (0=all)")
    m.add_argument("--condition-auto-topsnp", dest="condition_auto_topsnp", action="store_true")
    m.set_defaults(func=_cmd_run_multi)
    return p


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return int(args.func(args))


if __name__ == "__main__":
    sys.exit(main())
