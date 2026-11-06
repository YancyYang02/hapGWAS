"""Real broiler haplotypes: full-panel GRM versus region-excluded likelihood.

The GRM is the VanRaden matrix of the scanned panel. joint_lowo_p removes
the current window and --flank SNPs on each side. This is the real-data
counterpart of sim_results/lowo_bench.json.
"""

from __future__ import annotations

import json
import time
from pathlib import Path

import pandas as pd

from hapgwas.blocks import fixed_windows
from hapgwas.genotype import load_geno
from hapgwas.io import align_ids, write_assoc
from hapgwas.joint_test import run_genome
from hapgwas.lmm import fit_null_lmm
from hapgwas.omnibus import clump_overlapping, lambda_gc_from_p
from hapgwas.window_kinship import prepare_window_basis

PANELS = [
    ("chr4_1500", Path("/fs5/home/clwen/YancongYang/hapgwas_try/demo_data/phased_chr4_1p5k.vcf.gz")),
    ("first5000", Path("/fs5/home/clwen/YancongYang/hapgwas_try/data/phased_first5000.vcf.gz")),
]
PHENO = Path("/fs5/home/clwen/YancongYang/hapgwas_try/data/pheno_key.tsv")
OUT = Path("/fs5/home/clwen/YancongYang/hapgwas_try/runs/lowo_real")
TRAITS = ["BW6", "AFWR", "BW3_6", "EBW6", "GBW6", "ERFI"]
WINDOW, STEP, FLANK = 4, 2, 20
PCOLS = ("omni_p", "joint_p", "joint_lrt_p", "joint_lowo_p", "beyond_p", "snp_bonf_p", "allele_bonf_p")


def _lambda(p: pd.Series) -> float:
    arr = pd.to_numeric(p, errors="coerce").to_numpy(dtype=float)
    # Raw multi-df p-values are mapped to a 1-df chi-square, same as the CLI.
    return float(lambda_gc_from_p(arr))


def _one(geno, y, blocks, trait: str, panel: str) -> dict:
    t0 = time.time()
    from hapgwas.grm import vanraden_factors

    _Z, _keep, K = vanraden_factors(geno.snp_dosage_matrix(), maf_min=0.01)
    null = fit_null_lmm(y, K=K, method="lmm")
    basis = prepare_window_basis(y, geno.snp_dosage_matrix(), maf_min=0.01)
    block_df, allele_df = run_genome(
        null,
        blocks,
        geno.hap1,
        geno.hap2,
        min_freq=0.01,
        min_count=5,
        lowo_basis=basis,
        lowo_flank=FLANK,
        n_snps=geno.n_snps,
        apply_fdr=True,
    )
    out = OUT / panel / trait
    out.mkdir(parents=True, exist_ok=True)
    write_assoc(block_df, out / "assoc_blocks.tsv")
    write_assoc(allele_df, out / "assoc_alleles.tsv")
    loci = clump_overlapping(block_df, "joint_lowo_p") if len(block_df) else block_df
    write_assoc(loci, out / "assoc_loci_lowo.tsv")
    n_blocks = max(len(block_df), 1)
    thr = 0.05 / n_blocks
    summary = {
        "panel": panel,
        "trait": trait,
        "n": int(null.n),
        "n_snps": int(geno.n_snps),
        "n_blocks": int(len(block_df)),
        "n_loci_lowo": int(len(loci)),
        "h2": float(null.vg / (null.vg + null.ve)) if (null.vg + null.ve) > 0 else None,
        "vg": float(null.vg),
        "ve": float(null.ve),
        "threshold": thr,
        "flank": FLANK,
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
        top = loci.sort_values("joint_lowo_p").head(5)
        keep = [c for c in ("block_id", "chrom", "snp_start", "snp_end", "joint_lowo_p", "joint_p", "joint_lrt_p", "omni_p", "snp_bonf_p", "allele_bonf_p") if c in top.columns]
        summary["top_loci"] = top[keep].to_dict(orient="records")
    (out / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(trait, panel, "h2", summary["h2"], {c: summary[c]["hits_blocks"] for c in PCOLS if c in summary}, "lam", {c: round(summary[c]["lambda_gc"], 3) for c in ("joint_p", "joint_lowo_p") if c in summary}, flush=True)
    return summary


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    pheno = pd.read_csv(PHENO, sep=None, engine="python")
    pheno["IID"] = pheno["IID"].astype(str)
    rows = []
    for panel, path in PANELS:
        geno0 = load_geno(path)
        chrom = str(geno0.chrom[0]) if geno0.chrom is not None else "4"
        blocks = fixed_windows(geno0.n_snps, window=WINDOW, step=STEP, chrom=chrom, positions=geno0.positions)
        print(f"PANEL {panel} n={geno0.n_samples} snps={geno0.n_snps} blocks={len(blocks)}", flush=True)
        for trait in TRAITS:
            if trait not in pheno.columns:
                print("skip", trait, flush=True)
                continue
            sub = pheno[["IID", trait]].dropna().rename(columns={trait: "PHENO"})
            y, keep, _X, _ids = align_ids(sub, geno0.ids, None)
            geno = geno0.subset_samples(keep)
            rows.append(_one(geno, y, blocks, trait, panel))
    text = json.dumps(rows, indent=2)
    (OUT / "lowo_real.json").write_text(text, encoding="utf-8")
    print("WROTE", OUT / "lowo_real.json", flush=True)


if __name__ == "__main__":
    main()
