"""Synthetic end-to-end demo for hapGWAS."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from hapgwas.blocks import fixed_windows, summarize_blocks
from hapgwas.genotype import PhasedGeno, save_npz
from hapgwas.io import write_assoc, write_phased_vcf
from hapgwas.joint_test import run_genome
from hapgwas.lmm import fit_null_lmm


def simulate_phased(
    n: int = 400,
    n_snps: int = 60,
    causal_block: int = 5,
    window: int = 3,
    beta: float = 0.6,
    h2_kinship: float = 0.3,
    seed: int = 42,
) -> dict:
    """
    Simulate phased SNP haplotypes with one causal haplotype allele in a block.
    Returns dict with hap1, hap2, y, ids, positions, meta.
    """
    rng = np.random.default_rng(seed)
    # ancestral haplotype frequencies: create diversity via Bernoulli SNPs with LD in blocks
    maf = rng.uniform(0.1, 0.4, size=n_snps)
    # build hap1/hap2 with modest LD within windows of `window`
    hap1 = np.zeros((n, n_snps), dtype=np.int8)
    hap2 = np.zeros((n, n_snps), dtype=np.int8)
    for start in range(0, n_snps, window):
        end = min(start + window, n_snps)
        # latent haplotype type per chromosome copy
        for mat in (hap1, hap2):
            types = rng.integers(0, 4, size=n)
            for t in range(4):
                pattern = (rng.random(end - start) < maf[start:end]).astype(np.int8)
                # flip a bit of noise
                idx = np.where(types == t)[0]
                if len(idx) == 0:
                    continue
                base = np.tile(pattern, (len(idx), 1))
                flip = rng.random(base.shape) < 0.05
                base[flip] = 1 - base[flip]
                mat[idx, start:end] = base

    # causal: specific haplotype string in block causal_block
    s = causal_block * window
    e = s + window
    if e > n_snps:
        s, e = 0, window
    # define causal allele as all-zero haplotype in that window (may be rare-ish)
    # inject: make ~15% chromosomes carry 000.. and give them phenotypic effect
    causal_pattern = np.zeros(e - s, dtype=np.int8)
    carry = rng.random(n) < 0.25
    hap1[carry, s:e] = causal_pattern
    carry2 = rng.random(n) < 0.25
    hap2[carry2, s:e] = causal_pattern
    dosage_causal = carry.astype(float) + carry2.astype(float)  # 0/1/2
    causal_snp_start = int(s)
    causal_snp_end = int(e)

    # kinship-ish GRM from random SNPs for null LMM demo
    # simple: low-rank relatedness
    Z = rng.normal(size=(n, 20))
    K = Z @ Z.T / 20.0
    # standardize K
    K = K / np.mean(np.diag(K))

    # phenotype
    u = rng.multivariate_normal(np.zeros(n), h2_kinship * K + 1e-6 * np.eye(n))
    noise = rng.normal(0, np.sqrt(max(1 - h2_kinship, 1e-3)), size=n)
    y = beta * (dosage_causal - dosage_causal.mean()) + u + noise

    ids = [f"S{i:04d}" for i in range(n)]
    positions = np.arange(1, n_snps + 1) * 1000
    return {
        "hap1": hap1,
        "hap2": hap2,
        "y": y,
        "ids": ids,
        "K": K,
        "positions": positions,
        "causal_snp_start": causal_snp_start,
        "causal_snp_end": causal_snp_end,
        "window": window,
        "dosage_causal": dosage_causal,
    }


def run_demo(out_dir: str | Path, use_lmm: bool = True) -> dict:
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    sim = simulate_phased()
    hap1, hap2 = sim["hap1"], sim["hap2"]
    n, n_snps = hap1.shape
    y = sim["y"]
    X = np.ones((n, 1))
    K = sim["K"] if use_lmm else None
    null = fit_null_lmm(y, X=X, K=K, method="auto")

    blocks = fixed_windows(n_snps, window=sim["window"], step=sim["window"], chrom="1", positions=sim["positions"])
    summarize_blocks(blocks).to_csv(out_dir / "blocks.tsv", sep="\t", index=False)

    block_df, allele_df = run_genome(
        null,
        blocks,
        hap1,
        hap2,
        min_freq=0.01,
        min_count=5,
        collapse_mode="freq",
    )
    write_assoc(block_df, out_dir / "assoc_blocks.tsv")
    write_assoc(allele_df, out_dir / "assoc_alleles.tsv")

    # phenotype table for reuse
    pd.DataFrame({"IID": sim["ids"], "PHENO": y}).to_csv(out_dir / "pheno.tsv", sep="\t", index=False)

    # save phased arrays + VCF for import smoke tests
    snp_ids = np.asarray([f"snp{i}" for i in range(n_snps)], dtype=object)
    chrom = np.asarray(["1"] * n_snps, dtype=object)
    g = PhasedGeno(
        hap1=hap1,
        hap2=hap2,
        ids=list(sim["ids"]),
        chrom=chrom,
        positions=sim["positions"],
        snp_ids=snp_ids,
        K=sim["K"],
        source="demo",
    )
    save_npz(g, out_dir / "phased_demo.npz")
    write_phased_vcf(g, out_dir / "phased_demo.vcf")

    # also write conditioned results (auto top SNP per block)
    block_df_c, allele_df_c = run_genome(
        null,
        blocks,
        hap1,
        hap2,
        min_freq=0.01,
        min_count=5,
        collapse_mode="freq",
        condition_auto_topsnp=True,
        snp_ids=snp_ids,
    )
    write_assoc(block_df_c, out_dir / "assoc_blocks_cond_auto.tsv")
    write_assoc(allele_df_c, out_dir / "assoc_alleles_cond_auto.tsv")

    top = block_df.sort_values("joint_p").head(5)
    causal = block_df[
        (block_df["snp_start"] == sim["causal_snp_start"])
        & (block_df["snp_end"] == sim["causal_snp_end"])
    ]
    summary = {
        "n": n,
        "n_snps": n_snps,
        "n_blocks": len(blocks),
        "null_method": null.method,
        "vg": null.vg,
        "ve": null.ve,
        "causal_block_snp_start": sim["causal_snp_start"],
        "causal_block_joint_p": float(causal["joint_p"].iloc[0]) if len(causal) else None,
        "top_blocks": top[["block_id", "snp_start", "snp_end", "joint_p", "joint_df"]].to_dict("records"),
        "out_dir": str(out_dir.resolve()),
        "cond_causal_joint_p": float(block_df_c.loc[(block_df_c["snp_start"]==sim["causal_snp_start"]) & (block_df_c["snp_end"]==sim["causal_snp_end"]), "joint_p"].iloc[0]) if ((block_df_c["snp_start"]==sim["causal_snp_start"]) & (block_df_c["snp_end"]==sim["causal_snp_end"])).any() else None,
        "cond_causal_cond_joint_p": float(block_df_c.loc[(block_df_c["snp_start"]==sim["causal_snp_start"]) & (block_df_c["snp_end"]==sim["causal_snp_end"]), "cond_joint_p"].iloc[0]) if "cond_joint_p" in block_df_c.columns and ((block_df_c["snp_start"]==sim["causal_snp_start"]) & (block_df_c["snp_end"]==sim["causal_snp_end"])).any() else None,
    }
    pd.Series(summary).to_json(out_dir / "demo_summary.json")
    return summary
