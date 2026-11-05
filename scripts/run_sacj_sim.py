#!/usr/bin/env python3
"""Type I + power for SACJ vs ACJ / joint / allele-Bonf / SNP-Bonf."""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
import pandas as pd

from hapgwas.adaptive_score import multi_allelic_scores
from hapgwas.dosage import build_dosage, haplotypes_from_phased_snps
from hapgwas.joint_test import score_snp
from hapgwas.lmm import fit_null_lmm

OUT = Path(__file__).resolve().parents[1] / "sim_results"
OUT.mkdir(parents=True, exist_ok=True)

N_REP = 300
N = 500
SEED = 7
ALPHA = 0.05
BETA_SINGLE = 0.40
BETA_NEIGHBOR = 0.24
BETA_DIST2 = 0.28  # two alleles at Hamming distance 2 (SACJ target regime)


def make_K(n, rank=12, seed=0):
    rng = np.random.default_rng(seed)
    Z = rng.normal(size=(n, rank))
    K = Z @ Z.T / rank
    return K / np.mean(np.diag(K))


def phenotype(n, effect, K, h2, seed):
    rng = np.random.default_rng(seed)
    evals, evecs = np.linalg.eigh(K)
    evals = np.clip(evals, 1e-8, None)
    u = evecs @ (np.sqrt(evals) * rng.normal(size=n))
    u = u / np.std(u) * np.sqrt(h2)
    e = rng.normal(0, np.sqrt(max(1 - h2, 1e-6)), size=n)
    base = effect if effect is not None else 0.0
    return base + u + e


def sim_neighbor_block(n, seed):
    rng = np.random.default_rng(seed)
    patterns = [
        (0, 0, 0, 0),
        (0, 0, 0, 1),
        (0, 0, 1, 0),
        (0, 1, 0, 0),
        (1, 0, 0, 0),
        (1, 1, 0, 0),
    ]
    w = np.array([0.22, 0.18, 0.15, 0.15, 0.15, 0.15], float); w /= w.sum()
    i1 = rng.choice(len(patterns), size=n, p=w)
    i2 = rng.choice(len(patterns), size=n, p=w)
    hap1 = np.asarray([patterns[i] for i in i1], dtype=np.int8)
    hap2 = np.asarray([patterns[i] for i in i2], dtype=np.int8)
    return hap1, hap2, patterns


def sim_dist2_block(n, seed):
    """Causal alleles at Hamming distance 2: 0000 and 0011 — ACJ (h=1) cannot merge; SACJ may."""
    rng = np.random.default_rng(seed)
    patterns = [
        (0, 0, 0, 0),
        (0, 0, 1, 1),  # d=2 from first
        (0, 0, 0, 1),
        (0, 0, 1, 0),
        (1, 0, 0, 0),
        (0, 1, 0, 0),
    ]
    w = np.array([0.20, 0.20, 0.15, 0.15, 0.15, 0.15], float); w /= w.sum()
    i1 = rng.choice(len(patterns), size=n, p=w)
    i2 = rng.choice(len(patterns), size=n, p=w)
    hap1 = np.asarray([patterns[i] for i in i1], dtype=np.int8)
    hap2 = np.asarray([patterns[i] for i in i2], dtype=np.int8)
    return hap1, hap2, patterns


def allele_labels_from_patterns(patterns):
    return ["".join(str(int(x)) for x in p) for p in patterns]


def block_pvals(null, hap1, hap2):
    pairs = haplotypes_from_phased_snps(hap1, hap2)
    hap = build_dosage(pairs, min_freq=0.01, min_count=5, collapse_mode="freq", drop_ref=True)
    ms = multi_allelic_scores(null, hap.dosage, allele_labels=hap.alleles)
    # allele bonf
    from hapgwas.joint_test import score_joint_and_marginal
    _, _, _, _, marg_p = score_joint_and_marginal(null, hap.dosage)
    mp = marg_p[np.isfinite(marg_p)]
    allele_bon = float(min(1.0, mp.min() * len(mp))) if len(mp) else np.nan
    # snp bonf
    snp_ps = []
    for j in range(hap1.shape[1]):
        dos = hap1[:, j].astype(float) + hap2[:, j].astype(float)
        _, p = score_snp(null, dos)
        if np.isfinite(p):
            snp_ps.append(p)
    snp_bon = float(min(1.0, min(snp_ps) * len(snp_ps))) if snp_ps else np.nan
    return {
        "joint_p": ms.full_p,
        "acj_p": ms.acj.p,
        "sacj_p": ms.sacj.p,
        "pca_p": ms.pca.p,
        "ridge_p": ms.ridge.p,
        "allele_bon_p": allele_bon,
        "snp_bon_p": snp_bon,
        "sacj_k": ms.sacj.k_spectral,
        "sacj_ngroups": ms.sacj.n_groups,
        "acj_ngroups": ms.acj.n_groups,
    }


def dosage_effect(hap1, hap2, patterns, causal_idx, beta):
    """Additive effect from listed pattern indices (haplotype copies)."""
    lab = {tuple(p): i for i, p in enumerate(patterns)}
    n = hap1.shape[0]
    eff = np.zeros(n)
    for i in range(n):
        for h in (hap1[i], hap2[i]):
            t = tuple(int(x) for x in h)
            if t in lab and lab[t] in causal_idx:
                eff[i] += beta
    return eff


def run_suite(name, block_fn, causal_builder, n_rep=N_REP, beta=None, null_only=False):
    rows = []
    K = make_K(N, seed=SEED)
    for r in range(n_rep):
        hap1, hap2, patterns = block_fn(N, seed=SEED + 1000 + r)
        if null_only:
            effect = np.zeros(N)
        else:
            cidx = causal_builder(patterns)
            effect = dosage_effect(hap1, hap2, patterns, cidx, beta)
        y = phenotype(N, effect, K, h2=0.3, seed=SEED + 5000 + r)
        null = fit_null_lmm(y, X=None, K=K, method="lmm")
        st = block_pvals(null, hap1, hap2)
        st["rep"] = r
        rows.append(st)
    df = pd.DataFrame(rows)
    summary = {"scenario": name, "n_rep": n_rep, "n": N}
    for col in ["joint_p", "acj_p", "sacj_p", "pca_p", "ridge_p", "allele_bon_p", "snp_bon_p"]:
        p = df[col].to_numpy(float)
        summary[f"rate_{col}"] = float(np.mean(p < ALPHA))
        summary[f"median_{col}"] = float(np.median(p))
    summary["mean_sacj_k"] = float(df["sacj_k"].mean())
    summary["mean_sacj_ngroups"] = float(df["sacj_ngroups"].mean())
    summary["mean_acj_ngroups"] = float(df["acj_ngroups"].mean())
    return df, summary


def main():
    # null on neighbor design
    d0, s0 = run_suite("null_neighbor", sim_neighbor_block, lambda p: set(), null_only=True)
    # single allele: first pattern
    d1, s1 = run_suite(
        "power_single",
        sim_neighbor_block,
        lambda patterns: {0},
        beta=BETA_SINGLE,
    )
    # hamming-1 neighbors: patterns 0 and 1
    d2, s2 = run_suite(
        "power_neighbor_h1",
        sim_neighbor_block,
        lambda patterns: {0, 1},
        beta=BETA_NEIGHBOR,
    )
    # hamming-2 pair: patterns 0 and 1 in dist2 design
    d3, s3 = run_suite(
        "power_dist2",
        sim_dist2_block,
        lambda patterns: {0, 1},
        beta=BETA_DIST2,
    )
    # null dist2 design
    d4, s4 = run_suite("null_dist2", sim_dist2_block, lambda p: set(), null_only=True)

    summaries = [s0, s4, s1, s2, s3]
    sum_df = pd.DataFrame(summaries)
    sum_df.to_csv(OUT / "sacj_sim_summary.tsv", sep="\t", index=False)
    (OUT / "sacj_sim_summary.json").write_text(json.dumps(summaries, indent=2))
    d0.to_csv(OUT / "sacj_null_neighbor_reps.tsv", sep="\t", index=False)
    d3.to_csv(OUT / "sacj_power_dist2_reps.tsv", sep="\t", index=False)

    lines = ["# SACJ simulation (v0.4)\n\n", f"n={N}, n_rep={N_REP}, LMM h2=0.3, α={ALPHA}\n\n"]
    lines.append("| scenario | joint | ACJ | **SACJ** | alleleBon | snpBon |\n|---|---:|---:|---:|---:|---:|\n")
    for s in summaries:
        lines.append(
            f"| {s['scenario']} | {s['rate_joint_p']:.3f} | {s['rate_acj_p']:.3f} | "
            f"**{s['rate_sacj_p']:.3f}** | {s['rate_allele_bon_p']:.3f} | {s['rate_snp_bon_p']:.3f} |\n"
        )
    lines.append("\nSACJ vs ACJ: ACJ merges Hamming≤1 components; SACJ soft Hamming affinity + spectral clusters (eigengap k).\n")
    lines.append("Target regime for novelty: **power_dist2** (causal alleles at Hamming distance 2).\n")
    (OUT / "SACJ_RESULTS.md").write_text("".join(lines), encoding="utf-8")
    print("".join(lines))
    print(json.dumps(summaries, indent=2))


if __name__ == "__main__":
    main()
