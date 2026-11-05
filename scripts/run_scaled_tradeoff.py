#!/usr/bin/env python3
"""Scaled sims + single-allele vs neighbor-allele power tradeoff figure."""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats

from hapgwas.adaptive_score import multi_allelic_scores
from hapgwas.dosage import build_dosage, haplotypes_from_phased_snps
from hapgwas.joint_test import score_snp
from hapgwas.lmm import fit_null_lmm, project_features

OUT = Path("/workspace/hapgwas/sim_results")
OUT.mkdir(parents=True, exist_ok=True)

# Scaled defaults
N_REP = 500
N = 600
SEED = 42
BETA_SINGLE = 0.40
BETA_NEIGHBOR = 0.24


def make_K(n, rank=15, seed=0):
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


def sim_random_block(n, n_snps=4, n_haps=6, seed=0):
    rng = np.random.default_rng(seed)
    patterns = []
    while len(patterns) < n_haps:
        pat = tuple(int(x) for x in rng.integers(0, 2, size=n_snps))
        if pat not in patterns:
            patterns.append(pat)
    w = rng.dirichlet(np.ones(n_haps) * 1.5)
    idx1 = rng.choice(n_haps, size=n, p=w)
    idx2 = rng.choice(n_haps, size=n, p=w)
    hap1 = np.asarray([patterns[i] for i in idx1], dtype=np.int8)
    hap2 = np.asarray([patterns[i] for i in idx2], dtype=np.int8)
    return hap1, hap2, np.asarray(patterns, dtype=np.int8)


def sim_neighbor_block(n, seed):
    rng = np.random.default_rng(seed)
    patterns = [
        (0, 0, 0, 0),
        (0, 0, 0, 1),  # Hamming-1 from first
        (0, 0, 1, 0),
        (0, 1, 0, 0),
        (1, 0, 0, 0),
        (1, 1, 0, 0),
    ]
    w = np.array([0.22, 0.18, 0.15, 0.15, 0.15, 0.15], float)
    w /= w.sum()
    idx1 = rng.choice(len(patterns), size=n, p=w)
    idx2 = rng.choice(len(patterns), size=n, p=w)
    hap1 = np.asarray([patterns[i] for i in idx1], dtype=np.int8)
    hap2 = np.asarray([patterns[i] for i in idx2], dtype=np.int8)
    return hap1, hap2, patterns


def block_stats(null, hap1, hap2):
    pairs = haplotypes_from_phased_snps(hap1, hap2)
    hap = build_dosage(pairs, min_freq=0.01, min_count=5, collapse_mode="freq", drop_ref=True)
    ms = multi_allelic_scores(null, hap.dosage, allele_labels=hap.alleles, hamming_max=1)
    Z = project_features(null, hap.dosage)
    y = null.y_tilde
    allele_ps = []
    for j in range(Z.shape[1]):
        ss = float(Z[:, j] @ Z[:, j])
        if ss <= 1e-8:
            continue
        allele_ps.append(float(stats.chi2.sf((Z[:, j] @ y) ** 2 / ss, 1)))
    allele_min = float(np.min(allele_ps)) if allele_ps else np.nan
    allele_bon = min(1.0, allele_min * len(allele_ps)) if allele_ps else np.nan
    snp_ps = []
    for j in range(hap1.shape[1]):
        dos = hap1[:, j].astype(float) + hap2[:, j].astype(float)
        _c, p = score_snp(null, dos)
        if np.isfinite(p):
            snp_ps.append(p)
    snp_min = float(np.min(snp_ps)) if snp_ps else np.nan
    snp_bon = min(1.0, snp_min * len(snp_ps)) if snp_ps else np.nan
    return {
        "snp_bon_p": snp_bon,
        "allele_bon_p": allele_bon,
        "joint_p": ms.full_p,
        "joint_df": ms.full_df,
        "acj_p": ms.acj.p,
        "acj_df": ms.acj.df,
        "acj_merges": ms.acj.merge_events,
        "acj_ngroups": ms.acj.n_groups,
        "pca_p": ms.pca.p,
        "ridge_p": ms.ridge.p,
        "n_alleles": int(hap.dosage.shape[1]),
    }


def run_null(n_rep, n, seed):
    rows = []
    for r in range(n_rep):
        hap1, hap2, _ = sim_random_block(n, seed=seed + r)
        K = make_K(n, seed=seed + 10000 + r)
        y = phenotype(n, None, K, 0.3, seed + 20000 + r)
        null = fit_null_lmm(y, X=np.ones((n, 1)), K=K, method="auto")
        d = block_stats(null, hap1, hap2)
        d.update({"rep": r, "scenario": "null"})
        rows.append(d)
        if (r + 1) % 50 == 0:
            print(f"  null {r+1}/{n_rep}", flush=True)
    return pd.DataFrame(rows)


def run_single_allele(n_rep, n, beta, seed):
    rows = []
    for r in range(n_rep):
        hap1, hap2, patterns = sim_random_block(n, seed=seed + r)
        pairs = haplotypes_from_phased_snps(hap1, hap2)
        causal = "".join(str(int(x)) for x in patterns[0])
        dosage = np.array([(pairs[i, 0] == causal) + (pairs[i, 1] == causal) for i in range(n)], float)
        if dosage.std() < 1e-6 or dosage.mean() < 0.05 or dosage.mean() > 1.5:
            causal = "".join(str(int(x)) for x in patterns[1])
            dosage = np.array([(pairs[i, 0] == causal) + (pairs[i, 1] == causal) for i in range(n)], float)
        K = make_K(n, seed=seed + 10000 + r)
        effect = beta * (dosage - dosage.mean())
        y = phenotype(n, effect, K, 0.3, seed + 20000 + r)
        null = fit_null_lmm(y, X=np.ones((n, 1)), K=K, method="auto")
        d = block_stats(null, hap1, hap2)
        d.update({"rep": r, "scenario": "single_allele", "beta": beta})
        rows.append(d)
        if (r + 1) % 50 == 0:
            print(f"  single {r+1}/{n_rep}", flush=True)
    return pd.DataFrame(rows)


def run_neighbor(n_rep, n, beta, seed, causal=True):
    rows = []
    for r in range(n_rep):
        hap1, hap2, patterns = sim_neighbor_block(n, seed + r)
        pairs = haplotypes_from_phased_snps(hap1, hap2)
        c0 = "".join(map(str, patterns[0]))
        c1 = "".join(map(str, patterns[1]))
        d0 = np.array([(pairs[i, 0] == c0) + (pairs[i, 1] == c0) for i in range(n)], float)
        d1 = np.array([(pairs[i, 0] == c1) + (pairs[i, 1] == c1) for i in range(n)], float)
        dosage = d0 + d1
        K = make_K(n, seed=seed + 10000 + r)
        effect = beta * (dosage - dosage.mean()) if causal else None
        y = phenotype(n, effect, K, 0.3, seed + 20000 + r)
        null = fit_null_lmm(y, X=np.ones((n, 1)), K=K, method="auto")
        d = block_stats(null, hap1, hap2)
        d.update({"rep": r, "scenario": "neighbor" if causal else "neighbor_null", "beta": beta if causal else 0.0})
        rows.append(d)
        if (r + 1) % 50 == 0:
            print(f"  neighbor{'null' if not causal else ''} {r+1}/{n_rep}", flush=True)
    return pd.DataFrame(rows)


METHODS = ["snp_bon_p", "allele_bon_p", "joint_p", "acj_p", "pca_p", "ridge_p"]


def summarize(df, alpha=0.05, kind="power"):
    rows = []
    for m in METHODS:
        p = df[m].to_numpy(float)
        p = p[np.isfinite(p)]
        rate = float(np.mean(p < alpha)) if len(p) else np.nan
        chi = np.full_like(p, np.nan)
        ok = p > 0
        chi[ok] = stats.chi2.ppf(1 - p[ok], 1)
        med = float(np.nanmedian(chi)) if len(chi) else np.nan
        lam = med / 0.4549364 if np.isfinite(med) else np.nan
        rows.append({"method": m, "n": len(p), kind: rate, "lambda_gc": lam})
    return pd.DataFrame(rows)


def plot_tradeoff(sum_single, sum_neigh, sum_null, out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    focus = ["allele_bon_p", "joint_p", "acj_p", "pca_p", "ridge_p", "snp_bon_p"]
    labels = ["allele-Bonf", "joint", "ACJ", "PCA", "ridge", "SNP-Bonf"]
    fig, axes = plt.subplots(1, 3, figsize=(12.5, 4.2))

    # Type I
    ax = axes[0]
    sub = sum_null.set_index("method").loc[focus]
    ax.bar(labels, sub["typeI"].to_numpy(), color="#4C78A8")
    ax.axhline(0.05, color="red", ls="--", lw=1)
    ax.set_ylim(0, max(0.12, float(sub["typeI"].max()) * 1.2))
    ax.set_ylabel("Type I @ 0.05")
    ax.set_title(f"Null (n={N}, reps={N_REP})")
    ax.tick_params(axis="x", rotation=30)

    # Single allele power
    ax = axes[1]
    sub = sum_single.set_index("method").loc[focus]
    colors = ["#F58518" if m != "acj_p" else "#E45756" for m in focus]
    ax.bar(labels, sub["power"].to_numpy(), color=colors)
    ax.set_ylim(0, 1)
    ax.set_ylabel("Power @ 0.05")
    ax.set_title(f"Single HapAllele β={BETA_SINGLE}")
    ax.tick_params(axis="x", rotation=30)

    # Neighbor power
    ax = axes[2]
    sub = sum_neigh.set_index("method").loc[focus]
    colors = ["#54A24B" if m == "acj_p" else "#F58518" for m in focus]
    ax.bar(labels, sub["power"].to_numpy(), color=colors)
    ax.set_ylim(0, 1)
    ax.set_ylabel("Power @ 0.05")
    ax.set_title(f"Neighbor alleles β={BETA_NEIGHBOR}")
    ax.tick_params(axis="x", rotation=30)

    fig.suptitle("ACJ power tradeoff: sharp single allele vs Hamming-neighbor multi-allele", y=1.02, fontsize=11)
    fig.tight_layout()
    path = out / "fig_acj_power_tradeoff.png"
    fig.savefig(path, dpi=160, bbox_inches="tight")
    plt.close(fig)
    return path


def main():
    print(f"SCALED n={N} n_rep={N_REP}", flush=True)
    null_df = run_null(N_REP, N, SEED)
    null_df.to_csv(OUT / "scaled_null_replicates.tsv", sep="\t", index=False)
    sum_null = summarize(null_df, kind="typeI")
    sum_null.to_csv(OUT / "scaled_typeI_summary.tsv", sep="\t", index=False)
    print(sum_null.to_string(index=False), flush=True)

    single_df = run_single_allele(N_REP, N, BETA_SINGLE, SEED + 500)
    single_df.to_csv(OUT / "scaled_power_single_replicates.tsv", sep="\t", index=False)
    sum_single = summarize(single_df, kind="power")
    sum_single.to_csv(OUT / "scaled_power_single_summary.tsv", sep="\t", index=False)
    print(sum_single.to_string(index=False), flush=True)

    neigh_df = run_neighbor(N_REP, N, BETA_NEIGHBOR, SEED + 900, causal=True)
    neigh_df.to_csv(OUT / "scaled_power_neighbor_replicates.tsv", sep="\t", index=False)
    sum_neigh = summarize(neigh_df, kind="power")
    sum_neigh.to_csv(OUT / "scaled_power_neighbor_summary.tsv", sep="\t", index=False)
    print(sum_neigh.to_string(index=False), flush=True)

    neigh_null = run_neighbor(N_REP, N, BETA_NEIGHBOR, SEED + 1300, causal=False)
    neigh_null.to_csv(OUT / "scaled_null_neighbor_replicates.tsv", sep="\t", index=False)
    sum_neigh_null = summarize(neigh_null, kind="typeI")
    sum_neigh_null.to_csv(OUT / "scaled_typeI_neighbor_summary.tsv", sep="\t", index=False)
    print("neighbor-design Type I:\n", sum_neigh_null.to_string(index=False), flush=True)

    fig = plot_tradeoff(sum_single, sum_neigh, sum_null, OUT)
    summary = {
        "n": N,
        "n_rep": N_REP,
        "beta_single": BETA_SINGLE,
        "beta_neighbor": BETA_NEIGHBOR,
        "typeI": sum_null.to_dict(orient="records"),
        "power_single": sum_single.to_dict(orient="records"),
        "power_neighbor": sum_neigh.to_dict(orient="records"),
        "typeI_neighbor_design": sum_neigh_null.to_dict(orient="records"),
        "figure": str(fig),
    }
    (OUT / "scaled_summary.json").write_text(json.dumps(summary, indent=2))
    print("Wrote", fig, flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
