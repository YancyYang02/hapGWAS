#!/usr/bin/env python3
"""Chip-like LD haplotype simulation parameterized by D8_qc window r/MAF."""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import norm

from hapgwas.adaptive_score import multi_allelic_scores
from hapgwas.dosage import build_dosage, haplotypes_from_phased_snps
from hapgwas.joint_test import score_snp
from hapgwas.lmm import fit_null_lmm, project_features

OUT = Path("/workspace/hapgwas/sim_results")
LD_JSON = Path("/workspace/hapgwas/chip_ld/d8_ld_windows.json")
N_REP = 400
N = 800
SEED = 2026
BETA_SINGLE = 0.35
BETA_NEIGHBOR = 0.22


def make_K(n, rank=20, seed=0):
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
    base = 0.0 if effect is None else effect
    return base + u + e


def nearest_psd(R):
    R = 0.5 * (R + R.T)
    np.fill_diagonal(R, 1.0)
    w, V = np.linalg.eigh(R)
    w = np.clip(w, 1e-6, None)
    R2 = V @ np.diag(w) @ V.T
    d = np.sqrt(np.diag(R2))
    R2 = R2 / np.outer(d, d)
    np.fill_diagonal(R2, 1.0)
    return R2


def simulate_haps_from_ld(n, maf, R, seed):
    """Gaussian copula -> binary alleles per chromosome; two chroms -> diploid genotypes with chip-like LD."""
    rng = np.random.default_rng(seed)
    maf = np.asarray(maf, float)
    m = len(maf)
    R = nearest_psd(np.asarray(R, float))
    # thresholds for allele coding 1 with freq maf
    thr = norm.ppf(1 - maf)
    L = np.linalg.cholesky(R)
    def one_chrom():
        z = rng.normal(size=(n, m)) @ L.T
        return (z > thr).astype(np.int8)
    hap1 = one_chrom()
    hap2 = one_chrom()
    return hap1, hap2


def neighborize_patterns(hap1, hap2, seed):
    """Force a Hamming-neighbor allele pair into the population by flipping last SNP on a fraction of 000.. haplotypes."""
    # Already continuous alleles from LD; for neighbor causal we pick two alleles differing by 1 SNP among observed.
    return hap1, hap2


def block_stats(null, hap1, hap2):
    pairs = haplotypes_from_phased_snps(hap1, hap2)
    hap = build_dosage(pairs, min_freq=0.01, min_count=5, collapse_mode="freq", drop_ref=True)
    if hap.dosage.shape[1] < 1:
        return None
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
        "acj_p": ms.acj.p,
        "acj_df": ms.acj.df,
        "acj_merges": ms.acj.merge_events,
        "pca_p": ms.pca.p,
        "ridge_p": ms.ridge.p,
        "n_alleles": int(hap.dosage.shape[1]),
        "alleles": list(hap.alleles),
    }


def pick_single_causal(pairs, alleles_fallback=None):
    # most common non-ref-looking allele string among observed pairs
    from collections import Counter
    c = Counter()
    for a, b in pairs:
        c[a] += 1
        c[b] += 1
    # pick mid-frequency
    items = sorted(c.items(), key=lambda kv: kv[1], reverse=True)
    for lab, cnt in items[1:4] if len(items) > 1 else items:
        return lab
    return items[0][0]


def pick_neighbor_pair(pairs):
    from collections import Counter
    c = Counter()
    for a, b in pairs:
        c[a] += 1
        c[b] += 1
    labs = [lab for lab, _ in c.most_common(12) if lab != "OTHER"]
    for i, a in enumerate(labs):
        for b in labs[i + 1 :]:
            if len(a) == len(b) and sum(x != y for x, y in zip(a, b)) == 1:
                return a, b
    # fabricate by flipping last base of top allele
    if labs:
        a = labs[0]
        b = a[:-1] + ("1" if a[-1] == "0" else "0")
        return a, b
    return None, None


def run_scenario(windows, scenario, n_rep, n, beta, seed):
    rows = []
    rng = np.random.default_rng(seed)
    for r in range(n_rep):
        win = windows[int(rng.integers(0, len(windows)))]
        maf = np.clip(np.asarray(win["maf"], float), 0.05, 0.45)
        R = np.asarray(win["r"], float)
        hap1, hap2 = simulate_haps_from_ld(n, maf, R, seed + r)
        pairs = haplotypes_from_phased_snps(hap1, hap2)
        if scenario == "null":
            effect = None
        elif scenario == "single":
            causal = pick_single_causal(pairs)
            dosage = np.array([(pairs[i, 0] == causal) + (pairs[i, 1] == causal) for i in range(n)], float)
            if dosage.std() < 1e-6:
                continue
            effect = beta * (dosage - dosage.mean())
        elif scenario == "neighbor":
            a, b = pick_neighbor_pair(pairs)
            if a is None:
                continue
            d0 = np.array([(pairs[i, 0] == a) + (pairs[i, 1] == a) for i in range(n)], float)
            d1 = np.array([(pairs[i, 0] == b) + (pairs[i, 1] == b) for i in range(n)], float)
            dosage = d0 + d1
            if dosage.std() < 1e-6:
                continue
            effect = beta * (dosage - dosage.mean())
        else:
            raise ValueError(scenario)
        K = make_K(n, seed=seed + 10000 + r)
        y = phenotype(n, effect, K, 0.3, seed + 20000 + r)
        null = fit_null_lmm(y, X=np.ones((n, 1)), K=K, method="auto")
        d = block_stats(null, hap1, hap2)
        if d is None:
            continue
        d.update({"rep": r, "scenario": scenario, "chrom": win["chrom"], "bp0": win["bp"][0], "mean_abs_r_adj": win["mean_abs_r_adj"]})
        rows.append(d)
        if (r + 1) % 50 == 0:
            print(f"  {scenario} {r+1}/{n_rep}", flush=True)
    return pd.DataFrame(rows)


METHODS = ["snp_bon_p", "allele_bon_p", "joint_p", "acj_p", "pca_p", "ridge_p"]


def summarize(df, kind="power", alpha=0.05):
    rows = []
    for m in METHODS:
        p = df[m].to_numpy(float)
        p = p[np.isfinite(p)]
        rows.append({"method": m, "n": len(p), kind: float(np.mean(p < alpha)) if len(p) else np.nan})
    return pd.DataFrame(rows)


def main():
    if not LD_JSON.exists():
        print("MISSING", LD_JSON)
        return 1
    payload = json.loads(LD_JSON.read_text())
    windows = [w for w in payload["windows"] if w.get("mean_abs_r_adj", 0) >= 0.05]
    if len(windows) < 3:
        windows = payload["windows"]
    print(f"Using {len(windows)}/{payload['n_windows']} D8 windows; n={N} n_rep={N_REP}", flush=True)

    null = run_scenario(windows, "null", N_REP, N, 0.0, SEED)
    null.to_csv(OUT / "chipld_null_replicates.tsv", sep="\t", index=False)
    snull = summarize(null, "typeI")
    snull.to_csv(OUT / "chipld_typeI_summary.tsv", sep="\t", index=False)
    print(snull.to_string(index=False), flush=True)

    single = run_scenario(windows, "single", N_REP, N, BETA_SINGLE, SEED + 50)
    single.to_csv(OUT / "chipld_power_single_replicates.tsv", sep="\t", index=False)
    ss = summarize(single, "power")
    ss.to_csv(OUT / "chipld_power_single_summary.tsv", sep="\t", index=False)
    print(ss.to_string(index=False), flush=True)

    neigh = run_scenario(windows, "neighbor", N_REP, N, BETA_NEIGHBOR, SEED + 90)
    neigh.to_csv(OUT / "chipld_power_neighbor_replicates.tsv", sep="\t", index=False)
    sn = summarize(neigh, "power")
    sn.to_csv(OUT / "chipld_power_neighbor_summary.tsv", sep="\t", index=False)
    print(sn.to_string(index=False), flush=True)

    # figure panel
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    focus = ["allele_bon_p", "joint_p", "acj_p", "pca_p", "ridge_p", "snp_bon_p"]
    labels = ["allele-Bonf", "joint", "ACJ", "PCA", "ridge", "SNP-Bonf"]
    fig, axes = plt.subplots(1, 3, figsize=(12.5, 4.2))
    axes[0].bar(labels, snull.set_index("method").loc[focus, "typeI"], color="#4C78A8")
    axes[0].axhline(0.05, color="red", ls="--", lw=1)
    axes[0].set_title(f"Chip-LD null (n={N}, {N_REP} reps)")
    axes[0].tick_params(axis="x", rotation=30)
    axes[0].set_ylabel("Type I")
    axes[1].bar(labels, ss.set_index("method").loc[focus, "power"],
                color=["#E45756" if m == "acj_p" else "#F58518" for m in focus])
    axes[1].set_ylim(0, 1)
    axes[1].set_title(f"Single allele β={BETA_SINGLE}")
    axes[1].tick_params(axis="x", rotation=30)
    axes[2].bar(labels, sn.set_index("method").loc[focus, "power"],
                color=["#54A24B" if m == "acj_p" else "#F58518" for m in focus])
    axes[2].set_ylim(0, 1)
    axes[2].set_title(f"Neighbor alleles β={BETA_NEIGHBOR}")
    axes[2].tick_params(axis="x", rotation=30)
    fig.suptitle("D8 chip-LD parameterized sims (phased VCF unavailable locally)", y=1.02)
    fig.tight_layout()
    figpath = OUT / "fig_chipld_acj_tradeoff.png"
    fig.savefig(figpath, dpi=160, bbox_inches="tight")
    plt.close(fig)

    summary = {
        "source": payload.get("source"),
        "n_windows_used": len(windows),
        "n": N,
        "n_rep": N_REP,
        "beta_single": BETA_SINGLE,
        "beta_neighbor": BETA_NEIGHBOR,
        "typeI": snull.to_dict(orient="records"),
        "power_single": ss.to_dict(orient="records"),
        "power_neighbor": sn.to_dict(orient="records"),
        "figure": str(figpath),
        "data_gap": (
            "Breast-shape phased VCF/PED not on local machine; only on HPC "
            "/fs2/home/clwen/YancongYang/DoubleMuscling/variant_calling_results/"
            "step6.beagle_phased.vcf.gz and step7_phased_final.{bed,bim,fam}. "
            "Local unphased chip D8_qc used for LD parameterization."
        ),
    }
    (OUT / "chipld_summary.json").write_text(json.dumps(summary, indent=2))
    print("Wrote", figpath, flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
