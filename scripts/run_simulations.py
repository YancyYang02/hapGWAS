#!/usr/bin/env python3
"""hapGWAS simulation suite: Type I + power with calibrated methods."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from hapgwas.adaptive_score import multi_allelic_scores
from hapgwas.dosage import build_dosage, haplotypes_from_phased_snps
from hapgwas.joint_test import score_snp
from hapgwas.lmm import fit_null_lmm, project_features


def simulate_block_phased(n, n_snps=4, n_haps=6, seed=0):
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


def make_K(n, rank=15, seed=0):
    rng = np.random.default_rng(seed)
    Z = rng.normal(size=(n, rank))
    K = Z @ Z.T / rank
    return K / np.mean(np.diag(K))


def phenotype(n, effect, K, h2, seed):
    rng = np.random.default_rng(seed)
    if K is None:
        e = rng.normal(0, 1.0, size=n)
        return (effect if effect is not None else 0.0) + e
    evals, evecs = np.linalg.eigh(K)
    evals = np.clip(evals, 1e-8, None)
    z = rng.normal(size=n)
    u = evecs @ (np.sqrt(evals) * z)
    u = u / np.std(u) * np.sqrt(h2)
    e = rng.normal(0, np.sqrt(max(1 - h2, 1e-6)), size=n)
    base = effect if effect is not None else 0.0
    return base + u + e


def block_methods(null, hap1, hap2, min_freq=0.01, min_count=5, with_acj_pheno=False):
    pairs = haplotypes_from_phased_snps(hap1, hap2)
    hap = build_dosage(pairs, min_freq=min_freq, min_count=min_count, collapse_mode="freq", drop_ref=True)

    snp_ps = []
    for j in range(hap1.shape[1]):
        dos = hap1[:, j].astype(float) + hap2[:, j].astype(float)
        _c, p = score_snp(null, dos)
        if np.isfinite(p):
            snp_ps.append(p)
    snp_min = float(np.min(snp_ps)) if snp_ps else np.nan
    snp_bon = min(1.0, snp_min * len(snp_ps)) if snp_ps else np.nan

    ms = multi_allelic_scores(null, hap.dosage, allele_labels=hap.alleles, with_acj_pheno=with_acj_pheno, acj_pheno_perms=49)
    Z = project_features(null, hap.dosage)
    y = null.y_tilde
    allele_ps = []
    for j in range(Z.shape[1]):
        ss = float(Z[:, j] @ Z[:, j])
        if ss <= 1e-8:
            continue
        chisq = float((Z[:, j] @ y) ** 2 / ss)
        allele_ps.append(float(stats.chi2.sf(chisq, 1)))
    allele_min = float(np.min(allele_ps)) if allele_ps else np.nan
    allele_bon = min(1.0, allele_min * len(allele_ps)) if allele_ps else np.nan

    out = {
        "n_alleles": int(hap.dosage.shape[1]),
        "n_snps": int(hap1.shape[1]),
        "snp_min_p": snp_min,
        "snp_bon_p": snp_bon,
        "allele_min_p": allele_min,
        "allele_bon_p": allele_bon,
        "joint_p": ms.full_p,
        "joint_df": ms.full_df,
        "acj_p": ms.acj.p,
        "acj_df": ms.acj.df,
        "acj_ngroups": ms.acj.n_groups,
        "acj_merges": ms.acj.merge_events,
        "pca_p": ms.pca.p,
        "pca_df": ms.pca.df,
        "ridge_p": ms.ridge.p,
        "ridge_df": ms.ridge.df_satt,
    }
    if ms.acj_pheno is not None:
        out["acj_pheno_p"] = ms.acj_pheno.p
        out["acj_pheno_df"] = ms.acj_pheno.df
    return out


def run_null(n_rep, n, use_K, seed, with_acj_pheno=False):
    rows = []
    for r in range(n_rep):
        hap1, hap2, _ = simulate_block_phased(n, seed=seed + r)
        K = make_K(n, seed=seed + 10000 + r) if use_K else None
        y = phenotype(n, None, K, h2=0.3 if use_K else 0.0, seed=seed + 20000 + r)
        null = fit_null_lmm(y, X=np.ones((n, 1)), K=K, method="auto")
        stats_d = block_methods(null, hap1, hap2, with_acj_pheno=with_acj_pheno)
        stats_d.update({"rep": r, "scenario": "null", "use_K": use_K})
        rows.append(stats_d)
    return pd.DataFrame(rows)


def run_power_hap(n_rep, n, beta, use_K, seed, with_acj_pheno=False):
    rows = []
    for r in range(n_rep):
        hap1, hap2, patterns = simulate_block_phased(n, seed=seed + r)
        pairs = haplotypes_from_phased_snps(hap1, hap2)
        causal = "".join(str(int(x)) for x in patterns[0])
        dosage = np.array([(pairs[i, 0] == causal) + (pairs[i, 1] == causal) for i in range(n)], dtype=float)
        if dosage.std() < 1e-6 or dosage.mean() < 0.05 or dosage.mean() > 1.5:
            causal = "".join(str(int(x)) for x in patterns[1])
            dosage = np.array([(pairs[i, 0] == causal) + (pairs[i, 1] == causal) for i in range(n)], dtype=float)
        K = make_K(n, seed=seed + 10000 + r) if use_K else None
        effect = beta * (dosage - dosage.mean())
        y = phenotype(n, effect, K, h2=0.3 if use_K else 0.0, seed=seed + 20000 + r)
        null = fit_null_lmm(y, X=np.ones((n, 1)), K=K, method="auto")
        stats_d = block_methods(null, hap1, hap2, with_acj_pheno=with_acj_pheno)
        stats_d.update({"rep": r, "scenario": "causal_hap", "beta": beta, "use_K": use_K})
        rows.append(stats_d)
    return pd.DataFrame(rows)


def run_power_snp(n_rep, n, beta, use_K, seed, with_acj_pheno=False):
    rows = []
    for r in range(n_rep):
        hap1, hap2, _ = simulate_block_phased(n, seed=seed + r)
        dosage = hap1[:, 1].astype(float) + hap2[:, 1].astype(float)
        K = make_K(n, seed=seed + 10000 + r) if use_K else None
        effect = beta * (dosage - dosage.mean())
        y = phenotype(n, effect, K, h2=0.3 if use_K else 0.0, seed=seed + 20000 + r)
        null = fit_null_lmm(y, X=np.ones((n, 1)), K=K, method="auto")
        stats_d = block_methods(null, hap1, hap2, with_acj_pheno=with_acj_pheno)
        stats_d.update({"rep": r, "scenario": "causal_snp", "beta": beta, "use_K": use_K})
        rows.append(stats_d)
    return pd.DataFrame(rows)


METHODS = [
    "snp_min_p", "snp_bon_p",
    "allele_min_p", "allele_bon_p",
    "joint_p", "acj_p", "pca_p", "ridge_p",
]


def summarize_typeI(df, alpha=0.05):
    methods = [m for m in METHODS if m in df.columns]
    if "acj_pheno_p" in df.columns:
        methods.append("acj_pheno_p")
    rows = []
    for m in methods:
        p = df[m].to_numpy(dtype=float)
        p = p[np.isfinite(p)]
        rate = float(np.mean(p < alpha)) if len(p) else np.nan
        chi = np.full_like(p, np.nan)
        ok = p > 0
        chi[ok] = stats.chi2.ppf(1 - p[ok], 1)
        med = float(np.nanmedian(chi)) if len(chi) else np.nan
        lambda_gc = med / 0.4549364 if np.isfinite(med) else np.nan
        rows.append({"method": m, "n": len(p), "typeI_alpha0.05": rate, "lambda_gc": lambda_gc})
    return pd.DataFrame(rows)


def summarize_power(df, alpha=0.05):
    methods = [m for m in METHODS if m in df.columns]
    if "acj_pheno_p" in df.columns:
        methods.append("acj_pheno_p")
    rows = []
    for m in methods:
        p = df[m].to_numpy(dtype=float)
        p = p[np.isfinite(p)]
        rows.append({"method": m, "n": len(p), f"power_alpha{alpha}": float(np.mean(p < alpha)) if len(p) else np.nan})
    return pd.DataFrame(rows)


def maybe_plot(out, null_df, typeI_tab, pow_hap_tab, pow_snp_tab):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as e:
        print("skip figures:", e)
        return

    # Type I: calibrated methods focus
    focus = typeI_tab[typeI_tab["method"].isin(
        ["snp_bon_p", "allele_bon_p", "joint_p", "acj_p", "pca_p", "ridge_p"]
    )]
    fig, ax = plt.subplots(figsize=(7.5, 4))
    ax.bar(focus["method"], focus["typeI_alpha0.05"], color="#4C78A8")
    ax.axhline(0.05, color="red", ls="--", lw=1)
    ax.set_ylabel("Type I @ 0.05")
    ax.set_title("Null Type I (Bonferroni for min-p; ACJ corr-calibrated)")
    ax.tick_params(axis="x", rotation=25)
    fig.tight_layout()
    fig.savefig(out / "fig_typeI_rates.png", dpi=140)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4), sharey=True)
    for ax, tab, title in [
        (axes[0], pow_hap_tab, "Power: causal haplotype"),
        (axes[1], pow_snp_tab, "Power: causal SNP"),
    ]:
        col = [c for c in tab.columns if c.startswith("power_")][0]
        sub = tab[tab["method"].isin(
            ["snp_bon_p", "allele_bon_p", "joint_p", "acj_p", "pca_p", "ridge_p"]
        )]
        ax.bar(sub["method"], sub[col], color="#F58518")
        ax.set_title(title)
        ax.tick_params(axis="x", rotation=25)
        ax.set_ylim(0, 1)
        ax.set_ylabel("Power @ 0.05")
    fig.tight_layout()
    fig.savefig(out / "fig_power_compare.png", dpi=140)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(5, 5))
    for col, lab in [("joint_p", "joint"), ("acj_p", "ACJ"), ("allele_bon_p", "allele-Bonf")]:
        p = np.sort(null_df[col].dropna().to_numpy())
        if len(p) == 0:
            continue
        exp = -np.log10((np.arange(1, len(p) + 1) - 0.5) / len(p))
        obs = -np.log10(np.clip(p, 1e-300, 1))
        ax.plot(exp, obs, ".", label=lab, alpha=0.7)
    ax.plot([0, 4], [0, 4], "k--", lw=1)
    ax.set_xlabel(r"Expected $-\log_{10}p$")
    ax.set_ylabel(r"Observed $-\log_{10}p$")
    ax.set_title("Null QQ (calibrated)")
    ax.legend()
    fig.tight_layout()
    fig.savefig(out / "fig_null_qq.png", dpi=140)
    plt.close(fig)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="/workspace/hapgwas/sim_results")
    ap.add_argument("--n-rep", type=int, default=200)
    ap.add_argument("--n", type=int, default=400)
    ap.add_argument("--beta-hap", type=float, default=0.45)
    ap.add_argument("--beta-snp", type=float, default=0.35)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--fast", action="store_true")
    ap.add_argument("--acj-pheno", action="store_true", help="also run ACJ+ with permutations (slower)")
    args = ap.parse_args(argv)
    if args.fast:
        args.n_rep = 80
        args.n = 300

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    print(f"null n_rep={args.n_rep} n={args.n} ...")
    null_df = run_null(args.n_rep, args.n, True, args.seed, with_acj_pheno=args.acj_pheno)
    null_df.to_csv(out / "null_replicates.tsv", sep="\t", index=False)
    typeI = summarize_typeI(null_df)
    typeI.to_csv(out / "typeI_summary.tsv", sep="\t", index=False)
    print(typeI.to_string(index=False))

    print(f"power hap beta={args.beta_hap} ...")
    pow_hap = run_power_hap(args.n_rep, args.n, args.beta_hap, True, args.seed + 500, with_acj_pheno=args.acj_pheno)
    pow_hap.to_csv(out / "power_hap_replicates.tsv", sep="\t", index=False)
    pow_hap_tab = summarize_power(pow_hap)
    pow_hap_tab.to_csv(out / "power_hap_summary.tsv", sep="\t", index=False)
    print(pow_hap_tab.to_string(index=False))

    print(f"power snp beta={args.beta_snp} ...")
    pow_snp = run_power_snp(args.n_rep, args.n, args.beta_snp, True, args.seed + 900, with_acj_pheno=args.acj_pheno)
    pow_snp.to_csv(out / "power_snp_replicates.tsv", sep="\t", index=False)
    pow_snp_tab = summarize_power(pow_snp)
    pow_snp_tab.to_csv(out / "power_snp_summary.tsv", sep="\t", index=False)
    print(pow_snp_tab.to_string(index=False))

    maybe_plot(out, null_df, typeI, pow_hap_tab, pow_snp_tab)
    summary = {
        "n_rep": args.n_rep,
        "n": args.n,
        "beta_hap": args.beta_hap,
        "beta_snp": args.beta_snp,
        "typeI": typeI.to_dict(orient="records"),
        "power_hap": pow_hap_tab.to_dict(orient="records"),
        "power_snp": pow_snp_tab.to_dict(orient="records"),
        "notes": "ACJ uses genotype-correlation collapse (chi2-calibrated). Raw min-p is anti-conservative; use *_bon_p for fair Type I.",
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2))
    print("Wrote", out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
