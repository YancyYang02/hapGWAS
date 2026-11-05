#!/usr/bin/env python3
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
from hapgwas.adaptive_score import multi_allelic_scores
from hapgwas.dosage import build_dosage, haplotypes_from_phased_snps
from hapgwas.lmm import fit_null_lmm, project_features
from hapgwas.joint_test import score_snp

out = Path("/workspace/hapgwas/sim_results")
n_rep, n, beta = 200, 400, 0.28
seed = 123


def make_K(n, rank=15, seed=0):
    rng = np.random.default_rng(seed)
    Z = rng.normal(size=(n, rank))
    K = Z @ Z.T / rank
    return K / np.mean(np.diag(K))


def sim_collinear_haps(n, seed):
    rng = np.random.default_rng(seed)
    patterns = [
        (0, 0, 0, 0),
        (0, 0, 0, 1),
        (0, 0, 1, 0),
        (0, 1, 0, 0),
        (1, 0, 0, 0),
        (1, 1, 0, 0),
    ]
    w = np.array([0.22, 0.18, 0.15, 0.15, 0.15, 0.15])
    w = w / w.sum()
    idx1 = rng.choice(len(patterns), size=n, p=w)
    idx2 = rng.choice(len(patterns), size=n, p=w)
    hap1 = np.asarray([patterns[i] for i in idx1], dtype=np.int8)
    hap2 = np.asarray([patterns[i] for i in idx2], dtype=np.int8)
    return hap1, hap2, patterns


def one_rep(r, causal=True):
    hap1, hap2, patterns = sim_collinear_haps(n, seed + r)
    pairs = haplotypes_from_phased_snps(hap1, hap2)
    c0 = "".join(map(str, patterns[0]))
    c1 = "".join(map(str, patterns[1]))
    d0 = np.array([(pairs[i, 0] == c0) + (pairs[i, 1] == c0) for i in range(n)], float)
    d1 = np.array([(pairs[i, 0] == c1) + (pairs[i, 1] == c1) for i in range(n)], float)
    dosage = d0 + d1
    K = make_K(n, seed=seed + 10000 + r)
    evals, evecs = np.linalg.eigh(K)
    evals = np.clip(evals, 1e-8, None)
    rng = np.random.default_rng(seed + 20000 + r)
    u = evecs @ (np.sqrt(evals) * rng.normal(size=n))
    u = u / np.std(u) * np.sqrt(0.3)
    e = rng.normal(0, np.sqrt(0.7), size=n)
    if causal:
        y = beta * (dosage - dosage.mean()) + u + e
    else:
        y = u + e
    null = fit_null_lmm(y, X=np.ones((n, 1)), K=K, method="auto")
    hap = build_dosage(pairs, min_freq=0.01, min_count=5, collapse_mode="freq", drop_ref=True)
    ms = multi_allelic_scores(null, hap.dosage, allele_labels=hap.alleles, hamming_max=1)
    Z = project_features(null, hap.dosage)
    yy = null.y_tilde
    allele_ps = []
    for j in range(Z.shape[1]):
        ss = float(Z[:, j] @ Z[:, j])
        if ss <= 1e-8:
            continue
        allele_ps.append(float(stats.chi2.sf((Z[:, j] @ yy) ** 2 / ss, 1)))
    allele_min = min(allele_ps) if allele_ps else np.nan
    allele_bon = min(1.0, allele_min * len(allele_ps)) if allele_ps else np.nan
    snp_ps = []
    for j in range(hap1.shape[1]):
        dos = hap1[:, j].astype(float) + hap2[:, j].astype(float)
        _c, p = score_snp(null, dos)
        if np.isfinite(p):
            snp_ps.append(p)
    snp_min = min(snp_ps) if snp_ps else np.nan
    snp_bon = min(1.0, snp_min * len(snp_ps)) if snp_ps else np.nan
    return {
        "rep": r,
        "snp_bon_p": snp_bon,
        "allele_bon_p": allele_bon,
        "joint_p": ms.full_p,
        "joint_df": ms.full_df,
        "acj_p": ms.acj.p,
        "acj_df": ms.acj.df,
        "acj_merges": ms.acj.merge_events,
        "pca_p": ms.pca.p,
        "ridge_p": ms.ridge.p,
        "n_alleles": int(hap.dosage.shape[1]),
    }


def summarize(df, alpha=0.05):
    methods = ["snp_bon_p", "allele_bon_p", "joint_p", "acj_p", "pca_p", "ridge_p"]
    rows = []
    for m in methods:
        p = df[m].to_numpy(float)
        rows.append({"method": m, "n": len(p), f"rate_alpha{alpha}": float(np.mean(p < alpha))})
    return pd.DataFrame(rows)


print("power multi-allelic ...")
pow_rows = [one_rep(r, True) for r in range(n_rep)]
pdf = pd.DataFrame(pow_rows)
pdf.to_csv(out / "power_multi_allelic_replicates.tsv", sep="\t", index=False)
ptab = summarize(pdf).rename(columns={"rate_alpha0.05": "power_alpha0.05"})
ptab.to_csv(out / "power_multi_allelic_summary.tsv", sep="\t", index=False)
print(ptab.to_string(index=False))
print("merge>0", float((pdf.acj_merges > 0).mean()), "mean merges", float(pdf.acj_merges.mean()),
      "mean df joint/acj", float(pdf.joint_df.mean()), float(pdf.acj_df.mean()))

print("null collinear ...")
nrows = [one_rep(r, False) for r in range(n_rep)]
ndf = pd.DataFrame(nrows)
ndf.to_csv(out / "null_collinear_replicates.tsv", sep="\t", index=False)
ntab = summarize(ndf).rename(columns={"rate_alpha0.05": "typeI_alpha0.05"})
ntab.to_csv(out / "typeI_collinear_summary.tsv", sep="\t", index=False)
print(ntab.to_string(index=False))
print("null merge>0", float((ndf.acj_merges > 0).mean()))

summary = json.loads((out / "summary.json").read_text())
summary["power_multi_allelic"] = ptab.to_dict(orient="records")
summary["typeI_collinear"] = ntab.to_dict(orient="records")
summary["acj_defaults"] = {"corr_thresh": 0.5, "max_df": 3}
(out / "summary.json").write_text(json.dumps(summary, indent=2))
print("done")
