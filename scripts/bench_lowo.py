"""
Region-excluded haplotype likelihood versus the GRM that still contains the region.

The phenotype random effect is drawn from background SNPs only. The fitted
full GRM also contains the haplotype window and a flank built to tag it.
The region-excluded likelihood removes that flank with one Woodbury downdate.
Betas are fixed before the run.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from hapgwas.dosage import build_dosage, haplotypes_from_phased_snps
from hapgwas.exact_lrt import block_exact_lrts
from hapgwas.grm import vanraden_factors
from hapgwas.joint_test import score_joint_and_marginal
from hapgwas.lmm import fit_null_lmm
from hapgwas.window_kinship import lowo_joint_lrt, prepare_window_basis

HAPS = ["0000", "1111", "0110", "1001", "0011", "1100"]
ALPHAS = (0.05, 0.01, 0.001)


def _haps(rng: np.random.Generator, n: int):
    i1 = rng.integers(0, len(HAPS), size=n)
    i2 = rng.integers(0, len(HAPS), size=n)

    def pack(idx):
        return np.array([[int(ch) for ch in HAPS[i]] for i in idx], dtype=float)

    return pack(i1), pack(i2), i1, i2


def _dosage_and_pheno(rng, n, m_bg, n_flank, beta, herit):
    bg = rng.binomial(2, 0.3, size=(n, m_bg)).astype(float)
    h1, h2, i1, i2 = _haps(rng, n)
    window = h1 + h2
    flank = np.empty((n, n_flank), dtype=float)
    for j in range(n_flank):
        src = window[:, j % window.shape[1]]
        flank[:, j] = np.clip(np.rint(0.85 * src + 0.15 * rng.binomial(2, 0.3, size=n)), 0, 2)
    dosage = np.column_stack([bg, window, flank])
    _Z, _keep, K_bg = vanraden_factors(bg, maf_min=0.01)
    evals, evecs = np.linalg.eigh(K_bg)
    evals = np.clip(evals, 0.0, None)
    u = evecs @ (np.sqrt(evals) * rng.normal(size=n))
    causal = (i1 == 2).astype(float) + (i2 == 2).astype(float)
    y = beta * causal + np.sqrt(herit) * u + np.sqrt(1.0 - herit) * rng.normal(size=n)
    return dosage, h1, h2, y, m_bg


def _eval(dosage, h1, h2, y, m_bg, n_flank):
    basis = prepare_window_basis(y, dosage, maf_min=0.01)
    # Z columns follow MAF-kept original indices. Exclude window + flank.
    exclude_from = m_bg
    cols = np.flatnonzero(basis.keep >= exclude_from)
    hap = build_dosage(haplotypes_from_phased_snps(h1, h2), min_freq=0.01, min_count=5, drop_ref=True)
    lowo_stat, lowo_df, lowo_p = lowo_joint_lrt(basis, hap.dosage, cols)
    _Z, _keep, K_full = vanraden_factors(dosage, maf_min=0.01)
    null = fit_null_lmm(y, K=K_full, method="lmm")
    _s, _df, full_p, _bs, _bdf, _bp = block_exact_lrts(null, hap.dosage, h1 + h2)
    _c, _d, score_p, _mc, _mp = score_joint_and_marginal(null, hap.dosage)
    return {
        "lowo_p": float(lowo_p),
        "full_lrt_p": float(full_p),
        "full_score_p": float(score_p),
        "lowo_df": int(lowo_df),
        "excluded": int(cols.size),
    }


def _summarize(rows: list[dict]) -> dict:
    out = {"n": len(rows), "excluded_mean": float(np.mean([r["excluded"] for r in rows]))}
    for key in ("lowo_p", "full_lrt_p", "full_score_p"):
        p = np.array([r[key] for r in rows], dtype=float)
        ok = np.isfinite(p)
        entry = {"n_finite": int(ok.sum())}
        if ok.any():
            entry["median_neglog10"] = float(np.median(-np.log10(np.clip(p[ok], 1e-300, 1.0))))
        for alpha in ALPHAS:
            entry[f"a{alpha}"] = float(np.mean(p[ok] < alpha)) if ok.any() else None
        out[key] = entry
    return out


def main():
    rng = np.random.default_rng(19)
    n = 360
    m_bg = 200
    n_flank = 36
    h2 = 0.35
    summary = {
        "n": n,
        "background_snps": m_bg,
        "flank_snps": n_flank,
        "h2_background": h2,
        "note": "random effect from background SNPs only; full GRM also contains the window and its flank",
    }
    plan = (("null", 0.0, 400), ("hap_moderate", 0.45, 300), ("hap_large", 0.80, 250))
    for name, beta, reps in plan:
        rows = []
        for _ in range(reps):
            dosage, hap1, hap2, y, start = _dosage_and_pheno(rng, n, m_bg, n_flank, beta, h2)
            rows.append(_eval(dosage, hap1, hap2, y, start, n_flank))
        summary[name] = {"beta": beta, **_summarize(rows)}
        print(name, summary[name]["lowo_p"]["a0.05"], summary[name]["full_lrt_p"]["a0.05"], flush=True)
    text = json.dumps(summary, indent=2)
    out = Path(__file__).resolve().parents[1] / "sim_results" / "lowo_bench.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(text, encoding="utf-8")
    print(text)
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
