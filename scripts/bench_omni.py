"""
Calibration and power for the v0.5 default test (Cauchy omnibus).

Compares, at alpha=0.05, one block at a time:
  omni, joint, ACJ, allele-Bonferroni (GHap-like), SNP-Bonferroni

Also checks:
  - profile REML recovers a planted heritability
  - OLS vs LMM type I when family structure confounds the tested SNPs

Usage:
  python scripts/bench_omni.py --reps 200
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from hapgwas.blocks import Block
from hapgwas.grm import vanraden_grm
from hapgwas.joint_test import test_block
from hapgwas.lmm import fit_null_lmm


BASES = ["0000", "1110", "1100", "0011", "0101", "0111", "1000", "1111"]
FREQ = np.array([0.34, 0.16, 0.12, 0.10, 0.08, 0.08, 0.07, 0.05], dtype=float)
FREQ = FREQ / FREQ.sum()
ALPHA = 0.05
METHODS = ("omni_p", "joint_p", "acj_p", "allele_bonf_p", "snp_bonf_p")


def _haps(rng: np.random.Generator, n: int):
    i1 = rng.choice(len(BASES), size=n, p=FREQ)
    i2 = rng.choice(len(BASES), size=n, p=FREQ)

    def pack(idx):
        return np.array([[int(c) for c in BASES[i]] for i in idx], dtype=float)

    return pack(i1), pack(i2), i1, i2


def _dosage(code: str, i1, i2) -> np.ndarray:
    a = np.asarray(BASES)
    return (a[i1] == code).astype(float) + (a[i2] == code).astype(float)


def _pvalues(h1, h2, y, K=None) -> dict:
    null = fit_null_lmm(y, K=K, method="auto" if K is not None else "ols")
    res = test_block(
        null,
        Block("b", "1", 0, h1.shape[1]),
        h1,
        h2,
        min_freq=0.01,
        min_count=5,
    )
    return {name: float(getattr(res, name)) for name in METHODS}


def _rate(flags: list[bool]) -> float:
    return float(np.mean(flags)) if flags else float("nan")


def reml_check(seed: int = 11) -> dict:
    rng = np.random.default_rng(seed)
    n, m = 360, 180
    dosage = rng.binomial(2, 0.3, size=(n, m)).astype(float)
    K = vanraden_grm(dosage)
    evals, evecs = np.linalg.eigh(0.5 * (K + K.T))
    evals = np.clip(evals, 0, None)
    h2 = 0.40
    u = evecs @ (np.sqrt(evals) * rng.normal(size=n))
    u = u / np.std(u)
    y = np.sqrt(h2) * u + np.sqrt(1 - h2) * rng.normal(size=n)
    fit = fit_null_lmm(y, K=K, method="lmm")
    est = fit.vg / (fit.vg + fit.ve) if (fit.vg + fit.ve) > 0 else float("nan")
    return {"planted_h2": h2, "estimated_h2": est, "vg": fit.vg, "ve": fit.ve}


def correct_lmm_typeI(reps: int, seed: int) -> dict:
    """Type I when y is drawn from the same GRM the null model is given."""
    rng = np.random.default_rng(seed)
    n, m = 400, 160
    dosage = rng.binomial(2, 0.3, size=(n, m)).astype(float)
    K = vanraden_grm(dosage)
    evals, evecs = np.linalg.eigh(0.5 * (K + K.T))
    evals = np.clip(evals, 0.0, None)
    flags = {name: [] for name in METHODS}
    for _ in range(reps):
        u = evecs @ (np.sqrt(evals) * rng.normal(size=n))
        y = np.sqrt(0.4) * u + np.sqrt(0.6) * rng.normal(size=n)
        h1 = rng.integers(0, 2, size=(n, 4)).astype(float)
        h2 = rng.integers(0, 2, size=(n, 4)).astype(float)
        got = _pvalues(h1, h2, y, K=K)
        for name in METHODS:
            flags[name].append(np.isfinite(got[name]) and got[name] < ALPHA)
    return {"n": n, "reps": reps, "alpha": ALPHA, "lmm": {name: _rate(flags[name]) for name in METHODS}}


def structured_typeI(reps: int, seed: int) -> dict:
    rng = np.random.default_rng(seed)
    n_fam, n_per = 8, 30
    n = n_fam * n_per
    fam = np.repeat(np.arange(n_fam), n_per)
    ols = {m: [] for m in METHODS}
    lmm = {m: [] for m in METHODS}

    def panel(m):
        h1 = np.zeros((n, m))
        h2 = np.zeros((n, m))
        for j in range(m):
            p = np.clip(0.5 + 0.35 * np.sin(2 * np.pi * fam / n_fam + j), 0.05, 0.95)
            h1[:, j] = rng.random(n) < p
            h2[:, j] = rng.random(n) < p
        return h1, h2

    for _ in range(reps):
        bg1, bg2 = panel(80)
        t1, t2 = panel(4)
        K = vanraden_grm(bg1 + bg2)
        y = rng.normal(size=n_fam)[fam] + rng.normal(scale=0.55, size=n)
        a = _pvalues(t1, t2, y, K=None)
        b = _pvalues(t1, t2, y, K=K)
        for m in METHODS:
            ols[m].append(np.isfinite(a[m]) and a[m] < ALPHA)
            lmm[m].append(np.isfinite(b[m]) and b[m] < ALPHA)
    return {
        "n": n,
        "reps": reps,
        "alpha": ALPHA,
        "ols": {m: _rate(ols[m]) for m in METHODS},
        "lmm": {m: _rate(lmm[m]) for m in METHODS},
    }


def power_suite(reps: int, seed: int) -> dict:
    rng = np.random.default_rng(seed)
    n = 240
    scenarios = {
        "null": lambda i1, i2, h1, h2, rr: rr.normal(size=n),
        "single_allele": lambda i1, i2, h1, h2, rr: 0.55 * _dosage("0101", i1, i2) + rr.normal(size=n),
        "neighbor_alleles": lambda i1, i2, h1, h2, rr: 0.40 * (
            _dosage("0101", i1, i2) + _dosage("0111", i1, i2)
        ) + rr.normal(size=n),
        "single_snp": lambda i1, i2, h1, h2, rr: 0.35 * (h1[:, 0] + h2[:, 0]) + rr.normal(size=n),
    }
    out = {}
    for name, effect in scenarios.items():
        flags = {m: [] for m in METHODS}
        for _ in range(reps):
            h1, h2, i1, i2 = _haps(rng, n)
            y = effect(i1, i2, h1, h2, rng)
            got = _pvalues(h1, h2, y, K=None)
            for m in METHODS:
                flags[m].append(np.isfinite(got[m]) and got[m] < ALPHA)
        out[name] = {m: _rate(flags[m]) for m in METHODS}
    powers = [out[s] for s in ("single_allele", "neighbor_alleles", "single_snp")]
    out["mean_power_across_architectures"] = {
        m: float(np.mean([row[m] for row in powers])) for m in METHODS
    }
    out["n"] = n
    out["reps"] = reps
    out["alpha"] = ALPHA
    return out


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--reps", type=int, default=200)
    parser.add_argument("--seed", type=int, default=7)
    args = parser.parse_args()
    summary = {
        "reml": reml_check(args.seed),
        "typeI_lmm_correct": correct_lmm_typeI(args.reps, args.seed + 3),
        "typeI_structure_incomplete_grm": structured_typeI(args.reps, args.seed + 1),
        "iid_scores": power_suite(args.reps, args.seed + 2),
    }
    text = json.dumps(summary, indent=2)
    print(text)
    out = Path(__file__).resolve().parents[1] / "sim_results" / "omni_bench.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(text, encoding="utf-8")
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
