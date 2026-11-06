"""
Pre-specified check of the Brownian haplotype-tree score.

The target regime was fixed before looking at power:
four common haplotypes (frequency 0.14, so the rare-haplotype collapse
does not touch them) whose pairwise Hamming distance is 2, carrying the
same effect. A tree kernel should pool that clade. A neighbor collapse
should not. A 4-df joint test pays for the extra alleles.

Betas below were chosen in a 40-replicate pilot only so that power is not
stuck at 0 or 1. The kernel itself was not changed after that pilot.

Usage:
  python scripts/bench_tree.py --reps 1000
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy import stats

from hapgwas.blocks import Block
from hapgwas.grm import vanraden_grm
from hapgwas.joint_test import test_block
from hapgwas.lmm import fit_null_lmm
from hapgwas.tree_kernel import kuonen_sf

HAPS = ["000000", "000111", "111000", "111011", "111101", "111110"]
FREQ = np.array([0.30, 0.14, 0.14, 0.14, 0.14, 0.14])
CLADE = {2, 3, 4, 5}
ALPHA = 0.05
# Locked after an effect-size pilot. Not tuned to maximize the tree test.
BETA_CLADE = 0.10
BETA_SINGLE = 0.32
METHODS = ("tree_p", "joint_p", "acj_p", "allele_bonf_p", "snp_bonf_p", "omni_p")


def _check_saddlepoint() -> dict:
    """Weighted-chi-square tail should match a plain chi-square when every weight is 1."""
    df = 6
    lam = np.ones(df)
    qs = stats.chi2.ppf([0.5, 0.9, 0.95, 0.99, 0.999], df)
    got = [kuonen_sf(float(q), lam) for q in qs]
    target = [0.5, 0.1, 0.05, 0.01, 0.001]
    return {"target": target, "kuonen": got}


def _sample(rng: np.random.Generator, n: int):
    i1 = rng.choice(len(HAPS), size=n, p=FREQ)
    i2 = rng.choice(len(HAPS), size=n, p=FREQ)

    def pack(idx):
        return np.array([[int(c) for c in HAPS[i]] for i in idx], dtype=float)

    return pack(i1), pack(i2), i1, i2


def _dos(code_index: int, i1, i2) -> np.ndarray:
    return (i1 == code_index).astype(float) + (i2 == code_index).astype(float)


def _pvalues(h1, h2, y, K=None) -> dict:
    null = fit_null_lmm(y, K=K, method="lmm" if K is not None else "ols")
    res = test_block(null, Block("b", "1", 0, h1.shape[1]), h1, h2, min_freq=0.01, min_count=5)
    return {name: float(getattr(res, name)) for name in METHODS}


def _rate(flags) -> float:
    return float(np.mean(flags)) if flags else float("nan")


def _run_iid(reps: int, seed: int, n: int = 500) -> dict:
    rng = np.random.default_rng(seed)

    def y_null(i1, i2, rr):
        return rr.normal(size=n)

    def y_single(i1, i2, rr):
        return BETA_SINGLE * _dos(2, i1, i2) + rr.normal(size=n)

    def y_clade(i1, i2, rr):
        load = sum(_dos(j, i1, i2) for j in CLADE)
        return BETA_CLADE * load + rr.normal(size=n)

    def y_snp(i1, i2, h1, h2, rr):
        # Last SNP. It does not mark the clade (the first three SNPs do).
        return 0.22 * (h1[:, -1] + h2[:, -1]) + rr.normal(size=n)

    scenarios = {
        "null": lambda i1, i2, h1, h2, rr: y_null(i1, i2, rr),
        "single_common_allele": lambda i1, i2, h1, h2, rr: y_single(i1, i2, rr),
        "clade_distance2": lambda i1, i2, h1, h2, rr: y_clade(i1, i2, rr),
        "single_snp": lambda i1, i2, h1, h2, rr: y_snp(i1, i2, h1, h2, rr),
    }
    out = {"n": n, "reps": reps, "alpha": ALPHA, "beta_clade": BETA_CLADE, "beta_single": BETA_SINGLE}
    for name, effect in scenarios.items():
        flags = {m: [] for m in METHODS}
        for _ in range(reps):
            h1, h2, i1, i2 = _sample(rng, n)
            y = effect(i1, i2, h1, h2, rng)
            got = _pvalues(h1, h2, y)
            for m in METHODS:
                flags[m].append(np.isfinite(got[m]) and got[m] < ALPHA)
        out[name] = {m: _rate(flags[m]) for m in METHODS}
    return out


def _run_lmm_null(reps: int, seed: int, n: int = 400) -> dict:
    rng = np.random.default_rng(seed)
    dosage = rng.binomial(2, 0.3, size=(n, 180)).astype(float)
    K = vanraden_grm(dosage)
    evals, evecs = np.linalg.eigh(0.5 * (K + K.T))
    evals = np.clip(evals, 0.0, None)
    flags = {m: [] for m in METHODS}
    for _ in range(reps):
        u = evecs @ (np.sqrt(evals) * rng.normal(size=n))
        y = np.sqrt(0.4) * u + np.sqrt(0.6) * rng.normal(size=n)
        h1, h2, _, _ = _sample(rng, n)
        got = _pvalues(h1, h2, y, K=K)
        for m in METHODS:
            flags[m].append(np.isfinite(got[m]) and got[m] < ALPHA)
    return {"n": n, "reps": reps, "alpha": ALPHA, "rates": {m: _rate(flags[m]) for m in METHODS}}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--reps", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=19)
    args = parser.parse_args()
    summary = {
        "saddlepoint_chi2": _check_saddlepoint(),
        "iid": _run_iid(args.reps, args.seed),
        "lmm_null": _run_lmm_null(max(200, args.reps // 2), args.seed + 1),
    }
    text = json.dumps(summary, indent=2)
    print(text)
    out = Path(__file__).resolve().parents[1] / "sim_results" / "tree_bench.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(text, encoding="utf-8")
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
