"""
Score test versus the exact nested haplotype likelihood ratio.

The kinship is a sibship matrix (groups of 8, relatedness 0.5), and the same
matrix is used to simulate and to fit. Betas are fixed before the run.
The exact tests re-estimate delta = ve/vg inside the block. The score tests
keep the null-model delta.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from hapgwas.blocks import Block
from hapgwas.joint_test import test_block
from hapgwas.lmm import fit_null_lmm

HAPS = ["0000", "1111", "0110", "1001", "0011", "1100"]
METHODS = ("joint_p", "joint_lrt_p", "beyond_p", "beyond_lrt_p", "snp_bonf_p", "omni_p")
ALPHAS = (0.05, 0.01, 0.001)


def sib_kinship(n: int, family: int = 8) -> np.ndarray:
    K = np.eye(n)
    for start in range(0, n, family):
        stop = min(n, start + family)
        K[start:stop, start:stop] = 0.5
        np.fill_diagonal(K[start:stop, start:stop], 1.0)
    return K


def _sample(rng: np.random.Generator, n: int):
    i1 = rng.integers(0, len(HAPS), size=n)
    i2 = rng.integers(0, len(HAPS), size=n)

    def pack(idx):
        return np.array([[int(ch) for ch in HAPS[i]] for i in idx], dtype=float)

    return pack(i1), pack(i2), i1, i2


def _fixed(kind: str, beta: float, i1, i2, h1, h2) -> np.ndarray:
    if kind == "null" or beta == 0.0:
        return np.zeros(i1.shape[0])
    if kind == "hap":
        return beta * ((i1 == 2).astype(float) + (i2 == 2).astype(float))
    if kind == "snp":
        return beta * (h1[:, 0] + h2[:, 0])
    if kind == "pair":
        return beta * (
            ((i1 == 2) | (i1 == 3)).astype(float) + ((i2 == 2) | (i2 == 3)).astype(float)
        )
    raise ValueError(kind)


def _summarize(rows: list[dict]) -> dict:
    out = {}
    for method in METHODS:
        p = np.array([row[method] for row in rows], dtype=float)
        ok = np.isfinite(p)
        entry = {"n_finite": int(ok.sum()), "n": int(p.size)}
        if ok.any():
            entry["median_neglog10"] = float(np.median(-np.log10(np.clip(p[ok], 1e-300, 1.0))))
        for alpha in ALPHAS:
            entry[f"a{alpha}"] = float(np.mean(p[ok] < alpha)) if ok.any() else None
        out[method] = entry
    finite = [
        (row["joint_lrt_p"], row["joint_p"])
        for row in rows
        if np.isfinite(row["joint_lrt_p"]) and np.isfinite(row["joint_p"])
    ]
    if finite:
        out["joint_lrt_smaller_than_score"] = float(np.mean([a < b for a, b in finite]))
    return out


def _collect(reps: int, rng, kind: str, beta: float, n: int, K, evecs, evals, vg: float, ve: float, ols: bool):
    rows = []
    for _ in range(reps):
        h1, h2, i1, i2 = _sample(rng, n)
        fixed = _fixed(kind, beta, i1, i2, h1, h2)
        if ols:
            y = fixed + rng.normal(size=n)
            null = fit_null_lmm(y, method="ols")
        else:
            u = evecs @ (np.sqrt(evals) * rng.normal(size=n))
            y = fixed + np.sqrt(vg) * u + np.sqrt(ve) * rng.normal(size=n)
            null = fit_null_lmm(y, K=K, method="lmm")
        res = test_block(null, Block("b", "1", 0, h1.shape[1]), h1, h2, min_count=5)
        rows.append({m: float(getattr(res, m)) for m in METHODS})
    return _summarize(rows)


def main():
    rng = np.random.default_rng(11)
    n = 320
    K = sib_kinship(n, 8)
    evals, evecs = np.linalg.eigh(K)
    evals = np.clip(evals, 0.0, None)
    summary = {
        "n": n,
        "family": 8,
        "relatedness_within_family": 0.5,
        "null_lmm_vg_ve": [0.5, 0.5],
        "null_lmm": _collect(800, rng, "null", 0.0, n, K, evecs, evals, 0.5, 0.5, False),
    }
    print("null_lmm", flush=True)
    summary["null_ols"] = _collect(400, rng, "null", 0.0, n, K, evecs, evals, 0.5, 0.5, True)
    print("null_ols", flush=True)
    scenarios = (
        ("hap_moderate", "hap", 0.35, 0.5, 0.5),
        ("hap_large", "hap", 0.60, 0.5, 0.5),
        ("hap_dominates", "hap", 0.70, 0.25, 0.25),
        ("snp_large", "snp", 0.55, 0.5, 0.5),
        ("pair", "pair", 0.40, 0.5, 0.5),
    )
    for name, kind, beta, vg, ve in scenarios:
        summary[name] = {
            "beta": beta,
            "vg": vg,
            "ve": ve,
            **_collect(400, rng, kind, beta, n, K, evecs, evals, vg, ve, False),
        }
        print(name, flush=True)
    text = json.dumps(summary, indent=2)
    out = Path(__file__).resolve().parents[1] / "sim_results" / "lrt_bench.json"
    out.write_text(text, encoding="utf-8")
    print(text)
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
