"""
Publication-scale simulation for the hapGWAS default and the
haplotype-beyond-SNP score.

Pre-specified before the large run. Betas are effect sizes, not tuned
to crown a winner. Null replicates are used at three alpha levels so the
genome-wide tail can be checked, not only alpha = 0.05.

Scenarios
  null                  no effect
  single_haplotype      one haplotype allele
  additive_snp          one SNP dosage
  nonadditive_pair      two complementary haplotypes, 0110 and 1001

The beyond-SNP score must stay near alpha under null and under a pure
SNP effect. It is allowed to have power only when the haplotype pattern
is not linear in the SNP dosages.
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

HAPS = ["0000", "1111", "0110", "1001", "0011", "1100"]
FREQ = np.repeat(1.0 / len(HAPS), len(HAPS))
METHODS = (
    "omni_p",
    "joint_p",
    "allele_bonf_p",
    "snp_bonf_p",
    "tree_p",
    "beyond_p",
    "acj_p",
)
ALPHAS = (0.05, 0.01, 0.001)
BETA = {
    "single_haplotype": 0.35,
    "additive_snp": 0.22,
    "nonadditive_pair": 0.28,
}


def _sample(rng: np.random.Generator, n: int):
    i1 = rng.choice(len(HAPS), size=n, p=FREQ)
    i2 = rng.choice(len(HAPS), size=n, p=FREQ)

    def pack(idx):
        return np.array([[int(ch) for ch in HAPS[i]] for i in idx], dtype=float)

    return pack(i1), pack(i2), i1, i2


def _dos(index: int, i1, i2) -> np.ndarray:
    return (i1 == index).astype(float) + (i2 == index).astype(float)


def _phenotype(name: str, i1, i2, h1, h2, rng: np.random.Generator, n: int) -> np.ndarray:
    noise = rng.normal(size=n)
    if name == "null":
        return noise
    if name == "single_haplotype":
        return BETA[name] * _dos(2, i1, i2) + noise
    if name == "additive_snp":
        return BETA[name] * (h1[:, 0] + h2[:, 0]) + noise
    if name == "nonadditive_pair":
        return BETA[name] * (_dos(2, i1, i2) + _dos(3, i1, i2)) + noise
    raise ValueError(name)


def _eval(h1, h2, y, K=None) -> dict:
    null = fit_null_lmm(y, K=K, method="lmm" if K is not None else "ols")
    res = test_block(null, Block("b", "1", 0, h1.shape[1]), h1, h2, min_freq=0.01, min_count=5)
    return {m: float(getattr(res, m)) for m in METHODS}


def _rates(pvalues: dict[str, list[float]]) -> dict:
    out = {}
    for method, plist in pvalues.items():
        arr = np.asarray(plist, dtype=float)
        ok = np.isfinite(arr)
        entry = {"n_finite": int(ok.sum()), "n": int(arr.size)}
        for alpha in ALPHAS:
            if ok.sum() == 0:
                entry[f"a{alpha}"] = None
            else:
                entry[f"a{alpha}"] = float(np.mean(arr[ok] < alpha))
        out[method] = entry
    return out


def _collect(reps: int, rng: np.random.Generator, scenario: str, n: int, K=None, evecs=None, evals=None):
    store = {m: [] for m in METHODS}
    for _ in range(reps):
        h1, h2, i1, i2 = _sample(rng, n)
        if K is None:
            y = _phenotype(scenario, i1, i2, h1, h2, rng, n)
        else:
            u = evecs @ (np.sqrt(evals) * rng.normal(size=n))
            y = np.sqrt(0.4) * u + np.sqrt(0.6) * _phenotype(scenario, i1, i2, h1, h2, rng, n)
        got = _eval(h1, h2, y, K=K)
        for m in METHODS:
            store[m].append(got[m])
    return _rates(store)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--null-reps", type=int, default=4000)
    parser.add_argument("--power-reps", type=int, default=1000)
    parser.add_argument("--lmm-reps", type=int, default=2000)
    parser.add_argument("--n", type=int, default=600)
    parser.add_argument("--seed", type=int, default=23)
    args = parser.parse_args()
    rng = np.random.default_rng(args.seed)
    n = args.n

    summary = {
        "n": n,
        "seed": args.seed,
        "betas": BETA,
        "alphas": list(ALPHAS),
        "rank": "pivoted-qr",
        "null_iid": _collect(args.null_reps, rng, "null", n),
    }
    print("null_iid done", flush=True)
    for name in ("single_haplotype", "additive_snp", "nonadditive_pair"):
        summary[name] = _collect(args.power_reps, rng, name, n)
        print(name, "done", flush=True)

    dosage = rng.binomial(2, 0.3, size=(n, 200)).astype(float)
    K = vanraden_grm(dosage)
    evals, evecs = np.linalg.eigh(0.5 * (K + K.T))
    evals = np.clip(evals, 0.0, None)
    summary["null_lmm"] = _collect(args.lmm_reps, rng, "null", n, K=K, evecs=evecs, evals=evals)
    print("null_lmm done", flush=True)
    summary["additive_snp_lmm"] = _collect(
        args.power_reps, rng, "additive_snp", n, K=K, evecs=evecs, evals=evals
    )
    print("additive_snp_lmm done", flush=True)

    text = json.dumps(summary, indent=2)
    out = Path(__file__).resolve().parents[1] / "sim_results" / "paper_sim.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(text, encoding="utf-8")
    print(text)
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
