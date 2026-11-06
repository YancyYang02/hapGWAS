"""
Exact likelihood ratio for one haplotype block, after a single eigen-decomposition of K.

GEMMA re-estimates delta = ve/vg for each SNP. A haplotype block is not one column:
the additive SNP dosages in the block are linear combinations of the haplotype
dosages (a SNP count is the sum of the haplotypes that carry that allele). The
two models are therefore nested.

  joint:  covariates  versus  covariates + haplotype dosages
  beyond: covariates + SNP dosages  versus  covariates + SNP dosages + haplotype dosages

Both profiles use the ordinary Gaussian ML likelihood, not REML. REML changes the
contrast space when columns are added, so it cannot be compared across these models.
The degrees of freedom are the difference in column rank. One decomposition of K
is shared by every block; each block only profiles a one-dimensional delta on a
skinny design.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import linalg, optimize, stats

from hapgwas.lmm import NullFit


@dataclass
class Profile:
    loglik: float
    delta: float
    rss: float
    logdet: float


def matrix_rank(A: np.ndarray, tol: float = 1e-8) -> int:
    A = np.asarray(A, dtype=float)
    if A.ndim != 2 or A.shape[0] == 0 or A.shape[1] == 0:
        return 0
    norms = np.sqrt(np.sum(A * A, axis=0))
    keep = norms > 1e-10
    if not np.any(keep):
        return 0
    scaled = A[:, keep] / norms[keep]
    _q, r, _piv = linalg.qr(scaled, mode="economic", pivoting=True)
    diag = np.abs(np.diag(r)) if np.ndim(r) == 2 else np.array([abs(float(r))])
    if diag.size == 0:
        return 0
    return int(np.sum(diag > tol))


def _as_2d(Z: np.ndarray, n: int) -> np.ndarray:
    Z = np.asarray(Z, dtype=float)
    if Z.ndim == 1:
        Z = Z.reshape(n, 1)
    if Z.shape[0] != n:
        raise ValueError("design rows must match the sample size")
    return Z


def _weighted_rss(y: np.ndarray, X: np.ndarray, sw: np.ndarray) -> float:
    yw = sw * y
    ynorm2 = float(np.dot(yw, yw))
    if X.shape[1] == 0:
        return max(ynorm2, 1e-18)
    Xw = X * sw[:, None]
    q, r, _piv = linalg.qr(Xw, mode="economic", pivoting=True)
    if r.size == 0:
        return max(ynorm2, 1e-18)
    diag = np.abs(np.diag(r)) if np.ndim(r) == 2 else np.array([abs(float(r))])
    scale = float(diag[0]) if diag.size else 1.0
    rank = int(np.sum(diag > 1e-8 * max(scale, 1.0)))
    if rank == 0:
        return max(ynorm2, 1e-18)
    proj = q[:, :rank].T @ yw
    rss = ynorm2 - float(np.dot(proj, proj))
    return float(max(rss, 1e-18 * max(ynorm2, 1.0)))


def _pack(n: int, rss: float, logdet: float, delta: float) -> Profile:
    loglik = -0.5 * (n * np.log(rss / n) + logdet)
    return Profile(loglik=float(loglik), delta=float(delta), rss=float(rss), logdet=float(logdet))


def profile_ml(y: np.ndarray, X: np.ndarray, evals: np.ndarray | None) -> Profile:
    """Profile Gaussian ML over delta = ve/vg. evals is None for an OLS null."""
    y = np.asarray(y, dtype=float).reshape(-1)
    X = np.asarray(X, dtype=float)
    n = y.shape[0]
    if X.size == 0:
        X = np.zeros((n, 0), dtype=float)
    if evals is None:
        rss = _weighted_rss(y, X, np.ones(n))
        return _pack(n, rss, 0.0, np.inf)

    evals = np.maximum(np.asarray(evals, dtype=float), 0.0)

    def at(delta: float) -> Profile:
        sw = 1.0 / np.sqrt(evals + delta)
        rss = _weighted_rss(y, X, sw)
        logdet = float(np.sum(np.log(evals + delta)))
        return _pack(n, rss, logdet, delta)

    best = max((at(float(delta)) for delta in np.logspace(-8, 8, 33)), key=lambda p: p.loglik)
    lo = np.log(best.delta) - 1.2
    hi = np.log(best.delta) + 1.2

    def neg(log_delta: float) -> float:
        return -at(float(np.exp(log_delta))).loglik

    bracket = optimize.minimize_scalar(neg, bounds=(lo, hi), method="bounded", options={"xatol": 1e-3})
    if bracket.success:
        refined = at(float(np.exp(bracket.x)))
        if refined.loglik > best.loglik:
            best = refined
    return best


def _rotate_features(null: NullFit, Z: np.ndarray) -> np.ndarray:
    Z = _as_2d(Z, null.n)
    evecs = getattr(null, "_evecs", None)
    if evecs is None or null.method == "ols":
        return Z
    return evecs.T @ Z


def _basis(null: NullFit) -> tuple[np.ndarray, np.ndarray, np.ndarray | None]:
    cached = getattr(null, "_lrt_basis", None)
    if cached is not None:
        return cached
    evecs = getattr(null, "_evecs", None)
    if evecs is None or null.method == "ols":
        cached = (null.y, null.X, None)
    else:
        cached = (evecs.T @ null.y, evecs.T @ null.X, np.maximum(null._evals, 0.0))
    null._lrt_basis = cached  # type: ignore[attr-defined]
    return cached


def _covariate_profile(null: NullFit, y: np.ndarray, X: np.ndarray, evals: np.ndarray | None) -> Profile:
    cached = getattr(null, "_cov_ml", None)
    if cached is None:
        cached = profile_ml(y, X, evals)
        null._cov_ml = cached  # type: ignore[attr-defined]
    return cached


def _lrt(alt: Profile, null_profile: Profile, df: int) -> tuple[float, int, float]:
    if df <= 0:
        return 0.0, 0, float("nan")
    stat = 2.0 * (alt.loglik - null_profile.loglik)
    if not np.isfinite(stat):
        return float("nan"), df, float("nan")
    # A shortfall this small is the one-dimensional search stopping early.
    # The alternative design contains the null design, so the maximum cannot be lower.
    if stat < 0.0:
        stat = 0.0
    return float(stat), int(df), float(stats.chi2.sf(stat, df))


def block_exact_lrts(
    null: NullFit,
    hap_dosage: np.ndarray,
    snp_dosage: np.ndarray,
) -> tuple[float, int, float, float, int, float]:
    """
    Returns joint_stat, joint_df, joint_p, beyond_stat, beyond_df, beyond_p.

    joint compares covariates with covariates + haplotypes.
    beyond compares covariates + SNPs with covariates + SNPs + haplotypes.
    """
    y, X, evals = _basis(null)
    H = _rotate_features(null, hap_dosage)
    S = _rotate_features(null, snp_dosage)
    base = _covariate_profile(null, y, X, evals)
    rank_x = matrix_rank(X)

    if H.shape[1] == 0:
        joint = (0.0, 0, float("nan"))
    else:
        joint_df = matrix_rank(np.column_stack([X, H])) - rank_x
        joint = _lrt(profile_ml(y, np.column_stack([X, H]), evals), base, joint_df)

    snp_design = np.column_stack([X, S]) if S.shape[1] else X
    rank_snp = matrix_rank(snp_design)
    if H.shape[1] == 0:
        beyond = (0.0, 0, float("nan"))
    else:
        alt_design = np.column_stack([snp_design, H])
        beyond_df = matrix_rank(alt_design) - rank_snp
        beyond = _lrt(profile_ml(y, alt_design, evals), profile_ml(y, snp_design, evals), beyond_df)
    return (*joint, *beyond)
