"""
Null LMM fit and FaST-style projection for hapGWAS score tests.

We do NOT claim a new kinship LMM engine. We reuse the standard spectral
trick (Kang et al. / FaST-LMM / GEMMA-style): eigen-decompose K, transform
to independent residuals, then run block tests on the whitened scale.

Novelty lives in joint_test.py (unified block multi-df + allele decomposition
on the same projected residuals), not in reinventing REML.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import linalg


@dataclass
class NullFit:
    """Null-model quantities shared by all haplotype block tests."""

    y: np.ndarray                 # original phenotype (n,)
    X: np.ndarray                 # covariates with intercept (n, c)
    n: int
    # Whitened / projected space for score tests:
    y_tilde: np.ndarray           # residualized & variance-weighted phenotype
    X_qr_rank: int
    sigma2: float                 # residual variance on whitened scale (≈1 if standardized)
    vg: float                     # genetic variance component (0 if OLS-only)
    ve: float                     # environmental variance
    delta: float                  # ve/vg when vg>0 else inf
    method: str                   # 'ols' | 'lmm-spectral'
    # Optional full whitening matrices for advanced use
    weights: np.ndarray | None    # (n,) diagonal weight after eigen rotate; None if OLS


def _add_intercept(X: np.ndarray | None, n: int) -> np.ndarray:
    if X is None:
        return np.ones((n, 1), dtype=float)
    X = np.asarray(X, dtype=float)
    if X.ndim == 1:
        X = X.reshape(-1, 1)
    # ensure intercept
    if not np.allclose(X[:, 0], 1.0):
        X = np.column_stack([np.ones(n), X])
    return X


def fit_null_ols(y: np.ndarray, X: np.ndarray | None = None) -> NullFit:
    """OLS null: project out covariates; unit residual scale for score tests."""
    y = np.asarray(y, dtype=float).reshape(-1)
    n = y.shape[0]
    X = _add_intercept(X, n)
    # QR projection
    Q, R = np.linalg.qr(X, mode="reduced")
    rank = Q.shape[1]
    y_res = y - Q @ (Q.T @ y)
    df = max(n - rank, 1)
    sigma2 = float(np.dot(y_res, y_res) / df)
    # standardize residual for numerical stability
    y_tilde = y_res / np.sqrt(sigma2) if sigma2 > 0 else y_res
    return NullFit(
        y=y,
        X=X,
        n=n,
        y_tilde=y_tilde,
        X_qr_rank=rank,
        sigma2=sigma2,
        vg=0.0,
        ve=sigma2,
        delta=np.inf,
        method="ols",
        weights=None,
    )


def _reml_vg_ve(y_rot: np.ndarray, X_rot: np.ndarray, evals: np.ndarray) -> tuple[float, float]:
    """
    Profile REML on delta = ve/vg after rotating by the eigenvectors of K.

    V = vg * (K + delta I). The profiled criterion (constants dropped) is

        (n-c) log(RSS) + log|K + delta I| + log|X' (K+delta I)^{-1} X|

    with RSS = (y-Xβ)' (K+delta I)^{-1} (y-Xβ) and vg = RSS/(n-c).
    This is the FaST-LMM / EMMA profile. An earlier draft added an extra
    n*log(sigma^2) term, which biases the heritability grid.
    """
    n, c = X_rot.shape
    evals = np.maximum(np.asarray(evals, dtype=float), 0.0)

    def search(deltas: np.ndarray):
        best = None
        for delta in deltas:
            d = evals + float(delta)
            w = 1.0 / d
            xtwx = X_rot.T @ (X_rot * w[:, None])
            xtwy = X_rot.T @ (w * y_rot)
            try:
                beta = linalg.solve(xtwx, xtwy, assume_a="pos")
            except linalg.LinAlgError:
                continue
            resid = y_rot - X_rot @ beta
            rss = float(np.sum(w * resid * resid))
            if not np.isfinite(rss) or rss <= 1e-12:
                continue
            sign, logdet_x = np.linalg.slogdet(xtwx)
            if sign <= 0 or not np.isfinite(logdet_x):
                continue
            crit = (n - c) * np.log(rss) + float(np.sum(np.log(d))) + logdet_x
            if best is None or crit < best[0]:
                vg = rss / max(n - c, 1)
                best = (crit, float(vg), float(vg * delta), float(delta))
        return best

    best = search(np.logspace(-5, 5, 41))
    if best is not None:
        refined = search(np.logspace(np.log10(best[3]) - 0.6, np.log10(best[3]) + 0.6, 21))
        if refined is not None:
            best = refined
    if best is None:
        resid = y_rot - X_rot @ linalg.lstsq(X_rot, y_rot)[0]
        ve = float(np.dot(resid, resid) / max(n - c, 1))
        return 0.0, ve
    return best[1], best[2]


def fit_null_lmm(
    y: np.ndarray,
    X: np.ndarray | None = None,
    K: np.ndarray | None = None,
    method: str = "auto",
) -> NullFit:
    """
    Fit null model.

    method:
      - 'ols': ignore K
      - 'lmm': require K, spectral LMM
      - 'auto': use LMM if K provided else OLS
    """
    y = np.asarray(y, dtype=float).reshape(-1)
    n = y.shape[0]
    X = _add_intercept(X, n)

    if method == "auto":
        method = "lmm" if K is not None else "ols"
    if method == "ols" or K is None:
        return fit_null_ols(y, X)

    K = np.asarray(K, dtype=float)
    if K.shape != (n, n):
        raise ValueError(f"K shape {K.shape} != ({n},{n})")
    # Symmetrize
    K = 0.5 * (K + K.T)
    evals, evecs = linalg.eigh(K)
    evals = np.clip(evals, 0.0, None)
    y_rot = evecs.T @ y
    X_rot = evecs.T @ X
    vg, ve = _reml_vg_ve(y_rot, X_rot, evals)
    delta = ve / vg if vg > 1e-12 else np.inf

    # Whiten: W^{-1/2} with diag sqrt(vg*evals + ve)
    d = vg * evals + ve
    inv_sqrt = 1.0 / np.sqrt(d)
    y_w = inv_sqrt * y_rot
    X_w = inv_sqrt[:, None] * X_rot
    Q, _R = np.linalg.qr(X_w, mode="reduced")
    rank = Q.shape[1]
    y_tilde_rot = y_w - Q @ (Q.T @ y_w)
    # Map residual back note: tests can stay in rotated+whitened space if
    # dosages are rotated the same way. Store rotation for project_Z.
    sigma2 = float(np.dot(y_tilde_rot, y_tilde_rot) / max(n - rank, 1))
    # standardize
    if sigma2 > 0:
        y_tilde_rot = y_tilde_rot / np.sqrt(sigma2)

    fit = NullFit(
        y=y,
        X=X,
        n=n,
        y_tilde=y_tilde_rot,
        X_qr_rank=rank,
        sigma2=sigma2,
        vg=float(vg),
        ve=float(ve),
        delta=float(delta),
        method="lmm-spectral",
        weights=inv_sqrt,
    )
    # stash eigenbases for projecting haplotype dosages
    fit._evecs = evecs  # type: ignore[attr-defined]
    fit._evals = evals  # type: ignore[attr-defined]
    fit._Q = Q  # type: ignore[attr-defined]
    return fit


def project_features(null: NullFit, Z: np.ndarray) -> np.ndarray:
    """
    Project feature matrix Z (n x k) into the same residualized whitened space
    as null.y_tilde, so score tests are correct under the null LMM/OLS.
    """
    Z = np.asarray(Z, dtype=float)
    if Z.ndim == 1:
        Z = Z.reshape(-1, 1)
    if Z.shape[0] != null.n:
        raise ValueError("Z rows must match sample size")

    if null.method == "ols" or null.weights is None:
        Q, _ = np.linalg.qr(null.X, mode="reduced")
        Z_res = Z - Q @ (Q.T @ Z)
        # scale like y_tilde (divide by sqrt(sigma2) already applied to y)
        return Z_res / np.sqrt(null.sigma2) if null.sigma2 > 0 else Z_res

    evecs = null._evecs  # type: ignore[attr-defined]
    Q = null._Q  # type: ignore[attr-defined]
    inv_sqrt = null.weights
    Z_rot = evecs.T @ Z
    Z_w = inv_sqrt[:, None] * Z_rot
    Z_tilde = Z_w - Q @ (Q.T @ Z_w)
    return Z_tilde / np.sqrt(null.sigma2) if null.sigma2 > 0 else Z_tilde
