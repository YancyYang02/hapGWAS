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
    Simple 1-D line search REML for y ~ N(X b, vg * S + ve I) after rotating
    by eigenvectors of K, where S = diag(evals).
    """
    n, c = X_rot.shape
    # grid over delta = ve/vg (FaST-LMM parameterization), vg = s2/(delta scale)
    # Use grid on h2 = vg/(vg+ve) in (1e-5, 1-1e-5)
    h2_grid = np.linspace(1e-5, 1 - 1e-5, 40)
    best = None
    for h2 in h2_grid:
        vg = h2
        ve = 1.0 - h2
        w = 1.0 / (vg * evals + ve)
        # weighted least squares
        Xw = X_rot * w[:, None]
        # solve (X' W X) b = X' W y
        try:
            XtWX = X_rot.T @ Xw
            XtWy = X_rot.T @ (w * y_rot)
            beta = linalg.solve(XtWX, XtWy, assume_a="pos")
        except linalg.LinAlgError:
            continue
        r = y_rot - X_rot @ beta
        # REML-ish objective: -0.5 ( log|V| + log|X'WX| + n log(rss) )
        # with V = vg K + ve I on rotated scale -> diag(vg*evals+ve)
        rss = float(np.sum(w * r * r))  # r' W r under ve+vg scale if vg+ve=1
        # rescale: actual variance sigma2 * (h2 K + (1-h2)I); estimate sigma2 = rss / (n-c)
        sigma2 = rss / max(n - c, 1)
        logdet_V = np.sum(np.log(vg * evals + ve)) + n * np.log(sigma2 + 1e-300)
        sign, logdet_XtWX = np.linalg.slogdet(XtWX)
        if sign <= 0:
            continue
        # profiled REML criterion to minimize
        crit = logdet_V + logdet_XtWX + (n - c) * np.log(sigma2 + 1e-300)
        if best is None or crit < best[0]:
            best = (crit, float(sigma2 * vg), float(sigma2 * ve), float(ve / max(vg, 1e-12)))
    if best is None:
        # fallback OLS-like
        r = y_rot - X_rot @ linalg.lstsq(X_rot, y_rot)[0]
        ve = float(np.dot(r, r) / max(n - c, 1))
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
    evals = np.clip(evals, 1e-8, None)  # numerical floor
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
