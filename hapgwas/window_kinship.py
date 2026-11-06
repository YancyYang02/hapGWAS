"""
Leave-one-region-out kinship for a haplotype window.

Chromosome LOCO drops every marker on the tested chromosome. Leaving the
window inside the GRM lets those SNPs, and the markers in LD with them, absorb
the haplotype signal into the random effect.

This update removes only the window plus a flank. The full GRM is decomposed
once. For each region of r SNPs,

    K_region = a K - B B^T,   a = m / (m - r),

and the Gaussian likelihood at delta = ve/vg is evaluated with the matrix
determinant lemma and the Woodbury identity. The excluded rank is the region
size, not the sample size, so the genome is not decomposed again.

The likelihood ratio is the ordinary profile ML comparison of covariates
versus covariates plus haplotype dosages, under that region-excluded covariance.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import linalg, stats

from hapgwas.exact_lrt import matrix_rank
from hapgwas.grm import vanraden_factors


@dataclass
class WindowBasis:
    evals: np.ndarray
    evecs: np.ndarray
    Z_rot: np.ndarray
    y_rot: np.ndarray
    X_rot: np.ndarray
    m: int
    keep: np.ndarray


def _intercept(X: np.ndarray | None, n: int) -> np.ndarray:
    if X is None:
        return np.ones((n, 1), dtype=float)
    X = np.asarray(X, dtype=float)
    if X.ndim == 1:
        X = X.reshape(-1, 1)
    if X.shape[1] == 0 or not np.allclose(X[:, 0], 1.0):
        X = np.column_stack([np.ones(n), X])
    return X


def prepare_window_basis(
    y: np.ndarray,
    dosage: np.ndarray,
    X: np.ndarray | None = None,
    maf_min: float = 0.01,
) -> WindowBasis:
    """Eigen-decompose the VanRaden GRM of `dosage` once."""
    y = np.asarray(y, dtype=float).reshape(-1)
    Z, keep, K = vanraden_factors(dosage, maf_min=maf_min)
    evals, evecs = linalg.eigh(K)
    evals = np.clip(evals, 0.0, None)
    Xr = _intercept(X, y.shape[0])
    return WindowBasis(
        evals=evals,
        evecs=evecs,
        Z_rot=evecs.T @ Z,
        y_rot=evecs.T @ y,
        X_rot=evecs.T @ Xr,
        m=int(Z.shape[1]),
        keep=keep,
    )


def region_columns(basis: WindowBasis, start: int, end: int, flank: int, n_snps: int) -> np.ndarray:
    """Z-column indices for SNPs in [start - flank, end + flank)."""
    lo = max(0, int(start) - int(flank))
    hi = min(int(n_snps), int(end) + int(flank))
    keep = basis.keep
    # keep is sorted original indices
    left = int(np.searchsorted(keep, lo, side="left"))
    right = int(np.searchsorted(keep, hi, side="left"))
    return np.arange(left, right, dtype=int)


def _pieces(basis: WindowBasis, cols: np.ndarray, delta: float):
    cols = np.asarray(cols, dtype=int)
    rank = int(cols.size)
    m_keep = basis.m - rank
    if m_keep < 5:
        return None
    a = basis.m / m_keep
    h = a * basis.evals + float(delta)
    if np.any(h <= 1e-14):
        return None
    inv_h = 1.0 / h
    if rank:
        B = basis.Z_rot[:, cols] / np.sqrt(m_keep)
        S = (B * inv_h[:, None]).T @ B
        eye = np.eye(rank)
        sign, logdet_s = np.linalg.slogdet(eye - S)
        if sign <= 0 or not np.isfinite(logdet_s):
            return None
        try:
            small = linalg.inv(eye - S, check_finite=False)
        except linalg.LinAlgError:
            return None
    else:
        B = np.zeros((basis.evals.shape[0], 0))
        small = np.zeros((0, 0))
        logdet_s = 0.0
    logdet = float(np.sum(np.log(h)) + logdet_s)
    if not np.isfinite(logdet):
        return None
    return inv_h, B, small, logdet


def _apply(V: np.ndarray, inv_h: np.ndarray, B: np.ndarray, small: np.ndarray) -> np.ndarray:
    HV = inv_h[:, None] * V
    if B.shape[1] == 0:
        return HV
    return HV + inv_h[:, None] * (B @ (small @ (B.T @ HV)))


def _ml(basis: WindowBasis, design_rot: np.ndarray, cols: np.ndarray, delta: float) -> float | None:
    packed = _pieces(basis, cols, delta)
    if packed is None:
        return None
    inv_h, B, small, logdet = packed
    y = basis.y_rot.reshape(-1, 1)
    Vinv_y = _apply(y, inv_h, B, small)
    Vinv_X = _apply(design_rot, inv_h, B, small)
    xtvx = design_rot.T @ Vinv_X
    xtvy = design_rot.T @ Vinv_y
    try:
        beta = linalg.solve(xtvx, xtvy, assume_a="pos", check_finite=False)
    except linalg.LinAlgError:
        beta = linalg.lstsq(xtvx, xtvy)[0]
    resid = y - design_rot @ beta
    rss = float(np.sum(resid * _apply(resid, inv_h, B, small)))
    if not np.isfinite(rss) or rss <= 1e-18:
        return None
    n = basis.y_rot.shape[0]
    return float(-0.5 * (n * np.log(rss / n) + logdet))


def _profile(basis: WindowBasis, design_rot: np.ndarray, cols: np.ndarray) -> float | None:
    best = None
    for delta in np.logspace(-8, 8, 25):
        val = _ml(basis, design_rot, cols, float(delta))
        if val is not None and (best is None or val > best[0]):
            best = (val, float(delta))
    if best is None:
        return None
    lo = np.log(best[1]) - 1.0
    hi = np.log(best[1]) + 1.0
    for log_delta in np.linspace(lo, hi, 15):
        val = _ml(basis, design_rot, cols, float(np.exp(log_delta)))
        if val is not None and val > best[0]:
            best = (val, float(np.exp(log_delta)))
    return best[0]


def lowo_joint_lrt(
    basis: WindowBasis,
    hap_dosage: np.ndarray,
    cols: np.ndarray,
) -> tuple[float, int, float]:
    """
    Likelihood ratio of haplotype dosages under the region-excluded GRM.

    Returns stat, df, p. df is the extra column rank of the haplotypes.
    """
    H = np.asarray(hap_dosage, dtype=float)
    if H.ndim == 1:
        H = H.reshape(-1, 1)
    H_rot = basis.evecs.T @ H
    X = basis.X_rot
    alt = np.column_stack([X, H_rot]) if H_rot.shape[1] else X
    df = matrix_rank(alt) - matrix_rank(X)
    if df <= 0:
        return 0.0, 0, float("nan")
    null_ll = _profile(basis, X, cols)
    alt_ll = _profile(basis, alt, cols)
    if null_ll is None or alt_ll is None:
        return float("nan"), df, float("nan")
    stat = 2.0 * (alt_ll - null_ll)
    if not np.isfinite(stat):
        return float("nan"), df, float("nan")
    if stat < 0.0:
        stat = 0.0
    return float(stat), int(df), float(stats.chi2.sf(stat, df))
