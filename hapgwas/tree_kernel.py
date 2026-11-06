"""Local haplotype Brownian-kernel score.

Haplotypes in one block are leaves of a UPGMA tree built from Hamming
distance. A Brownian motion on that tree gives a covariance among
haplotype alleles: alleles that share a recent branch move together.
The score is the SKAT quadratic form with this covariance, and the
p-value is the tail of the corresponding weighted chi-square mixture
(Kuonen saddlepoint, with a Satterthwaite fallback).

The tree is built from the allele strings only. Phenotypes are not used
to decide which haplotypes are similar, so the null remains a mixture
of chi-squares after the LMM projection.

This is a phylogenetic kernel applied to short phased haplotypes. It is
not a new quadratic-form test, and it is not claimed to dominate SNP,
joint, or marginal haplotype tests outside the regime it targets:
several common haplotypes that sit in one clade and carry the same effect,
with pairwise Hamming distance greater than 1 so a neighbor-collapse rule
leaves them separate.
"""

from __future__ import annotations

import numpy as np
from scipy import stats


def hamming(a: str, b: str) -> float:
    if a == "OTHER" or b == "OTHER":
        return float(max(len(a), len(b), 1))
    if len(a) != len(b):
        width = max(len(a), len(b))
        a = a.ljust(width, "?")
        b = b.ljust(width, "?")
    return float(sum(x != y for x, y in zip(a, b)))


def upgma_covariance(labels: list[str]) -> np.ndarray:
    """
    Brownian covariance of haplotype alleles.

    cov_ii = 1 after scaling. cov_ij is the fraction of tree height shared
    by i and j, so recently diverged haplotypes have covariance near 1 and
    a distant OTHER allele has covariance near 0.
    """
    k = len(labels)
    if k == 0:
        return np.zeros((0, 0))
    if k == 1:
        return np.ones((1, 1))

    dist = np.zeros((k, k))
    for i in range(k):
        for j in range(i + 1, k):
            d = hamming(labels[i], labels[j])
            dist[i, j] = dist[j, i] = d

    clusters: list[list[int]] = [[i] for i in range(k)]
    active = list(range(k))
    # mrca[i, j] = UPGMA node height (from the present) of the pair
    mrca = np.zeros((k, k))
    while len(active) > 1:
        best_d = None
        best_pair = (active[0], active[1])
        for ia, i in enumerate(active):
            for j in active[ia + 1 :]:
                d = dist[i, j]
                if best_d is None or d < best_d:
                    best_d = d
                    best_pair = (i, j)
        i, j = best_pair
        height = float(best_d) / 2.0 if best_d is not None else 0.0
        for a in clusters[i]:
            for b in clusters[j]:
                mrca[a, b] = mrca[b, a] = height
        new_id = len(clusters)
        clusters.append(clusters[i] + clusters[j])
        # average-linkage distance to remaining clusters
        dist = np.pad(dist, ((0, 1), (0, 1)))
        ni, nj = len(clusters[i]), len(clusters[j])
        for t in active:
            if t == i or t == j:
                continue
            dist[new_id, t] = dist[t, new_id] = (ni * dist[i, t] + nj * dist[j, t]) / (ni + nj)
        active = [t for t in active if t != i and t != j]
        active.append(new_id)

    root = float(np.max(mrca)) if np.max(mrca) > 0 else 1.0
    cov = root - mrca
    np.fill_diagonal(cov, root)
    cov = cov / root
    cov = 0.5 * (cov + cov.T)
    # Numerical floor so a singular tree still yields a PSD matrix.
    evals, evecs = np.linalg.eigh(cov)
    evals = np.clip(evals, 0.0, None)
    return (evecs * evals) @ evecs.T


def _satterthwaite_sf(q: float, evals: np.ndarray) -> float:
    c1 = float(np.sum(evals))
    c2 = float(np.sum(evals**2))
    if c1 <= 1e-12 or c2 <= 1e-12 or not np.isfinite(q):
        return float("nan")
    df = c1 * c1 / c2
    scale = c2 / c1
    return float(stats.chi2.sf(max(q, 0.0) / scale, df))


def kuonen_sf(q: float, evals: np.ndarray) -> float:
    """
    Tail probability of sum_i λ_i χ²_1 by Kuonen's saddlepoint.
    Falls back to Satterthwaite if the saddlepoint equation is unstable.
    """
    lam = np.asarray(evals, dtype=float)
    lam = lam[np.isfinite(lam) & (lam > 1e-10)]
    if lam.size == 0 or not np.isfinite(q):
        return float("nan")
    q = float(q)
    mean = float(np.sum(lam))
    if q <= 0:
        return 1.0
    # At the mean the saddlepoint is 0 and the formula is singular; use the center directly.
    if abs(q - mean) <= 1e-8 * max(mean, 1.0):
        return 0.5

    lmax = float(np.max(lam))
    tmax = 0.5 / lmax * (1.0 - 1e-8)

    def k1(t: float) -> float:
        return float(np.sum(lam / (1.0 - 2.0 * t * lam)))

    lo, hi = -50.0 / lmax, tmax
    # Bisection: K'(t) is increasing in t.
    if k1(hi) < q:
        return _satterthwaite_sf(q, lam)
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if k1(mid) < q:
            lo = mid
        else:
            hi = mid
    t = 0.5 * (lo + hi)
    denom = 1.0 - 2.0 * t * lam
    if np.any(denom <= 1e-12):
        return _satterthwaite_sf(q, lam)
    k = float(-0.5 * np.sum(np.log(denom)))
    k2 = float(np.sum(2.0 * lam**2 / denom**2))
    if k2 <= 0 or not np.isfinite(k):
        return _satterthwaite_sf(q, lam)
    inner = t * q - k
    if inner <= 0 or not np.isfinite(inner):
        return _satterthwaite_sf(q, lam)
    w = np.sign(t) * np.sqrt(2.0 * inner)
    u = t * np.sqrt(k2)
    if not np.isfinite(w) or abs(w) < 1e-8:
        return _satterthwaite_sf(q, lam)
    z = w + np.log(u / w) / w
    if not np.isfinite(z):
        return _satterthwaite_sf(q, lam)
    return float(stats.norm.sf(z))


def tree_kernel_score(Z: np.ndarray, y: np.ndarray, labels: list[str]) -> tuple[float, float, float]:
    """
    Return (Q, effective_df, p) for the Brownian haplotype kernel.

    Z is the dosage matrix already projected into the null residual space.
    y is the matching whitened phenotype.
    """
    Z = np.asarray(Z, dtype=float)
    y = np.asarray(y, dtype=float).reshape(-1)
    if Z.ndim != 2 or Z.shape[1] == 0 or len(labels) != Z.shape[1]:
        return 0.0, 0.0, float("nan")
    keep = np.sum(Z * Z, axis=0) > 1e-8
    if not np.any(keep):
        return 0.0, 0.0, float("nan")
    Z = Z[:, keep]
    labels = [labels[i] for i, flag in enumerate(keep) if flag]
    cov = upgma_covariance(labels)
    v = Z.T @ y
    q = float(v @ cov @ v)
    gram = Z.T @ Z
    # Nonzero eigenvalues of Z cov Z' equal those of cov @ gram.
    mixed = cov @ gram
    mixed = 0.5 * (mixed + mixed.T)
    evals = np.linalg.eigvalsh(mixed)
    evals = evals[evals > 1e-8]
    if evals.size == 0:
        return q, 0.0, float("nan")
    c1 = float(np.sum(evals))
    c2 = float(np.sum(evals**2))
    df = (c1 * c1 / c2) if c2 > 0 else 0.0
    return q, float(df), kuonen_sf(q, evals)
