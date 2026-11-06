"""Phase-pair (compound-het) scores for consecutive SNP pairs.

Optional block test. For every pair of neighbouring SNPs in the block, the
indicator is 1 when one haplotype is alternate/reference and the other is
reference/alternate (Ab in trans with aB). Each indicator is residualized on
the additive haplotype dosages in the same block and scored against the
whitened phenotype, the same projection the joint and beyond scores use.
The primary column is an equal-weight ACAT of those 1-df p-values.

A single pair is what an earlier trial scored at the first two SNPs only.
That p-value is not the column this module returns. ``pair_minp`` keeps the
smallest pair p-value for debugging. Multi-df scores that put an indicator
beside the haplotype or SNP columns are not part of this test.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import linalg

from hapgwas.adaptive_score import _qr_joint_score
from hapgwas.lmm import NullFit, project_features
from hapgwas.omnibus import acat


def trans_compound_indicator(hap1: np.ndarray, hap2: np.ndarray, i: int = 0, j: int = 1) -> np.ndarray:
    """Trans compound heterozygote at SNP columns ``i`` and ``j``.

    hap1, hap2: (n_samples, n_snps), alleles coded 0/1.
    A sample is 1 only when one haplotype is (1, 0) and the other is (0, 1)
    at those two columns.
    """
    hap1 = np.asarray(hap1)
    hap2 = np.asarray(hap2)
    if hap1.shape != hap2.shape or hap1.ndim != 2:
        raise ValueError("trans indicator needs two phased matrices of the same shape")
    n_snps = hap1.shape[1]
    if not (0 <= i < n_snps and 0 <= j < n_snps) or i == j:
        raise ValueError(f"SNP columns {i}, {j} are outside a block of {n_snps}")
    a1 = hap1[:, i]
    b1 = hap1[:, j]
    a2 = hap2[:, i]
    b2 = hap2[:, j]
    trans = ((a1 == 1) & (b1 == 0) & (a2 == 0) & (b2 == 1)) | (
        (a1 == 0) & (b1 == 1) & (a2 == 1) & (b2 == 0)
    )
    return trans.astype(float)


def consecutive_pairs(n_snps: int) -> list[tuple[int, int]]:
    """Neighbouring SNP index pairs inside one block: (0, 1), (1, 2), ..."""
    if n_snps < 2:
        return []
    return [(i, i + 1) for i in range(n_snps - 1)]


def pair_indicator_score(
    null: NullFit,
    hap_dosage: np.ndarray,
    indicator: np.ndarray,
) -> tuple[float, int, float]:
    """1-df score of `indicator` after removing additive haplotype dosages.

    Both designs are mapped with ``project_features`` into the null's
    whitened residual space. The indicator is then residualized on the
    haplotype columns (pivoted QR). The score uses ``null.y_tilde`` with the
    null residual variance, the same approximation as ``beyond_snp_score``.
    Returns (chisq, df, p). df is 0 and p is NaN when the residual indicator
    has no variation.
    """
    indicator = np.asarray(indicator, dtype=float).reshape(-1)
    if indicator.shape[0] != null.n:
        raise ValueError("indicator length must match the null sample size")
    if not np.isfinite(indicator).all() or float(np.std(indicator)) <= 1e-12:
        return 0.0, 0, float("nan")

    w = project_features(null, indicator.reshape(-1, 1))
    dosage = np.asarray(hap_dosage, dtype=float)
    if dosage.size and dosage.ndim == 2 and dosage.shape[1] > 0:
        projected = project_features(null, dosage)
        keep = np.sum(projected * projected, axis=0) > 1e-8
        if np.any(keep):
            q_factor, r_factor, _piv = linalg.qr(projected[:, keep], mode="economic", pivoting=True)
            diag = np.abs(np.diag(r_factor)) if np.ndim(r_factor) == 2 else np.asarray([abs(float(r_factor))])
            rank = int(np.sum(diag > 1e-8))
            if rank:
                q_keep = q_factor[:, :rank]
                w = w - q_keep @ (q_keep.T @ w)
    return _qr_joint_score(w, null.y_tilde)


@dataclass
class PairBlockScore:
    """ACAT of consecutive trans-pair scores, plus the smallest pair p-value."""

    pair_p: float
    pair_n: int
    pair_minp: float
    pair_snp_i: int | None
    pair_snp_j: int | None


def score_consecutive_pairs(
    null: NullFit,
    hap_dosage: np.ndarray,
    hap1: np.ndarray,
    hap2: np.ndarray,
    snp_start: int = 0,
) -> PairBlockScore:
    """Score every consecutive SNP pair and combine the finite p-values with ACAT.

    ``pair_n`` is the number of pairs with a finite 1-df score. ``pair_minp``
    is the smallest of those p-values, not a Bonferroni adjustment. ``pair_snp_i``
    and ``pair_snp_j`` are the genome indices of that smallest pair.
    ``snp_start`` is the block's first SNP index.
    """
    hap1 = np.asarray(hap1)
    hap2 = np.asarray(hap2)
    n_snps = int(hap1.shape[1]) if hap1.ndim == 2 else 0
    empty = PairBlockScore(float("nan"), 0, float("nan"), None, None)
    if n_snps < 2:
        return empty

    pvals: list[float] = []
    sites: list[tuple[int, int]] = []
    for i, j in consecutive_pairs(n_snps):
        indicator = trans_compound_indicator(hap1, hap2, i, j)
        _chisq, _df, p = pair_indicator_score(null, hap_dosage, indicator)
        if np.isfinite(p):
            pvals.append(float(p))
            sites.append((int(snp_start) + i, int(snp_start) + j))
    if not pvals:
        return empty
    best = int(np.argmin(pvals))
    return PairBlockScore(
        pair_p=float(acat(pvals)),
        pair_n=len(pvals),
        pair_minp=float(pvals[best]),
        pair_snp_i=sites[best][0],
        pair_snp_j=sites[best][1],
    )
