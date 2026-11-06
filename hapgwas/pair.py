"""Phase-pair (compound-het) indicator score.

Optional 1-df block test. At the first two SNPs of a block, the indicator is
1 when the two haplotypes are alternate/reference and reference/alternate
(Ab in trans with aB). That indicator is residualized on the additive
haplotype dosages in the same block, then scored against the whitened
phenotype. The projection is the same one used by the joint and beyond scores,
so an LMM or OLS null carries over unchanged.

The column to report is this indicator p-value. Multi-df scores that put the
indicator beside the haplotype or SNP columns called far more windows on the
broiler463 BW6 scan and are not part of this test.
"""

from __future__ import annotations

import numpy as np
from scipy import linalg

from hapgwas.adaptive_score import _qr_joint_score
from hapgwas.lmm import NullFit, project_features


def trans_compound_indicator(hap1: np.ndarray, hap2: np.ndarray) -> np.ndarray:
    """Trans compound heterozygote at the first two SNPs.

    hap1, hap2: (n_samples, n_snps), alleles coded 0/1.
    Returns a length-n float indicator. Requires at least two SNP columns.
    A sample is 1 only when one haplotype is (1, 0) and the other is (0, 1).
    """
    hap1 = np.asarray(hap1)
    hap2 = np.asarray(hap2)
    if hap1.shape != hap2.shape or hap1.ndim != 2 or hap1.shape[1] < 2:
        raise ValueError("trans indicator needs two phased SNP columns of the same shape")
    a1 = hap1[:, 0]
    b1 = hap1[:, 1]
    a2 = hap2[:, 0]
    b2 = hap2[:, 1]
    trans = ((a1 == 1) & (b1 == 0) & (a2 == 0) & (b2 == 1)) | (
        (a1 == 0) & (b1 == 1) & (a2 == 1) & (b2 == 0)
    )
    return trans.astype(float)


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
