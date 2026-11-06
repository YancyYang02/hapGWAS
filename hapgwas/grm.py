"""VanRaden SNP-GRM used as the default kinship in hapGWAS.

This is the same estimator GEMMA/GCTA users expect. It is not a new
relatedness model. Leave-one-chromosome-out is done by the CLI, which
passes dosages from the chromosomes that are not being tested.
"""

from __future__ import annotations

import numpy as np


def vanraden_factors(dosage: np.ndarray, maf_min: float = 0.01) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Return Z, keep, K for a VanRaden GRM.

    Z columns are the centered and scaled dosages that survived the MAF filter.
    keep gives their positions in the original dosage matrix. K = Z Z^T / m,
    in float64, so removing a set of columns from Z reproduces the
    corresponding leave-those-SNPs-out GRM exactly.
    """
    X = np.asarray(dosage, dtype=np.float64)
    if X.ndim != 2:
        raise ValueError("dosage must be (n, m)")
    if X.shape[0] < 2 or X.shape[1] < 1:
        raise ValueError("need at least 2 samples and 1 SNP to build a GRM")
    af = np.nanmean(X, axis=0) / 2.0
    keep_mask = np.isfinite(af) & (af > maf_min) & (af < 1.0 - maf_min)
    if int(keep_mask.sum()) < 5:
        raise ValueError(
            f"only {int(keep_mask.sum())} SNPs pass MAF>{maf_min}; cannot build a GRM"
        )
    keep = np.flatnonzero(keep_mask)
    X = X[:, keep]
    p = af[keep]
    Z = (X - 2.0 * p) / np.sqrt(2.0 * p * (1.0 - p))
    m = Z.shape[1]
    K = (Z @ Z.T) / m
    K = 0.5 * (K + K.T)
    return Z, keep, K


def load_gcta_grm(prefix: str) -> tuple[list[str], np.ndarray]:
    """Read a GCTA ``--make-grm`` prefix.

    Returns the IID column of ``.grm.id`` and the symmetric GRM. The
    ``.grm.bin`` file is the lower triangle, diagonal included, stored as
    float32 in row-major order. This does not invent SNP factors: a later
    region downdate needs the dosages that built the matrix, which this
    file does not contain.
    """
    id_path = f"{prefix}.grm.id"
    bin_path = f"{prefix}.grm.bin"
    ids: list[str] = []
    with open(id_path, encoding="utf-8") as handle:
        for line in handle:
            parts = line.split()
            if len(parts) >= 2:
                ids.append(parts[1])
    n = len(ids)
    if n < 2:
        raise ValueError(f"{id_path} has {n} samples")
    raw = np.fromfile(bin_path, dtype=np.float32)
    expect = n * (n + 1) // 2
    if int(raw.size) != expect:
        raise ValueError(f"{bin_path} has {raw.size} values, expected {expect} for {n} samples")
    kinship = np.zeros((n, n), dtype=np.float64)
    kinship[np.tril_indices(n)] = raw
    kinship = kinship + kinship.T - np.diag(np.diag(kinship))
    return ids, kinship


def vanraden_grm(dosage: np.ndarray, maf_min: float = 0.01) -> np.ndarray:
    """
    dosage: (n_samples, n_snps) coded 0/1/2.
    Returns a symmetric GRM with mean diagonal near 1 when SNPs are centered
    by their own allele frequencies.
    """
    _Z, _keep, K = vanraden_factors(dosage, maf_min=maf_min)
    if _Z.shape[0] * _Z.shape[1] > 2_000_000:
        Z32 = np.asarray(_Z, dtype=np.float32)
        K32 = (Z32 @ Z32.T).astype(np.float64) / _Z.shape[1]
        return 0.5 * (K32 + K32.T)
    return K
