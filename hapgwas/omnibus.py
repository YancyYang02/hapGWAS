"""Cauchy omnibus for haplotype block tests.

Liu & Xie (JASA 2020) ACAT stays valid when the input p-values are
correlated, which joint / collapsed-joint / marginal haplotype tests are.
The default block p-value is ACAT of three design-only scores:

- full multi-df joint (signal spread across alleles)
- ACAT of the 1-df HapAllele marginals (one sharp allele)
- ACAT of the 1-df SNP scores inside the window

ACJ is still computed and written to the table. It is left out of the
default combination because Hamming collapse on a short window often dilutes.

No component looks at the phenotype when it decides the design matrix.
"""

from __future__ import annotations

import numpy as np
from scipy import stats


def _as_finite(pvalues) -> np.ndarray:
    p = np.asarray(pvalues, dtype=float).reshape(-1)
    return p[np.isfinite(p)]


def acat(pvalues) -> float:
    """Equal-weight Cauchy combination. Returns NaN if nothing is finite."""
    p = _as_finite(pvalues)
    if p.size == 0:
        return float("nan")
    # tan((0.5-p)*pi) = cot(p*pi); the ratio form does not overflow at tiny p.
    p = np.clip(p, 1e-300, 1.0 - 1e-15)
    t = np.sum(np.cos(p * np.pi) / np.sin(p * np.pi)) / p.size
    if not np.isfinite(t):
        return float(0.0 if t > 0 else 1.0)
    return float(0.5 - np.arctan(t) / np.pi)


def bonf_min(pvalues) -> tuple[float, float]:
    """Return (minimum raw p, Bonferroni p capped at 1)."""
    p = _as_finite(pvalues)
    if p.size == 0:
        return float("nan"), float("nan")
    raw = float(np.min(p))
    return raw, float(min(1.0, raw * p.size))


def lambda_gc_from_p(pvalues) -> float:
    """
    Genomic inflation from raw p-values, mapped onto a 1-df chi-square.
    Do not pass Bonferroni-adjusted p-values: that pushes lambda far below 1
    even when the raw test is calibrated.
    """
    p = _as_finite(pvalues)
    p = p[(p > 0.0) & (p < 1.0)]
    if p.size < 20:
        return float("nan")
    chi = stats.chi2.isf(p, 1)
    return float(np.median(chi) / stats.chi2.ppf(0.5, 1))


def clump_overlapping(df, pcol: str = "omni_p"):
    """
    Greedy locus list. Walk blocks from smallest p and drop any later block
    whose SNP index interval overlaps a kept block on the same chromosome.
    Overlapping windows of one QTL then count as one locus.
    """
    import pandas as pd

    if df is None or len(df) == 0 or pcol not in df.columns:
        return df
    work = df.sort_values(pcol, kind="mergesort").reset_index(drop=True)
    kept: list[tuple[str, int, int]] = []
    keep_rows: list[int] = []
    for i, row in work.iterrows():
        chrom = str(row["chrom"])
        start, end = int(row["snp_start"]), int(row["snp_end"])
        overlaps = any(c == chrom and start < e and s < end for c, s, e in kept)
        if overlaps:
            continue
        kept.append((chrom, start, end))
        keep_rows.append(int(i))
    out = work.loc[keep_rows].copy()
    out.insert(0, "locus_id", [f"L{i + 1}" for i in range(len(out))])
    return out
