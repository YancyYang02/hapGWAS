"""Haplotype allele dosage coding and rare-haplotype collapse."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass
class HapDosage:
    """Dosage matrix for one block after rare-allele collapse."""

    alleles: list[str]          # column labels (kept + optional OTHER)
    dosage: np.ndarray          # (n_samples, n_alleles) float 0/1/2
    frequencies: np.ndarray     # allele frequencies (count / 2n)
    n_raw_alleles: int
    n_collapsed: int
    collapse_rule: str


def haplotypes_from_phased_snps(geno_hap1: np.ndarray, geno_hap2: np.ndarray) -> np.ndarray:
    """
    Build per-sample haplotype allele IDs from phased SNP alleles.

    geno_hap1, geno_hap2: (n_samples, n_snps) with alleles 0/1 (or ints).
    Returns hap_pair: (n_samples, 2) object array of haplotype strings.
    """
    if geno_hap1.shape != geno_hap2.shape:
        raise ValueError("haplotype matrices must match")
    n, p = geno_hap1.shape
    pairs = np.empty((n, 2), dtype=object)
    for i in range(n):
        h1 = "".join(str(int(x)) for x in geno_hap1[i])
        h2 = "".join(str(int(x)) for x in geno_hap2[i])
        pairs[i, 0] = h1
        pairs[i, 1] = h2
    return pairs


def adaptive_collapse_threshold(
    counts: np.ndarray,
    min_freq: float = 0.01,
    min_count: int = 5,
    mode: str = "freq",
) -> float:
    """
    Adaptive rare threshold.

    mode='freq': use max(min_freq, min_count/(2n)) conceptually via count floor.
    mode='keep_mass': keep alleles from most common until cum freq >= 1 - min_freq,
                      collapse the tail (implemented in build_dosage).
    Returns a frequency cutoff for mode='freq'.
    """
    if mode == "freq":
        return float(min_freq)
    return float(min_freq)


def build_dosage(
    hap_pairs: np.ndarray,
    min_freq: float = 0.01,
    min_count: int = 5,
    collapse_mode: str = "freq",
    drop_ref: bool = True,
) -> HapDosage:
    """
    Code haplotype alleles as 0/1/2 dosages.

    collapse_mode:
      - 'freq': alleles with freq < min_freq OR count < min_count -> OTHER
      - 'keep_mass': keep common alleles until cumulative freq >= 1-min_freq; rest -> OTHER
      - 'none': no collapse

    If drop_ref, omit the most common allele column (baseline) for joint tests
    that include an intercept in the null; dosages for remaining alleles still
    sum toward 2. Joint test uses the kept columns.
    """
    if hap_pairs.ndim != 2 or hap_pairs.shape[1] != 2:
        raise ValueError("hap_pairs must be (n, 2)")
    n = hap_pairs.shape[0]
    flat = hap_pairs.reshape(-1)
    alleles, counts = np.unique(flat, return_counts=True)
    order = np.argsort(-counts)
    alleles = alleles[order]
    counts = counts[order]
    freqs = counts / (2.0 * n)
    n_raw = len(alleles)

    if collapse_mode == "none":
        keep_mask = np.ones(n_raw, dtype=bool)
        rule = "none"
    elif collapse_mode == "keep_mass":
        cum = np.cumsum(freqs)
        # keep while previous cum < 1-min_freq; always keep at least 1
        keep_mask = np.zeros(n_raw, dtype=bool)
        keep_mask[0] = True
        for i in range(1, n_raw):
            if cum[i - 1] < 1.0 - min_freq:
                keep_mask[i] = True
            else:
                break
        # also enforce min_count on kept
        keep_mask &= counts >= min_count
        if not keep_mask.any():
            keep_mask[0] = True
        rule = f"keep_mass(min_freq={min_freq},min_count={min_count})"
    else:  # freq
        thr = adaptive_collapse_threshold(counts, min_freq, min_count, "freq")
        keep_mask = (freqs >= thr) & (counts >= min_count)
        if not keep_mask.any():
            keep_mask[np.argmax(freqs)] = True
        rule = f"freq(min_freq={thr},min_count={min_count})"

    kept = list(alleles[keep_mask].astype(str))
    collapsed = list(alleles[~keep_mask].astype(str))
    n_collapsed = len(collapsed)

    labels = list(kept)
    if n_collapsed:
        labels.append("OTHER")

    # map allele -> column
    col_index = {a: i for i, a in enumerate(kept)}
    other_i = len(kept) if n_collapsed else None

    dosage = np.zeros((n, len(labels)), dtype=float)
    for i in range(n):
        for j in range(2):
            a = str(hap_pairs[i, j])
            if a in col_index:
                dosage[i, col_index[a]] += 1.0
            elif other_i is not None:
                dosage[i, other_i] += 1.0
            else:
                # should not happen
                pass

    # frequencies of columns
    col_freq = dosage.sum(axis=0) / (2.0 * n)

    if drop_ref and dosage.shape[1] >= 2:
        # drop most frequent column as baseline
        ref = int(np.argmax(col_freq))
        keep_cols = [i for i in range(dosage.shape[1]) if i != ref]
        labels = [labels[i] for i in keep_cols]
        dosage = dosage[:, keep_cols]
        col_freq = col_freq[keep_cols]

    return HapDosage(
        alleles=labels,
        dosage=dosage,
        frequencies=col_freq,
        n_raw_alleles=n_raw,
        n_collapsed=n_collapsed,
        collapse_rule=rule,
    )
