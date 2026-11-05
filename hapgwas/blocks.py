"""LD / analysis block definitions for hapGWAS."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd


@dataclass
class Block:
    block_id: str
    chrom: str
    snp_start: int  # inclusive index into SNP array
    snp_end: int    # exclusive index
    start_bp: int | None = None
    end_bp: int | None = None

    @property
    def n_snps(self) -> int:
        return self.snp_end - self.snp_start


def fixed_windows(
    n_snps: int,
    window: int = 3,
    step: int | None = None,
    chrom: str = "1",
    positions: np.ndarray | None = None,
) -> list[Block]:
    """Sliding fixed-SNP windows (thesis Hap-assoc style windowing, clean reimpl)."""
    if window < 1:
        raise ValueError("window must be >= 1")
    if step is None:
        step = window  # non-overlapping by default for joint-block tests
    blocks: list[Block] = []
    i = 0
    k = 0
    while i + window <= n_snps:
        j = i + window
        start_bp = int(positions[i]) if positions is not None else None
        end_bp = int(positions[j - 1]) if positions is not None else None
        blocks.append(
            Block(
                block_id=f"{chrom}:w{k}",
                chrom=str(chrom),
                snp_start=i,
                snp_end=j,
                start_bp=start_bp,
                end_bp=end_bp,
            )
        )
        k += 1
        i += step
    return blocks


def blocks_from_frame(df: pd.DataFrame, n_snps: int) -> list[Block]:
    """Accept snp_start/snp_end (0-based, end exclusive) block tables."""
    cols = {c.lower(): c for c in df.columns}

    def col(*names: str) -> str | None:
        for n in names:
            if n in df.columns:
                return n
            if n.lower() in cols:
                return cols[n.lower()]
        return None

    id_c = col("block_id", "id", "block")
    chrom_c = col("chrom", "chr", "chromosome")
    s_c = col("snp_start", "start_idx", "snp_start_idx")
    e_c = col("snp_end", "end_idx", "snp_end_idx")
    if s_c is None or e_c is None:
        raise ValueError("Block file needs snp_start and snp_end (0-based, end exclusive)")

    blocks: list[Block] = []
    for i, row in df.iterrows():
        s = int(row[s_c])
        e = int(row[e_c])
        if e <= s:
            raise ValueError(f"Invalid block range {s}:{e}")
        if e > n_snps:
            raise ValueError(f"snp_end {e} > n_snps {n_snps}")
        bid = str(row[id_c]) if id_c else f"b{i}"
        chrom = str(row[chrom_c]) if chrom_c else "1"
        bp_s = int(row[col("start_bp")]) if col("start_bp") else None
        bp_e = int(row[col("end_bp")]) if col("end_bp") else None
        blocks.append(Block(bid, chrom, s, e, bp_s, bp_e))
    return blocks


def read_plink_blocks(
    path: str | Path,
    snp_ids: np.ndarray | None = None,
    positions: np.ndarray | None = None,
    chrom: np.ndarray | None = None,
) -> list[Block]:
    """
    Read PLINK 1.9 ``--blocks`` table:

        CHR  BP1  BP2  KB  NSNPS  SNPS
        1    100  500  ...  5     rs1|rs2|rs3|rs4|rs5

    Maps SNPS names (or BP span) onto indices using snp_ids / positions.
    """
    path = Path(path)
    df = pd.read_csv(path, sep=r"\s+", engine="python")
    cols = {c.upper(): c for c in df.columns}
    if "SNPS" not in cols and "snp_start" not in {c.lower() for c in df.columns}:
        # maybe already our format
        return blocks_from_frame(df, n_snps=len(snp_ids) if snp_ids is not None else int(df["snp_end"].max()))

    blocks: list[Block] = []
    id_index = {}
    if snp_ids is not None:
        for i, s in enumerate(snp_ids.astype(str)):
            id_index.setdefault(s, i)

    for i, row in df.iterrows():
        snps_field = str(row[cols["SNPS"]]) if "SNPS" in cols else ""
        names = [x for x in snps_field.split("|") if x]
        if names and id_index:
            idxs = [id_index[n] for n in names if n in id_index]
            if len(idxs) < 2:
                continue
            s, e = min(idxs), max(idxs) + 1
        elif positions is not None and "BP1" in cols and "BP2" in cols:
            bp1, bp2 = int(row[cols["BP1"]]), int(row[cols["BP2"]])
            mask = (positions.astype(int) >= bp1) & (positions.astype(int) <= bp2)
            if chrom is not None and "CHR" in cols:
                mask &= chrom.astype(str) == str(row[cols["CHR"]])
            idxs = np.where(mask)[0]
            if len(idxs) < 1:
                continue
            s, e = int(idxs[0]), int(idxs[-1]) + 1
        else:
            raise ValueError("PLINK blocks file needs SNPS matched to snp_ids, or BP1/BP2+positions")
        chrom_s = str(row[cols["CHR"]]) if "CHR" in cols else "1"
        bp_s = int(row[cols["BP1"]]) if "BP1" in cols else None
        bp_e = int(row[cols["BP2"]]) if "BP2" in cols else None
        blocks.append(Block(f"{chrom_s}:plink{i}", chrom_s, s, e, bp_s, bp_e))
    return blocks


def _pairwise_dp_from_phased(hap1: np.ndarray, hap2: np.ndarray, i: int, j: int) -> tuple[float, float, float]:
    """
    Estimate D' and a rough 95% CI using haplotype counts from phased data.

    Returns (Dprime, ci_low, ci_high). CI uses a normal approx on D with MAF floors —
    sufficient for Gabriel-style partitioning MVP (not a full Haploview port).
    """
    # 2n chromosomes
    a = np.concatenate([hap1[:, i], hap2[:, i]]).astype(int)
    b = np.concatenate([hap1[:, j], hap2[:, j]]).astype(int)
    n = len(a)
    # haplotype counts for 00,01,10,11
    c00 = int(np.sum((a == 0) & (b == 0)))
    c01 = int(np.sum((a == 0) & (b == 1)))
    c10 = int(np.sum((a == 1) & (b == 0)))
    c11 = int(np.sum((a == 1) & (b == 1)))
    p_a = (c10 + c11) / n
    p_b = (c01 + c11) / n
    if p_a < 1e-6 or p_a > 1 - 1e-6 or p_b < 1e-6 or p_b > 1 - 1e-6:
        return np.nan, np.nan, np.nan
    p_ab = c11 / n
    d = p_ab - p_a * p_b
    if d >= 0:
        dmax = min(p_a * (1 - p_b), (1 - p_a) * p_b)
    else:
        dmax = min(p_a * p_b, (1 - p_a) * (1 - p_b))
    if dmax <= 0:
        return np.nan, np.nan, np.nan
    dp = d / dmax
    # variance approx for D (Weir)
    var_d = (
        p_a * (1 - p_a) * p_b * (1 - p_b) / n
        + (1 - 2 * p_a) * (1 - 2 * p_b) * d / n
        - d * d / n
    )
    var_d = max(float(var_d), 1e-12)
    se = np.sqrt(var_d) / dmax
    z = 1.96
    lo = max(-1.0, dp - z * se)
    hi = min(1.0, dp + z * se)
    return float(dp), float(lo), float(hi)


def gabriel_blocks(
    hap1: np.ndarray,
    hap2: np.ndarray,
    *,
    positions: np.ndarray | None = None,
    chrom: str = "1",
    max_block_snps: int = 40,
    max_span_bp: int | None = 1_000_000,
    strong_lower: float = 0.70,
    strong_upper: float = 0.98,
    recomb_upper: float = 0.90,
    min_strong_frac: float = 0.95,
    min_snps: int = 2,
) -> list[Block]:
    """
    Gabriel et al.–style CI block partition on phased genotypes (MVP).

    Informative pairs classified as:
      - strong LD: |D'| CI lower >= strong_lower and upper >= strong_upper
      - strong recombination: CI upper < recomb_upper
    A candidate span is a block if among informative pairs, the fraction that
    are strong-LD is >= min_strong_frac (Haploview/Gabriel default ~0.95).

    This is a practical reimplementation for hapGWAS block input — not a
    bit-for-bit Haploview/PLINK clone. Prefer ``read_plink_blocks`` when you
    already have PLINK ``--blocks`` output.
    """
    n, p = hap1.shape
    if p == 0:
        return []

    # Precompute pairwise classifications for j>i within a sliding max span
    # To keep O(p * L^2) with L=max_block_snps
    strong = np.zeros((p, p), dtype=np.int8)   # 1 strong LD, -1 recomb, 0 uninformative/nan
    for i in range(p):
        j_max = min(p, i + max_block_snps)
        for j in range(i + 1, j_max):
            if positions is not None and max_span_bp is not None:
                if int(positions[j]) - int(positions[i]) > max_span_bp:
                    break
            _dp, lo, hi = _pairwise_dp_from_phased(hap1, hap2, i, j)
            if not np.isfinite(lo):
                continue
            # use abs for CI on signed D' by folding: compare on |D'| via hi/lo carefully
            # Gabriel uses D' magnitude CI; with signed approx, take abs bounds roughly
            alo, ahi = abs(lo), abs(hi)
            if lo < 0 < hi:
                # CI crosses 0 -> weak
                strong[i, j] = 0
            elif ahi < recomb_upper and alo < strong_lower:
                strong[i, j] = -1
            elif alo >= strong_lower and max(alo, ahi) >= strong_upper:
                strong[i, j] = 1
            elif ahi < recomb_upper:
                strong[i, j] = -1
            else:
                strong[i, j] = 0

    blocks: list[Block] = []
    i = 0
    b_i = 0
    while i < p:
        best_end = i + 1
        j_lim = min(p, i + max_block_snps)
        for end in range(i + min_snps, j_lim + 1):
            if positions is not None and max_span_bp is not None:
                if int(positions[end - 1]) - int(positions[i]) > max_span_bp:
                    break
            n_strong = 0
            n_recomb = 0
            n_info = 0
            for a in range(i, end):
                for b in range(a + 1, end):
                    v = strong[a, b]
                    if v == 1:
                        n_strong += 1
                        n_info += 1
                    elif v == -1:
                        n_recomb += 1
                        n_info += 1
            if n_info == 0:
                continue
            frac = n_strong / n_info
            if frac >= min_strong_frac and n_recomb == 0:
                best_end = end
        if best_end - i >= min_snps:
            start_bp = int(positions[i]) if positions is not None else None
            end_bp = int(positions[best_end - 1]) if positions is not None else None
            blocks.append(
                Block(
                    block_id=f"{chrom}:g{b_i}",
                    chrom=str(chrom),
                    snp_start=i,
                    snp_end=best_end,
                    start_bp=start_bp,
                    end_bp=end_bp,
                )
            )
            b_i += 1
            i = best_end
        else:
            # singleton / no block — skip one SNP (optionally emit 1-SNP blocks; we skip)
            i += 1
    return blocks


def summarize_blocks(blocks: list[Block]) -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "block_id": b.block_id,
                "chrom": b.chrom,
                "snp_start": b.snp_start,
                "snp_end": b.snp_end,
                "n_snps": b.n_snps,
                "start_bp": b.start_bp,
                "end_bp": b.end_bp,
            }
            for b in blocks
        ]
    )

