"""
Unified block-level multi-df joint test + per-HapAllele dosage tests.

This is the statistical core claimed as hapGWAS's methods contribution:
both tests share one null-LMM projection (see lmm.project_features), so
block evidence and allele decomposition are commensurate.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy import linalg, stats

from hapgwas.dosage import HapDosage, build_dosage, haplotypes_from_phased_snps
from hapgwas.window_kinship import WindowBasis, lowo_joint_lrt, region_columns
from hapgwas.blocks import Block
from hapgwas.lmm import NullFit, project_features
from hapgwas.adaptive_score import _qr_joint_score, multi_allelic_scores
from hapgwas.exact_lrt import block_exact_lrts
from hapgwas.omnibus import acat, bonf_min
from hapgwas.pair import pair_indicator_score, trans_compound_indicator
from hapgwas.tree_kernel import tree_kernel_score


@dataclass
class BlockResult:
    block_id: str
    chrom: str
    snp_start: int
    snp_end: int
    n_snps: int
    n_alleles_tested: int
    n_raw_alleles: int
    n_collapsed: int
    collapse_rule: str
    joint_chisq: float
    joint_df: int
    joint_p: float
    allele_rows: list[dict]
    cond_snp_index: int | None = None
    cond_snp_id: str | None = None
    cond_snp_p: float | None = None
    cond_joint_chisq: float | None = None
    cond_joint_df: int | None = None
    cond_joint_p: float | None = None
    acj_chisq: float | None = None
    acj_df: int | None = None
    acj_p: float | None = None
    acj_ngroups: int | None = None
    acj_merges: int | None = None
    pca_chisq: float | None = None
    pca_df: int | None = None
    pca_p: float | None = None
    ridge_stat: float | None = None
    ridge_df: float | None = None
    ridge_p: float | None = None
    allele_min_p: float | None = None
    allele_bonf_p: float | None = None
    allele_acat_p: float | None = None
    snp_min_p: float | None = None
    snp_bonf_p: float | None = None
    snp_acat_p: float | None = None
    omni_p: float | None = None
    tree_qstat: float | None = None
    tree_df: float | None = None
    tree_p: float | None = None
    beyond_chisq: float | None = None
    beyond_df: int | None = None
    beyond_p: float | None = None
    joint_lrt_stat: float | None = None
    joint_lrt_df: int | None = None
    joint_lrt_p: float | None = None
    beyond_lrt_stat: float | None = None
    beyond_lrt_df: int | None = None
    beyond_lrt_p: float | None = None
    pair_chisq: float | None = None
    pair_df: int | None = None
    pair_p: float | None = None
    pair_n: int | None = None
    pair_snp_i: int | None = None
    pair_snp_j: int | None = None


def _safe_p_chi2(chisq: float, df: int) -> float:
    if df <= 0 or not np.isfinite(chisq):
        return np.nan
    return float(stats.chi2.sf(max(chisq, 0.0), df))


def score_joint_and_marginal(null: NullFit, dosage: np.ndarray) -> tuple[float, int, float, np.ndarray, np.ndarray]:
    """
    Score tests on projected scale.

    Under null, y_tilde ~ approx N(0, I) after whitening+residualizing.
    Joint: chi-square with rank(Z) df; marginal: 1-df per column.
    """
    Z = project_features(null, dosage)
    y = null.y_tilde
    col_ss = np.sum(Z * Z, axis=0)
    keep = col_ss > 1e-8
    if not np.any(keep):
        return 0.0, 0, np.nan, np.full(dosage.shape[1], np.nan), np.full(dosage.shape[1], np.nan)

    joint_chisq, joint_df, joint_p = _qr_joint_score(Z, y)

    marg_chisq = np.full(dosage.shape[1], np.nan)
    marg_p = np.full(dosage.shape[1], np.nan)
    Zy = Z.T @ y
    for j in range(dosage.shape[1]):
        ss = col_ss[j]
        if ss <= 1e-8:
            continue
        c = float((Zy[j] ** 2) / ss)
        marg_chisq[j] = c
        marg_p[j] = _safe_p_chi2(c, 1)
    return joint_chisq, joint_df, joint_p, marg_chisq, marg_p


def beyond_snp_score(null: NullFit, hap_dosage: np.ndarray, snp_dosage: np.ndarray) -> tuple[float, int, float]:
    """
    Score test for haplotype signal orthogonal to the additive SNP dosages
    in the same block.

    Haplotype dosages span the SNP dosages, plus extra contrasts whenever
    more than one haplotype is present. Residualize the haplotype design on
    the SNP design in the whitened space, then compute the joint score.
    Under a purely additive SNP effect those residuals are null, so the
    p-value stays uniform. A haplotype pattern that is not a linear function
    of the SNP dosages remains in the residual and can be detected.
    """
    H = project_features(null, hap_dosage)
    S = project_features(null, snp_dosage)
    keep = np.sum(S * S, axis=0) > 1e-8
    if np.any(keep):
        # Pivoted QR so a near-zero direction in the middle is not treated as a real SNP.
        Qs, Rs, _piv = linalg.qr(S[:, keep], mode="economic", pivoting=True)
        diag = np.abs(np.diag(Rs)) if Rs.ndim == 2 else np.array([abs(float(Rs))])
        rank = int(np.sum(diag > 1e-8))
        if rank:
            H = H - Qs[:, :rank] @ (Qs[:, :rank].T @ H)
    return _qr_joint_score(H, null.y_tilde)


def score_snp_matrix(null: NullFit, dosage: np.ndarray, chunk: int = 1024) -> tuple[np.ndarray, np.ndarray]:
    """1-df score tests for many SNP dosages. dosage is (n,) or (n, m)."""
    dosage = np.asarray(dosage, dtype=float)
    if dosage.ndim == 1:
        dosage = dosage.reshape(-1, 1)
    n_snps = dosage.shape[1]
    chisq = np.full(n_snps, np.nan)
    y = null.y_tilde
    for start in range(0, n_snps, chunk):
        stop = min(n_snps, start + chunk)
        Z = project_features(null, dosage[:, start:stop])
        ss = np.sum(Z * Z, axis=0)
        zy = Z.T @ y
        ok = ss > 1e-8
        part = np.full(stop - start, np.nan)
        part[ok] = (zy[ok] ** 2) / ss[ok]
        chisq[start:stop] = part
    p = np.full(n_snps, np.nan)
    ok = np.isfinite(chisq)
    p[ok] = stats.chi2.sf(np.clip(chisq[ok], 0.0, None), 1)
    return chisq, p


def score_snp(null: NullFit, snp_dosage: np.ndarray) -> tuple[float, float]:
    """1-df score test for a SNP dosage vector. Returns (chisq, p)."""
    chisq, p = score_snp_matrix(null, snp_dosage)
    return float(chisq[0]), float(p[0])


def condition_on_snp(
    null: NullFit,
    dosage: np.ndarray | HapDosage,
    snp_dosage: np.ndarray,
) -> tuple[float, int, float]:
    """
    SNP-conditioning: residualize hap dosages on a SNP dosage in the same
    whitened space, then re-run joint score test.
    """
    Z = dosage.dosage if isinstance(dosage, HapDosage) else np.asarray(dosage, dtype=float)
    S = snp_dosage.reshape(-1, 1)
    Sp = project_features(null, S)
    Zp = project_features(null, Z)
    ss = float(np.sum(Sp * Sp))
    if ss > 1e-8:
        coef = (Sp.T @ Zp) / ss
        Zc = Zp - Sp @ coef
    else:
        Zc = Zp
    return _qr_joint_score(Zc, null.y_tilde)


def _pick_top_snp_in_block(
    null: NullFit,
    hap1: np.ndarray,
    hap2: np.ndarray,
    block: Block,
) -> tuple[int, float, float]:
    """Return (global_snp_index, chisq, p) for strongest SNP inside the block."""
    best_j, best_c, best_p = block.snp_start, -1.0, np.nan
    for j in range(block.snp_start, block.snp_end):
        dos = hap1[:, j].astype(float) + hap2[:, j].astype(float)
        c, p = score_snp(null, dos)
        if np.isfinite(c) and c > best_c:
            best_j, best_c, best_p = j, c, p
    return best_j, best_c, best_p


def test_block(
    null: NullFit,
    block: Block,
    hap1: np.ndarray,
    hap2: np.ndarray,
    min_freq: float = 0.01,
    min_count: int = 5,
    collapse_mode: str = "freq",
    condition_snp: int | None = None,
    condition_auto_topsnp: bool = False,
    snp_ids: np.ndarray | None = None,
    compute_pair: bool = False,
) -> BlockResult:
    """Run unified tests for one block, optionally with SNP conditioning."""
    h1 = hap1[:, block.snp_start : block.snp_end]
    h2 = hap2[:, block.snp_start : block.snp_end]
    pairs = haplotypes_from_phased_snps(h1, h2)
    hap = build_dosage(
        pairs,
        min_freq=min_freq,
        min_count=min_count,
        collapse_mode=collapse_mode,
        drop_ref=True,
    )
    joint_chisq, joint_df, joint_p, marg_c, marg_p = score_joint_and_marginal(null, hap.dosage)
    ms = multi_allelic_scores(
        null, hap.dosage, allele_labels=hap.alleles, frequencies=hap.frequencies
    )
    _snp_chi, snp_p = score_snp_matrix(null, h1.astype(float) + h2.astype(float))
    allele_min_p, allele_bonf_p = bonf_min(marg_p)
    snp_min_p, snp_bonf_p = bonf_min(snp_p)
    allele_acat_p = acat(marg_p)
    snp_acat_p = acat(snp_p)
    # Default block p-value. ACJ is reported beside it, but it is not in the
    # combination: on short windows it often dilutes. ACAT of the joint score,
    # the allele-marginal ACAT, and the in-window SNP ACAT tracks whichever
    # architecture is actually present (Liu & Xie 2020).
    omni_p = acat([joint_p, allele_acat_p, snp_acat_p])
    Z_hap = project_features(null, hap.dosage)
    tree_qstat, tree_df, tree_p = tree_kernel_score(Z_hap, null.y_tilde, hap.alleles)
    snp_dosage = h1.astype(float) + h2.astype(float)
    beyond_chisq, beyond_df, beyond_p = beyond_snp_score(null, hap.dosage, snp_dosage)
    pair_chisq = pair_p = None
    pair_df = pair_n = pair_snp_i = pair_snp_j = None
    if compute_pair:
        pair_chisq = float("nan")
        pair_df = 0
        pair_p = float("nan")
        pair_n = 0
        if h1.shape[1] >= 2:
            indicator = trans_compound_indicator(h1, h2)
            pair_n = int(np.sum(indicator))
            pair_snp_i = int(block.snp_start)
            pair_snp_j = int(block.snp_start + 1)
            pair_chisq, pair_df, pair_p = pair_indicator_score(null, hap.dosage, indicator)
    (
        joint_lrt_stat,
        joint_lrt_df,
        joint_lrt_p,
        beyond_lrt_stat,
        beyond_lrt_df,
        beyond_lrt_p,
    ) = block_exact_lrts(null, hap.dosage, snp_dosage)

    cond_idx = condition_snp
    cond_snp_p = None
    cond_joint = (None, None, None)
    if condition_auto_topsnp and cond_idx is None:
        cond_idx, _c, cond_snp_p = _pick_top_snp_in_block(null, hap1, hap2, block)
    elif cond_idx is not None:
        dos = hap1[:, cond_idx].astype(float) + hap2[:, cond_idx].astype(float)
        _c, cond_snp_p = score_snp(null, dos)

    if cond_idx is not None:
        # allow conditioning on SNP outside block, but auto mode stays inside
        snp_dos = hap1[:, cond_idx].astype(float) + hap2[:, cond_idx].astype(float)
        cond_joint = condition_on_snp(null, hap.dosage, snp_dos)

    allele_rows = []
    for j, allele in enumerate(hap.alleles):
        row = {
            "block_id": block.block_id,
            "chrom": block.chrom,
            "allele": allele,
            "freq": float(hap.frequencies[j]) if j < len(hap.frequencies) else np.nan,
            "marg_chisq": float(marg_c[j]) if np.isfinite(marg_c[j]) else np.nan,
            "marg_p": float(marg_p[j]) if np.isfinite(marg_p[j]) else np.nan,
            "joint_p_block": joint_p,
            "joint_df_block": joint_df,
        }
        if cond_idx is not None:
            row["cond_snp_index"] = cond_idx
            row["cond_joint_p_block"] = cond_joint[2]
        allele_rows.append(row)

    cond_snp_id = None
    if cond_idx is not None and snp_ids is not None and 0 <= cond_idx < len(snp_ids):
        cond_snp_id = str(snp_ids[cond_idx])

    return BlockResult(
        block_id=block.block_id,
        chrom=block.chrom,
        snp_start=block.snp_start,
        snp_end=block.snp_end,
        n_snps=block.n_snps,
        n_alleles_tested=len(hap.alleles),
        n_raw_alleles=hap.n_raw_alleles,
        n_collapsed=hap.n_collapsed,
        collapse_rule=hap.collapse_rule,
        joint_chisq=joint_chisq,
        joint_df=joint_df,
        joint_p=joint_p,
        allele_rows=allele_rows,
        cond_snp_index=cond_idx,
        cond_snp_id=cond_snp_id,
        cond_snp_p=float(cond_snp_p) if cond_snp_p is not None and np.isfinite(cond_snp_p) else None,
        cond_joint_chisq=cond_joint[0],
        cond_joint_df=cond_joint[1],
        cond_joint_p=cond_joint[2],
        acj_chisq=ms.acj.chisq,
        acj_df=ms.acj.df,
        acj_p=ms.acj.p,
        acj_ngroups=ms.acj.n_groups,
        acj_merges=ms.acj.merge_events,
        pca_chisq=ms.pca.chisq,
        pca_df=ms.pca.df,
        pca_p=ms.pca.p,
        ridge_stat=ms.ridge.stat,
        ridge_df=ms.ridge.df_satt,
        ridge_p=ms.ridge.p,
        allele_min_p=allele_min_p,
        allele_bonf_p=allele_bonf_p,
        allele_acat_p=allele_acat_p,
        snp_min_p=snp_min_p,
        snp_bonf_p=snp_bonf_p,
        snp_acat_p=snp_acat_p,
        omni_p=omni_p,
        tree_qstat=tree_qstat,
        tree_df=tree_df,
        tree_p=tree_p,
        beyond_chisq=beyond_chisq,
        beyond_df=int(beyond_df),
        beyond_p=beyond_p,
        joint_lrt_stat=joint_lrt_stat,
        joint_lrt_df=joint_lrt_df,
        joint_lrt_p=joint_lrt_p,
        beyond_lrt_stat=beyond_lrt_stat,
        beyond_lrt_df=beyond_lrt_df,
        beyond_lrt_p=beyond_lrt_p,
        pair_chisq=None if pair_chisq is None else float(pair_chisq),
        pair_df=None if pair_df is None else int(pair_df),
        pair_p=None if pair_p is None else float(pair_p),
        pair_n=None if pair_n is None else int(pair_n),
        pair_snp_i=pair_snp_i,
        pair_snp_j=pair_snp_j,
    )


def run_genome(
    null: NullFit,
    blocks: list[Block],
    hap1: np.ndarray,
    hap2: np.ndarray,
    min_freq: float = 0.01,
    min_count: int = 5,
    collapse_mode: str = "freq",
    condition_snp: int | None = None,
    condition_auto_topsnp: bool = False,
    snp_ids: np.ndarray | None = None,
    apply_fdr: bool = True,
    lowo_basis: WindowBasis | None = None,
    lowo_flank: int = 0,
    n_snps: int | None = None,
    compute_pair: bool = False,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Returns (block_level_table, allele_level_table)."""
    block_rows = []
    allele_rows: list[dict] = []
    for b in blocks:
        res = test_block(
            null, b, hap1, hap2,
            min_freq=min_freq,
            min_count=min_count,
            collapse_mode=collapse_mode,
            condition_snp=condition_snp,
            condition_auto_topsnp=condition_auto_topsnp,
            snp_ids=snp_ids,
            compute_pair=compute_pair,
        )
        row = {
            "block_id": res.block_id,
            "chrom": res.chrom,
            "snp_start": res.snp_start,
            "snp_end": res.snp_end,
            "n_snps": res.n_snps,
            "n_alleles_tested": res.n_alleles_tested,
            "n_raw_alleles": res.n_raw_alleles,
            "n_collapsed": res.n_collapsed,
            "collapse_rule": res.collapse_rule,
            "joint_chisq": res.joint_chisq,
            "joint_df": res.joint_df,
            "joint_p": res.joint_p,
            "acj_chisq": res.acj_chisq,
            "acj_df": res.acj_df,
            "acj_p": res.acj_p,
            "acj_ngroups": res.acj_ngroups,
            "acj_merges": res.acj_merges,
            "pca_chisq": res.pca_chisq,
            "pca_df": res.pca_df,
            "pca_p": res.pca_p,
            "ridge_stat": res.ridge_stat,
            "ridge_df": res.ridge_df,
            "ridge_p": res.ridge_p,
            "allele_min_p": res.allele_min_p,
            "allele_bonf_p": res.allele_bonf_p,
            "allele_acat_p": res.allele_acat_p,
            "snp_min_p": res.snp_min_p,
            "snp_bonf_p": res.snp_bonf_p,
            "snp_acat_p": res.snp_acat_p,
            "omni_p": res.omni_p,
            "tree_qstat": res.tree_qstat,
            "tree_df": res.tree_df,
            "tree_p": res.tree_p,
            "beyond_chisq": res.beyond_chisq,
            "beyond_df": res.beyond_df,
            "beyond_p": res.beyond_p,
            "joint_lrt_stat": res.joint_lrt_stat,
            "joint_lrt_df": res.joint_lrt_df,
            "joint_lrt_p": res.joint_lrt_p,
            "beyond_lrt_stat": res.beyond_lrt_stat,
            "beyond_lrt_df": res.beyond_lrt_df,
            "beyond_lrt_p": res.beyond_lrt_p,
        }
        if compute_pair:
            row.update(
                {
                    "pair_chisq": res.pair_chisq,
                    "pair_df": res.pair_df,
                    "pair_p": res.pair_p,
                    "pair_n": res.pair_n,
                    "pair_snp_i": res.pair_snp_i,
                    "pair_snp_j": res.pair_snp_j,
                }
            )
        if lowo_basis is not None:
            h1 = hap1[:, b.snp_start : b.snp_end]
            h2 = hap2[:, b.snp_start : b.snp_end]
            pairs = haplotypes_from_phased_snps(h1, h2)
            hap = build_dosage(pairs, min_freq=min_freq, min_count=min_count, collapse_mode=collapse_mode, drop_ref=True)
            cols = region_columns(
                lowo_basis, b.snp_start, b.snp_end, lowo_flank, n_snps if n_snps is not None else hap1.shape[1]
            )
            stat, df, p = lowo_joint_lrt(lowo_basis, hap.dosage, cols)
            row["joint_lowo_stat"] = stat
            row["joint_lowo_df"] = df
            row["joint_lowo_p"] = p
            row["lowo_nsnps"] = int(cols.size)
        if res.cond_snp_index is not None:
            row.update(
                {
                    "cond_snp_index": res.cond_snp_index,
                    "cond_snp_id": res.cond_snp_id,
                    "cond_snp_p": res.cond_snp_p,
                    "cond_joint_chisq": res.cond_joint_chisq,
                    "cond_joint_df": res.cond_joint_df,
                    "cond_joint_p": res.cond_joint_p,
                }
            )
        block_rows.append(row)
        allele_rows.extend(res.allele_rows)

    block_df = pd.DataFrame(block_rows)
    allele_df = pd.DataFrame(allele_rows)
    if apply_fdr:
        block_df, allele_df = apply_block_fdr(block_df, allele_df)
    return block_df, allele_df


def apply_block_fdr(block_df: pd.DataFrame, allele_df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """BH on raw block and allele p-values. Safe to call once after concatenating chromosomes."""
    if len(block_df):
        for src, dst in (
            ("joint_p", "joint_q"),
            ("acj_p", "acj_q"),
            ("pca_p", "pca_q"),
            ("ridge_p", "ridge_q"),
            ("omni_p", "omni_q"),
            ("tree_p", "tree_q"),
            ("beyond_p", "beyond_q"),
            ("joint_lrt_p", "joint_lrt_q"),
            ("joint_lowo_p", "joint_lowo_q"),
            ("beyond_lrt_p", "beyond_lrt_q"),
            ("pair_p", "pair_q"),
            ("cond_joint_p", "cond_joint_q"),
        ):
            if src in block_df.columns:
                block_df = _add_bh(block_df, src, dst)
    if len(allele_df) and "marg_p" in allele_df.columns:
        allele_df = _add_bh(allele_df, "marg_p", "marg_q")
    return block_df, allele_df


def _add_bh(df: pd.DataFrame, pcol: str, qcol: str) -> pd.DataFrame:
    out = df.copy()
    p = out[pcol].to_numpy(dtype=float)
    q = np.full_like(p, np.nan, dtype=float)
    ok = np.isfinite(p)
    if ok.any():
        pv = p[ok]
        m = len(pv)
        order = np.argsort(pv)
        ranked = pv[order]
        qv = ranked * m / (np.arange(1, m + 1))
        qv = np.minimum.accumulate(qv[::-1])[::-1]
        qv = np.clip(qv, 0, 1)
        q_ok = np.empty(m)
        q_ok[order] = qv
        q[ok] = q_ok
    out[qcol] = q
    return out
