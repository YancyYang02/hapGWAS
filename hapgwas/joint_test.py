"""
Unified block-level multi-df joint test + per-HapAllele dosage tests.

Statistical core (v0.4): shared null projection; joint / ACJ / SACJ / PCA / ridge / allele.
both tests share one null-LMM projection (see lmm.project_features), so
block evidence and allele decomposition are commensurate.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy import stats

from hapgwas.dosage import HapDosage, build_dosage, haplotypes_from_phased_snps
from hapgwas.blocks import Block
from hapgwas.lmm import NullFit, project_features
from hapgwas.adaptive_score import multi_allelic_scores


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
    sacj_chisq: float | None = None
    sacj_df: int | None = None
    sacj_p: float | None = None
    sacj_ngroups: int | None = None
    sacj_k: int | None = None


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

    Zk = Z[:, keep]
    Qz, Rz = np.linalg.qr(Zk, mode="reduced")
    diag = np.abs(np.diag(Rz)) if Rz.ndim == 2 else np.array([abs(Rz)])
    rank_mask = diag > 1e-8
    r = int(rank_mask.sum())
    if r == 0:
        joint_chisq, joint_df, joint_p = 0.0, 0, np.nan
    else:
        Qr = Qz[:, :r]
        joint_chisq = float(np.sum((Qr.T @ y) ** 2))
        joint_df = r
        joint_p = _safe_p_chi2(joint_chisq, joint_df)

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


def score_snp(null: NullFit, snp_dosage: np.ndarray) -> tuple[float, float]:
    """1-df score test for a SNP dosage vector. Returns (chisq, p)."""
    Z = project_features(null, snp_dosage.reshape(-1, 1))
    y = null.y_tilde
    ss = float(np.sum(Z * Z))
    if ss <= 1e-8:
        return np.nan, np.nan
    chisq = float((float(np.asarray(Z.T @ y).ravel()[0]) ** 2) / ss)
    return chisq, _safe_p_chi2(chisq, 1)


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
    y = null.y_tilde
    col_ss = np.sum(Zc * Zc, axis=0)
    keep = col_ss > 1e-8
    if not np.any(keep):
        return 0.0, 0, np.nan
    Qz, Rz = np.linalg.qr(Zc[:, keep], mode="reduced")
    diag = np.abs(np.diag(Rz)) if Rz.ndim == 2 else np.array([abs(float(Rz))])
    r = int((diag > 1e-8).sum())
    if r == 0:
        return 0.0, 0, np.nan
    chisq = float(np.sum((Qz[:, :r].T @ y) ** 2))
    return chisq, r, _safe_p_chi2(chisq, r)


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
    ms = multi_allelic_scores(null, hap.dosage, allele_labels=hap.alleles)

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
        sacj_chisq=ms.sacj.chisq,
        sacj_df=ms.sacj.df,
        sacj_p=ms.sacj.p,
        sacj_ngroups=ms.sacj.n_groups,
        sacj_k=ms.sacj.k_spectral,
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
            "sacj_chisq": res.sacj_chisq,
            "sacj_df": res.sacj_df,
            "sacj_p": res.sacj_p,
            "sacj_ngroups": res.sacj_ngroups,
            "sacj_k": res.sacj_k,
        }
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
    if len(block_df):
        block_df = _add_bh(block_df, "joint_p", "joint_q")
        if "acj_p" in block_df.columns:
            block_df = _add_bh(block_df, "acj_p", "acj_q")
        if "pca_p" in block_df.columns:
            block_df = _add_bh(block_df, "pca_p", "pca_q")
        if "ridge_p" in block_df.columns:
            block_df = _add_bh(block_df, "ridge_p", "ridge_q")
        if "sacj_p" in block_df.columns:
            block_df = _add_bh(block_df, "sacj_p", "sacj_q")
        if "cond_joint_p" in block_df.columns:
            block_df = _add_bh(block_df, "cond_joint_p", "cond_joint_q")
    if len(allele_df):
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
