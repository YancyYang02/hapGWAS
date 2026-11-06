"""
Sparse multi-allelic block scores — hapGWAS v0.3 methods core.

(a) Full multi-df joint score
(b) ACJ: Adaptive Collapsed Joint via Hamming-neighborhood allele merge
    (genotype-driven; χ² calibrated) + optional corr cleanup / df cap
(c) PCA joint + Ridge (Satterthwaite)

Optional ACJ+ : phenotype-informed same-sign corr merge with permutation p.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import linalg, stats

from hapgwas.lmm import NullFit, project_features


def _safe_chi2_sf(chisq: float, df: float) -> float:
    if not np.isfinite(chisq) or df <= 0:
        return float("nan")
    return float(stats.chi2.sf(max(chisq, 0.0), df))


def _qr_joint_score(Z: np.ndarray, y: np.ndarray, col_tol: float = 1e-8) -> tuple[float, int, float]:
    if Z.size == 0 or Z.ndim != 2 or Z.shape[1] == 0:
        return 0.0, 0, float("nan")
    col_ss = np.sum(Z * Z, axis=0)
    keep = col_ss > col_tol
    if not np.any(keep):
        return 0.0, 0, float("nan")
    Zk = Z[:, keep]
    # Pivoting keeps the nonzero directions in the leading columns.
    # Without it, a zero in the middle makes Q[:, :r] include an arbitrary
    # direction outside the column space.
    Q, R, _piv = linalg.qr(Zk, mode="economic", pivoting=True)
    diag = np.abs(np.diag(R)) if np.ndim(R) == 2 else np.array([abs(float(R))])
    r = int(np.sum(diag > col_tol))
    if r == 0:
        return 0.0, 0, float("nan")
    chisq = float(np.sum((Q[:, :r].T @ y) ** 2))
    return chisq, r, _safe_chi2_sf(chisq, r)


def full_joint_score(Z: np.ndarray, y: np.ndarray) -> tuple[float, int, float]:
    return _qr_joint_score(Z, y)


def marginal_z(Z: np.ndarray, y: np.ndarray) -> np.ndarray:
    k = Z.shape[1]
    u = np.full(k, np.nan)
    ss = np.sum(Z * Z, axis=0)
    Zy = Z.T @ y
    for j in range(k):
        if ss[j] > 1e-8:
            u[j] = float(Zy[j] / np.sqrt(ss[j]))
    return u


def _hamming(a: str, b: str) -> int:
    if a == "OTHER" or b == "OTHER":
        return 10**9
    if len(a) != len(b):
        # pad / mismatch → large distance
        m = max(len(a), len(b))
        a = a.ljust(m, "?")
        b = b.ljust(m, "?")
    return sum(ch1 != ch2 for ch1, ch2 in zip(a, b))


def _union_find_merge(n: int, should_merge) -> list[int]:
    parent = list(range(n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra

    for i in range(n):
        for j in range(i + 1, n):
            if should_merge(i, j):
                union(i, j)
    return [find(i) for i in range(n)]


@dataclass
class ACJResult:
    chisq: float
    df: int
    p: float
    n_alleles_in: int
    n_groups: int
    merge_events: int
    soft_dropped: int
    rule: str
    group_labels: list[str]
    p_method: str = "chi2"


def adaptive_collapsed_joint(
    Z: np.ndarray,
    y: np.ndarray,
    *,
    allele_labels: list[str] | None = None,
    hamming_max: int = 1,
    corr_thresh: float = 1.01,
    max_df: int | None = None,
    frequencies: np.ndarray | None = None,
    common_freq: float = 0.10,
) -> ACJResult:
    """
    ACJ (calibrated): collapse only uncommon HapAlleles, then joint χ².

    A 4-SNP window is a hypercube. Union-find on every Hamming-1 pair connects
    almost every haplotype, and the test collapses to a single burden. v0.5
    therefore merges two alleles only when both have frequency < common_freq
    (default 0.10) and Hamming distance ≤ h. Common haplotypes stay separate.

    Uses only genotype/design info (not y), so the χ² null still applies.
    """
    Z = np.asarray(Z, dtype=float)
    y = np.asarray(y, dtype=float).reshape(-1)
    n, k = Z.shape
    if k == 0:
        return ACJResult(0.0, 0, float("nan"), 0, 0, 0, 0, "acj", [], "chi2")

    labels = list(allele_labels) if allele_labels is not None else [f"a{j}" for j in range(k)]
    if len(labels) != k:
        labels = [f"a{j}" for j in range(k)]

    if frequencies is None:
        rare = np.ones(k, dtype=bool)
    else:
        freq = np.asarray(frequencies, dtype=float).reshape(-1)
        if freq.shape[0] != k:
            freq = np.zeros(k)
        rare = freq < common_freq
    # Only uncommon alleles may merge. Common haplotypes are never chained together.
    roots = _union_find_merge(
        k,
        lambda i, j: bool(
            rare[i] and rare[j] and _hamming(str(labels[i]), str(labels[j])) <= hamming_max
        ),
    )
    groups: dict[int, list[int]] = {}
    for i, r in enumerate(roots):
        groups.setdefault(r, []).append(i)
    merge_events = sum(len(v) - 1 for v in groups.values() if len(v) > 1)

    cols = []
    glabels = []
    for _r, idxs in groups.items():
        z = np.sum(Z[:, idxs], axis=1)
        cols.append(z)
        glabels.append("+".join(labels[i] for i in idxs))

    Zg = np.column_stack(cols) if cols else np.zeros((n, 0))

    # Step 2/3: optional high positive-corr merge among groups
    def corr_force_merge(Zcur, labs, thresh, force_to: int | None):
        nonlocal merge_events
        cols_l = [Zcur[:, j].copy() for j in range(Zcur.shape[1])]
        labs_l = list(labs)
        changed = True
        while changed:
            changed = False
            m = len(cols_l)
            if m <= 1:
                break
            if force_to is not None and m <= force_to and thresh is not None:
                # still allow thresh merges
                pass
            best = None
            for i in range(m):
                ss_i = float(cols_l[i] @ cols_l[i])
                if ss_i <= 1e-8:
                    continue
                for j in range(i + 1, m):
                    ss_j = float(cols_l[j] @ cols_l[j])
                    if ss_j <= 1e-8:
                        continue
                    c = float(cols_l[i] @ cols_l[j] / np.sqrt(ss_i * ss_j))
                    if thresh is not None and c < thresh:
                        continue
                    key = c
                    if best is None or key > best[0]:
                        best = (key, i, j)
            do_force = force_to is not None and len(cols_l) > force_to
            if best is None and do_force:
                # pick max corr even if below thresh
                best = None
                for i in range(m):
                    ss_i = float(cols_l[i] @ cols_l[i])
                    if ss_i <= 1e-8:
                        continue
                    for j in range(i + 1, m):
                        ss_j = float(cols_l[j] @ cols_l[j])
                        if ss_j <= 1e-8:
                            continue
                        c = float(cols_l[i] @ cols_l[j] / np.sqrt(ss_i * ss_j))
                        if best is None or c > best[0]:
                            best = (c, i, j)
            if best is None:
                break
            if (thresh is not None and best[0] >= thresh) or do_force:
                _c, i, j = best
                znew = cols_l[i] + cols_l[j]
                lab = f"({labs_l[i]}+{labs_l[j]})"
                cols_l[i] = znew
                labs_l[i] = lab
                cols_l.pop(j)
                labs_l.pop(j)
                merge_events += 1
                changed = True
            else:
                break
        if not cols_l:
            return np.zeros((n, 0)), []
        return np.column_stack(cols_l), labs_l

    Zstar, glabels = corr_force_merge(Zg, glabels, corr_thresh, max_df)
    chisq, df, p = _qr_joint_score(Zstar, y)
    rule = (
        f"acj_rare_hamming(h<={hamming_max},common_freq<{common_freq},"
        f"corr>={corr_thresh},max_df={max_df})"
    )
    return ACJResult(
        chisq=chisq,
        df=df,
        p=p,
        n_alleles_in=k,
        n_groups=Zstar.shape[1] if Zstar.ndim == 2 else 0,
        merge_events=merge_events,
        soft_dropped=0,
        rule=rule,
        group_labels=glabels,
        p_method="chi2",
    )


def adaptive_collapsed_joint_pheno(
    Z: np.ndarray,
    y: np.ndarray,
    *,
    allele_labels: list[str] | None = None,
    hamming_max: int = 1,
    corr_thresh: float = 0.7,
    max_df: int | None = 5,
    soft_z: float | None = None,
    n_perm: int = 99,
    rng: np.random.Generator | None = None,
) -> ACJResult:
    """ACJ+: start from Hamming groups, then same-sign corr merge using y; perm p."""
    Z = np.asarray(Z, dtype=float)
    y = np.asarray(y, dtype=float).reshape(-1)
    rng = rng or np.random.default_rng(0)
    # First Hamming collapse (design-only)
    base = adaptive_collapsed_joint(
        Z, y, allele_labels=allele_labels, hamming_max=hamming_max, corr_thresh=1.01, max_df=None
    )
    # rebuild Zstar from labels... easier: recompute Hamming groups dosages
    # Use full Z with phenotype signed merges on top of Hamming-collapsed design
    # Reconstruct collapsed design via ACJ structure:
    acj0 = adaptive_collapsed_joint(
        Z, np.zeros(Z.shape[0]), allele_labels=allele_labels, hamming_max=hamming_max, corr_thresh=1.01, max_df=None
    )
    # get design by calling internal path again
    tmp = adaptive_collapsed_joint(
        Z, y * 0.0, allele_labels=allele_labels, hamming_max=hamming_max, corr_thresh=1.01, max_df=None
    )
    # Build Zhamm by re-running with y and extracting from a helper:
    # simplest: compute Hamming design here
    labels = list(allele_labels) if allele_labels is not None else [f"a{j}" for j in range(Z.shape[1])]
    k = Z.shape[1]
    roots = _union_find_merge(k, lambda i, j: _hamming(str(labels[i]), str(labels[j])) <= hamming_max)
    groups: dict[int, list[int]] = {}
    for i, r in enumerate(roots):
        groups.setdefault(r, []).append(i)
    cols = [np.sum(Z[:, idxs], axis=1) for idxs in groups.values()]
    Zhamm = np.column_stack(cols) if cols else np.zeros((Z.shape[0], 0))

    if soft_z is None:
        soft_z = float(np.sqrt(2.0 * np.log(max(Zhamm.shape[1], 2))))

    def signed_stat(y_vec):
        if Zhamm.shape[1] == 0:
            return 0.0, 0, 0, []
        u = marginal_z(Zhamm, y_vec)
        keep = [j for j in range(Zhamm.shape[1]) if np.isfinite(u[j]) and abs(u[j]) >= soft_z]
        soft_dropped = Zhamm.shape[1] - len(keep)
        if not keep:
            keep = [j for j in range(Zhamm.shape[1]) if np.isfinite(u[j])]
            soft_dropped = 0
        Zk = Zhamm[:, keep]
        uk = u[keep]
        # signed corr merge
        cols_l = [Zk[:, j].copy() for j in range(Zk.shape[1])]
        labs = [f"g{j}" for j in range(Zk.shape[1])]
        u_l = uk.copy()
        merges = 0
        def pass_once(force=False):
            nonlocal cols_l, labs, u_l, merges
            m = len(cols_l)
            if m <= 1:
                return False
            best = None
            for i in range(m):
                for j in range(i + 1, m):
                    if np.sign(u_l[i]) != np.sign(u_l[j]) or np.sign(u_l[i]) == 0:
                        continue
                    ss_i = float(cols_l[i] @ cols_l[i])
                    ss_j = float(cols_l[j] @ cols_l[j])
                    if ss_i <= 1e-8 or ss_j <= 1e-8:
                        continue
                    c = float(cols_l[i] @ cols_l[j] / np.sqrt(ss_i * ss_j))
                    if not force and c < corr_thresh:
                        continue
                    key = abs(c) * abs(u_l[i]) * abs(u_l[j]) if not force else abs(c)
                    if best is None or key > best[0]:
                        best = (key, i, j)
            if best is None:
                return False
            _, i, j = best
            znew = cols_l[i] + cols_l[j]
            ss = float(znew @ znew)
            unew = float((znew @ y_vec) / np.sqrt(ss)) if ss > 1e-8 else 0.0
            cols_l[i] = znew
            labs[i] = f"({labs[i]}+{labs[j]})"
            u_l[i] = unew
            cols_l.pop(j)
            labs.pop(j)
            u_l = np.delete(u_l, j)
            merges += 1
            return True
        while pass_once(False):
            pass
        if max_df is not None:
            guard = 0
            while len(cols_l) > max_df and guard < 1000:
                if not pass_once(True):
                    break
                guard += 1
        if not cols_l:
            return 0.0, 0, merges, []
        Zs = np.column_stack(cols_l)
        chisq, df, _ = _qr_joint_score(Zs, y_vec)
        return chisq, df, merges, labs

    chisq, df, merges, labs = signed_stat(y)
    ge = 0
    for _ in range(n_perm):
        cp, _, _, _ = signed_stat(rng.permutation(y))
        if cp >= chisq - 1e-12:
            ge += 1
    p = (ge + 1) / (n_perm + 1)
    return ACJResult(
        chisq=chisq,
        df=df,
        p=float(p),
        n_alleles_in=k,
        n_groups=len(labs),
        merge_events=merges,
        soft_dropped=0,
        rule=f"acj_pheno_hamming(h<={hamming_max},B={n_perm})",
        group_labels=labs,
        p_method="perm",
    )


@dataclass
class PCAResult:
    chisq: float
    df: int
    p: float
    n_comp: int
    var_explained: float
    rule: str


def pca_joint_score(Z: np.ndarray, y: np.ndarray, *, var_frac: float = 0.9, max_comp: int | None = None) -> PCAResult:
    Z = np.asarray(Z, dtype=float)
    y = np.asarray(y, dtype=float).reshape(-1)
    if Z.shape[1] == 0:
        return PCAResult(0.0, 0, float("nan"), 0, 0.0, "pca")
    Zc = Z - Z.mean(axis=0, keepdims=True)
    try:
        U, S, _ = np.linalg.svd(Zc, full_matrices=False)
    except np.linalg.LinAlgError:
        return PCAResult(0.0, 0, float("nan"), 0, 0.0, "pca_fail")
    if S.size == 0 or float(np.sum(S**2)) <= 1e-12:
        return PCAResult(0.0, 0, float("nan"), 0, 0.0, "pca_zero")
    ev = S**2
    frac = np.cumsum(ev) / np.sum(ev)
    r = int(np.searchsorted(frac, var_frac) + 1)
    r = max(1, min(r, len(S)))
    if max_comp is not None:
        r = min(r, max_comp)
    chisq = float(np.sum((U[:, :r].T @ y) ** 2))
    return PCAResult(chisq, r, _safe_chi2_sf(chisq, r), r, float(frac[r - 1]), f"pca(var_frac={var_frac},r={r})")


@dataclass
class RidgeResult:
    stat: float
    df_satt: float
    p: float
    lam: float
    rule: str


def ridge_score(Z: np.ndarray, y: np.ndarray, *, lam: float | None = None) -> RidgeResult:
    Z = np.asarray(Z, dtype=float)
    y = np.asarray(y, dtype=float).reshape(-1)
    k = Z.shape[1]
    if k == 0:
        return RidgeResult(0.0, 0.0, float("nan"), 0.0, "ridge")
    G = Z.T @ Z
    evals = np.clip(np.linalg.eigvalsh(G), 0.0, None)
    if lam is None:
        pos = evals[evals > 1e-10]
        lam = float(np.median(pos)) if len(pos) else 1.0
        lam = max(lam, 1e-3)
    A = G + lam * np.eye(k)
    try:
        Zty = Z.T @ y
        sol = np.linalg.solve(A, Zty)
        stat = float(Zty @ sol)
    except np.linalg.LinAlgError:
        return RidgeResult(0.0, 0.0, float("nan"), lam, "ridge_fail")
    w = evals / (evals + lam)
    e1 = float(np.sum(w))
    e2 = float(np.sum(w**2))
    if e1 <= 1e-12 or e2 <= 1e-12:
        return RidgeResult(stat, 0.0, float("nan"), lam, "ridge_degen")
    df = e1 * e1 / e2
    scale = e2 / e1
    return RidgeResult(stat, df, _safe_chi2_sf(stat / scale, df), lam, f"ridge(lam={lam:.4g})")


@dataclass
class MultiScoreResult:
    full_chisq: float
    full_df: int
    full_p: float
    acj: ACJResult
    pca: PCAResult
    ridge: RidgeResult
    acj_pheno: ACJResult | None = None


def multi_allelic_scores(
    null: NullFit,
    dosage: np.ndarray,
    *,
    allele_labels: list[str] | None = None,
    frequencies: np.ndarray | None = None,
    common_freq: float = 0.10,
    hamming_max: int = 1,
    acj_corr: float = 1.01,
    acj_max_df: int | None = None,
    pca_var_frac: float = 0.9,
    ridge_lam: float | None = None,
    with_acj_pheno: bool = False,
    acj_pheno_perms: int = 49,
) -> MultiScoreResult:
    Z = project_features(null, np.asarray(dosage, dtype=float))
    y = null.y_tilde
    full_chisq, full_df, full_p = full_joint_score(Z, y)
    acj = adaptive_collapsed_joint(
        Z, y,
        allele_labels=allele_labels,
        hamming_max=hamming_max,
        corr_thresh=acj_corr,
        max_df=acj_max_df,
        frequencies=frequencies,
        common_freq=common_freq,
    )
    pca = pca_joint_score(Z, y, var_frac=pca_var_frac)
    ridge = ridge_score(Z, y, lam=ridge_lam)
    acj_pheno = None
    if with_acj_pheno:
        acj_pheno = adaptive_collapsed_joint_pheno(
            Z, y, allele_labels=allele_labels, hamming_max=hamming_max, n_perm=acj_pheno_perms
        )
    return MultiScoreResult(
        full_chisq=full_chisq,
        full_df=full_df,
        full_p=full_p,
        acj=acj,
        pca=pca,
        ridge=ridge,
        acj_pheno=acj_pheno,
    )
