"""
Sparse multi-allelic block scores — hapGWAS v0.4 methods core.

(a) Full multi-df joint score
(b) ACJ: Adaptive Collapsed Joint via Hamming≤h hard neighborhood merge
(c) SACJ (NEW): Spectral Adaptive Collapsed Joint — phenotype-blind
    Hamming-affinity graph → normalized Laplacian → eigengap k →
    hard cluster allele dosages → joint χ² (calibrated)
(d) PCA joint + Ridge (Satterthwaite)

Optional ACJ+: phenotype-informed same-sign corr merge with permutation p.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import stats

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
    Q, R = np.linalg.qr(Zk, mode="reduced")
    diag = np.abs(np.diag(R)) if R.ndim == 2 else np.array([abs(float(R))])
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
    corr_thresh: float = 0.9,
    max_df: int | None = None,
) -> ACJResult:
    """
    ACJ (calibrated): Hamming-neighborhood collapse of HapAlleles, then joint χ².

    Formula:
      1. For allele strings s_i, s_j (excluding OTHER), if Hamming(s_i,s_j) ≤ h,
         put them in the same group (union-find).
      2. Dosage of a group G: z_G = Σ_{j∈G} z_j.
      3. Optional: further merge groups with dosage corr ≥ corr_thresh.
      4. Optional: if #groups > max_df, greedily merge highest-corr groups.
      5. T_ACJ = ||P_{Z*} y||² ~ χ²_rank(Z*).

    Uses only genotype/design info (not y) ⇒ χ² null calibration holds under
    the LMM projection model.
    """
    Z = np.asarray(Z, dtype=float)
    y = np.asarray(y, dtype=float).reshape(-1)
    n, k = Z.shape
    if k == 0:
        return ACJResult(0.0, 0, float("nan"), 0, 0, 0, 0, "acj", [], "chi2")

    labels = list(allele_labels) if allele_labels is not None else [f"a{j}" for j in range(k)]
    if len(labels) != k:
        labels = [f"a{j}" for j in range(k)]

    # Step 1: Hamming union-find
    roots = _union_find_merge(
        k, lambda i, j: _hamming(str(labels[i]), str(labels[j])) <= hamming_max
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
    rule = f"acj_hamming(h<={hamming_max},corr>={corr_thresh},max_df={max_df})"
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
class SACJResult:
    chisq: float
    df: int
    p: float
    n_alleles_in: int
    n_groups: int
    k_spectral: int
    eigengap: float
    tau: float
    rule: str
    group_labels: list[str]
    p_method: str = "chi2"


def _hamming_matrix(labels: list[str]) -> np.ndarray:
    k = len(labels)
    H = np.zeros((k, k), dtype=float)
    for i in range(k):
        for j in range(i + 1, k):
            d = float(_hamming(str(labels[i]), str(labels[j])))
            H[i, j] = H[j, i] = d
    return H


def _kmeans(X: np.ndarray, k: int, rng: np.random.Generator, n_iter: int = 40) -> np.ndarray:
    """Simple k-means; returns cluster id per row."""
    n, d = X.shape
    if k <= 1 or n <= 1:
        return np.zeros(n, dtype=int)
    k = min(k, n)
    # init: farthest-point + random
    centers = np.empty((k, d), dtype=float)
    centers[0] = X[int(rng.integers(0, n))]
    for c in range(1, k):
        dist = np.min(((X[:, None, :] - centers[None, :c, :]) ** 2).sum(-1), axis=1)
        dist = np.maximum(dist, 0)
        s = dist.sum()
        if s <= 1e-18:
            centers[c] = X[int(rng.integers(0, n))]
        else:
            centers[c] = X[int(rng.choice(n, p=dist / s))]
    labels = np.zeros(n, dtype=int)
    for _ in range(n_iter):
        dist = ((X[:, None, :] - centers[None, :, :]) ** 2).sum(-1)
        new = np.argmin(dist, axis=1).astype(int)
        if np.array_equal(new, labels):
            break
        labels = new
        for c in range(k):
            m = labels == c
            if np.any(m):
                centers[c] = X[m].mean(axis=0)
            else:
                centers[c] = X[int(rng.integers(0, n))]
    return labels


def spectral_adaptive_collapsed_joint(
    Z: np.ndarray,
    y: np.ndarray,
    *,
    allele_labels: list[str] | None = None,
    tau: float = 1.0,
    h_cap: int = 2,
    k_max: int | None = None,
    min_k: int | None = None,
    rng: np.random.Generator | None = None,
) -> SACJResult:
    """
    SACJ — Spectral Adaptive Collapsed Joint (phenotype-blind, χ²-calibrated).

    Design-only (no y):
      1. Hamming matrix H on HapAllele strings.
      2. Sparsified affinity: A_ij = exp(-H_ij/τ) if 0 < H_ij ≤ h_cap (default 2),
         else 0. OTHER alleles are isolated.
         (ACJ uses only the hard Hamming≤1 graph; SACJ uses soft weights up to h_cap.)
      3. Normalized Laplacian L = I - D^{-1/2} A D^{-1/2}.
      4. Eigengap chooses k ∈ [min_k, k_max] (default min_k=2 when ≥3 alleles).
      5. Embed with first k eigenvectors (row-normalized); k-means → clusters.
      6. If clustering yields a single group while ≥3 alleles exist, fall back to a
         soft spectral basis Z* = Z V_{:,1:m} (drop trivial eigenvector), m=min(k_max,k-1).
      7. Group/basis dosages → joint χ² on projected ỹ.

    Calibration: steps 1–6 are genotype/design only ⇒ T ~ χ²_rank under the null projection.
    """
    Z = np.asarray(Z, dtype=float)
    y = np.asarray(y, dtype=float).reshape(-1)
    rng = rng or np.random.default_rng(0)
    n, k = Z.shape
    if k == 0:
        return SACJResult(0.0, 0, float("nan"), 0, 0, 0, 0.0, tau, "sacj", [], "chi2")

    labels = list(allele_labels) if allele_labels is not None else [f"a{j}" for j in range(k)]
    if len(labels) != k:
        labels = [f"a{j}" for j in range(k)]

    if k == 1:
        chisq, df, p = _qr_joint_score(Z, y)
        return SACJResult(chisq, df, p, 1, 1 if df else 0, 1, 0.0, tau, "sacj_k1", labels[:], "chi2")

    H = _hamming_matrix(labels)
    other = np.array([str(lb) == "OTHER" for lb in labels], dtype=bool)
    A = np.exp(-H / max(tau, 1e-6))
    np.fill_diagonal(A, 0.0)
    A[H > float(h_cap)] = 0.0
    A[H <= 0] = 0.0
    if other.any():
        A[other, :] = 0.0
        A[:, other] = 0.0

    deg = A.sum(axis=1)
    deg_safe = np.where(deg > 1e-12, deg, 1.0)
    Dmh = 1.0 / np.sqrt(deg_safe)
    L = np.eye(k) - (Dmh[:, None] * A * Dmh[None, :])
    L = 0.5 * (L + L.T)
    try:
        evals, evecs = np.linalg.eigh(L)
    except np.linalg.LinAlgError:
        chisq, df, p = _qr_joint_score(Z, y)
        return SACJResult(chisq, df, p, k, df, k, 0.0, tau, "sacj_eigh_fail", labels[:], "chi2")

    evals = np.clip(evals, 0.0, None)
    k_cap = k_max if k_max is not None else min(6, k)
    k_cap = max(1, min(k_cap, k))
    mk = min_k if min_k is not None else (2 if k >= 3 else 1)
    mk = max(1, min(mk, k_cap))

    # eigengap among first k_cap eigenvalues; prefer k>=mk
    gaps = np.diff(evals[: min(k_cap + 1, k)])
    k_star = mk
    gap = 0.0
    if gaps.size:
        # candidates: gap after index i ⇒ k=i+1
        best = None
        for i, g in enumerate(gaps):
            kk = i + 1
            if mk <= kk <= k_cap:
                if best is None or g > best[0]:
                    best = (float(g), kk)
        if best is None:
            # fall back: largest gap overall, clamped
            i = int(np.argmax(gaps))
            k_star = max(mk, min(i + 1, k_cap))
            gap = float(gaps[i])
        else:
            gap, k_star = best

    def soft_basis(m: int):
        # drop trivial (first) eigenvector when possible
        m = max(1, min(m, k - 1 if k > 1 else 1))
        start = 1 if k > 1 and m <= k - 1 else 0
        V = evecs[:, start : start + m]
        return Z @ V

    X = evecs[:, :k_star].copy()
    rn = np.linalg.norm(X, axis=1, keepdims=True)
    rn = np.where(rn > 1e-12, rn, 1.0)
    X = X / rn
    cl = _kmeans(X, k_star, rng)

    if other.any():
        next_id = int(cl.max()) + 1
        for i in np.where(other)[0]:
            cl[i] = next_id
            next_id += 1

    groups: dict[int, list[int]] = {}
    for i, c in enumerate(cl):
        groups.setdefault(int(c), []).append(i)

    mode = "hard"
    if len(groups) <= 1 and k >= 3:
        # soft spectral fallback — still phenotype-blind
        Zstar = soft_basis(min(k_cap, max(mk, 2)))
        glabels = [f"spec{j}" for j in range(Zstar.shape[1])]
        mode = "soft"
        k_star = Zstar.shape[1]
    else:
        cols = [np.sum(Z[:, idxs], axis=1) for idxs in groups.values()]
        glabels = ["+".join(labels[i] for i in idxs) for idxs in groups.values()]
        Zstar = np.column_stack(cols) if cols else np.zeros((n, 0))

    chisq, df, p = _qr_joint_score(Zstar, y)
    rule = f"sacj({mode},tau={tau},h_cap={h_cap},k={k_star},eigengap={gap:.4g})"
    return SACJResult(
        chisq=chisq,
        df=df,
        p=p,
        n_alleles_in=k,
        n_groups=Zstar.shape[1] if Zstar.ndim == 2 else 0,
        k_spectral=int(k_star),
        eigengap=float(gap),
        tau=float(tau),
        rule=rule,
        group_labels=glabels,
        p_method="chi2",
    )



@dataclass
class MultiScoreResult:
    full_chisq: float
    full_df: int
    full_p: float
    acj: ACJResult
    sacj: SACJResult
    pca: PCAResult
    ridge: RidgeResult
    acj_pheno: ACJResult | None = None


def multi_allelic_scores(
    null: NullFit,
    dosage: np.ndarray,
    *,
    allele_labels: list[str] | None = None,
    hamming_max: int = 1,
    acj_corr: float = 0.9,
    acj_max_df: int | None = None,
    pca_var_frac: float = 0.9,
    ridge_lam: float | None = None,
    with_acj_pheno: bool = False,
    acj_pheno_perms: int = 49,
    sacj_tau: float = 1.0,
    sacj_k_max: int | None = None,
) -> MultiScoreResult:
    Z = project_features(null, np.asarray(dosage, dtype=float))
    y = null.y_tilde
    full_chisq, full_df, full_p = full_joint_score(Z, y)
    k = Z.shape[1]
    md = acj_max_df
    if md is None and k >= 10:
        md = 6
    acj = adaptive_collapsed_joint(
        Z, y, allele_labels=allele_labels, hamming_max=hamming_max, corr_thresh=acj_corr, max_df=md
    )
    sacj = spectral_adaptive_collapsed_joint(
        Z, y, allele_labels=allele_labels, tau=sacj_tau, k_max=sacj_k_max
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
        sacj=sacj,
        pca=pca,
        ridge=ridge,
        acj_pheno=acj_pheno,
    )
