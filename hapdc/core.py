
"""hapdc.core -- mathematical core.

All routines operate on haplotype matrices in (n_haplotype, n_marker) layout,
where the two haplotypes of individual i occupy rows 2i and 2i+1.
"""
import numpy as np
from scipy.special import betainc
from scipy import optimize

__all__ = ["windows", "window_design", "reml_delta", "acat", "t2p", "f2p",
           "lmm_scan_additive", "lambda_gc", "ld_blocks", "jitter_blocks",
           "check_fixed_effects", "FixedEffectRankError"]


def windows(chrom, w_size, w_step):
    """Sliding window start indices over one chromosome, never crossing a boundary."""
    out = []
    i, nm = 0, len(chrom)
    while i < nm:
        j = i
        while j < nm and chrom[j] == chrom[i]:
            j += 1
        out.extend((s, s + w_size) for s in range(i, j - w_size + 1, w_step))
        i = j
    return out


def window_design(H, minf=0.05, minc=5):
    """Classify haplotypes and diplotypes inside one window.

    H : (n_hap, n_marker) uint8 0/1.
    Returns dos (n, nk) additive haplotype dosages, C (n, ncell) diplotype
    cell memberships, nk haplotype classes retained, ncell diplotype cells
    retained, and nrare (haplotype copies reassigned) / nsmall (individuals
    reassigned, summed membership).

    Haplotypes below `minf` are reassigned to the nearest retained haplotype by
    Hamming distance; among equally near retained haplotypes the more frequent
    one is taken. Individuals in diplotype cells below `minc` are moved to the
    most frequent retained cell sharing one of their haplotypes (else to the most
    frequent retained cell). Every individual's haplotype-class dosages sum to
    2 and its cell memberships sum to 1: no implicit group with a free mean is
    created (leaving rare carriers out of all columns would, next to an
    intercept, give them their own parameter).

    Exact ties (equal distance AND equal frequency, or equally frequent
    candidate cells) carry no information for choosing one class. They are not
    broken by the byte order of the patterns, which depends on how each allele
    is coded, but shared equally among the tied candidates: a rare haplotype
    copy tied between t retained classes contributes 1/t to each of their
    dosages, and an individual in a small cell tied between t target cells
    contributes 1/t of its membership to each. Such an individual then has
    fractional dosage / cell membership (the masses above still sum to 2 and 1),
    so with ties it does not belong to exactly one class pair or one cell. The
    design is invariant to recoding the alleles of any marker and to the order
    of the samples. Without ties the result is the 0/1 design.
    """
    H = np.ascontiguousarray(H, dtype=np.uint8)
    nh, nm = H.shape
    n = nh // 2
    if nm == 0:
        return dict(dos=np.zeros((n, 0)), C=np.zeros((n, 0)), nk=0, ncell=0, nrare=0, nsmall=0)
    packed = np.ascontiguousarray(np.packbits(H, axis=1))
    vv = packed.view(np.dtype((np.void, packed.shape[1]))).ravel()
    _, first, inv, cnt = np.unique(vv, return_index=True, return_inverse=True, return_counts=True)
    inv = np.asarray(inv).ravel()
    keep = np.where(cnt / nh >= minf)[0]
    # class order only names the columns; it does not change any column space
    keep = keep[np.argsort(-cnt[keep], kind="stable")]
    nk = len(keep)
    if nk == 0:
        return dict(dos=np.zeros((n, 0)), C=np.zeros((n, 0)), nk=0, ncell=0, nrare=0, nsmall=0)
    A = np.zeros((len(cnt), nk))          # pattern -> class membership
    A[keep, np.arange(nk)] = 1.0
    rare = np.setdiff1d(np.arange(len(cnt)), keep)
    if len(rare):
        Kh = H[first[keep]].astype(np.float64)
        Hr = H[first[rare]].astype(np.float64)
        dist = np.rint(Hr @ (1.0 - Kh).T + (1.0 - Hr) @ Kh.T)
        ck = cnt[keep].astype(np.float64)
        near = dist == dist.min(1, keepdims=True)
        score = np.where(near, ck[None, :], -1.0)
        tied = score == score.max(1, keepdims=True)
        A[rare] = tied / tied.sum(1, keepdims=True)
    nrare = int(cnt[rare].sum()) if len(rare) else 0
    Ah = A[inv]                            # (nh, nk) haplotype-copy membership
    a1, a2 = Ah[0::2], Ah[1::2]
    dos = a1 + a2
    # unordered diplotype cells (k <= l), coded k * nk + l; membership mass
    # a1_k a2_l + a1_l a2_k (k < l) or a1_k a2_k. Rows without ties are one-hot.
    k1, k2 = a1.argmax(1), a2.argmax(1)
    lo, hi = np.minimum(k1, k2), np.maximum(k1, k2)
    rows = [np.arange(n)]; codes = [lo.astype(np.int64) * nk + hi]; mass = [np.ones(n)]
    fr = np.where((a1.max(1) < 1.0) | (a2.max(1) < 1.0))[0]
    if len(fr):
        keepr = np.ones(n, bool); keepr[fr] = False
        rows, codes, mass = [rows[0][keepr]], [codes[0][keepr]], [mass[0][keepr]]
        for i_ in fr:
            ks, ls = np.flatnonzero(a1[i_] > 0), np.flatnonzero(a2[i_] > 0)
            for k_ in ks:
                for l_ in ls:
                    rows.append(np.array([i_])); mass.append(np.array([a1[i_, k_] * a2[i_, l_]]))
                    codes.append(np.array([min(k_, l_) * nk + max(k_, l_)], np.int64))
    rows, codes, mass = np.concatenate(rows), np.concatenate(codes), np.concatenate(mass)
    ucode, ci = np.unique(codes, return_inverse=True)
    M = np.zeros((n, len(ucode)))
    np.add.at(M, (rows, np.asarray(ci).ravel()), mass)
    iu, ju = ucode // nk, ucode % nk
    cc = M.sum(0)
    tol = 1e-9
    live = cc > tol
    big = np.where(live & (cc >= minc - tol))[0]
    small = np.where(live & (cc < minc - tol))[0]
    nsmall = 0
    if len(big) == 0:
        return dict(dos=dos, C=np.zeros((n, 0)), nk=nk, ncell=0, nrare=nrare, nsmall=0)
    for c in small:
        a, b = iu[c], ju[c]
        share = big[(iu[big] == a) | (ju[big] == a) | (iu[big] == b) | (ju[big] == b)]
        cand = share if len(share) else big
        best = cand[np.abs(cc[cand] - cc[cand].max()) <= tol]
        mv = M[:, c].copy()
        M[:, best] += mv[:, None] / len(best)
        M[:, c] = 0.0
        nsmall += float(mv.sum())
    order = big[np.argsort(-cc[big], kind="stable")]
    C = M[:, order]
    return dict(dos=dos, C=C, nk=nk, ncell=len(order), nrare=nrare,
                nsmall=int(round(nsmall)))

class FixedEffectRankError(ValueError):
    """Raised when the fixed-effect design is not of full column rank and the
    caller asked for an error instead of dropping the redundant columns."""


def _seq_resid(X):
    """Sequential Gram-Schmidt residual norms of the unit-scaled columns of X.

    Entry j is the norm of column j / |column j| after projecting out the
    columns before it that were themselves kept (residual > 1e-8). A residual
    of 0 marks a zero column. Invariant to an orthogonal rotation of the rows,
    so the same rule gives the same answer on the original and the
    eigen-rotated design."""
    X = np.asarray(X, float)
    out = np.zeros(X.shape[1])
    basis = []
    for j in range(X.shape[1]):
        x = X[:, j]
        nx = np.linalg.norm(x)
        if nx == 0:
            continue
        r_ = x / nx
        if basis:
            Q = np.column_stack(basis)
            r_ = r_ - Q @ (Q.T @ r_)
            r_ = r_ - Q @ (Q.T @ r_)               # second pass for orthogonality
        out[j] = float(np.linalg.norm(r_))
        if out[j] > 1e-8:
            basis.append(r_ / out[j])
    return out


def check_fixed_effects(covar, n, names=None, policy="drop", tol=1e-8):
    """Check the fixed-effect design [1, covar] before any REML fit.

    covar : (n, c) array or None. The intercept is always included and always
        kept; it is not part of `covar`.
    names : optional column names (length c), used in messages.
    policy : "drop" removes redundant covariate columns, "error" raises
        FixedEffectRankError naming them.
    tol : relative rank tolerance on unit-norm columns (the same rule as the
        column-pivoted basis used by the scan).

    Columns are tested in the order given: a covariate is redundant if, after
    scaling to unit norm, its residual on the intercept and on the covariates
    kept before it has norm <= tol. The result therefore depends only on the
    values and the column order the user supplied, not on rounding of a
    pivoted factorisation. Non-finite values always raise ValueError.

    Returns (covar_kept or None, kept_names, dropped), where dropped is a list
    of dicts {name, index, reason} with reason "constant" (collinear with the
    intercept), "duplicate of <name>" or "linear combination of earlier columns".
    """
    if covar is None:
        return None, [], []
    C = np.asarray(covar, float)
    if C.size == 0:
        return None, [], []
    if C.ndim == 1:
        C = C[:, None]
    if C.shape[0] != n:
        raise ValueError("covariate matrix has %d rows, expected %d" % (C.shape[0], n))
    c = C.shape[1]
    names = list(names) if names is not None else ["covar%d" % (j + 1) for j in range(c)]
    if len(names) != c:
        raise ValueError("%d covariate names for %d columns" % (len(names), c))
    bad = [names[j] for j in range(c) if not np.isfinite(C[:, j]).all()]
    if bad:
        raise ValueError("non-finite values in covariate column(s) %s; remove or impute the "
                         "affected samples before the scan" % ", ".join(bad))
    rr = _seq_resid(np.column_stack([np.ones(n), C]))[1:]
    kept, dropped = [], []
    for j in range(c):
        x = C[:, j]
        r = rr[j]
        if r <= tol:
            if np.ptp(x) == 0:
                reason = "constant"
            else:
                dup = [names[k] for k in kept if np.array_equal(C[:, k], x)]
                reason = ("duplicate of %s" % dup[0]) if dup else "linear combination of earlier columns"
            dropped.append(dict(name=names[j], index=j, reason=reason))
            continue
        kept.append(j)
    if policy not in ("drop", "error"):
        raise ValueError("policy must be 'drop' or 'error'")
    if dropped and policy == "error":
        raise FixedEffectRankError(
            "fixed-effect design [intercept, covariates] is rank deficient: "
            + "; ".join("%s (%s)" % (d["name"], d["reason"]) for d in dropped))
    if not kept:
        return None, [], dropped
    return C[:, kept], [names[j] for j in kept], dropped


def _reml_nll(yr, S, Xr):
    # profile REML negative log-likelihood minimised over delta; returns (nll, delta)
    n = len(yr); Xr = np.asarray(Xr, float)
    X2 = Xr[:, None] if Xr.ndim == 1 else Xr
    q = X2.shape[1]
    def nll(ld):
        d = np.exp(ld); w = 1.0 / (S + d)
        Xw = X2 * w[:, None]; XtX = X2.T @ Xw
        b = np.linalg.solve(XtX, Xw.T @ yr); r = yr - X2 @ b
        s2 = (w * r * r).sum() / (n - q)
        return 0.5 * ((n - q) * np.log(s2) + np.log(S + d).sum() + np.linalg.slogdet(XtX)[1])
    r = optimize.minimize_scalar(nll, bounds=(-10.0, 10.0), method="bounded")
    return float(r.fun), float(np.exp(r.x))


KPOW_GRID = (0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 2.5, 3.0, 4.0)


def kernel_power(yr, S, Xr, grid=KPOW_GRID, crit=3.841458820694124):
    # Spectral-power kernel selection. K = U diag(S) U' is replaced by
    # K_g = U diag(S**g / mean(S**g)) U'; g = 1 is the standard GRM. g is chosen by
    # profile REML on a grid and departs from 1 only if 2*(l(g_hat) - l(1)) > crit
    # (chi-square 1 df, 0.05). Returns (g, lrt, table).
    tab = []
    for g in grid:
        Sg = S ** g; Sg = Sg / Sg.mean()
        f, d = _reml_nll(yr, Sg, Xr); tab.append((float(g), -f, d))
    l1 = [t for t in tab if t[0] == 1.0][0][1]
    gb, lb, _ = max(tab, key=lambda t: t[1])
    lrt = 2.0 * (lb - l1)
    return (gb if lrt > crit else 1.0), float(lrt), tab


def reml_delta(yr, S, Xr):
    """REML variance ratio for the GRM random effect, on the eigen scale.

    Xr is the rotated fixed-effect design: a vector (intercept only) or an
    (n, q) matrix (intercept plus covariates). The design must have full column
    rank; check it with `check_fixed_effects` first. A rank-deficient matrix
    raises FixedEffectRankError here instead of a bare LinAlgError."""
    n = len(yr)
    Xr = np.asarray(Xr, float)
    if Xr.ndim == 2 and Xr.shape[1] > 1 and (_seq_resid(Xr) <= 1e-8).any():
            raise FixedEffectRankError(
                "fixed-effect design passed to REML is rank deficient (%d columns); "
                "run core.check_fixed_effects on the covariates first" % Xr.shape[1])
    if Xr.ndim == 1:
        def nll(ld):
            d = np.exp(ld)
            w = 1.0 / (S + d)
            Xw = Xr * w
            b = (Xw @ yr) / (Xw @ Xr)
            r = yr - Xr * b
            s2 = (w * r * r).sum() / (n - 1)
            return 0.5 * ((n - 1) * np.log(s2) + np.log(S + d).sum() + np.log(Xw @ Xr))
    else:
        q = Xr.shape[1]

        def nll(ld):
            d = np.exp(ld)
            w = 1.0 / (S + d)
            Xw = Xr * w[:, None]
            XtX = Xr.T @ Xw
            b = np.linalg.solve(XtX, Xw.T @ yr)
            r = yr - Xr @ b
            s2 = (w * r * r).sum() / (n - q)
            return 0.5 * ((n - q) * np.log(s2) + np.log(S + d).sum() + np.linalg.slogdet(XtX)[1])

    r = optimize.minimize_scalar(nll, bounds=(-10.0, 10.0), method="bounded")
    return float(np.exp(r.x))


def _T(p):
    p = np.clip(p, 1e-300, 1 - 1e-15)
    return np.where(p < 1e-10, 1.0 / (p * np.pi), np.tan((0.5 - p) * np.pi))


def acat(P, w):
    """ACAT: P (m,B) p-values along axis 0, w (m,) weights -> (B,) combined p.

    The closed form 0.5 - arctan(S)/pi loses every significant digit once S is
    large, because arctan saturates at pi/2 in float64. That used to clip real
    extreme p-values to the 1e-300 floor. The exact asymptotic tail 1/(pi*S) is
    used instead wherever it applies.
    """
    P = np.atleast_2d(P)
    w = np.asarray(w, float)
    w = w / w.sum()
    S = (w[:, None] * _T(P)).sum(0)
    small = 0.5 - np.arctan(S) / np.pi
    large = 1.0 / (np.pi * np.maximum(S, 1e-300))
    return np.clip(np.where(S > 1e7, large, small), 1e-300, 1.0)


def t2p(t, df):
    return betainc(0.5 * df, 0.5, df / (df + t * t))


def f2p(F, d1, d2):
    return betainc(0.5 * d2, 0.5 * d1, d2 / (d2 + d1 * np.clip(F, 0, None)))


def lmm_scan_additive(Gt, Yt, l0, dlt, one_t):
    """Vectorised additive scan on the eigen scale.

    Gt : (n, m) eigen-rotated dosages (U.T @ dosage).
    Yt : (n, B) eigen-rotated phenotypes, column 0 = observed.
    one_t : (n,) eigen-rotated intercept, U.T @ 1. Required: in the rotated
        space the intercept is not a constant vector, and projecting out the
        weights alone (as earlier versions did) leaves the mean in the residual
        and can produce spuriously small p-values.
    Returns t statistics (m, B) and their degrees of freedom (n - 2).
    """
    n = Gt.shape[0]
    w = 1.0 / np.sqrt(l0 + dlt)
    X = Gt * w[:, None]
    Yw = Yt * w[:, None]
    o = np.asarray(one_t, float).ravel() * w
    q0 = o / np.linalg.norm(o)
    Xr = X - np.outer(q0, q0 @ X)
    Yr = Yw - np.outer(q0, q0 @ Yw)
    nrm = (Xr * Xr).sum(0)
    ok = nrm > 1e-10
    b = np.zeros((Gt.shape[1], Yt.shape[1]))
    b[ok] = (Xr[:, ok].T @ Yr) / nrm[ok][:, None]
    rss0 = (Yr * Yr).sum(0)
    rss = np.clip(rss0[None, :] - b * b * nrm[:, None], 0, None)
    df = n - 2
    t = np.zeros_like(b)
    good = ok[:, None] & (rss > 0)
    t[good] = b[good] / np.sqrt(rss[good] / df / np.broadcast_to(nrm[:, None], b.shape)[good])
    return t, df

def lambda_gc(P):
    """Genomic control factor from a vector of p-values."""
    from scipy.stats import chi2
    P = np.asarray(P, float)
    P = P[np.isfinite(P)]
    if P.size == 0:
        return float("nan")
    x = chi2.isf(np.clip(P, 1e-300, 1), 1)
    return float(np.median(x) / 0.4549364)


def _adj_stats(H):
    """Vectorised LD between each marker and its right neighbour.

    Returns (|D'|, r2), each of length m-1, computed on the haplotype matrix
    itself so no genotype calling or phasing assumptions are introduced here.
    D' is returned as its absolute value: the sign of D depends only on which
    allele is coded 1, so strong LD in either phase counts as strong LD.
    Both statistics are formed from integer haplotype counts (N*n11 - nA*nB and
    integer D_max), so recoding the alleles of any marker or reordering the
    samples gives bit-identical values and therefore identical block boundaries
    (including which of several equal |D'| values is the weakest link).
    """
    H = np.asarray(H)
    N = np.int64(H.shape[0])
    A = H[:, :-1].astype(np.int64)
    B = H[:, 1:].astype(np.int64)
    nA = A.sum(0)
    nB = B.sum(0)
    n11 = (A * B).sum(0)
    Dn = N * n11 - nA * nB                       # N^2 * D, exact
    up = np.minimum(nA * (N - nB), (N - nA) * nB)
    dn = np.minimum(nA * nB, (N - nA) * (N - nB))
    dmax = np.where(Dn > 0, up, dn)
    aD = np.abs(Dn).astype(np.float64)
    with np.errstate(divide="ignore", invalid="ignore"):
        dp = np.where(dmax > 0, aD / np.maximum(dmax, 1).astype(np.float64), 0.0)
        va = (nA * (N - nA)).astype(np.float64)
        vb = (nB * (N - nB)).astype(np.float64)
        den = va * vb
        r2 = np.where(den > 0, (aD * aD) / np.where(den > 0, den, 1.0), 0.0)
    return np.clip(dp, 0.0, 1.0), np.clip(r2, 0.0, 1.0)

def _hap_cover(Hb, minf):
    packed = np.packbits(np.ascontiguousarray(Hb, dtype=np.uint8), axis=1)
    _, cnt = np.unique(packed, axis=0, return_counts=True)
    f = cnt / Hb.shape[0]
    # integer sum of counts: independent of the (coding-dependent) order of patterns
    return float(int(cnt[f >= minf].sum()) / Hb.shape[0])


def ld_blocks(H, pos, chrom=None, dmin=0.8, r2min=0.0, max_span=200000,
              min_mk=2, hap_cover=0.8, minf=0.05):
    """Haplotype blocks cut where linkage disequilibrium collapses.

    A boundary is placed between markers j and j+1 when |D'| (or r^2, if r2min is
    given) falls below the threshold, when the gap to the next marker exceeds
    max_span, or at a chromosome change. The block is therefore a unit of
    measured LD rather than a unit of marker count, so a cut cannot land inside
    a stretch of markers that are in LD.

    Two further rules keep every block testable, fixed in advance and applied
    identically to every panel:
      * short-block merging never crosses a hard cut (chromosome change or a
        marker gap > max_span), so no block spans a chromosome boundary;
      * a block whose common haplotypes (frequency >= minf) cover less than
        hap_cover of the chromosomes is split at its weakest internal D', and
        the rule is re-applied to the pieces. On sparse arrays in long-LD
        populations D' alone produces long blocks in which every haplotype is
        rare, leaving the haplotype tests with nothing to test.

    Returns a dict with blocks (list of (start, end) marker indices), the
    adjacency dprime / r2 arrays, and counts of each cut reason.
    """
    H = np.asarray(H)
    m = H.shape[1]
    if m < 2:
        return dict(blocks=[(0, m)] if m else [], dprime=np.zeros(0),
                    r2=np.zeros(0), n_cut=0, n_span_cut=0, n_chrom_cut=0)
    pos = np.asarray(pos, np.int64)
    dp, r2v = _adj_stats(H)
    rule = (r2v < r2min) if r2min > 0 else (dp < dmin)

    span_cut = (pos[1:] - pos[:-1]) > max_span
    if chrom is None:
        chrom_cut = np.zeros(m - 1, bool)
    else:
        ch = np.asarray(chrom)
        chrom_cut = ch[1:] != ch[:-1]
    hard = span_cut | chrom_cut
    cut = rule | hard
    hard_after = set((np.where(hard)[0] + 1).tolist())

    bounds = [0] + [int(j) + 1 for j in np.where(cut)[0]] + [m]
    blocks = [(bounds[i], bounds[i + 1])
              for i in range(len(bounds) - 1) if bounds[i + 1] > bounds[i]]

    if min_mk > 1:
        merged = []
        for s, e in blocks:
            if merged and (e - s) < min_mk and s not in hard_after:
                merged[-1] = (merged[-1][0], e)
            else:
                merged.append((s, e))
        blocks = merged

    n_div = 0
    if hap_cover and hap_cover > 0:
        out, stack = [], list(reversed(blocks))
        while stack:
            s, e = stack.pop()
            if e - s >= 2 and _hap_cover(H[:, s:e], minf) < hap_cover:
                j = s + int(np.argmin(dp[s:e - 1]))          # weakest adjacent |D'| (leftmost if tied)
                stack.append((j + 1, e)); stack.append((s, j + 1)); n_div += 1
            else:
                out.append((s, e))
        blocks = out

    return dict(blocks=blocks, dprime=dp, r2=r2v,
                n_cut=int(rule.sum()), n_span_cut=int(span_cut.sum()),
                n_chrom_cut=int(chrom_cut.sum()), n_div_split=n_div)


def jitter_blocks(blocks, m, step=1):
    """Boundary-perturbed copies of a block set, for ACAT over cut placement.

    Returns the original partition, one with every block widened by `step`
    markers, and one with every block narrowed by `step`. Combining the three at
    equal weights makes the reported statistic insensitive to the exact cut
    site, which is the honest way to report a data-driven partition.
    """
    wide, narrow = [], []
    for s, e in blocks:
        wide.append((max(0, s - step), min(m, e + step)))
        ns, ne = min(m - 1, s + step), max(min(m, e - step), min(m, s + step))
        narrow.append((ns, max(ne, ns + 1)))
    return [list(blocks), wide, narrow]
