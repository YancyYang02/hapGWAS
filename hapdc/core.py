
"""hapdc.core -- mathematical core.

All routines operate on haplotype matrices in (n_haplotype, n_marker) layout,
where the two haplotypes of individual i occupy rows 2i and 2i+1.
"""
import numpy as np
from scipy.special import betainc
from scipy import optimize

__all__ = ["windows", "window_design", "reml_delta", "acat", "t2p", "f2p",
           "lmm_scan_additive", "lambda_gc", "ld_blocks", "jitter_blocks"]


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
    indicators, nk haplotype classes retained, ncell diplotype cells retained,
    and nrare (haplotype copies reassigned) / nsmall (individuals reassigned).

    Haplotypes below `minf` are reassigned to the nearest retained haplotype by
    Hamming distance (ties -> the more frequent class), and individuals in
    diplotype cells below `minc` to the most frequent retained cell sharing one
    of their haplotypes (else the most frequent cell). Every individual is
    therefore in exactly one class pair and one cell: no implicit group with a
    free mean is created (leaving rare carriers out of all columns would, next
    to an intercept, give them their own parameter).
    Rows are compared as packed bytes, not as integers.
    """
    H = np.ascontiguousarray(H, dtype=np.uint8)
    nh, nm = H.shape
    n = nh // 2
    if nm == 0:
        return dict(dos=np.zeros((n, 0)), C=np.zeros((n, 0)), nk=0, ncell=0, nrare=0, nsmall=0)
    packed = np.ascontiguousarray(np.packbits(H, axis=1))
    # rows compared as raw bytes through a void view: same lexicographic order as
    # np.unique(axis=0) on uint8 rows, but a 1-D sort (exact)
    vv = packed.view(np.dtype((np.void, packed.shape[1]))).ravel()
    _, first, inv, cnt = np.unique(vv, return_index=True, return_inverse=True, return_counts=True)
    inv = np.asarray(inv).ravel()
    upat = packed[first]
    keep = np.where(cnt / nh >= minf)[0]
    keep = keep[np.argsort(-cnt[keep], kind="stable")]
    nk = len(keep)
    if nk == 0:
        return dict(dos=np.zeros((n, 0)), C=np.zeros((n, 0)), nk=0, ncell=0, nrare=0, nsmall=0)
    lab = np.full(len(upat), -1, np.int64)
    lab[keep] = np.arange(nk)
    rare = np.where(lab < 0)[0]
    if len(rare):
        Kh = H[first[keep]].astype(np.float64)
        Hr = H[first[rare]].astype(np.float64)
        # Hamming distances rare x kept in one product; argmin returns the first
        # minimum, and keep is frequency-sorted: ties -> more frequent class
        dist = Hr @ (1.0 - Kh).T + (1.0 - Hr) @ Kh.T
        lab[rare] = np.argmin(np.rint(dist), axis=1)
    cls = lab[inv]
    nrare = int(cnt[rare].sum()) if len(rare) else 0
    h1, h2 = cls[0::2], cls[1::2]
    dos = np.zeros((n, nk))
    for k in range(nk):
        dos[:, k] = (h1 == k).astype(float) + (h2 == k).astype(float)
    lo = np.minimum(h1, h2)
    hi = np.maximum(h1, h2)
    cell = lo.astype(np.int64) * nk + hi
    cu, cc = np.unique(cell, return_counts=True)
    order = np.argsort(-cc, kind="stable")
    cu, cc = cu[order], cc[order]
    big = cu[cc >= minc]
    nsmall = 0
    if len(big) == 0:
        return dict(dos=dos, C=np.zeros((n, 0)), nk=nk, ncell=0, nrare=nrare, nsmall=0)
    bl, bh = big // nk, big % nk
    small = cu[cc < minc]
    for c in small:
        m = cell == c
        a, b = c // nk, c % nk
        share = np.flatnonzero((bl == a) | (bh == a) | (bl == b) | (bh == b))
        cell[m] = big[share[0]] if len(share) else big[0]
        nsmall += int(m.sum())
    C = np.stack([(cell == c).astype(float) for c in big], 1)
    return dict(dos=dos, C=C, nk=nk, ncell=len(big), nrare=nrare, nsmall=nsmall, h1=h1, h2=h2)

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
    (n, q) matrix (intercept plus covariates)."""
    n = len(yr)
    Xr = np.asarray(Xr, float)
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


def lmm_scan_additive(Gt, Yt, l0, dlt):
    """Vectorised genome-wide additive scan on the eigen scale.

    Gt : (n, m) already eigen-rotated dosages (U.T @ dosage).
    Yt : (n, B) already rotated phenotypes, column 0 = observed.
    """
    n = Gt.shape[0]
    w = 1.0 / np.sqrt(l0 + dlt)
    X = Gt * w[:, None]
    Yw = Yt * w[:, None]
    q0 = w / np.linalg.norm(w)
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

    Returns (dprime, r2), each of length m-1, computed on the haplotype matrix
    itself so no genotype calling or phasing assumptions are introduced here.
    """
    H = np.asarray(H)
    A = H[:, :-1].astype(np.float64)
    B = H[:, 1:].astype(np.float64)
    pA = A.mean(0)
    pB = B.mean(0)
    p11 = (A * B).mean(0)
    D = p11 - pA * pB
    up = np.minimum(pA * (1.0 - pB), (1.0 - pA) * pB)
    dn = np.minimum(pA * pB, (1.0 - pA) * (1.0 - pB))
    dmax = np.where(D > 0, up, dn)
    with np.errstate(divide="ignore", invalid="ignore"):
        dp = np.where(dmax > 1e-12, D / dmax, 0.0)
        den = pA * (1.0 - pA) * pB * (1.0 - pB)
        r2 = np.where(den > 1e-12, (D * D) / den, 0.0)
    return np.clip(dp, -1.0, 1.0), np.clip(r2, 0.0, 1.0)


def _hap_cover(Hb, minf):
    packed = np.packbits(np.ascontiguousarray(Hb, dtype=np.uint8), axis=1)
    _, cnt = np.unique(packed, axis=0, return_counts=True)
    f = cnt / Hb.shape[0]
    return float(f[f >= minf].sum())


def ld_blocks(H, pos, chrom=None, dmin=0.8, r2min=0.0, max_span=200000,
              min_mk=2, hap_cover=0.8, minf=0.05):
    """Haplotype blocks cut where linkage disequilibrium collapses.

    A boundary is placed between markers j and j+1 when D' (or r^2, if r2min is
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
                j = s + int(np.argmin(dp[s:e - 1]))          # weakest adjacent D'
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
