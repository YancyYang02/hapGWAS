"""hapdc.nadda -- anchor-wide additive-closure non-additivity test (addendum NF2).

For each haplotype unit inside an anchor [a, e) the reduced model holds the fixed
effects, the unit's common-haplotype dosages and the single-marker dosages of ALL
markers of the anchor; the full model adds the unit's diplotype-cell indicators.
The cell F test is ACAT-combined over the anchor's units with the same
scale weights as P_NADD (ms.wacat)."""
import numpy as np
from scipy.linalg import qr as _qr
from .core import window_design, f2p
from .ms import anchor_units, wacat


def _extend(Q, X, tol=1e-8):
    """Orthonormal basis of span([Q, X]) as [Q, Qx]; rank-pruned by pivoted QR."""
    if X.shape[1] == 0:
        return Q
    nx = np.sqrt((X * X).sum(0)).max()
    R = X - Q @ (Q.T @ X)
    R = R - Q @ (Q.T @ R)
    Qx, Rx, piv = _qr(R, mode="economic", pivoting=True)
    d = np.abs(np.diag(Rx))
    r = int((d > tol * max(nx, 1e-30)).sum())
    return np.column_stack([Q, Qx[:, :r]])


def nadd_anchor(sc, anc_i, YS=None):
    """P_NADDA (B,) for one anchor anc_i = (a, e, lds) of a hapdc Scan `sc`."""
    a, e, lds = anc_i
    if YS is None:
        YS = sc._design_matrix(None, 0)
    Ut = sc.U.T
    n = sc.n
    Q0 = np.linalg.qr(sc.one if sc.one.ndim == 2 else sc.one[:, None])[0]
    Xg = (Ut @ sc.G[:, a:e].astype(np.float64)) * sc.w[:, None]
    Qg = _extend(Q0, Xg)
    pl, sl = [], []
    for scl, s, t in anchor_units(a, e, lds):
        d = window_design(sc.H[:, s:t], sc.minf, sc.minc)
        if d["nk"] < 2 or d["ncell"] < 2:
            continue
        Xd = (Ut @ d["dos"]) * sc.w[:, None]
        Xc = (Ut @ d["C"]) * sc.w[:, None]
        Qa = _extend(Qg, Xd)
        Qb = _extend(Qa, Xc)
        ka, kb = Qa.shape[1], Qb.shape[1]
        if kb <= ka or n - kb <= 0:
            continue
        Ra = YS - Qa @ (Qa.T @ YS)
        Rb = YS - Qb @ (Qb.T @ YS)
        rsa = (Ra * Ra).sum(0); rsb = (Rb * Rb).sum(0)
        F = ((rsa - rsb) / (kb - ka)) / (rsb / (n - kb))
        pl.append(np.clip(f2p(F, kb - ka, n - kb), 1e-300, 1.0)); sl.append(scl)
    if not pl:
        return np.ones(YS.shape[1])
    return wacat(pl, sl)[0]
