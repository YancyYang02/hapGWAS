"""hapdc.ms -- multi-scale haplotype scan.

Existing haplotype GWAS tools ask the user to choose one unit in advance: a
fixed window width (GHap, RAINBOWR) or an LD-block partition. The resolution
that carries the signal depends on the unknown genetic architecture -- a single
causal SNP is best seen at the single-marker scale, an allelic series or a
recessive haplotype at the scale of the haplotype that carries it.

The multi-scale scan removes that choice. The genome is cut into anchors
(consecutive LD blocks merged to >= `min_mk` markers). Inside each anchor every
scale of a prespecified ladder is tested with the same mixed model and the same
variance ratio -- single markers, sliding windows of 2/4/8/16 markers, the LD
blocks and the anchor itself -- and all unit p-values are combined by the
Cauchy combination (ACAT), which stays valid under arbitrary dependence among
the units. Each scale carries total weight 1, split equally among its units, so
no scale dominates merely because it has many units. Weights are fixed before
any data are analysed.

The anchor is the test unit: the genome-wide threshold is 0.05 / n_anchors.
The scale contributing most to the anchor statistic is reported as the
resolution at which the signal is carried, and the same combination applied to
the non-additivity statistic NADD_F gives an anchor-level test of the mode of
inheritance.
"""
import numpy as np

from .core import _T, acat, t2p

SCALES = ("SNP", "W2", "W4", "W8", "W16", "LD", "ANCHOR")
WSPEC = {"W2": (2, 1), "W4": (4, 2), "W8": (8, 4), "W16": (16, 8)}

__all__ = ["SCALES", "anchors", "anchor_units", "wacat", "ms_scan"]


def anchors(blocks, min_mk=16, chrom=None):
    """Merge consecutive LD blocks into anchors of >= min_mk markers.

    blocks : sorted list of (s, e) half-open marker intervals covering the markers.
    Returns a list of (a, e, [ld blocks inside]). LD blocks are never cut; an
    anchor never crosses a chromosome change; a short remainder at the end of a
    chromosome is merged into the previous anchor of the same chromosome.
    """
    out, cur = [], None

    def same(i, j):
        return chrom is None or chrom[i] == chrom[j]

    for s, e in blocks:
        if cur is not None and not same(cur[0], s):
            if out and same(out[-1][0], cur[0]) and out[-1][1] == cur[0]:
                out[-1] = (out[-1][0], cur[1], out[-1][2] + cur[2])
            else:
                out.append(tuple(cur))
            cur = None
        if cur is None:
            cur = [s, e, [(s, e)]]
        else:
            cur[1] = e
            cur[2].append((s, e))
        if cur[1] - cur[0] >= min_mk:
            out.append(tuple(cur))
            cur = None
    if cur is not None:
        if out and same(out[-1][0], cur[0]) and out[-1][1] == cur[0]:
            out[-1] = (out[-1][0], cur[1], out[-1][2] + cur[2])
        else:
            out.append(tuple(cur))
    return out


def anchor_units(a, e, lds):
    """All haplotype units of one anchor as (scale, s, e)."""
    units = []
    L = e - a
    for sc, (w, st) in WSPEC.items():
        if L < w:
            continue
        starts = list(range(a, e - w + 1, st))
        if starts[-1] + w < e:
            starts.append(e - w)
        units += [(sc, s, s + w) for s in starts]
    units += [("LD", s, t) for s, t in lds]
    units.append(("ANCHOR", a, e))
    return units


def wacat(plist, slist):
    """Scale-weighted ACAT.

    plist : list of (B,) p-value arrays, one per unit; slist : scale label per unit.
    Every scale present gets total weight 1, split equally among its units.
    Callers pass testable units only (ms_scan drops units whose design admits
    no test, using genotype-only validity flags). A p-value equal to 1 (to
    within 1e-12) from a testable unit is replaced by 1 - 1/d (d = number of
    p-values combined): the Cauchy transform of p = 1 is -infinity, so one such
    unit would otherwise veto the rest. All other p-values enter the standard
    Cauchy combination unchanged.
    Returns combined p (B,) and the per-scale contribution matrix (nscale, B).
    """
    P = np.vstack(plist)
    d = P.shape[0]
    P = np.where(P >= 1.0 - 1e-12, 1.0 - 1.0 / max(d, 2), P)
    P = np.clip(P, 1e-300, 1.0)
    slist = np.asarray(slist)
    present = [s for s in SCALES if (slist == s).any()]
    w = np.zeros(d)
    for s in present:
        k = slist == s
        w[k] = 1.0 / k.sum()
    p = acat(P, w)
    T = _T(P) * (w / w.sum())[:, None]
    contrib = np.vstack([T[slist == s].sum(0) for s in present])
    return p, contrib, present


def _snp_p(sc, idx):
    """Single-marker LMM p-values (len(idx), B) under the scan's variance ratio,
    and a (len(idx),) flag: marker testable (non-constant after the fixed effects)."""
    YS = sc._design_matrix(None, 0)
    R0 = sc._p0(YS)
    rss0 = (R0 * R0).sum(0)
    n = sc.n
    df = n - sc.nfix - 1
    out = np.ones((len(idx), YS.shape[1]))
    valid = np.zeros(len(idx), bool)
    for c0 in range(0, len(idx), 1000):
        jj = idx[c0:c0 + 1000]
        # project only the marker range this chunk touches (not the whole
        # chromosome)
        j0, j1 = int(jj.min()), int(jj.max()) + 1
        Xa = sc._p0(sc._projG(j0, j1) * sc.w[:, None])
        ok = sc.testable(j0, j1, Xa)[jj - j0]
        X = Xa[:, jj - j0]
        nr = (X * X).sum(0)
        if ok.any():
            b = (X[:, ok].T @ R0) / nr[ok][:, None]
            rss = np.clip(rss0[None, :] - b * b * nr[ok][:, None], 1e-300, None)
            out[c0 + np.where(ok)[0]] = t2p(b / np.sqrt(rss / df / nr[ok][:, None]), df)
        valid[c0:c0 + len(jj)] = ok
    return out, valid


def ms_scan(sc, anc, which=None, variants=("O", "D"), nadd=True, scales=False, nadds=True):
    """Multi-scale scan of the anchors `anc[i]` for i in `which` (default all).

    sc : a hapdc.scan.Scan (or subclass with several phenotype columns).
    Unit statistic per variant: "O" = OMNI (A_joint, CELL_F, SACAT);
    "D" = equal-weight ACAT of A_joint, CELL_F, SACAT and NADD_F.
    Returns dict with P_MS[v] (na, B), P_NADD (na, B), best_scale[v] (na,) for
    column 0, best_unit[v] (na, 2) of the smallest-p unit (column 0), unit counts.
    """
    which = list(range(len(anc))) if which is None else list(which)
    # unique haplotype intervals across all requested anchors
    plan, uniq = [], {}
    for i in which:
        a, e, lds = anc[i]
        us = anchor_units(a, e, lds)
        plan.append(us)
        for _, s, t in us:
            if (s, t) not in uniq:
                uniq[(s, t)] = len(uniq)
    ivals = sorted(uniq, key=uniq.get)
    need = {"A_joint", "SACAT", "CELL_F", "NADD_F", "OMNI"}
    if "S" in variants:
        need.add("P1_acat")
    if nadd and nadds:
        need.add("NADD_FS")
    res = sc.run(blocks=ivals, verbose=False, need=need) if ivals else None
    snp_idx = np.concatenate([np.arange(anc[i][0], anc[i][1]) for i in which]) if which else np.zeros(0, int)
    pS, vS = _snp_p(sc, snp_idx)
    spos = {int(j): r for r, j in enumerate(snp_idx)}
    B = pS.shape[1]
    na = len(which)
    P = res["P"] if res is not None else {}
    V = res["V"] if res is not None else {}
    # Components are selected with genotype-only validity flags V (never with
    # the observed p-values), so a phenotype gets the same combined p whether it
    # is scanned alone or in a batch with other phenotype columns.
    unitp, unitv = {}, {}
    if "O" in variants:
        unitp["O"] = P.get("OMNI"); unitv["O"] = V.get("OMNI")
    if "D" in variants and res is not None:
        D = np.ones_like(P["OMNI"]); DV = np.zeros(D.shape[0], bool)
        for k in range(D.shape[0]):
            parts = [P[q][k] for q in ("A_joint", "CELL_F", "SACAT", "NADD_F") if V[q][k]]
            if parts:
                D[k] = acat(np.clip(np.vstack(parts), 1e-300, 1 - 1e-12), np.ones(len(parts)))
                DV[k] = True
        unitp["D"] = D; unitv["D"] = DV
    if "S" in variants and res is not None:
        S5 = np.ones_like(P["OMNI"]); SV = np.zeros(S5.shape[0], bool)
        for k in range(S5.shape[0]):
            parts = [P[q][k] for q in ("A_joint", "CELL_F", "SACAT", "P1_acat", "NADD_F") if V[q][k]]
            if parts:
                S5[k] = acat(np.clip(np.vstack(parts), 1e-300, 1 - 1e-12), np.ones(len(parts)))
                SV[k] = True
        unitp["S"] = S5; unitv["S"] = SV
    PMS = {v: np.ones((na, B)) for v in variants}
    best = {v: np.empty(na, dtype=object) for v in variants}
    bunit = {v: np.zeros((na, 2), int) for v in variants}
    PNA = np.ones((na, B))
    PNS = np.ones((na, B))
    PSC = {v: {q: np.ones((na, B)) for q in SCALES} for v in variants} if scales else None
    nunits = np.zeros(na, int)
    for r, (i, us) in enumerate(zip(which, plan)):
        a, e, _ = anc[i]
        sl = ["SNP"] * (e - a) + [u[0] for u in us]
        iv = [(j, j + 1) for j in range(a, e)] + [(s, t) for _, s, t in us]
        ks = [uniq[(s, t)] for _, s, t in us]
        snp = [pS[spos[j]] for j in range(a, e)]
        snpv = [bool(vS[spos[j]]) for j in range(a, e)]
        for v in variants:
            pl_all = snp + [unitp[v][k] for k in ks]
            ok_u = snpv + [bool(unitv[v][k]) for k in ks]
            sel = [u for u, o in enumerate(ok_u) if o]
            if not sel:
                best[v][r] = "none"
                bunit[v][r] = iv[0]
                continue
            pl = [pl_all[u] for u in sel]
            slv = [sl[u] for u in sel]
            p, contrib, present = wacat(pl, slv)
            PMS[v][r] = p
            best[v][r] = present[int(np.argmax(contrib[:, 0]))]
            bunit[v][r] = iv[sel[int(np.argmin([x[0] for x in pl]))]]
            if scales:
                Pv = np.vstack(pl)
                Pv = np.clip(np.where(Pv >= 1.0 - 1e-12, 1.0 - 1.0 / max(len(pl), 2), Pv), 1e-300, 1.0)
                sla = np.asarray(slv)
                for q in present:
                    kq = sla == q
                    PSC[v][q][r] = acat(Pv[kq], np.ones(int(kq.sum())))
        nunits[r] = len(sl)
        if nadd:
            for key, dst in (("NADD_F", PNA), ("NADD_FS", PNS)):
                if key == "NADD_FS" and not nadds:
                    continue
                pn, sn = [], []
                for (scl, _, _), k in zip(us, ks):
                    if V[key][k]:
                        pn.append(P[key][k]); sn.append(scl)
                if pn:
                    dst[r] = wacat(pn, sn)[0]

    return dict(which=which, P_MS=PMS, P_NADD=PNA, best_scale=best, best_unit=bunit,
                nunits=nunits, n_unique=len(ivals), P_SC=PSC, P_NADDS=PNS)
