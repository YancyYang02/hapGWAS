
"""hapdc.cli -- command line entry points: prep, scan, report, selftest."""
import argparse
import json
import os
import sys

import numpy as np

from . import core, io, report, scan


def _mkout(p):
    os.makedirs(p, exist_ok=True)
    return p


# ---------------------------------------------------------------- prep
def cmd_prep(a):
    _mkout(a.out)
    ids, H, meta, dropped = io.read_vcf(a.vcf, a.nmax)
    if a.grm == "auto":
        # VanRaden GRM from every --grm-thin-th marker (MAF > 0.01)
        G = (H[0::2, ::a.grm_thin] + H[1::2, ::a.grm_thin]).astype(np.float32)
        pf = G.mean(0) / 2; ok_ = (pf > 0.01) & (pf < 0.99)
        Z = (G[:, ok_] - 2 * pf[ok_]) / np.sqrt(2 * pf[ok_] * (1 - pf[ok_]))
        gids, K = list(ids), (Z @ Z.T / Z.shape[1]).astype(float)
        print("GRM computed from %d markers" % Z.shape[1], flush=True)
    else:
        gids, K = io.read_grm(a.grm)
    vids = list(ids)
    keep_ids = [g for g in gids if g in vids]
    order = [vids.index(g) for g in keep_ids]
    m = H.shape[1]
    # H is (n_hap, m); rows 2j, 2j+1 are the two haplotypes of sample j
    H = np.ascontiguousarray(H.reshape(len(vids), 2, m)[order, :, :].reshape(2 * len(order), m))
    K = K[np.ix_([gids.index(g) for g in keep_ids], [gids.index(g) for g in keep_ids])]
    n = len(keep_ids)
    S, U = np.linalg.eigh(K)
    S = np.clip(S, 1e-6, None)
    np.savez(os.path.join(a.out, "eig.npz"), S=S, U=U)
    np.save(os.path.join(a.out, "K_keep.npy"), K.astype(np.float32))
    np.save(os.path.join(a.out, "hap.npy"), H.astype(np.int8))
    miss_codes = tuple(float(x) for x in str(a.miss).split(","))
    Y = []
    names = []
    for p in a.pheno.split(","):
        y = io.read_pheno(p, keep_ids, miss=miss_codes)
        names.append(os.path.splitext(os.path.basename(p))[0])
        Y.append(y)
    Y = np.column_stack(Y)
    ok = np.isfinite(Y).any(1)
    np.save(os.path.join(a.out, "traits.npy"), Y[ok])
    np.save(os.path.join(a.out, "keep.npy"), np.where(ok)[0])
    json.dump(dict(ids=[keep_ids[i] for i in np.where(ok)[0]], traits=names,
                   n=int(ok.sum()), m=int(m), dropped=dropped,
                   vcf=os.path.abspath(a.vcf), grm=os.path.abspath(a.grm)),
              open(os.path.join(a.out, "meta.json"), "w"), indent=1)
    json.dump([[list(x) for x in meta]], open(os.path.join(a.out, "markers.json"), "w"))
    print("prep done: n=%d m=%d traits=%s" % (ok.sum(), m, names))
    return 0


# ---------------------------------------------------------------- scan
def cmd_scan(a):
    out = _mkout(a.out)
    H = np.load(os.path.join(a.data, "hap.npy"))
    Y = np.load(os.path.join(a.data, "traits.npy"))
    K = np.load(os.path.join(a.data, "K_keep.npy"))
    meta = json.load(open(os.path.join(a.data, "markers.json")))[0]
    meta = [tuple(x) for x in meta]
    widths = [int(x) for x in str(a.wsize).split(",")]
    ti = a.trait
    y = np.asarray(Y[:, ti], float)
    ok = np.isfinite(y)
    if not ok.all():
        # Missing individuals are dropped for this trait only. H holds two
        # rows per individual, so the reduced GRM is re-decomposed rather
        # than subsetting the stored eigenvectors.
        m = H.shape[1]
        H = np.ascontiguousarray(H.reshape(len(y), 2, m)[ok].reshape(-1, m))
        Kr = K[np.ix_(ok, ok)].astype(float)
        S, U = np.linalg.eigh(Kr)
        S = np.clip(S, 1e-6, None)
        y = y[ok]
        print("trait %d: %d of %d kept after dropping missing"
              % (ti, int(ok.sum()), len(ok)), flush=True)
    else:
        E = np.load(os.path.join(a.data, "eig.npz"))
        S, U = E["S"], E["U"]
    sc = scan.Scan(H, meta, U, S, y, minf=a.minf, minc=a.minc, seed=a.seed)
    if a.nperm:
        sc.permute(a.nperm)
    results = []
    if a.units == "ld":
        # Default: LD-defined haplotype blocks (D-prime adjacency rule),
        # cut at chromosome changes and at gaps larger than --ld-maxspan.
        pos = np.array([int(x[2]) for x in meta]); chrom = [x[0] for x in meta]
        bl = core.ld_blocks(H, pos, chrom=chrom, dmin=a.ld_dmin, max_span=a.ld_maxspan, min_mk=2, hap_cover=a.ld_hapcover, minf=a.minf)
        units = [(int(s_), int(e_)) for s_, e_ in bl["blocks"]]
        if a.lo is not None or a.hi is not None:
            units = units[(a.lo or 0):(a.hi if a.hi is not None else len(units))]
        print("LD blocks: %d (D-prime<%.2f cuts %d, span cuts %d, chrom cuts %d)"
              % (len(units), a.ld_dmin, bl["n_cut"], bl["n_span_cut"], bl["n_chrom_cut"]), flush=True)
        r = sc.run(blocks=units, verbose=not a.quiet); r["w_size"] = "ld"
        results.append(r)
        np.savez_compressed(os.path.join(out, "scan_ld.npz"),
                            **{"P_" + k: r["P"][k] for k in r["P"]},
                            nk=r["nk"], ncell=r["ncell"], starts=np.array(r["starts"]),
                            ends=np.array([e_ for _, e_ in units]), B=r["B"])
        print("LD units done: %d units %.0fs" % (r["nwin"], r["elapsed"]), flush=True)
    else:
      for ws in widths:
        r = sc.run(ws, a.wstep, lo=a.lo, hi=a.hi, verbose=not a.quiet)
        results.append(r)
        np.savez_compressed(os.path.join(out, "scan_w%d.npz" % ws),
                            **{"P_" + k: r["P"][k] for k in r["P"]},
                            nk=r["nk"], ncell=r["ncell"], starts=np.array(r["starts"]),
                            w_size=ws, w_step=a.wstep, B=r["B"])
        print("width %d done: %d windows %.0fs" % (ws, r["nwin"], r["elapsed"]), flush=True)
    if a.nperm:
        rows = []
        for r in results:
            cal = sc.calibrate(r)
            for k, v in cal.items():
                rows.append(dict(trait=a.trait_name or str(ti), w_size=r["w_size"], stat=k,
                                 fpr_null=v["fpr"], thr_scan=v["thr"],
                                 n_below=v["n_below"], n_perm=v.get("n_perm", 0),
                                 delta=r["delta"], n=int(sc.n), nwin=int(r["nwin"])))
        report.write_summary(os.path.join(out, "calibration.csv"), rows)
        for row in rows:
            print("  %-9s w=%s fpr_null=%.3f thr=%.4g n_below=%d"
                  % (row["stat"], row["w_size"], row["fpr_null"], row["thr_scan"], row["n_below"]))
    report.write_windows_table(os.path.join(out, "windows_%s.csv" % (results[0]["w_size"] if a.units == "ld" else "w%d" % results[0]["w_size"])),
                               results[0], str(ti))
    json.dump(dict(trait=ti, units=a.units, ld_dmin=a.ld_dmin, ld_maxspan=a.ld_maxspan, ld_hapcover=a.ld_hapcover, widths=widths, wstep=a.wstep, minf=a.minf, minc=a.minc,
                   nperm=a.nperm, n=int(sc.n), delta=[r["delta"] for r in results]),
              open(os.path.join(out, "scan_meta.json"), "w"), indent=1)
    print("scan done ->", out)
    return 0


# ---------------------------------------------------------------- ms (primary analysis)
def _read_covar(path, ids):
    """Tab-separated table, header line, first column IID, remaining columns numeric."""
    tab = {}
    with open(path) as f:
        hdr = f.readline().rstrip("\n").split("\t")
        for ln in f:
            p = ln.rstrip("\n").split("\t")
            try:
                tab[p[0]] = [float(x) for x in p[1:len(hdr)]]
            except ValueError:
                tab[p[0]] = [np.nan] * (len(hdr) - 1)
    k = len(hdr) - 1
    return np.array([tab.get(i, [np.nan] * k) for i in ids], float), hdr[1:]


_MS_CTX = {}
_BLAS_LIMIT = None


def _limit_blas(n):
    try:
        from threadpoolctl import threadpool_limits
        return threadpool_limits(max(1, int(n)))
    except Exception:
        return None


def _ms_chrom(ch):
    """Scan one chromosome; reads shared inputs from _MS_CTX (fork-inherited)."""
    from . import ms as M
    g = _MS_CTX; a = g["a"]
    lim = None
    if g.get("blas"):
        try:
            from threadpoolctl import threadpool_limits
            lim = threadpool_limits(g["blas"])
        except Exception:
            lim = None
    cols = np.where(g["chrom"] == ch)[0]
    # read only this chromosome's columns from the memory-mapped haplotypes
    blk = np.asarray(g["Hall"][:, cols[0]:cols[-1] + 1])
    H = np.ascontiguousarray(blk[g["hrow"]][:, cols - cols[0]])
    p = g["pos"][cols]
    bl = core.ld_blocks(H, p, dmin=a.ld_dmin, max_span=a.ld_maxspan, min_mk=2)["blocks"]
    cov = np.zeros(len(cols), bool); units = []
    for s_, e_ in bl:
        units.append((int(s_), int(e_))); cov[s_:e_] = True
    units += [(int(j), int(j) + 1) for j in np.where(~cov)[0]]; units.sort()
    anc = M.anchors(units, a.anchor_min)
    meta = [g["markers"][j] for j in cols]
    sc = scan.Scan(H, meta, g["U"], g["S"], g["yt"], minf=a.minf, minc=a.minc, covar=g["C"])
    rows = []
    for lo in range(0, len(anc), a.chunk):
        r = M.ms_scan(sc, anc, which=range(lo, min(lo + a.chunk, len(anc))), variants=("D",), nadd=True, nadds=False)
        for k, i in enumerate(r["which"]):
            a0, e0, _ = anc[i]; bu = r["best_unit"]["D"][k]
            rows.append((ch, int(p[a0]), int(p[e0 - 1]), e0 - a0, float(r["P_MS"]["D"][k, 0]),
                         str(r["best_scale"]["D"][k]), int(p[bu[0]]), int(p[max(bu[1] - 1, 0)]),
                         float(r["P_NADD"][k, 0])))
    return ch, rows, len(cols), len(anc)


def cmd_ms(a):
    """Multi-scale anchor scan with unit statistic D; Bonferroni over anchors."""
    from . import ms as M
    import time
    out = _mkout(a.out)
    meta_j = json.load(open(os.path.join(a.data, "meta.json")))
    ids = meta_j["ids"]
    Hall = np.load(os.path.join(a.data, "hap.npy"), mmap_mode="r")
    Y = np.load(os.path.join(a.data, "traits.npy"))
    K = np.load(os.path.join(a.data, "K_keep.npy")).astype(float)
    keep_rows = np.load(os.path.join(a.data, "keep.npy"))
    markers = [tuple(x) for x in json.load(open(os.path.join(a.data, "markers.json")))[0]]
    y = np.asarray(Y[:, a.trait], float)
    C = None
    if a.covar:
        C, cnames = _read_covar(a.covar, ids)
        print("covariates:", cnames, flush=True)
    ok = np.isfinite(y) & (np.isfinite(C).all(1) if C is not None else True)
    idx = np.where(ok)[0]
    K = K[np.ix_(keep_rows, keep_rows)][np.ix_(idx, idx)]
    S, U = np.linalg.eigh(K)
    S = np.clip(S, 1e-6, None)
    hrow = np.sort(np.concatenate([2 * keep_rows[idx], 2 * keep_rows[idx] + 1]))
    chrom = np.array([str(x[0]) for x in markers]); pos = np.array([int(x[2]) for x in markers])
    chroms = [str(c) for c in dict.fromkeys(chrom.tolist()) if a.chrom is None or c in a.chrom.split(",")]
    tname = a.trait_name or meta_j["traits"][a.trait]
    print("trait %s: n=%d of %d, chromosomes %s" % (tname, len(idx), len(ids), chroms), flush=True)
    yt = y[idx].copy()
    if a.transform == "int" and len(np.unique(yt)) > 2:
        # Rank-based inverse normal transformation (Blom offsets). The
        # diplotype-cell F test compares means of small cells, and a single
        # extreme phenotype in a small cell inflates false positives on the
        # raw scale; INT is therefore the default for quantitative traits.
        from scipy.stats import rankdata, norm
        yt = norm.ppf((rankdata(yt) - 0.375) / (len(yt) + 0.25))
        print("phenotype: rank-based inverse normal transformation applied", flush=True)
    else:
        print("phenotype: analysed on the original scale", flush=True)
    # Kernel diagnostic: REML profile over spectral powers K^g. A large
    # likelihood-ratio statistic flags a relationship matrix that does not fit
    # the phenotypic covariance (e.g. near-Mendelian loci or breed structure);
    # the scan's anchor p-values are then not calibrated. With --kernel-power auto
    # the selected K^g is used (sensitivity analysis; it can absorb
    # population-differentiated loci into the random effect).
    Xk = np.ones(len(idx)) if C is None else np.column_stack([np.ones(len(idx)), C[idx]])
    kp_sel, kp_lrt, _tab = core.kernel_power(U.T @ yt, S, U.T @ Xk)
    kp_g = 1.0
    print("kernel diagnostic: REML-best g=%.2f, LRT vs g=1 = %.2f" % (max(_tab, key=lambda t: t[1])[0], kp_lrt), flush=True)
    if getattr(a, "kernel_power", "1") == "auto" and kp_sel != 1.0:
        kp_g = kp_sel
        S = S ** kp_g; S = S / S.mean()
    rows = []
    t0 = time.time()
    ctx = dict(Hall=Hall, hrow=hrow, chrom=chrom, pos=pos, markers=markers, U=U, S=S, yt=yt,
               C=C[idx] if C is not None else None, a=a)
    jobs = max(1, min(int(getattr(a, "jobs", 1) or 1), len(chroms)))
    if jobs == 1:
        _MS_CTX.update(ctx)
        for ch in chroms:
            ch_, rr, nm, na = _ms_chrom(ch)
            rows += rr
            print("chr %s: %d markers, %d anchors, %.0fs" % (ch_, nm, na, time.time() - t0), flush=True)
    else:
        # chromosomes in parallel (fork); 1 BLAS thread per worker unless --threads
        import multiprocessing as mp
        _MS_CTX.update(ctx)
        _MS_CTX["blas"] = a.threads if a.threads else 1
        order = sorted(chroms, key=lambda c: -int((chrom == c).sum()))
        res = {}
        with mp.get_context("fork").Pool(jobs) as pool:
            for ch_, rr, nm, na in pool.imap_unordered(_ms_chrom, order):
                res[ch_] = rr
                print("chr %s: %d markers, %d anchors, %.0fs" % (ch_, nm, na, time.time() - t0), flush=True)
        for ch in chroms:
            rows += res[ch]
    NA = len(rows); thr = 0.05 / NA
    _lam = float(core.lambda_gc(np.array([r_[4] for r_ in rows])))
    print("lambda_gc of anchor p-values = %.3f (null reference 0.61-0.68: the Cauchy combination is "
          "conservative away from the tail)" % _lam, flush=True)
    if _lam < 0.45:
        print("WARNING: lambda_gc = %.3f is well below the null reference; the single-GRM null model may fit "
              "poorly (kernel LRT = %.1f), as in samples made of discrete breeds or strata. Do not use the "
              "Bonferroni threshold; take the threshold from a permutation that preserves that structure."
              % (_lam, kp_lrt), flush=True)
    path = os.path.join(out, "%s.anchors.tsv" % tname)
    with open(path, "w") as f:
        f.write("chrom\tstart_bp\tend_bp\tn_markers\tP\tlead_scale\tlead_unit_start\tlead_unit_end\tP_nonadditive\tsig_bonf\n")
        for rw in rows:
            f.write("%s\t%d\t%d\t%d\t%.4g\t%s\t%d\t%d\t%.4g\t%d\n" % (rw + (int(rw[4] < thr),)))
    json.dump(dict(trait=tname, n=int(len(idx)), n_anchors=NA, threshold_bonf=thr, n_sig=int(sum(r_[4] < thr for r_ in rows)),
                   lambda_gc=float(core.lambda_gc(np.array([r_[4] for r_ in rows]))), anchor_min=a.anchor_min,
                   kernel_power=kp_g, kernel_power_lrt=kp_lrt, kernel_power_best=float(max(_tab, key=lambda t: t[1])[0]),
                   ld_dmin=a.ld_dmin, ld_maxspan=a.ld_maxspan, minf=a.minf, minc=a.minc, covar=a.covar,
                   elapsed_s=round(time.time() - t0, 1)),
              open(os.path.join(out, "%s.summary.json" % tname), "w"), indent=1)
    print("ms done: %d anchors, threshold %.3g, %d significant -> %s" % (NA, thr, sum(r_[4] < thr for r_ in rows), path))
    return 0


# ---------------------------------------------------------------- report
def cmd_report(a):
    out = _mkout(a.out)
    metas = []
    for d in a.scans:
        mp = os.path.join(d, "scan_meta.json")
        if not os.path.exists(mp):
            continue
        sm = json.load(open(mp))
        for f in sorted(os.listdir(d)):
            if f.startswith("scan_w") and f.endswith(".npz"):
                z = np.load(os.path.join(d, f))
                metas.append((d, f, z, sm))
    if not metas:
        print("nothing to report")
        return 1
    for d, f, z, sm in metas:
        ws = int(z["w_size"])
        for k in ("A_joint", "P1_sidak", "D_acat", "HDCb"):
            key = "P_" + k
            if key not in z:
                continue
            p = z[key][:, 0]
            report.write_manhattan_svg(os.path.join(out, "manhattan_%s_w%d.svg" % (k, ws)),
                                       p, z["starts"], None, label=k)
            report.write_qq_svg(os.path.join(out, "qq_%s_w%d.svg" % (k, ws)),
                                p[z["nk"] >= 2], label=k)
    print("report done ->", out)
    return 0


# ---------------------------------------------------------------- selftest
def cmd_selftest(a):
    """Known-architecture check on founder-haplotype blocks.

    Random 0/1 haplotypes are useless here: a 6-marker window yields up to 64
    patterns, none reaching MINF, so every window is skipped and all p-values
    are 1. Haplotypes are therefore drawn from a small founder set, which is what
    real phased data looks like inside a block.
    """
    rng = np.random.default_rng(7)
    n, nblk, blk = 300, 3, 10
    m = nblk * blk
    # founder haplotypes per block: 3 founders, distinguishable patterns
    H = np.zeros((2 * n, m), np.int8)
    for b in range(nblk):
        F = np.zeros((3, blk), np.int8)
        F[1] = rng.integers(0, 2, blk)
        F[2] = rng.integers(0, 2, blk)
        draws = rng.integers(0, 3, size=(2 * n,))
        sl = slice(b * blk, (b + 1) * blk)
        H[:, sl] = F[draws]
    S, U = np.linalg.eigh(np.eye(n))
    meta = [("1", "m%d" % i, i + 1, "A", "C") for i in range(m)]
    causal = slice(0, blk)

    def run_case(y, tag):
        y = (y - y.mean()) / y.std()
        sc = scan.Scan(H, meta, U, S, y)
        sc.permute(200)
        r = sc.run(4, 2, verbose=False)
        p1 = {k: float(r["P"][k][0, 0]) for k in ("A_joint", "P1_sidak", "D_acat", "HDCb")}
        print(tag, {k: float("%.3g" % v) for k, v in p1.items()})
        return p1

    # additive: dosage of founder-1 haplotype in the causal block
    dos = (H[:, causal] == 0).sum(1).astype(float)
    yadd = dos[0::2] + dos[1::2] + rng.normal(0, 1.0, n)
    print("=== window 0 sits on the causal block ===")
    run_case(yadd, "ADDITIVE   ")

    # overdominance: heterozygosity at the causal block
    h1, h2 = H[0::2, causal], H[1::2, causal]
    het = (np.abs(h1.astype(int) - h2.astype(int)).sum(1) > 0).astype(float)
    yhet = het + rng.normal(0, 1.0, n)
    run_case(yhet, "OVERDOM    ")

    # pure noise: both should be null
    run_case(rng.normal(0, 1.0, n), "NULL       ")
    return 0


def main(argv=None):
    p = argparse.ArgumentParser(prog="hapdc", description="Haplotype + diplotype GWAS with a permutation null")
    sub = p.add_subparsers(dest="cmd", required=True)

    q = sub.add_parser("prep", help="intersect VCF, GRM and phenotypes")
    q.add_argument("vcf"); q.add_argument("grm", help="GCTA .grm prefix, or 'auto' to compute a VanRaden GRM from the VCF"); q.add_argument("pheno", help="comma separated")
    q.add_argument("--grm-thin", type=int, default=10, help="with grm=auto: use every k-th marker")
    q.add_argument("out"); q.add_argument("--nmax", type=int, default=None)
    q.add_argument("--miss", default="-9", help="comma separated missing codes")
    q.set_defaults(func=cmd_prep)

    q = sub.add_parser("scan", help="window scan for one trait")
    q.add_argument("data"); q.add_argument("out"); q.add_argument("trait", type=int)
    q.add_argument("--units", choices=["ld", "fixed"], default="ld",
                   help="test units: LD-defined haplotype blocks (default) or fixed marker windows")
    q.add_argument("--ld-dmin", type=float, default=0.8, help="D-prime threshold for a block cut")
    q.add_argument("--ld-hapcover", type=float, default=0.8, help="split an LD block until haplotypes with frequency >= --minf cover this fraction of chromosomes")
    q.add_argument("--ld-maxspan", type=int, default=200000, help="cut blocks at gaps larger than this (bp)")
    q.add_argument("--wsize", default="4", help="--units fixed: comma separated ladder, e.g. 2,4,8")
    q.add_argument("--wstep", type=int, default=2)
    q.add_argument("--minf", type=float, default=0.05)
    q.add_argument("--minc", type=int, default=5)
    q.add_argument("--nperm", type=int, default=1000)
    q.add_argument("--seed", type=int, default=20261110)
    q.add_argument("--lo", type=int, default=None); q.add_argument("--hi", type=int, default=None)
    q.add_argument("--quiet", action="store_true")
    q.add_argument("--trait-name", default=None)
    q.add_argument("--threads", type=int, default=0, help="BLAS threads (default 1)")
    q.set_defaults(func=cmd_scan)

    q = sub.add_parser("ms", help="multi-scale anchor scan (primary analysis, statistic D)")
    q.add_argument("data", help="directory written by `hapdc prep`"); q.add_argument("out")
    q.add_argument("trait", type=int, help="0-based trait column from prep")
    q.add_argument("--covar", default=None, help="tab-separated covariate table: IID then numeric columns, header line")
    q.add_argument("--transform", choices=["int", "none"], default="int",
                   help="phenotype transformation: rank-based inverse normal (default; skipped for binary traits) or none")
    q.add_argument("--chrom", default=None, help="comma separated subset of chromosomes")
    q.add_argument("--anchor-min", type=int, default=16, help="minimum markers per anchor")
    q.add_argument("--ld-dmin", type=float, default=0.8); q.add_argument("--ld-maxspan", type=int, default=200000)
    q.add_argument("--minf", type=float, default=0.05); q.add_argument("--minc", type=int, default=5)
    q.add_argument("--chunk", type=int, default=500, help="anchors per batch (memory)")
    q.add_argument("--jobs", type=int, default=1, help="chromosomes scanned in parallel (POSIX fork); results identical to --jobs 1")
    q.add_argument("--threads", type=int, default=0, help="BLAS threads per worker (default 1; more threads give no speed-up for hapdc's small per-unit matrices)")
    q.add_argument("--trait-name", default=None)
    q.add_argument("--kernel-power", choices=["auto", "1"], default="1",
                   help="relationship kernel: 1 = standard GRM (default); auto = spectral power K^g selected by REML (LRT, 1 df, 0.05), sensitivity analysis")
    q.set_defaults(func=cmd_ms)

    q = sub.add_parser("report", help="tables and SVG figures")
    q.add_argument("scans", nargs="+"); q.add_argument("--out", default="report")
    q.set_defaults(func=cmd_report)

    q = sub.add_parser("selftest", help="known-architecture sanity check")
    q.set_defaults(func=cmd_selftest)

    a = p.parse_args(argv)
    # BLAS threading: the per-unit matrices are small, so multithreaded BLAS adds
    # spin-wait overhead without speed-up (and oversubscribes Slurm allocations,
    # where os.cpu_count() reports the whole node). Default 1 thread; --threads overrides.
    global _BLAS_LIMIT
    _BLAS_LIMIT = _limit_blas(int(getattr(a, "threads", 0) or 1))
    return a.func(a)


if __name__ == "__main__":
    sys.exit(main())
