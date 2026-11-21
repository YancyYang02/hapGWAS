
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


def _sha256(path, bs=1 << 22):
    import hashlib
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(bs), b""):
            h.update(b)
    return h.hexdigest()


def _md5(path, bs=1 << 22):
    import hashlib
    h = hashlib.md5()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(bs), b""):
            h.update(b)
    return h.hexdigest()


def _input_key(path, md5):
    """Unique provenance key of one input file: normalised absolute path + MD5."""
    return "%s#md5=%s" % (os.path.normcase(os.path.realpath(os.path.abspath(path))), md5)


def _input_record(path, role=None):
    md5 = _md5(path)
    rec = dict(path=os.path.realpath(os.path.abspath(path)), bytes=os.path.getsize(path),
               md5=md5, sha256=_sha256(path))
    if role:
        rec["role"] = role
    return _input_key(path, md5), rec


def _provenance(a, files, roles=None):
    """Software version, command line, environment and input hashes for summary.json.

    `inputs` is keyed by the normalised absolute path plus the MD5 of the file,
    so two inputs with the same file name in different directories (for
    example per-chromosome tables of one trait) are both recorded. A file passed
    twice is recorded once, with every role it was passed under."""
    import platform
    import scipy
    from . import __version__
    roles = list(roles) if roles is not None else [None] * len(files)
    inputs = {}
    for f, ro in zip(files, roles):
        if not f or not os.path.exists(f):
            continue
        k, rec = _input_record(f, ro)
        if k in inputs:
            prev = inputs[k].get("role")
            if ro and prev and ro not in prev.split(","):
                inputs[k]["role"] = prev + "," + ro
            continue
        inputs[k] = rec
    return dict(hapdc_version=__version__, argv=sys.argv, command=a.cmd,
                options={k: v for k, v in vars(a).items() if k != "func"},
                python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__,
                inputs=inputs, n_inputs=len(inputs))


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
        # every trait column of every file; a two-column file keeps its file stem as the trait name
        Yp, hn = io.read_pheno_table(p, keep_ids, miss=miss_codes)
        names += [os.path.splitext(os.path.basename(p))[0]] if Yp.shape[1] == 1 else list(hn)
        Y.append(Yp)
    Y = np.column_stack(Y)
    if len(set(names)) < len(names):
        raise SystemExit("duplicate trait names across phenotype files: %s" % names)
    print("markers dropped while reading the VCF: %s" % dropped, flush=True)
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
    # traits.npy holds only the prep samples with at least one trait value
    # (keep.npy indexes them in hap.npy / K_keep.npy); this trait further drops
    # its own missing individuals. The GRM is re-decomposed whenever the sample
    # set differs from the one eig.npz was computed for.
    kp = os.path.join(a.data, "keep.npy")
    keep_rows = np.load(kp) if os.path.exists(kp) else np.arange(K.shape[0])
    y = np.asarray(Y[:, ti], float)
    ok = np.isfinite(y)
    rows = keep_rows[ok]
    if len(rows) < K.shape[0]:
        m = H.shape[1]
        H = np.ascontiguousarray(H.reshape(K.shape[0], 2, m)[rows].reshape(-1, m))
        S, U = np.linalg.eigh(K[np.ix_(rows, rows)].astype(float))
        S = np.clip(S, 1e-6, None)
        y = y[ok]
        print("trait %d: %d of %d samples used" % (ti, len(rows), K.shape[0]), flush=True)
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


# lambda_gc of `ms` anchor p-values is below 1 under the null (the Cauchy
# combination is conservative away from the tail) by an amount that depends on
# the reference panel, the marker density, the anchor count and the transform.
# No generic null range is built in; a condition-matched range can be supplied
# with --lambda-ref LO,HI (e.g. the central 95% of lambda_gc over null
# phenotypes scanned with the same panel and options).

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
    bl = core.ld_blocks(H, p, dmin=a.ld_dmin, max_span=a.ld_maxspan, min_mk=2,
                        hap_cover=a.ld_hapcover, minf=a.minf)["blocks"]
    cov = np.zeros(len(cols), bool); units = []
    for s_, e_ in bl:
        units.append((int(s_), int(e_))); cov[s_:e_] = True
    units += [(int(j), int(j) + 1) for j in np.where(~cov)[0]]; units.sort()
    anc = M.anchors(units, a.anchor_min)
    meta = [g["markers"][j] for j in cols]
    sc = scan.Scan(H, meta, g["U"], g["S"], g["yt"], minf=a.minf, minc=a.minc, covar=g["C"],
                   covar_policy="error")
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
    cnames = []
    if a.covar:
        C, cnames = _read_covar(a.covar, ids)
        print("covariates:", cnames, flush=True)
    ok = np.isfinite(y) & (np.isfinite(C).all(1) if C is not None else True)
    idx = np.where(ok)[0]
    covar_dropped = []
    if C is not None:
        # Rank check of [intercept, covariates] on the samples actually analysed,
        # before any REML fit: a constant column (e.g. a level that became
        # constant after removing missing phenotypes, or a user-supplied
        # intercept) or a duplicated / linearly dependent column is dropped
        # (--covar-rank drop, default) or rejected (--covar-rank error).
        Ck, kn, covar_dropped = core.check_fixed_effects(C[idx], len(idx), names=cnames,
                                                         policy=a.covar_rank)
        for d_ in covar_dropped:
            print("covariate %s dropped before REML: %s" % (d_["name"], d_["reason"]), flush=True)
        if Ck is None:
            print("covariates: none left after the rank check; intercept only", flush=True)
            C = None
        else:
            C = np.full((len(ids), Ck.shape[1]), np.nan); C[idx] = Ck
        cnames = kn
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
    NA = len(rows)
    n_tests = int(a.n_tests) if a.n_tests else NA
    if a.chrom is not None and not a.n_tests:
        print("WARNING: --chrom scans a subset; the Bonferroni threshold below uses only its %d anchors. "
              "For a genome-wide threshold pass --n-tests <genome-wide anchor count>, or combine the "
              "per-chromosome tables with `hapdc merge`." % NA, flush=True)
    thr = 0.05 / n_tests
    _lam = float(core.lambda_gc(np.array([r_[4] for r_ in rows])))
    lam_ref = None
    if a.lambda_ref:
        lam_ref = tuple(float(x) for x in a.lambda_ref.split(","))
        if len(lam_ref) != 2 or not lam_ref[0] < lam_ref[1]:
            raise SystemExit("--lambda-ref expects LO,HI with LO < HI")
    print("lambda_gc of anchor p-values = %.3f. Under the null this is below 1 by an amount that depends "
          "on the panel, marker density, anchor count and transform; it is descriptive, not a calibration "
          "test (use a permutation threshold to check calibration)." % _lam, flush=True)
    if lam_ref is not None:
        print("lambda_gc %.3f is %s the user-supplied reference range %.3f-%.3f"
              % (_lam, "inside" if lam_ref[0] <= _lam <= lam_ref[1] else "OUTSIDE", lam_ref[0], lam_ref[1]),
              flush=True)
    path = os.path.join(out, "%s.anchors.tsv" % tname)
    with open(path, "w") as f:
        f.write("chrom\tstart_bp\tend_bp\tn_markers\tP\tlead_scale\tlead_unit_start\tlead_unit_end\tP_nonadditive\tsig_bonf\n")
        for rw in rows:
            f.write("%s\t%d\t%d\t%d\t%.4g\t%s\t%d\t%d\t%.4g\t%d\n" % (rw + (int(rw[4] < thr),)))
    json.dump(dict(trait=tname, n=int(len(idx)), n_anchors=NA, n_tests_bonf=n_tests, chrom=a.chrom, threshold_bonf=thr, n_sig=int(sum(r_[4] < thr for r_ in rows)),
                   lambda_gc=_lam, lambda_ref=(list(lam_ref) if lam_ref is not None else None),
                   anchor_min=a.anchor_min,
                   kernel_power=kp_g, kernel_power_lrt=kp_lrt, kernel_power_best=float(max(_tab, key=lambda t: t[1])[0]),
                   ld_dmin=a.ld_dmin, ld_maxspan=a.ld_maxspan, ld_hapcover=a.ld_hapcover, minf=a.minf, minc=a.minc,
                   covar=a.covar, covariates_used=cnames, covariates_dropped=covar_dropped,
                   covar_rank_policy=a.covar_rank, transform=("int" if (a.transform == "int" and len(np.unique(y[idx])) > 2) else "none"),
                   kernel_power_option=a.kernel_power, n_tests_source=("--n-tests" if a.n_tests else "anchors scanned"),
                   elapsed_s=round(time.time() - t0, 1),
                   provenance=_provenance(a, [os.path.join(a.data, f) for f in
                                              ("meta.json", "markers.json", "hap.npy", "traits.npy", "K_keep.npy", "keep.npy")]
                                          + ([a.covar] if a.covar else []),
                                          roles=["prep"] * 6 + (["covar"] if a.covar else []))),
              open(os.path.join(out, "%s.summary.json" % tname), "w"), indent=1)
    print("ms done: %d anchors, threshold %.3g, %d significant -> %s" % (NA, thr, sum(r_[4] < thr for r_ in rows), path))
    return 0


# ---------------------------------------------------------------- merge
def cmd_merge(a):
    """Combine per-chromosome `ms` anchor tables of one trait; genome-wide Bonferroni."""
    hdr, rows = None, []
    for p in a.tables:
        with open(p) as f:
            h = f.readline().rstrip("\n").split("\t")
            if hdr is None:
                hdr = h
            elif h != hdr:
                raise SystemExit("header of %s differs from %s" % (p, a.tables[0]))
            rows += [ln.rstrip("\n").split("\t") for ln in f if ln.strip()]
    key = [(r[0], int(r[1]), int(r[2])) for r in rows]
    if len(set(key)) < len(key):
        raise SystemExit("the same anchor appears in more than one table; pass each chromosome once")
    ip, isg = hdr.index("P"), hdr.index("sig_bonf")
    P = np.array([float(r[ip]) for r in rows])
    thr = a.alpha / len(rows)
    for r, pv in zip(rows, P):
        r[isg] = str(int(pv < thr))
    out = a.out if a.out.endswith(".tsv") else a.out + ".anchors.tsv"
    with open(out, "w") as f:
        f.write("\t".join(hdr) + "\n")
        for r in rows:
            f.write("\t".join(r) + "\n")
    lam = float(core.lambda_gc(P))
    # per-shard record: the shard's own summary.json (trait, n, transform,
    # software version, key options), when it sits next to the table
    shards, sfiles, sroles = [], [], []
    keys = ("trait", "n", "n_anchors", "transform", "covar", "covariates_used", "covariates_dropped",
            "anchor_min", "ld_dmin", "ld_maxspan", "ld_hapcover", "minf", "minc", "kernel_power", "chrom")
    for p in a.tables:
        sj = p[:-len(".anchors.tsv")] + ".summary.json" if p.endswith(".anchors.tsv") else None
        rec = dict(table=os.path.realpath(os.path.abspath(p)), table_md5=_md5(p))
        if sj and os.path.exists(sj):
            m_ = json.load(open(sj))
            rec.update({k: m_.get(k) for k in keys})
            rec["hapdc_version"] = (m_.get("provenance") or {}).get("hapdc_version")
            rec["summary"] = os.path.realpath(os.path.abspath(sj))
            sfiles.append(sj); sroles.append("shard_summary")
        shards.append(rec)
    traits = sorted({str(r["trait"]) for r in shards if r.get("trait") is not None})
    if len(traits) > 1:
        print("WARNING: merging tables of different traits: %s" % traits, flush=True)
    js = dict(tables=[os.path.realpath(os.path.abspath(p)) for p in a.tables], n_anchors=len(rows),
              threshold_bonf=thr, n_sig=int((P < thr).sum()), lambda_gc=lam, alpha=a.alpha, shards=shards,
              provenance=_provenance(a, list(a.tables) + sfiles, roles=["table"] * len(a.tables) + sroles))
    json.dump(js, open(out[:-len(".anchors.tsv")] + ".summary.json", "w"), indent=1)
    print("merge done: %d anchors, threshold %.3g, %d significant, lambda_gc %.3f -> %s"
          % (len(rows), thr, js["n_sig"], lam, out))
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
    """End-to-end check of the primary `ms` path (statistic D) and its invariances.

    Founder-haplotype blocks with known effects are scanned with the same
    functions `hapdc ms` calls (ld_blocks -> anchors -> ms_scan, variant D).
    Checks: (1) an additive and (2) a diplotype-specific (compound genotype)
    signal are each detected at their causal anchor below the Bonferroni level;
    (3) over 200 null phenotypes the family-wise rate at that level is <= 0.10; (4) the combined p of a
    phenotype is identical scanned alone or in a batch with another phenotype;
    (5) recoding the alleles of every marker and (6) reordering the samples
    leave every anchor p-value unchanged. Exits non-zero on any failure.
    """
    from . import ms as M
    rng = np.random.default_rng(7)
    n, nblk, blk = 600, 12, 8
    m = nblk * blk
    H = np.zeros((2 * n, m), np.uint8)
    for b in range(nblk):
        F = rng.integers(0, 2, (4, blk)).astype(np.uint8)
        F[0] = 0
        H[:, b * blk:(b + 1) * blk] = F[rng.choice(4, 2 * n, p=[.4, .3, .2, .1])]
    meta = [("1", "m%d" % j, 1000 * (j + 1), "A", "C") for j in range(m)]
    pos = np.array([x[2] for x in meta])
    S, U = np.ones(n), np.eye(n)

    class _Batch(scan.Scan):
        """Several phenotype columns scanned together (variance ratio of column 0)."""
        def __init__(self, Hh, Y):
            super().__init__(Hh, meta, U, S, Y[:, 0])
            self._Y = (U.T @ Y) * self.w[:, None]

        def _design_matrix(self, rows, B):
            return self._Y

    def scan_all(Hh, Y):
        bl = core.ld_blocks(Hh, pos, dmin=0.8, max_span=200000, min_mk=2)["blocks"]
        anc = M.anchors([(int(s), int(e)) for s, e in bl], 16)
        sc = _Batch(Hh, Y)
        r = M.ms_scan(sc, anc, variants=("D",), nadd=True, nadds=False)
        return anc, r["P_MS"]["D"]

    c0 = slice(2 * blk, 3 * blk)                         # causal block 2
    def hap_of(Hh, sl):
        F0 = (Hh[:, sl] == 0).all(1)
        return F0[0::2].astype(float), F0[1::2].astype(float)
    h1, h2 = hap_of(H, c0)
    fails = []
    tests = {"additive": h1 + h2, "diplotype": h1 * h2 * 1.0}
    for tag, g in tests.items():
        y = (0.45 if tag == "additive" else 1.2) * g + rng.normal(0, 1.0, n)
        y = (y - y.mean()) / y.std()
        anc, P = scan_all(H, y[:, None])
        thr = 0.05 / len(anc)
        hit = [k for k, (s, e, _) in enumerate(anc) if s <= c0.start and e >= c0.stop]
        lead = int(np.argmin(P[:, 0]))
        ok = bool(P[lead, 0] < thr and lead in hit)
        print("%-10s lead anchor %d (causal %s) p=%.3g  threshold %.3g  %s"
              % (tag, lead, hit, P[lead, 0], thr, "PASS" if ok else "FAIL"))
        if not ok:
            fails.append(tag)
    # (3) 200 null phenotypes scanned as one batch: family-wise rate at the
    # Bonferroni level (expected <= 0.05; fail above 0.10, i.e. > 3 binomial SE)
    Y0 = rng.normal(0, 1.0, (n, 200))
    anc, PN = scan_all(H, Y0)
    fw = float((PN.min(0) < 0.05 / len(anc)).mean())
    ok = fw <= 0.10
    print("null       family-wise rate %.3f over 200 null phenotypes (Bonferroni 0.05)  %s"
          % (fw, "PASS" if ok else "FAIL"))
    if not ok:
        fails.append("null")
    y0 = Y0[:, 0]
    P0 = PN[:, :1]
    # (4) batch invariance: same phenotype alone and next to a pure-noise column
    yb = np.column_stack([y0, rng.normal(0, 1.0, n) * 1e-3 + (H[0::2, 0] == H[1::2, 0])])
    _, Pb = scan_all(H, yb)
    d = float(np.max(np.abs(np.log(Pb[:, 0]) - np.log(P0[:, 0]))))
    ok = d < 1e-9
    print("batched    max |log p - log p0| = %.2e  %s" % (d, "PASS" if ok else "FAIL"))
    if not ok:
        fails.append("batch")
    # (5) allele recoding of every marker, (6) sample reordering
    perm = rng.permutation(n)
    Hf = (1 - H).astype(np.uint8)
    _, Pf = scan_all(Hf, y0[:, None])
    Hp = np.ascontiguousarray(H.reshape(n, 2, m)[perm].reshape(2 * n, m))
    _, Pp = scan_all(Hp, y0[perm][:, None])
    for tag, Q in (("recoded", Pf), ("reordered", Pp)):
        d = float(np.max(np.abs(np.log(Q[:, 0]) - np.log(P0[:, 0]))))
        ok = d < 1e-6
        print("%-10s max |log p - log p0| = %.2e  %s" % (tag, d, "PASS" if ok else "FAIL"))
        if not ok:
            fails.append(tag)
    print("selftest:", "all checks passed" if not fails else "FAILED: " + ", ".join(fails))
    return 0 if not fails else 1

def main(argv=None):
    p = argparse.ArgumentParser(prog="hapdc", description="Multi-scale haplotype and diplotype association scans in a linear mixed model")
    sub = p.add_subparsers(dest="cmd", required=True)

    q = sub.add_parser("prep", help="intersect VCF, GRM and phenotypes")
    q.add_argument("vcf"); q.add_argument("grm", help="GCTA .grm prefix, or 'auto' to compute a VanRaden GRM from the VCF"); q.add_argument("pheno", help="comma separated phenotype tables (header; IID then one column per trait)")
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
    q.add_argument("--nperm", type=int, default=1000,
                   help="null replicates by permuting whitened LMM residuals (variance ratio fixed at its estimate)")
    q.add_argument("--seed", type=int, default=20261110)
    q.add_argument("--lo", type=int, default=None); q.add_argument("--hi", type=int, default=None)
    q.add_argument("--quiet", action="store_true")
    q.add_argument("--trait-name", default=None)
    q.add_argument("--threads", type=int, default=0, help="BLAS threads (default 1)")
    q.set_defaults(func=cmd_scan)

    q = sub.add_parser("ms", help="multi-scale anchor scan (primary analysis, statistic D)")
    q.add_argument("data", help="directory written by `hapdc prep`"); q.add_argument("out")
    q.add_argument("trait", type=int, help="0-based trait column from prep")
    q.add_argument("--covar", default=None, help="tab-separated covariate table: IID then numeric columns, header line; the intercept is always added, do not include it")
    q.add_argument("--covar-rank", choices=["drop", "error"], default="drop",
                   help="covariate columns that are constant or linearly dependent on the intercept and earlier columns (checked on the analysed samples, before REML): drop and log them (default) or stop with an error")
    q.add_argument("--lambda-ref", default=None,
                   help="optional LO,HI condition-matched null range of lambda_gc (from null phenotypes scanned with the same panel and options); only reported, never used to change the threshold")
    q.add_argument("--transform", choices=["int", "none"], default="int",
                   help="phenotype transformation: rank-based inverse normal (default; skipped for binary traits) or none")
    q.add_argument("--chrom", default=None, help="comma separated subset of chromosomes")
    q.add_argument("--n-tests", type=int, default=None,
                   help="genome-wide anchor count for the Bonferroni threshold when --chrom scans a subset")
    q.add_argument("--anchor-min", type=int, default=16, help="minimum markers per anchor")
    q.add_argument("--ld-dmin", type=float, default=0.8); q.add_argument("--ld-maxspan", type=int, default=200000)
    q.add_argument("--ld-hapcover", type=float, default=0.8, help="split an LD block until haplotypes with frequency >= --minf cover this fraction of chromosomes")
    q.add_argument("--minf", type=float, default=0.05, help="haplotype frequency floor, used both for LD-block splitting and for haplotype classes")
    q.add_argument("--minc", type=int, default=5)
    q.add_argument("--chunk", type=int, default=500, help="anchors per batch (memory)")
    q.add_argument("--jobs", type=int, default=1, help="chromosomes scanned in parallel (POSIX fork); results identical to --jobs 1")
    q.add_argument("--threads", type=int, default=0, help="BLAS threads per worker (default 1; more threads give no speed-up for hapdc's small per-unit matrices)")
    q.add_argument("--trait-name", default=None)
    q.add_argument("--kernel-power", choices=["auto", "1"], default="1",
                   help="relationship kernel: 1 = standard GRM (default); auto = spectral power K^g selected by REML (LRT, 1 df, 0.05), sensitivity analysis")
    q.set_defaults(func=cmd_ms)

    q = sub.add_parser("merge", help="combine per-chromosome `ms` anchor tables with a genome-wide Bonferroni threshold")
    q.add_argument("out", help="output prefix"); q.add_argument("tables", nargs="+")
    q.add_argument("--alpha", type=float, default=0.05)
    q.set_defaults(func=cmd_merge)

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
