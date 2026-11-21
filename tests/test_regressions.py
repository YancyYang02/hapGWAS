"""Regression tests for hapdc. Run: python -m pytest tests  (or python tests/test_regressions.py)."""
import gzip
import json
import os
import subprocess
import sys
import tempfile

import numpy as np
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from hapdc import core, io, ms, scan  # noqa: E402

N = 120
PATS = np.array([[0, 0, 0, 0], [1, 1, 0, 0], [0, 1, 1, 0], [1, 0, 1, 1]], np.int8)
PAIRS = {0: (0, 1), 1: (2, 2), 2: (2, 3), 3: (3, 3)}
META = [("1", "m%d" % j, j + 1, "A", "C") for j in range(4)]


def _grm(seed=0, n=N):
    rng = np.random.default_rng(seed)
    X = rng.normal(size=(n, n // 4))
    K = X @ X.T / X.shape[1]
    S, U = np.linalg.eigh(K)
    return K, U, np.clip(S, 1e-6, None)


def _data(seed, p_ab=0.6):
    """Haplotypes A and B are always carried together, so dos_A == dos_B: an
    early linear dependency among the haplotype-dosage columns."""
    rng = np.random.default_rng(seed)
    typ = rng.choice(4, N, p=[p_ab, (1 - p_ab) * .4, (1 - p_ab) * .4, (1 - p_ab) * .2])
    H = np.zeros((2 * N, 4), np.int8)
    for i, t in enumerate(typ):
        a, b = PAIRS[t]
        H[2 * i], H[2 * i + 1] = PATS[a], PATS[b]
    y = 0.5 * (typ == 2) + 0.4 * (typ == 1) + rng.normal(size=N)
    return H, y, rng


def _ref(H, U, S, y):
    """Nested-model F tests from an SVD column basis (reference)."""
    sc = scan.Scan(H, META, U, S, y)
    d = core.window_design(H, 0.05, 5)
    Xd = (U.T @ d["dos"]) * sc.w[:, None]
    Xc = (U.T @ d["C"]) * sc.w[:, None]
    ys = sc.yr * sc.w

    def rss(X):
        u, s, _ = np.linalg.svd(X, full_matrices=False)
        Q = u[:, s > 1e-9 * s[0]]
        r = ys - Q @ (Q.T @ ys)
        return r @ r, Q.shape[1]

    one = sc.one[:, None]
    r0, k0 = rss(one)
    r1, k1 = rss(np.column_stack([one, Xd]))
    rc, kc = rss(np.column_stack([one, Xc]))
    rf, kf = rss(np.column_stack([one, Xd, Xc]))
    f = lambda a, b, ka, kb: stats.f.sf(((a - b) / (kb - ka)) / (b / (N - kb)), kb - ka, N - kb)
    return dict(A_joint=f(r0, r1, k0, k1), CELL_F=f(r0, rc, k0, kc), NADD_F=f(r1, rf, k1, kf))


def _p(H, U, S, y):
    r = scan.Scan(H, META, U, S, y).run(blocks=[(0, 4)], verbose=False)
    return {k: float(r["P"][k][0, 0]) for k in ("A_joint", "CELL_F", "NADD_F")}


def test_f_tests_match_reference_and_sample_order():
    K, U, S = _grm()
    for seed in range(20):
        H, y, rng = _data(seed, 0.6 if seed % 2 else 0.3)
        perm = rng.permutation(N)
        hp = np.concatenate([[2 * i, 2 * i + 1] for i in perm])
        Sp, Up = np.linalg.eigh(K[np.ix_(perm, perm)])
        a = _p(H, U, S, y)
        b = _p(H[hp], Up, np.clip(Sp, 1e-6, None), y[perm])
        r = _ref(H, U, S, y)
        for k in a:
            assert abs(np.log10(a[k]) - np.log10(r[k])) < 1e-6, (seed, k, a[k], r[k])
            assert abs(np.log10(b[k]) - np.log10(r[k])) < 1e-6, (seed, k, b[k], r[k])


def test_permutation_null_is_calibrated():
    """Replicates from permuted whitened residuals: p-values ~ U(0,1) under H0."""
    K, U, S = _grm(1, 200)
    rng = np.random.default_rng(5)
    L = np.linalg.cholesky(K + 1e-6 * np.eye(200))
    H = rng.integers(0, 2, size=(400, 1)).astype(np.int8)
    y = 2.0 * (L @ rng.normal(size=200)) + rng.normal(size=200)
    sc = scan.Scan(H, [("1", "m", 1, "A", "C")], U, S, y).permute(2000)
    YS = sc._design_matrix(None, 0)[:, 1:]
    g = sc._p0(((U.T @ (H[0::2, 0] + H[1::2, 0]).astype(float)) * sc.w)[:, None])[:, 0]
    r0 = sc._p0(YS)
    b = g @ r0 / (g @ g)
    t = b / np.sqrt(((r0 * r0).sum(0) - b * b * (g @ g)) / (200 - 2) / (g @ g))
    p = 2 * stats.t.sf(np.abs(t), 198)
    assert stats.kstest(p, "uniform").pvalue > 1e-3
    assert abs((p < 0.05).mean() - 0.05) < 0.02


def test_wacat_is_standard_acat_except_p_one():
    p, _, _ = ms.wacat([np.array([0.8]), np.array([0.9])], ["SNP", "SNP"])
    t = (np.tan((0.5 - 0.8) * np.pi) + np.tan((0.5 - 0.9) * np.pi)) / 2
    assert abs(p[0] - (0.5 - np.arctan(t) / np.pi)) < 1e-12
    p1, _, _ = ms.wacat([np.array([1.0]), np.array([1e-6])], ["SNP", "SNP"])
    assert p1[0] < 1e-5


def _write_vcf(path, rows, ids):
    with gzip.open(path, "wt") as f:
        f.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(ids) + "\n")
        for chrom, pos, ref, alt, gts in rows:
            f.write("%s\t%d\t.\t%s\t%s\t.\tPASS\t.\tGT\t%s\n" % (chrom, pos, ref, alt, "\t".join(gts)))


def test_vcf_reader_rejects_multiallelic_unphased_and_counts_nmax():
    d = tempfile.mkdtemp()
    ids = ["s1", "s2", "s3"]
    rows = [("1", 10, "A", "C,G", ["1|2", "0|0", "0|1"]),   # multi-allelic ALT
            ("1", 20, "A", "C", ["0/1", "0|0", "1|1"]),     # unphased heterozygote
            ("1", 30, "A", "C", ["0|2", "0|0", "1|1"]),     # allele index 2 with one ALT
            ("1", 40, "A", "C", ["1", "0|0", "1|1"]),       # haploid call
            ("1", 50, "A", "C", ["1|0", "0|0", "1|1"]),
            ("1", 60, "A", "C", ["1/1", "0|1", "0|0"])]     # unphased homozygote: fine
    p = os.path.join(d, "x.vcf.gz")
    _write_vcf(p, rows, ids)
    _, H, meta, dr = io.read_vcf(p)
    assert [m[2] for m in meta] == [50, 60]
    assert H[:, 0].tolist() == [1, 0, 0, 0, 1, 1]
    assert dr["multiallelic"] == 2 and dr["unphased"] == 1 and dr["nondiploid"] == 1
    rows = [("1", 100 + j, "A", "C", ["1|0", "0|0", "0|1"]) for j in range(5000)]
    _write_vcf(p, rows, ids)
    _, H, meta, _ = io.read_vcf(p, nmax=4500)
    assert H.shape == (6, 4500)


def _toy_prep(d, n=20, m=60, missing_all=("i7",), missing_t2=("i3",)):
    rng = np.random.default_rng(3)
    ids = ["i%d" % i for i in range(n)]
    F = rng.integers(0, 2, size=(4, m))
    rows = []
    hap = F[rng.integers(0, 4, size=2 * n)]
    for j in range(m):
        if hap[:, j].min() == hap[:, j].max():
            hap[0, j] = 1 - hap[0, j]
        rows.append(("1", 1000 * (j + 1), "A", "C", ["%d|%d" % (hap[2 * i, j], hap[2 * i + 1, j]) for i in range(n)]))
    _write_vcf(os.path.join(d, "g.vcf.gz"), rows, ids)
    with open(os.path.join(d, "ph.txt"), "w") as f:
        f.write("IID\tt1\tt2\n")
        for i in ids:
            v1 = "-9" if i in missing_all else "%.3f" % rng.normal()
            v2 = "-9" if (i in missing_all or i in missing_t2) else "%.3f" % rng.normal()
            f.write("%s\t%s\t%s\n" % (i, v1, v2))
    return ids


def _cli(*args):
    env = dict(os.environ, PYTHONPATH=os.path.dirname(HERE) + os.pathsep + os.environ.get("PYTHONPATH", ""))
    r = subprocess.run([sys.executable, "-c", "import sys; from hapdc.cli import main; sys.exit(main(sys.argv[1:]))"]
                       + list(args), capture_output=True, text=True, env=env)
    assert r.returncode == 0, r.stderr[-2000:]
    return r.stdout


def test_prep_multitrait_and_scan_after_missing():
    d = tempfile.mkdtemp()
    _toy_prep(d)
    _cli("prep", os.path.join(d, "g.vcf.gz"), "auto", os.path.join(d, "ph.txt"), os.path.join(d, "prep"), "--grm-thin", "1")
    mj = json.load(open(os.path.join(d, "prep", "meta.json")))
    assert mj["traits"] == ["t1", "t2"] and mj["n"] == 19
    Y = np.load(os.path.join(d, "prep", "traits.npy"))
    assert Y.shape == (19, 2) and np.isnan(Y[:, 1]).sum() == 1
    for t in ("0", "1"):
        _cli("scan", os.path.join(d, "prep"), os.path.join(d, "scan" + t), t, "--units", "fixed",
             "--wsize", "4", "--nperm", "0", "--minc", "1", "--minf", "0.05", "--quiet")


def test_merge_genome_wide_threshold():
    d = tempfile.mkdtemp()
    hdr = "chrom\tstart_bp\tend_bp\tn_markers\tP\tlead_scale\tlead_unit_start\tlead_unit_end\tP_nonadditive\tsig_bonf\n"
    for c, ps in (("1", [1e-3, 0.5]), ("2", [0.02, 0.2])):  # 0.02: significant per chromosome (0.025), not genome-wide (0.0125)
        with open(os.path.join(d, "c%s.tsv" % c), "w") as f:
            f.write(hdr)
            for k, p in enumerate(ps):
                f.write("%s\t%d\t%d\t20\t%g\tSNP\t1\t2\t1\t%d\n" % (c, 10 * k, 10 * k + 5, p, int(p < 0.05 / 2)))
    _cli("merge", os.path.join(d, "all"), os.path.join(d, "c1.tsv"), os.path.join(d, "c2.tsv"))
    js = json.load(open(os.path.join(d, "all.summary.json")))
    assert js["n_anchors"] == 4 and abs(js["threshold_bonf"] - 0.0125) < 1e-12 and js["n_sig"] == 1


def test_projection_cache_is_bounded():
    rng = np.random.default_rng(0)
    n, m = 50, 2048 * 20
    H = rng.integers(0, 2, size=(2 * n, m)).astype(np.int8)
    K, U, S = _grm(2, n)
    sc = scan.Scan(H, [("1", "m%d" % j, j, "A", "C") for j in range(m)], U, S, rng.normal(size=n))
    for s in range(0, m - 100, 1500):
        sc._projG(s, s + 100)
    assert len(sc._UG) <= sc.cache_blocks
    X = sc._projG(5000, 5100)
    assert np.allclose(X, U.T @ (H[0::2, 5000:5100] + H[1::2, 5000:5100]).astype(float))


# ---------------------------------------------------------------- audit items

def _founders(seed, n=300, nblk=6, blk=6):
    rng = np.random.default_rng(seed)
    m = nblk * blk
    H = np.zeros((2 * n, m), np.uint8)
    for b in range(nblk):
        F = rng.integers(0, 2, (4, blk)).astype(np.uint8)
        H[:, b * blk:(b + 1) * blk] = F[rng.choice(4, 2 * n, p=rng.dirichlet(np.ones(4) * 2))]
    H = (H ^ (rng.random(H.shape) < 0.01)).astype(np.uint8)
    pos = np.cumsum(rng.integers(100, 3000, m))
    return H, pos, rng


def _proj(X):
    from scipy.linalg import orth
    Q = orth(np.column_stack([np.ones(X.shape[0]), X]))
    return Q @ Q.T


def test_window_design_is_invariant_to_allele_coding_and_sample_order():
    """S02: rare-haplotype and small-cell ties are shared, not broken by byte order."""
    # auditor's example: two equally frequent retained haplotypes, rare 01 equidistant
    H = np.repeat(np.array([[0, 0], [1, 1], [0, 1]], np.uint8), [98, 98, 4], axis=0)
    Hf = H.copy(); Hf[:, 0] = 1 - Hf[:, 0]
    a, b = core.window_design(H), core.window_design(Hf)
    assert np.allclose(_proj(a["dos"]), _proj(b["dos"]))
    assert np.allclose(_proj(a["C"]), _proj(b["C"]))
    rng = np.random.default_rng(1)
    for _ in range(150):
        nf = rng.integers(2, 6)
        F = rng.integers(0, 2, (nf, rng.integers(2, 8))).astype(np.uint8)
        n = int(rng.integers(30, 200))
        Hh = F[rng.choice(nf, 2 * n, p=rng.dirichlet(np.ones(nf)))]
        Hh = (Hh ^ (rng.random(Hh.shape) < 0.01)).astype(np.uint8)
        fl = rng.random(Hh.shape[1]) < .5
        perm = rng.permutation(n)
        Hp = Hh.copy(); Hp[:, fl] = 1 - Hp[:, fl]
        Hp = Hp.reshape(n, 2, -1)[perm][:, ::-1].reshape(2 * n, -1)
        a, b = core.window_design(Hh), core.window_design(Hp)
        # each individual sits in exactly one class pair and one cell in total
        assert np.allclose(a["dos"].sum(1), 2) and (a["ncell"] == 0 or np.allclose(a["C"].sum(1), 1))
        assert np.allclose(_proj(a["dos"])[np.ix_(perm, perm)], _proj(b["dos"]))
        assert np.allclose(_proj(a["C"])[np.ix_(perm, perm)], _proj(b["C"]))


def test_ld_blocks_are_invariant_to_allele_coding_and_sample_order():
    """S01/S12: |D'| and coverage from integer counts; minf reaches the splitter."""
    for seed in range(40):
        H, pos, rng = _founders(seed)
        n = H.shape[0] // 2
        fl = rng.random(H.shape[1]) < .5
        perm = rng.permutation(n)
        Hp = H.copy(); Hp[:, fl] = 1 - Hp[:, fl]
        Hp = Hp.reshape(n, 2, -1)[perm].reshape(2 * n, -1)
        for kw in (dict(), dict(hap_cover=0.9, minf=0.1)):
            a = core.ld_blocks(H, pos, dmin=0.8, max_span=200000, min_mk=2, **kw)["blocks"]
            b = core.ld_blocks(Hp, pos, dmin=0.8, max_span=200000, min_mk=2, **kw)["blocks"]
            assert [tuple(x) for x in a] == [tuple(x) for x in b]
    # minf changes the splitting (it used to be ignored by the ms worker)
    src = open(os.path.join(os.path.dirname(HERE), "hapdc", "cli.py"), encoding="utf-8").read()
    assert "hap_cover=a.ld_hapcover, minf=a.minf" in src


class _Batch(scan.Scan):
    def __init__(self, H, meta, U, S, Y):
        super().__init__(H, meta, U, S, Y[:, 0])
        self._Y = (U.T @ Y) * self.w[:, None]

    def _design_matrix(self, rows, B):
        return self._Y


def test_ms_scan_does_not_depend_on_the_batch():
    """S04: components are selected by genotype-only validity, not by batch p-values."""
    H, pos, rng = _founders(3)
    n = H.shape[0] // 2
    H[:, 5] = 0                                    # a monomorphic marker inside a block
    meta = [("1", "m%d" % j, int(pos[j]), "A", "C") for j in range(H.shape[1])]
    K, U, S = _grm(4, n)
    bl = core.ld_blocks(H, pos, dmin=0.8, max_span=200000, min_mk=2)["blocks"]
    anc = ms.anchors([(int(s), int(e)) for s, e in bl], 16)
    y = rng.normal(size=n)
    alone = ms.ms_scan(_Batch(H, meta, U, S, y[:, None]), anc, variants=("D", "O"), nadd=True, nadds=True)
    Y = np.column_stack([y, rng.normal(size=(n, 7)), np.where(np.arange(n) % 2, 1.0, 0.0)])
    batch = ms.ms_scan(_Batch(H, meta, U, S, Y), anc, variants=("D", "O"), nadd=True, nadds=True)
    for v in ("D", "O"):
        assert np.allclose(np.log(alone["P_MS"][v][:, 0]), np.log(batch["P_MS"][v][:, 0]), atol=1e-9)
    assert np.allclose(alone["P_NADD"][:, 0], batch["P_NADD"][:, 0], atol=1e-12)


def test_ms_scan_is_invariant_to_allele_coding():
    H, pos, rng = _founders(5)
    n = H.shape[0] // 2
    meta = [("1", "m%d" % j, int(pos[j]), "A", "C") for j in range(H.shape[1])]
    K, U, S = _grm(6, n)
    y = rng.normal(size=n)
    out = []
    for Hh in (H, (1 - H).astype(np.uint8)):
        bl = core.ld_blocks(Hh, pos, dmin=0.8, max_span=200000, min_mk=2)["blocks"]
        anc = ms.anchors([(int(s), int(e)) for s, e in bl], 16)
        out.append(ms.ms_scan(scan.Scan(Hh, meta, U, S, y), anc, variants=("D",))["P_MS"]["D"][:, 0])
    assert np.allclose(np.log(out[0]), np.log(out[1]), atol=1e-8)


def test_lmm_scan_additive_uses_rotated_intercept():
    """S14: the intercept is U'1 on the eigen scale, not the weights alone."""
    rng = np.random.default_rng(2)
    n = 150
    K, U, S = _grm(2, n)
    G = rng.integers(0, 3, (n, 4)).astype(float)
    y = 5.0 + rng.normal(size=n)
    t, df = core.lmm_scan_additive(U.T @ G, (U.T @ y)[:, None], S, 0.7, U.T @ np.ones(n))
    L = np.linalg.cholesky(K + 0.7 * np.eye(n) + 1e-9 * np.eye(n))
    Li = np.linalg.inv(L)
    for j in range(4):
        X = Li @ np.column_stack([np.ones(n), G[:, j]])
        b, *_ = np.linalg.lstsq(X, Li @ y, rcond=None)
        r = Li @ y - X @ b
        se = np.sqrt(r @ r / (n - 2) * np.linalg.inv(X.T @ X)[1, 1])
        assert abs(t[j, 0] - b[1] / se) < 1e-6


def test_ms_scan_batch_invariance_at_p_equal_one():
    """S04, auditor's construction: a valid A_joint whose observed p rounds to 1."""
    from scipy import linalg
    rng = np.random.default_rng(704); n = 160
    pats = np.array([[0, 0], [0, 1], [1, 0], [1, 1]], np.uint8)
    H = pats[rng.integers(0, 4, 2 * n)]
    d = core.window_design(H)
    Qa = linalg.orth(np.column_stack([np.ones(n), d["dos"]]))
    y = d["C"][:, 1] + .08 * rng.normal(size=n)
    y = y - Qa @ (Qa.T @ y)
    other = d["dos"][:, 0] + rng.normal(size=n)
    meta = [("1", str(j), j, "A", "G") for j in range(2)]
    U, S = np.eye(n), np.ones(n)
    anc = [(0, 2, [(0, 2)])]
    pa = ms.ms_scan(_Batch(H, meta, U, S, y[:, None]), anc, variants=("D",), nadds=False)["P_MS"]["D"][0, 0]
    pb = ms.ms_scan(_Batch(H, meta, U, S, np.column_stack([y, other])), anc, variants=("D",), nadds=False)["P_MS"]["D"][0, 0]
    assert abs(np.log(pa) - np.log(pb)) < 1e-9, (pa, pb)


# --------------------------------------------------------------- v0.6.6 fixes
def _v01_counterexample():
    """Reviewer's construction (independent_checks, seed 20261121): n = 120, one
    marker with MAF 0.129, covariate g + 5.416e-5 * z, nearly collinear with the
    marker dosage. Before the fix the testability reference norm used the whitening
    vector w instead of the whitened intercept (U'1) * w, and the marker was
    testable under one allele coding and not under the other."""
    rng = np.random.default_rng(20261121)
    n = 120
    H = rng.binomial(1, .15, (2 * n, 1)).astype(np.int8)
    g = (H[0::2, 0] + H[1::2, 0]).astype(float)
    U, _ = np.linalg.qr(rng.normal(size=(n, n)))
    S = np.linspace(.3, 2, n)
    z = rng.normal(size=n)
    y = z + .2 * rng.normal(size=n)
    eps = None
    for e in np.geomspace(1e-6, .001, 200):          # same grid as the review
        if abs(e - 5.416e-5) / 5.416e-5 < 1e-3:
            eps = e
            break
    assert eps is not None
    C = (g + eps * z)[:, None]
    return H, g, U, S, y, C, z


def test_v01_testable_is_invariant_to_allele_coding_with_collinear_covariate():
    H, g, U, S, y, C, _ = _v01_counterexample()
    assert abs(g.mean() / 2 - 0.129) < 5e-4, g.mean() / 2
    meta = [("1", "m1", 100, "A", "C")]
    res = []
    for HH in (H, 1 - H):
        sc = scan.Scan(HH, meta, U, S, y, covar=C)
        p, v = ms._snp_p(sc, np.array([0]))
        pd = ms.ms_scan(sc, [(0, 1, [(0, 1)])], variants=("D",))["P_MS"]["D"][0, 0]
        res.append((bool(v[0]), p[0, 0], pd))
    assert res[0][0] == res[1][0], res
    assert abs(np.log(res[0][1]) - np.log(res[1][1])) < 1e-6, res
    assert abs(np.log(res[0][2]) - np.log(res[1][2])) < 1e-6, res


def test_v01_testable_reference_is_whitened_intercept():
    H, g, U, S, y, C, _ = _v01_counterexample()
    sc = scan.Scan(H, [("1", "m1", 100, "A", "C")], U, S, y, covar=C)
    o = (U.T @ np.ones(len(y))) * sc.w
    assert np.allclose(sc.o, o)
    # Q0 is always the basis used by _p0, also with intercept only
    sc1 = scan.Scan(H, [("1", "m1", 100, "A", "C")], U, S, y)
    X = np.random.default_rng(3).normal(size=(len(y), 3))
    o1 = sc1.o / np.linalg.norm(sc1.o)
    assert sc1.nfix == 1 and np.allclose(sc1.Q0[:, 0], o1)
    assert np.allclose(sc1._p0(X), X - np.outer(o1, o1 @ X))


def test_v01_eigenvector_sign_flip_and_allele_flip():
    """Flipping the sign of GRM eigenvectors, alone or together with recoding
    the alleles, changes neither testability nor the p-values."""
    H, g, U, S, y, C, _ = _v01_counterexample()
    meta = [("1", "m1", 100, "A", "C")]
    Uf = U.copy()
    Uf[:, ::3] *= -1
    base = None
    for UU in (U, Uf):
        for HH in (H, 1 - H):
            for CC in (None, C):
                sc = scan.Scan(HH, meta, UU, S, y, covar=CC)
                p, v = ms._snp_p(sc, np.array([0]))
                key = CC is None
                if base is None:
                    base = {}
                if key not in base:
                    base[key] = (bool(v[0]), p[0, 0], sc.delta)
                b = base[key]
                assert bool(v[0]) == b[0]
                assert abs(np.log(p[0, 0]) - np.log(b[1])) < 1e-6, (p[0, 0], b[1])
                assert abs(sc.delta - b[2]) < 1e-6 * max(1, b[2])
    # and the multi-marker window statistics
    Hm, ym, _ = _data(11)
    K, U2, S2 = _grm(5)
    U2f = U2.copy()
    U2f[:, 1::2] *= -1
    cov = np.random.default_rng(9).normal(size=(N, 1))
    ra = scan.Scan(Hm, META, U2, S2, ym, covar=cov).run(blocks=[(0, 4)], verbose=False)
    rb = scan.Scan(1 - Hm, META, U2f, S2, ym, covar=cov).run(blocks=[(0, 4)], verbose=False)
    for q in ("A_joint", "CELL_F", "NADD_F", "SACAT"):
        a_, b_ = float(ra["P"][q][0, 0]), float(rb["P"][q][0, 0])
        assert np.isfinite(a_) == np.isfinite(b_), q
        if np.isfinite(a_):
            assert abs(np.log(a_) - np.log(b_)) < 1e-6, (q, a_, b_)


def test_v02_covariate_rank_check_before_reml():
    H, g, U, S, y, C, z = _v01_counterexample()
    n = len(y)
    meta = [("1", "m1", 100, "A", "C")]
    ref = scan.Scan(H, meta, U, S, y, covar=z[:, None])
    p_ref = ms._snp_p(ref, np.array([0]))[0][0, 0]
    p_int = ms._snp_p(scan.Scan(H, meta, U, S, y), np.array([0]))[0][0, 0]
    cases = {
        "constant": (np.ones((n, 1)), ["covar1"], p_int, 1),
        "constant_5": (np.full((n, 1), 5.0), ["covar1"], p_int, 1),
        "duplicated": (np.column_stack([z, z]), ["covar2"], p_ref, 2),
        "scaled_duplicate": (np.column_stack([z, 3 * z + 1]), ["covar2"], p_ref, 2),
        "constant_then_valid": (np.column_stack([np.ones(n), z]), ["covar1"], p_ref, 2),
    }
    for lab, (CC, drop, pexp, nfix) in cases.items():
        sc = scan.Scan(H, meta, U, S, y, covar=CC)
        assert [d_["name"] for d_ in sc.covar_dropped] == drop, (lab, sc.covar_dropped)
        assert sc.nfix == nfix, lab
        p = ms._snp_p(sc, np.array([0]))[0][0, 0]
        assert abs(np.log(p) - np.log(pexp)) < 1e-8, (lab, p, pexp)
        try:
            scan.Scan(H, meta, U, S, y, covar=CC, covar_policy="error")
            raise AssertionError("no error for " + lab)
        except core.FixedEffectRankError as e:
            assert drop[0] in str(e)
    assert core.check_fixed_effects(np.ones((n, 1)), n)[2][0]["reason"] == "constant"
    assert core.check_fixed_effects(np.column_stack([z, z]), n, names=["pc1", "pc1b"])[2][0]["reason"] == "duplicate of pc1"
    # a nearly collinear but full-rank covariate (the V01 case) is kept
    assert core.check_fixed_effects(C, n)[2] == []
    # REML itself refuses a singular design instead of raising LinAlgError
    try:
        core.reml_delta(U.T @ y, S, U.T @ np.column_stack([np.ones(n), np.ones(n)]))
        raise AssertionError("reml_delta accepted a singular design")
    except core.FixedEffectRankError:
        pass
    try:
        core.check_fixed_effects(np.column_stack([z, np.r_[np.nan, z[1:]]]), n)
        raise AssertionError("NaN accepted")
    except ValueError:
        pass


def _cli_raw(*args):
    env = dict(os.environ, PYTHONPATH=os.path.dirname(HERE) + os.pathsep + os.environ.get("PYTHONPATH", ""))
    return subprocess.run([sys.executable, "-c", "import sys; from hapdc.cli import main; sys.exit(main(sys.argv[1:]))"]
                          + list(args), capture_output=True, text=True, env=env)


def test_v02_cli_drops_and_logs_rank_deficient_covariates():
    rng = np.random.default_rng(31)
    n, m = 80, 40
    H = rng.binomial(1, .3, (2 * n, m)).astype(np.int8)
    d = tempfile.mkdtemp()
    ids = ["s%d" % i for i in range(n)]
    rows = [("1", 1000 * (j + 1), "A", "C", ["%d|%d" % (H[2 * i, j], H[2 * i + 1, j]) for i in range(n)])
            for j in range(m)]
    _write_vcf(os.path.join(d, "g.vcf.gz"), rows, ids)
    y = H[0::2, 3] + H[1::2, 3] + rng.normal(size=n)
    with open(os.path.join(d, "p.txt"), "w") as f:
        f.write("IID\tt\n" + "".join("%s\t%.6f\n" % (i, v) for i, v in zip(ids, y)))
    z = rng.normal(size=n)
    with open(os.path.join(d, "c.txt"), "w") as f:
        f.write("IID\tsex\tpc1\tpc1copy\n" + "".join("%s\t1\t%.6f\t%.6f\n" % (i, v, v) for i, v in zip(ids, z)))
    with open(os.path.join(d, "c1.txt"), "w") as f:
        f.write("IID\tpc1\n" + "".join("%s\t%.6f\n" % (i, v) for i, v in zip(ids, z)))
    pr = os.path.join(d, "prep")
    _cli("prep", os.path.join(d, "g.vcf.gz"), "auto", os.path.join(d, "p.txt"), pr, "--grm-thin", "1")
    o1 = _cli("ms", pr, os.path.join(d, "o1"), "0", "--anchor-min", "8", "--covar", os.path.join(d, "c.txt"))
    assert "covariate sex dropped before REML: constant" in o1, o1[-1500:]
    assert "covariate pc1copy dropped before REML: duplicate of pc1" in o1
    _cli("ms", pr, os.path.join(d, "o2"), "0", "--anchor-min", "8", "--covar", os.path.join(d, "c1.txt"))
    j1 = json.load(open(os.path.join(d, "o1", "p.summary.json")))
    assert j1["covariates_used"] == ["pc1"]
    assert [x["name"] for x in j1["covariates_dropped"]] == ["sex", "pc1copy"]
    t1 = open(os.path.join(d, "o1", "p.anchors.tsv")).read()
    t2 = open(os.path.join(d, "o2", "p.anchors.tsv")).read()
    assert t1 == t2
    r3 = _cli_raw("ms", pr, os.path.join(d, "o3"), "0", "--anchor-min", "8", "--covar", os.path.join(d, "c.txt"),
                  "--covar-rank", "error")
    assert r3.returncode != 0 and "sex (constant)" in (r3.stderr + r3.stdout), r3.stderr[-800:]
    # V15: no hardcoded lambda reference range; an optional user range is only reported
    assert "0.61" not in o1 and "WARNING: lambda_gc" not in o1
    assert j1["lambda_ref"] is None
    o4 = _cli("ms", pr, os.path.join(d, "o4"), "0", "--anchor-min", "8", "--lambda-ref", "0.5,0.9")
    assert "user-supplied reference range 0.500-0.900" in o4
    assert json.load(open(os.path.join(d, "o4", "p.summary.json")))["lambda_ref"] == [0.5, 0.9]


def test_v03_provenance_keys_are_unique_per_path_and_content():
    from types import SimpleNamespace
    from hapdc.cli import _provenance
    with tempfile.TemporaryDirectory() as d:
        paths = []
        for chrom in ("chr1", "chr2", "chr3"):
            os.mkdir(os.path.join(d, chrom))
            p = os.path.join(d, chrom, "height.anchors.tsv")
            open(p, "w").write(chrom if chrom != "chr3" else "chr1")   # chr3 has the same content as chr1
            paths.append(p)
        pr = _provenance(SimpleNamespace(cmd="merge"), paths + [paths[0]])
        assert pr["n_inputs"] == 3 and len(pr["inputs"]) == 3, pr["inputs"]
        recs = list(pr["inputs"].values())
        assert sorted(r["path"] for r in recs) == sorted(os.path.realpath(p) for p in paths)
        for k, r in pr["inputs"].items():
            assert k.endswith("#md5=" + r["md5"]) and len(r["md5"]) == 32 and len(r["sha256"]) == 64
        assert recs[0]["md5"] == recs[2]["md5"] or len({r["md5"] for r in recs}) == 2


def test_v03_merge_records_every_shard():
    with tempfile.TemporaryDirectory() as d:
        tabs = []
        hdr = "chrom\tstart_bp\tend_bp\tn_markers\tP\tlead_scale\tlead_unit_start\tlead_unit_end\tP_nonadditive\tsig_bonf\n"
        for k, chrom in enumerate(("1", "2")):
            sub_ = os.path.join(d, "chr" + chrom)
            os.mkdir(sub_)
            p = os.path.join(sub_, "height.anchors.tsv")
            open(p, "w").write(hdr + "%s\t1\t100\t16\t%.3g\tSNP\t1\t1\t0.5\t0\n" % (chrom, 0.01 * (k + 1)))
            json.dump(dict(trait="height", n=100 + k, transform="rank_inverse_normal", covar=None,
                           provenance=dict(hapdc_version="test")), open(os.path.join(sub_, "height.summary.json"), "w"))
            tabs.append(p)
        out = os.path.join(d, "merged")
        _cli("merge", out, *tabs)
        m = json.load(open(out + ".summary.json"))
        assert len(m["shards"]) == 2 and [s_["n"] for s_ in m["shards"]] == [100, 101]
        assert m["provenance"]["n_inputs"] == 4


if __name__ == "__main__":
    for k, f in list(globals().items()):
        if k.startswith("test_"):
            f()
            print("ok", k)
