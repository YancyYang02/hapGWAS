
"""hapdc.io -- readers and writers for phased genotypes, GRM and phenotypes."""
import gzip
import numpy as np

__all__ = ["iter_vcf", "read_vcf", "read_grm", "read_pheno", "write_matrix"]


def iter_vcf(path, chunk=1):
    """Yield (sample_ids, hap_matrix, meta) per chunk of markers.

    hap_matrix is (n_hap, n_marker) int8 with rows 2i, 2i+1 the two haplotypes of
    sample i, matching the convention in this project. Markers with missing
    genotypes or no variation are dropped, so chunks are not fixed length.
    """
    ids = None
    G, pos, mk, nmiss, nmono = [], [], [], 0, 0
    with gzip.open(path, "rt") as f:
        for ln in f:
            if ln.startswith("##"):
                continue
            p = ln.rstrip("\n").split("\t")
            if ln.startswith("#CHROM"):
                ids = p[9:]
                continue
            gts = [g.split(":", 1)[0] for g in p[9:]]
            if any("." in g for g in gts):
                nmiss += 1
                continue
            row = np.empty(2 * len(ids), np.int8)
            for j, gt in enumerate(gts):
                a = gt.replace("/", "|").split("|")
                row[2 * j] = 1 if a[0] == "1" else 0
                row[2 * j + 1] = 1 if (len(a) > 1 and a[1] == "1") else 0
            if row.min() == row.max():
                nmono += 1
                continue
            G.append(row)
            pos.append(int(p[1]))
            mk.append((p[0], p[2] if p[2] not in (".", "") else "%s:%s" % (p[0], p[1]),
                       int(p[1]), p[3], p[4].split(",")[0]))
            if len(G) >= chunk:
                yield ids, np.array(G), mk, (nmiss, nmono)
                G, pos, mk, nmiss, nmono = [], [], [], 0, 0
    if G:
        yield ids, np.array(G), mk, (nmiss, nmono)


def read_vcf(path, nmax=None):
    """Read a whole phased VCF. Returns (ids, H (n_hap, m) int8, meta, dropped)."""
    ids = None
    outs = []
    nmiss = nmono = 0
    for cids, G, mk, dr in iter_vcf(path, chunk=2000):
        if ids is None:
            ids = cids
        outs.append((G, mk))
        nmiss += dr[0]
        nmono += dr[1]
        if nmax and sum(o[0].shape[1] for o in outs) >= nmax:
            break
    if not outs:
        raise SystemExit("no usable markers in %s" % path)
    # chunk buffers are (marker, haplotype): concatenate on the marker axis
    H = np.concatenate([o[0] for o in outs], axis=0)
    meta = [m for o in outs for m in o[1]]
    if nmax:
        H, meta = H[:nmax, :], meta[:nmax]
    # internal buffers are (marker, haplotype); the public API is (haplotype, marker)
    return ids, np.ascontiguousarray(H.T), meta, dict(missing=nmiss, monomorphic=nmono)


def read_grm(pref):
    """Read a GCTA .grm.{id,bin,n} triple into (ids, K)."""
    ids = [l.split()[1] for l in open(pref + ".grm.id")]
    n = len(ids)
    v = np.fromfile(pref + ".grm.bin", dtype="<f4")
    K = np.zeros((n, n))
    K[np.tril_indices(n)] = v
    K = K + K.T - np.diag(np.diag(K))
    return ids, K


def read_pheno(path, ids, col=None, miss=(-9.0,)):
    """Read a two-column (IID, PHENO) phenotype file aligned to ids.

    Values in `miss` (default -9, the study missing code) become NaN so that
    each trait drops only its own missing individuals.
    """
    bad = set(float(v) for v in miss)
    tab = {}
    with open(path) as f:
        head = f.readline().rstrip("\n").split("\t")
        idx = 1 if col is None else col
        for ln in f:
            p = ln.rstrip("\n").split("\t")
            if len(p) <= idx:
                continue
            try:
                v = float(p[idx])
            except ValueError:
                continue
            if v in bad:
                tab[p[0]] = np.nan
            elif np.isfinite(v):
                tab[p[0]] = v
    y = np.array([tab.get(i, np.nan) for i in ids], float)
    return y


def write_matrix(path, M, ids=None, renames=None):
    np.save(path, M)
    return path
