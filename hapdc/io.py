
"""hapdc.io -- readers and writers for phased genotypes, GRM and phenotypes."""
import gzip
import numpy as np

__all__ = ["iter_vcf", "read_vcf", "read_grm", "read_pheno", "read_pheno_table", "write_matrix"]


def iter_vcf(path, chunk=1):
    """Yield (sample_ids, hap_matrix, meta, dropped) per chunk of markers.

    hap_matrix is (n_marker, n_hap) int8 inside a chunk; read_vcf returns the
    public (n_hap, n_marker) layout with rows 2i, 2i+1 the two haplotypes of
    sample i. Only biallelic, fully phased, diploid, non-missing markers are
    used. Markers are dropped (and counted in `dropped`) when they are
    multi-allelic (ALT lists more than one allele, or a genotype carries an
    allele index other than 0/1), carry a missing call, an unphased
    heterozygote or a haploid call, or show no variation.
    """
    ids = None
    keys = ("missing", "monomorphic", "multiallelic", "unphased", "nondiploid")
    G, mk = [], []
    dr = dict.fromkeys(keys, 0)
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as f:
        for ln in f:
            if ln.startswith("##"):
                continue
            p = ln.rstrip("\n").split("\t")
            if ln.startswith("#CHROM"):
                ids = p[9:]
                continue
            if "," in p[4]:
                dr["multiallelic"] += 1
                continue
            gts = [g.split(":", 1)[0] for g in p[9:]]
            row = np.empty(2 * len(ids), np.int8)
            bad = None
            for j, gt in enumerate(gts):
                if "|" in gt:
                    a = gt.split("|")
                elif "/" in gt:
                    a = gt.split("/")
                    if len(a) == 2 and a[0] != a[1] and "." not in a:
                        bad = "unphased"
                        break
                else:
                    a = [gt]
                if "." in a:
                    bad = "missing"
                    break
                if len(a) != 2:
                    bad = "nondiploid"
                    break
                if a[0] not in ("0", "1") or a[1] not in ("0", "1"):
                    bad = "multiallelic"
                    break
                row[2 * j] = a[0] == "1"
                row[2 * j + 1] = a[1] == "1"
            if bad is not None:
                dr[bad] += 1
                continue
            if row.min() == row.max():
                dr["monomorphic"] += 1
                continue
            G.append(row)
            mk.append((p[0], p[2] if p[2] not in (".", "") else "%s:%s" % (p[0], p[1]),
                       int(p[1]), p[3], p[4]))
            if len(G) >= chunk:
                yield ids, np.array(G), mk, dr
                G, mk = [], []
                dr = dict.fromkeys(keys, 0)
    if G or any(dr.values()):
        yield ids, (np.array(G) if G else np.zeros((0, 2 * len(ids or [])), np.int8)), mk, dr


def read_vcf(path, nmax=None):
    """Read a whole phased VCF. Returns (ids, H (n_hap, m) int8, meta, dropped).

    `dropped` counts markers skipped per reason (see iter_vcf). nmax, if given,
    caps the number of markers read.
    """
    ids = None
    outs = []
    dropped = {}
    nread = 0
    for cids, G, mk, dr in iter_vcf(path, chunk=2000):
        if ids is None:
            ids = cids
        if len(mk):
            outs.append((G, mk))
            nread += G.shape[0]          # chunk buffers are (marker, haplotype)
        for k_, v_ in dr.items():
            dropped[k_] = dropped.get(k_, 0) + v_
        if nmax and nread >= nmax:
            break
    if not outs:
        raise SystemExit("no usable markers in %s (dropped: %s)" % (path, dropped))
    H = np.concatenate([o[0] for o in outs], axis=0)
    meta = [m for o in outs for m in o[1]]
    if nmax:
        H, meta = H[:nmax, :], meta[:nmax]
    # internal buffers are (marker, haplotype); the public API is (haplotype, marker)
    return ids, np.ascontiguousarray(H.T), meta, dropped


def read_grm(pref):
    """Read a GCTA .grm.{id,bin,n} triple into (ids, K)."""
    ids = [l.split()[1] for l in open(pref + ".grm.id")]
    n = len(ids)
    v = np.fromfile(pref + ".grm.bin", dtype="<f4")
    K = np.zeros((n, n))
    K[np.tril_indices(n)] = v
    K = K + K.T - np.diag(np.diag(K))
    return ids, K


def _split(ln):
    ln = ln.rstrip("\r\n")
    return ln.split("\t") if "\t" in ln else ln.split()


def read_pheno_table(path, ids, miss=(-9.0,)):
    """Read a phenotype table: header line, first column IID, one column per trait.

    Returns (Y (len(ids), k) float with NaN for missing, trait names). Values in
    `miss` and non-numeric entries become NaN, so each trait drops only its own
    missing individuals. A two-column file gives one trait.
    """
    bad = set(float(v) for v in miss)
    with open(path) as f:
        head = _split(f.readline())
        k = len(head) - 1
        if k < 1:
            raise SystemExit("phenotype file %s needs an IID column and at least one trait" % path)
        tab = {}
        for ln in f:
            p = _split(ln)
            if not p or not p[0]:
                continue
            v = np.full(k, np.nan)
            for t in range(k):
                if t + 1 < len(p):
                    try:
                        x = float(p[t + 1])
                    except ValueError:
                        continue
                    if x not in bad and np.isfinite(x):
                        v[t] = x
            tab[p[0]] = v
    Y = np.array([tab.get(i, np.full(k, np.nan)) for i in ids], float).reshape(len(ids), k)
    return Y, head[1:]


def read_pheno(path, ids, col=None, miss=(-9.0,)):
    """One trait from a phenotype table (column `col`, 1-based after IID; default 1)."""
    Y, _ = read_pheno_table(path, ids, miss=miss)
    return Y[:, (1 if col is None else col) - 1]


def write_matrix(path, M, ids=None, renames=None):
    np.save(path, M)
    return path
