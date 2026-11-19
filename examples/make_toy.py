"""Generate the toy dataset used in the README (deterministic).

Writes, next to this script:
  toy/panel.vcf.gz   phased VCF, 400 individuals, 2 chromosomes x 300 markers
  toy/height.txt     phenotype (IID, value), tab-separated, header line
  toy/covar.txt      two numeric covariates (IID, sex, batch), tab-separated
  toy/truth.txt      the planted signal

Haplotypes are mosaics of a few founder haplotypes per 20-marker segment, so
the panel has block-like LD. The trait has a polygenic background plus one
haplotype-specific effect on chromosome 1: carrying two copies of one
founder haplotype of the segment at markers 140-159 (a recessive haplotype
effect that no single marker tags well). Run with: python make_toy.py
"""
import gzip
import os

import numpy as np

SEED = 20261117
N, M_PER_CHR, CHROMS, SEG, NFOUND = 400, 300, ("1", "2"), 20, 5


def main():
    rng = np.random.default_rng(SEED)
    here = os.path.dirname(os.path.abspath(__file__))
    out = os.path.join(here, "toy")
    os.makedirs(out, exist_ok=True)
    haps, meta = [], []
    causal = None
    for ch in CHROMS:
        H = np.zeros((2 * N, M_PER_CHR), np.int8)
        for s in range(0, M_PER_CHR, SEG):
            w = min(SEG, M_PER_CHR - s)
            founders = (rng.random((NFOUND, w)) < rng.uniform(0.2, 0.8, w)).astype(np.int8)
            freq = rng.dirichlet(np.full(NFOUND, 2.0))
            pick = rng.choice(NFOUND, size=2 * N, p=freq)
            H[:, s:s + w] = founders[pick]
            if ch == "1" and s == 140:
                # the most frequent founder that is not the commonest one
                order = np.argsort(-np.bincount(pick, minlength=NFOUND))
                causal = (pick == order[1]).astype(int)
        flip = rng.random(H.shape) < 0.01          # mutation noise
        H = np.where(flip, 1 - H, H).astype(np.int8)
        haps.append(H)
        meta += [(ch, 10000 + 2000 * j) for j in range(M_PER_CHR)]
    H = np.concatenate(haps, 1)
    ids = ["ind%03d" % i for i in range(N)]

    G = (H[0::2] + H[1::2]).astype(float)
    p = G.mean(0) / 2
    ok = (p > 0.01) & (p < 0.99)
    Z = (G[:, ok] - 2 * p[ok]) / np.sqrt(2 * p[ok] * (1 - p[ok]))
    poly = Z @ rng.normal(0, 1, Z.shape[1])
    poly = poly / poly.std() * np.sqrt(0.30)
    two = (causal[0::2] + causal[1::2] == 2).astype(float)
    hap = (two - two.mean()) / two.std() * np.sqrt(0.12)
    sex = rng.integers(0, 2, N).astype(float)
    batch = rng.normal(0, 1, N)
    y = poly + hap + 0.3 * sex + 0.2 * batch + rng.normal(0, np.sqrt(0.58), N)

    with gzip.open(os.path.join(out, "panel.vcf.gz"), "wt", compresslevel=9) as f:
        f.write("##fileformat=VCFv4.2\n")
        f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(ids) + "\n")
        for j, (ch, pos) in enumerate(meta):
            gt = ["%d|%d" % (H[2 * i, j], H[2 * i + 1, j]) for i in range(N)]
            f.write("%s\t%d\t%s_%d\tA\tG\t.\tPASS\t.\tGT\t%s\n" % (ch, pos, ch, pos, "\t".join(gt)))
    with open(os.path.join(out, "height.txt"), "w") as f:
        f.write("IID\theight\n")
        f.writelines("%s\t%.5f\n" % (i, v) for i, v in zip(ids, y))
    with open(os.path.join(out, "covar.txt"), "w") as f:
        f.write("IID\tsex\tbatch\n")
        f.writelines("%s\t%d\t%.5f\n" % (i, s, b) for i, s, b in zip(ids, sex, batch))
    with open(os.path.join(out, "truth.txt"), "w") as f:
        f.write("chrom\tstart_bp\tend_bp\tmodel\tvariance_explained\n")
        f.write("1\t%d\t%d\trecessive_haplotype\t0.12\n" % (meta[140][1], meta[159][1]))
    print("toy data written to", out, "| carriers of two copies:", int(two.sum()))


if __name__ == "__main__":
    main()
