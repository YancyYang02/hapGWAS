#!/usr/bin/env python3
"""Extract chip LD windows from D8_qc BED (run on user Windows machine)."""
import json
from pathlib import Path
from collections import defaultdict
import numpy as np

BFILE = Path(r"D:\12性状\D8_qc")
OUT = Path(r"C:\Users\杨颜聪\Desktop\hapgwas\chip_ld\d8_ld_windows.json")


def main():
    bim = []
    with open(str(BFILE) + ".bim", encoding="utf-8", errors="ignore") as f:
        for i, line in enumerate(f):
            chrom, snp, _, bp, a1, a2 = line.split()
            bim.append((chrom, snp, int(bp), i))
    fam_n = sum(1 for _ in open(str(BFILE) + ".fam", encoding="utf-8", errors="ignore"))
    print("n_indiv", fam_n, "n_snp", len(bim), flush=True)

    byc = defaultdict(list)
    for chrom, snp, bp, i in bim:
        byc[chrom].append((bp, i, snp))

    bytes_per_snp = (fam_n + 3) // 4

    def read_snps(indices):
        X = np.zeros((fam_n, len(indices)), dtype=np.float32)
        with open(str(BFILE) + ".bed", "rb") as f:
            magic = f.read(3)
            assert magic == b"\x6c\x1b\x01"
            for col, idx in enumerate(indices):
                f.seek(3 + idx * bytes_per_snp)
                raw = f.read(bytes_per_snp)
                codes = np.frombuffer(raw, dtype=np.uint8)
                g = np.full(fam_n, np.nan, dtype=np.float32)
                for j in range(fam_n):
                    byte = int(codes[j // 4])
                    code = (byte >> (2 * (j % 4))) & 0b11
                    if code == 0:
                        g[j] = 0.0
                    elif code == 2:
                        g[j] = 1.0
                    elif code == 3:
                        g[j] = 2.0
                X[:, col] = g
        return X

    windows = []
    rng = np.random.default_rng(7)
    for chrom in ["1", "2", "3", "4", "5", "Z"]:
        arr = byc.get(chrom)
        if not arr or len(arr) < 20:
            continue
        arr = sorted(arr)
        n_starts = min(4, max(1, len(arr) // 50))
        starts = rng.choice(max(1, len(arr) - 8), size=n_starts, replace=False)
        for s in starts:
            s = int(s)
            block = arr[s : s + 5]
            if len(block) < 4:
                continue
            block2 = [block[0]]
            for b in block[1:]:
                if b[0] - block2[0][0] <= 100_000 and len(block2) < 5:
                    block2.append(b)
            block = block2
            if len(block) < 4:
                continue
            idxs = [b[1] for b in block]
            X = read_snps(idxs)
            ok = np.all(np.isfinite(X), axis=1)
            X = X[ok]
            if X.shape[0] < 200:
                continue
            maf = X.mean(0) / 2.0
            maf = np.minimum(maf, 1 - maf)
            Xc = X - X.mean(0)
            sd = Xc.std(0)
            sd[sd < 1e-8] = np.nan
            R = (Xc.T @ Xc) / (X.shape[0] * np.outer(sd, sd))
            win = {
                "chrom": chrom,
                "snps": [b[2] for b in block],
                "bp": [int(b[0]) for b in block],
                "n": int(X.shape[0]),
                "maf": [float(x) for x in maf],
                "r": [
                    [float(R[i, j]) if np.isfinite(R[i, j]) else 0.0 for j in range(R.shape[1])]
                    for i in range(R.shape[0])
                ],
                "mean_abs_r_adj": float(np.nanmean(np.abs(np.diag(R, k=1)))),
                "mean_abs_r": float(np.nanmean(np.abs(R[np.triu_indices_from(R, 1)]))),
            }
            windows.append(win)
            print("win", chrom, block[0][0], len(block), "mean|r_adj|", win["mean_abs_r_adj"], flush=True)

    OUT.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "source": str(BFILE),
        "n_indiv": fam_n,
        "n_snp": len(bim),
        "n_windows": len(windows),
        "windows": windows,
        "note": "Unphased chip LD from D8_qc; phased DoubleMuscling VCF not on this machine",
    }
    OUT.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print("wrote", OUT, "n_windows", len(windows), flush=True)


if __name__ == "__main__":
    main()
