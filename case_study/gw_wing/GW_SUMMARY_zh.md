# Genome-wide wing hapGWAS (HPC)

- trait=`wing` (gemma_final col 4), n=2351, chroms=31, blocks=16361, alleles=97766, snps=32802
- Bonferroni thr = 3.06e-06; walltime = 383.4s

## Hits @ block-Bonferroni

| method | n_hits | min_p |
|---|---:|---:|
| ACJ | 5 | 6.19e-09 |
| joint | 5 | 1.81e-08 |
| GHap_like_bon | 6 | 2.64e-08 |
| SNP_bon | 0 | 8.27e-06 |
| PCA | 0 | 4.15e-06 |
| ridge | 2 | 3.4e-07 |

## ACJ vs GHap-like

- shared: **2**
- ACJ only: **3** → 1:1:w2323, 4:4:w105, 6:6:w398
- GHap only: **4** → 1:1:w2328, 2:2:w399, 2:2:w401, 2:2:w414
- ACJ not in SNP-Bonf: **5**
