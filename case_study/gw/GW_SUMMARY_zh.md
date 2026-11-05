# Genome-wide claw/BW hapGWAS (HPC)

- n=2351, chromosomes=31, blocks=16361, alleles=97766, snps=32802
- Bonferroni thr (0.05/n_blocks) = 3.06e-06
- Walltime total sec = 381.3

## Hits @ block-Bonferroni

| method | n_hits | min_p |
|---|---:|---:|
| ACJ | 8 | 1.31e-08 |
| joint | 8 | 2.55e-10 |
| GHap_like_bon | 34 | 2.5e-12 |
| SNP_bon | 0 | 1.24e-05 |
| PCA | 5 | 6.99e-11 |
| ridge | 2 | 6.15e-08 |

## ACJ vs GHap-like (allele-Bonf within block, then genome Bonf)

- shared: **6**
- ACJ only: **2** → 2:2:w611, 3:3:w228
- GHap only: **28** → 17:17:w67, 1:1:w1338, 1:1:w1611, 1:1:w1635, 1:1:w1965, 1:1:w1982, 1:1:w2835, 20:20:w10, 20:20:w151, 20:20:w5, 20:20:w6, 20:20:w9, 23:23:w24, 2:2:w1029, 2:2:w1830, 3:3:w1506, 3:3:w1507, 3:3:w1864, 3:3:w1865, 4:4:w1289
- ACJ hits not in SNP-Bonf: **8**

Proxy note: GHap-like = hapGWAS allele marginal min-p with within-block Bonferroni (same dosages/null as ACJ). Full R GHap not required for head-to-head on identical blocks.
