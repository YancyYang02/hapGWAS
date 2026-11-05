# 重叠窗口收成位点 — BW6 (lmm)

- 输入：`BW6_lmm/assoc_blocks.tsv`
- n_blocks=749, Bonf thr=6.6756e-05, max_gap=0 (SNP index)
- 合并规则：窗口 SNP 区间重叠（或间距 ≤ max_gap）则并入同一 locus；lead = min p

| method | raw hits | loci | compression |
|---|---:|---:|---:|
| joint | 1 | 1 | 1.00x |
| acj | 1 | 1 | 1.00x |
| allele_bonf | 1 | 1 | 1.00x |
