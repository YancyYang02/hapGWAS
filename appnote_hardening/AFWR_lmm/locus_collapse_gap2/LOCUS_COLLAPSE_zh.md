# 重叠窗口收成位点 — AFWR (lmm)

- 输入：`AFWR_lmm/assoc_blocks.tsv`
- n_blocks=749, Bonf thr=6.6756e-05, max_gap=2 (SNP index)
- 合并规则：窗口 SNP 区间重叠（或间距 ≤ max_gap）则并入同一 locus；lead = min p

| method | raw hits | loci | compression |
|---|---:|---:|---:|
| joint | 42 | 28 | 1.50x |
| acj | 25 | 18 | 1.39x |
| allele_bonf | 73 | 42 | 1.74x |
