# 重叠窗口收成位点 — AFWR（LMM）

- n_blocks=749, Bonf thr=6.6756e-05
- 窗口 4/2；`snp_end` 为开区间
- **主口径 touch**：区间重叠或首尾相接即合并（覆盖滑动窗相邻 hit）
- 另报 overlap（仅严格重叠）与 gap2（允许间隔 ≤2 个 SNP）

| mode | joint raw→loci | acj | allele_bonf |
|---|---|---|---|
| overlap | 42→33 | 25→20 | 73→47 |
| touch | 42→28 | 25→18 | 73→42 |
| gap2 | 42→28 | 25→18 | 73→42 |

## 解读
- AFWR LMM 下 joint 42 个窗口 hit **不是**单一 QTL 被步长撕成几十个；touch 合并后仍有 **28** 个分散 locus。
- 说明项 4 必要，但合并后仍需项 5（条件分析）和更强多重检验/全基因组验证，不能直接当 28 个生物学位点。
- BW6：1 hit → 1 locus（无需合并）。

## 文件
- `loci_<method>.tsv`：主口径 touch
- `loci_<method>_<mode>.tsv`：三种合并规则
- `raw_hits_*.tsv`：合并前窗口
