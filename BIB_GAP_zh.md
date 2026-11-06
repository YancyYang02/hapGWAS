# Bioinformatics Application Note 差距清单（相对 BIB App Note）

**0.7.1：** 包版本已到 0.7.1。默认排序是 `omni_p`；`beyond_p` / `beyond_lrt_p` 与 `tree_p` 是次要列；CLI 有 `--lowo`、`--kin-geno`、`--gcta-grm` 和精确 LRT 列。方法口径以 `README.md` 为准，不把 LOCO / RHM 写成新方法，也不把 SACJ 当主结果。下面 Availability 里的 DOI 数字仍是 **v0.4.1** 的 Zenodo 记录，重新存档之前不要改。



对照典型 *Bioinformatics* Application Note：短文、软件可用、1–2 真实数据集、与现有工具比较、清晰 novelty。

| 要件 | 状态 | 证据 / 缺口 |
|------|------|-------------|
| **Novelty claim** | ✅ 诚实 | ACJ + 同投影多 score；明确不虚报 Gabriel/SKAT/LMM 发明 |
| **Availability** | ✅ | MIT + CITATION + README；https://github.com/YancyYang02/hapGWAS （`v0.4.1`）；版本 DOI [10.5281/zenodo.23984394](https://doi.org/10.5281/zenodo.23984394)；概念 DOI [10.5281/zenodo.23984393](https://doi.org/10.5281/zenodo.23984393) |
| **Docs** | ✅ | `APPLICATION_NOTE_en.md` 近投稿英文稿；`METHOD_NOTE_zh.md`；`README.md` |
| **Simulation** | ✅ | Type I / 单等位 / 双近邻数字写入 App Note |
| **真实数据 #1** | ✅ | claw/BW，n=2351，31 chr |
| **真实数据 #2** | ✅ | wing，同队列 GW |
| **与现有工具比较** | ✅ proxy | GHap-like allele-Bonf 同块对照；两性状 unique/shared |
| **Runtime** | ✅ | claw ≈381 s；wing ≈383 s |
| **英文稿** | ✅ 近投稿 | `APPLICATION_NOTE_en.md`（~2k words；正式排版/图仍缺；软件 DOI 10.5281/zenodo.23984394） |
| **参考文献** | ⚠️ | 仅 `[REF:…]` 占位，**未编造** GHap/GEMMA/Beagle/FaST-LMM 条目 |
| **Web server** | 不需要 | 可选 |

## 两性状基因组 hits（Bonferroni 0.05/16361）

| 性状 | ACJ | joint | GHap-like | SNP | ∩ | ACJ-only | GHap-only | top ACJ |
|------|----:|------:|----------:|----:|--:|---------:|----------:|---------|
| claw/BW | 8 | 8 | 34 | 0 | 6 | 2 | 28 | `4:w1290` 1.31e-8 |
| wing | 5 | 5 | 6 | 0 | 2 | 3 | 4 | `6:w398` 6.19e-9 |

## 投稿前仍缺

0. （文档已做）主表口径 joint/ACJ/GHap/SNP；**SACJ 降级**；硬化数字见 `appnote_hardening/HARDENING_STATUS_zh.md`  
1. 公共仓库与 Zenodo 已写入 Availability：https://github.com/YancyYang02/hapGWAS ；版本 DOI https://doi.org/10.5281/zenodo.23984394 ；概念 DOI https://doi.org/10.5281/zenodo.23984393  
2. 将 `[REF:…]` 换成核对过的正式文献  
3. 期刊双栏 2 页排版 + 1–2 主图（sim tradeoff + GW Manhattan）  
4.（可选）官方 R GHap 同块复测以加强 “vs GHap” 表述  
