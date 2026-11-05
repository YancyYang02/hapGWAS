# GHap 安装记录（已解决安装；assoc 收尾中）

更新：2026-11-05 20:43 Asia/Shanghai

## 已解决
1. 旧脚本 `ganoninc/GHap` → GitHub **404**。
2. 改装 `ytutsunomiya/GHap` **4.0.0.14**，`library(GHap)` OK。
3. phase 格式：行=位点、列=2n、编码 0/1。
4. GRM：`gcta_grm.grm.id/.bin`。

## 收尾中
- `ghap.subset` variants 必须用等位基因 **名字**；修正后重跑 `ghap.assoc`（BW6、AFWR）。
- 成功后写入 `ghap_h2h/H2H_table.tsv`；在此之前主文对照列勿把 hapGWAS `allele_bonf` 写成 GHap。

## 与项4并行交付
- AFWR LMM joint 42 → **28 loci**（touch）；见 `AFWR_lmm/locus_collapse/`。
