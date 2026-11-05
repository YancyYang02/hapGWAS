# Bioinformatics Application Note 硬证据清单（进行中）

更新时间：2026-11-05 21:29（Asia/Shanghai）
目标：`Bioinformatics` Application Note（先钉死证据；BIB 另议）

## 总原则
- 纯生产/屠体表型（GPE 61 性状），**不含代谢物**。
- 主表四列：**joint / ACJ / allele（正式 GHap）/ SNP**；**SACJ 降级**，不作主结果。
- 肉鸡 + D8 **禁止只报 OLS**；必须有亲缘校正 LMM 与基因组膨胀系数 λ。

## 清单状态

| # | 项 | 状态 | 备注 |
|---|---|---|---|
| 1 | GRM-LMM + λ | **完成（chr4≈1500 SNP 子集）** | EBW6 OLS 膨胀后 LMM 清零；AFWR LMM 仍多 hap hit |
| 2 | 正式 GHap 对打 | **完成（chr4≈1500）** | GHap **4.0.0.14**；BW6 allele/block=11/20；AFWR allele/block=98/117；PD-regularized 同套 GCTA GRM；见 `ghap_h2h/H2H_table.tsv` |
| 3 | 全位点 SNP Bonf（按 #SNP） | **完成（关键性状已制表）** | LMM 下 BW6/BW3_6/EBW6/GBW6/ERFI=0 hit；AFWR=1（40:345244）；见 `item5_conditional/results/snp_gw_summary.tsv` |
| 4 | 重叠窗口收成位点 | **本轮完成（AFWR/BW6 LMM）** | 见下表；主口径 **touch** |
| 5 | hap>SNP 条件分析+基因 | **本轮完成（2 loci）** | BW6/AFWR lead 控最强 SNP 后 joint 仍显著；映射 chr4≈RXFP2；见 `item5_conditional/ITEM5_zh.md` |
| 6 | 公开仓库版本 + Zenodo DOI | 已完成 | https://github.com/YancyYang02/hapGWAS ，代码 **0.4.1**（`--grm`）；版本 DOI 10.5281/zenodo.23984394 |

## 软件改动
- `hapgwas` **0.4.1**：`--grm`/`--kinship` → spectral LMM（GCTA `.grm.bin+.id` / npy / npz / 文本）。
- HPC：`/fs5/home/clwen/YancongYang/hapgwas_try`；结果 `runs/lmm_grm/`、`runs/ghap_h2h/`。
- 本机：`Desktop/hapgwas/appnote_hardening/`。

## 1) GRM-LMM + λ

### 设定
- 基因型：`phased_chr4_1p5k.vcf.gz`（n=465, p≈1500）；表型 n=463
- GRM：`Genetic_Parameter_Estimation/gcta_grm`（与 `a*` ID 全对齐；diag_mean≈1.079）
- 窗口 4/2 → 749 blocks；块 Bonf≈6.68e-5；SNP-gw Bonf=`0.05/1500≈3.33e-5`
- λ：median χ²₁(p) / 0.4549362344
- **口径**：单染色体子集上的 λ ≠ 全基因组 λ；本轮用于证明亲缘校正必要并打通流水线。

### 结果表（精确值来自 `key_traits_lmm_vs_ols.json`）

| trait | null | vg | ve | λ_joint | λ_acj | λ_allele | λ_snp_gw | hits j/acj/allele | snp_gw hits |
|---|---|---:|---:|---:|---:|---:|---:|---|---:|
| BW6 | **lmm** | 2.703e+04 | 1.047e+05 | 1.159 | 1.561 | 0.625 | 1.851 | 1/1/1 | 0 |
| BW6 | **ols** | 0 | 1.341e+05 | 1.089 | 1.544 | 0.565 | 2.189 | 1/1/1 | 0 |
| AFWR | **lmm** | 0.01785 | 0.1562 | 0.005 | 0.035 | 0.000 | 0.834 | 42/25/73 | 1 |
| AFWR | **ols** | 0 | 0.1754 | 0.016 | 0.072 | 0.000 | 1.275 | 41/24/76 | 1 |
| BW3_6 | **lmm** | 1.807e+04 | 7.002e+04 | 1.049 | 1.403 | 0.451 | 2.153 | 0/1/1 | 0 |
| BW3_6 | **ols** | 0 | 8.974e+04 | 1.290 | 1.635 | 0.553 | 2.924 | 0/2/1 | 0 |
| EBW6 | **lmm** | 3066 | 0.03066 | 1.124 | 1.398 | 0.573 | 2.399 | 0/0/0 | 0 |
| EBW6 | **ols** | 0 | 4182 | 2.214 | 2.129 | 1.715 | 2.447 | 1/3/0 | 8 |
| GBW6 | **lmm** | 2410 | 0.0241 | 0.241 | 0.299 | 0.000 | 0.350 | 0/0/0 | 0 |
| GBW6 | **ols** | 0 | 3688 | 1.028 | 0.872 | 0.302 | 1.171 | 0/0/0 | 0 |
| ERFI | **lmm** | 2.375 | 0.1284 | 0.311 | 0.399 | 0.000 | 0.517 | 0/0/0 | 0 |
| ERFI | **ols** | 0 | 3.699 | 1.062 | 0.937 | 0.415 | 1.057 | 0/0/0 | 0 |

### 解读
1. **EBW6（最强证据）**：OLS λ_joint=2.214 且 j/acj/snp_gw=1/3/8；LMM 后 λ_joint=1.124，**全部 hit 清零** → OLS-only 会把亲缘结构当关联。
2. **BW6 / BW3_6**：LMM 降低 SNP-gw λ（仍>1，子集局限）；haplotype 保留少量 hit，SNP-gw Bonf=0。
3. **AFWR**：LMM 后 haplotype 仍大量 hit（42/25/73），SNP-gw 仅 1；λ≪1 与尾部 hit 并存 → **已做项4 overlapping-window 合并**（见下）；生物学解释仍待项5。
4. 投稿主表应报 **LMM**；OLS 仅作膨胀对照。SACJ 不进主表。

## 2) 正式 GHap（**完成**）
- **包**：ytutsunomiya/GHap **4.0.0.14**；`R_LIBS_USER=/fs5/home/clwen/YancongYang/hapgwas_try/R_libs`。
- **设计**：同一批块（window=4/slide=2 → 749 块 / 8449 等位基因）、同一 GCTA GRM、`ghap.assoc` LMM（`pheno ~ 1 + (1|id)` + covmat）。
- **subset**：haplo 的 `ids`/`variants` 在 `index=TRUE` 下必须用 **numeric indices**（等位基因名跨块不唯一，`index=FALSE` 会失败）。
- **GRM PD**：GHap/cholmod 对非 PD 的 GCTA GRM 会报 `Unable to solve LHS`；本轮对 K 做对称化、特征值下截断 1e-4、对角 +1e-3。与 hapGWAS 同一源文件，求解器不同。
- **结果摘要**：

| trait | hapGWAS joint/ACJ/allele/snp_gw | GHap allele hits | GHap block-minp hits | GHap min_p |
|---|---|---:|---:|---:|
| BW6 | 1/1/1/0 | 11 | 20 | 1.0040496155035e-11 |
| AFWR | 42/25/73/1 | 98 | 117 | 1.11223335667036e-103 |

- 产物：`runs/ghap_h2h/`；本机 `Desktop/hapgwas/appnote_hardening/ghap_h2h/`（`H2H_table.tsv`、`H2H_zh.md`、`*_ghap_assoc.tsv`）。
- **禁止**再用流水线 `allele_bonf` 冒充 GHap。
- AFWR 正式 GHap 等位基因 hit 很多（98），且 min_p 极小；需与 hapGWAS LMM 的 locus collapse / 条件分析一起解读，不能直接当 98 个独立 QTL。

## 3) 全位点 SNP Bonf（按 #SNP）— 已制表

- 路径：`runs/lmm_grm/<TRAIT>_{lmm,ols}/snp_gw.tsv`；汇总 `item5_conditional/results/snp_gw_summary.tsv`。
- 阈值 **0.05 / n_SNP**（本子集 n=1500 → 3.33e-5），与窗口内 `snp_bonf`（按块）分列。

| trait | null | hits_snp_gw | min_p | min SNP (chrom:pos) |
|---|---|---:|---:|---|
| BW6 | lmm | **0** | 1.69e-3 | 40:360033 |
| AFWR | lmm | **1** | 1.80e-5 | **40:345244** |
| BW3_6 | lmm | 0 | 1.29e-3 | 40:374236 |
| EBW6 | lmm | **0** | 1.86e-3 | 40:388244 |
| EBW6 | ols | **8** | 5.18e-6 | （对照：无亲缘时虚高） |
| GBW6 | lmm | 0 | 3.97e-2 | 40:369537 |
| ERFI | lmm | 0 | 3.33e-3 | 40:127643 |

解读：AFWR 的唯一 snp_gw hit **不在** haplotype lead（~699 kb）附近；BW6 haplotype hit 对应的 snp_gw=0。不能把「窗口内 snp_bonf」和「全位点 snp_gw」混用。

## 4) 重叠窗口收成位点（完成）

设定：window=4 / step=2；`snp_end` 为开区间；Bonf = 0.05/749 ≈ 6.68e-5。

| 性状 | null | method | raw window hits | loci (touch) | loci (overlap) | loci (gap2) |
|---|---|---|---:|---:|---:|---:|
| AFWR | LMM | joint | 42 | **28** | 33 | 28 |
| AFWR | LMM | ACJ | 25 | **18** | 20 | 18 |
| AFWR | LMM | allele_bonf | 73 | **42** | 47 | 42 |
| BW6 | LMM | joint/ACJ/allele | 1 | **1** | 1 | 1 |

- **主口径 touch**：区间重叠或首尾相接即合并（覆盖滑动窗相邻 hit）。
- 关键发现：AFWR 的 42 个 joint hit **不是**「一个 QTL 被步长撕成几十个」——相邻双 hit 仅约 9 对；touch 后仍有 **28** 个分散 locus。项4 必要，但合并后仍需条件分析/全基因组，不能直接当 28 个生物学 QTL。
- 产物：`AFWR_lmm/locus_collapse/`（`loci_*.tsv`、`LOCUS_COLLAPSE_zh.md`、`collapse_summary.json`）；脚本 `collapse_loci.py`。

## 5) 条件分析 + 基因注释（完成）

细节见 `item5_conditional/ITEM5_zh.md`。一句话：

| locus | uncond joint | 块内最强 SNP p | cond joint（控 SNP） | 最近基因（按 chr4 / GRCg7b） |
|---|---:|---:|---:|---|
| BW6 `40:w543` (625414–627113) | 2.93e-6 | 0.178 | **1.99e-6（存活）** | RXFP2 ~57 kb |
| AFWR `40:w743` (699040–699111) | 4.07e-16 | 0.012 | **5.36e-18（存活）** | RXFP2 ~130 kb |

- BW3_6 在 BW6 lead：ACJ 过块 Bonf（2.19e-5），joint 临界未过。
- VCF `#CHROM=40`，文件名 chr4 → 按 **chr4** 注释；组装 Ensembl 当前 **GRCg7b**。
- 最强 hap 等位基因均为 **OTHER**（低频塌缩类）——主文必须写明。

产物：`item5_conditional/results/`（本机 Desktop 同步目录同名）。

## 下一步
1. 用 `phased_first5000.vcf.gz` 或全 chr4 / 全基因组 + 同 GRM 扩项1/3/5。
2. 核对芯片原始坐标对应的正式组装（Galgal6 vs GRCg7b），必要时 liftOver 后再钉基因。
3. D8 GRM-LMM + λ。
4. 公开仓库 https://github.com/YancyYang02/hapGWAS ，VERSION 0.4.1；版本 DOI https://doi.org/10.5281/zenodo.23984394 。
5. Application Note 初稿按 LMM 主表重写（joint/ACJ/GHap allele/snp_gw），去掉 OLS 主结果与 SACJ 主表。
