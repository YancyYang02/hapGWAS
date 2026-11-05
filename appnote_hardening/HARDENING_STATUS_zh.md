# App Note 硬化状态摘要（可粘贴 Methods/Results）

更新时间：2026-11-05 21:40（Asia/Shanghai）  
范围：肉鸡 GPE 纯生产表型；**chr4≈1500 SNP 子集**（`phased_chr4_1p5k.vcf.gz`，n_geno=465，n_pheno=463）；窗 4/2 → 749 blocks；GCTA GRM + spectral LMM。  
主表口径：**joint / ACJ / 正式 GHap allele / SNP**。**SACJ 仅附录/探索**，不作主宣称。  
产物根目录：`Desktop/hapgwas/appnote_hardening/`。

## 1. GRM-LMM 与 λ（项1，完成）

- GRM：`Genetic_Parameter_Estimation/gcta_grm`（与 `a*` ID 对齐）。
- λ = median χ²₁(p) / 0.4549362344；**单染色体子集 λ ≠ 全基因组 λ**。
- 投稿主结果报 **LMM**；OLS 仅作膨胀对照。

| trait | null | λ_joint | λ_acj | λ_allele | λ_snp_gw | hits j/acj/allele | snp_gw |
|---|---|---:|---:|---:|---:|---|---:|
| BW6 | LMM | 1.159 | 1.561 | 0.625 | 1.851 | 1/1/1 | 0 |
| AFWR | LMM | 0.005 | 0.035 | 0.000 | 0.834 | 42/25/73 | 1 |
| EBW6 | LMM | 1.124 | 1.398 | 0.573 | 2.399 | 0/0/0 | 0 |
| EBW6 | OLS | 2.214 | 2.129 | 1.715 | 2.447 | 1/3/0 | 8 |

要点：EBW6 在 OLS 下严重膨胀且有 hit，LMM 后 **全部清零** → 禁止只报 OLS。

## 2. 正式 GHap 头对头（项2，完成）

- 包：ytutsunomiya/GHap **4.0.0.14**；`ghap.assoc` LMM；同一批块与同套 GCTA GRM（GHap 侧 PD 正则：eig clip 1e-4 + diag 1e-3）。
- **禁止**把流水线 `allele_bonf` 写成正式 GHap。

| trait | hapGWAS joint/ACJ/allele/snp_gw | GHap allele hits | GHap block-minp hits | GHap min_p |
|---|---|---:|---:|---:|
| BW6 | 1/1/1/0 | 11 | 20 | 1.00e-11 |
| AFWR | 42/25/73/1 | 98 | 117 | 1.11e-103 |

解读：**本子集上正式 GHap 等位基因/块命中往往多于 hapGWAS**；不能写成 hapGWAS 全面压过 GHap。GHap 多 hit 也不等于独立 QTL。

## 3. 全位点 SNP Bonf（按 #SNP）（项3，完成）

阈值 `0.05/1500 ≈ 3.33e-5`（与窗口内 `snp_bonf` 分列）。

| trait | null | snp_gw hits | min_p | lead |
|---|---|---:|---:|---|
| BW6 | LMM | **0** | 1.69e-3 | 40:360033 |
| AFWR | LMM | **1** | 1.80e-5 | **40:345244** |
| EBW6 | LMM | **0** | 1.86e-3 | 40:388244 |

AFWR 唯一 snp_gw 约在 345 kb，**不在** haplotype lead（~699 kb）。

## 4. 重叠窗口收成位点（项4，完成）

主口径 **touch**（区间重叠或首尾相接合并）。

| 性状 | method | raw window hits | loci (touch) |
|---|---|---:|---:|
| AFWR LMM | joint | 42 | **28** |
| AFWR LMM | ACJ | 25 | 18 |
| AFWR LMM | allele_bonf | 73 | 42 |
| BW6 LMM | joint/ACJ/allele | 1 | 1 |

42→28 后仍分散，不能直接当 28 个生物学 QTL。

## 5. 条件分析 + 基因（项5，完成）

VCF `#CHROM=40`，文件名 chr4 → 按 **chr4 / GRCg7b** 注释。

| locus | uncond joint | 块内最强 SNP p | cond joint（控 SNP） | 最近编码基因 |
|---|---:|---:|---:|---|
| BW6 `40:w543` (625414–627113) | 2.93e-6 | 0.178 | **1.99e-6（存活）** | **RXFP2** ~57 kb |
| AFWR `40:w743` (699040–699111) | 4.07e-16 | 0.012 | **5.36e-18（存活）** | **RXFP2** ~130 kb |

两处 top hap 等位基因均为 **OTHER**（低频塌缩类）——主文必须披露，不可写成具体四位点序列。

## 6. 公开仓库 + Zenodo（项6：仓库已公开；DOI 仍待 Zenodo）

- 本地版本已对齐：**0.4.1**（`VERSION` / `pyproject.toml` / `__version__` / `CITATION.cff` / README 标题 / `CHANGELOG.md`；含 `--grm`/`--kinship`）。
- 公开仓库：**https://github.com/YancyYang02/hapGWAS**（tag `v0.4.1`）。
- Zenodo DOI **未申请、未编造**。需所有者在 Zenodo 启用 GitHub↔Zenodo，并归档 `v0.4.1` GitHub Release 后，才把正式 DOI 回填到 `CITATION.cff` 与 Availability。
- **完整勾选清单**见 `ITEM6_RELEASE_zh.md`。

## 披露要点（主文必写）

1. 本子集上 **正式 GHap 常比 hapGWAS 报更多 allele/block hit**（BW6 11/20 vs 1；AFWR 98/117 vs 42/25/73/1）。
2. Lead hap 等位基因为 **OTHER**；条件分析后 joint 仍显著，但不等于已解析具体单倍型序列。
3. 两 lead 最近蛋白编码基因均为 **RXFP2**（距离约 57 kb / 130 kb）；邻近非因果不能排除。
4. AFWR genome snp_gw 唯一 hit 在 ~345 kb，与 ~699 kb haplotype lead 分离。
5. 结果基于 **~1500 SNP chr4 子集**，不是全基因组；λ 口径同限。
6. **SACJ 降级为探索/附录**；主表四列 joint / ACJ / GHap allele / SNP。

## 建议粘贴句（英文 Methods 骨架，可改写）

> Association used a spectral LMM with the GCTA GRM. Primary haplotype scores were the full joint and ACJ tests; allele-wise GHap (v4.0.0.14) and SNP genome-wide Bonferroni (0.05/n_SNP) were co-reported. SACJ was retained in software output but treated as exploratory. On a chr4 ~1,500-SNP subset, conditioning on the strongest local SNP left BW6 and AFWR lead haplotype joints significant (1.99×10⁻⁶ and 5.36×10⁻¹⁸). Nearest coding gene (GRCg7b chr4) was RXFP2; top alleles were OTHER collapses. Official GHap yielded more allele/block hits than hapGWAS on this subset.
