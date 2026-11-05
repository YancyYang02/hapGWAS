# hapGWAS 真实连续表型案例（HPC 全流程）

> 计算全部在 `192.168.10.50` 完成；本机仅回传小结果表/图/md，**未下载** phased VCF。

## 设计

| 项 | 值 |
|----|----|
| 性状 | **爪重相对体重残差**（`gemma_claw_per_bw.pheno`，连续） |
| n | **2351**（VCF 2384 中去掉 NA 后对齐） |
| 染色体 | **chr33**（既有翅峰线索）、**chr4**（芯片 SNPs 较多） |
| 基因型 | `.../work_hybrid/phased_chip_chr{33,4}.vcf.gz` |
| 分块 | 固定窗口 window=4, step=2 |
| Null | OLS（本案例未拉 GRM；后续可接 GCTA/GEMMA kinship） |
| 工作目录 | `/fs2/home/clwen/YancongYang/hapgwas_case/results` |

表型：mean=34.78, sd=9.93, nunique=2322（明确非有序胸型）。

## 染色体摘要

| chrom | n | snps | blocks | min joint_p | min acj_p |
|------:|--:|-----:|-------:|------------:|----------:|
| 33 | 2351 | 78 | 38 | 2.18e-2 | 8.6e-3 |
| 4 | 2351 | 3165 | 1581 | **3.76e-7** | **1.31e-8** |

## Top 块方法对比（按 ACJ p）

块内 allele/SNP 报 Bonferroni；ACJ/joint 为设计驱动 χ²。

| chrom | block | joint_p | acj_p | allele_bon | snp_bon | winner |
|-------|-------|--------:|------:|-----------:|--------:|--------|
| 4 | 4:w1290 | 3.76e-7 | **1.31e-8** | 2.01e-8 | 0.557 | **ACJ** |
| 4 | 4:w1289 | 3.77e-6 | 4.22e-6 | **8.84e-7** | 0.219 | allele-Bonf |
| 4 | 4:w1268 | 1.93e-4 | 1.68e-4 | **2.32e-5** | 0.172 | allele-Bonf |
| 4 | 4:w1343 | 6.21e-4 | **1.87e-4** | 4.58e-4 | 6.04e-3 | **ACJ** |
| 4 | 4:w178 | 6.23e-4 | **2.22e-4** | 3.35e-3 | 1.42e-2 | **ACJ** |
| 4 | 4:w599 | 9.82e-3 | **3.58e-4** | 6.30e-4 | 0.591 | **ACJ** |
| 4 | 4:w1303 | **3.94e-4** | 4.95e-4 | 4.30e-4 | 0.125 | joint |

要点：

1. **chr4:w1290**：ACJ（1.3e-8）强于 joint（3.8e-7）；块内 SNP-Bonf≈0.56 → **单 SNP 基线几乎看不到、单倍型/ACJ 才检出**的典型场景。
2. 邻近块 **w1289/w1268**：allele-Bonf 更尖 → 存在“尖锐单等位”成分，**不能只报 ACJ**。
3. chr33 在本窗口设定下仅达 suggestive（min joint≈0.022）；翅峰需更密标记或全基因组扫描再议。

## 何时不要单独用 ACJ

与 `METHOD_NOTE_zh.md` 一致：尖锐单等位优先 allele-Bonf/joint；探索性扫描必须并列 joint + ACJ + allele-Bonf。

## 产物路径

**HPC：** `/fs2/home/clwen/YancongYang/hapgwas_case/results/`  
（`assoc_blocks_all.tsv`, `assoc_alleles_all.tsv`, `assoc_snps_all.tsv`, `top_blocks_method_compare.tsv`, `fig_case_manhattan_snp_vs_hap.png`, `case_meta.json`）

**Desktop：** `C:\Users\杨颜聪\Desktop\hapgwas\case_study\`

## 基因组扩展（31 chr，HPC）

| 项 | 值 |
|----|----|
| chroms | 31（缺 29/32 等无 phased 文件） |
| blocks / alleles / SNPs | 16361 / 97766 / 32802 |
| walltime | **381.3 s** |
| Bonf thr | 3.06×10⁻⁶ |

Hits @ Bonf：ACJ **8**，joint **8**，GHap-like **34**，SNP-Bonf **0**。  
ACJ vs GHap-like：shared 6，ACJ-only 2（`2:w611`,`3:w228`），GHap-only 28。  
详表：`case_study/gw/hits_unique_shared.tsv`、`GW_SUMMARY_zh.md`。

## 第二性状：翅重 wing（绝对）

| 项 | 值 |
|----|----|
| 来源 | `gemma_final.pheno` 第 4 列（脚本标注 Wing） |
| n / mean±sd | 2351 / 238.5±33.3 g |
| Hits @ Bonf | ACJ **5** · joint **5** · GHap-like **6** · SNP **0** |
| ACJ vs GHap | shared 2 · ACJ-only 3 · GHap-only 4 |
| top ACJ | `6:w398` **6.19e-9** |

Desktop：`case_study/gw_wing/` · HPC：`hapgwas_case/gw_results_wing/`
