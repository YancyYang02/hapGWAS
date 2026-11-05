# 项5：条件分析 + 基因注释（LMM）

更新时间：2026-11-05 21:29（Asia/Shanghai）

## 设定
- 数据：`phased_chr4_1p5k.vcf.gz`（n=465, p=1500），纯生产表型 BW6 / AFWR / BW3_6
- null：GCTA GRM + spectral LMM（与项1同）
- 窗口 Bonf：0.05/749 ≈ 6.68e-5；SNP-gw Bonf：0.05/1500 ≈ 3.33e-5
- 主方法：joint / ACJ；SACJ 不进主结论
- 条件：把块内最强 SNP 剂量当作协变量，对单倍型 dosage 做残差联合检验；反向把最强 hap 等位基因剂量当作协变量再测该 SNP

## 选点
| locus | 性状 | lead block | 坐标（VCF 字段） | 为何选 |
|---|---|---|---|---|
| BW6_40_L1 | BW6 | 40:w543 | chrom=`40`, 625414–627113 | LMM 下唯一 joint/ACJ hit；snp_gw hits=0 |
| AFWR_40_L1 | AFWR | 40:w743 | chrom=`40`, 699040–699111 | touch 合并后 top locus（joint_p 最小）；局部 snp_gw 未过 Bonf |

## 条件分析结果（摘要）

| locus | uncond joint_p | uncond ACJ_p | 块内最强 SNP p | cond joint_p（控 SNP） | 最强 hap 等位基因 | allele\|SNP p | SNP\|allele p | 解读 |
|---|---:|---:|---:|---:|---|---:|---:|---|
| BW6_40_L1 | 2.93e-6 | 5.08e-7 | 0.178（40:625414） | **1.99e-6** | OTHER (freq≈0.0097) | 2.53e-7 | 0.775 | **hap 存活**；SNP 几乎无信号，控 SNP 后 joint 几乎不变 |
| AFWR_40_L1 | 4.07e-16 | 1.03e-16 | 0.0117（40:699051） | **5.36e-18** | OTHER (freq≈0.0065) | 2.70e-17 | 0.444 | **hap 存活且更强**；局部最强 SNP 远未过 snp_gw Bonf |

文件：`results/conditional_summary.tsv`、`*_conditional.json`、`*_local_snps.tsv`

### 诚实边界
1. 两个 locus 的最强 hap 等位基因标签都是 **OTHER**（低频等位基因被 `freq` 塌缩合并）。这说明信号来自稀有/合并类，**不能**直接写成「某个具体四位点单倍型序列」。投稿需写明 collapse 规则与 OTHER 频率。
2. AFWR 的 cond_joint_p 比 uncond 更小，符合「SNP 解释不了联合信号」；也提示多自由度/稀有类在小样本上可能偏乐观——应用项1的 λ 与亲缘校正一起报。
3. 仍是 **chr 子集（~1500 SNP）**，不是全基因组。

## 相关性状（BW3_6 @ BW6 lead）
| 检验 | 数值 | 是否过块 Bonf |
|---|---:|---|
| BW3_6 joint @ BW6_40_L1 | 1.17e-4 | 否（阈 6.68e-5） |
| BW3_6 ACJ @ BW6_40_L1 | **2.19e-5** | **是** |
| BW3_6 局部最强 SNP | 0.105 | 否 |
| 控局部 SNP 后 joint | 2.15e-6 | — |

→ 体重增量性状在同一区间有 ACJ 信号，方向与 BW6 一致，支持「非孤立假阳性」的弱旁证。  
交叉：AFWR 在 BW6 位点、BW6/BW3_6 在 AFWR 位点均无信号（`cross_trait_at_leads.tsv`）。

## 基因注释
- VCF `#CHROM` 字面量是 **`40`**（无 `##contig`）；文件名 `phased_chr4_1p5k` → **按鸡 chr4 解释**。Ensembl 当前默认组装 **bGalGal1.mat.broiler.GRCg7b** 无 chrom 40；chrom=40 查询 400。
- 注释：把坐标映射到 **chr4**，用 Ensembl REST overlap（中点 ±250 kb）。

| locus | 解释坐标 | 最近蛋白编码基因 | 距离 | 关系（相对 locus） | 其他邻近 |
|---|---|---|---:|---|---|
| BW6_40_L1 | 4:625414–627113 | **RXFP2** (ENSGALG00010015843) | ~57.4 kb | 基因在更低坐标侧 | STARD8 ~70 kb；YIPF6；OPHN1；AR |
| AFWR_40_L1 | 4:699040–699111 | **RXFP2** | ~130 kb | 同上 | STARD8 ~143 kb；区间内无基因重叠 |

**注意**：这是 GRCg7b 上的就近注释；若芯片坐标其实是 Galgal6/GRCg6a，基因距离可能偏移。投稿应写清「VCF 标 40、按 chr4+GRCg7b 解释」，并在全基因组重跑后用正式 build 再钉一次。

生物学一句话（谨慎）：RXFP2（relaxin-family peptide receptor 2）在哺乳动物繁殖/隐睾相关，禽类生长性状关联需文献核对，**本轮只报到基因符号与距离，不做过度机制推断**。

文件：`gene_annotation.json`

## 项3补强：snp_gw（Bonf by n_SNP）主表

| trait | null | n_SNP | thr=0.05/n | hits | min_p | min SNP |
|---|---|---:|---:|---:|---:|---|
| BW6 | **lmm** | 1500 | 3.33e-5 | **0** | 1.69e-3 | 40:360033 |
| BW6 | ols | 1500 | 3.33e-5 | 0 | 5.66e-4 | 40:374236 |
| AFWR | **lmm** | 1500 | 3.33e-5 | **1** | 1.80e-5 | **40:345244**（≠ AFWR hap locus） |
| AFWR | ols | 1500 | 3.33e-5 | 1 | 2.83e-5 | 40:345244 |
| BW3_6 | **lmm** | 1500 | 3.33e-5 | 0 | 1.29e-3 | 40:374236 |
| EBW6 | **lmm** | 1500 | 3.33e-5 | **0** | 1.86e-3 | 40:388244 |
| EBW6 | ols | 1500 | 3.33e-5 | **8** | 5.18e-6 | （亲缘未校正假阳性） |
| GBW6 | **lmm** | 1500 | 3.33e-5 | 0 | 3.97e-2 | 40:369537 |
| ERFI | **lmm** | 1500 | 3.33e-5 | 0 | 3.33e-3 | 40:127643 |

要点：AFWR 唯一 snp_gw hit 在 ~345 kb，与 hap lead（~699 kb）不是同一位点；BW6/BW3_6 在 LMM 下 snp_gw=0，而 haplotype 有 hit。

文件：`results/snp_gw_summary.tsv`、`*_lmm_snp_gw.tsv`、`snp_gw_hits_*.tsv`

## 产物路径
- HPC：`/fs5/home/clwen/YancongYang/hapgwas_try/runs/lmm_grm/item5_conditional/`
- 本机：`Desktop/hapgwas/appnote_hardening/item5_conditional/`
- 脚本：`run_conditional_loci.py`（venv + `PYTHONPATH=src_new`）
