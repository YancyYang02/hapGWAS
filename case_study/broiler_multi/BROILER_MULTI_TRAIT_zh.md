# Broiler multi-trait hapGWAS（纯生产/屠体表型）

## 更正说明
- **上一版错误**：用了 9 个生产性状 + top-40 高 h² **代谢物**（`phenotype2.txt` / `Id#####` 等），不符合课题设定。
- **本版更正**：仅使用纯物理/生产/屠体表型面板（体重、增重、采食与饲料效率、胸腿翅爪器官重量与比率、肌纤维等），**不含任何代谢物/脂质**。

## 一句话汇总
- 队列：Genetic_Parameter_Estimation / zhouyan broiler **n=463**（fam 463），表型 **61** 个纯生产性状。
- 表型来源：`/fs2/home/broiler/YancongYang/hapgwas_broiler_multi/gcta_pheno_final.txt`（列序 = `run_all_traits.sh` 的 `ALL_TRAITS`）。
- 基因型：phased VCF 前 5000 位点子集（`2.Filter_SNPs.vcf`），与 gemma_463 同队列 IID。
- 工具：`hapgwas run-multi --methods all` → joint / ACJ / SACJ / allele_bonf / snp_bonf / pca / ridge。
- **SACJ 命中数 > SNP-bonf：6/61 性状**（`BW3_6`, `LWR`, `AFWR`, `EBW6`, `GBW6`, `ERFI`）；ACJ 同指标：**5/61**（同上除 `EBW6`）。
- 有任一方法命中的性状共 **12/61**：BW3, BW3_6, BW3_zq, RFI2, LWR, HWR, AFWR, WCLM, PM_CSA, EBW6, GBW6, ERFI。
- 位点子集上 SACJ 总 hits=25，SNP-bonf 总 hits=22（窗口前 5000 位点、2499 blocks，Bonf thr≈2e-5）。

## 数据与方法
- 原始表型目录：`/fs2/home/clwen/YancongYang/Genetic_Parameter_Estimation/`（`gcta_pheno_final.txt` / `phe_gemma_all.txt` / `run_all_traits.sh`）
- cohort geno fam: `/fs2/home/broiler/zhouyan/GWAS/1.1Filter_SNPs.fam`（IID=`a1`…，与 phased VCF / gcta_pheno 一致；**不用** `gemma_463_gwas_updated.fam` 的系谱号）
- pheno staged copy → `pheno_pure61.tsv`，**n_traits=61**（纯生产/屠体，无代谢物；含 BFCf/LFCf 等肌肉成分，属 ALL_TRAITS 物理表型）
- geno subset: `phased_first5000.vcf.gz` (first 5000 phased sites)
- methods: joint / ACJ / SACJ / allele_bonf (GHap-like) / snp_bonf / pca / ridge
- CLI: `hapgwas run-multi --methods all --block-mode window --window 4 --step 2`

## 61 个性状列表（ALL_TRAITS）

| # | trait | 中文/含义 |
|---:|---|---|
| 1 | `BW3` | 21日龄体重 |
| 2 | `BW6` | 42日龄体重 |
| 3 | `BW3_6` | 3-6周增重 |
| 4 | `BW` | 宰前体重 |
| 5 | `BW3_zq` | 3周-宰前增重 |
| 6 | `BW6_zq` | 6周-宰前增重 |
| 7 | `FI` | 采食量 |
| 8 | `AFI` | 平均日采食量 |
| 9 | `RFI1` | 剩余采食量1 |
| 10 | `RFI2` | 剩余采食量2 |
| 11 | `RFI3` | 剩余采食量3 |
| 12 | `RFI4` | 剩余采食量4 |
| 13 | `FRI` | 饲料转化率FCR |
| 14 | `PMA` | 胸大肌重 |
| 15 | `PMAR` | 胸大肌率 |
| 16 | `PMI` | 胸小肌重 |
| 17 | `PMIR` | 胸小肌率 |
| 18 | `PM` | 胸肌重 |
| 19 | `PMR` | 胸肌率 |
| 20 | `LW` | 肝脏重 |
| 21 | `LWR` | 肝脏比 |
| 22 | `HW` | 心脏重 |
| 23 | `HWR` | 心脏比 |
| 24 | `AFW` | 腹脂重 |
| 25 | `AFWR` | 腹脂率 |
| 26 | `MW` | 肌胃重 |
| 27 | `MWR` | 肌胃比 |
| 28 | `GAS` | 腓肠肌重 |
| 29 | `GASR` | 腓肠肌比 |
| 30 | `TM2` | 腿重(去毛皮) |
| 31 | `TMR2` | 腿重比2 |
| 32 | `RV` | 右心室重 |
| 33 | `RVR` | 右心室比 |
| 34 | `TV` | 心室总重 |
| 35 | `TVR` | 心室比 |
| 36 | `AHI` | 腹水心脏指数 |
| 37 | `TM1` | 腿重 |
| 38 | `TMR1` | 腿重比1 |
| 39 | `CW` | 爪重 |
| 40 | `CWR` | 爪比 |
| 41 | `DW` | 翅根重 |
| 42 | `DWR` | 翅根比 |
| 43 | `SW` | 翅中重 |
| 44 | `SWR` | 翅中比 |
| 45 | `TW` | 翅尖重 |
| 46 | `TWR` | 翅尖比 |
| 47 | `WW` | 翅重 |
| 48 | `WWR` | 翅重比 |
| 49 | `BFCf` | 胸肌脂肪(鲜) |
| 50 | `LFCf` | 腿肌脂肪(鲜) |
| 51 | `MCCM` | 胸肌含水量 |
| 52 | `WCLM` | 腿肌含水量 |
| 53 | `TNM` | 腓肠肌纤维总数 |
| 54 | `TM_CSA` | 腓肠肌纤维面积 |
| 55 | `PM_CSA` | 胸大肌纤维面积 |
| 56 | `EBW3` | EBW3 |
| 57 | `EBW6` | EBW6 |
| 58 | `GBW3` | GBW3 |
| 59 | `GBW6` | GBW6 |
| 60 | `GRFI` | GRFI |
| 61 | `ERFI` | ERFI |

## 各性状命中汇总

| trait | n | joint | ACJ | SACJ | alleleBon | snpBon | SACJ>SNP | ACJ>SNP |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| BW3 | 463 | 1 | 1 | 1 | 1 | 2 | 0 | 0 |
| BW6 | 463 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| BW3_6 | 463 | 1 | 1 | 1 | 0 | 0 | 1 | 1 |
| BW | 463 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| BW3_zq | 463 | 0 | 0 | 0 | 0 | 2 | 0 | 0 |
| BW6_zq | 463 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| FI | 463 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| AFI | 463 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| RFI1 | 463 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| RFI2 | 463 | 2 | 2 | 2 | 3 | 8 | 0 | 0 |
| RFI3 | 463 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| RFI4 | 463 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| FRI | 463 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| PMA | 463 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| PMAR | 463 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| PMI | 463 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| PMIR | 463 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| PM | 463 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| PMR | 463 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| LW | 463 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| LWR | 463 | 1 | 1 | 1 | 2 | 0 | 1 | 1 |
| HW | 463 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| HWR | 463 | 2 | 2 | 4 | 9 | 4 | 0 | 0 |
| AFW | 463 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| AFWR | 463 | 1 | 1 | 1 | 4 | 0 | 1 | 1 |
| MW | 463 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| MWR | 463 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| GAS | 463 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| GASR | 463 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| TM2 | 463 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| TMR2 | 463 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| RV | 463 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| RVR | 463 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| TV | 463 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| TVR | 463 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| AHI | 463 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| TM1 | 463 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| TMR1 | 463 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| CW | 463 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| CWR | 463 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| DW | 463 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| DWR | 463 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| SW | 463 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| SWR | 463 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| TW | 463 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| TWR | 463 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| WW | 463 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| WWR | 463 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| BFCf | 463 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| LFCf | 463 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| MCCM | 463 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| WCLM | 463 | 0 | 0 | 0 | 1 | 0 | 0 | 0 |
| TNM | 463 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| TM_CSA | 463 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| PM_CSA | 463 | 0 | 0 | 0 | 1 | 0 | 0 | 0 |
| EBW3 | 463 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| EBW6 | 463 | 2 | 2 | 3 | 1 | 2 | 1 | 0 |
| GBW3 | 463 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| GBW6 | 463 | 11 | 14 | 11 | 5 | 4 | 1 | 1 |
| GRFI | 463 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| ERFI | 463 | 1 | 1 | 1 | 1 | 0 | 1 | 1 |

SACJ hits>SNP: **6**/61; ACJ hits>SNP: **5**/61
