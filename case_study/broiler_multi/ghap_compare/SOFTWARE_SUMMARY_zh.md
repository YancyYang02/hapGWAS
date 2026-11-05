# hapGWAS 软件汇总（中文，BIB / PPT 用）

## 一句话
在同一 null（LMM/OLS）投影与同一块划分下，一次跑完 **joint / ACJ / SACJ / GHap-like allele-Bonf / SNP-Bonf / PCA / Ridge**，输出统一表。

## 方法（CLI 名）
| 方法 | 含义 | 角色 |
|---|---|---|
| `joint` | 块内多 HapAllele 全联合 χ² | 多等位联合基线 |
| `acj` | Hamming≤1 近邻塌缩后再联合 | **工程化新颖组合**（设计不看表型） |
| `sacj` | 软 Hamming 亲和 + 谱聚类再联合 | ACJ 扩展（允许 d=2 共簇） |
| `allele_bonf` / `ghap` | 块内等位边际 min-p × Bonferroni | **GHap-like** 基线 |
| `snp_bonf` | 块内 SNP min-p × Bonf | SNP-GWAS 对照 |
| `pca` / `ridge` | 高维截断 / 岭收缩 | 高等位块备选 |

## CLI 用法
```bash
# 单性状
hapgwas run --geno phased.vcf.gz --pheno pheno.tsv --ols \
  --methods all --block-mode window --window 4 --step 2 --out out1

# 多性状宽表
hapgwas run-multi --geno phased.vcf.gz --pheno traits_wide.tsv --ols \
  --methods all --window 4 --step 2 --out out_multi

# 仅 GHap-like + SACJ
hapgwas run ... --methods sacj,ghap,snp_bonf
```

### 输入
- 基因型：phased VCF / PED（相位假定需注明）
- 表型：`IID` + 性状列；可选协变量
- 块：固定窗口 / Gabriel / 外部 PLINK blocks

### 输出
- `assoc_blocks.tsv`：每块各方法 p / χ² / df
- `assoc_alleles.tsv`：HapAllele 边际
- `multi_trait_hits.tsv`：多性状命中汇总
- JSON meta（样本数、窗口、方法列表）

## 与现有软件对比
| | GHap (R) | PLINK hap-assoc/linear | **hapGWAS** |
|---|---|---|---|
| 单倍型边际 | ✓ | ✓（窗口） | ✓（allele-Bonf） |
| 块内联合 / 自适应塌缩 | ✗ / 弱 | ✗ | ✓ ACJ/SACJ |
| 多方法同表 | 需自拼 | 有限 | ✓ 一次 CLI |
| LMM 投影 | 有限 | 无（线性为主） | ✓ 可选 GRM/OLS |
| 社区成熟度 | 高 | 极高 | 新工具 |

## Broiler 纯生产案例数字（本轮）
- n=463；**61** 纯生产/屠体性状；phased **前 5000 位点**；2499 blocks；thr≈2e-5
- 总 hits：GHap-like **28** · SACJ **25** · ACJ **25** · joint **22** · SNP **22**
- SACJ>GHap：**3**/61（BW3_6, EBW6, GBW6）；GHap>SACJ：**6**/61；并列 **52**/61
- SACJ>SNP：**6**/61；任一方法有 hit：**12**/61

## 诚实新颖性声明
- **算创新（BIB 应用笔记级）**：把 Laplacian/近邻塌缩/联合检验拧成可调、可校准的 ACJ/SACJ，并做成多方法对照 CLI。
- **不算“前所未有新数学”**：组件别处已有；实证是局部补强，不是全面碾压 GHap/SNP。
- 发 BIB：靠方法组合 + 软件 + 对照；冲方法学顶刊不够。

## 路径
- 对照与图：`case_study/broiler_multi/ghap_compare/`
- HPC：`/fs2/home/broiler/YancongYang/hapgwas_broiler_multi/run_chr4/`
