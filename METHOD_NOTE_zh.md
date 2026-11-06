# 状态

不要按这一稿投稿。这一稿把 ACJ 当成方法创新。后来的模拟里，单倍型都比较常见（频率至少 0.10）时，ACJ 和普通联合检验是同一个统计量。和当前代码一致的数字在 `README.md` 与 `sim_results/paper_sim.json`（v0.5.0 起的带主元 QR；0.7.1 沿用）。发现用 Cauchy 综合 p 值（`omni_p`）。`beyond_p` 是扣掉块内加性 SNP 剂量之后的单倍型联合得分：纯 SNP 时维持在显著性水平，互补单倍型 `0110` 与 `1001` 同时有效时有功效。`tree_p` 与 SACJ 不是 0.7.1 的主结果。`--lowo` 是已有亲缘矩阵上的窗口降更新，不是新的 LOCO / RHM。下面的鸡数据扫描（claw/wing）是 OLS，不是亲缘 LMM；0.7.1 的 GCTA-GRM 扫描见 README「Real broiler data」。

# hapGWAS — Application Note 风格方法说明（v0.3.1 稿，勿作 0.7.1 主文）

## Background

单倍型关联可捕获 SNP 边际难以表达的等位基因组合效应。现有工具（如 **GHap**）侧重单倍型等位基因的**边际**检验；块内多等位基因的**联合**与**自适应塌缩**仍缺少统一、可校准的 score 框架。本文/本软件在同一 null（LMM 或 OLS）投影下，给出全联合、**ACJ**、PCA/Ridge 与边际/SNP 基线的并列实现。

## Methods

### 模型与投影
\[
y=X\beta+u+e,\quad u\sim\mathcal N(0,\sigma_g^2 K),\ e\sim\mathcal N(0,\sigma_e^2 I).
\]
谱白化残差化得 \(\tilde y\)；HapAllele 剂量经同一变换得 \(Z\)（丢参考等位基因）。无 GRM 时退化为 OLS 投影。

### 块定义
默认固定 SNP 窗口（`--window/--step`）；可选 Gabriel 风格或外部 PLINK `--blocks`。

### Score 族
1. **Full joint：** \(T=\|Q_Z^\top\tilde y\|^2\sim\chi^2_{\mathrm{rank}(Z)}\)。  
2. **ACJ（本稿旧口径，0.7.1 不再当主创新）：** HapAllele 字符串 Hamming\(\le h\)（默认 1）并查集合并 → 组剂量求和 → 再联合；**不看表型**，故 χ² 可校准。  
3. **PCA / Ridge：** 方差截断 / 岭惩罚 + Satterthwaite。  
4. **Allele marginal：** 各 HapAllele 1-df score（GHap-like）。  
5. **SNP：** 块内/全染色体 1-df；比较时对 min-p 做 Bonferroni。

### 何时不要单独用 ACJ
尖锐单一 HapAllele、或需定位具体等位基因串时：优先 allele-Bonf / joint，ACJ 仅并列。Hamming 距离均 ≫1 时 ACJ 近似 joint。

## Results（摘要数字）

### 模拟（n=600, 500 reps, LMM h²=0.3）
- Null Type I：ACJ 0.044（λ_GC≈1.09）。  
- 单因果 HapAllele β=0.40：joint 0.942 > allele-Bonf 0.932 > **ACJ 0.700**（稀释）。  
- 双近邻等位 β=0.24：**ACJ 0.980** > joint 0.964 > allele-Bonf 0.722。

### 真实数据（鸡 D8 phased chip，连续 claw/BW）
- n=2351；**31** 条染色体 phased chip；**16361** blocks / **97766** alleles / **32802** SNPs；HPC walltime **381 s**。  
- 块级 Bonferroni thr = 0.05/16361 = **3.06×10⁻⁶**。  
- Hits：ACJ **8** · joint **8** · GHap-like（块内 allele-Bonf）**34** · SNP-Bonf **0** · PCA 5 · ridge 2。  
- ACJ vs GHap-like：shared **6**，ACJ-only **2**（`2:w611`,`3:w228`），GHap-only **28**。  
- 顶块仍为 `4:w1290`：ACJ **1.31×10⁻⁸**，joint 3.76×10⁻⁷，GHap-like 2.01×10⁻⁸，SNP-Bonf **0.557**（单 SNP 基因组无显著）。  
- 表/图：`case_study/gw/`。


### 真实数据 #2（翅重 wing，绝对重量）
- gemma_final 第 4 列；n=2351；mean≈238.5 g；同 31 chr / 16361 blocks；wall **383 s**。  
- Hits：ACJ **5** · joint **5** · GHap-like **6** · SNP-Bonf **0**。  
- ACJ∩GHap **2**，ACJ-only **3**，GHap-only **4**；顶块 `6:w398` ACJ **6.19×10⁻⁹**。  
- 小结：`case_study/gw_wing/`。

### v0.4 补充（不是 0.7.1 主结果）
SACJ 与 broiler 多性状记录仍在 `NOVELTY_SACJ_zh.md`、`sim_results/SACJ_RESULTS.md`、`case_study/sacj_chr4/`、`case_study/broiler_multi/`。0.7.1 扫描不输出 SACJ。

## Availability
- 许可：MIT（`LICENSE`）。  
- 引用：`CITATION.cff`（作者 Yancong Yang）。  
- 安装：`README.md`（`pip install -e .`）。  
- 差距清单：`BIB_GAP_zh.md`。  
- 公共仓库：https://github.com/YancyYang02/hapGWAS 。包版本 **0.7.1**。版本 DOI 仍是 v0.4.1 存档：[10.5281/zenodo.23984394](https://doi.org/10.5281/zenodo.23984394)。概念 DOI（全版本）：[10.5281/zenodo.23984393](https://doi.org/10.5281/zenodo.23984393)。重新存档前不要把该 DOI 写成 0.7.1。
