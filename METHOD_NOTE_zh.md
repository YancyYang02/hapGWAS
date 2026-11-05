# hapGWAS — Application Note 风格方法说明（v0.3.1）

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
2. **ACJ（novelty）：** HapAllele 字符串 Hamming\(\le h\)（默认 1）并查集合并 → 组剂量求和 → 再联合；**不看表型**，故 χ² 可校准。  
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
- 表/图：`case_study/gw/`（Desktop）与 HPC `/fs2/home/clwen/YancongYang/hapgwas_case/gw_results/`。


### 真实数据 #2（翅重 wing，绝对重量）
- gemma_final 第 4 列；n=2351；mean≈238.5 g；同 31 chr / 16361 blocks；wall **383 s**。  
- Hits：ACJ **5** · joint **5** · GHap-like **6** · SNP-Bonf **0**。  
- ACJ∩GHap **2**，ACJ-only **3**，GHap-only **4**；顶块 `6:w398` ACJ **6.19×10⁻⁹**。  
- 小结：`case_study/gw_wing/`；HPC `.../hapgwas_case/gw_results_wing/`。



## SACJ（v0.4 新增，相对 ACJ）

**Spectral Adaptive Collapsed Joint**：在 HapAllele 软 Hamming 亲和图（\(A_{ij}=e^{-d/τ}\) 且 \(d\le h_{\mathrm{cap}}\)，默认 \(h_{\mathrm{cap}}=2\)）上构造规范化拉普拉斯，按 eigengap 选 \(k\)，谱嵌入 + k-means 得到簇剂量，再联合 χ²。设计步不看表型 → 与 ACJ 同属可校准 χ² 族。

- 相对 ACJ：ACJ 仅合并 Hamming≤1 连通分量；SACJ 允许距离 2 的等位基因共簇。  
- 模拟（n=500, 300 reps, LMM h²=0.3）：Type I SACJ≈0.03–0.05；**Hamming=2 双因果** power SACJ **0.910** > ACJ 0.870（joint 0.987）；Hamming=1 近邻仍 ACJ 更强（0.950 vs SACJ 0.870）。  
- 实数据 chr4 claw/BW（n=2351, 1581 blocks, thr≈3.16e-5）：hits SACJ **3** / ACJ 2 / joint 3 / GHap-like 4 / SNP 0；SACJ 独有 `4:w1343`（sacj_p=1.79e-5 < thr；acj/joint/GHap 均未过阈）。  
- 详见 `NOVELTY_SACJ_zh.md`、`sim_results/SACJ_RESULTS.md`、`case_study/sacj_chr4/`。



## Broiler 多性状整合测试（v0.4）
- 队列：`/fs2/home/broiler/zhouyan/GWAS`，n=463；表型 49（生产+代谢物 top-h²）。
- CLI：`hapgwas run` / `hapgwas run-multi --methods all` → 单表含 joint/ACJ/SACJ/allele_bonf/snp_bonf/pca/ridge。
- 结果：`case_study/broiler_multi/BROILER_MULTI_TRAIT_zh.md`；SACJ>SNP 11/49 性状。

## Availability
- 许可：MIT（`LICENSE`）。  
- 引用：`CITATION.cff`。  
- 安装：`README.md`（HPC `module load python/3.10_mambaforge_*` + `pip install -e .`）。  
- 差距清单：`BIB_GAP_zh.md`。  
- 公共仓库：https://github.com/YancyYang02/hapGWAS （`v0.4.1`）。Zenodo DOI 尚未铸造；需所有者启用 GitHub↔Zenodo 并归档 Release，勿编造 DOI。
