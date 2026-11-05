# Broiler 纯生产面板：hapGWAS vs GHap-like 对照
## 设置
- 队列 n=463；表型 **61** 纯生产/屠体（BW3…ERFI），**无代谢物**。
- 基因型：phased VCF **前 5000 位点**；`--block-mode window --window 4 --step 2` → **2499 blocks**。
- Bonf thr = 0.05/2499 ≈ **2.00e-5**；OLS null（`--ols`）。
- **GHap 列**：本软件 `allele_bonf` / `ghap_like_bon_p`（块内 HapAllele 边际 min-p × 等位基因数），与 GHap 典型“单倍型剂量/等位关联”同口径、同块划分。
- 说明：HPC 上 CRAN **GHap R 包未预装**；已启动用户库安装（`ghap_install.log`）。在安装可用前，以同流水线 allele-Bonf 作为 head-to-head 基线（CLI 别名 `--methods ghap`）。

## 总体命中数（61 性状求和）
| joint | ACJ | SACJ | GHap-like | SNP |
|---:|---:|---:|---:|---:|
| 22 | 25 | 25 | 28 | 22 |

## 性状级优劣计数
- SACJ hits > GHap：**3**/61
- ACJ hits > GHap：**3**/61
- GHap hits > SACJ：**6**/61
- SACJ hits = GHap：**52**/61
- SACJ hits > SNP：**6**/61
- ACJ hits > SNP：**5**/61
- 任一方法有 hit 的性状：**12**/61
- 有 hit 性状中按命中数胜出：{'GHap': 5, 'SNP': 3, 'joint': 2, 'SACJ': 1, 'ACJ': 1}

## 诚实解读
- 总 hits：GHap-like(**28**) ≈ SACJ(**25**) > SNP(**22**)；GHap 在**尖锐单等位**块上更易报多 hit（如 HWR 9 vs SACJ 4）。
- SACJ/ACJ 的优势不在全面碾压，而在少数多等位/近邻同向信号（如 BW3_6、LWR、AFWR、EBW6、GBW6、ERFI 相对 SNP）。
- 相对 GHap：SACJ>GHap 的性状少；多数并列或 GHap 边际更灵敏。卖点是**同一 CLI 内联合/塌缩+基线并列**，而非处处比 GHap 强。

## 文件
- `head_to_head_by_trait.tsv` — 逐性状 hits/min_p/胜出
- `fig_a_hit_counts.png` / `fig_b_*.png` / `fig_c_architecture.png`
- `SOFTWARE_SUMMARY_zh.md`
