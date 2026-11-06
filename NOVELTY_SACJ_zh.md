# SACJ 相对 ACJ 的新意（v0.4，历史说明）

> **0.7.1：** 当前包不再计算 SACJ，CLI 也没有 `run-multi`。默认排序是 `omni_p`。本文件只保留 v0.4 的补充记录，不是 0.7.1 的主结果。
>
> **App Note 口径（2026-11-05）**：SACJ 保留在当时的软件与本说明中，但 **Application Note 主表不报 SACJ**；主宣称列为 joint / ACJ / allele（正式 GHap）/ SNP。详见 `appnote_hardening/HARDENING_STATUS_zh.md`。


## 一句话
**SACJ（Spectral Adaptive Collapsed Joint）**：在 HapAllele 的 **软 Hamming 亲和图**（距离 ≤ h_cap，默认 2，权重 \(e^{-d/τ}\)）上做规范化拉普拉斯谱聚类，再对簇剂量做联合 χ²；全程 **不看表型**，故与 ACJ 一样可 χ² 校准。

## 和 ACJ 差在哪
| | ACJ | SACJ |
|--|--|--|
| 图 | 硬边：Hamming≤1 | 软边：Hamming≤h_cap（默认 2）+ 指数权重 |
| 合并 | 连通分量（并查集） | 谱嵌入 + eigengap 选 k + k-means |
| 距离 2 的等位基因 | **不会**直接合并 | **可以**共簇（经中间亲和/谱结构） |
| 退化时 | — | 单簇时退回软谱基 \(Z V\)（仍设计驱动） |

## 公式（设计步）
1. \(H_{ij}=\mathrm{Hamming}(s_i,s_j)\)  
2. \(A_{ij}=e^{-H_{ij}/τ}\cdot 1\{0<H_{ij}\le h_{\mathrm{cap}}\}\)  
3. \(L=I-D^{-1/2}AD^{-1/2}\)  
4. eigengap → \(k\)；特征向量嵌入 → 簇；\(z_G=\sum_{j\in G}z_j\)  
5. \(T=\|P_{Z^*}\tilde y\|^2\sim\chi^2_{\mathrm{rank}(Z^*)}\)

## 诚实边界
- 不是表型自适应 SKAT 权重；不是深度学习。  
- 相对 ACJ 的增益预期在 **Hamming=2 共信号**；尖锐单等位仍可能稀释（与 ACJ 同类风险）。  
- 与 GHap-like 边际检验互补，应并列报告。

## 实数据快检（chr4 claw/BW）
n=2351, snps=3165, blocks=1581, Bonf thr≈3.16×10⁻⁵，OLS null。

| method | hits | min_p |
|--|--:|--:|
| joint | 3 | 3.76e-7 |
| ACJ | 2 | 1.31e-8 |
| **SACJ** | **3** | 5.23e-8 |
| GHap-like | 4 | 2.01e-8 |
| SNP-bon | 0 | 2.85e-3 |

SACJ∩ACJ=2；**SACJ\\ACJ=1**（`4:w1343`，sacj_p=1.79e-5，k=2）；ACJ\\SACJ=0。
