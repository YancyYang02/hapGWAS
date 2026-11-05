# 正式 GHap 头对头（chr4≈1500 SNP；窗4/2；同 GCTA GRM）

GHap **4.0.0.14**（ytutsunomiya）；`ghap.assoc` + 同一套 `gcta_grm`（对角正则后 PD）。

| trait | null | n | n_blocks_hapgwas | joint_hits | acj_hits | allele_bonf_hits_hapgwas | snp_gw_hits | ghap_version | ghap_n_hap_alleles | ghap_allele_hits | ghap_block_minp_hits | ghap_min_p | hap_joint_min_p | AFWR_joint_loci_touch | notes |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| BW6 | LMM+GRM | 463 | 749 | 1 | 1 | 1 | 0 | 4.0.0.14 | 8449 | 11 | 20 | 1.0040496155035e-11 | 2.9295660500064728e-06 | 1 | GHap LMM used PD-regularized GCTA GRM (eig clip 1e-4 + diag 1e-3) |
| AFWR | LMM+GRM | 463 | 749 | 42 | 25 | 73 | 1 | 4.0.0.14 | 8449 | 98 | 117 | 1.11223335667036e-103 | 4.0654901866518853e-16 | 28 | GHap LMM used PD-regularized GCTA GRM (eig clip 1e-4 + diag 1e-3) |

## 口径
- 主表：joint / ACJ / **正式 GHap** / SNP；**SACJ 降级**。
- hapGWAS 流水线 `allele_bonf` ≠ GHap；本表为正式 `ghap.assoc`。
- GHap `ghap_allele_hits`：等位基因 Bonferroni（0.05/#hap-alleles）。
- GHap `ghap_block_minp_hits`：块内 min-p 后再按 #blocks Bonferroni。
- GRM：与 hapGWAS LMM 同一 GCTA 文件；GHap/cholmod 要求 PD，故 eig clip + diag ridge（见 notes）。
- AFWR joint 窗口 hit：locus_collapse touch 42→28。

## 解读要点
- **BW6**：hapGWAS joint/ACJ/allele/snp_gw = 1/1/1/0；正式 GHap allele/block-minp = 11/20 (min_p=1.0040496155035e-11).
- **AFWR**：hapGWAS joint/ACJ/allele/snp_gw = 42/25/73/1；正式 GHap allele/block-minp = 98/117 (min_p=1.11223335667036e-103).
