# Simulation results (v0.3.1)

## A. Scaled suite — n=600, n_rep=500, relatedness LMM (h²=0.3)

Script: `scripts/run_scaled_tradeoff.py`  
Figure: `fig_acj_power_tradeoff.png`  
JSON: `scaled_summary.json`

### Type I @ α=0.05 (main null)

| method | Type I | λ_GC |
|--------|--------|------|
| snp_bon_p | 0.050 | 0.11 |
| allele_bon_p | 0.052 | 0.34 |
| joint_p | 0.040 | 0.89 |
| **acj_p** | **0.044** | **1.09** |
| pca_p | 0.056 | 0.94 |
| ridge_p | 0.040 | 0.91 |

### Power — single causal HapAllele (β=0.40)

| method | power |
|--------|-------|
| joint_p | **0.942** |
| allele_bon_p | 0.932 |
| ridge_p | 0.924 |
| pca_p | 0.884 |
| snp_bon_p | 0.766 |
| **acj_p** | **0.700** |

### Power — two Hamming-neighbor alleles (β=0.24)

| method | power |
|--------|-------|
| **acj_p** | **0.980** |
| joint_p | 0.964 |
| ridge_p | 0.928 |
| snp_bon_p | 0.852 |
| pca_p | 0.770 |
| allele_bon_p | 0.722 |

Neighbor-design Type I: ACJ 0.052, joint 0.050.

### Takeaway

- ACJ wins on **neighbor multi-allele**; loses on **sharp single allele** vs joint/allele-Bonf.
- Do **not** report ACJ alone — always with joint + allele-Bonf.

## B. Earlier suite (n=400, 200 reps) — retained

See prior tables in git history / `power_*_summary.tsv`. Qualitative same tradeoff (ACJ 0.965 vs joint 0.940 on neighbors; ACJ 0.625 on single).

## C. Chip-LD fallback (D8_qc) — n=800, n_rep=400

Phased breast-shape VCF not local — see `DATA_GAP_chicken_phased.md`.  
LD: `../chip_ld/d8_ld_windows.json` (20 windows; adjacent |r| ≈ 0.19–0.87).  
Figure: `fig_chipld_acj_tradeoff.png` · JSON: `chipld_summary.json`

### Type I
ACJ 0.0475 · joint 0.0475 · snp-Bonf 0.0525 (calibrated).

### Power single HapAllele β=0.35
allele-Bonf **0.995** · ridge 0.983 · PCA 0.960 · joint 0.930 · snp-Bonf 0.838 · **ACJ 0.473** (strong dilution under chip LD).

### Power neighbor alleles β=0.22
snp-Bonf 0.933 · **ACJ 0.903** · joint 0.900 · ridge 0.875 · PCA 0.665 · allele-Bonf 0.493.

Same qualitative rule: ACJ for neighbor multi-allele; not alone for sharp single allele.
