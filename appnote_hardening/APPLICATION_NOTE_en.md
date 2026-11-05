# hapGWAS: Adaptive Collapsed Joint haplotype association under a shared LMM projection

**Bioinformatics Application Note (near-submission draft)** · software v0.3.1 · MIT License

## Abstract

Haplotype tests can capture multi-allelic architecture that single-SNP genome-wide association studies (GWAS) miss, yet widely used haplotype software often emphasizes allele-wise marginal association. We introduce **hapGWAS**, an open-source Python package that evaluates multiple block-haplotype score tests under **one** null linear (mixed) model projection: a full multi-allelic joint score, an **Adaptive Collapsed Joint (ACJ)** score obtained by phenotype-blind merging of Hamming-neighbor HapAlleles, PCA/ridge truncations, and allele- and SNP-level baselines. Simulations show that ACJ is approximately calibrated under the null and gains power when two neighboring HapAlleles share causal signal, while diluting power for a single sharp HapAllele relative to joint or allele-Bonferroni tests. In a chicken discovery cohort (n=2351) with phased chip genotypes, genome-wide scans of claw weight residualized on body weight and of absolute wing weight show complementary discoveries versus a GHap-like allele-Bonferroni baseline; block-level SNP Bonferroni recovered no genome-wide hits for either trait. hapGWAS is intended for exploratory haplotype scanning alongside—not instead of—marginal haplotype and SNP analyses.

## Introduction

Association studies in livestock and other structured populations frequently use linear mixed models (LMMs) to control relatedness Zhou & Stephens, Nat Genet 44:821–824 (2012); Lippert et al., Nat Methods 8:833–835 (2011). When causal variants are poorly tagged by individual SNPs, haplotype or local multi-allelic encodings can increase power. Packages such as **GHap** Utsunomiya et al., Bioinformatics 32:2861–2862 (2016), doi:10.1093/bioinformatics/btw356 provide flexible haplotype dosage construction and largely **marginal** tests of individual haplotype alleles. Full joint tests of many HapAlleles in a block can be more powerful for diffuse allelic architectures, but rare alleles inflate degrees of freedom and can reduce power. Phenotype-dependent adaptive procedures risk invalid null distributions; phenotype-blind design choices (for example merging alleles that differ by at most one SNP) can reduce dimensionality while remaining χ²-calibrated under the projected null.

Our **novelty claim is deliberately narrow**: hapGWAS implements **ACJ**—union-find collapse of HapAllele strings within Hamming distance ≤ h (default 1), followed by a joint score on merged dosages—inside a **shared residualized projection** that also returns full joint, PCA/ridge, allele-marginal (GHap-like), and SNP scores for the same blocks. We do **not** claim to invent dosage-based LMM association, LD-block discovery, spectral LMM projection, or SKAT-style kernel tests. The practical contribution is a single, reproducible scan that makes the tradeoff between joint, collapsed-joint, and marginal haplotype tests explicit.

## Methods

### Null model and projection

For phenotype \(y\) and optional covariates \(X\), we fit
\[
y = X\beta + u + e,\qquad u\sim\mathcal N(0,\sigma_g^2 K),\quad e\sim\mathcal N(0,\sigma_e^2 I),
\]
with kinship \(K\) when available. Following spectral LMM / FaST-style practice Lippert et al., Nat Methods 8:833–835 (2011), doi:10.1038/nmeth.1681, residuals are whitened to \(\tilde y\). HapAllele (or SNP) dosages are transformed by the same operator. When no GRM is supplied, the implementation falls back to OLS residualization after covariates (used for the chicken genome-wide scans reported below). Score statistics are then computed on \(\tilde y\).

### Blocks and haplotype dosages

Blocks may be fixed SNP windows (`--window`, `--step`) or user-supplied intervals (e.g. external PLINK blocks). Within each block, phased genotypes (VCF) define HapAllele strings; dosages count copies of each allele (0–2 per individual). A reference allele is dropped to avoid collinearity. Phasing of real data may use Beagle or equivalent Browning et al., Am J Hum Genet 103:338–348 (2018), doi:10.1016/j.ajhg.2018.07.015; hapGWAS itself consumes already-phased calls.

### Score family

1. **Full joint.** Let \(Z\) be the projected HapAllele dosage matrix. The joint score is \(T=\|Q_Z^\top\tilde y\|^2\) compared to \(\chi^2_{\mathrm{rank}(Z)}\).
2. **ACJ (Adaptive Collapsed Joint).** HapAlleles whose string labels differ by Hamming distance ≤ h (default 1) are merged with a union-find; dosages are summed within merges; the joint score is reapplied. Merging depends only on genotype design, **not** on \(y\), so the χ² reference remains appropriate under the projected null.
3. **PCA / ridge.** Truncation or ridge shrinkage of the dosage Gram matrix with Satterthwaite degrees of freedom for approximate p-values.
4. **Allele marginal (GHap-like).** One-degree-of-freedom scores for each HapAllele dosage; within-block Bonferroni on the minimum p-value summarizes the block for discovery tables.
5. **SNP.** Analogous one-df SNP scores; within-block SNP Bonferroni for comparison.

A CLI (`hapgwas run`) writes per-block and per-allele tables (`joint_p`, `acj_p`, `pca_p`, `ridge_p`, allele `marg_p`, SNP p-values). Genome-wide discovery in this note uses Bonferroni on the number of blocks (\(\alpha=0.05/n_{\mathrm{blocks}}\)); GHap-like and SNP block summaries apply within-block Bonferroni before the same genome threshold.

### Honest use of ACJ

ACJ is recommended when neighboring HapAlleles may share causal signal. For a single sharp HapAllele, allele-Bonferroni or full joint typically dominate; ACJ can dilute signal by merging causal with non-causal neighbors. Exploratory scans should report **joint + ACJ + allele-Bonferroni** (and SNP) together. A spectral soft-collapse score (**SACJ**) is implemented in the software for sensitivity analyses but is **not** a primary claim of this Application Note.

## Simulation results

We summarize a scaled suite (n=600 individuals, 500 replicates, relatedness LMM with heritability 0.3; details and figures in package `sim_results/`).

**Type I error** at \(\alpha=0.05\): ACJ 0.044 (\(\lambda_{\mathrm{GC}}\approx1.09\)), joint 0.040, allele-Bonferroni 0.052, SNP-Bonferroni 0.050—approximately calibrated.

**Power, single causal HapAllele** (\(\beta=0.40\)): joint 0.942 > allele-Bonferroni 0.932 > ridge 0.924 > PCA 0.884 > SNP-Bonferroni 0.766 > **ACJ 0.700**. ACJ loses here by dilution.

**Power, two Hamming-neighbor causal alleles** (\(\beta=0.24\) each): **ACJ 0.980** > joint 0.964 > ridge 0.928 > SNP-Bonferroni 0.852 > PCA 0.770 > allele-Bonferroni 0.722. Neighbor-design Type I remained near nominal (ACJ 0.052, joint 0.050).

A chip-LD fallback suite built from empirical D8 LD windows (n=800, 400 replicates) reproduced the same qualitative tradeoff: strong ACJ dilution under sharp single-allele architecture (ACJ power 0.473 vs allele-Bonferroni 0.995 at \(\beta=0.35\)) and competitive ACJ power for neighbor alleles (0.903 vs joint 0.900; allele-Bonferroni 0.493 at \(\beta=0.22\)).

## Real-data results

### Cohort and scan design

Chicken discovery line D8: continuous phenotypes aligned to phased chip VCFs (31 chromosomes with available `phased_chip_chr*.vcf.gz`). After removing missing phenotypes, **n=2351**. Scans used fixed windows of 4 SNPs with step 2, OLS null (no GRM in the reported runs), and genome-wide Bonferroni on **16 361** blocks (threshold \(3.06\times10^{-6}\)). The marker set comprised **32 802** SNPs and **97 766** HapAlleles after filtering. Wall-clock per trait was ≈6.4 min on institutional Linux HPC (Python 3.10 mambaforge). Phasing and prior SNP GWAS pipelines used standard tools Browning et al., Am J Hum Genet 103:338–348 (2018); Zhou & Stephens, Nat Genet 44:821–824 (2012); genotypes are not redistributed with the software.

### Trait 1: claw weight residualized on body weight

Hits at block Bonferroni:

| Method | Hits | Minimum p |
|--------|-----:|----------:|
| ACJ | 8 | \(1.31\times10^{-8}\) |
| Full joint | 8 | \(2.55\times10^{-10}\) |
| GHap-like (allele Bonf.) | 34 | \(2.50\times10^{-12}\) |
| SNP Bonf. | 0 | \(1.24\times10^{-5}\) |
| PCA | 5 | \(6.99\times10^{-11}\) |
| Ridge | 2 | \(6.15\times10^{-8}\) |

ACJ∩GHap-like: **6** blocks; ACJ-only: **2**; GHap-like-only: **28**. The lead ACJ block `4:w1290` had ACJ \(p=1.31\times10^{-8}\), joint \(3.76\times10^{-7}\), GHap-like \(2.01\times10^{-8}\), and **SNP Bonferroni 0.557**—a clear haplotype-only discovery relative to the SNP block summary.

### Trait 2: absolute wing weight

Wing weight (mean 238.5 g, SD 33.3) used the same individuals, windows, and multiple-testing rule.

| Method | Hits | Minimum p |
|--------|-----:|----------:|
| ACJ | 5 | \(6.19\times10^{-9}\) |
| Full joint | 5 | \(1.81\times10^{-8}\) |
| GHap-like (allele Bonf.) | 6 | \(2.64\times10^{-8}\) |
| SNP Bonf. | 0 | \(8.27\times10^{-6}\) |
| PCA | 0 | \(4.15\times10^{-6}\) |
| Ridge | 2 | \(3.40\times10^{-7}\) |

ACJ∩GHap-like: **2**; ACJ-only: **3** (`1:w2323`, `4:w105`, `6:w398`); GHap-like-only: **4**. Lead ACJ block `6:w398`: ACJ \(6.19\times10^{-9}\), joint \(1.81\times10^{-8}\), GHap-like Bonferroni \(2.40\times10^{-5}\), SNP Bonferroni 0.0025.

## Comparison to GHap-like analysis

For a fair head-to-head on **identical blocks, dosages, and null projection**, we treat hapGWAS allele-marginal scores with within-block Bonferroni as a **GHap-like** baseline. This is a proxy for the marginal haplotype-allele workflow commonly used with GHap Utsunomiya et al., Bioinformatics 32:2861–2862 (2016), doi:10.1093/bioinformatics/btw356; it is not a reimplementation of the full R package API. Across both traits, GHap-like recovered more (claw/BW) or similar (wing) Bonferroni hits than ACJ, consistent with many sharp single-allele peaks. ACJ nevertheless contributed unique blocks (2 and 3 respectively) and, with the full joint, highlighted loci where SNP Bonferroni was null (all 8 claw/BW ACJ hits; all 5 wing ACJ hits). We therefore recommend reporting ACJ **with** GHap-like marginal and SNP summaries rather than claiming universal superiority.

## Availability and Implementation

- **Project name:** hapGWAS  
- **Version:** 0.4.1 (GRM/LMM CLI; `--grm` / `--kinship`)  
- **License:** MIT  
- **Language:** Python 3 (≥3.9 recommended); dependencies NumPy, SciPy, pandas  
- **Operating system:** platform-independent; tested on Linux HPC with `module load python/3.10_mambaforge_*`  
- **Install and one-liner:**
```bash
pip install -e .
hapgwas run --geno phased.vcf.gz --pheno pheno.tsv \
  --window 4 --step 2 --ols --out hapgwas_out
```
- **Citation metadata:** `CITATION.cff`  
- **Documentation:** `README.md`, `METHOD_NOTE_zh.md`, simulation summaries under `sim_results/`  
- **Public source repository:** https://github.com/YancyYang02/hapGWAS (tag `v0.4.1`)  
- **Zenodo DOI (v0.4.1):** https://doi.org/10.5281/zenodo.23984394 ([record](https://zenodo.org/records/23984394))  
- **Zenodo concept DOI (all versions):** https://doi.org/10.5281/zenodo.23984393  
- **Hardening note (chr4 ~1.5k subset, Chinese summary):** `appnote_hardening/HARDENING_STATUS_zh.md` — LMM λ, official GHap v4.0.0.14 H2H (GHap often more hits than hapGWAS), locus collapse 42→28, conditional survival at RXFP2-proximal leads with OTHER alleles.  
- **Restrictions:** Example chicken genotypes/phenotypes remain under institutional data governance and are not bundled; scripts in `scripts/` document HPC analysis layouts used for this note.

## Acknowledgements

[Institutional acknowledgements, funding, and HPC support to be completed by the authors.] We thank collaborators who generated the chicken D8 phenotypes and phased chip resources analysed computationally for this note. Genotype data are not redistributed with the software.

## References

1. Utsunomiya et al., Bioinformatics 32:2861–2862 (2016), doi:10.1093/bioinformatics/btw356 — GHap: haplotype analysis package (citation to be completed).  
2. Zhou & Stephens, Nat Genet 44:821–824 (2012), doi:10.1038/ng.2310 — GEMMA linear mixed model association (citation to be completed).  
3. Lippert et al., Nat Methods 8:833–835 (2011), doi:10.1038/nmeth.1681 — FaST-LMM / spectral LMM projection (citation to be completed).  
4. Browning et al., Am J Hum Genet 103:338–348 (2018), doi:10.1016/j.ajhg.2018.07.015 — Beagle genotype phasing (citation to be completed).  

*(Do not treat see References tokens as bibliographic entries until replaced with verified citations.)*
