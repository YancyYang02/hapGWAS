# hapGWAS 0.7.1

One command for a phased haplotype scan under a relatedness model, in the role GEMMA plays for SNP GWAS: a kinship null, a default test that does not require you to guess the allelic architecture, and a SNP scan in the same projection so the comparison is fair.

```bash
git clone https://github.com/YancyYang02/hapGWAS
cd hapGWAS
pip install -e .
hapgwas run \
  --geno phased.vcf.gz \
  --pheno pheno.tsv \
  --kin-geno other_chromosomes.vcf.gz \
  --block-mode window --window 4 --step 2 \
  --out hapgwas_out
```

`--kin-geno` is the leave-one-chromosome-out kinship (test markers in `--geno`, GRM from everywhere else). If the genotype file itself contains two or more chromosomes, `--grm auto` builds that GRM per chromosome. `--gcta-grm prefix` reads a GCTA `--make-grm` file (`.grm.id` and `.grm.bin`) and aligns it on IID. `--ols` turns relatedness off. `--lowo` downdates SNP factors from the scanned panel or from `--kin-geno`; a GCTA binary has no SNP factors, so it cannot be combined with `--lowo`.

Leaving markers that tag the tested locus inside the kinship is the proximal-contamination problem in mixed-model GWAS. Leave-one-chromosome-out (LOCO) kinship is the usual remedy and is already how GEMMA and GCTA are run (Listgarten et al., *Nature Methods* 2012; Yang et al., *Nature Genetics* 2014). `--kin-geno`, and `--grm auto` when one file contains two or more chromosomes, follow that LOCO pattern. `--lowo` is narrower: after one decomposition of K it removes SNP factors for the window plus a flank, using the matrix-determinant lemma and the Woodbury identity, on the same VanRaden estimator. This package does not claim to have introduced LOCO or regional heritability mapping.

## What is in the table

Block ranking and `assoc_loci.tsv` use `omni_p`. `beyond_p` and `beyond_lrt_p` are a follow-up: whether any haplotype signal remains after the additive SNPs in the block. `tree_p` is a secondary kernel column and is not in the omnibus. The v0.4 spectral collapse (SACJ) is not computed in 0.7.1; the earlier note stays in `NOVELTY_SACJ_zh.md` as a historical supplement.

| Column | What it is |
|---|---|
| `omni_p` | Default. Cauchy combination (ACAT) of the joint haplotype score, the allele-marginal scores, and the SNP scores in the same window |
| `tree_p` | Brownian kernel on the UPGMA tree of haplotype Hamming distances. Higher power only when several common haplotypes in one clade share an effect; worse for a single haplotype or a single SNP. Not part of `omni_p` |
| `beyond_p` | Joint score of the haplotype dosages after they have been residualized on the additive SNP dosages in the same block. Calibrated when the truth is a SNP. Power when the haplotype pattern is outside that span. Not part of `omni_p` |
| `joint_p` | Multi-df score on HapAllele dosages. Uses the variance ratio from the null model |
| `joint_lrt_p` | Exact likelihood ratio for the same haplotypes. Profiles δ = ve/vg again with the block in the model. One eigen-decomposition of K is shared |
| `joint_lowo_p` | Same likelihood ratio after the window and a flank are removed from K. One decomposition, then a Woodbury downdate per region. Use `--lowo` |
| `beyond_lrt_p` | Exact likelihood ratio for the haplotype directions outside the block's SNP dosages. The SNP model is nested in the haplotype model, and each profile has its own δ |
| `acj_p` | Joint score after merging only haplotypes below frequency 0.10 at Hamming distance ≤ 1. Common haplotypes are not chained together |
| `allele_bonf_p` | Smallest allele p × number of alleles. This is the GHap-like baseline, not the R package GHap |
| `snp_bonf_p` | Smallest SNP p in the window × number of SNPs in the window |
| `assoc_snps.tsv` | Ordinary 1-df SNP scan. Genome-wide threshold is 0.05 / n_SNPs |
| `assoc_loci.tsv` | Overlapping windows collapsed to one row per locus |

ACAT is Liu & Xie, *JASA* 2020. The kinship estimator is VanRaden. The scan uses one eigen-decomposition of K. Score columns keep the null-model δ, which is the EMMAX approximation. `joint_lrt_p` and `beyond_lrt_p` profile the Gaussian ML likelihood again inside the block, which is the GEMMA calculation, extended to a multi-column design. A SNP dosage is the sum of the haplotype dosages that carry that allele, so the SNP columns sit inside the haplotype columns. The beyond likelihood ratio uses that nesting: its degrees of freedom are rank(haplotypes) − rank(SNPs), not the number of alleles and not the number of SNPs. REML is not used for this comparison, because adding columns changes the REML contrast.

## What the simulations actually say

`sim_results/omni_bench.json`, 150 replicates, α = 0.05.

Independent samples, null: omni 0.053, joint 0.047, allele-Bonferroni 0.047, SNP-Bonferroni 0.040.

Correct LMM (phenotype drawn from the GRM that was fitted), n = 400: omni 0.040, joint 0.020, allele-Bonferroni 0.040, SNP-Bonferroni 0.033. Profile REML on a planted heritability of 0.40 returned 0.39.

Power when the architecture is not known in advance (n = 240):

| Architecture | omni | joint | allele-Bonf | SNP-Bonf |
|---|---:|---:|---:|---:|
| one HapAllele | 0.71 | 0.65 | 0.71 | 0.33 |
| two neighboring HapAlleles | 0.67 | 0.65 | 0.60 | 0.45 |
| one SNP | 0.85 | 0.74 | 0.44 | 0.89 |
| mean of the three | **0.74** | 0.68 | 0.58 | 0.56 |

Omni stays next to the best test in each row. It does not beat that best test. Mean power is higher than any one of the single tests because those single tests each fail on the architecture they were not built for.

`sim_results/tree_bench.json` is a separate, pre-specified check of `tree_p` (800 replicates, n = 500). Null rejection was 0.054. When four common haplotypes at Hamming distance 2 shared one effect, tree power was 0.27, joint 0.17, allele-Bonferroni 0.11, SNP-Bonferroni 0.14. When only one of those haplotypes was causal, allele-Bonferroni was 0.82 and tree was 0.57. When one SNP was causal, SNP-Bonferroni was 0.79 and tree was 0.53.

`sim_results/paper_sim.json` is the confirmatory run (n = 600, seed 23, pivoted QR rank). Six haplotypes at equal frequency: `0000`, `1111`, `0110`, `1001`, `0011`, `1100`. Null uses 4000 independent replicates and 2000 LMM replicates (planted h² = 0.40). Each alternative uses 1000 replicates. Effect sizes were fixed before the run: one haplotype β = 0.35, one SNP β = 0.22, complementary pair β = 0.28.

Null rejection at α = 0.05 / 0.01 / 0.001, independent samples: omni 0.054 / 0.0098 / 0.0005, beyond 0.050 / 0.0078 / 0.0008. Under the LMM null, beyond was 0.058 / 0.010 / 0.001.

Power at α = 0.05:

| Architecture | omni | joint | allele-Bonf | SNP-Bonf | beyond |
|---|---:|---:|---:|---:|---:|
| one HapAllele (`0110`) | 0.958 | 0.955 | 0.855 | 0.699 | 0.734 |
| one SNP (first site) | 0.901 | 0.867 | 0.636 | 0.922 | 0.047 |
| `0110` and `1001` together | 0.936 | 0.964 | 0.845 | 0.039 | 0.989 |

The middle row is the calibration check for `beyond_p`: the phenotype is an additive SNP, and the residual haplotype score rejects at 0.047 / 0.009 / 0.001. The same three rates under the LMM version of that SNP are 0.050 / 0.009 / 0.001. The last row is the complementary pair. Those two haplotypes differ at every SNP, so one copy of each gives SNP dosage `(1,1,1,1)`, the same dosage as `1111` plus `0000`. The SNP model cannot separate them. SNP-Bonferroni stays at 0.039; `beyond_p` reaches 0.989 at α = 0.05 and 0.811 at α = 0.001.

Use `omni_p` to find blocks. Use `beyond_p` to ask whether a block still has signal after the additive SNPs in that block are removed.

`sim_results/lrt_bench.json` compares the score columns with the exact block likelihood ratios. n = 320, sibships of 8 with relatedness 0.5, and the fitted kinship is the one used to simulate. Null LMM, 800 replicates, α = 0.05: joint score 0.040, joint LRT 0.044, beyond score 0.050, beyond LRT 0.051.

With background vg = ve = 0.5 and one haplotype at β = 0.35, power at α = 0.05 is 0.765 (score) and 0.770 (exact). At β = 0.60 both reach 1. The median −log10 p moves from 6.40 to 6.96. When that haplotype is large relative to the background (β = 0.70, vg = ve = 0.25), the median −log10 p moves from 15.1 to 17.9. A pure SNP at β = 0.55 leaves the beyond LRT at 0.055. The complementary pair at β = 0.40 leaves SNP-Bonferroni at 0.038 and puts both beyond tests at 0.995. The exact ratio is the calculation to trust when a block explains a large share of the variance and the sample contains close relatives. For a small effect it stays next to the score.

`sim_results/lowo_bench.json` records the package LOWO check. n = 360. The random effect uses 200 background SNPs only. The fitted full GRM also contains the 4-SNP window and 36 flank SNPs built as noisy copies of that window, so the region can absorb a haplotype signal. `--lowo` removes those 40 SNPs with the Woodbury update. Null, 400 replicates, α = 0.05: region-excluded likelihood ratio 0.043, full-GRM likelihood ratio 0.0025, full-GRM score 0.0025. One haplotype at β = 0.45, 300 replicates: region-excluded power 0.963, full likelihood ratio 0.587, full score 0.567. At β = 0.80 the median −log10 p is 12.9 against 5.7. The gain is the region leaving the GRM. A flank that does not tag the haplotype has nothing to remove.

Family structure that is only partly captured by the GRM is still badly inflated (OLS omni type I 0.84; LMM 0.31). A haplotype scan is not trustworthy until `run_meta.json` shows raw-p λ_GC near 1. Do not compute λ on Bonferroni p-values.

## Real broiler data

These scans used the phased chromosome whose VCF chromosome field is `40` (465 birds) and the GCTA GRM built from 8,251,083 autosomal SNPs. After the phenotype merge, n = 463. Windows are 4 SNPs, step 2. A window is counted when its raw p-value is below 0.05 / n_windows. Windows within 250 kb are then one region: a step of 2 leaves non-overlapping copies of a single peak, so an overlap-only clump over-counts. λ_GC is the median of a 1-df chi-square transformed from the raw p-value, divided by 0.4549. The REML search reports heritability on a grid; the values below are that grid, not a continuous MLE to many decimals. EBW6, GBW6 and ERFI pinned δ = ve/vg to the bottom of the grid, so those three traits are not interpreted.

1500 SNPs, 749 windows, threshold 6.68×10⁻⁵. Cluster output: `sim_results/gcta_hap.json`.

| Trait | h² | λ joint | λ LRT | Regions | Lead window | joint LRT p | in-window SNP-Bonferroni |
|---|---:|---:|---:|---:|---|---:|---:|
| BW6 | 0.24 | 1.12 | 1.18 | 1 | 625414–627113 | 1.35×10⁻⁶ | 0.66 |
| AFWR | 0.08 | 0.006 | 0.006 | 1 | 699040–699111 | 9.83×10⁻¹⁸ | 0.047 |
| BW3_6 | 0.24 | 1.03 | 1.06 | 0 | 625414–627113 is the minimum, 7.47×10⁻⁵ | — | — |

AFWR's 49 windows that pass the threshold are one interval, 320282–699691. Its λ is not a calibration result: most windows on this chromosome are conservative, and one long interval is not. BW6 is the calibration row.

5000 SNPs, 2499 windows, threshold 2.00×10⁻⁵. BW6 λ for the joint LRT is 1.54. Two regions: 625414–801545, led by 801457–801545 at 5.80×10⁻²¹ (in-window SNP-Bonferroni 5.6×10⁻⁴), and a single window 1220763–1221282 at 1.53×10⁻⁵. AFWR is one interval, 320282–1060476, led by 973499–973644 at 1.97×10⁻¹⁸. BW3_6 calls the same 801457–801545 window at 1.39×10⁻¹⁶. The 625414 window is unchanged from the 1500-SNP scan (1.35×10⁻⁶); the longer panel adds the stronger window downstream of it.

GHap 4.0.0.14 on the same 749 windows, without this GRM: BW6 has 20 block-minimum hits and a smallest p of 1.00×10⁻¹¹, including windows the LMM does not call. The shared window 625414–627113 is GHap p = 1.79×10⁻¹⁰ and joint LRT p = 1.35×10⁻⁶. AFWR has 117 GHap hits and a smallest p of 1.11×10⁻¹⁰³. The LMM keeps one interval. That is the difference relatedness makes on this chromosome. It is not a claim that one program is more powerful at every locus.

A VanRaden kinship built only from the scanned 1500 or 5000 SNPs has h² near 0 for BW6 and AFWR. `joint_lowo_p` then matches `joint_lrt_p`, because that matrix is not the genome-wide relationship matrix. Those tables are the panel control, not the association results above.

Region exclusion was checked on the chip itself. At 40:4921 the PLINK A1 dosage matches the phased ALT dosage with correlation 1 (424 birds). For one window the Woodbury downdate and the GRM rebuilt without those SNPs differ by at most 6.7×10⁻¹⁶. Over 749 windows, with kinship from 300 SNPs on the other chromosomes, the median absolute difference in log10 p between the downdate and that rebuilt kinship is 2.6×10⁻⁶.

The p-value moves when the flanking markers are a real part of the kinship. BW6 window 365062–368158, using 2000 other-chromosome SNPs (h² about 0.19) plus up to 40 chip SNPs within 20 kb: p = 0.075 with those SNPs left in, p = 0.0083 after the downdate, and p = 0.0083 if they were never included. Neighbouring windows move in the same direction by a smaller factor. This comparison uses the 428 birds present in both the phased file and the chip. It is not the genome-wide significant lead; it is the real-genotype counterpart of the simulation, where the flank is inside K.

Removing the same 40 SNPs from the 8,251,083-SNP GCTA matrix changes the AFWR lead 699040–699111 from 4.5750×10⁻²⁰ to 4.5747×10⁻²⁰. Forty markers are not a meaningful fraction of that GRM. Use `--lowo` when the kinship still contains the local tags. Use `--gcta-grm` for the dense genome-wide null.

The earlier Ensembl GRCg7b lookup, which reads VCF chromosome 40 as chromosome 4, places RXFP2 57 kb from 625414–627113 and 130 kb from 699040–699111. Neither window overlaps the gene. The strongest allele class in the earlier conditional table was the collapsed OTHER class, not one 4-SNP string. The 801 kb and 973 kb windows were not part of that lookup, and this note does not assign them a gene symbol.

## Outputs

- `assoc_blocks.tsv`, `assoc_alleles.tsv`, `assoc_snps.tsv`, `assoc_loci.tsv`, `run_meta.json`
- `run_meta.json` records heritability, raw-p λ_GC, block hits, and locus hits after merging overlaps

## Demo

```bash
hapgwas demo --out demo_out
hapgwas run --geno demo_out/phased_demo.npz --pheno demo_out/pheno.tsv --out demo_run
```

Benchmark:

```bash
python scripts/bench_omni.py --reps 150
```

## Availability / cite

Source code: https://github.com/YancyYang02/hapGWAS

Author metadata: `CITATION.cff` (Yancong Yang, China Agricultural University). License: MIT (`LICENSE`).

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23984394.svg)](https://doi.org/10.5281/zenodo.23984394)

This tree is package version **0.7.1**. The Zenodo version DOI [10.5281/zenodo.23984394](https://doi.org/10.5281/zenodo.23984394) ([record](https://zenodo.org/records/23984394)) and the concept DOI [10.5281/zenodo.23984393](https://doi.org/10.5281/zenodo.23984393) are unchanged and still record the **v0.4.1** deposit. They will be refreshed when a new Zenodo version is minted. Cite the software version you ran, and do not describe the current DOI as a 0.7.1 deposit.
