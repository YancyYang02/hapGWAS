# Chicken breast-shape real-data case study — data gap

## Searched (machine `53c7843e-7642-4a2b-aaa3-ad5ecb15fa86`)

| Location | Found | Usable for hapGWAS? |
|----------|-------|---------------------|
| `Desktop/毕设/` | phenotypes, scripts, docs, pedigree; **no** VCF/BED/PED genotypes | No genotypes |
| `Desktop/单倍型分析.docx` | Pipeline pointing to HPC Beagle-phased VCF + PLINK hap-assoc | Scripts only |
| `D:\phased_hap_gwas.sh` | Same HPC paths | Script only |
| `D:\12性状\D8_qc.{bed,bim,fam}` | **Unphased** chip: n=2384, 40126 SNPs | LD patterns yes; **not phased** |
| `D:\Y4_*.bed`, `D:\D8\*.bed` | Other unphased beds | Not phased |
| HPC paths in docs | `/fs2/home/clwen/YancongYang/DoubleMuscling/variant_calling_results/step6.beagle_phased.vcf.gz` and `step7_phased_final` | **Not on this PC** |

## Missing for a real breast-shape hapGWAS case study

1. Phased genotypes: `step6.beagle_phased.vcf.gz` **or** `step7_phased_final.{bed,bim,fam}` (phase-aware PED/Oxford also OK).
2. Breast-shape / 胸型 phenotype file keyed to the same sample IDs.
3. Optional: kinship / GRM matched to those samples.

## Fallback used

Empirical LD from `D8_qc` → `chip_ld/d8_ld_windows.json` (20 windows, mean adjacent |r| ≈ 0.19–0.87) → `scripts/run_chipld_sim.py` Gaussian-copula haplotype sims.
