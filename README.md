# hapGWAS (v0.4.1)

Block **haplotype GWAS** under one null-LMM (or OLS) projection, with a sparse multi-allelic score family:

| Score | Role |
|-------|------|
| **Joint full** | multi-df score on HapAllele dosages |
| **ACJ** *(novel)* | Hamming-neighborhood allele collapse → joint χ² |
| Allele-Bonf / SNP-Bonf | GHap-like + SNP baselines in same table |
| **PCA / Ridge** | high-allelic truncation / shrinkage |
| Allele marginals | GHap-like 1-df decomposition |
| SNP scores | optional conditioning / baseline |

## Novelty (honest)

**Claim:** Hamming-driven **Adaptive Collapsed Joint (ACJ)** inside a projected multi-score block pipeline.

**Do not claim:** inventing dosage→LMM, Gabriel blocks, spectral LMM, or SKAT.

**When ACJ helps:** same-direction Hamming-neighbor HapAlleles (sim power ACJ 0.980 > joint 0.964 > allele-Bonf 0.722 at n=600, 500 reps).  

**SACJ (exploratory / supplement only):** soft Hamming affinity + spectral clusters remain available via `--methods sacj` and in default `all` output, but App Note **main tables report joint / ACJ / allele (GHap) / SNP**—not SACJ as a primary claim.

**When not to use ACJ alone:** sharp single HapAllele — prefer allele-Bonf / joint. See `METHOD_NOTE_zh.md`, `sim_results/RESULTS.md`, `case_study/CASE_STUDY_zh.md`.

## Install (Linux / HPC)

```bash
# example: module Python with numpy/scipy/pandas
module load python/3.10_mambaforge_4.14.0-0   # site-specific
python -m pip install --user -i https://pypi.tuna.tsinghua.edu.cn/simple pandas  # if needed
git clone https://github.com/YancyYang02/hapGWAS hapgwas && cd hapgwas
python -m pip install --user -e .
# or without install:
export PYTHONPATH=/path/to/hapgwas:$PYTHONPATH
```

Windows / laptop: same `pip install -e .` in a Python ≥3.10 env.

## One-liner example

```bash
# phenotype TSV columns: IID  PHENO
hapgwas run \
  --geno phased_chip_chr4.vcf.gz \
  --pheno pheno_continuous.tsv \
  --block-mode window --window 4 --step 2 \
  --ols \
  --out hapgwas_out_chr4
```

Demo (no data):

```bash
hapgwas demo --out out_demo --ols
```

## Real-data case (chicken D8, HPC)

Continuous **claw/BW residual**, n=2351, phased chip VCFs.  
Results live on HPC under `.../hapgwas_case/`; see `case_study/CASE_STUDY_zh.md`.


## Methods in one run

```bash
# single trait — primary App Note columns: joint/acj + allele_bonf + snp; sacj/pca/ridge also in table (SACJ exploratory)
hapgwas run --geno phased.vcf.gz --pheno pheno.tsv --ols \
  --methods all --block-mode window --window 4 --step 2 --out out1

# multi-trait wide table (IID + many columns)
hapgwas run-multi --geno phased.vcf.gz --pheno traits_wide.tsv --ols \
  --methods all --window 4 --step 2 --out out_multi
```

Flags: `--methods all|joint,acj,sacj,allele_bonf,snp_bonf,pca,ridge` (default all; all scores still share one null).

## Availability / cite

Source code: https://github.com/YancyYang02/hapGWAS

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23984394.svg)](https://doi.org/10.5281/zenodo.23984394)

MIT — `LICENSE`. Cite release **v0.4.1** with version DOI [10.5281/zenodo.23984394](https://doi.org/10.5281/zenodo.23984394) ([Zenodo record](https://zenodo.org/records/23984394)). Concept DOI for all versions: [10.5281/zenodo.23984393](https://doi.org/10.5281/zenodo.23984393). Metadata: `CITATION.cff`. Gap checklist vs BIB App Note: `BIB_GAP_zh.md`.
