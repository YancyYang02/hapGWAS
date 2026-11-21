# hapdc

Multi-scale haplotype and diplotype scan for genome-wide association.

`hapdc` tests the joint effect of the haplotype and diplotype configuration of
a genomic region on a quantitative or binary trait, in a linear mixed model
that accounts for population structure through a genomic relationship matrix.
Each anchor region is tested at seven scales at once -- single markers,
sliding windows of 2 to 16 markers, LD blocks and the whole anchor -- and the
per-scale p-values are combined by ACAT into one p-value per anchor. The
genome-wide threshold is a Bonferroni correction over anchors; an empirical
threshold from phenotype permutation can be used instead.

## Install

```bash
git clone https://github.com/YancyYang02/hapdc.git
cd hapdc
pip install -e .
hapdc --help
```

Requires Python 3.9 or newer, numpy, scipy and threadpoolctl. No compiler,
no external binary, no R.

## Quick start

```bash
# 1. intersect a phased VCF, a GRM and a phenotype table
hapdc prep genotypes.vcf.gz auto pheno.txt out/prep

# 2. multi-scale anchor scan of trait column 0
hapdc ms out/prep out/ms0 0 --anchor-min 16 --jobs 8
```

`prep` writes the haplotype arrays, the relationship matrix and the phenotype
matrix into one directory. `ms` reads that directory and writes, for the trait
named in the phenotype header, `<trait>.anchors.tsv` and
`<trait>.summary.json`. `--jobs N` scans N chromosomes in parallel and returns
the same result as `--jobs 1` (parallel mode uses `fork`, so Linux or macOS).

## Output

`<trait>.anchors.tsv`, one row per anchor:

| column | meaning |
| --- | --- |
| `chrom`, `start_bp`, `end_bp`, `n_markers` | anchor position and size |
| `P` | anchor p-value (statistic D, combined over the seven scales) |
| `lead_scale`, `lead_unit_start`, `lead_unit_end` | scale and position of the unit with the smallest p-value |
| `P_nonadditive` | scale-weighted ACAT of the `NADD_F` p-values in the anchor: departure from additivity *of the grouped haplotype classes* (after rare-haplotype reassignment), not a test against every additive SNP model (see "What is tested") |
| `sig_bonf` | 1 if `P` < 0.05 / number of anchors |

`<trait>.summary.json` records the sample size, number of anchors, Bonferroni
threshold, number of significant anchors, the genomic-control lambda of the
anchor p-values, the phenotype transform actually applied, the covariates used
and any covariates dropped by the rank check (with the reason), every option,
and a provenance block (hapdc, Python, NumPy and SciPy versions, the command
line, and for every input file its absolute path, size, MD5 and SHA-256; inputs
are keyed by absolute path plus MD5, so files with the same name in different
directories are all recorded). `hapdc merge` also lists, per input table, the
trait, sample size, transform, covariates, hapdc version and key options read
from the table's own `summary.json`.

### Reading lambda_gc

The lambda of anchor p-values is printed and stored as a description of the
scan, not as a calibration test. Under the null it is below 1, because the
Cauchy combination is conservative away from the tail, and how far below 1
depends on the reference panel, marker density, anchor count and phenotype
transform; a low lambda neither shows nor rules out a miscalibrated tail. No
generic reference range is built in. If you have a condition-matched null range
(lambda over null phenotypes scanned with the same panel and options), pass it
as `--lambda-ref LO,HI` and the scan reports whether the observed value lies in
it; this never changes the threshold. To check calibration, use the permutation
threshold under "Thresholds".

## Inputs

* **Genotypes** -- a phased, biallelic VCF (`0|1`). A marker is dropped, and
  counted under its reason in the `prep` log, if any sample has a missing call,
  an unphased heterozygote (`0/1`), a haploid call, or an allele index other
  than 0/1 (split multi-allelic sites with `bcftools norm -m-` first), or if it
  does not vary. Unphased homozygotes (`0/0`, `1/1`) carry no phase information
  and are accepted. Genotypes are read as 0/1 haplotypes, two rows per individual.
* **Relationship matrix** -- either a GCTA `.grm` prefix, or `auto` to compute
  a VanRaden GRM from the VCF itself (`--grm-thin` markers are used).
* **Phenotypes** -- a table whose first column is the individual ID and whose
  remaining columns are traits, one column per trait, missing values coded
  `-9` by default (`--miss`). Numeric covariates, for example leading principal
  components or batch indicators, go in a separate file passed to `hapdc ms
  --covar`, one header line, IID followed by numeric columns. The intercept is
  always added; do not include a column of ones. For a categorical covariate
  with L levels supply L - 1 indicator columns.
* **Covariate rank check** -- before the variance components are fitted, the
  design [intercept, covariates] is checked on the samples actually analysed
  (after removing samples with a missing phenotype or covariate). Columns are
  examined in file order; a column that is constant, duplicates an earlier
  column, or is a linear combination of the intercept and earlier columns
  (relative tolerance 1e-8 on unit-norm columns) is dropped and logged
  (`--covar-rank drop`, default) or stops the run with the offending column
  names (`--covar-rank error`). Dropped columns are listed in `summary.json`.
  A covariate that is nearly, but not exactly, collinear is kept.

## What is tested

Markers are grouped into LD blocks, and consecutive blocks are merged into
anchors of at least `--anchor-min` markers (default 16). Blocks are never cut
and anchors never cross a chromosome boundary. Inside each anchor seven scales
are tested:

| scale | units |
| --- | --- |
| SNP | single markers |
| W2, W4, W8, W16 | sliding windows of 2, 4, 8 and 16 markers, 50% overlap |
| LD | the LD blocks inside the anchor |
| ANCHOR | the anchor itself |

Inside one unit the haplotypes present above `--minf` are retained; rarer
haplotypes are reassigned to the nearest retained haplotype by Hamming
distance (the more frequent one when several are equally near). Individuals are
grouped into diplotype cells; individuals in cells below `--minc` are moved to
the most frequent retained cell sharing one of their haplotypes (else to the
most frequent retained cell). Exact ties (equal distance and equal frequency,
or equally frequent target cells) are not broken by allele coding but shared
equally: a haplotype copy tied between t classes adds 1/t to each class dosage,
and an individual tied between t target cells has membership 1/t in each. Each
individual's class dosages therefore sum to 2 and its cell memberships to 1,
but with ties they are fractional rather than one class pair and one cell.
Without ties the design is 0/1. The results do not change when the alleles of
any marker are recoded, the signs of GRM eigenvectors are flipped, or the
samples are reordered. A unit or marker that admits no test (for example a
monomorphic marker, or one explained by the covariates) is left out of the
combination; this is decided from the genotypes and covariates alone, so a
phenotype gets the same p-value whether it is scanned alone or together with
other traits. Four statistics are computed per unit from the same mixed-model
fit (variance ratio estimated once per trait by REML under the null model):

* `A_joint` -- joint F test of the haplotype-class dosages against the null
  model (intercept and covariates);
* `CELL_F` -- joint F test of the diplotype-cell memberships against the null
  model (dominance, overdominance, recessive and compound-heterozygous effects);
* `SACAT` -- Cauchy combination of the single-marker additive t tests in the unit;
* `NADD_F` -- F test of the diplotype cells conditional on the haplotype-class
  dosages: departure from additivity of the grouped haplotype classes. A small
  `NADD_F` does not by itself exclude an additive effect of a marker that the
  grouping merges or splits; the scan can also compute `NADD_FS`, which
  conditions on the class dosages and on every single-marker dosage of the unit.

Degrees of freedom are the actual rank differences of the nested designs, from a
column-pivoted QR with relative tolerance 1e-8: for `A_joint` the numerator df
is rank([X0, dosages]) - rank(X0), where X0 is the null design, and similarly
for `CELL_F` and `NADD_F`; the denominator df is n - rank of the larger design.
With nk retained classes and no further dependency `A_joint` has nk - 1 df, but
that is not assumed.

The unit p-value is the equal-weight ACAT of the four. The anchor p-value is
the ACAT of all unit p-values, with every scale given total weight 1 split
equally among its units, so that no scale dominates because it has many
units. Every anchor contributes exactly one p-value, which is what the
Bonferroni correction below counts.

## Thresholds

The default is Bonferroni over the number of anchors. For a genome-wide
empirical threshold, permute the phenotype, rescan, and take the 5% quantile of
the per-permutation minimum anchor p-value.

## Options worth knowing

| option | meaning |
| --- | --- |
| `--transform {int,none}` | rank-based inverse normal transform of the phenotype; on by default, skipped for binary traits |
| `--anchor-min` | minimum markers per anchor; controls the anchor count and thus the correction |
| `--ld-dmin`, `--ld-maxspan` | LD block definition (D-prime cut, maximum span in bp) |
| `--minf`, `--minc` | haplotype frequency and diplotype cell size floors |
| `--kernel-power {1,auto}` | `1` (default) is the standard GRM; `auto` selects a spectral power of the GRM by REML and is a sensitivity analysis, not a default |

## Self-test

```bash
hapdc selftest
```

Runs the same path as `hapdc ms` (LD blocks, anchors, the combined statistic D)
on founder haplotypes with a known architecture and checks that (1) an additive
and (2) a diplotype-specific effect are each detected at the causal anchor
below the Bonferroni level; (3) over 200 null phenotypes the family-wise rate at
that level is at most 0.10 (expected 0.05); (4) a phenotype gets the same
p-values scanned alone or in a batch; (5) recoding the alleles of every marker
and (6) reordering the samples leave every anchor p-value unchanged. It exits
non-zero if any check fails. Takes about ten seconds. The test suite in `tests/`
covers the remaining code paths (`pytest tests`).

## Example

```bash
cd examples
bash run_example.sh
```

`make_toy.py` writes a deterministic toy panel (400 individuals, two
chromosomes of 300 phased markers) with one planted recessive haplotype effect
on chromosome 1 (`toy/truth.txt`); the script then runs `prep` and `ms` and
prints the top anchors.

## Citation

If you use hapdc, please cite the accompanying paper (in preparation; see
`CITATION.cff`).

## License

MIT. See `LICENSE`.
