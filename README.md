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
| `P_nonadditive` | scale-weighted ACAT of the `NADD_F` p-values in the anchor |
| `sig_bonf` | 1 if `P` < 0.05 / number of anchors |

`<trait>.summary.json` records the sample size, number of anchors, Bonferroni
threshold, number of significant anchors, the genomic-control lambda of the
anchor p-values, and the options used. The lambda is also printed during the scan; a value well below the
null reference (0.61-0.68) triggers a warning that the Bonferroni threshold
should be replaced by a permutation threshold that preserves the sample structure.

## Inputs

* **Genotypes** -- a VCF, phased (unphased genotypes are accepted but every
  window is then built from phase-unknown pairs and power drops). Genotypes are
  read as 0/1 haplotypes, two rows per individual.
* **Relationship matrix** -- either a GCTA `.grm` prefix, or `auto` to compute
  a VanRaden GRM from the VCF itself (`--grm-thin` markers are used).
* **Phenotypes** -- a table whose first column is the individual ID and whose
  remaining columns are traits, one column per trait, missing values coded
  `-9` by default (`--miss`). Numeric covariates, for example leading principal
  components or batch indicators, go in a separate file passed to `hapdc ms
  --covar`, one header line, IID followed by numeric columns.

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
distance. Individuals are grouped into diplotype cells; cells below `--minc`
individuals are pooled. Four statistics are computed per unit from the same
mixed-model fit:

* `A_joint` -- joint F test of the additive haplotype dosages;
* `CELL_F` -- joint F test of the diplotype-cell indicators (dominance,
  overdominance, recessive and compound-heterozygous effects);
* `SACAT` -- Cauchy combination of the single-marker additive tests in the unit;
* `NADD_F` -- F test of the diplotype cells conditional on the additive
  haplotype dosages.

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

Builds founder haplotypes with a known architecture and prints p-values for an
additive, an overdominant and a null region. The additive region is detected by
the additive statistic and the overdominant region by the diplotype statistic;
the null region is not detected. Takes a few seconds.

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
