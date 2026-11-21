# Changelog

## 0.6.6 (2026-11-21)

Current state of the package.

- `Scan.testable` measures each marker's reference norm after removing the
  whitened rotated intercept `o = (U'1) * w`. The set of testable markers does
  not depend on allele coding, on the signs of the GRM eigenvectors or on
  the other markers, including when a covariate is nearly collinear with a
  marker.
- Before REML, the fixed-effect design `[intercept, covariates]` is checked
  for full column rank on the samples that are analysed
  (`core.check_fixed_effects`). Constant columns, duplicated columns and
  linear combinations of earlier columns are dropped and logged
  (`--covar-rank drop`, the default) or rejected with an explicit error
  (`--covar-rank error`). `summary.json` records `covariates_used`,
  `covariates_dropped` and `covar_rank_policy`. `reml_delta` raises
  `FixedEffectRankError` on a rank-deficient design.
- Provenance in `summary.json` is keyed by the full normalised path plus the
  MD5 of each input file, so inputs that share a basename are recorded
  separately. `merge` records the parameters and input hashes of every
  shard.
- There is no built-in reference range for `lambda_gc`. `--lambda-ref LO,HI`
  only reports whether the observed value falls in a range the user
  supplies.
- The README describes fractional haplotype-class dosages for tied
  candidates, F-test degrees of freedom computed as rank differences, the
  meaning of `NADD_F` and `NADD_FS`, and per-cohort covariate handling.
