# Changelog

All notable changes to hapGWAS are documented in this file.

## [0.4.1] — 2026-11-05

### Added
- CLI `--grm` / `--kinship` for spectral LMM using an external GRM (GCTA `.grm.bin`+`.id`, npy/npz, or text square matrix), aligned to sample IDs.
- App Note hardening notes under `appnote_hardening/` (LMM λ, official GHap head-to-head, SNP genome-wide Bonferroni, locus collapse, conditional analysis).

### Changed
- Package / `VERSION` / `CITATION.cff` aligned to **0.4.1**.
- Documentation: SACJ demoted to exploratory / supplement; App Note main tables use **joint / ACJ / official GHap allele / SNP**.
- README title bumped to v0.4.1. Install instructions use https://github.com/YancyYang02/hapGWAS.
- Zenodo version DOI for this release: https://doi.org/10.5281/zenodo.23984394 ([record](https://zenodo.org/records/23984394)).
- Zenodo concept DOI (all versions): https://doi.org/10.5281/zenodo.23984393.

### Notes
- Release packaging checklist: `appnote_hardening/ITEM6_RELEASE_zh.md`.
- Public source: https://github.com/YancyYang02/hapGWAS (tag `v0.4.1`).
