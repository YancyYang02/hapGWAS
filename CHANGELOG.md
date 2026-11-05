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
- Zenodo DOI is **not** included. Enable the GitHub–Zenodo integration on this repository and archive the `v0.4.1` GitHub Release so Zenodo can mint the DOI.

### Notes
- Release packaging checklist: `appnote_hardening/ITEM6_RELEASE_zh.md`.
- Public source: https://github.com/YancyYang02/hapGWAS (tag `v0.4.1`). Do not invent a Zenodo DOI.
