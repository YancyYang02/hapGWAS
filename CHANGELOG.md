# Changelog

All notable changes to hapGWAS are documented in this file.

## [0.7.1] — 2026-11-06

### Added
- Leave-one-window-out kinship (`--lowo`, `--lowo-flank`): one eigen-decomposition of a VanRaden GRM, then a Woodbury downdate that drops the scanned window and a flank. The column is `joint_lowo_p`. This follows proximal-contamination / LOCO practice. It is the same kinship estimator, and it is not a claim to have introduced LOCO or regional heritability mapping.
- `--kin-geno` / `--kin-map` / `--kin-sample` to build kinship from a file that leaves the test markers out, and `--grm auto|vanraden|loco` when one genotype file spans more than one chromosome.
- GCTA `--make-grm` reader (`--gcta-grm`, `.grm.id` + `.grm.bin`, aligned on IID). A binary GRM has no SNP factors, so it is not combined with `--lowo`.
- Default block ranking `omni_p`: ACAT (Liu & Xie, JASA 2020) of the joint haplotype score, the allele-marginal scores, and the SNP scores in the window. Locus collapse uses `omni_p`.
- Exact block likelihood-ratio columns `joint_lrt_p` and `beyond_lrt_p` (profile ML on the shared eigen-decomposition of K).
- Secondary `beyond_p`: joint score of haplotype dosages after residualizing the additive SNP dosages in the same block. It asks whether signal remains after those SNPs. It is not the default rank.
- `tree_p`: Brownian kernel on a UPGMA tree of haplotype Hamming distances. Secondary column only. It is not part of `omni_p`.
- Summary JSON under `sim_results/` for the omnibus, confirmatory, LRT, LOWO, tagged-flank, GCTA-GRM, and tree checks (`omni_bench.json`, `paper_sim.json`, `lrt_bench.json`, `lowo_bench.json`, `lowo_tagged.json`, `lowo_real_lowo_real.json`, `gcta_hap.json`, `tree_bench.json`, and the smaller summaries next to them). Package scripts: `scripts/bench_omni.py`, `scripts/bench_paper.py`, `scripts/bench_lrt.py`, `scripts/bench_lowo.py`, `scripts/bench_tree.py`, `scripts/run_lowo_real.py`, `scripts/run_lowo_tagged.py`, `scripts/run_gcta_hap.py`.

### Changed
- Package version is **0.7.1** in `pyproject.toml`, `VERSION`, `hapgwas.__version__`, and `CITATION.cff`.
- README describes the 0.7.1 scan: `omni_p` ranks blocks; beyond and tree stay secondary; broiler GCTA-GRM numbers are reported with the calibration limits in the tables (including traits pinned to the heritability grid, and the chromosome-40 λ values).
- ACJ remains a written column (`acj_p`). It is not the default rank. SACJ and `hapgwas run-multi` from 0.4.1 are not in the 0.7.1 CLI. Historical scripts that call them (`scripts/run_sacj_sim.py`, `scripts/run_chr4_sacj_hpc.py`, `scripts/run_broiler_multi_hpc.py`) are records of those earlier analyses.

### Real-data notes
- Chromosome field `40`, GCTA GRM from 8,251,083 autosomal SNPs, n = 463 after the phenotype merge, windows of 4 SNPs with step 2. Counts and λ are in the README and in `sim_results/gcta_hap.json`. BW6 is the calibration row (λ near 1 on the 1500-SNP panel). AFWR's λ is not a calibration result. A 40-SNP downdate of that dense GRM does not move the AFWR lead; `--lowo` matters when the kinship still contains the local tags.
- Chip check of the Woodbury update: one window matches a GRM rebuilt without those SNPs to about 1e-15, and across 749 windows the median |Δ log10 p| versus rebuilt kinship is about 2.6×10⁻⁶.
- GHap 4.0.0.14 on the same 749 windows, without this GRM, calls more windows. That comparison is the effect of the relatedness model on this chromosome, not a claim that one program is more powerful at every locus.

### Notes
- **DOI.** `CITATION.cff` version is 0.7.1. The `doi` and `identifiers` fields are unchanged from v0.4.1. They still point at the Zenodo deposit for **v0.4.1** (version DOI `10.5281/zenodo.23984394`, concept DOI `10.5281/zenodo.23984393`, record https://zenodo.org/records/23984394). A new tag and a new Zenodo version will be cut later. Until that re-deposit, do not treat the current DOI as a 0.7.1 record.
- Private-host SSH helpers from the working fork were not added to this public tree. Passwords were never in those files; the host layout stays off the public branch.
- Family structure that the GRM does not capture can still inflate type I error. Trust a scan when `run_meta.json` shows raw-p λ_GC near 1.

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
