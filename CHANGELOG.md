# Changelog

All notable changes to hapGWAS are documented in this file.

## [0.7.3] — 2026-11-06

### Changed
- `--pair` now scores every consecutive SNP pair in the block. Each trans compound-het indicator is residualized on the block's additive haplotype dosages in the same whitened space as `joint_p` and `beyond_p`. `pair_p` is the ACAT of those 1-df p-values. A two-SNP block still has one pair, so the ACAT is that pair. The first-pair-only p-value is no longer the primary column.
- Debug columns are `pair_n` (how many pairs entered the ACAT), `pair_minp` (smallest pair p-value, not Bonferroni-adjusted), and `pair_snp_i` / `pair_snp_j` (that pair). The single-pair `pair_chisq` and `pair_df` columns are not written.
- Package version is **0.7.3**. `CITATION.cff` version is 0.7.3. The `doi` and `identifiers` fields are unchanged and still point at the Zenodo deposit for **v0.4.1**.

### Notes
- Experimental. A simulation with the phase effect on a pair other than the first SNP pair gained 55–59 percentage points of rejection rate over the first pair alone. Null and additive-only stayed near the nominal level in that check.
- One broiler463 BW6 spectral LMM (GCTA GRM, h² ≈ 0.19, n = 463, first 20,000 phased SNPs, 9,999 windows of 4 SNPs with step 2, Bonferroni 5.0×10⁻⁶) compared the two scores on the same run. First pair only: 18 windows below 1×10⁻⁴, 16 Bonferroni hits, λ_GC 0.76. Consecutive-pair ACAT: 61 and 45, λ_GC 0.67. The in-window SNP joint score was 25 and 16, λ_GC 0.96. Window `40:w600` still passes (ACAT 6.68×10⁻⁹, first pair 3.34×10⁻⁹, SNP score 0.071). A Bonferroni minimum of the pair p-values matched the 61 / 45 counts with λ_GC 0.36 and is not the shipped column. This panel is not a genome-wide scan. The flag does not claim a new leave-one-chromosome-out method.

## [0.7.2] — 2026-11-06

### Added
- Optional `--pair` on `hapgwas run`. The first cut of this flag scored only the first two SNPs. **0.7.3 replaces that primary column** with an ACAT over every consecutive pair; the notes below describe the scan that motivated adding the flag, not the 0.7.3 statistic. The score uses the same whitened phenotype as `joint_p` and `beyond_p`. It is not part of `omni_p`. SACJ stays out of the CLI.
- `tests/test_pair.py`: the indicator rejects more often when the phenotype is that phase configuration than under a null or an additive haplotype effect, under OLS and under an LMM.

### Notes
- Experimental. The flag follows a broiler463 BW6 spectral-LMM scan of the first 20,000 phased SNPs (4-SNP windows, step 2, 9,999 windows, GCTA GRM h² ≈ 0.19). The indicator had 10 Bonferroni hits (threshold 5.0×10⁻⁶) and 15 windows below 1×10⁻⁴, and all 15 of those had a two-SNP score above 1×10⁻². λ_GC was 0.45. Window `40:w600` was 1.40×10⁻⁶ on both the 5,000-SNP and the 20,000-SNP LMM scans. Joint and beyond companions of the same indicator were not kept (210 and 184 Bonferroni windows). This release does not claim a new leave-one-chromosome-out method.
- **DOI.** `CITATION.cff` version is 0.7.2. The `doi` and `identifiers` fields are unchanged. They still point at the Zenodo deposit for **v0.4.1** (version DOI `10.5281/zenodo.23984394`, concept DOI `10.5281/zenodo.23984393`, record https://zenodo.org/records/23984394).

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
