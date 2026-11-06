"""PAIR indicator: phase-discordant signal, additive null, CLI column."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from hapgwas.blocks import fixed_windows
from hapgwas.cli import main
from hapgwas.dosage import build_dosage, haplotypes_from_phased_snps
from hapgwas.genotype import PhasedGeno, save_npz
from hapgwas.joint_test import run_genome
from hapgwas.lmm import fit_null_lmm
from hapgwas.pair import pair_indicator_score, score_consecutive_pairs, trans_compound_indicator


def _two_site(mode: str, rng: np.random.Generator, n: int = 500, beta: float = 0.7):
    """Two phased SNPs. Matches the trial codes: 0=ab, 1=aB, 2=Ab, 3=AB."""
    p_a, p_b = 0.3, 0.3
    probs = np.array(
        [(1 - p_a) * (1 - p_b), (1 - p_a) * p_b, p_a * (1 - p_b), p_a * p_b],
        dtype=float,
    )
    h1 = rng.choice(4, size=n, p=probs)
    h2 = rng.choice(4, size=n, p=probs)

    def allele(h, which):
        if which == 0:
            return (h >= 2).astype(np.int8)
        return ((h == 1) | (h == 3)).astype(np.int8)

    hap1 = np.column_stack([allele(h1, 0), allele(h1, 1)])
    hap2 = np.column_stack([allele(h2, 0), allele(h2, 1)])
    ab = ((h1 == 3).astype(float) + (h2 == 3).astype(float))
    pair = (((h1 == 2) & (h2 == 1)) | ((h1 == 1) & (h2 == 2))).astype(float)
    noise = rng.normal(0.0, 1.0, n)
    if mode == "null":
        y = noise
    elif mode == "additive":
        y = beta * ab + noise
    elif mode == "phase_het":
        y = beta * pair + noise
    else:
        raise ValueError(mode)
    return hap1, hap2, y, pair


def _dosage(hap1, hap2):
    return build_dosage(
        haplotypes_from_phased_snps(hap1, hap2),
        min_freq=0.0,
        min_count=1,
        collapse_mode="none",
        drop_ref=True,
    )


def _pair_summary(hap1, hap2, y, K=None):
    method = "lmm" if K is not None else "ols"
    null = fit_null_lmm(y, K=K, method=method)
    dosage = _dosage(hap1, hap2)
    return null, dosage.dosage, score_consecutive_pairs(null, dosage.dosage, hap1, hap2)


def _pair_p(hap1, hap2, y, K=None) -> float:
    _null, _dosage, scored = _pair_summary(hap1, hap2, y, K=K)
    return float(scored.pair_p)


def _first_pair_p(hap1, hap2, y) -> float:
    null = fit_null_lmm(y, method="ols")
    dosage = _dosage(hap1, hap2)
    _chisq, _df, p = pair_indicator_score(null, dosage.dosage, trans_compound_indicator(hap1, hap2, 0, 1))
    return float(p)


def _four_site(mode: str, rng: np.random.Generator, n: int = 420, beta: float = 0.9):
    """4 SNPs. The phase configuration sits on the last consecutive pair, not the first."""
    front1, front2, _y, _pair = _two_site("null", rng, n=n, beta=0.0)
    back1, back2, y, pair = _two_site(mode, rng, n=n, beta=beta)
    hap1 = np.column_stack([front1, back1])
    hap2 = np.column_stack([front2, back2])
    return hap1, hap2, y, pair


def _reject_rate(mode: str, rng: np.random.Generator, reps: int, n: int, beta: float, lmm: bool) -> float:
    hits = 0
    finite = 0
    for _ in range(reps):
        hap1, hap2, y, _pair = _two_site(mode, rng, n=n, beta=beta)
        K = None
        if lmm:
            factors = rng.normal(size=(n, 10))
            raw = factors @ factors.T / 10.0
            scale = float(np.mean(np.diag(raw)))
            K = raw / scale
            genetic = np.sqrt(0.35 / scale) * (factors @ rng.normal(size=10)) / np.sqrt(10.0)
            y = y + genetic
        p = _pair_p(hap1, hap2, y, K=K)
        if np.isfinite(p):
            finite += 1
            hits += int(p < 0.05)
    assert finite == reps
    return hits / reps


def test_indicator_is_trans_only():
    hap1 = np.array([[1, 0], [1, 1], [1, 0], [0, 0]], dtype=np.int8)
    hap2 = np.array([[0, 1], [0, 0], [1, 0], [0, 1]], dtype=np.int8)
    got = trans_compound_indicator(hap1, hap2)
    # row0 Ab/aB trans; row1 cis AB/ab; row2 two copies of Ab; row3 not a compound het
    np.testing.assert_array_equal(got, np.array([1.0, 0.0, 0.0, 0.0]))


def test_pair_rejects_more_under_phase_het_than_additive_or_null():
    rng = np.random.default_rng(9)
    n, beta, reps = 480, 0.75, 80
    rates = {mode: _reject_rate(mode, rng, reps, n, beta, lmm=False) for mode in ("null", "additive", "phase_het")}
    assert rates["phase_het"] > 0.40, rates
    assert rates["phase_het"] > rates["additive"] + 0.20, rates
    assert rates["phase_het"] > rates["null"] + 0.20, rates
    assert rates["additive"] < 0.15, rates
    assert rates["null"] < 0.15, rates


def test_pair_lmm_keeps_the_same_ordering():
    rng = np.random.default_rng(11)
    n, beta, reps = 220, 0.90, 36
    rates = {mode: _reject_rate(mode, rng, reps, n, beta, lmm=True) for mode in ("null", "additive", "phase_het")}
    assert rates["phase_het"] > 0.30, rates
    assert rates["phase_het"] > rates["additive"] + 0.15, rates
    assert rates["phase_het"] > rates["null"] + 0.15, rates
    assert rates["additive"] < 0.20, rates
    assert rates["null"] < 0.20, rates


def test_two_snp_window_matches_the_single_pair():
    rng = np.random.default_rng(4)
    hap1, hap2, y, _pair = _two_site("phase_het", rng, n=180, beta=0.8)
    null, dosage, scored = _pair_summary(hap1, hap2, y)
    _chisq, _df, single = pair_indicator_score(null, dosage, trans_compound_indicator(hap1, hap2, 0, 1))
    assert scored.pair_n == 1
    assert np.isfinite(scored.pair_p) and np.isfinite(single)
    assert abs(scored.pair_p - single) < 1e-8
    assert scored.pair_snp_i == 0 and scored.pair_snp_j == 1


def test_later_pair_is_recovered_by_acat_and_missed_by_the_first_pair():
    rng = np.random.default_rng(21)
    n, beta, reps = 400, 0.95, 40
    rates = {mode: {"v2": 0, "v1": 0, "finite_v2": 0, "finite_v1": 0} for mode in ("null", "additive", "phase_het")}
    for mode in rates:
        for _ in range(reps):
            hap1, hap2, y, _pair = _four_site(mode, rng, n=n, beta=beta)
            p_v2 = _pair_p(hap1, hap2, y)
            p_v1 = _first_pair_p(hap1, hap2, y)
            if np.isfinite(p_v2):
                rates[mode]["finite_v2"] += 1
                rates[mode]["v2"] += int(p_v2 < 0.05)
            if np.isfinite(p_v1):
                rates[mode]["finite_v1"] += 1
                rates[mode]["v1"] += int(p_v1 < 0.05)
    for mode in rates:
        assert rates[mode]["finite_v2"] == reps, rates
        assert rates[mode]["finite_v1"] > reps * 0.8, rates
        rates[mode]["v2"] /= reps
        rates[mode]["v1"] /= rates[mode]["finite_v1"]
    phase = rates["phase_het"]
    assert phase["v2"] > 0.45, rates
    assert phase["v2"] > phase["v1"] + 0.35, rates
    assert phase["v1"] < 0.20, rates
    assert rates["null"]["v2"] < 0.15, rates
    assert rates["additive"]["v2"] < 0.20, rates

    hap1, hap2, y, _pair = _four_site("phase_het", rng, n=360, beta=1.25)
    _null, _dosage, scored = _pair_summary(hap1, hap2, y)
    first = _first_pair_p(hap1, hap2, y)
    assert np.isfinite(scored.pair_p) and scored.pair_p < 1e-4
    assert scored.pair_n == 3
    assert np.isfinite(first) and first > scored.pair_p * 50
    assert (scored.pair_snp_i, scored.pair_snp_j) != (0, 1)


def test_run_genome_uses_every_consecutive_pair_and_omits_column_by_default():
    rng = np.random.default_rng(3)
    hap1, hap2, y, _pair = _four_site("phase_het", rng, n=280, beta=1.2)
    null = fit_null_lmm(y, method="ols")
    blocks = fixed_windows(4, window=4, chrom="1")
    plain, _alleles = run_genome(null, blocks, hap1, hap2, apply_fdr=False, min_freq=0.0, min_count=1, collapse_mode="none")
    assert "pair_p" not in plain.columns
    scored, _alleles = run_genome(
        null, blocks, hap1, hap2, apply_fdr=True, compute_pair=True, min_freq=0.0, min_count=1, collapse_mode="none"
    )
    row = scored.iloc[0]
    assert int(row["pair_n"]) == 3
    assert (int(row["pair_snp_i"]), int(row["pair_snp_j"])) != (0, 1)
    assert np.isfinite(row["pair_p"]) and row["pair_p"] < 1e-4
    assert np.isfinite(row["pair_minp"]) and row["pair_minp"] <= row["pair_p"] + 1e-8
    assert "pair_q" in scored.columns
    assert np.isfinite(row["pair_q"])
    assert "pair_chisq" not in scored.columns


def test_run_help_lists_pair(capsys):
    with pytest.raises(SystemExit) as exc:
        main(["run", "--help"])
    assert exc.value.code == 0
    text = capsys.readouterr().out
    assert "--pair" in text
    assert "ACAT" in text
    assert "consecutive" in text


def test_cli_fixture_writes_finite_pair_p(tmp_path):
    rng = np.random.default_rng(7)
    n = 240
    _front1, _front2, _noise, _pair = _two_site("null", rng, n=n, beta=0.0)
    back1, back2, y, _causal = _two_site("phase_het", rng, n=n, beta=1.2)
    hap1 = np.column_stack([_front1, back1])
    hap2 = np.column_stack([_front2, back2])
    ids = [f"S{i:04d}" for i in range(n)]
    geno = PhasedGeno(
        hap1=hap1.astype(np.int8),
        hap2=hap2.astype(np.int8),
        ids=ids,
        chrom=np.asarray(["1", "1", "1", "1"], dtype=object),
        positions=np.arange(1, 5) * 1000,
        snp_ids=np.asarray(["s0", "s1", "s2", "s3"], dtype=object),
    )
    geno_path = tmp_path / "phased.npz"
    pheno_path = tmp_path / "pheno.tsv"
    save_npz(geno, geno_path)
    pd.DataFrame({"IID": ids, "PHENO": y}).to_csv(pheno_path, sep="\t", index=False)

    out = tmp_path / "pair_out"
    code = main(
        [
            "run",
            "--geno",
            str(geno_path),
            "--pheno",
            str(pheno_path),
            "--ols",
            "--pair",
            "--window",
            "4",
            "--step",
            "4",
            "--min-count",
            "1",
            "--min-freq",
            "0",
            "--collapse",
            "none",
            "--chrom",
            "1",
            "--out",
            str(out),
        ]
    )
    assert code == 0
    blocks = pd.read_csv(out / "assoc_blocks.tsv", sep="\t")
    loci = pd.read_csv(out / "assoc_loci.tsv", sep="\t")
    assert "pair_p" in blocks.columns
    assert np.isfinite(blocks.loc[0, "pair_p"])
    assert blocks.loc[0, "pair_p"] < 1e-3
    assert int(blocks.loc[0, "pair_n"]) == 3
    assert np.isfinite(blocks.loc[0, "pair_minp"])
    assert "pair_p" in loci.columns
    assert np.isfinite(loci["pair_p"]).any()

    plain = tmp_path / "plain_out"
    code = main(
        [
            "run",
            "--geno",
            str(geno_path),
            "--pheno",
            str(pheno_path),
            "--ols",
            "--window",
            "4",
            "--collapse",
            "none",
            "--out",
            str(plain),
        ]
    )
    assert code == 0
    plain_blocks = pd.read_csv(plain / "assoc_blocks.tsv", sep="\t")
    assert "pair_p" not in plain_blocks.columns
