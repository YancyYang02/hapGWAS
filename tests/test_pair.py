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
from hapgwas.pair import pair_indicator_score, trans_compound_indicator


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


def _pair_p(hap1, hap2, y, K=None) -> float:
    method = "lmm" if K is not None else "ols"
    null = fit_null_lmm(y, K=K, method=method)
    dosage = build_dosage(
        haplotypes_from_phased_snps(hap1, hap2),
        min_freq=0.01,
        min_count=5,
        collapse_mode="freq",
        drop_ref=True,
    )
    indicator = trans_compound_indicator(hap1, hap2)
    _chisq, _df, p = pair_indicator_score(null, dosage.dosage, indicator)
    return float(p)


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


def test_run_genome_records_first_two_sites_and_omits_column_by_default():
    rng = np.random.default_rng(3)
    hap1, hap2, y, pair = _two_site("phase_het", rng, n=240, beta=1.1)
    extra_h1 = rng.integers(0, 2, size=(240, 2), dtype=np.int8)
    extra_h2 = rng.integers(0, 2, size=(240, 2), dtype=np.int8)
    hap1 = np.column_stack([hap1, extra_h1])
    hap2 = np.column_stack([hap2, extra_h2])
    null = fit_null_lmm(y, method="ols")
    blocks = fixed_windows(4, window=4, chrom="1")
    plain, _alleles = run_genome(null, blocks, hap1, hap2, apply_fdr=False)
    assert "pair_p" not in plain.columns
    scored, _alleles = run_genome(null, blocks, hap1, hap2, apply_fdr=True, compute_pair=True)
    row = scored.iloc[0]
    assert int(row["pair_snp_i"]) == 0
    assert int(row["pair_snp_j"]) == 1
    assert int(row["pair_n"]) == int(pair.sum())
    assert np.isfinite(row["pair_p"])
    assert row["pair_p"] < 1e-4
    assert "pair_q" in scored.columns
    assert np.isfinite(row["pair_q"])


def test_run_help_lists_pair(capsys):
    with pytest.raises(SystemExit) as exc:
        main(["run", "--help"])
    assert exc.value.code == 0
    assert "--pair" in capsys.readouterr().out


def test_cli_fixture_writes_finite_pair_p(tmp_path):
    rng = np.random.default_rng(7)
    n = 240
    hap1, hap2, y, _pair = _two_site("phase_het", rng, n=n, beta=1.15)
    extra_h1 = rng.integers(0, 2, size=(n, 2), dtype=np.int8)
    extra_h2 = rng.integers(0, 2, size=(n, 2), dtype=np.int8)
    hap1 = np.column_stack([hap1, extra_h1])
    hap2 = np.column_stack([hap2, extra_h2])
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
            "2",
            "--step",
            "2",
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
            "2",
            "--step",
            "2",
            "--out",
            str(plain),
        ]
    )
    assert code == 0
    plain_blocks = pd.read_csv(plain / "assoc_blocks.tsv", sep="\t")
    assert "pair_p" not in plain_blocks.columns
