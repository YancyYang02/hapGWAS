"""Input / output helpers for hapGWAS."""

from __future__ import annotations

from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

# Re-export genotype API for a single import surface
from hapgwas.genotype import (  # noqa: F401
    PhasedGeno,
    load_geno,
    load_npz,
    read_oxford_haps,
    read_phased_vcf,
    read_plink_pedmap,
    save_npz,
)


def read_phenotype(
    path: str | Path,
    id_col: str = "IID",
    pheno_col: str = "PHENO",
) -> pd.DataFrame:
    """Read a whitespace or CSV phenotype table with IID and PHENO columns."""
    path = Path(path)
    sep = "," if path.suffix.lower() == ".csv" else r"\s+"
    df = pd.read_csv(path, sep=sep, engine="python")
    cols = {c.lower(): c for c in df.columns}
    if id_col not in df.columns and id_col.lower() in cols:
        id_col = cols[id_col.lower()]
    if pheno_col not in df.columns and pheno_col.lower() in cols:
        pheno_col = cols[pheno_col.lower()]
    if id_col not in df.columns or pheno_col not in df.columns:
        raise ValueError(f"Need columns {id_col!r} and {pheno_col!r} in {path}")
    out = df[[id_col, pheno_col]].copy()
    out.columns = ["IID", "PHENO"]
    out["IID"] = out["IID"].astype(str)
    return out.dropna(subset=["PHENO"]).reset_index(drop=True)


def read_covariates(path: str | Path | None, id_col: str = "IID") -> pd.DataFrame | None:
    if path is None:
        return None
    path = Path(path)
    sep = "," if path.suffix.lower() == ".csv" else r"\s+"
    df = pd.read_csv(path, sep=sep, engine="python")
    cols = {c.lower(): c for c in df.columns}
    if id_col not in df.columns and id_col.lower() in cols:
        id_col = cols[id_col.lower()]
    df = df.rename(columns={id_col: "IID"})
    df["IID"] = df["IID"].astype(str)
    return df


def read_blocks(path: str | Path) -> pd.DataFrame:
    """Block file reader (generic TSV/CSV)."""
    path = Path(path)
    sep = "," if path.suffix.lower() == ".csv" else r"\s+"
    return pd.read_csv(path, sep=sep, engine="python")


def write_assoc(df: pd.DataFrame, path: str | Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False, float_format="%.6g")


def align_ids(
    pheno: pd.DataFrame,
    sample_ids: Iterable[str],
    cov: pd.DataFrame | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray | None, list[str]]:
    """Return y, keep_mask into genotype order, X covariates, kept IIDs."""
    sample_ids = [str(s) for s in sample_ids]
    pheno_map = dict(zip(pheno["IID"], pheno["PHENO"]))
    keep = [i for i, s in enumerate(sample_ids) if s in pheno_map]
    ids = [sample_ids[i] for i in keep]
    y = np.asarray([pheno_map[s] for s in ids], dtype=float)
    if cov is not None:
        cov_idx = cov.set_index("IID")
        missing = [s for s in ids if s not in cov_idx.index]
        if missing:
            raise ValueError(f"{len(missing)} samples missing covariates, e.g. {missing[:3]}")
        cov_mat = cov_idx.loc[ids].select_dtypes(include=[np.number]).to_numpy(dtype=float)
        X = np.column_stack([np.ones(len(ids)), cov_mat])
    else:
        X = np.ones((len(ids), 1), dtype=float)
    return y, np.asarray(keep, dtype=int), X, ids


def write_phased_vcf(g, path: str | Path) -> None:
    """Write a minimal phased biallelic VCF from PhasedGeno (for smoke tests)."""
    path = Path(path)
    chrom = g.chrom if g.chrom is not None else np.asarray(["1"] * g.n_snps, dtype=object)
    pos = g.positions if g.positions is not None else np.arange(1, g.n_snps + 1, dtype=np.int64)
    snp_ids = g.snp_ids if g.snp_ids is not None else np.asarray([f"snp{i}" for i in range(g.n_snps)], dtype=object)
    if g.alleles is not None:
        refs = g.alleles[:, 0]
        alts = g.alleles[:, 1]
    else:
        refs = np.asarray(["A"] * g.n_snps)
        alts = np.asarray(["G"] * g.n_snps)

    with open(path, "wt") as fh:
        fh.write("##fileformat=VCFv4.2\n")
        fh.write("##source=hapGWAS\n")
        fh.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(g.ids) + "\n")
        for j in range(g.n_snps):
            gts = [f"{int(g.hap1[i, j])}|{int(g.hap2[i, j])}" for i in range(g.n_samples)]
            fh.write(
                f"{chrom[j]}\t{int(pos[j])}\t{snp_ids[j]}\t{refs[j]}\t{alts[j]}\t.\tPASS\t.\tGT\t"
                + "\t".join(gts)
                + "\n"
            )
