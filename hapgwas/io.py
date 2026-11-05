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


def load_grm(path: str | Path, sample_ids: list[str] | None = None) -> np.ndarray:
    """Load a square kinship/GRM matrix and optionally reorder to sample_ids.

    Accepted formats:
      - .npy binary matrix
      - whitespace text matrix (n x n), optional leading ID column/row
      - GCTA-style .grm.bin + .grm.id (pass path without suffix or .grm.bin)
      - npz with key 'K' or 'grm'
    """
    path = Path(path)
    ids_from_file: list[str] | None = None
    K: np.ndarray | None = None

    # GCTA binary GRM: stem.grm.bin + stem.grm.id (or stem without suffix)
    stem = path
    if path.name.endswith(".grm.bin"):
        stem = path.with_suffix("").with_suffix("")  # drop .bin then .grm
    elif path.suffix == ".id" and path.name.endswith(".grm.id"):
        stem = Path(str(path)[: -len(".grm.id")])
    grm_bin = Path(str(stem) + ".grm.bin")
    grm_id = Path(str(stem) + ".grm.id")
    if grm_bin.exists() and grm_id.exists():
        ids_from_file = [ln.split()[1] if len(ln.split()) > 1 else ln.split()[0]
                         for ln in grm_id.read_text().splitlines() if ln.strip()]
        n = len(ids_from_file)
        # lower triangle including diagonal, float32
        n_elem = n * (n + 1) // 2
        raw = np.fromfile(grm_bin, dtype=np.float32)
        if raw.size < n_elem:
            raise ValueError(f"GRM bin too short: {raw.size} < {n_elem}")
        K = np.zeros((n, n), dtype=float)
        idx = 0
        for i in range(n):
            for j in range(i + 1):
                K[i, j] = raw[idx]
                K[j, i] = raw[idx]
                idx += 1
    elif path.suffix == ".npy":
        K = np.load(path)
    elif path.suffix == ".npz":
        z = np.load(path, allow_pickle=True)
        key = "K" if "K" in z.files else ("grm" if "grm" in z.files else z.files[0])
        K = np.asarray(z[key], dtype=float)
        if "ids" in z.files:
            ids_from_file = [str(x) for x in z["ids"]]
    else:
        # text: try ID-aware then plain
        rows = []
        with path.open() as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                rows.append(line.split())
        if not rows:
            raise ValueError(f"Empty GRM file {path}")
        # If first cell non-numeric -> row IDs
        def _is_num(s: str) -> bool:
            try:
                float(s)
                return True
            except ValueError:
                return False
        if not _is_num(rows[0][0]):
            has_header = all(not _is_num(x) for x in rows[0])
            body = rows[1:] if has_header else rows
            # leading ID column?
            if body and not _is_num(body[0][0]):
                ids_from_file = [r[0] for r in body]
                data = [r[1:] for r in body]
            else:
                ids_from_file = rows[0] if has_header else None
                data = body
            K = np.asarray([[float(x) for x in r] for r in data], dtype=float)
        else:
            K = np.asarray([[float(x) for x in r] for r in rows], dtype=float)

    if K is None:
        raise ValueError(f"Could not load GRM from {path}")
    K = np.asarray(K, dtype=float)
    if K.ndim != 2 or K.shape[0] != K.shape[1]:
        raise ValueError(f"GRM must be square, got {K.shape}")

    if sample_ids is not None:
        if ids_from_file is None:
            if K.shape[0] != len(sample_ids):
                raise ValueError(
                    f"GRM n={K.shape[0]} != n_samples={len(sample_ids)} and no IDs to reorder"
                )
        else:
            pos = {s: i for i, s in enumerate(ids_from_file)}
            missing = [s for s in sample_ids if s not in pos]
            if missing:
                raise ValueError(f"GRM missing {len(missing)} sample IDs e.g. {missing[:5]}")
            idx = [pos[s] for s in sample_ids]
            K = K[np.ix_(idx, idx)]
    return 0.5 * (K + K.T)
