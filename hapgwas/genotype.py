"""Phased genotype container and loaders (VCF / Oxford haps / PLINK ped-map / npz)."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd


@dataclass
class PhasedGeno:
    """Internal phased genotype representation used by hapGWAS."""

    hap1: np.ndarray          # (n, p) int8/int alleles coded 0/1 (biallelic) or small ints
    hap2: np.ndarray
    ids: list[str]
    chrom: np.ndarray | None = None       # (p,) object/str
    positions: np.ndarray | None = None   # (p,) int bp
    snp_ids: np.ndarray | None = None     # (p,)
    alleles: np.ndarray | None = None     # (p, 2) REF/ALT or A1/A2
    K: np.ndarray | None = None           # optional GRM (n, n)
    source: str = "unknown"
    meta: dict = field(default_factory=dict)

    @property
    def n_samples(self) -> int:
        return int(self.hap1.shape[0])

    @property
    def n_snps(self) -> int:
        return int(self.hap1.shape[1])

    def snp_dosage(self, index: int) -> np.ndarray:
        """0/1/2 dosage for one SNP from phased alleles."""
        return self.hap1[:, index].astype(float) + self.hap2[:, index].astype(float)

    def snp_dosage_matrix(self, start: int = 0, end: int | None = None) -> np.ndarray:
        end = self.n_snps if end is None else end
        return self.hap1[:, start:end].astype(float) + self.hap2[:, start:end].astype(float)

    def subset_samples(self, keep: np.ndarray) -> "PhasedGeno":
        keep = np.asarray(keep, dtype=int)
        K = self.K[np.ix_(keep, keep)] if self.K is not None else None
        return PhasedGeno(
            hap1=self.hap1[keep],
            hap2=self.hap2[keep],
            ids=[self.ids[i] for i in keep],
            chrom=self.chrom,
            positions=self.positions,
            snp_ids=self.snp_ids,
            alleles=self.alleles,
            K=K,
            source=self.source,
            meta=dict(self.meta),
        )

    def find_snp(self, key: str) -> int:
        """Resolve SNP by 0-based index string, snp_id, or chrom:pos."""
        key = str(key)
        if key.isdigit():
            idx = int(key)
            if idx < 0 or idx >= self.n_snps:
                raise ValueError(f"SNP index {idx} out of range 0..{self.n_snps-1}")
            return idx
        if self.snp_ids is not None:
            hits = np.where(self.snp_ids.astype(str) == key)[0]
            if len(hits) == 1:
                return int(hits[0])
            if len(hits) > 1:
                raise ValueError(f"SNP id {key!r} is not unique")
        if ":" in key and self.chrom is not None and self.positions is not None:
            chrom, _, rest = key.partition(":")
            pos = int(rest)
            hits = np.where((self.chrom.astype(str) == chrom) & (self.positions.astype(int) == pos))[0]
            if len(hits) == 1:
                return int(hits[0])
        raise ValueError(f"Could not resolve SNP {key!r}")


def save_npz(g: PhasedGeno, path: str | Path) -> None:
    path = Path(path)
    payload = {
        "hap1": g.hap1,
        "hap2": g.hap2,
        "ids": np.asarray(g.ids, dtype=object),
    }
    if g.chrom is not None:
        payload["chrom"] = g.chrom
    if g.positions is not None:
        payload["positions"] = g.positions
    if g.snp_ids is not None:
        payload["snp_ids"] = g.snp_ids
    if g.alleles is not None:
        payload["alleles"] = g.alleles
    if g.K is not None:
        payload["K"] = g.K
    np.savez_compressed(path, **payload)


def load_npz(path: str | Path) -> PhasedGeno:
    npz = np.load(path, allow_pickle=True)
    return PhasedGeno(
        hap1=npz["hap1"],
        hap2=npz["hap2"],
        ids=[str(x) for x in npz["ids"]],
        chrom=npz["chrom"] if "chrom" in npz.files else None,
        positions=npz["positions"] if "positions" in npz.files else None,
        snp_ids=npz["snp_ids"] if "snp_ids" in npz.files else None,
        alleles=npz["alleles"] if "alleles" in npz.files else None,
        K=npz["K"] if "K" in npz.files else None,
        source=str(path),
    )


def _allele_to_code(a: str, ref: str, alt: str) -> int:
    if a in {".", ""}:
        return -1
    if a == ref:
        return 0
    if a == alt:
        return 1
    # multi-allelic / mismatch
    return -1


def read_phased_vcf(
    path: str | Path,
    *,
    max_snps: int | None = None,
    require_phased: bool = True,
) -> PhasedGeno:
    """
    Read a biallelic phased VCF/VCF.gz-like text VCF (GT with '|').

    gzip: if path ends with .gz, use gzip.open. No pysam dependency.
    Unphased GT ('/') raises unless require_phased=False (then treats order as hap order with warning in meta).
    """
    import gzip

    path = Path(path)
    opener = gzip.open if path.suffix == ".gz" or str(path).endswith(".vcf.gz") else open
    sample_ids: list[str] = []
    chroms: list[str] = []
    positions: list[int] = []
    snp_ids: list[str] = []
    alleles: list[tuple[str, str]] = []
    # store columns as lists then stack — for small/medium VCFs
    hap1_cols: list[np.ndarray] = []
    hap2_cols: list[np.ndarray] = []
    n_unphased = 0

    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                parts = line.rstrip("\n").split("\t")
                sample_ids = parts[9:]
                continue
            if not sample_ids:
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 10:
                continue
            chrom, pos, vid, ref, alt, _qual, _filt, _info, fmt = parts[:9]
            alt0 = alt.split(",")[0]
            if "," in alt:
                # skip multi-allelic for MVP
                continue
            fmt_keys = fmt.split(":")
            if "GT" not in fmt_keys:
                continue
            gt_i = fmt_keys.index("GT")
            h1 = np.empty(len(sample_ids), dtype=np.int8)
            h2 = np.empty(len(sample_ids), dtype=np.int8)
            ok = True
            for s_i, cell in enumerate(parts[9:]):
                gt = cell.split(":")[gt_i]
                if "|" in gt:
                    a, b = gt.split("|", 1)
                elif "/" in gt:
                    n_unphased += 1
                    if require_phased:
                        ok = False
                        break
                    a, b = gt.split("/", 1)
                else:
                    ok = False
                    break
                # numeric alleles 0/1 relative to REF/ALT
                try:
                    ca = int(a) if a != "." else -1
                    cb = int(b) if b != "." else -1
                except ValueError:
                    ca = _allele_to_code(a, ref, alt0)
                    cb = _allele_to_code(b, ref, alt0)
                if ca < 0 or cb < 0 or ca > 1 or cb > 1:
                    ok = False
                    break
                h1[s_i] = ca
                h2[s_i] = cb
            if not ok:
                if require_phased and n_unphased:
                    raise ValueError(
                        f"Unphased GT found near {chrom}:{pos}. "
                        "Provide a phased VCF (Beagle/SHAPEIT) or set require_phased=False."
                    )
                continue
            chroms.append(str(chrom))
            positions.append(int(pos))
            snp_ids.append(vid if vid != "." else f"{chrom}:{pos}")
            alleles.append((ref, alt0))
            hap1_cols.append(h1)
            hap2_cols.append(h2)
            if max_snps is not None and len(hap1_cols) >= max_snps:
                break

    if not hap1_cols:
        raise ValueError(f"No usable phased biallelic variants in {path}")

    hap1 = np.column_stack(hap1_cols)
    hap2 = np.column_stack(hap2_cols)
    return PhasedGeno(
        hap1=hap1,
        hap2=hap2,
        ids=sample_ids,
        chrom=np.asarray(chroms, dtype=object),
        positions=np.asarray(positions, dtype=np.int64),
        snp_ids=np.asarray(snp_ids, dtype=object),
        alleles=np.asarray(alleles, dtype=object),
        source=str(path),
        meta={"n_unphased_sites_seen": n_unphased, "format": "vcf"},
    )


def read_oxford_haps(haps_path: str | Path, sample_path: str | Path | None = None) -> PhasedGeno:
    """
    Read SHAPEIT/IMPUTE2 .haps (+ optional .sample).

    .haps columns: chrom id pos a1 a2 then 2*n hap alleles (0/1).
    .sample: ID_1 ID_2 missing [covariates...]; data rows start after header+0 line.
    """
    haps_path = Path(haps_path)
    rows_chrom, rows_id, rows_pos, rows_a1, rows_a2 = [], [], [], [], []
    hap_rows: list[list[int]] = []
    with open(haps_path, "rt") as fh:
        for line in fh:
            parts = line.split()
            if len(parts) < 7:
                continue
            rows_chrom.append(parts[0])
            rows_id.append(parts[1])
            rows_pos.append(int(parts[2]))
            rows_a1.append(parts[3])
            rows_a2.append(parts[4])
            hap_rows.append([int(x) for x in parts[5:]])
    if not hap_rows:
        raise ValueError(f"Empty haps file {haps_path}")
    H = np.asarray(hap_rows, dtype=np.int8)  # (p, 2n)
    if H.shape[1] % 2 != 0:
        raise ValueError("haps allele columns must be even (2 per sample)")
    n = H.shape[1] // 2
    hap1 = H[:, 0::2].T.copy()
    hap2 = H[:, 1::2].T.copy()

    if sample_path is not None:
        sample_path = Path(sample_path)
        ids: list[str] = []
        with open(sample_path, "rt") as fh:
            header = fh.readline()
            _ = fh.readline()  # types line
            for line in fh:
                parts = line.split()
                if not parts:
                    continue
                # ID_1 usually family, ID_2 individual — use ID_2 if present
                ids.append(parts[1] if len(parts) > 1 else parts[0])
        if len(ids) != n:
            raise ValueError(f"sample file has {len(ids)} IDs but haps has {n} samples")
    else:
        ids = [f"S{i}" for i in range(n)]

    return PhasedGeno(
        hap1=hap1,
        hap2=hap2,
        ids=ids,
        chrom=np.asarray(rows_chrom, dtype=object),
        positions=np.asarray(rows_pos, dtype=np.int64),
        snp_ids=np.asarray(rows_id, dtype=object),
        alleles=np.column_stack([rows_a1, rows_a2]),
        source=str(haps_path),
        meta={"format": "oxford_haps"},
    )


def read_plink_pedmap(
    ped_path: str | Path,
    map_path: str | Path,
    *,
    treat_allele_order_as_phase: bool = True,
) -> PhasedGeno:
    """
    Read PLINK text PED/MAP.

    PED allele pairs are NOT guaranteed phased. By default we treat the written
    order as haplotype order (common after Beagle re-export workflows) and set
    meta['phase_assumed']=True. This is NOT a claim of true phase.
    """
    if not treat_allele_order_as_phase:
        raise ValueError("Unphased PED cannot supply haplotypes; export phased VCF/haps instead.")

    map_path = Path(map_path)
    ped_path = Path(ped_path)
    mdf = pd.read_csv(
        map_path,
        sep=r"\s+",
        engine="python",
        header=None,
        names=["chrom", "snp_id", "gm", "pos"],
    )
    p = len(mdf)

    ids: list[str] = []
    hap1_rows: list[list[int]] = []
    hap2_rows: list[list[int]] = []
    # allele coding per SNP: first-seen non-missing allele -> 0, other -> 1
    a0 = [None] * p
    a1 = [None] * p

    with open(ped_path, "rt") as fh:
        for line in fh:
            parts = line.split()
            if len(parts) < 6 + 2 * p:
                raise ValueError(f"PED line has {len(parts)} fields; expected >= {6+2*p}")
            ids.append(parts[1])  # IID
            alleles = parts[6:]
            h1 = []
            h2 = []
            for j in range(p):
                x, y = alleles[2 * j], alleles[2 * j + 1]
                for raw in (x, y):
                    if raw in {"0", ".", "N"}:
                        continue
                    if a0[j] is None:
                        a0[j] = raw
                    elif raw != a0[j] and a1[j] is None:
                        a1[j] = raw
                def code(raw: str, j: int = j) -> int:
                    if raw in {"0", ".", "N"}:
                        return 0  # naive impute to major for MVP
                    if raw == a0[j]:
                        return 0
                    return 1
                h1.append(code(x))
                h2.append(code(y))
            hap1_rows.append(h1)
            hap2_rows.append(h2)

    alleles_arr = np.column_stack([
        [x if x is not None else "A" for x in a0],
        [y if y is not None else "B" for y in a1],
    ])
    return PhasedGeno(
        hap1=np.asarray(hap1_rows, dtype=np.int8),
        hap2=np.asarray(hap2_rows, dtype=np.int8),
        ids=ids,
        chrom=mdf["chrom"].astype(str).to_numpy(),
        positions=mdf["pos"].to_numpy(dtype=np.int64),
        snp_ids=mdf["snp_id"].astype(str).to_numpy(),
        alleles=alleles_arr,
        source=str(ped_path),
        meta={"format": "plink_pedmap", "phase_assumed": True},
    )


def load_geno(path: str | Path, *, companion: str | Path | None = None) -> PhasedGeno:
    """
    Autodetect loader from suffix.

    - .npz / .npz.npz
    - .vcf / .vcf.gz
    - .haps (+ optional companion .sample)
    - .ped (+ companion .map required)
    """
    path = Path(path)
    name = path.name.lower()
    if name.endswith(".npz") or name.endswith(".npz.npz"):
        return load_npz(path)
    if name.endswith(".vcf") or name.endswith(".vcf.gz"):
        return read_phased_vcf(path)
    if name.endswith(".haps"):
        return read_oxford_haps(path, companion)
    if name.endswith(".ped"):
        if companion is None:
            guess = path.with_suffix(".map")
            if not guess.exists():
                raise ValueError("PED input requires --map companion")
            companion = guess
        return read_plink_pedmap(path, companion)
    raise ValueError(f"Unrecognized genotype format: {path}")
