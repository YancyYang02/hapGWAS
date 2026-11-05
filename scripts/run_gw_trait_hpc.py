#!/usr/bin/env python3
"""Genome-wide hapGWAS + GHap-like baseline for one continuous trait on HPC.
Env:
  HAPGWAS_WORK, HAPGWAS_HYB, HAPGWAS_D8
  HAPGWAS_TRAIT_NAME (e.g. wing)
  HAPGWAS_TRAIT_COL  (1-based column in gemma_final.pheno; default 4 = Wing)
  HAPGWAS_OUT_SUBDIR (default gw_results_<trait>)
"""
from __future__ import annotations
import json, os, sys, time, gzip
from pathlib import Path
import numpy as np
import pandas as pd

WORK = Path(os.environ["HAPGWAS_WORK"])
HYB = Path(os.environ["HAPGWAS_HYB"])
D8 = Path(os.environ["HAPGWAS_D8"])
TRAIT = os.environ.get("HAPGWAS_TRAIT_NAME", "wing")
COL = int(os.environ.get("HAPGWAS_TRAIT_COL", "4"))  # 1-based
OUT = WORK / os.environ.get("HAPGWAS_OUT_SUBDIR", f"gw_results_{TRAIT}")
OUT.mkdir(parents=True, exist_ok=True)

sys.path.insert(0, str(WORK / "hapgwas_src" / "hapgwas"))
from hapgwas.genotype import load_geno
from hapgwas.io import read_phenotype, align_ids
from hapgwas.lmm import fit_null_lmm
from hapgwas.blocks import fixed_windows
from hapgwas.joint_test import run_genome, score_snp

# Build phenotype from gemma_final column
PHENO = OUT / f"pheno_{TRAIT}.tsv"
vcf0 = sorted(HYB.glob("phased_chip_chr*.vcf.gz"))[0]
with gzip.open(vcf0, "rt") as f:
    for line in f:
        if line.startswith("#CHROM"):
            vids = line.rstrip().split("\t")[9:]
            break
final = D8 / "gemma_final.pheno"
vals = []
for line in final.read_text().splitlines():
    if not line.strip():
        continue
    parts = line.split()
    x = parts[COL - 1]
    vals.append(float("nan") if x.upper() in {"NA", "NAN", "."} else float(x))
assert len(vals) == len(vids), (len(vals), len(vids))
df = pd.DataFrame({"IID": vids, "PHENO": vals}).dropna()
df.to_csv(PHENO, sep="\t", index=False)
print(f"trait={TRAIT} col={COL} n={len(df)} mean={df.PHENO.mean():.4g} sd={df.PHENO.std():.4g}", flush=True)

pheno = read_phenotype(PHENO)

def chrom_key(p: Path):
    t = p.name.replace("phased_chip_chr", "").replace(".vcf.gz", "")
    try:
        return (0, int(t))
    except ValueError:
        return (1, t)

vcfs = sorted(HYB.glob("phased_chip_chr*.vcf.gz"), key=chrom_key)
print("n_vcfs", len(vcfs), flush=True)

all_blocks, all_alleles, all_snps, timing = [], [], [], []

for vcf in vcfs:
    chrom = vcf.name.replace("phased_chip_chr", "").replace(".vcf.gz", "")
    t0 = time.time()
    print(f"=== chr{chrom} ===", flush=True)
    geno = load_geno(vcf)
    y, keep, X, kept_ids = align_ids(pheno, geno.ids, None)
    geno = geno.subset_samples(keep)
    null = fit_null_lmm(y, X=X, K=None, method="ols")
    blocks = fixed_windows(geno.n_snps, window=4, step=2, chrom=str(chrom), positions=geno.positions)
    block_df, allele_df = run_genome(
        null, blocks, geno.hap1, geno.hap2,
        min_freq=0.02, min_count=5, collapse_mode="freq",
        snp_ids=geno.snp_ids,
    )
    block_df["chrom"] = chrom
    allele_df["chrom"] = chrom
    snp_rows = []
    for j in range(geno.n_snps):
        dos = geno.hap1[:, j].astype(float) + geno.hap2[:, j].astype(float)
        chisq, p = score_snp(null, dos)
        if np.isfinite(p):
            snp_rows.append({
                "chrom": chrom, "snp_index": j, "snp_id": str(geno.snp_ids[j]),
                "pos": int(geno.positions[j]) if geno.positions is not None else j,
                "snp_p": float(p),
            })
    snp_df = pd.DataFrame(snp_rows)
    gh = []
    for bid, sub in allele_df.groupby("block_id"):
        ps = pd.to_numeric(sub["marg_p"], errors="coerce").dropna()
        if len(ps) == 0:
            continue
        m = float(ps.min())
        gh.append({
            "chrom": chrom, "block_id": bid,
            "ghap_like_min_p": m,
            "ghap_like_bon_p": min(1.0, m * len(ps)),
            "n_alleles": int(len(ps)),
        })
    gh_df = pd.DataFrame(gh)
    block_df = block_df.merge(gh_df, on=["chrom", "block_id"], how="left")
    snp_bon = []
    for _, r in block_df.iterrows():
        s0, s1 = int(r["snp_start"]), int(r["snp_end"])
        sub = snp_df[(snp_df.snp_index >= s0) & (snp_df.snp_index < s1)]
        if len(sub) == 0:
            snp_bon.append((np.nan, np.nan)); continue
        m = float(sub["snp_p"].min())
        snp_bon.append((m, min(1.0, m * len(sub))))
    block_df["snp_min_p"] = [a for a, b in snp_bon]
    block_df["snp_bon_p"] = [b for a, b in snp_bon]
    block_df.to_csv(OUT / f"blocks_chr{chrom}.tsv", sep="\t", index=False)
    allele_df.to_csv(OUT / f"alleles_chr{chrom}.tsv", sep="\t", index=False)
    snp_df.to_csv(OUT / f"snps_chr{chrom}.tsv", sep="\t", index=False)
    all_blocks.append(block_df); all_alleles.append(allele_df); all_snps.append(snp_df)
    dt = time.time() - t0
    timing.append({"chrom": chrom, "n": int(geno.n_samples), "n_snps": int(geno.n_snps),
                   "n_blocks": int(len(block_df)), "sec": round(dt, 2)})
    print(f"chr{chrom} {dt:.1f}s snps={geno.n_snps} blocks={len(block_df)} min_acj={block_df['acj_p'].min():.3g}", flush=True)

blocks = pd.concat(all_blocks, ignore_index=True)
alleles = pd.concat(all_alleles, ignore_index=True)
snps = pd.concat(all_snps, ignore_index=True)
blocks.to_csv(OUT / "blocks_all.tsv", sep="\t", index=False)
alleles.to_csv(OUT / "alleles_all.tsv", sep="\t", index=False)
snps.to_csv(OUT / "snps_all.tsv", sep="\t", index=False)
pd.DataFrame(timing).to_csv(OUT / "timing_by_chrom.tsv", sep="\t", index=False)

n_blocks, n_alleles, n_snps = len(blocks), len(alleles), len(snps)
thr_block = 0.05 / max(n_blocks, 1)
blocks["block_uid"] = blocks["chrom"].astype(str) + ":" + blocks["block_id"].astype(str)

methods = {
    "ACJ": "acj_p", "joint": "joint_p", "GHap_like_bon": "ghap_like_bon_p",
    "SNP_bon": "snp_bon_p", "PCA": "pca_p", "ridge": "ridge_p",
}

def hits_uid(col, thr):
    return set(blocks.loc[pd.to_numeric(blocks[col], errors="coerce") < thr, "block_uid"])

venn = {name: hits_uid(col, thr_block) for name, col in methods.items()}
summary_thr = {
    name: {
        "n_hits_bonf_blocks": len(venn[name]),
        "thr": thr_block,
        "min_p": float(pd.to_numeric(blocks[col], errors="coerce").min()),
    }
    for name, col in methods.items()
}
acj, gh, joint, snp = set(venn["ACJ"]), set(venn["GHap_like_bon"]), set(venn["joint"]), set(venn["SNP_bon"])
compare = {
    "trait": TRAIT, "pheno_col_1based": COL, "n": int(len(df)),
    "n_blocks": n_blocks, "n_alleles": n_alleles, "n_snps": n_snps,
    "thr_block_bonferroni": thr_block,
    "hits_bonferroni": {k: v["n_hits_bonf_blocks"] for k, v in summary_thr.items()},
    "min_p": {k: v["min_p"] for k, v in summary_thr.items()},
    "ACJ_only_vs_GHap": sorted(acj - gh),
    "GHap_only_vs_ACJ": sorted(gh - acj),
    "shared_ACJ_GHap": sorted(acj & gh),
    "ACJ_not_SNP": sorted(acj - snp),
    "suggestive_1e-5_counts": {k: len(hits_uid(col, 1e-5)) for k, col in methods.items()},
    "walltime_sec": float(sum(t["sec"] for t in timing)),
}
(OUT / "gw_compare_summary.json").write_text(json.dumps(compare, indent=2))

cols = ["block_uid", "chrom", "block_id", "snp_start", "snp_end", "joint_p", "acj_p", "pca_p", "ridge_p",
        "ghap_like_min_p", "ghap_like_bon_p", "snp_min_p", "snp_bon_p", "acj_df", "acj_merges"]
blocks.sort_values("acj_p").head(30)[cols].to_csv(OUT / "top30_acj_with_baselines.tsv", sep="\t", index=False)

rows = []
for uid in sorted(acj | gh | joint):
    r = blocks.loc[blocks.block_uid == uid].iloc[0]
    rows.append({
        "block_uid": uid,
        "joint_p": float(r["joint_p"]), "acj_p": float(r["acj_p"]),
        "ghap_like_bon_p": float(r["ghap_like_bon_p"]) if pd.notna(r["ghap_like_bon_p"]) else np.nan,
        "snp_bon_p": float(r["snp_bon_p"]) if pd.notna(r["snp_bon_p"]) else np.nan,
        "in_ACJ_bonf": uid in acj, "in_GHap_bonf": uid in gh,
        "in_joint_bonf": uid in joint, "in_SNP_bonf": uid in snp,
        "class": (
            "shared_ACJ_GHap" if uid in acj & gh else
            "ACJ_only" if uid in acj - gh else
            "GHap_only" if uid in gh - acj else
            "joint_only" if uid in joint - acj - gh else "other"
        ),
    })
pd.DataFrame(rows).to_csv(OUT / "hits_unique_shared.tsv", sep="\t", index=False)

lines = [
    f"# Genome-wide {TRAIT} hapGWAS (HPC)\n\n",
    f"- trait=`{TRAIT}` (gemma_final col {COL}), n={len(df)}, chroms={len(vcfs)}, "
    f"blocks={n_blocks}, alleles={n_alleles}, snps={n_snps}\n",
    f"- Bonferroni thr = {thr_block:.3g}; walltime = {sum(t['sec'] for t in timing):.1f}s\n\n",
    "## Hits @ block-Bonferroni\n\n| method | n_hits | min_p |\n|---|---:|---:|\n",
]
for name, col in methods.items():
    lines.append(f"| {name} | {summary_thr[name]['n_hits_bonf_blocks']} | {summary_thr[name]['min_p']:.3g} |\n")
lines += [
    "\n## ACJ vs GHap-like\n\n",
    f"- shared: **{len(acj & gh)}**\n",
    f"- ACJ only: **{len(acj - gh)}** → {', '.join(sorted(acj - gh)[:15]) or '—'}\n",
    f"- GHap only: **{len(gh - acj)}** → {', '.join(sorted(gh - acj)[:15]) or '—'}\n",
    f"- ACJ not in SNP-Bonf: **{len(acj - snp)}**\n",
]
(OUT / "GW_SUMMARY_zh.md").write_text("".join(lines), encoding="utf-8")

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(5, 4))
    ax.bar(["shared", "ACJ only", "GHap only"],
           [len(acj & gh), len(acj - gh), len(gh - acj)],
           color=["#54A24B", "#E45756", "#4C78A8"])
    ax.set_ylabel("n blocks (Bonferroni)")
    ax.set_title(f"{TRAIT}: ACJ vs GHap-like")
    fig.tight_layout(); fig.savefig(OUT / "fig_acj_vs_ghap_hits.png", dpi=140); plt.close()
    fig, ax = plt.subplots(figsize=(10, 3.5))
    chroms = [t["chrom"] for t in timing]
    x = 0; xticks, xlabels = [], []
    for ch in chroms:
        sub = blocks[blocks.chrom.astype(str) == str(ch)]
        if len(sub) == 0: continue
        xpos = np.arange(len(sub)) + x
        ax.scatter(xpos, -np.log10(np.clip(sub["acj_p"], 1e-300, 1)), s=4, c="#54A24B",
                   label="ACJ" if ch == chroms[0] else None)
        ax.scatter(xpos, -np.log10(np.clip(sub["ghap_like_bon_p"], 1e-300, 1)), s=3, c="#4C78A8", alpha=0.5,
                   label="GHap-like" if ch == chroms[0] else None)
        xticks.append(x + len(sub) / 2); xlabels.append(str(ch)); x += len(sub) + 10
    ax.axhline(-np.log10(thr_block), color="red", ls="--", lw=0.8)
    ax.set_xticks(xticks); ax.set_xticklabels(xlabels, fontsize=7)
    ax.set_ylabel("-log10 p"); ax.set_title(f"Genome-wide {TRAIT}: ACJ vs GHap-like")
    ax.legend(loc="upper right", fontsize=8)
    fig.tight_layout(); fig.savefig(OUT / "fig_gw_manhattan_acj_ghap.png", dpi=140); plt.close()
except Exception as e:
    print("plot skip", e, flush=True)

print("DONE", OUT, flush=True)
print(json.dumps(compare, indent=2)[:2500], flush=True)
