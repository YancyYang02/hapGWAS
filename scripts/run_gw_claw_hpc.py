#!/usr/bin/env python3
"""Genome-wide claw/BW hapGWAS + GHap-like allele baseline on HPC.
Outputs only under HAPGWAS_WORK/gw_results/.
"""
from __future__ import annotations
import json, os, sys, time
from pathlib import Path
import numpy as np
import pandas as pd

WORK = Path(os.environ["HAPGWAS_WORK"])
HYB = Path(os.environ["HAPGWAS_HYB"])
D8 = Path(os.environ["HAPGWAS_D8"])
OUT = WORK / "gw_results"
OUT.mkdir(parents=True, exist_ok=True)

sys.path.insert(0, str(WORK / "hapgwas_src" / "hapgwas"))
from hapgwas.genotype import load_geno
from hapgwas.io import read_phenotype, align_ids
from hapgwas.lmm import fit_null_lmm
from hapgwas.blocks import fixed_windows, summarize_blocks
from hapgwas.joint_test import run_genome, score_snp

PHENO = WORK / "results" / "pheno_continuous.tsv"
if not PHENO.exists():
    # rebuild from claw
    claw = D8 / "gemma_claw_per_bw.pheno"
    # need VCF ids from any chrom
    import gzip
    vcf = sorted(HYB.glob("phased_chip_chr*.vcf.gz"))[0]
    with gzip.open(vcf, "rt") as f:
        for line in f:
            if line.startswith("#CHROM"):
                vids = line.rstrip().split("\t")[9:]
                break
    vals = []
    for line in claw.read_text().splitlines():
        x = line.strip()
        if not x:
            continue
        vals.append(float("nan") if x.upper() in {"NA", "NAN", "."} else float(x))
    assert len(vals) == len(vids)
    df = pd.DataFrame({"IID": vids, "PHENO": vals}).dropna()
    PHENO.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(PHENO, sep="\t", index=False)

pheno = read_phenotype(PHENO)
print("pheno n", len(pheno), flush=True)

vcfs = sorted(HYB.glob("phased_chip_chr*.vcf.gz"))
# numeric sort by chrom token
def chrom_key(p: Path):
    t = p.name.replace("phased_chip_chr", "").replace(".vcf.gz", "")
    try:
        return (0, int(t))
    except ValueError:
        return (1, t)

vcfs = sorted(vcfs, key=chrom_key)
print("n_vcfs", len(vcfs), flush=True)

all_blocks = []
all_alleles = []
all_snps = []
timing = []

for vcf in vcfs:
    chrom = vcf.name.replace("phased_chip_chr", "").replace(".vcf.gz", "")
    t0 = time.time()
    print(f"=== chr{chrom} {vcf.name} ===", flush=True)
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
    # SNP scores
    snp_rows = []
    for j in range(geno.n_snps):
        dos = geno.hap1[:, j].astype(float) + geno.hap2[:, j].astype(float)
        chisq, p = score_snp(null, dos)
        if np.isfinite(p):
            snp_rows.append({
                "chrom": chrom,
                "snp_index": j,
                "snp_id": str(geno.snp_ids[j]),
                "pos": int(geno.positions[j]) if geno.positions is not None else j,
                "snp_p": float(p),
            })
    snp_df = pd.DataFrame(snp_rows)
    # GHap-like: per-block min allele marg_p + Bonferroni within block
    gh = []
    for bid, sub in allele_df.groupby("block_id"):
        ps = pd.to_numeric(sub["marg_p"], errors="coerce").dropna()
        if len(ps) == 0:
            continue
        m = float(ps.min())
        gh.append({
            "chrom": chrom,
            "block_id": bid,
            "ghap_like_min_p": m,
            "ghap_like_bon_p": min(1.0, m * len(ps)),
            "n_alleles": int(len(ps)),
        })
    gh_df = pd.DataFrame(gh)
    # merge onto blocks
    block_df = block_df.merge(gh_df, on=["chrom", "block_id"], how="left")
    # snp bonferroni within block
    snp_bon = []
    for _, r in block_df.iterrows():
        s0, s1 = int(r["snp_start"]), int(r["snp_end"])
        sub = snp_df[(snp_df.snp_index >= s0) & (snp_df.snp_index < s1)]
        if len(sub) == 0:
            snp_bon.append((np.nan, np.nan))
            continue
        m = float(sub["snp_p"].min())
        snp_bon.append((m, min(1.0, m * len(sub))))
    block_df["snp_min_p"] = [a for a, b in snp_bon]
    block_df["snp_bon_p"] = [b for a, b in snp_bon]

    block_df.to_csv(OUT / f"blocks_chr{chrom}.tsv", sep="\t", index=False)
    allele_df.to_csv(OUT / f"alleles_chr{chrom}.tsv", sep="\t", index=False)
    snp_df.to_csv(OUT / f"snps_chr{chrom}.tsv", sep="\t", index=False)
    all_blocks.append(block_df)
    all_alleles.append(allele_df)
    all_snps.append(snp_df)
    dt = time.time() - t0
    timing.append({"chrom": chrom, "n": int(geno.n_samples), "n_snps": int(geno.n_snps), "n_blocks": int(len(block_df)), "sec": round(dt, 2)})
    print(f"chr{chrom} done in {dt:.1f}s  snps={geno.n_snps} blocks={len(block_df)} min_acj={block_df['acj_p'].min():.3g}", flush=True)

blocks = pd.concat(all_blocks, ignore_index=True)
alleles = pd.concat(all_alleles, ignore_index=True)
snps = pd.concat(all_snps, ignore_index=True)
blocks.to_csv(OUT / "blocks_all.tsv", sep="\t", index=False)
alleles.to_csv(OUT / "alleles_all.tsv", sep="\t", index=False)
snps.to_csv(OUT / "snps_all.tsv", sep="\t", index=False)
pd.DataFrame(timing).to_csv(OUT / "timing_by_chrom.tsv", sep="\t", index=False)

# Genome-wide thresholds: Bonferroni on number of blocks / alleles / snps
n_blocks = len(blocks)
n_alleles = len(alleles)
n_snps = len(snps)
alpha = 0.05
thr_block = alpha / max(n_blocks, 1)
thr_allele = alpha / max(n_alleles, 1)
thr_snp = alpha / max(n_snps, 1)

# Hit sets at (1) block-Bonferroni on each method's p, (2) lenient suggestive 1e-4
def hits(df, col, thr):
    return set(df.loc[pd.to_numeric(df[col], errors="coerce") < thr, "block_id"].astype(str))

# Ensure block_id unique across chrom: prefix
blocks["block_uid"] = blocks["chrom"].astype(str) + ":" + blocks["block_id"].astype(str)

def hits_uid(col, thr):
    return set(blocks.loc[pd.to_numeric(blocks[col], errors="coerce") < thr, "block_uid"])

methods = {
    "ACJ": "acj_p",
    "joint": "joint_p",
    "GHap_like_bon": "ghap_like_bon_p",
    "SNP_bon": "snp_bon_p",
    "PCA": "pca_p",
    "ridge": "ridge_p",
}

# Primary: Bonferroni by n_blocks for block-level scores; GHap_like already block-Bonf within then compare at thr_block
# For GHap_like_bon_p and snp_bon_p which are within-block corrected, still apply genome thr_block for discovery table
summary_thr = {}
venn = {}
for name, col in methods.items():
    H = hits_uid(col, thr_block)
    venn[name] = sorted(H)
    summary_thr[name] = {"n_hits_bonf_blocks": len(H), "thr": thr_block, "min_p": float(pd.to_numeric(blocks[col], errors="coerce").min())}

# Suggestive 1e-5
suggestive = {}
for name, col in methods.items():
    suggestive[name] = sorted(hits_uid(col, 1e-5))

# Unique / shared vs GHap-like
acj = set(venn["ACJ"])
gh = set(venn["GHap_like_bon"])
joint = set(venn["joint"])
snp = set(venn["SNP_bon"])

compare = {
    "n_blocks": n_blocks,
    "n_alleles": n_alleles,
    "n_snps": n_snps,
    "thr_block_bonferroni": thr_block,
    "hits_bonferroni": {k: v["n_hits_bonf_blocks"] for k, v in summary_thr.items()},
    "min_p": {k: v["min_p"] for k, v in summary_thr.items()},
    "ACJ_only_vs_GHap": sorted(acj - gh),
    "GHap_only_vs_ACJ": sorted(gh - acj),
    "shared_ACJ_GHap": sorted(acj & gh),
    "ACJ_only_vs_joint": sorted(acj - joint),
    "joint_only_vs_ACJ": sorted(joint - acj),
    "shared_ACJ_joint": sorted(acj & joint),
    "ACJ_not_SNP": sorted(acj - snp),
    "suggestive_1e-5_counts": {k: len(v) for k, v in suggestive.items()},
}

# Top 30 by ACJ with method columns
top = blocks.sort_values("acj_p").head(30)
cols = ["block_uid", "chrom", "block_id", "snp_start", "snp_end", "joint_p", "acj_p", "pca_p", "ridge_p",
        "ghap_like_min_p", "ghap_like_bon_p", "snp_min_p", "snp_bon_p", "acj_df", "acj_merges"]
top[cols].to_csv(OUT / "top30_acj_with_baselines.tsv", sep="\t", index=False)

# Unique hit table
rows = []
for uid in sorted(acj | gh | joint):
    r = blocks.loc[blocks.block_uid == uid].iloc[0]
    rows.append({
        "block_uid": uid,
        "joint_p": float(r["joint_p"]),
        "acj_p": float(r["acj_p"]),
        "ghap_like_bon_p": float(r["ghap_like_bon_p"]) if pd.notna(r["ghap_like_bon_p"]) else np.nan,
        "snp_bon_p": float(r["snp_bon_p"]) if pd.notna(r["snp_bon_p"]) else np.nan,
        "in_ACJ_bonf": uid in acj,
        "in_GHap_bonf": uid in gh,
        "in_joint_bonf": uid in joint,
        "in_SNP_bonf": uid in snp,
        "class": (
            "shared_ACJ_GHap" if uid in acj & gh else
            "ACJ_only" if uid in acj - gh else
            "GHap_only" if uid in gh - acj else
            "joint_only" if uid in joint - acj - gh else
            "other"
        ),
    })
pd.DataFrame(rows).to_csv(OUT / "hits_unique_shared.tsv", sep="\t", index=False)
(OUT / "gw_compare_summary.json").write_text(json.dumps(compare, indent=2))
pd.DataFrame(timing).to_csv(OUT / "timing_by_chrom.tsv", sep="\t", index=False)

# Markdown summary (small)
lines = []
lines.append("# Genome-wide claw/BW hapGWAS (HPC)\n\n")
lines.append(f"- n={len(pheno)}, chromosomes={len(vcfs)}, blocks={n_blocks}, alleles={n_alleles}, snps={n_snps}\n")
lines.append(f"- Bonferroni thr (0.05/n_blocks) = {thr_block:.3g}\n")
lines.append(f"- Walltime total sec = {sum(t['sec'] for t in timing):.1f}\n\n")
lines.append("## Hits @ block-Bonferroni\n\n")
lines.append("| method | n_hits | min_p |\n|---|---:|---:|\n")
for name, col in methods.items():
    lines.append(f"| {name} | {summary_thr[name]['n_hits_bonf_blocks']} | {summary_thr[name]['min_p']:.3g} |\n")
lines.append("\n## ACJ vs GHap-like (allele-Bonf within block, then genome Bonf)\n\n")
lines.append(f"- shared: **{len(acj & gh)}**\n")
lines.append(f"- ACJ only: **{len(acj - gh)}** → {', '.join(sorted(acj - gh)[:20]) or '—'}\n")
lines.append(f"- GHap only: **{len(gh - acj)}** → {', '.join(sorted(gh - acj)[:20]) or '—'}\n")
lines.append(f"- ACJ hits not in SNP-Bonf: **{len(acj - snp)}**\n")
lines.append("\nProxy note: GHap-like = hapGWAS allele marginal min-p with within-block Bonferroni (same dosages/null as ACJ). Full R GHap not required for head-to-head on identical blocks.\n")
(OUT / "GW_SUMMARY_zh.md").write_text("".join(lines), encoding="utf-8")

# optional figure
try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(5, 4))
    cats = ["shared", "ACJ only", "GHap only"]
    vals = [len(acj & gh), len(acj - gh), len(gh - acj)]
    ax.bar(cats, vals, color=["#54A24B", "#E45756", "#4C78A8"])
    ax.set_ylabel("n blocks (Bonferroni)")
    ax.set_title("ACJ vs GHap-like hits")
    fig.tight_layout()
    fig.savefig(OUT / "fig_acj_vs_ghap_hits.png", dpi=140)
    plt.close()
    # manhattan ACJ
    fig, ax = plt.subplots(figsize=(10, 3.5))
    # map chrom order
    chroms = [t["chrom"] for t in timing]
    x = 0
    xticks, xlabels = [], []
    for ch in chroms:
        sub = blocks[blocks.chrom.astype(str) == str(ch)]
        if len(sub) == 0:
            continue
        xpos = np.arange(len(sub)) + x
        ax.scatter(xpos, -np.log10(np.clip(sub["acj_p"], 1e-300, 1)), s=4, c="#54A24B", label="ACJ" if ch == chroms[0] else None)
        ax.scatter(xpos, -np.log10(np.clip(sub["ghap_like_bon_p"], 1e-300, 1)), s=3, c="#4C78A8", alpha=0.5, label="GHap-like" if ch == chroms[0] else None)
        xticks.append(x + len(sub) / 2)
        xlabels.append(str(ch))
        x += len(sub) + 10
    ax.axhline(-np.log10(thr_block), color="red", ls="--", lw=0.8)
    ax.set_xticks(xticks)
    ax.set_xticklabels(xlabels, fontsize=7)
    ax.set_ylabel("-log10 p")
    ax.set_title("Genome-wide claw/BW: ACJ vs GHap-like (block)")
    ax.legend(loc="upper right", fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT / "fig_gw_manhattan_acj_ghap.png", dpi=140)
    plt.close()
except Exception as e:
    print("plot skip", e, flush=True)

print("DONE", OUT, flush=True)
print(json.dumps(compare, indent=2)[:2000], flush=True)
