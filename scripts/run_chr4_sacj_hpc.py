#!/usr/bin/env python3
"""chr4 claw/BW: SACJ vs ACJ/joint/GHap-like/SNP on HPC."""
from __future__ import annotations
import json, os, sys, time
from pathlib import Path
import numpy as np
import pandas as pd

WORK = Path(os.environ["HAPGWAS_WORK"])
HYB = Path(os.environ["HAPGWAS_HYB"])
OUT = WORK / "sacj_chr4"
OUT.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(WORK / "hapgwas_src" / "hapgwas"))

from hapgwas.genotype import load_geno
from hapgwas.io import read_phenotype, align_ids
from hapgwas.lmm import fit_null_lmm
from hapgwas.blocks import fixed_windows
from hapgwas.joint_test import run_genome, score_snp

PHENO = WORK / "results" / "pheno_continuous.tsv"
vcf = HYB / "phased_chip_chr4.vcf.gz"
assert vcf.exists(), vcf
assert PHENO.exists(), PHENO

pheno = read_phenotype(PHENO)
t0 = time.time()
geno = load_geno(vcf)
y, keep, X, _ = align_ids(pheno, geno.ids, None)
geno = geno.subset_samples(keep)
null = fit_null_lmm(y, X=X, K=None, method="ols")
blocks = fixed_windows(geno.n_snps, window=4, step=2, chrom="4", positions=geno.positions)
block_df, allele_df = run_genome(
    null, blocks, geno.hap1, geno.hap2,
    min_freq=0.02, min_count=5, collapse_mode="freq", snp_ids=geno.snp_ids,
)
# GHap-like + SNP bon within block
gh = []
for bid, sub in allele_df.groupby("block_id"):
    ps = pd.to_numeric(sub["marg_p"], errors="coerce").dropna()
    if len(ps)==0: continue
    m=float(ps.min()); gh.append({"block_id":bid,"ghap_like_bon_p":min(1.0,m*len(ps)),"ghap_like_min_p":m})
gh_df=pd.DataFrame(gh)
block_df=block_df.merge(gh_df,on="block_id",how="left")
snp_rows=[]
for j in range(geno.n_snps):
    dos=geno.hap1[:,j].astype(float)+geno.hap2[:,j].astype(float)
    _,p=score_snp(null,dos)
    if np.isfinite(p):
        snp_rows.append({"snp_index":j,"snp_p":float(p)})
snp_df=pd.DataFrame(snp_rows)
sb=[]
for _,r in block_df.iterrows():
    s0,s1=int(r.snp_start),int(r.snp_end)
    sub=snp_df[(snp_df.snp_index>=s0)&(snp_df.snp_index<s1)]
    if len(sub)==0: sb.append(np.nan); continue
    m=float(sub.snp_p.min()); sb.append(min(1.0,m*len(sub)))
block_df["snp_bon_p"]=sb
dt=time.time()-t0

n_blocks=len(block_df)
thr=0.05/max(n_blocks,1)
methods={"joint":"joint_p","acj":"acj_p","sacj":"sacj_p","ghap_like":"ghap_like_bon_p","snp":"snp_bon_p","pca":"pca_p","ridge":"ridge_p"}
hits={name:int((pd.to_numeric(block_df[col],errors="coerce")<thr).sum()) for name,col in methods.items()}
mins={name:float(pd.to_numeric(block_df[col],errors="coerce").min()) for name,col in methods.items()}

# unique/shared SACJ vs ACJ / GHap
def hitset(col):
    return set(block_df.loc[pd.to_numeric(block_df[col],errors="coerce")<thr,"block_id"].astype(str))
S,A,G=hitset("sacj_p"),hitset("acj_p"),hitset("ghap_like_bon_p")
compare={
    "n":int(len(y)),"n_snps":int(geno.n_snps),"n_blocks":n_blocks,"thr":thr,"sec":round(dt,2),
    "hits":hits,"min_p":mins,
    "shared_SACJ_ACJ":sorted(S&A),"SACJ_only_vs_ACJ":sorted(S-A),"ACJ_only_vs_SACJ":sorted(A-S),
    "shared_SACJ_GHap":sorted(S&G),"SACJ_only_vs_GHap":sorted(S-G),"GHap_only_vs_SACJ":sorted(G-S),
}
block_df.to_csv(OUT/"blocks_chr4_sacj.tsv",sep="\t",index=False)
top=block_df.sort_values("sacj_p").head(25)
cols=[c for c in ["block_id","joint_p","acj_p","sacj_p","sacj_k","sacj_ngroups","ghap_like_bon_p","snp_bon_p","pca_p","ridge_p"] if c in top.columns]
top[cols].to_csv(OUT/"top25_sacj.tsv",sep="\t",index=False)
(OUT/"chr4_sacj_compare.json").write_text(json.dumps(compare,indent=2))
md=[f"# chr4 claw/BW SACJ compare\n\nn={len(y)} snps={geno.n_snps} blocks={n_blocks} thr={thr:.3g} sec={dt:.1f}\n\n"]
md.append("| method | hits | min_p |\n|---|---:|---:|\n")
for name in methods:
    md.append(f"| {name} | {hits[name]} | {mins[name]:.3g} |\n")
md.append(f"\nSACJ∩ACJ={len(S&A)} SACJ\\\\ACJ={len(S-A)} ACJ\\\\SACJ={len(A-S)}\n")
md.append(f"SACJ∩GHap={len(S&G)} SACJ\\\\GHap={len(S-G)} GHap\\\\SACJ={len(G-S)}\n")
(OUT/"CHR4_SACJ_zh.md").write_text("".join(md),encoding="utf-8")
print("DONE", OUT)
print(json.dumps(compare,indent=2)[:2000])
