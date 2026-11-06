#!/usr/bin/env python3
"""On HPC (broiler): pure physical/production multi-trait panel (NO metabolites).

Historical v0.4 script. `hapgwas run-multi` and `--methods` are not in the 0.7.1 CLI.
"""
from __future__ import annotations
import os, sys, time, shutil
from pathlib import Path
import numpy as np
import pandas as pd

WORK = Path(os.environ.get("HAPGWAS_WORK", "/fs2/home/broiler/YancongYang/hapgwas_broiler_multi"))
GPE = Path("/fs2/home/clwen/YancongYang/Genetic_Parameter_Estimation")
ZHOU = Path("/fs2/home/broiler/zhouyan/GWAS")
VCF_PHASED = ZHOU / "imputed" / "2.Filter_SNPs.vcf"
OUT = WORK / "run_chr4"
OUT.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(WORK / "hapgwas_src"))
print("WORK", WORK, flush=True)

# Exact ALL_TRAITS order from GPE/run_all_traits.sh (61 pure production/carcass traits)
ALL_TRAITS = [
    "BW3", "BW6", "BW3_6", "BW", "BW3_zq", "BW6_zq", "FI", "AFI",
    "RFI1", "RFI2", "RFI3", "RFI4", "FRI",
    "PMA", "PMAR", "PMI", "PMIR", "PM", "PMR",
    "LW", "LWR", "HW", "HWR", "AFW", "AFWR", "MW", "MWR",
    "GAS", "GASR", "TM2", "TMR2", "RV", "RVR", "TV", "TVR", "AHI",
    "TM1", "TMR1", "CW", "CWR", "DW", "DWR", "SW", "SWR", "TW", "TWR",
    "WW", "WWR", "BFCf", "LFCf", "MCCM", "WCLM", "TNM", "TM_CSA", "PM_CSA",
    "EBW3", "EBW6", "GBW3", "GBW6", "GRFI", "ERFI",
]
assert len(ALL_TRAITS) == 61, len(ALL_TRAITS)

TRAIT_ZH = {
    "BW3": "21日龄体重", "BW6": "42日龄体重", "BW3_6": "3-6周增重", "BW": "宰前体重",
    "BW3_zq": "3周-宰前增重", "BW6_zq": "6周-宰前增重", "FI": "采食量", "AFI": "平均日采食量",
    "RFI1": "剩余采食量1", "RFI2": "剩余采食量2", "RFI3": "剩余采食量3", "RFI4": "剩余采食量4",
    "FRI": "饲料转化率FCR", "PMA": "胸大肌重", "PMAR": "胸大肌率", "PMI": "胸小肌重",
    "PMIR": "胸小肌率", "PM": "胸肌重", "PMR": "胸肌率", "LW": "肝脏重", "LWR": "肝脏比",
    "HW": "心脏重", "HWR": "心脏比", "AFW": "腹脂重", "AFWR": "腹脂率", "MW": "肌胃重",
    "MWR": "肌胃比", "GAS": "腓肠肌重", "GASR": "腓肠肌比", "TM2": "腿重(去毛皮)", "TMR2": "腿重比2",
    "RV": "右心室重", "RVR": "右心室比", "TV": "心室总重", "TVR": "心室比", "AHI": "腹水心脏指数",
    "TM1": "腿重", "TMR1": "腿重比1", "CW": "爪重", "CWR": "爪比", "DW": "翅根重", "DWR": "翅根比",
    "SW": "翅中重", "SWR": "翅中比", "TW": "翅尖重", "TWR": "翅尖比", "WW": "翅重", "WWR": "翅重比",
    "BFCf": "胸肌脂肪(鲜)", "LFCf": "腿肌脂肪(鲜)", "MCCM": "胸肌含水量", "WCLM": "腿肌含水量",
    "TNM": "腓肠肌纤维总数", "TM_CSA": "腓肠肌纤维面积", "PM_CSA": "胸大肌纤维面积",
    "EBW3": "EBW3", "EBW6": "EBW6", "GBW3": "GBW3", "GBW6": "GBW6", "GRFI": "GRFI", "ERFI": "ERFI",
}

# Prefer gcta_pheno_final (FID IID + 61 traits); fallback phe_gemma_all (61 cols, fam order)
pheno_src = Path(os.environ.get("HAPGWAS_PHENO", str(WORK / "gcta_pheno_final.txt")))
if not pheno_src.exists():
    pheno_src = GPE / "gcta_pheno_final.txt"
if not pheno_src.exists():
    pheno_src = WORK / "gcta_pheno_final.txt"
# Prefer zhouyan/step1 fam (a1/a10 IDs matching phased VCF + gcta_pheno_final).
# gemma_463_gwas_updated.fam uses pedigree IIDs (1401...) — do NOT use for VCF merge.
fam_path = ZHOU / "1.1Filter_SNPs.fam"
if not fam_path.exists():
    fam_path = GPE / "gemma_463_gwas_step1.fam"

fam = pd.read_csv(fam_path, sep=r"\s+", header=None, names=["FID", "IID", "F", "M", "SEX", "PHENO"])
fam["IID"] = fam["IID"].astype(str)
print("fam", fam_path, "n", len(fam), "sample", fam["IID"].head(3).tolist(), flush=True)

raw = pd.read_csv(pheno_src, sep=r"\s+", header=None)
print("pheno_src", pheno_src, "shape", raw.shape, flush=True)
if raw.shape[1] == 63:
    raw.columns = ["FID", "IID"] + ALL_TRAITS
elif raw.shape[1] == 61:
    # phe_gemma_all: row order matches step1/zhouyan fam
    raw.columns = ALL_TRAITS
    raw.insert(0, "IID", fam["IID"].values[: len(raw)])
else:
    raise SystemExit(f"unexpected ncols {raw.shape[1]}")

raw["IID"] = raw["IID"].astype(str)
# Keep pheno birds; left-join onto fam so geno-missing birds stay but traits NA
wide = fam[["IID"]].merge(raw[["IID"] + ALL_TRAITS], on="IID", how="inner")
if len(wide) < 100:
    # fallback: use pheno IIDs alone (VCF overlap checked later)
    print("WARN inner-merge small", len(wide), "fallback to pheno IIDs", flush=True)
    wide = raw[["IID"] + ALL_TRAITS].copy()
# coerce numeric
for t in ALL_TRAITS:
    wide[t] = pd.to_numeric(wide[t], errors="coerce")
n_nonmiss = wide[ALL_TRAITS].notna().any(axis=1).sum()
print("wide", wide.shape, "birds_with_any_pheno", int(n_nonmiss), flush=True)
print("traits", ALL_TRAITS, flush=True)

pheno_path = WORK / "pheno_pure61.tsv"
wide.to_csv(pheno_path, sep="\t", index=False)
(WORK / "trait_list.txt").write_text("\n".join(ALL_TRAITS), encoding="utf-8")
(WORK / "trait_zh.tsv").write_text(
    "trait\tzh\n" + "\n".join(f"{t}\t{TRAIT_ZH.get(t,t)}" for t in ALL_TRAITS),
    encoding="utf-8",
)
print("wrote", pheno_path, flush=True)

n_var = int(os.environ.get("HAPGWAS_NVAR", "5000"))
vcf_out = WORK / f"phased_first{n_var}.vcf.gz"
src_existing = WORK / "phased_first20000.vcf.gz"
if src_existing.exists() and n_var <= 20000:
    import gzip as gz
    print("subset existing vcf to", n_var, flush=True)
    n_keep = 0
    with gz.open(src_existing, "rt") as fin, gz.open(vcf_out, "wt") as fout:
        for line in fin:
            if line.startswith("#"):
                fout.write(line); continue
            fout.write(line); n_keep += 1
            if n_keep >= n_var:
                break
    print("kept", n_keep, "size", vcf_out.stat().st_size, flush=True)
elif (not vcf_out.exists()) or vcf_out.stat().st_size < 1000:
    import gzip as gz
    print("extracting first", n_var, flush=True)
    t0 = time.time(); n_keep = 0
    with open(VCF_PHASED, "rt", errors="replace") as fin, gz.open(vcf_out, "wt") as fout:
        for line in fin:
            if line.startswith("#"):
                fout.write(line); continue
            fout.write(line); n_keep += 1
            if n_keep >= n_var:
                break
    print("kept", n_keep, "sec", round(time.time() - t0, 1), flush=True)
else:
    print("reuse", vcf_out, "size", vcf_out.stat().st_size, flush=True)

# Verify IID overlap geno vs pheno
import gzip as gz
with gz.open(vcf_out, "rt") as fh:
    for line in fh:
        if line.startswith("#CHROM"):
            geno_ids = line.rstrip("\n").split("\t")[9:]
            break
overlap = set(geno_ids) & set(wide["IID"].astype(str))
print("geno_n", len(geno_ids), "pheno_n", len(wide), "overlap", len(overlap), flush=True)
if len(overlap) < 100:
    raise SystemExit(f"too few overlapping IDs: {len(overlap)}")

from hapgwas.cli import main as cli_main
out_dir = OUT / "multi"
if out_dir.exists():
    shutil.rmtree(out_dir)
argv = [
    "run-multi", "--geno", str(vcf_out), "--pheno", str(pheno_path), "--ols",
    "--methods", "all", "--block-mode", "window", "--window", "4", "--step", "2",
    "--min-freq", "0.05", "--min-count", "5", "--out", str(out_dir),
]
print("CLI", argv, flush=True)
t0 = time.time()
rc = cli_main(argv)
print("cli_rc", rc, "sec", round(time.time() - t0, 1), flush=True)
hits = pd.read_csv(out_dir / "multi_trait_hits.tsv", sep="\t")

sacj_win = int(hits["sacj_beats_snp"].sum()) if "sacj_beats_snp" in hits else 0
acj_win = int(hits["acj_beats_snp"].sum()) if "acj_beats_snp" in hits else 0
n_traits = len(hits)

md = []
md.append("# Broiler multi-trait hapGWAS（纯生产/屠体表型）\n\n")
md.append("## 更正说明\n")
md.append("- **上一版错误**：用了 9 个生产性状 + top-40 高 h² **代谢物**（`phenotype2.txt` / `Id#####` 等），不符合课题设定。\n")
md.append("- **本版更正**：仅使用纯物理/生产/屠体表型面板（体重、增重、采食与饲料效率、胸腿翅爪器官重量与比率、肌纤维等），**不含任何代谢物/脂质**。\n\n")
md.append("## 一句话汇总\n")
md.append(f"- 队列：Genetic_Parameter_Estimation / zhouyan broiler **n={len(overlap)}**（fam {len(fam)}），表型 **{n_traits}** 个纯生产性状。\n")
md.append(f"- 表型来源：`{pheno_src}`（列序 = `run_all_traits.sh` 的 `ALL_TRAITS`）。\n")
md.append(f"- 基因型：phased VCF 前 {n_var} 位点子集（`{VCF_PHASED.name}`），与 gemma_463 同队列 IID。\n")
md.append("- 工具：`hapgwas run-multi --methods all` → joint / ACJ / SACJ / allele_bonf / snp_bonf / pca / ridge。\n")
md.append(f"- **SACJ 命中数 > SNP-bonf：{sacj_win}/{n_traits} 性状**；ACJ 同指标：**{acj_win}/{n_traits}**。\n\n")
md.append("## 数据与方法\n")
md.append(f"- cohort geno fam: `{fam_path}`\n")
md.append(f"- pheno: `{pheno_src}` → `{pheno_path.name}`，**n_traits={n_traits}**（纯生产，无代谢物）\n")
md.append(f"- geno subset: `{vcf_out.name}` (first {n_var} phased sites)\n")
md.append("- methods: joint / ACJ / SACJ / allele_bonf (GHap-like) / snp_bonf / pca / ridge\n")
md.append("- CLI: `hapgwas run-multi --methods all --block-mode window --window 4 --step 2`\n\n")
md.append("## 61 个性状列表（ALL_TRAITS）\n\n")
md.append("| # | trait | 中文/含义 |\n|---:|---|---|\n")
for i, t in enumerate(ALL_TRAITS, 1):
    md.append(f"| {i} | `{t}` | {TRAIT_ZH.get(t, t)} |\n")
md.append("\n## 各性状命中汇总\n\n")
md.append("| trait | n | joint | ACJ | SACJ | alleleBon | snpBon | SACJ>SNP | ACJ>SNP |\n")
md.append("|---|---:|---:|---:|---:|---:|---:|---:|---:|\n")
for _, r in hits.iterrows():
    md.append(
        f"| {r['trait']} | {int(r['n'])} | {int(r.get('hits_joint', 0))} | {int(r.get('hits_acj', 0))} | "
        f"{int(r.get('hits_sacj', 0))} | {int(r.get('hits_allele_bonf', 0))} | {int(r.get('hits_snp_bonf', 0))} | "
        f"{int(r.get('sacj_beats_snp', 0))} | {int(r.get('acj_beats_snp', 0))} |\n"
    )
md.append(f"\nSACJ hits>SNP: **{sacj_win}**/{n_traits}; ACJ hits>SNP: **{acj_win}**/{n_traits}\n")
(OUT / "BROILER_MULTI_TRAIT_zh.md").write_text("".join(md), encoding="utf-8")
hits.to_csv(OUT / "multi_trait_hits.tsv", sep="\t", index=False)
print("DONE", OUT, "sacj_win", sacj_win, "acj_win", acj_win, flush=True)
