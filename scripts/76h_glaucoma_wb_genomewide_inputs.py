#!/usr/bin/env python3
# ============================================================
# Step 76h-a (29 Oct 2026; 10 Nov 2026): per-chromosome FUSION inputs for the genome-wide whole-blood
#   glaucoma analysis (step 76h), restricted to the 1000G EUR LD-reference SNPs of each chromosome
#   (as step 76d does for chr11).
#   allele_corrected/  data/gwas/glaucoma_gwas_fusion.tsv, the allele-matched input of step 60b
#                      (the input reported from release v14)
#   as_reported/       data/gwas/glaucoma_gwas_fusion_firstrecord.tsv.gz, the first FinnGen record per
#                      rsID (step 60), the input reported up to release v13
#   Up to release v13 this step built allele_corrected/ itself with step 60b's rule. Since step 60b
#   also drops the 8 records that stay allele-inconsistent, chromosomes holding one of them can differ
#   from the v13 allele_corrected/ input in those rows only.
#   Also copies analysis/glaucoma_fusion_allele_matching.tsv (step 60b) to replaced_records.tsv.
# Usage: python3 76h_glaucoma_wb_genomewide_inputs.py <out dir>
# ============================================================
import gzip, os, shutil, sys

BASE = os.environ.get("FECD_BASE", os.path.expanduser("~/FECD_GWAS"))
os.chdir(BASE)
OUT = sys.argv[1] if len(sys.argv) > 1 else "results/twas_wb_genomewide/inputs"
LDREF = "data/ld_ref/LDREF/1000G.EUR"
SETS = {"allele_corrected": "data/gwas/glaucoma_gwas_fusion.tsv",
        "as_reported": "data/gwas/glaucoma_gwas_fusion_firstrecord.tsv.gz"}

CHR = {}
for c in range(1, 23):
    for line in open(f"{LDREF}.{c}.bim"):
        CHR[line.split()[1]] = c

for d, src in SETS.items():
    os.makedirs(f"{OUT}/{d}", exist_ok=True)
    fh = {c: open(f"{OUT}/{d}/glaucoma_fusion_chr{c}.tsv", "w") for c in range(1, 23)}
    n = 0
    with (gzip.open(src, "rt") if src.endswith(".gz") else open(src)) as f:
        hdr = f.readline()
        for o in fh.values():
            o.write(hdr)
        for line in f:
            c = CHR.get(line.split("\t", 1)[0])
            if c is not None:
                fh[c].write(line)
                n += 1
    for o in fh.values():
        o.close()
    print(f"{d}: {n} rows on the LD reference <- {src}")
shutil.copyfile("analysis/glaucoma_fusion_allele_matching.tsv", f"{OUT}/replaced_records.tsv")
