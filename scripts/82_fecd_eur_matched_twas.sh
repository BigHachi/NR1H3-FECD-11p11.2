#!/bin/bash
# ============================================================
# Step 82: Ancestry-matched FECD FUSION TWAS
#   GWAS:    EUR-only FECD (GCST90475865, 3,473 EUR cases /
#            445,239 controls) — matches the EUR FUSION weights
#            and 1000G EUR LD, removing the multi-ancestry mismatch.
#   Weights: GTExv8.EUR (fibroblast, esophagus, skin)
#   Output:  results/twas_fecd_eurGWAS/{tissue}/...
#            analysis/fecd_eurGWAS_twas_{all,significant}.tsv
#   beta = ln(OR); SE reconstructed from CI (SE col is #NA in source)
# ============================================================
source "$(dirname "${BASH_SOURCE[0]}")/config.sh"   # FECD_BASE, FUSION_DIR, METAXCAN_DIR, METAXCAN_PY, PLINK, SMR, TABIX
set -euo pipefail
cd "$FECD_BASE"
LDREF=data/ld_ref/LDREF/1000G.EUR
WDIR=data/weights_eur/
GWASN=448712
RAW=data/gwas/fecd_eur_GCST90475865.tsv.gz
FUS=data/gwas/fecd_eur_fusion.tsv
OUTBASE=results/twas_fecd_eurGWAS
LOG=logs/82_fecd_eur_matched_twas.log
mkdir -p "$OUTBASE" logs
exec > >(tee "$LOG") 2>&1

# ---- 1. Build FUSION input (SNP A1 A2 Z N) ----
echo "=== Building FUSION input from EUR FECD sumstats ==="
if [ ! -f "$FUS" ]; then
Rscript - << 'EOF'
suppressPackageStartupMessages(library(data.table))
g <- fread("data/gwas/fecd_eur_GCST90475865.tsv.gz",
           select=c("rsid","effect_allele","other_allele","odds_ratio",
                    "ci_lower","ci_upper","p_value"))
setnames(g, c("rsid","effect_allele","other_allele","odds_ratio","ci_lower","ci_upper","p_value"),
            c("SNP","A1","A2","OR","CIL","CIU","P"))
g <- g[grepl("^rs", SNP)]
num <- function(x) suppressWarnings(as.numeric(x))
g[, OR:=num(OR)][, CIL:=num(CIL)][, CIU:=num(CIU)]
g <- g[!is.na(OR) & OR>0 & !is.na(CIL) & CIL>0 & !is.na(CIU) & CIU>0]
g[, beta := log(OR)]
g[, se   := (log(CIU) - log(CIL)) / (2*1.959964)]
g <- g[is.finite(beta) & is.finite(se) & se>0]
g[, Z := beta/se]
g <- g[is.finite(Z)][!duplicated(SNP)]
out <- g[, .(SNP, A1=toupper(A1), A2=toupper(A2), Z, N=448712L)]
fwrite(out, "data/gwas/fecd_eur_fusion.tsv", sep="\t", quote=FALSE)
cat(sprintf("  FUSION SNPs: %d\n", nrow(out)))
EOF
else
  echo "  (already built)"
fi
echo "  SNPs: $(($(wc -l < $FUS)-1))"

# ---- 1b. Split by chromosome (FUSION reads one chromosome at a time; keeps memory low) ----
# The FUSION file is an ordered subsequence of the raw GWAS, so chromosomes are assigned by
# streaming both files together.
if [ ! -s data/gwas/fecd_eur_fusion_bychr/chr22.tsv ]; then
python3 - "$RAW" "$FUS" data/gwas/fecd_eur_fusion_bychr <<'PYEOF'
import gzip, os, sys
raw, fus, outdir = sys.argv[1:4]
os.makedirs(outdir, exist_ok=True)
fr = gzip.open(raw, 'rt'); hr = fr.readline().rstrip('\n').split('\t'); iC, iR = hr.index('chromosome'), hr.index('rsid')
ff = open(fus); hdr = ff.readline()
outs = {}
def out(c):
    if c not in outs:
        outs[c] = open(os.path.join(outdir, f'chr{c}.tsv'), 'w'); outs[c].write(hdr)
    return outs[c]
line = ff.readline()
for r in fr:                      # the FUSION file is an ordered subsequence of the raw GWAS
    if not line: break
    x = r.split('\t')
    if x[iR] == line.split('\t', 1)[0]:
        out(x[iC]).write(line); line = ff.readline()
assert not line, f'unmatched FUSION row: {line[:60]}'
for f in outs.values(): f.close()
print('chromosomes written:', len(outs))
PYEOF
fi
# ---- 2. FUSION TWAS over 3 EUR tissues ----
for TISSUE in Cells_Cultured_fibroblasts Esophagus_Mucosa Skin_Not_Sun_Exposed_Suprapubic; do
  WPOS="${WDIR}GTExv8.EUR.${TISSUE}.pos"
  ODIR="${OUTBASE}/${TISSUE}"; mkdir -p "$ODIR"
  echo ""
  echo "--- ${TISSUE} ---"
  for CHR in $(seq 1 22); do
    OUTF="${ODIR}/fecd_eurGWAS_${TISSUE}_chr${CHR}"
    [ -s "${OUTF}" ] && { echo "  chr${CHR} exists, skip"; continue; }   # resume
    CHRSUM="data/gwas/fecd_eur_fusion_bychr/chr${CHR}.tsv"               # per-chr (low mem)
    Rscript ${FUSION_DIR}/FUSION.assoc_test.R \
      --sumstats ${CHRSUM} --weights ${WPOS} --weights_dir ${WDIR} \
      --ref_ld_chr ${LDREF}. --chr ${CHR} --out "${OUTF}" --GWASN ${GWASN} \
      > "logs/fecd_eurGWAS_${TISSUE}_chr${CHR}.log" 2>&1 || echo "  chr${CHR} warn"
  done
  COMB="${ODIR}/fecd_eurGWAS_${TISSUE}_all_chr.dat"
  head -1 "${ODIR}/fecd_eurGWAS_${TISSUE}_chr1" > "${COMB}"
  for CHR in $(seq 1 22); do tail -n +2 "${ODIR}/fecd_eurGWAS_${TISSUE}_chr${CHR}" >> "${COMB}" 2>/dev/null || true; done
  echo "  genes: $(($(wc -l < ${COMB})-1))"
done
echo ""
echo "=== Done: $OUTBASE ==="
