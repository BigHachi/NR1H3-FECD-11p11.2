#!/bin/bash
# ============================================================
# Step 82: Ancestry-matched FECD FUSION TWAS
#   GWAS:    EUR-only FECD (GCST90475865, 3,473 EUR cases /
#            445,239 controls) — matches the EUR FUSION weights
#            and 1000G EUR LD, removing the multi-ancestry mismatch.
#   Weights: GTExv8.EUR (fibroblast, esophagus, skin)
#   12 Nov 2026: the TWAS runs through scripts/fusion_genomewide.sh, the procedure used for every trait (steps 76h,
#   76i); the all-chromosome results equal the earlier per-chromosome loop's (models in the MHC were not in them
#   before either), and the genome-wide conditional/joint analysis of the Bonferroni-significant models is added.
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

# ---- 2. FUSION TWAS over 3 EUR tissues, with the genome-wide procedure shared by all traits ----
# (scripts/fusion_genomewide.sh: per-chromosome inputs restricted to the LD-reference SNPs, FUSION.assoc_test,
#  models in the MHC kept aside, genome-wide Bonferroni threshold, conditional/joint analysis of the significant models)
for TISSUE in Cells_Cultured_fibroblasts Esophagus_Mucosa Skin_Not_Sun_Exposed_Suprapubic; do
  echo ""
  echo "--- ${TISSUE} ---"
  bash "$SCRIPTS_DIR/fusion_genomewide.sh" "fecd_eurGWAS_${TISSUE}" "$FUS" "${WDIR}GTExv8.EUR.${TISSUE}.pos" "${OUTBASE}/${TISSUE}"
  echo "  genes: $(($(wc -l < ${OUTBASE}/${TISSUE}/fecd_eurGWAS_${TISSUE}_all_chr.dat)-1))"
done
echo ""
echo "=== Done: $OUTBASE ==="
