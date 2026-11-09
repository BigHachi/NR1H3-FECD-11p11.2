#!/bin/bash
# ============================================================
# Glaucoma Analysis - Step 45: Download FinnGen R12 H7_GLAUCOMA
# and convert to standard analysis format
# ============================================================
source "$(dirname "${BASH_SOURCE[0]}")/config.sh"   # FECD_BASE, FUSION_DIR, METAXCAN_DIR, METAXCAN_PY, PLINK, SMR, TABIX
set -e
cd "$FECD_BASE"

GWAS_RAW=data/gwas/finngen_R12_H7_GLAUCOMA.gz
GWAS_FULL=data/gwas/glaucoma_gwas_full.tsv

echo "=== Downloading FinnGen R12 Glaucoma GWAS ==="
if [ ! -f "$GWAS_RAW" ]; then
  curl -L --fail --retry 3 -o "$GWAS_RAW" \
    "https://storage.googleapis.com/finngen-public-data-r12/summary_stats/release/finngen_R12_H7_GLAUCOMA.gz"
  echo "  Downloaded: $(du -sh $GWAS_RAW | cut -f1)"
else
  echo "  Already downloaded: $GWAS_RAW"
fi

echo ""
echo "=== Converting FinnGen R12 GWAS to analysis format ==="
Rscript - <<'EOF'
suppressPackageStartupMessages(library(data.table))

cat("Reading FinnGen R12 Glaucoma GWAS...\n")
gwas <- fread("data/gwas/finngen_R12_H7_GLAUCOMA.gz")
cat(sprintf("  Raw rows: %d\n", nrow(gwas)))
cat("  Columns:", paste(names(gwas), collapse=", "), "\n")

# FinnGen columns: #chrom, pos, ref, alt, rsids, nearest_genes,
#                  pval, mlogp, beta, sebeta, af_alt, ...
gwas <- gwas[!is.na(pval) & !is.na(beta) & !is.na(sebeta) & sebeta > 0]
gwas <- gwas[grepl("^rs", rsids)]
gwas[, CHR := as.integer(sub("^chr", "", `#chrom`))]
gwas <- gwas[!is.na(CHR) & CHR %in% 1:22]

gwas_out <- gwas[, .(
  SNP  = rsids,
  CHR  = CHR,
  BP   = pos,
  A1   = toupper(alt),   # effect allele (FinnGen convention)
  A2   = toupper(ref),   # non-effect allele
  BETA = beta,
  SE   = sebeta,
  P    = pval,
  FREQ = af_alt
)]
gwas_out <- gwas_out[!duplicated(SNP)]
gwas_out <- gwas_out[is.finite(BETA) & is.finite(SE)]

cat(sprintf("  SNPs after filtering: %d\n", nrow(gwas_out)))
fwrite(gwas_out, "data/gwas/glaucoma_gwas_full.tsv", sep="\t", quote=FALSE)
cat("  Saved: data/gwas/glaucoma_gwas_full.tsv\n")
EOF

echo ""
echo "=== Done ==="
echo "Output: $GWAS_FULL"
echo "Lines:  $(wc -l < $GWAS_FULL)"
