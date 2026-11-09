#!/bin/bash
# ============================================================
# FECD Analysis - Step 15: Prepare FECD GWAS for FUSION
# ============================================================
source "$(dirname "${BASH_SOURCE[0]}")/config.sh"   # FECD_BASE, FUSION_DIR, METAXCAN_DIR, METAXCAN_PY, PLINK, SMR, TABIX
set -e
cd "$FECD_BASE"

echo "=== Preparing FECD GWAS for FUSION ==="

Rscript - <<'EOF'
suppressPackageStartupMessages(library(data.table))

BASE_DIR <- Sys.getenv("FECD_BASE", path.expand("~/FECD_GWAS"))

cat("Loading FECD GWAS (fecd_gwas_full.tsv.gz)...\n")
fecd <- fread(file.path(BASE_DIR, "data/gwas/fecd_gwas_full.tsv.gz"),
              select = c("SNP", "A1", "A2", "ZSCORE"))
cat(sprintf("  Raw rows: %d\n", nrow(fecd)))

# FUSION format: SNP, A1, A2, Z, N
# A1 = effect allele, A2 = ref allele; ZSCORE already computed in script 11
fecd <- fecd[!is.na(ZSCORE) & is.finite(ZSCORE) & grepl("^rs", SNP)]
fecd <- fecd[!duplicated(SNP)]

N_EFF <- 570321L
fusion <- data.table(
  SNP = fecd$SNP,
  A1  = toupper(fecd$A1),
  A2  = toupper(fecd$A2),
  Z   = fecd$ZSCORE,
  N   = N_EFF
)
fusion <- fusion[is.finite(Z)]

out <- file.path(BASE_DIR, "data/gwas/fecd_gwas_fusion.tsv")
fwrite(fusion, out, sep = "\t", quote = FALSE)
cat(sprintf("  FUSION format: %d SNPs saved to %s\n", nrow(fusion), out))
EOF

echo ""
echo "=== FECD FUSION prep complete ==="
echo "Output: $(wc -l < data/gwas/fecd_gwas_fusion.tsv) lines"
