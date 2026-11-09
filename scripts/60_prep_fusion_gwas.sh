#!/bin/bash
# ============================================================
# Step 60: Prepare GWAS files for FUSION TWAS
# Converts glaucoma, IOP, and myopia GWAS to FUSION format:
#   SNP  A1  A2  Z  N
# ============================================================
source "$(dirname "${BASH_SOURCE[0]}")/config.sh"   # FECD_BASE, FUSION_DIR, METAXCAN_DIR, METAXCAN_PY, PLINK, SMR, TABIX
set -e
cd "$FECD_BASE"

echo "=== Preparing GWAS files for FUSION TWAS ==="

Rscript - << 'EOF'
suppressPackageStartupMessages(library(data.table))

prep_fusion <- function(in_file, out_file, N, label) {
  cat(sprintf("\n--- %s ---\n", label))
  gwas <- fread(in_file, select = c("SNP","A1","A2","BETA","SE"))
  cat(sprintf("  Loaded: %d rows\n", nrow(gwas)))

  gwas <- gwas[grepl("^rs", SNP) & SE > 0 & is.finite(BETA)]
  gwas[, Z := BETA / SE]
  gwas <- gwas[is.finite(Z)][!duplicated(SNP)]

  out <- gwas[, .(SNP, A1 = toupper(A1), A2 = toupper(A2), Z, N = as.integer(N))]
  fwrite(out, out_file, sep = "\t", quote = FALSE)
  cat(sprintf("  Written: %d SNPs -> %s\n", nrow(out), out_file))
}

# Glaucoma: first FinnGen record per rsID, as used up to release v13. Step 60b replaces the records
# that do not match the LD-reference alleles and writes data/gwas/glaucoma_gwas_fusion.tsv (10 Nov 2026).
prep_fusion(
  "data/gwas/glaucoma_gwas_full.tsv",
  "data/gwas/glaucoma_gwas_fusion_firstrecord.tsv.gz",   # fwrite gzips by extension
  N     = 500348,   # 26,591 cases + 473,757 controls (FinnGen R12)
  label = "Glaucoma (FinnGen R12)"
)

prep_fusion(
  "data/gwas/iop_gwas_full.tsv",
  "data/gwas/iop_gwas_fusion.tsv",
  N     = 31269,    # Bonnemaijer 2019 European
  label = "IOP (Bonnemaijer 2019)"
)

# Myopia: BETA = Zscore, SE = 1 (see step 55/56); Z = BETA/SE = Zscore
prep_fusion(
  "data/gwas/myopia_gwas_full.tsv",
  "data/gwas/myopia_gwas_fusion.tsv",
  N     = 351091,   # Hysi 2020 European excl. 23andMe
  label = "Myopia (Hysi 2020)"
)
EOF

echo ""
echo "=== FUSION GWAS prep complete ==="
for f in data/gwas/glaucoma_gwas_fusion_firstrecord.tsv.gz data/gwas/iop_gwas_fusion.tsv data/gwas/myopia_gwas_fusion.tsv; do
  echo "  $(basename $f): $(gzip -cdf $f | wc -l) SNPs"
done
