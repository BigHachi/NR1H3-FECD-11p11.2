#!/usr/bin/env Rscript
# ============================================================
# BK Analysis - Step 11: Prepare FECD GWAS (GCST90480062)
# Multi-ancestry meta-analysis: 3,928 cases / 566,393 controls
# ============================================================
suppressPackageStartupMessages(library(data.table))

BASE_DIR <- Sys.getenv("FECD_BASE", path.expand("~/FECD_GWAS"))
IN_FILE  <- file.path(BASE_DIR, "data/gwas/raw/GCST90480062.tsv.gz")
OUT_FILE <- file.path(BASE_DIR, "data/gwas/fecd_gwas_full.tsv.gz")

cat("Loading FECD GWAS:", IN_FILE, "\n")
fecd <- fread(IN_FILE, na.strings = "#NA")
cat(sprintf("  Raw SNPs: %d\n", nrow(fecd)))

# Rename columns to consistent internal names
setnames(fecd,
  old = c("chromosome","base_pair_location","effect_allele","other_allele",
          "odds_ratio","standard_error","effect_allele_frequency","p_value",
          "rsid","ci_upper","ci_lower"),
  new = c("CHR","BP","A1","A2","OR","SE","FREQ","P",
          "SNP","CI_UPPER","CI_LOWER"),
  skip_absent = TRUE)

# Compute BETA = log(OR); SE from 95% CI (all SE are NA in this file)
fecd[, BETA := log(OR)]
fecd[, SE := as.numeric(SE)]   # coerce logical NA column to numeric
fecd[is.na(SE) & !is.na(CI_UPPER) & !is.na(CI_LOWER),
     SE := (log(as.numeric(CI_UPPER)) - log(as.numeric(CI_LOWER))) / (2 * 1.96)]

# Filter
n_before <- nrow(fecd)
fecd <- fecd[
  !is.na(CHR) & !is.na(BP) &
  !is.na(BETA) & !is.na(SE) & SE > 0 &
  !is.na(P) & P > 0 & P <= 1 &
  !is.na(A1) & !is.na(A2) &
  CHR %in% 1:22
]
cat(sprintf("  After filtering: %d SNPs (removed %d)\n", nrow(fecd), n_before - nrow(fecd)))

fecd[, CHR := as.integer(CHR)]
fecd[, BP  := as.integer(BP)]

# panel_variant_id: chr{CHR}_{BP}_{A2}_{A1}_b38  (A2=ref, A1=effect)
fecd[, panel_variant_id := paste0("chr", CHR, "_", BP, "_", A2, "_", A1, "_b38")]

# Zscore
fecd[, ZSCORE := BETA / SE]

# Keep useful columns
fecd <- fecd[, .(SNP, CHR, BP, A1, A2, BETA, SE, OR, FREQ, P, ZSCORE,
                 panel_variant_id, n, num_cases, num_controls,
                 CI_UPPER, CI_LOWER, i2, q_pval, direction)]

setkey(fecd, CHR, BP)

# GW-significant summary
gw <- fecd[P < 5e-8]
cat(sprintf("\nGenome-wide significant SNPs (P<5e-8): %d\n", nrow(gw)))

fwrite(fecd, OUT_FILE, sep = "\t", quote = FALSE, compress = "gzip")
cat(sprintf("Saved: %s  (%d SNPs)\n", OUT_FILE, nrow(fecd)))
cat("\nTop 10 FECD hits:\n")
print(fecd[order(P)[1:10], .(SNP, CHR, BP, A1, A2, BETA, SE, P, FREQ)])
