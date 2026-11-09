#!/usr/bin/env Rscript
# ============================================================
# Glaucoma Analysis - Step 49: Analyse S-PrediXcan + S-MultiXcan
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
})

BASE_DIR <- Sys.getenv("FECD_BASE", path.expand("~/FECD_GWAS"))
ANA_DIR  <- file.path(BASE_DIR, "analysis")
dir.create(ANA_DIR, showWarnings = FALSE, recursive = TRUE)

ensg_to_symbol <- function(ensg_vec) {
  ensg_clean <- sub("[.].*", "", ensg_vec)
  m <- suppressMessages(
    select(org.Hs.eg.db, keys = unique(ensg_clean),
           columns = "SYMBOL", keytype = "ENSEMBL")
  )
  m <- as.data.table(m)[!duplicated(ENSEMBL)]
  res <- m[match(ensg_clean, ENSEMBL), SYMBOL]
  res[is.na(res)] <- ensg_vec[is.na(res)]
  res
}

# ============================================================
# PART A: S-MultiXcan
# ============================================================
cat("======================================================\n")
cat("PART A: S-MultiXcan Analysis\n")
cat("======================================================\n\n")

smx_file <- file.path(BASE_DIR, "results/multixcan/glaucoma_smultixcan.txt")
if (!file.exists(smx_file)) stop("Run 48_run_glaucoma_smultixcan.sh first.")

smx <- fread(smx_file)
cat(sprintf("S-MultiXcan results: %d genes\n", nrow(smx)))

smx[, ENSG   := sub("[.].*", "", gene)]
smx[, SYMBOL := ensg_to_symbol(ENSG)]

n_tests_smx <- nrow(smx)
bonf_p_smx  <- 0.05 / n_tests_smx
cat(sprintf("Bonferroni threshold: P < %.2e (%d tests)\n", bonf_p_smx, n_tests_smx))

smx[, FDR := p.adjust(pvalue, method = "BH")]
sig_smx <- smx[pvalue < bonf_p_smx][order(pvalue)]
cat(sprintf("Significant hits: %d\n\n", nrow(sig_smx)))
if (nrow(sig_smx) > 0) print(sig_smx[, .(SYMBOL, ENSG, pvalue, FDR, n_models_used = n)])

fwrite(smx[order(pvalue)],
       file.path(ANA_DIR, "glaucoma_smultixcan_all.tsv"),    sep = "\t", quote = FALSE)
fwrite(sig_smx,
       file.path(ANA_DIR, "glaucoma_smultixcan_significant.tsv"), sep = "\t", quote = FALSE)

# Volcano
smx[, log10p      := -log10(pvalue)]
smx[, significant := pvalue < bonf_p_smx]

p_volcano <- ggplot(smx, aes(x = z_mean, y = log10p, color = significant)) +
  geom_point(size = 0.8, alpha = 0.6) +
  geom_hline(yintercept = -log10(bonf_p_smx), linetype = "dashed", color = "red") +
  scale_color_manual(values = c("FALSE" = "grey60", "TRUE" = "firebrick"),
                     name = "Bonferroni sig.") +
  labs(x = "Mean Z-score", y = "-log10(P)",
       title = "S-MultiXcan Volcano — Glaucoma (FinnGen R12)",
       subtitle = sprintf("GTEx v8 MASHR | 3 tissues | Bonferroni P < %.2e", bonf_p_smx)) +
  theme_bw(base_size = 11)

if (nrow(sig_smx) > 0) {
  p_volcano <- p_volcano +
    geom_text_repel(seed = 1, max.time = Inf, data = smx[significant == TRUE],
                    aes(label = SYMBOL), size = 2.5, color = "firebrick", max.overlaps = 20)
}
ggsave(file.path(ANA_DIR, "glaucoma_smultixcan_volcano.pdf"), p_volcano, width = 8, height = 6)
cat("Volcano: ", file.path(ANA_DIR, "glaucoma_smultixcan_volcano.pdf"), "\n\n")

# ============================================================
# PART B: S-PrediXcan per tissue
# ============================================================
cat("======================================================\n")
cat("PART B: S-PrediXcan Per-Tissue Analysis\n")
cat("======================================================\n\n")

tissue_files <- list(
  Esophagus_Mucosa                = file.path(BASE_DIR, "results/sprediXcan/glaucoma_esophagus_mucosa.csv"),
  Skin_Not_Sun_Exposed_Suprapubic = file.path(BASE_DIR, "results/sprediXcan/glaucoma_skin_not_sun_exposed_suprapubic.csv"),
  Cells_Cultured_fibroblasts      = file.path(BASE_DIR, "results/sprediXcan/glaucoma_cells_cultured_fibroblasts.csv")
)

spx_all <- list()
spx_sig <- list()

for (tissue in names(tissue_files)) {
  fpath <- tissue_files[[tissue]]
  if (!file.exists(fpath)) { cat(sprintf("WARNING: %s not found — skipping\n", fpath)); next }

  spx <- fread(fpath)
  spx[, TISSUE := tissue]
  spx[, ENSG   := sub("[.].*", "", gene)]
  spx[, SYMBOL := ensg_to_symbol(ENSG)]

  bonf_p <- 0.05 / nrow(spx)
  cat(sprintf("[%s]\n  Genes: %d | Bonferroni P < %.2e\n", tissue, nrow(spx), bonf_p))

  spx[, FDR := p.adjust(pvalue, method = "BH")]
  sig <- spx[pvalue < bonf_p][order(pvalue)]
  cat(sprintf("  Significant: %d\n", nrow(sig)))
  if (nrow(sig) > 0) print(sig[, .(SYMBOL, ENSG, zscore = zscore, pvalue, FDR)])
  cat("\n")

  spx_all[[tissue]] <- spx
  spx_sig[[tissue]] <- sig
}

if (length(spx_all) > 0) {
  all_spx_dt <- rbindlist(spx_all, fill = TRUE)
  sig_spx_dt <- rbindlist(spx_sig, fill = TRUE)
  fwrite(all_spx_dt[order(TISSUE, pvalue)],
         file.path(ANA_DIR, "glaucoma_sprediXcan_all.tsv"),         sep = "\t", quote = FALSE)
  fwrite(sig_spx_dt,
         file.path(ANA_DIR, "glaucoma_sprediXcan_significant.tsv"), sep = "\t", quote = FALSE)
  cat("All S-PrediXcan:  ", file.path(ANA_DIR, "glaucoma_sprediXcan_all.tsv"), "\n")
  cat("Sig S-PrediXcan:  ", file.path(ANA_DIR, "glaucoma_sprediXcan_significant.tsv"), "\n\n")
}

# ============================================================
# PART C: Manhattan plot (S-MultiXcan)
# ============================================================
cat("======================================================\n")
cat("PART C: Manhattan Plot\n")
cat("======================================================\n\n")

# Need gene coordinates — pull from S-PrediXcan fibroblasts output (has CHR/pos)
fib_file <- file.path(BASE_DIR, "results/sprediXcan/glaucoma_cells_cultured_fibroblasts.csv")
if (file.exists(fib_file)) {
  # Use the per-tissue file that has gene_cov/chromosome columns
  fib  <- fread(fib_file)
  # Merge chromosome position from S-PrediXcan output into S-MultiXcan
  if ("chromosome" %in% names(fib)) {
    smx_plot <- merge(smx,
                      fib[, .(gene, CHR = as.integer(chromosome), P0 = start)],
                      by = "gene", all.x = TRUE)
  } else {
    smx_plot <- smx
    smx_plot[, CHR := NA_integer_]
    smx_plot[, P0  := NA_real_]
  }
} else {
  smx_plot <- smx
  smx_plot[, CHR := NA_integer_]
  smx_plot[, P0  := NA_real_]
}

# If no positional info, skip Manhattan
if (all(is.na(smx_plot$CHR))) {
  cat("No positional data available for Manhattan — skipping\n\n")
} else {
  smx_plot <- smx_plot[!is.na(CHR) & !is.na(P0)]
  setorder(smx_plot, CHR, P0)
  chr_offsets <- smx_plot[, .(max_pos = max(as.numeric(P0), na.rm=TRUE)), by=CHR]
  chr_offsets[, offset := c(0, cumsum(as.numeric(head(max_pos, -1))))]
  smx_plot <- merge(smx_plot, chr_offsets[, .(CHR, offset)], by="CHR")
  smx_plot[, cum_pos := as.numeric(P0) + offset]
  chr_mids <- smx_plot[, .(mid = mean(cum_pos)), by=CHR][order(CHR)]

  p_man <- ggplot(smx_plot, aes(x = cum_pos, y = log10p, color = factor(CHR %% 2))) +
    geom_point(size = 0.8, alpha = 0.6) +
    geom_hline(yintercept = -log10(bonf_p_smx), linetype = "dashed", color = "red") +
    scale_color_manual(values = c("0" = "steelblue", "1" = "navy"), guide = "none") +
    scale_x_continuous(breaks = chr_mids$mid, labels = chr_mids$CHR,
                       guide = guide_axis(check.overlap = FALSE)) +
    labs(x = "Chromosome", y = "-log10(P)",
         title = "Glaucoma S-MultiXcan Manhattan Plot (GTEx v8 MASHR)",
         subtitle = sprintf("FinnGen R12 H7_GLAUCOMA | Bonferroni P < %.2e | %d genes | %d significant",
                            bonf_p_smx, n_tests_smx, nrow(sig_smx))) +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_text(size = 7))

  if (nrow(sig_smx) > 0) {
    sig_plot <- smx_plot[pvalue < bonf_p_smx]
    p_man <- p_man +
      geom_text_repel(seed = 1, max.time = Inf, data = sig_plot, aes(label = SYMBOL),
                      size = 2.5, color = "red", max.overlaps = 20)
  }

  man_pdf <- file.path(ANA_DIR, "glaucoma_smultixcan_manhattan.pdf")
  ggsave(man_pdf, p_man, width = 14, height = 5)
  cat("Manhattan: ", man_pdf, "\n\n")
}

cat("=== Analysis complete ===\n")
cat(sprintf("S-MultiXcan significant genes: %d\n", nrow(sig_smx)))
if (nrow(sig_smx) > 0) cat(paste(unique(sig_smx$SYMBOL), collapse = ", "), "\n")
