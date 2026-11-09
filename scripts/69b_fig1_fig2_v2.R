#!/usr/bin/env Rscript
# ============================================================
# Step 69b (Oct 9 2026; gene labels harmonised Oct 10 2026 via data/gene_symbol_display.tsv): Figures 1-2 regenerated from a single pipeline: FUSION on the European FECD GWAS (twas_fecd_eurGWAS) and the remapped 3-tissue S-MultiXcan (step 19b).
# Step 69: NR1H3-focused replots
#   (1) FECD FUSION TWAS Manhattan (fibroblasts) — NR1H3 labelled
#   (2) FECD S-MultiXcan volcano — NR1H3 labelled
#   (3) Cross-trait NR1H3 forest plot — myopia row uses JOINT Z/P
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
})

BASE <- Sys.getenv("FECD_BASE", path.expand("~/FECD_GWAS"))
ANA  <- file.path(BASE, "analysis")

NR1H3_ENSG <- "ENSG00000025434"

ensg_to_symbol <- function(v) {
  e <- sub("[.].*", "", v)
  m <- suppressMessages(select(org.Hs.eg.db, keys=unique(e), columns="SYMBOL", keytype="ENSEMBL"))
  m <- as.data.table(m)[!duplicated(ENSEMBL)]
  r <- m[match(e, ENSEMBL), SYMBOL]; r[is.na(r)] <- v[is.na(r)]; display_symbol(e, r)
}

# Oct 10 2026: one gene-name convention for Figs 1-3, Table 2 and S6
# (current HGNC symbol where one exists, else GENCODE v26 name, else Ensembl ID)
SYM_OVR <- fread(file.path(BASE, "data/gene_symbol_display.tsv"))
display_symbol <- function(ensg, sym) {
  e <- sub("[.].*", "", ensg)
  i <- match(e, SYM_OVR$ENSG)
  ifelse(is.na(i), sym, SYM_OVR$display_symbol[i])
}

# ============================================================
# (1) FECD FUSION TWAS Manhattan — EUR fibroblasts, NR1H3 labelled
# ============================================================
cat("=== (1) FECD FUSION TWAS Manhattan (EUR fibroblasts) ===\n")
twas <- fread(file.path(BASE, "results/twas_fecd_eurGWAS/Cells_Cultured_fibroblasts",
                        "fecd_eurGWAS_Cells_Cultured_fibroblasts_all_chr.dat"), fill=TRUE)
twas <- twas[!is.na(TWAS.P) & TWAS.P > 0]
twas[, ENSG   := sub("[.].*", "", ID)]
twas[, SYMBOL := ensg_to_symbol(ID)]

n_tests <- nrow(twas)
bonf_p  <- 0.05 / n_tests
bonf_y  <- -log10(bonf_p)

twas[, LOG10P := -log10(TWAS.P)]
twas[, CHR    := as.integer(CHR)]
setorder(twas, CHR, P0)
chr_off <- twas[, .(maxp = max(as.numeric(P0), na.rm=TRUE)), by=CHR]
chr_off[, offset := c(0, cumsum(as.numeric(head(maxp, -1))))]
twas <- merge(twas, chr_off[, .(CHR, offset)], by="CHR")
twas[, cum_pos := as.numeric(P0) + offset]
chr_mids <- twas[, .(mid = mean(cum_pos)), by=CHR][order(CHR)]

# Full (non-truncated) y-axis — original style, no clipping
nr <- twas[ENSG == NR1H3_ENSG]
cat(sprintf("  NR1H3: Z=%.2f  P=%.2e  -log10P=%.2f\n",
            nr$TWAS.Z[1], nr$TWAS.P[1], nr$LOG10P[1]))

sig    <- twas[TWAS.P < bonf_p]
y_top  <- max(twas$LOG10P, na.rm=TRUE) * 1.05
# 10 Nov 2026 (review): a broken y-axis (fig_common.R broken_y) at the largest gap in -log10 P keeps TCF4 and LINC01415 on the
# plot while giving the remaining loci room; gene symbols in italics; shared style and export (fig_common.R)
source(file.path(Sys.getenv("SCRIPTS_DIR", file.path(BASE, "scripts")), "fig_common.R"))
axis_break <- function(y) {        # (top of lower segment, bottom of upper segment) or NULL
  v <- sort(unique(y), decreasing = TRUE)[1:6]; g <- -diff(v); i <- which.max(g)
  if (g[i] < 40) return(NULL)
  c(ceiling((v[i + 1] + 4) / 5) * 5, floor((v[i] - 6) / 5) * 5)
}
brk <- axis_break(twas$LOG10P)
cat(sprintf("  y-axis break: %s\n", paste(brk, collapse = "-")))
sig[, face := gface(SYMBOL, bold = TRUE)]
lab1 <- rbind(sig[ENSG != NR1H3_ENSG, .(cum_pos, LOG10P, SYMBOL, face, key = "sig")],
              nr[, .(cum_pos, LOG10P, SYMBOL = "NR1H3", face = "bold.italic", key = "NR1H3")])
# background points rasterised at 600 dpi (ggrastr; text, labels and highlighted genes stay vector): ~9,500 vector
# points made the SVG too large for a reliable PDF
p_man <- ggplot(twas, aes(cum_pos, LOG10P, color=factor(CHR %% 2))) +
  ggrastr::rasterise(geom_point(size=0.7, alpha=0.6), dpi=600, dev="ragg") +
  geom_hline(yintercept=bonf_y, linetype="dashed", color="red", linewidth=0.5) +
  geom_point(data=sig, aes(cum_pos, LOG10P), color=COL_SIG, size=1.7) +
  geom_point(data=nr, aes(cum_pos, LOG10P), color=COL_NR1H3, size=2.6) +
  geom_text_repel(data=lab1, aes(cum_pos, LOG10P, label=SYMBOL, fontface=face, color=key),
                  size=2.8, family=FIG_FONT, min.segment.length=0, box.padding=0.5, point.padding=0.3,
                  max.overlaps=Inf, seed=1, max.time=Inf, force=3) +
  scale_color_manual(values=c("0"="steelblue","1"="navy", sig=COL_SIG, NR1H3=COL_NR1H3), guide="none") +
  scale_x_continuous(breaks=chr_mids$mid, labels=chr_mids$CHR, expand=expansion(mult=0.01),
                     guide=guide_axis(check.overlap=TRUE)) +
  labs(x="Chromosome", y="\u2212log<sub>10</sub>(*P*)",
       title="FUSION TWAS of FECD in cultured fibroblasts",
       subtitle=sprintf(paste0("GTEx v8 EUR weights; European-ancestry FECD GWAS. Dashed line: Bonferroni threshold, ",
                               "*P* < %s (%s genes)<br>%s highlighted",
                               if (!is.null(brk)) "; y-axis broken between %d and %d" else ""),
                        fmt_p(bonf_p, 3), format(n_tests, big.mark=","), gi("NR1H3"), brk[1], brk[2])) +
  theme_paper(9) + theme(axis.text.x=element_text(size=7.5), panel.grid.major.x=element_blank())
p_man <- if (!is.null(brk)) broken_y(p_man, brk, y_top, rel=0.35) else p_man + scale_y_continuous(limits=c(0, y_top))
save_fig(p_man, file.path(ANA, "fecd_twas_manhattan_NR1H3_v2"), width=FIG_W, height=4.0, fig="Figure1")
cat("  Saved: analysis/fecd_twas_manhattan_NR1H3.pdf\n\n")

# ============================================================
# (2) FECD S-MultiXcan volcano — NR1H3 labelled
# ============================================================
cat("=== (2) FECD S-MultiXcan volcano ===\n")
smx <- fread(file.path(BASE, "results/multixcan/fecd_smultixcan_remapped.txt"))[!is.na(pvalue)]; smx[, ENSG := sub("[.].*", "", gene)]; smx[, SYMBOL := display_symbol(ENSG, gene_name)]
n_smx  <- nrow(smx)
bonf_s <- 0.05 / n_smx
smx[, log10p := -log10(pvalue)]
smx[, significant := !is.na(pvalue) & pvalue < bonf_s]   # !is.na avoids an NA legend entry

nrs <- smx[ENSG == NR1H3_ENSG]
cat(sprintf("  NR1H3: z_mean=%.2f  P=%.2e  -log10P=%.2f\n",
            nrs$z_mean[1], nrs$pvalue[1], nrs$log10p[1]))
# Oct 21 2026 (review 18): ZNF408, the other 11p11.2 candidate (stronger S-MultiXcan P), highlighted too,
# and the x-axis labelled explicitly as the mean of the per-tissue S-PrediXcan Z-scores
ZNF408_ENSG <- "ENSG00000175213"
zns <- smx[ENSG == ZNF408_ENSG]
cat(sprintf("  ZNF408: z_mean=%.2f  P=%.2e\n", zns$z_mean[1], zns$pvalue[1]))

smx[, sig_lab := factor(ifelse(significant, "Bonferroni-significant", "Not significant"),
                        levels=c("Bonferroni-significant", "Not significant"))]
lab_s <- rbind(smx[significant==TRUE & !ENSG %in% c(NR1H3_ENSG, ZNF408_ENSG),
                   .(z_mean, log10p, SYMBOL, face = gface(SYMBOL), key = "sig")],
               nrs[, .(z_mean, log10p, SYMBOL = "NR1H3", face = "bold.italic", key = "NR1H3")],
               zns[, .(z_mean, log10p, SYMBOL = "ZNF408", face = "bold.italic", key = "ZNF408")])
brk2 <- axis_break(smx$log10p)
cat(sprintf("  y-axis break: %s\n", paste(brk2, collapse = "-")))
p_volc <- ggplot(smx, aes(z_mean, log10p)) +
  ggrastr::rasterise(geom_point(aes(color=sig_lab), size=0.8, alpha=0.6), dpi=600, dev="ragg") +   # 17,285 points; see Figure 1
  geom_hline(yintercept=-log10(bonf_s), linetype="dashed", color="red", linewidth=0.5) +
  geom_point(data=nrs, aes(z_mean, log10p), color=COL_NR1H3, size=2.6) +
  geom_point(data=zns, aes(z_mean, log10p), color=COL_ZNF408, size=2.6) +
  geom_text_repel(seed=1, max.time=Inf, data=lab_s, aes(z_mean, log10p, label=SYMBOL, fontface=face, color=key),
                  size=3, family=FIG_FONT, max.overlaps=Inf, min.segment.length=0, box.padding=0.6,
                  point.padding=0.3, force=3, show.legend=FALSE) +
  scale_color_manual(values=c("Not significant"="grey65", "Bonferroni-significant"=COL_SIG, sig=COL_SIG,
                              NR1H3=COL_NR1H3, ZNF408=COL_ZNF408),
                     breaks=c("Bonferroni-significant", "Not significant"), name=NULL) +
  scale_x_continuous(labels=lab_minus) +
  labs(x="Mean per-tissue S-PrediXcan Z-score (three tissues)", y="\u2212log<sub>10</sub>(*P*)",
       title="S-MultiXcan of FECD across three proxy tissues",
       subtitle=sprintf(paste0("GTEx v8 MASHR models; multi-ancestry FECD GWAS. Dashed line: Bonferroni threshold, ",
                               "*P* < %s (%s genes)<br>11p11.2: %s (mean Z = %s, *P* = %s), ",
                               "%s (mean Z = %s, *P* = %s)", if (!is.null(brk2)) "<br>y-axis broken between %d and %d" else ""),
                        fmt_p(bonf_s, 3), format(n_smx, big.mark=","), gi("ZNF408"), fmt_num(zns$z_mean[1]),
                        fmt_p(zns$pvalue[1]), gi("NR1H3"), fmt_num(nrs$z_mean[1]), fmt_p(nrs$pvalue[1]), brk2[1], brk2[2])) +
  theme_paper(10) + theme(legend.position="bottom")
if (!is.null(brk2)) p_volc <- broken_y(p_volc, brk2, max(smx$log10p) * 1.02, rel=0.22,
                                       upper_breaks=seq(ceiling(brk2[2] / 5) * 5, 200, 5))
save_fig(p_volc, file.path(ANA, "fecd_smultixcan_volcano_NR1H3_v2"), width=FIG_W, height=6.4, fig="Figure2")
cat("  Saved: analysis/fecd_smultixcan_volcano_NR1H3.pdf\n\n")
