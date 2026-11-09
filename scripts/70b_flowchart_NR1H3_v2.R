#!/usr/bin/env Rscript
# ============================================================
# Step 70b: Figure 7, flow diagram of the FECD analyses (blue) and the glaucoma cross-trait analyses (yellow).
# History: Oct 9-10 2026 statistics and wording updated to the review reruns; Oct 21 2026 glaucoma cis-MR read
#   from the Table 3 output (step 89c).
# 10 Nov 2026 (review: "increase the font size"): rebuilt on a 10 x 12 in canvas with 12.5-14 pt text (about
#   9-10 pt when printed at a 7-in page width; before, 7-10 pt text on a 10 x 11 in canvas printed at 5-7 pt), shorter
#   lines in taller boxes, one arrow per box, gene symbols in italics and P-values as a x 10^b (fig_common.R).
#   Every statistic is read from the pipeline outputs: S-MultiXcan (steps 19b, 75), FECD cis-MR (step 91) and
#   glaucoma cis-MR (step 89c).
#   Output: analysis/fig7_flowchart_NR1H3_v2.{pdf,svg,png} (+ Figure7.* when FIG_DIR is set)
# ============================================================
suppressPackageStartupMessages({ library(ggplot2); library(data.table) })
BASE <- Sys.getenv("FECD_BASE", path.expand("~/FECD_GWAS")); ANA <- file.path(BASE, "analysis")
source(file.path(Sys.getenv("SCRIPTS_DIR", file.path(BASE, "scripts")), "fig_common.R"))

ens <- function(x) sub("[.].*", "", x)
smx <- fread(file.path(BASE, "results/multixcan/fecd_smultixcan_remapped.txt"))
P_NR1H3  <- smx[ens(gene) == "ENSG00000025434", pvalue]
P_ZNF408 <- smx[ens(gene) == "ENSG00000175213", pvalue]
gsm <- fread(file.path(BASE, "results/multixcan/glaucoma_smultixcan_brain.txt"))
P_GL_SMX <- gsm[ens(gene) == "ENSG00000025434", pvalue]
fe <- fread(file.path(ANA, "nr1h3_eqtlgen_fecd_mr_EUR.tsv"))[method == "Inverse variance weighted"]
fd <- fread(file.path(ANA, "nr1h3_eqtlgen_fecd_mr_EUR_diagnostics.tsv"))
t3 <- fread(file.path(ANA, "table3_nr1h3_crosstrait_cismr.tsv"))[outcome == "Glaucoma" & method == "Inverse variance weighted"]
stopifnot(length(P_NR1H3) == 1, length(P_ZNF408) == 1, length(P_GL_SMX) == 1, nrow(fe) == 1, nrow(t3) == 1)
GL_OR <- sprintf("%.2f", exp(t3$b)); FE_OR <- sprintf("%.2f", exp(fe$b))

BLUE <- c("#D9E2F3", "#3B5C9F"); YEL <- c("#FBF1D0", "#BF9B30"); RED <- c("#F8DDDD", "#B23A3A")
TITLE_PT <- 14; BODY_PT <- 12.5
it <- function(g) paste0("*", g, "*")
box <- function(x0, x1, y1, h, col, title, body) data.table(xmin = x0, xmax = x1, ymax = y1, ymin = y1 - h,
  fill = col[1], border = col[2], label = paste0("<span style='font-size:", TITLE_PT, "pt'>**", title, "**</span><br>",
                                                 paste(body, collapse = "<br>")))
# rows (inches from the bottom of the canvas): section label at T[k], box top at T[k] - 0.3, box height h[k]
H <- c(1.25, 1.0, 1.25, 1.25, 1.5, 1.0); CANVAS_H <- 0.25 + sum(0.3 + H + 0.55) - 0.55 + 0.1
R <- numeric(6); R[1] <- CANVAS_H - 0.25
for (k in 2:6) R[k] <- R[k - 1] - 0.3 - H[k - 1] - 0.55
top <- R - 0.3; bot <- top - H
b <- rbind(
  # (i) data
  box(0.10, 2.10, top[1], H[1], BLUE, "MVP FECD GWAS", c("Verma et al. 2024", "3,928 cases")),
  box(2.25, 4.25, top[1], H[1], BLUE, "GTEx v8", c("MASHR and FUSION", "weights, 3 proxy tissues", "~500 donors each")),
  box(4.40, 6.40, top[1], H[1], BLUE, "1000 Genomes", c("LD reference (EUR)", "phase 3; n = 503")),
  box(6.65, 8.225, top[1], H[1], YEL,  "FinnGen R12", c("H7_GLAUCOMA", "26,591 cases")),
  box(8.325, 9.90, top[1], H[1], YEL,  "GTEx v8", c("MASHR models:", "brain cortex,", "fibroblasts, skin")),
  # (ii) preprocessing
  box(1.00, 5.50, top[2], H[2], BLUE, "GWAS preprocessing", c("custom scripts; allele matching to", "model and LD-reference variants")),
  box(6.65, 9.90, top[2], H[2], YEL,  "Glaucoma TWAS", c("S-PrediXcan + S-MultiXcan", "brain cortex, fibroblasts, skin")),
  # (iii) TWAS
  box(0.10, 3.15, top[3], H[3], BLUE, "S-PrediXcan + S-MultiXcan", c("MetaXcan v0.8.1; MASHR models",
      paste0(it("NR1H3"), " *P* = ", fmt_p(P_NR1H3)), paste0(it("ZNF408"), " *P* = ", fmt_p(P_ZNF408)))),
  box(3.35, 6.40, top[3], H[3], BLUE, "FUSION TWAS", c("Gusev et al. 2016", "cultured fibroblasts", "(GTEx v8 EUR weights)")),
  box(6.65, 9.90, top[3], H[3], YEL,  paste(it("NR1H3"), "in glaucoma"), c("S-MultiXcan (cross-tissue)",
      paste0("*P* = ", fmt_p(P_GL_SMX), " with brain cortex;"), "not significant without")),
  # (iv) prioritisation
  box(0.60, 5.90, top[4], H[4], BLUE, "Conditional analysis and fine-mapping",
      c(paste0("FUSION: ", it("NR1H3"), " independent in fibroblasts only"),
        paste0("(no ", it("ZNF408"), " model); SuSiE: no credible set;"), "FOCUS: effector gene unresolved")),
  # (v) colocalisation and causal inference
  box(0.10, 3.15, top[5], H[5], BLUE, "Colocalisation", c("coloc.abf, full cis-eQTLs", "(eQTL Catalogue)",
      paste0(it("NR1H3"), " PP4 \u2264 0.23"), paste0(it("ZNF408"), " PP4 \u2264 0.60"))),
  box(3.35, 6.40, top[5], H[5], BLUE, "FECD cis-MR", c("eQTLGen blood", sprintf("N = 31,684; F = %.0f", fd$F_set),
      "European FECD GWAS:", paste0("OR ", FE_OR, ", *P* = ", fmt_p(fe$pval), " (null)"))),
  box(6.65, 9.90, top[5], H[5], YEL,  "Glaucoma cis-MR", c("eQTLGen blood, 4 instruments",
      paste0("IVW OR ", GL_OR, ", *P* = ", fmt_p(t3$pval)))),
  # (vi) result
  box(0.10, 9.90, top[6], H[6], RED, "11p11.2: suggestive FECD locus, effector gene unresolved",
      c("FECD: sub-threshold TWAS, no colocalisation, fine-mapping unresolved, MR null, no association in FinnGen H18.5",
        paste0("Glaucoma: cis-MR OR ", GL_OR, "; little support for a causal variant shared with FECD")))
)
b[, `:=`(cx = (xmin + xmax) / 2, ty = ymax - 0.12)]

sect <- data.table(x = 0.1, y = R, lab = c("(i) Study design and data acquisition", "(ii) Preprocessing and harmonisation",
  "(iii) Statistical analysis: TWAS discovery", "(iv) Locus prioritisation", "(v) Colocalisation and causal inference",
  "(vi) Result"))
head_r <- data.table(x = 8.275, y = R[1], lab = "Cross-trait analysis: glaucoma")

# arrows: x, from y, to y (one per target box)
L1 <- 1.625; L2 <- 4.875; RC <- 8.275
# the two glaucoma data boxes (equal widths either side of the column's centre line) join in a bracket that feeds
# one arrow on the centre line, so every arrow in the glaucoma column is centred on its boxes
GX <- c((6.65 + 8.225) / 2, (8.325 + 9.90) / 2); JY <- bot[1] - 0.25
join <- data.table(x = c(GX, GX[1]), xend = c(GX, GX[2]), y = c(bot[1], bot[1], JY), yend = c(JY, JY, JY))
arr <- data.table(
  x    = c(3.25, L1, L2, L1, L2, L1, L2, L1, L2, RC, RC, RC, RC),
  y    = c(bot[1], bot[2], bot[2], bot[3], bot[3], bot[4], bot[4], bot[5], bot[5], JY, bot[2], bot[3], bot[5]),
  yend = c(top[2], top[3], top[3], top[4], top[4], top[5], top[5], top[6], top[6], top[2], top[3], top[5], top[6]))
arr[, yend := yend + 0.03]

p <- ggplot() +
  geom_rect(data = b, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = b$fill, colour = b$border,
            linewidth = 0.6) +
  geom_segment(data = join, aes(x = x, y = y, xend = xend, yend = yend), linewidth = 0.5, colour = "grey25",
               lineend = "square") +
  geom_segment(data = arr, aes(x = x, y = y, xend = x, yend = yend), linewidth = 0.5, colour = "grey25",
               arrow = arrow(length = unit(0.16, "cm"), type = "closed")) +
  geom_richtext(data = b, aes(cx, ty, label = label), vjust = 1, size = BODY_PT / .pt, family = FIG_FONT,
                lineheight = 1.25, fill = NA, label.colour = NA, label.padding = unit(0, "pt"), colour = "grey10") +
  geom_label(data = sect, aes(x, y, label = lab), hjust = 0, vjust = 0.5, fontface = "bold", size = TITLE_PT / .pt,
             family = FIG_FONT, fill = "white", linewidth = 0, label.padding = unit(2, "pt")) +   # arrows pass behind
  geom_text(data = head_r, aes(x, y, label = lab), vjust = 0.5, fontface = "bold", size = TITLE_PT / .pt,
            family = FIG_FONT) +
  coord_cartesian(xlim = c(0, 10), ylim = c(0, CANVAS_H), expand = FALSE) +
  theme_void(base_family = FIG_FONT) +
  theme(plot.margin = margin(4, 4, 4, 4), plot.background = element_rect(fill = "white", colour = NA))
save_fig(p, file.path(ANA, "fig7_flowchart_NR1H3_v2"), width = 10, height = CANVAS_H, fig = "Figure7", png_dpi = 200)
cat("Saved: analysis/fig7_flowchart_NR1H3_v2.{pdf,svg,png} (Figure 7)\n")
