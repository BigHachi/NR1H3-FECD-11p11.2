#!/usr/bin/env Rscript
# ============================================================
# Step 89b: Figure 6 — NR1H3 eQTLGen cis-MR on glaucoma, IOP and refractive error.
#   10 Nov 2026 (review): each outcome is drawn on its own scale, from the Table 3 output of step 89c, so the
#   figure and Table 3 share instruments, estimators and values:
#     glaucoma          odds ratio (exp of the Table 3 log-OR), log axis
#     IOP               mmHg (GWAS beta/SE as published)
#     refractive error  SD units: the Hysi et al. statistics are Z-scores only, so step 89c standardises the
#                       outcome from Z with the per-variant N and the outcome-GWAS allele frequency
#   all per SD of genetically predicted NR1H3 expression (eQTLGen Z standardised with the eQTLGen allele
#   frequency and per-variant N, step 89c). 95% CIs are normal-based except MR-Egger (t, k - 2 df), as in
#   Figure 4 and Supplementary Table S7.
#   Up to release v13 this step re-estimated all three outcomes on one "standardised beta" axis. For glaucoma
#   that was an effect on the observed 0/1 scale, which depends on the case fraction of FinnGen and is not
#   comparable with the continuous traits, so the common axis was dropped.
#   Outputs: analysis/nr1h3_cismr_forest_v3.{pdf,png} (Figure 6), analysis/nr1h3_cismr_fig6_values.tsv
# ============================================================
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork) })
BASE <- Sys.getenv("FECD_BASE", path.expand("~/FECD_GWAS")); ANA <- file.path(BASE, "analysis")
source(file.path(Sys.getenv("SCRIPTS_DIR", file.path(BASE, "scripts")), "fig_common.R"))

t3 <- fread(file.path(ANA, "table3_nr1h3_crosstrait_cismr.tsv"))
lev <- c("Inverse variance weighted", "MR Egger", "Weighted median", "Weighted mode")
lab <- c("Inverse variance weighted" = "IVW (primary)", "MR Egger" = "MR-Egger",
         "Weighted median" = "Weighted median", "Weighted mode" = "Weighted mode")
t3[, crit := ifelse(method == "MR Egger", qt(0.975, nsnp - 2), qnorm(0.975))]
t3[, `:=`(lo = b - crit * se, hi = b + crit * se)]

outc <- list(
  Glaucoma = list(title = "Glaucoma (FinnGen R12 H7_GLAUCOMA; 26,591 cases, 473,757 controls)",
                  x = "Odds ratio for glaucoma (log scale)", or = TRUE, unit = "OR"),
  IOP      = list(title = "Intraocular pressure (Bonnemaijer et al.; n = 31,269)",
                  x = "Difference in IOP (mmHg)", or = FALSE, unit = "mmHg"),
  Myopia   = list(title = "Refractive error (Hysi et al.; n = 351,091)",
                  x = "Difference in refractive error (SD; positive = more hyperopic)", or = FALSE, unit = "SD"))

vals <- list(); rows <- list()
for (nm in names(outc)) {
  o <- outc[[nm]]
  d <- t3[outcome == nm][match(lev, method)]
  stopifnot(nrow(d) == 4, !anyNA(d$b))
  if (o$or) d[, `:=`(est = exp(b), elo = exp(lo), ehi = exp(hi))] else d[, `:=`(est = b, elo = lo, ehi = hi)]
  dig <- if (nm == "Myopia") 3 else 2
  d[, txt := sprintf("%s%s (%s to %s); *P* = %s", if (o$or) "OR " else "",
                     fmt_num(est, dig), fmt_num(elo, dig), fmt_num(ehi, dig), fmt_p(pval))]
  d[, ylab := factor(lab[method], levels = rev(lab[lev]))]
  vals[[nm]] <- d[, .(outcome = nm, scale = o$unit, method, nsnp, b, se, pval, estimate = est, ci_low = elo, ci_high = ehi)]

  pf <- ggplot(d, aes(est, ylab)) +
    geom_vline(xintercept = if (o$or) 1 else 0, linetype = "dashed", colour = "grey45") +
    geom_errorbarh(aes(xmin = elo, xmax = ehi), height = 0.2, linewidth = 0.6, colour = COL_EST) +
    geom_point(aes(size = method == "Inverse variance weighted"), colour = COL_EST) +
    scale_size_manual(values = c(`FALSE` = 2.2, `TRUE` = 3.6), guide = "none") +
    labs(x = o$x, y = NULL, title = o$title) +
    theme_paper(10) +
    theme(panel.grid.major.y = element_blank(), plot.title = element_markdown(size = rel(1), face = "bold"))
  pf <- pf + if (o$or) scale_x_log10(breaks = c(0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2), labels = lab_minus)
              else scale_x_continuous(labels = lab_minus)
  pt <- ggplot(d, aes(0, ylab)) +
    geom_richtext(aes(label = txt), hjust = 0, size = 3.3, family = FIG_FONT, fill = NA, label.colour = NA,
                  label.padding = unit(0, "pt"), colour = "grey15") +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    labs(title = " ") + theme_void(base_family = FIG_FONT) +
    theme(plot.title = element_text(size = 10))
  rows[[nm]] <- pf + pt + plot_layout(widths = c(1.35, 1))
}
fwrite(rbindlist(vals), file.path(ANA, "nr1h3_cismr_fig6_values.tsv"), sep = "\t", quote = FALSE)

F_SET <- 597   # set-level F of the four instruments (steps 89c/91; independent of the allele frequency used)
p <- wrap_plots(rows, ncol = 1) +
  plot_annotation(
    title = paste0("Cross-trait cis-Mendelian randomisation of ", gi("NR1H3"), " expression"),
    subtitle = paste0("Four independent eQTLGen whole-blood cis-eQTL instruments (*P* < 5 × 10<sup>",
                      MINUS, "8</sup>; r<sup>2</sup> < 0.01 within 5,000 kb); set-level F = ", F_SET, ".<br>",
                      "Effects per SD higher genetically predicted expression; 95% confidence intervals<br>",
                      "(MR-Egger: t distribution, k ", MINUS, " 2 df)."),
    theme = theme(plot.title = element_markdown(face = "bold", size = 11, family = FIG_FONT),
                  plot.subtitle = element_markdown(size = 9, colour = "grey25", family = FIG_FONT, lineheight = 1.25)))
save_fig(p, file.path(ANA, "nr1h3_cismr_forest_v3"), width = FIG_W, height = 7.2, fig = "Figure6")
print(rbindlist(vals)[, .(outcome, method, estimate = signif(estimate, 3), ci_low = signif(ci_low, 3),
                         ci_high = signif(ci_high, 3), pval = signif(pval, 3))])
cat("Saved: analysis/nr1h3_cismr_forest_v3.{pdf,png} (Figure 6)\n")
