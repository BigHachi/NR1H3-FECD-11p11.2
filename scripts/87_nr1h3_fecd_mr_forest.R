#!/usr/bin/env Rscript
# ============================================================
# Step 87: Forest plot — NR1H3 (eQTLGen blood) -> FECD cis-MR
#   Four estimators (IVW primary + 3 sensitivity), OR scale.
#   Usage: Rscript 87_nr1h3_fecd_mr_forest.R [multiancestry|EUR]
#   Ancestry is explicit in the filename, title and subtitle so the
#   figure can never be mistaken for the other cohort.
# ============================================================
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
BASE <- Sys.getenv("FECD_BASE", path.expand("~/FECD_GWAS"))

ANC <- commandArgs(trailingOnly=TRUE)[1]
if (is.na(ANC)) ANC <- "multiancestry"
if (!ANC %in% c("multiancestry","EUR"))
  stop("ancestry must be 'multiancestry' or 'EUR', got: ", ANC)

d  <- fread(file.path(BASE, sprintf("analysis/nr1h3_eqtlgen_fecd_mr_%s.tsv", ANC)))
dg <- fread(file.path(BASE, sprintf("analysis/nr1h3_eqtlgen_fecd_mr_%s_diagnostics.tsv", ANC)))
anc_label <- if (ANC=="EUR") "European subset (GCST90475865)" else "multi-ancestry (GCST90480062)"
# Method-aware 95% CI multiplier: MR-Egger P is t(n-2)-based (intercept+slope),
# so its CI must use the t quantile, NOT the normal 1.96 (else CI<->P disagree).
d[, crit := ifelse(method=="MR Egger", qt(0.975, nsnp-2), qnorm(0.975))]
d[, OR    := exp(b)]
d[, OR_lo := exp(b - crit*se)]
d[, OR_hi := exp(b + crit*se)]

# order top->bottom: IVW (primary) first
lev <- c("Weighted mode","Weighted median","MR Egger","Inverse variance weighted")
lab <- c("Inverse variance weighted"="IVW (primary)", "MR Egger"="MR-Egger",
         "Weighted median"="Weighted median", "Weighted mode"="Weighted mode")
d[, method := factor(method, levels=lev)]
d[, ylab := lab[as.character(method)]]
d[, ylab := factor(ylab, levels=lab[lev])]
d[, txt  := sprintf("OR %.2f (%.2f-%.2f),  P = %s",
                    OR, OR_lo, OR_hi, ifelse(pval<0.001, formatC(pval,format="e",digits=1),
                                             sprintf("%.3f", pval)))]

# 10 Nov 2026 (review): same layout and style as Figure 6 (forest + estimate column, italic gene symbol,
# P-values as a x 10^b, vector PDF), via fig_common.R
source(file.path(Sys.getenv("SCRIPTS_DIR", file.path(BASE, "scripts")), "fig_common.R"))
suppressPackageStartupMessages(library(patchwork))
d[, txt := sprintf("OR %s (%s to %s); *P* = %s", fmt_num(OR), fmt_num(OR_lo), fmt_num(OR_hi), fmt_p(pval))]
pf <- ggplot(d, aes(x=OR, y=ylab)) +
  geom_vline(xintercept=1, linetype="dashed", colour="grey45") +
  geom_errorbarh(aes(xmin=OR_lo, xmax=OR_hi), height=0.2, linewidth=0.6, colour=COL_EST) +
  geom_point(aes(size=(method=="Inverse variance weighted")), colour=COL_EST) +
  scale_size_manual(values=c(`FALSE`=2.2, `TRUE`=3.6), guide="none") +
  scale_x_log10(breaks=c(0.6,0.8,1,1.25,1.5,2,2.5), labels=lab_minus) +
  labs(x=paste0("Odds ratio for FECD per SD higher genetically predicted ", gi("NR1H3"), " expression (log scale)"),
       y=NULL) +
  theme_paper(10) + theme(panel.grid.major.y=element_blank())
pt <- ggplot(d, aes(0, ylab)) +
  geom_richtext(aes(label=txt), hjust=0, size=3.3, family=FIG_FONT, fill=NA, label.colour=NA,
                label.padding=unit(0, "pt"), colour="grey15") +
  scale_x_continuous(limits=c(0, 1), expand=c(0, 0)) + theme_void(base_family=FIG_FONT)
anc_title <- if (ANC=="EUR") "European-ancestry FECD GWAS (GCST90475865; 3,473 cases, 445,239 controls)" else
  "Multi-ancestry FECD GWAS (GCST90480062; 3,928 cases, 566,393 controls)"
p <- pf + pt + plot_layout(widths=c(1.35, 1)) +
  plot_annotation(
    title=paste0("cis-Mendelian randomisation of genetically predicted ", gi("NR1H3"), " expression on FECD"),
    subtitle=sprintf(paste0("%s<br>Four independent eQTLGen whole-blood cis-eQTL instruments; set-level F = %.0f<br>",
                            "Cochran's Q *P* = %s; MR-Egger intercept *P* = %s"),
                     anc_title, dg$F_set, fmt_p(dg$Q_pval), fmt_p(dg$egger_int_p)),
    theme=theme(plot.title=element_markdown(face="bold", size=11, family=FIG_FONT),
                plot.subtitle=element_markdown(size=9, colour="grey25", family=FIG_FONT, lineheight=1.25)))

stem <- sprintf("analysis/nr1h3_fecd_mr_forest_%s", ANC)
save_fig(p, file.path(BASE, stem), width=FIG_W, height=3.5, fig=if (ANC=="EUR") "Figure4" else NULL, png_dpi=200)
cat(sprintf("Wrote %s.{pdf,png}\n", stem))
print(d[, .(method, OR=round(OR,3), OR_lo=round(OR_lo,3), OR_hi=round(OR_hi,3), P=signif(pval,3))])
