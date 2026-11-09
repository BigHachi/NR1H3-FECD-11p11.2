#!/usr/bin/env Rscript
# ============================================================
# Step 76c: Figure 5 - cross-trait NR1H3 FUSION TWAS forest
#   Single framework (FUSION TWAS), tissue varies by trait:
#     FECD      -> disease-proxy tissues (fibroblast, skin; GTEx EUR weights)
#     glaucoma/IOP/ -> whole blood (GTEx EUR weights)
#     Myopia
#   Cross-trait rows also carry FUSION conditional/joint status:
#     conditionally independent (survives joint analysis) vs conditioned away.
#   Oct 21 2026 (review 13): signed-Z point plot. A TWAS Z is a test statistic, not an effect per SD of
#   expression, so the axis is "FUSION TWAS Z-score" and the old Z +/- 1.96 bars (no inferential target) are gone.
#   "Refractive error" replaces "myopia" (Hysi sign convention: positive = more hyperopic refraction).
#   Values READ from FUSION outputs (fully reproducible; steps 60/76*/post-process):
#     analysis/nr1h3_crosstrait_fusion_conditional.tsv
# ============================================================
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
BASE <- Sys.getenv("FECD_BASE", path.expand("~/FECD_GWAS")); ANA <- file.path(BASE,"analysis")
d <- fread(file.path(ANA,"nr1h3_crosstrait_fusion_conditional.tsv"))

lab <- c(FECD_fibroblast="FECD \u2013 cultured fibroblasts",
         FECD_skin      ="FECD \u2013 skin (not sun exposed)",
         glaucoma       ="Glaucoma \u2013 whole blood",
         myopia         ="Refractive error \u2013 whole blood",
         iop            ="IOP \u2013 whole blood")
grp <- c(FECD_fibroblast="FECD", FECD_skin="FECD", glaucoma="Glaucoma",
         myopia="Myopia", iop="IOP")
d <- d[trait %in% names(lab)]
d[, label := factor(lab[trait], levels=rev(lab[c("FECD_fibroblast","FECD_skin","glaucoma","myopia","iop")]))]
d[, Z := as.numeric(TWAS_Z)][, P := as.numeric(TWAS_P)]
# 10 Nov 2026 (review): shared figure style (fig_common.R): italic gene symbol, P as a x 10^b, vector PDF + TIFF
source(file.path(Sys.getenv("SCRIPTS_DIR", file.path(BASE, "scripts")), "fig_common.R"))
d[, txt := sprintf("Z = %s; *P* = %s", sub("^\\+", "+", fmt_num(Z)), fmt_p(P))]
d[Z > 0, txt := paste0("+", txt)]
d[, txt := sub("^\\+Z = ", "Z = +", txt)]
p <- ggplot(d, aes(Z, label)) +
  geom_vline(xintercept=0, linetype=2, colour="grey55") +
  geom_point(size=3.2, colour=COL_EST) +
  geom_richtext(aes(label=txt), vjust=-0.6, size=3.2, family=FIG_FONT, fill=NA, label.colour=NA, colour="grey15") +
  scale_x_continuous(limits=c(-9.5, 8), labels=lab_minus) +
  scale_y_discrete(expand=expansion(add=c(0.5, 0.9))) +
  labs(x=paste0("FUSION TWAS Z-score for genetically predicted ", gi("NR1H3"), " expression"), y=NULL,
       title=paste0("Cross-trait TWAS of ", gi("NR1H3"), " (FUSION, GTEx v8 EUR weights)"),
       subtitle=paste0("FECD in disease-proxy tissues; glaucoma (FinnGen R12, all subtypes), refractive error and IOP in whole blood<br>",
                       "Negative Z: higher predicted expression with lower FECD or glaucoma risk or lower IOP<br>",
                       "Positive Z for refractive error: higher predicted expression with more hyperopic refraction")) +
  theme_paper(10) + theme(panel.grid.major.y=element_blank())
save_fig(p, file.path(ANA,"nr1h3_cross_trait_twas_forest_v2"), width=FIG_W, height=4.1, fig="Figure5", png_dpi=200)
cat("Wrote analysis/nr1h3_cross_trait_twas_forest_v2.{pdf,png}\n")
print(d[, .(label, Z, P=signif(P,2))])
