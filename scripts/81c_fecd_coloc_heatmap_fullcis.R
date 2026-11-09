#!/usr/bin/env Rscript
# ============================================================
# Step 81c (Oct 13 2026): Figure 3 + Supplementary Table S6 from FULL cis eQTL summary statistics.
#   Replaces step 81b (GTEx Portal significant variant-gene pairs only, which coloc.abf is not designed for).
#   Source : step 102c (analysis/fecd_coloc_allpairs_v2.tsv): coloc.abf (p1 = p2 = 1e-4, p12 = 1e-5; PP4 at
#            p12 = 5e-6 also reported), eQTL Catalogue r6 GTEx v8 all variant-gene pairs in the gene's cis
#            window (steps 101/101c), FECD European GWAS (primary) and multi-ancestry GWAS (sensitivity).
#   Genes  : data/tophits21/genes_coloc.tsv = European S-MultiXcan ranks 1-23 (down to NR1H3) + four genes from
#            earlier versions of S6; S_MultiXcan_rank = European rank written by step 19f (NA = no European P).
#   Output : analysis/fecd_coloc_S6_fullcis.tsv, analysis/fecd_coloc_heatmap_fullcis.{pdf,png}
# ============================================================
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
BASE <- Sys.getenv("FECD_BASE", path.expand("~/FECD_GWAS"))
ANA  <- file.path(BASE, "analysis")

SYM_OVR <- fread(file.path(BASE, "data/gene_symbol_display.tsv"))
display_symbol <- function(ensg, sym) {
  i <- match(sub("[.].*", "", ensg), SYM_OVR$ENSG)
  ifelse(is.na(i), sym, SYM_OVR$display_symbol[i])
}
tissue_labels <- c(fibroblast = "Fibroblasts", skin_not_sun_exposed = "Skin\n(Not Sun Exposed)",
                   esophagus_mucosa = "Oesophagus\nMucosa")

co <- fread(file.path(ANA, "fecd_coloc_allpairs_v2.tsv"))
co[, SYMBOL := display_symbol(ensg, gene)]
co[, ancestry := fifelse(ancestry == "EUR", "European", "Multi-ancestry")]
stopifnot(!anyDuplicated(co[, .(ensg, tissue, ancestry)]))

# ---- S6: every gene x tissue x GWAS, all posterior probabilities ----
s6 <- co[, .(SYMBOL, ENSG = sub("[.].*", "", ensg), S_MultiXcan_rank = rank, tissue = gsub("\n", " ", tissue_labels[tissue]),
             GWAS = ancestry, n_variants = n_snps, min_P_GWAS = min_p_gwas, min_P_eQTL = min_p_eqtl,
             PP0, PP1, PP2, PP3, PP4, PP4_p12_5e6, PP4_over_PP3_plus_PP4 = PP4_given_H3H4,
             top_variant = top_snp, top_variant_PP4 = top_snp_PP4, status)]
setorder(s6, GWAS, -PP4, na.last = TRUE)
fwrite(s6, file.path(ANA, "fecd_coloc_S6_fullcis.tsv"), sep = "\t", quote = FALSE, na = "NA")

# ---- Figure 3: PP4, European GWAS ----
pd <- co[ancestry == "European", .(SYMBOL, tissue, PP4)]
pd[, tissue_label := factor(tissue_labels[tissue], levels = tissue_labels)]
gene_order <- pd[, .(m = max(c(PP4, -1), na.rm = TRUE)), by = SYMBOL][order(m), SYMBOL]
pd[, SYMBOL := factor(SYMBOL, levels = gene_order)]
# 10 Nov 2026 (review): shared figure style (fig_common.R): italic gene symbols, vector PDF + TIFF (Figure 3)
source(file.path(Sys.getenv("SCRIPTS_DIR", file.path(BASE, "scripts")), "fig_common.R"))
p <- ggplot(pd, aes(x = tissue_label, y = SYMBOL, fill = PP4)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = ifelse(is.na(PP4), "\u2014", sprintf("%.2f", PP4)),
                fontface = ifelse(!is.na(PP4) & PP4 >= 0.8, "bold", "plain")), size = 3.1, family = FIG_FONT,
            color = ifelse(is.na(pd$PP4), "grey55", "black")) +
  scale_fill_gradient2(low = "white", mid = "steelblue", high = "firebrick", midpoint = 0.5,
                       limits = c(0, 1), na.value = "grey92", name = "PP4") +
  scale_y_discrete(labels = function(g) gi(as.character(g))) +
  labs(x = NULL, y = NULL, title = "Colocalisation of FECD TWAS genes with GTEx v8 cis-eQTLs",
       subtitle = paste0("Posterior probability of a shared causal variant (PP4; coloc.abf), European-ancestry FECD GWAS<br>",
                         "and full cis-eQTL summary statistics (eQTL Catalogue); bold, PP4 \u2265 0.8; grey (\u2014), no eQTL data")) +
  theme_paper(10) + theme(axis.text.x = element_text(size = 10, colour = "grey15"), panel.grid = element_blank())
h <- max(4.5, 1.4 + 0.28 * length(gene_order))   # 27 genes: 9.0 in, within the journal page height
save_fig(p, file.path(ANA, "fecd_coloc_heatmap_fullcis"), width = 7, height = h, fig = "Figure3")

eu <- co[ancestry == "European"]
cat(sprintf("European GWAS: %d genes x 3 tissues; pairs with PP4 >= 0.8: %s\n", uniqueN(eu$SYMBOL),
            paste(eu[PP4 >= 0.8, sprintf("%s (%s %.2f)", SYMBOL, tissue, PP4)], collapse = "; ")))
print(dcast(eu, SYMBOL ~ tissue, value.var = "PP4")[order(SYMBOL)], digits = 2)
