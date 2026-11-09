#!/usr/bin/env Rscript
# ============================================================
# Step 79b: NR1H3 eQTL x glaucoma colocalisation on FULL cis summary statistics.
#   Replaces step 79, which used the GTEx Portal API (significant variant-gene pairs only).
#   eQTL   : eQTL Catalogue r6 GTEx v8, all variant-gene pairs for NR1H3 (step 100):
#            fibroblasts, whole blood (the two tissues in the Results), skin, oesophagus mucosa.
#   GWAS   : FinnGen R12 H7_GLAUCOMA (26,591 cases / 500,348), chr11 region (step 78 file).
#   Match  : position + alleles (either orientation). coloc.abf, p1 = p2 = 1e-4, p12 = 1e-5
#            (PP4 at p12 = 5e-6 also reported).
#   Output : results/mr_nr1h3_glaucoma/nr1h3_glaucoma_coloc_fullcis.tsv
# ============================================================
suppressPackageStartupMessages({ library(data.table); library(coloc) })
BASE <- Sys.getenv("FECD_BASE", path.expand("~/FECD_GWAS"))
OUT  <- file.path(BASE, "results/mr_nr1h3_glaucoma"); dir.create(OUT, FALSE, TRUE)
GLA_N <- 500348L; GLA_S <- 26591 / 500348
TISSUES <- c(fibroblast = "Fibroblast", blood = "Whole blood", skin_not_sun_exposed = "Skin",
             esophagus_mucosa = "Oesophagus mucosa")

reg <- file.path(BASE, "data/gwas/glaucoma_nr1h3_region.tsv")
if (!file.exists(reg))
  system(sprintf("awk -F'\t' 'NR==1||($2==11&&$3>=46000000&&$3<=48500000)' %s > %s",
                 file.path(BASE, "data/gwas/glaucoma_gwas_full.tsv"), reg))
gw <- fread(reg)[is.finite(BETA) & is.finite(SE) & SE > 0]

res <- rbindlist(lapply(names(TISSUES), function(t) {
  f <- file.path(BASE, sprintf("data/eqtl_catalogue/ENSG00000025434_%s_v8_allpairs.tsv", t))
  if (!file.exists(f)) return(data.table(tissue = TISSUES[[t]], status = "eQTL file missing (run step 100)"))
  e <- fread(f, colClasses = list(character = c("chromosome", "ref", "alt")))
  e <- e[is.finite(beta) & is.finite(se) & se > 0 & maf > 0 & maf < 1 & !duplicated(variant)]
  m <- merge(e[, .(BP = position, ref, alt, b_e = beta, se_e = se, p_e = pvalue, maf, N = an / 2, rsid)],
             gw[, .(BP, A1, A2, BETA, SE, P)], by = "BP", allow.cartesian = TRUE)
  m <- m[(A1 == alt & A2 == ref) | (A1 == ref & A2 == alt)]
  m[, id := paste(BP, ref, alt, sep = ":")]; m <- m[!duplicated(id)]
  d1 <- list(beta = m$BETA, varbeta = m$SE^2, snp = m$id, type = "cc", N = GLA_N, s = GLA_S)
  d2 <- list(beta = m$b_e, varbeta = m$se_e^2, snp = m$id, type = "quant", N = round(median(m$N)), MAF = m$maf)
  r  <- suppressMessages(coloc.abf(d1, d2, p1 = 1e-4, p2 = 1e-4, p12 = 1e-5)); pp <- r$summary
  r5 <- suppressMessages(coloc.abf(d1, d2, p1 = 1e-4, p2 = 1e-4, p12 = 5e-6))$summary
  top <- as.data.table(r$results)[which.max(SNP.PP.H4)]
  data.table(tissue = TISSUES[[t]], n_snps = nrow(m), eqtl_N = round(median(m$N)),
             min_p_glaucoma = min(m$P), min_p_eqtl = min(m$p_e),
             PP0 = pp[["PP.H0.abf"]], PP1 = pp[["PP.H1.abf"]], PP2 = pp[["PP.H2.abf"]],
             PP3 = pp[["PP.H3.abf"]], PP4 = pp[["PP.H4.abf"]], PP4_p12_5e6 = r5[["PP.H4.abf"]],
             top_snp = m[id == top$snp, rsid][1], top_snp_PP4 = top$SNP.PP.H4, status = "OK")
}), fill = TRUE)
fwrite(res, file.path(OUT, "nr1h3_glaucoma_coloc_fullcis.tsv"), sep = "\t")
print(res[, .(tissue, n_snps, min_p_glaucoma = signif(min_p_glaucoma, 2), min_p_eqtl = signif(min_p_eqtl, 2),
              PP3 = round(PP3, 3), PP4 = round(PP4, 3), PP4_p12_5e6 = round(PP4_p12_5e6, 3), top_snp)])
