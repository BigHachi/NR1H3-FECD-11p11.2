#!/usr/bin/env Rscript
# ============================================================
# Step 108: Full-cis colocalisation at 11p11.2 for the glaucoma arm
#   (review point 2 / 3). Replaces script 79, which used significant-only
#   GTEx Portal eQTLs (fibroblast PP4 0.68, whole blood 0.61 in v2.7).
#   A. NR1H3 eQTL (eQTL Catalogue r6 = GTEx v8, ALL pairs) x glaucoma
#      (FinnGen R12 H7_GLAUCOMA, 26,591 cases / 473,757 controls, GRCh38)
#      in fibroblast and whole blood.
#   B. FECD GWAS x glaucoma GWAS at 11p11.2 (NR1H3 TSS +/- 1 Mb): do the two
#      diseases share a causal variant? FECD = EUR (primary) and multi-ancestry.
#   Matching: chr + pos + allele set. Priors p1=p2=1e-4, p12=1e-5 (+ p12=5e-6).
# ============================================================
suppressPackageStartupMessages({ library(data.table); library(coloc) })
BASE <- Sys.getenv("FECD_BASE", path.expand("~/FECD_GWAS"))
CHR <- "11"; TSS <- 47248300; S <- TSS - 1.05e6; E <- TSS + 1.05e6

gl <- fread(file.path(BASE, "data/gwas/glaucoma_gwas_full.tsv"), colClasses = list(character = "CHR"))
gl <- gl[CHR == "11" & BP >= S & BP <= E & is.finite(BETA) & SE > 0]
GL <- list(N = 26591 + 473757, s = 26591 / (26591 + 473757))

win <- sprintf("awk -F'\\t' 'NR==1 || ($%d==\"11\" && $%d>=%d && $%d<=%d)'", 1, 2, as.integer(S), 2, as.integer(E))
eur <- fread(cmd = paste("gzcat", file.path(BASE, "data/gwas/fecd_eur_GCST90475865.tsv.gz"), "|", win),
             select = c("chromosome","base_pair_location","effect_allele","other_allele","odds_ratio","ci_upper","ci_lower","p_value"),
             colClasses = list(character = "chromosome"))
setnames(eur, c("CHR","BP","A1","A2","OR","CIU","CIL","P"))
eur[, `:=`(BETA = log(as.numeric(OR)), SE = (log(as.numeric(CIU)) - log(as.numeric(CIL))) / (2*qnorm(0.975)))]
eur <- eur[is.finite(BETA) & is.finite(SE) & SE > 0]
win2 <- sprintf("awk -F'\\t' 'NR==1 || ($%d==\"11\" && $%d>=%d && $%d<=%d)'", 2, 3, as.integer(S), 3, as.integer(E))
mul <- fread(cmd = paste("gzcat", file.path(BASE, "data/gwas/fecd_gwas_full.tsv.gz"), "|", win2),
             select = c("CHR","BP","A1","A2","BETA","SE","P"), colClasses = list(character = "CHR"))
mul <- mul[is.finite(BETA) & SE > 0]
FE <- list(EUR = list(d = eur, N = 448712, s = 3473/448712),
           multiancestry = list(d = mul, N = 570321, s = 3928/570321))

key <- function(d, bp, a, b) d[, id := paste(get(bp), pmin(get(a), get(b)), pmax(get(a), get(b)), sep = ":")][!duplicated(id)]
gl <- key(copy(gl), "BP", "A1", "A2")

summ <- function(lab1, lab2, r, n, extra = list()) {
  p <- r$summary; r5 <- extra$r5
  data.table(trait1 = lab1, trait2 = lab2, n_snps = n,
             PP0 = p[["PP.H0.abf"]], PP1 = p[["PP.H1.abf"]], PP2 = p[["PP.H2.abf"]],
             PP3 = p[["PP.H3.abf"]], PP4 = p[["PP.H4.abf"]],
             PP4_p12_5e6 = if (!is.null(r5)) r5[["PP.H4.abf"]] else NA_real_,
             top_snp = r$results$snp[which.max(r$results$SNP.PP.H4)],
             top_snp_PP4 = max(r$results$SNP.PP.H4))
}
out <- list()
# ---- A: NR1H3 eQTL x glaucoma ----
for (t in c("fibroblast", "blood")) {
  e <- fread(file.path(BASE, sprintf("data/eqtl_catalogue/ENSG00000025434_%s_v8_allpairs.tsv", t)),
             colClasses = list(character = c("chromosome","ref","alt")))
  e <- e[is.finite(beta) & se > 0 & maf > 0 & maf < 1]
  e <- key(e, "position", "ref", "alt")
  m <- merge(e[, .(id, b_e = beta, se_e = se, maf, N = an/2)], gl[, .(id, BETA, SE)], by = "id")
  d1 <- list(beta = m$b_e, varbeta = m$se_e^2, snp = m$id, type = "quant", N = round(median(m$N)), MAF = m$maf)
  d2 <- list(beta = m$BETA, varbeta = m$SE^2, snp = m$id, type = "cc", s = GL$s, N = GL$N)
  r  <- suppressMessages(coloc.abf(d1, d2, p12 = 1e-5)); r5 <- suppressMessages(coloc.abf(d1, d2, p12 = 5e-6))$summary
  out[[length(out)+1]] <- summ(paste0("NR1H3 eQTL (", t, ")"), "Glaucoma (FinnGen R12)", r, nrow(m), list(r5 = r5))
}
# ---- B: FECD x glaucoma ----
for (anc in names(FE)) {
  f <- key(copy(FE[[anc]]$d), "BP", "A1", "A2")
  m <- merge(f[, .(id, b_f = BETA, se_f = SE)], gl[, .(id, BETA, SE)], by = "id")
  d1 <- list(beta = m$b_f, varbeta = m$se_f^2, snp = m$id, type = "cc", s = FE[[anc]]$s, N = FE[[anc]]$N)
  d2 <- list(beta = m$BETA, varbeta = m$SE^2, snp = m$id, type = "cc", s = GL$s, N = GL$N)
  r  <- suppressMessages(coloc.abf(d1, d2, p12 = 1e-5)); r5 <- suppressMessages(coloc.abf(d1, d2, p12 = 5e-6))$summary
  out[[length(out)+1]] <- summ(paste0("FECD (", anc, ")"), "Glaucoma (FinnGen R12)", r, nrow(m), list(r5 = r5))
}
res <- rbindlist(out)
cat(sprintf("glaucoma window min P = %.2e | FECD EUR min P = %.2e | FECD multi min P = %.2e\n",
            min(gl$P), min(eur$P), min(mul$P)))
print(res[, lapply(.SD, function(x) if (is.numeric(x)) signif(x, 3) else x)])
fwrite(res, file.path(BASE, "analysis/coloc_11p11_glaucoma_fullcis.tsv"), sep = "\t")
