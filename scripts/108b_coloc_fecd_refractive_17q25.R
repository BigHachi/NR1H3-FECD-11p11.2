#!/usr/bin/env Rscript
# ============================================================
# Step 108b (13 Nov 2026, review; release v15.0): FECD x refractive error at 17q25.3 (TSPAN10 TSS +/- 1.05 Mb, GRCh38), the same
#   coloc.abf set-up as step 108 part B (FECD x glaucoma at 11p11.2). Question: does the TSPAN10 FECD signal share a
#   causal variant with refractive error (corneal shape / refraction, or ascertainment), and which way does the FECD risk
#   allele shift refraction?
#   FECD: European GWAS (GCST90475865; primary) and multi-ancestry GWAS. Refractive error: Hysi et al. (publicly released
#   subset; Z-scores with per-variant N and allele frequency; positive Z = more hyperopic refraction), step 55/55b file
#   data/gwas/myopia_gwas_full.tsv (BETA = Z, SE = 1). Standardised effects for coloc (quantitative trait, sdY = 1):
#   beta = Z / sqrt(2p(1-p)(N + Z^2)), varbeta = 1 / (2p(1-p)(N + Z^2)).
#   Matching: chr + pos + allele set. Priors p1 = p2 = 1e-4, p12 = 1e-5 (+ p12 = 5e-6).
#   Outputs: analysis/coloc_17q25_fecd_refractive.tsv (coloc summaries)
#            analysis/coloc_17q25_fecd_refractive_direction.tsv (FECD lead and colocalising variants: FECD risk allele,
#            its refractive-error Z and direction)
# ============================================================
suppressPackageStartupMessages({ library(data.table); library(coloc) })
BASE <- Sys.getenv("FECD_BASE", path.expand("~/FECD_GWAS"))
CHR <- "17"; TSS <- 81637171; S <- TSS - 1.05e6; E <- TSS + 1.05e6

re <- fread(file.path(BASE, "data/gwas/myopia_gwas_full.tsv"), colClasses = list(character = "CHR"),
            select = c("SNP", "CHR", "BP", "A1", "A2", "BETA", "P", "FREQ", "N"))
re <- re[CHR == CHR & BP >= S & BP <= E & is.finite(BETA) & is.finite(FREQ) & FREQ > 0 & FREQ < 1 & N > 0]
setnames(re, "BETA", "Z")
re[, maf := pmin(FREQ, 1 - FREQ)]
re[, d := 2 * FREQ * (1 - FREQ) * (N + Z^2)]
re[, `:=`(b_r = Z / sqrt(d), v_r = 1 / d)]

win <- sprintf("awk -F'\\t' 'NR==1 || ($%d==\"%s\" && $%d>=%d && $%d<=%d)'", 1, CHR, 2, as.integer(S), 2, as.integer(E))
eur <- fread(cmd = paste("gzcat", file.path(BASE, "data/gwas/fecd_eur_GCST90475865.tsv.gz"), "|", win),
             select = c("chromosome", "base_pair_location", "effect_allele", "other_allele", "odds_ratio", "ci_upper",
                        "ci_lower", "p_value", "rsid"), colClasses = list(character = "chromosome"))
setnames(eur, c("CHR", "BP", "A1", "A2", "OR", "CIU", "CIL", "P", "SNP"))
eur[, `:=`(BETA = log(as.numeric(OR)), SE = (log(as.numeric(CIU)) - log(as.numeric(CIL))) / (2 * qnorm(0.975)))]
eur <- eur[is.finite(BETA) & is.finite(SE) & SE > 0]
win2 <- sprintf("awk -F'\\t' 'NR==1 || ($%d==\"%s\" && $%d>=%d && $%d<=%d)'", 2, CHR, 3, as.integer(S), 3, as.integer(E))
mul <- fread(cmd = paste("gzcat", file.path(BASE, "data/gwas/fecd_gwas_full.tsv.gz"), "|", win2),
             select = c("SNP", "CHR", "BP", "A1", "A2", "BETA", "SE", "P"), colClasses = list(character = "CHR"))
mul <- mul[is.finite(BETA) & SE > 0]
FE <- list(EUR = list(d = eur, N = 448712, s = 3473 / 448712),
           multiancestry = list(d = mul, N = 570321, s = 3928 / 570321))

key <- function(d, bp, a, b) d[, id := paste(get(bp), pmin(toupper(get(a)), toupper(get(b))),
                                              pmax(toupper(get(a)), toupper(get(b))), sep = ":")][!duplicated(id)]
re <- key(copy(re), "BP", "A1", "A2")

out <- list(); dirs <- list()
for (anc in names(FE)) {
  f <- key(copy(FE[[anc]]$d), "BP", "A1", "A2")
  m <- merge(f[, .(id, SNP_f = SNP, A1_f = toupper(A1), b_f = BETA, se_f = SE, P_f = P)],
             re[, .(id, SNP, A1_r = toupper(A1), b_r, v_r, Z, P_r = P, maf, N)], by = "id")
  m[, b_r_aligned := ifelse(A1_r == A1_f, b_r, -b_r)]          # refractive-error effect of the FECD effect allele
  m[, Z_aligned := ifelse(A1_r == A1_f, Z, -Z)]
  d1 <- list(beta = m$b_f, varbeta = m$se_f^2, snp = m$id, type = "cc", s = FE[[anc]]$s, N = FE[[anc]]$N)
  d2 <- list(beta = m$b_r_aligned, varbeta = m$v_r, snp = m$id, type = "quant", sdY = 1, N = round(median(m$N)))
  r  <- suppressMessages(coloc.abf(d1, d2, p12 = 1e-5)); r5 <- suppressMessages(coloc.abf(d1, d2, p12 = 5e-6))$summary
  p  <- r$summary; top <- r$results$snp[which.max(r$results$SNP.PP.H4)]
  out[[anc]] <- data.table(trait1 = paste0("FECD (", anc, ")"), trait2 = "Refractive error (Hysi et al.)",
                           locus = "17q25.3 (TSPAN10 TSS +/- 1.05 Mb)", n_snps = nrow(m),
                           min_P_FECD = min(m$P_f), min_P_refractive = min(m$P_r),
                           PP0 = p[["PP.H0.abf"]], PP1 = p[["PP.H1.abf"]], PP2 = p[["PP.H2.abf"]],
                           PP3 = p[["PP.H3.abf"]], PP4 = p[["PP.H4.abf"]], PP4_p12_5e6 = r5[["PP.H4.abf"]],
                           top_snp = m$SNP[m$id == top], top_snp_PP4 = max(r$results$SNP.PP.H4))
  # direction: FECD lead variant, the eQTL-colocalising variants of S6, the coloc top variant
  pick <- unique(c(m$SNP[which.min(m$P_f)], "rs7405790", "rs7405453", m$SNP[m$id == top]))
  for (s in pick) {
    x <- m[SNP == s]
    if (!nrow(x)) next
    risk_sign <- sign(x$b_f)                                    # FECD risk allele = effect allele if b_f > 0
    zr <- x$Z_aligned * risk_sign                               # refractive-error Z of the FECD RISK allele
    dirs[[length(dirs) + 1]] <- data.table(
      gwas = anc, variant = s, role = c(if (s == m$SNP[which.min(m$P_f)]) "FECD lead in window",
                                         if (s %in% c("rs7405790", "rs7405453")) "TSPAN10 eQTL colocalising variant (S6)",
                                         if (s == m$SNP[m$id == top]) "FECD x refractive-error coloc top variant")[1],
      FECD_risk_allele = ifelse(risk_sign > 0, x$A1_f, setdiff(strsplit(x$id, ":")[[1]][2:3], x$A1_f)[1]),
      FECD_OR_risk = exp(abs(x$b_f)), FECD_P = x$P_f, refractive_Z_of_risk_allele = zr, refractive_P = x$P_r,
      direction = ifelse(zr > 0, "more hyperopic", "more myopic"))
  }
}
res <- rbindlist(out); dd <- rbindlist(dirs)
print(res[, lapply(.SD, function(x) if (is.numeric(x)) signif(x, 3) else x)])
print(dd[, lapply(.SD, function(x) if (is.numeric(x)) signif(x, 3) else x)])
fwrite(res, file.path(BASE, "analysis/coloc_17q25_fecd_refractive.tsv"), sep = "\t")
fwrite(dd, file.path(BASE, "analysis/coloc_17q25_fecd_refractive_direction.tsv"), sep = "\t")
