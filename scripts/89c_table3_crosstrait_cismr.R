#!/usr/bin/env Rscript
# ============================================================
# Step 89c (Oct 10 2026): Table 3 — NR1H3 eQTLGen cis-MR on glaucoma, IOP and
#   refractive error, all four estimators + diagnostics.
#   Same instruments and outcome scaling as step 89b (Figure 6): eQTLGen NR1H3,
#   |Z| > 5.45 (P < 5e-8), MAF > 0.01, PLINK clump r2 < 0.01 / 5,000 kb (k = 4).
#   Effect units: glaucoma log-OR and IOP mmHg (GWAS beta/SE as published);
#   refractive error standardised from Z with per-variant N (Z-only GWAS), as in 89b.
#   Exposure betas are standardised (per SD of expression). Weighted median /
#   mode use a bootstrap, so the seed is fixed for reproducibility.
#   Replaces the inline July 2026 table (analysis/nr1h3_cismr_long.tsv), whose
#   refractive-error rows used a fixed N and whose bootstrap was unseeded.
#   Output: analysis/table3_nr1h3_crosstrait_cismr.tsv
#   Exposure N (Oct 20 2026, review): the eQTLGen Z is standardised with each instrument's NrSamples, as in
#   the FECD cis-MR (steps 83/88/91): 31,676 / 31,684 / 31,346 / 21,502 (rs4647760, 21 of 37 cohorts).
#   Until Oct 20 a fixed N = 31,684 was used (previous Table 3: results/bak_exposureN_fixed_20261020/). The
#   per-instrument exposure beta/SE under both conventions are written to
#   analysis/nr1h3_cismr_instruments_exposureN.tsv (the comparison reported to the reviewer), with the
#   eQTLGen and outcome-GWAS MAF and the v13 exposure beta/SE (outcome-GWAS MAF) for comparison.
# ============================================================
suppressPackageStartupMessages({ library(data.table); library(TwoSampleMR) })
set.seed(1)
BASE  <- Sys.getenv("FECD_BASE", path.expand("~/FECD_GWAS")); ANA <- file.path(BASE, "analysis")
PLINK <- file.path(BASE, "software/plink/plink"); LD <- file.path(BASE, "data/ld_ref/LDREF/1000G.EUR.11")
N_EQTL <- 31684L
z2b <- function(z, maf, n) { d <- sqrt(2 * maf * (1 - maf) * (n + z^2)); list(b = z / d, se = 1 / d) }

e <- fread(file.path(BASE, "data/eqtl_eqtlgen/nr1h3_eqtlgen.txt"))[
       , .(SNP, A1e = AssessedAllele, A2e = OtherAllele, Z = Zscore, Ns = NrSamples)]
# eQTLGen frequency of the assessed allele (exposure source; 10 Nov 2026, review). It standardises the exposure Z
# and defines the MAF filter; up to release v13 the outcome GWAS frequency was used for both. The outcome GWAS
# frequency still standardises the Z-only refractive-error outcome.
af <- fread(file.path(BASE, "data/eqtl_eqtlgen/nr1h3_eqtlgen_af.tsv"))
e <- merge(e, af[, .(SNP, AlleleA, AlleleB, AlleleB_all)], by = "SNP")
e[, eaf_eqtl := fifelse(A1e == AlleleB, AlleleB_all, fifelse(A1e == AlleleA, 1 - AlleleB_all, NA_real_))]
e <- e[is.finite(eaf_eqtl), .(SNP, A1e, A2e, Z, Ns, eaf_eqtl)]
# Refractive-error N = METAL Weight per variant, the N column of the step-55 merge (Oct 19 2026;
# replaces values hard-coded for the four instruments, which it reproduces).
outc <- list(
  Glaucoma = list(f = "data/gwas/glaucoma_gwas_full.tsv", z_only = FALSE),
  IOP      = list(f = "data/gwas/iop_gwas_full.tsv",      z_only = FALSE),
  Myopia   = list(f = "data/gwas/myopia_gwas_full.tsv",   z_only = TRUE)   # BETA = Z, SE = 1
)
res <- list(); inst <- list()
for (nm in names(outc)) {
  o <- outc[[nm]]; g <- fread(file.path(BASE, o$f))
  g <- g[SNP %in% e$SNP, .(SNP, A1_g = A1, A2_g = A2, BETA, SE, FREQ, Nv = if (o$z_only) as.numeric(N) else NA_real_)]
  m <- merge(e, g, by = "SNP"); m[, `:=`(maf_e = pmin(eaf_eqtl, 1 - eaf_eqtl), maf = pmin(FREQ, 1 - FREQ))]
  m[, Zo := if (o$z_only) BETA else BETA / SE]
  m <- m[maf_e > 0.01 & is.finite(Z) & abs(Z) > 5.45 & is.finite(Zo)]
  be <- z2b(m$Z, m$maf_e, m$Ns); m[, `:=`(be = be$b, se_e = be$se)]
  if (o$z_only) { bo <- z2b(m$Zo, m$maf, m$Nv); m[, `:=`(bo = bo$b, se_o = bo$se)] }
  else m[, `:=`(bo = BETA, se_o = SE)]
  cf <- tempfile(fileext = ".tsv"); fwrite(data.table(SNP = m$SNP, P = 2 * pnorm(-abs(m$Z))), cf, sep = "\t"); co <- tempfile()
  system(sprintf("%s --bfile %s --clump %s --clump-p1 5e-8 --clump-p2 1 --clump-r2 0.01 --clump-kb 5000 --out %s --memory 4000 --silent",
                 PLINK, LD, cf, co), ignore.stdout = TRUE, ignore.stderr = TRUE)
  ins <- m[SNP %in% fread(paste0(co, ".clumped"))$SNP]
  bf <- z2b(ins$Z, ins$maf_e, N_EQTL); bv <- z2b(ins$Z, ins$maf_e, ins$Ns); b13 <- z2b(ins$Z, ins$maf, ins$Ns)
  inst[[nm]] <- ins[, .(outcome = nm, SNP, exposure_Z = Z, NrSamples = Ns, maf_eqtlgen = maf_e, maf_outcome_gwas = maf,
                        beta_fixedN = bf$b, se_fixedN = bf$se, beta_variantN = bv$b, se_variantN = bv$se,
                        beta_v13_outcome_maf = b13$b, se_v13_outcome_maf = b13$se)]
  ex <- format_data(as.data.frame(ins), type = "exposure", snp_col = "SNP", beta_col = "be", se_col = "se_e",
                    effect_allele_col = "A1e", other_allele_col = "A2e", eaf_col = "eaf_eqtl"); ex$exposure <- "NR1H3"
  ou <- format_data(as.data.frame(ins), type = "outcome", snp_col = "SNP", beta_col = "bo", se_col = "se_o",
                    effect_allele_col = "A1_g", other_allele_col = "A2_g", eaf_col = "FREQ"); ou$outcome <- nm
  dat <- harmonise_data(ex, ou, action = 2); dat <- dat[dat$mr_keep, ]
  r <- as.data.table(mr(dat, method_list = c("mr_ivw", "mr_weighted_median", "mr_weighted_mode", "mr_egger_regression")))
  het <- as.data.table(mr_heterogeneity(dat, method_list = "mr_ivw")); ple <- mr_pleiotropy_test(dat)
  bx <- abs(dat$beta.exposure); sx <- dat$se.exposure
  isq <- { w <- 1 / sx^2; mu <- sum(w * bx) / sum(w); Q <- sum(w * (bx - mu)^2); max(0, (Q - (length(bx) - 1)) / Q) }
  res[[nm]] <- r[, .(outcome = nm, method, nsnp, b, se, pval, Isq_GX = isq,
                     egger_intercept_pval = ple$pval, cochran_Q_pval = het$Q_pval)]
}
d <- rbindlist(res)
fwrite(d, file.path(ANA, "table3_nr1h3_crosstrait_cismr.tsv"), sep = "\t", quote = FALSE)
fwrite(rbindlist(inst), file.path(ANA, "nr1h3_cismr_instruments_exposureN.tsv"), sep = "\t", quote = FALSE)
print(d[, .(outcome, method, nsnp, b = signif(b, 3), se = signif(se, 2), pval = signif(pval, 3),
            Isq_GX = round(Isq_GX, 3), egger_int_p = signif(egger_intercept_pval, 3), Q_p = signif(cochran_Q_pval, 3))])
