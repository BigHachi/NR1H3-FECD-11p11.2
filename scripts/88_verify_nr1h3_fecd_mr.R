#!/usr/bin/env Rscript
# ============================================================
# Step 88: VERIFICATION of NR1H3 (eQTLGen) -> FECD cis-MR
#   Full audit of clumping + allele harmonisation.
# ============================================================
suppressPackageStartupMessages({ library(data.table); library(TwoSampleMR) })
BASE  <- Sys.getenv("FECD_BASE", path.expand("~/FECD_GWAS"))
PLINK <- file.path(BASE,"software/plink/plink")
LD    <- file.path(BASE,"data/ld_ref/LDREF/1000G.EUR.11")
N_EQTL <- 31684L

e <- fread(file.path(BASE,"data/eqtl_eqtlgen/nr1h3_eqtlgen.txt"))
e <- e[, .(SNP, A1e=AssessedAllele, A2e=OtherAllele, Z=Zscore, N=NrSamples)]
# eQTLGen frequency of the assessed allele (exposure source; 10 Nov 2026, review: up to release v13 the
# outcome GWAS frequency was used). Standardises the exposure Z and defines the MAF filter.
af <- fread(file.path(BASE,"data/eqtl_eqtlgen/nr1h3_eqtlgen_af.tsv"))
e <- merge(e, af[, .(SNP, AlleleA, AlleleB, AlleleB_all)], by="SNP")
e[, eaf_eqtl := fifelse(A1e==AlleleB, AlleleB_all, fifelse(A1e==AlleleA, 1-AlleleB_all, NA_real_))]
e <- e[is.finite(eaf_eqtl), .(SNP, A1e, A2e, Z, N, eaf_eqtl)]
g <- fread(file.path(BASE,"data/gwas/fecd_nr1h3_eqtlgen_region.tsv"))
setnames(g, c("A1","A2","BETA","SE","P","FREQ"),
            c("A1_gwas","A2_gwas","BETA_gwas","SE_gwas","P_gwas","FREQ_gwas"))
m <- merge(e, g, by="SNP")
# MAF of the instrument in eQTLGen (the exposure sample)
m[, maf := pmin(eaf_eqtl, 1-eaf_eqtl)]
m <- m[maf>0.01 & is.finite(Z)]
m[, denom := sqrt(2*maf*(1-maf)*(N + Z^2))]
m[, beta_eqtl := Z/denom]
m[, se_eqtl   := 1/denom]

cands <- m[abs(Z) > 5.45]
cat(sprintf("SNPs in both eQTLGen & FECD: %d | genome-wide-sig (|Z|>5.45): %d\n", nrow(m), nrow(cands)))
cf <- tempfile(fileext=".tsv"); fwrite(data.table(SNP=cands$SNP, P=2*pnorm(-abs(cands$Z))), cf, sep="\t")
co <- tempfile()
system(sprintf("%s --bfile %s --clump %s --clump-p1 5e-8 --clump-p2 1 --clump-r2 0.01 --clump-kb 5000 --out %s --memory 4000 --silent",
               PLINK, LD, cf, co), ignore.stdout=TRUE, ignore.stderr=TRUE)
clumped <- fread(paste0(co,".clumped"))$SNP
ins <- m[SNP %in% clumped]
cat(sprintf("\n>>> CLUMPED independent instruments: %d\n", nrow(ins)))
cat("   SNPs:", paste(ins$SNP, collapse=", "), "\n\n")

# ---- palindrome check ----
comp <- c(A="T",T="A",C="G",G="C")
ins[, palindromic := (A1e==comp[A2e])]
cat("=== Per-SNP allele table (exposure eQTLGen vs outcome FECD) ===\n")
print(ins[, .(SNP, eQTL_A1=A1e, eQTL_A2=A2e, Z=round(Z,2), b_eqtl=round(beta_eqtl,4),
              FECD_A1=A1_gwas, FECD_A2=A2_gwas, FECD_beta=round(BETA_gwas,4),
              FECD_FRQ=round(FREQ_gwas,3), palindromic)])

# ---- independent LD check among clumped SNPs ----
sf <- tempfile(fileext=".txt"); writeLines(ins$SNP, sf); ldo <- tempfile()
system(sprintf("%s --bfile %s --extract %s --r2 --ld-window-r2 0 --ld-window 99999 --ld-window-kb 1000000 --out %s --memory 4000 --silent",
               PLINK, LD, sf, ldo), ignore.stdout=TRUE, ignore.stderr=TRUE)
cat("\n=== Pairwise LD (r2) among the clumped instruments (should all be < 0.01) ===\n")
if (file.exists(paste0(ldo,".ld"))) {
  ld <- fread(paste0(ldo,".ld"))
  print(ld[, .(SNP_A, SNP_B, R2=round(R2,4))])
  cat(sprintf("MAX pairwise r2 = %.4f\n", max(ld$R2)))
} else cat("(no pairs reported = all beyond window / no LD)\n")

# ---- harmonise & MR ----
exp_dat <- format_data(as.data.frame(ins), type="exposure", snp_col="SNP",
  beta_col="beta_eqtl", se_col="se_eqtl", effect_allele_col="A1e",
  other_allele_col="A2e", eaf_col="eaf_eqtl")
exp_dat$exposure <- "NR1H3"; exp_dat$pval.exposure <- 2*pnorm(-abs(ins$Z))
out_dat <- format_data(as.data.frame(ins), type="outcome", snp_col="SNP",
  beta_col="BETA_gwas", se_col="SE_gwas", pval_col="P_gwas",
  effect_allele_col="A1_gwas", other_allele_col="A2_gwas", eaf_col="FREQ_gwas")
out_dat$outcome <- "FECD"
dat <- harmonise_data(exp_dat, out_dat, action=2)
cat("\n=== HARMONISED data (TwoSampleMR action=2) ===\n")
print(as.data.table(dat)[, .(SNP, ea=effect_allele.exposure, oa=other_allele.exposure,
    b.exp=round(beta.exposure,4), b.out=round(beta.outcome,4),
    eaf=round(eaf.exposure,3), remove, palindromic, ambiguous, mr_keep)])
cat(sprintf("\nInstruments kept for MR: %d of %d\n", sum(dat$mr_keep), nrow(dat)))
set.seed(2024)
res <- mr(dat[dat$mr_keep,], method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median","mr_weighted_mode"))
res$OR <- exp(res$b)
cat("\n=== MR RESULT ===\n")
print(as.data.table(res)[, .(method, nsnp, b=round(b,4), se=round(se,4), OR=round(OR,3), P=signif(pval,3))])
