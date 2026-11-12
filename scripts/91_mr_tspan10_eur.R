#!/usr/bin/env Rscript
# ============================================================
# Step 91-TSPAN10 (13 Nov 2026: copy of the NR1H3 script with only the gene changed to TSPAN10, ENSG00000182612, 17q25.3,
#   and the LD reference to chromosome 17; run by run_full_pipeline.sh from release v14.8; S7 TSPAN10 rows): cis-MR — TSPAN10 -> FECD (EUR subset)
#   Exposure: eQTLGen blood cis-eQTL (N=31,684) — well-powered
#             instrument source (vs underpowered GTEx N~500).
#   Outcome:  FECD GWAS, EUROPEAN subset (GCST90475865;
#             3,473 cases / 445,239 controls)
#   beta/SE derived from Z, eQTLGen MAF, N (Zhu 2016 standardised formula)
#   Clumped (PLINK 1000G EUR chr17); set-level F (Pierce/Burgess).
# ============================================================
suppressPackageStartupMessages({ library(data.table); library(TwoSampleMR) })
BASE  <- Sys.getenv("FECD_BASE", path.expand("~/FECD_GWAS"))
PLINK <- file.path(BASE,"software/plink/plink")
LD    <- file.path(BASE,"data/ld_ref/LDREF/1000G.EUR.17")
N_EQTL <- 31684L

# ---- eQTLGen TSPAN10 cis-eQTLs ----
e <- fread(file.path(BASE,"data/eqtl_eqtlgen/tspan10_eqtlgen.txt"))
e <- e[, .(SNP, A1e=AssessedAllele, A2e=OtherAllele, Z=Zscore, N=NrSamples)]
# eQTLGen frequency of the assessed allele (exposure source; 10 Nov 2026, review: up to release v13 the
# outcome GWAS frequency was used). Standardises the exposure Z and defines the MAF filter.
af <- fread(file.path(BASE,"data/eqtl_eqtlgen/tspan10_eqtlgen_af.tsv"))
e <- merge(e, af[, .(SNP, AlleleA, AlleleB, AlleleB_all)], by="SNP")
e[, eaf_eqtl := fifelse(A1e==AlleleB, AlleleB_all, fifelse(A1e==AlleleA, 1-AlleleB_all, NA_real_))]
e <- e[is.finite(eaf_eqtl), .(SNP, A1e, A2e, Z, N, eaf_eqtl)]

# ---- FECD GWAS (outcome) matched rows ----
g <- fread(file.path(BASE,"data/gwas/fecd_eur_tspan10_eqtlgen_region.tsv"))
setnames(g, c("A1","A2","BETA","SE","P","FREQ"),
            c("A1_gwas","A2_gwas","BETA_gwas","SE_gwas","P_gwas","FREQ_gwas"))

m <- merge(e, g, by="SNP")
# MAF of the instrument in eQTLGen (the exposure sample)
m[, maf := pmin(eaf_eqtl, 1-eaf_eqtl)]
m <- m[maf>0.01 & is.finite(Z)]

# ---- beta/se of eQTL from Z, MAF, N (standardised; Zhu et al. 2016) ----
m[, denom := sqrt(2*maf*(1-maf)*(N + Z^2))]
m[, beta_eqtl := Z/denom]      # per AssessedAllele (A1e)
m[, se_eqtl   := 1/denom]

# restrict to genome-wide-significant eQTL instruments, then clump
cands <- m[abs(Z) > 5.45]      # p<5e-8
cf <- tempfile(fileext=".tsv"); fwrite(data.table(SNP=cands$SNP, P=2*pnorm(-abs(cands$Z))), cf, sep="\t")
co <- tempfile()
# clump-kb 5000 (not 1000): cis window spans >1Mb, so a 1000kb window left a long-range
# LD pair (rs7395581/rs1316604, r2=0.166, 1.03Mb apart) through. 5000 -> truly independent set.
system(sprintf("%s --bfile %s --clump %s --clump-p1 5e-8 --clump-p2 1 --clump-r2 0.01 --clump-kb 5000 --out %s --memory 4000 --silent",
               PLINK, LD, cf, co), ignore.stdout=TRUE, ignore.stderr=TRUE)
clumped <- fread(paste0(co,".clumped"))$SNP
unlink(c(cf, paste0(co,c(".clumped",".log",".nosex"))))
ins <- m[SNP %in% clumped]
cat(sprintf("eQTLGen TSPAN10 cis-eQTLs: %d in FECD GWAS | genome-wide-sig: %d | clumped independent: %d\n",
            nrow(m), nrow(cands), nrow(ins)))

# ---- set-level F (Pierce/Burgess) ----
k  <- nrow(ins)
ins[, r2 := 2*maf*(1-maf)*beta_eqtl^2]
R2 <- min(sum(ins$r2), 0.999)
Fset <- (R2*(N_EQTL - k - 1))/((1-R2)*k)
cat(sprintf("Instrument SET: k=%d  R2=%.4f  N=%d  ->  F_set=%.1f  [%s]\n",
            k, R2, N_EQTL, Fset, ifelse(Fset>10,"PASS","EXCLUDE")))
if (Fset <= 10) { cat("Excluded: weak instruments\n"); quit(status=0) }

# ---- MR ----
exp_dat <- format_data(as.data.frame(ins), type="exposure", snp_col="SNP",
  beta_col="beta_eqtl", se_col="se_eqtl",
  effect_allele_col="A1e", other_allele_col="A2e", eaf_col="eaf_eqtl")
exp_dat$exposure <- "TSPAN10 (eQTLGen blood)"
exp_dat$id.exposure <- "TSPAN10_eQTLGen"   # fixed IDs (format_data draws random ones) -> byte-identical reruns
exp_dat$pval.exposure <- 2*pnorm(-abs(ins$Z))
out_dat <- format_data(as.data.frame(ins), type="outcome", snp_col="SNP",
  beta_col="BETA_gwas", se_col="SE_gwas", pval_col="P_gwas",
  effect_allele_col="A1_gwas", other_allele_col="A2_gwas", eaf_col="FREQ_gwas")
out_dat$outcome <- "FECD"
out_dat$id.outcome <- "FECD_EUR"
dat <- harmonise_data(exp_dat, out_dat, action=2); dat <- dat[dat$mr_keep,]
cat(sprintf("Valid instruments after harmonisation: %d\n\n", nrow(dat)))

set.seed(2024)  # weighted median/mode SEs are bootstrapped -> seed for reproducible P
res <- mr(dat, method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median","mr_weighted_mode"))
res$OR <- exp(res$b); res$OR_lo <- exp(res$b-1.96*res$se); res$OR_hi <- exp(res$b+1.96*res$se)
cat("=== TSPAN10 (eQTLGen) -> FECD: MR estimates ===\n")
print(as.data.table(res)[,.(method, nsnp, b=round(b,3), se=round(se,3),
       OR=round(OR,3), CI=sprintf("%.2f-%.2f",OR_lo,OR_hi), P=signif(pval,3))])
iv <- res[res$method=="Inverse variance weighted",]
cat(sprintf("\n>>> IVW: OR=%.3f (%.2f-%.2f), P=%.2e -> higher TSPAN10 %s FECD risk\n",
            exp(iv$b), exp(iv$b-1.96*iv$se), exp(iv$b+1.96*iv$se), iv$pval,
            ifelse(iv$b>0,"INCREASES","DECREASES")))
het <- mr_heterogeneity(dat); ple <- mr_pleiotropy_test(dat)
cat(sprintf("Heterogeneity IVW Q P=%.3g | Egger intercept=%.4f (P=%.3g)\n",
            het$Q_pval[het$method=="Inverse variance weighted"], ple$egger_intercept, ple$pval))
fwrite(as.data.table(res), file.path(BASE,"analysis/tspan10_eqtlgen_fecd_mr_EUR.tsv"), sep="\t", quote=FALSE)

# ---- diagnostics sidecar: so downstream figures never hardcode these ----
diag <- data.table(
  ancestry      = "EUR",
  outcome_gwas  = "FECD GWAS European subset (GCST90475865; 3,473/445,239)",
  k             = k,
  R2            = R2,
  N_eqtl        = N_EQTL,
  F_set         = Fset,
  clump_kb      = 5000L,
  Q_pval        = het$Q_pval[het$method=="Inverse variance weighted"],
  egger_int     = ple$egger_intercept,
  egger_int_p   = ple$pval)
fwrite(diag, file.path(BASE,"analysis/tspan10_eqtlgen_fecd_mr_EUR_diagnostics.tsv"), sep="\t", quote=FALSE)
cat(sprintf("Wrote analysis/tspan10_eqtlgen_fecd_mr_EUR{,_diagnostics}.tsv\n"))
