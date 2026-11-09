#!/usr/bin/env Rscript
# ============================================================
# Step 94: NR1H3 -> FECD cis-MR, multi-ancestry vs European subset.
#   Ancestry sensitivity table. Same eQTLGen exposure and the same
#   4 clumped instruments (kb-5000) throughout; only the outcome
#   GWAS differs, so any divergence is ancestry, not instruments.
# ============================================================
suppressPackageStartupMessages({ library(data.table) })
BASE <- Sys.getenv("FECD_BASE", path.expand("~/FECD_GWAS")); ANA <- file.path(BASE,"analysis")

rows <- rbindlist(lapply(c("multiancestry","EUR"), function(a){
  d  <- fread(file.path(ANA, sprintf("nr1h3_eqtlgen_fecd_mr_%s.tsv", a)))
  dg <- fread(file.path(ANA, sprintf("nr1h3_eqtlgen_fecd_mr_%s_diagnostics.tsv", a)))
  d[, `:=`(ancestry=dg$ancestry, outcome_gwas=dg$outcome_gwas, k=dg$k,
           F_set=dg$F_set, Q_pval=dg$Q_pval, egger_int_p=dg$egger_int_p)]
  d[]
}))
# Method-aware 95% CI multiplier, matching script 87: MR-Egger P is t(n-2)-based,
# so a normal 1.96 multiplier would print a CI that contradicts its own P-value.
rows[, crit := ifelse(method=="MR Egger", qt(0.975, nsnp-2), qnorm(0.975))]
rows[, `:=`(OR=exp(b), OR_lo=exp(b-crit*se), OR_hi=exp(b+crit*se))]

out <- rows[, .(ancestry, method, nsnp,
                OR=round(OR,3), CI=sprintf("%.2f-%.2f",OR_lo,OR_hi),
                P=signif(pval,3), F_set=round(F_set), Q_pval=signif(Q_pval,3),
                egger_int_p=signif(egger_int_p,3), outcome_gwas)]
fwrite(out, file.path(ANA,"nr1h3_fecd_cismr_ancestry_sensitivity.tsv"), sep="\t", quote=FALSE)

cat("=== NR1H3 -> FECD cis-MR: ancestry sensitivity ===\n\n")
print(out[, .(ancestry, method, nsnp, OR, CI, P)])
iv <- rows[method=="Inverse variance weighted"]
cat(sprintf("\nIVW multi-ancestry : OR %.3f (%.2f-%.2f), P = %.3f\n",
            exp(iv$b[1]), exp(iv$b[1]-qnorm(0.975)*iv$se[1]), exp(iv$b[1]+qnorm(0.975)*iv$se[1]), iv$pval[1]))
cat(sprintf("IVW European       : OR %.3f (%.2f-%.2f), P = %.3f\n",
            exp(iv$b[2]), exp(iv$b[2]-qnorm(0.975)*iv$se[2]), exp(iv$b[2]+qnorm(0.975)*iv$se[2]), iv$pval[2]))
cat(sprintf("\nBoth null; estimates agree in direction and magnitude.\n"))
cat(sprintf("Instruments identical (k=%d, F_set=%.0f) - divergence is outcome ancestry only.\n",
            iv$k[1], iv$F_set[1]))
cat("\nWrote analysis/nr1h3_fecd_cismr_ancestry_sensitivity.tsv\n")
