#!/usr/bin/env Rscript
# ============================================================
# Step 112: statistical fine-mapping of the FECD association at 11p11.2 (chr11:46.2-48.3 Mb, GRCh38)
#   (a) SuSiE-RSS (multiple causal variants, L = 10) on the European GWAS (GCST90475865) with an
#       LD matrix from 1000 Genomes phase 3 Europeans (n = 503; GRCh38 biallelic release; MAF >= 1%),
#       with the susieR LD-mismatch diagnostic (estimate_s_rss) and allele-flip check (kriging_rss).
#   (b) Single-causal-variant fine-mapping (Wakefield ABF; coloc::finemap.abf) of the European and
#       multi-ancestry GWAS, which needs no LD matrix: 95% credible sets.
#   Inputs : analysis/fecd_eur_11p11_region.tsv, data/ld_1000g_grch38/eur_chr11_46.2-48.3Mb.{bed,bim,fam}
#            (prepared by step 112a), data/gwas/fecd_gwas_full.tsv.gz (multi-ancestry, step 11)
#   Outputs: results/finemap_11p11/{susie_summary.tsv, susie_pip.tsv, abf_credible_sets.tsv, diagnostics.txt}
# ============================================================
suppressPackageStartupMessages({ library(data.table); library(susieR); library(coloc) })
set.seed(2024)
BASE  <- Sys.getenv("FECD_BASE", path.expand("~/FECD_GWAS"))
PLINK <- Sys.getenv("PLINK", file.path(BASE, "software/plink/plink"))
LDP   <- file.path(BASE, "data/ld_1000g_grch38/eur_chr11_46.2-48.3Mb")
OUT   <- file.path(BASE, "results/finemap_11p11"); dir.create(OUT, FALSE, TRUE)
NCASE <- 3473; NCTRL <- 445239; NEFF <- 4 / (1 / NCASE + 1 / NCTRL)
diag  <- file(file.path(OUT, "diagnostics.txt"), "w"); say <- function(...) { m <- sprintf(...); cat(m, "\n"); writeLines(m, diag) }

# ---- European GWAS; SE from the 95% CI where the file gives none ----
g <- fread(file.path(BASE, "analysis/fecd_eur_11p11_region.tsv"), na.strings = c("#NA", "NA", ""))
g <- g[, .(pos = base_pair_location, ea = effect_allele, oa = other_allele, or = odds_ratio, se = standard_error,
           lo = ci_lower, hi = ci_upper, eaf = effect_allele_frequency, p = p_value, rsid)]
g[, c("se", "lo", "hi", "or") := lapply(.SD, as.numeric), .SDcols = c("se", "lo", "hi", "or")]   # an all-blank column is read as logical
g[, beta := log(or)]
g[is.na(se), se := (log(hi) - log(lo)) / (2 * qnorm(0.975))]
g <- g[is.finite(beta) & is.finite(se) & se > 0]
g[, z := beta / se]
say("European GWAS variants in region: %d; min P %.2e (%s); cor(log10 P, log10 P from z) = %.4f",
    nrow(g), min(g$p), g$rsid[which.min(g$p)], cor(log10(g$p), log10(2 * pnorm(-abs(g$z)))))

# ---- (b) ABF credible sets (single causal variant), European and multi-ancestry ----
abf_cs <- function(d, label, s, N) {
  d <- d[!duplicated(id)]
  fm <- as.data.table(finemap.abf(list(beta = d$beta, varbeta = d$se^2, type = "cc", s = s, N = N, snp = d$id)))
  fm <- fm[snp != "null"][order(-SNP.PP)]
  fm[, cum := cumsum(SNP.PP)]
  cs <- fm[seq_len(which(cum >= 0.95)[1])]
  cs <- merge(cs[, .(id = snp, PIP = SNP.PP)], d[, .(id, rsid, pos, p)], by = "id")[order(-PIP)]
  say("ABF %s: 95%% credible set = %d variants spanning %.3f-%.3f Mb; top %s (PIP %.3f, P %.2e); P(no causal variant in region | prior) not modelled",
      label, nrow(cs), min(cs$pos) / 1e6, max(cs$pos) / 1e6, cs$rsid[1], cs$PIP[1], cs$p[1])
  cs[, set := label]
}
g[, id := sprintf("11_%d_%s_%s", pos, pmin(ea, oa), pmax(ea, oa))]
cs_eur <- abf_cs(g, "European", NCASE / (NCASE + NCTRL), NCASE + NCTRL)

mm <- fread(cmd = sprintf("gzip -cd %s | awk -F'\\t' 'NR==1 || ($2==11 && $3>=46200000 && $3<=48300000)'",
                          file.path(BASE, "data/gwas/fecd_gwas_full.tsv.gz")))
mm <- mm[is.finite(BETA) & is.finite(SE) & SE > 0, .(pos = BP, ea = A1, oa = A2, beta = BETA, se = SE, p = P, rsid = SNP)]
mm[, id := sprintf("11_%d_%s_%s", pos, pmin(ea, oa), pmax(ea, oa))]
say("Multi-ancestry GWAS variants in region: %d; min P %.2e (%s)", nrow(mm), min(mm$p), mm$rsid[which.min(mm$p)])
cs_multi <- abf_cs(mm, "Multi-ancestry", 3928 / 570321, 570321)
fwrite(rbind(cs_eur, cs_multi), file.path(OUT, "abf_credible_sets.tsv"), sep = "\t")

# ---- (a) SuSiE-RSS on the European GWAS with 1000G European LD ----
bim <- fread(paste0(LDP, ".bim"), header = FALSE, col.names = c("chr", "ldid", "cm", "pos", "a1", "a2"))
m <- merge(g, bim[, .(pos, ldid, a1, a2)], by = "pos")
m <- m[(ea == a1 & oa == a2) | (ea == a2 & oa == a1)]
m <- m[!duplicated(ldid) & !duplicated(pos)][order(pos)]
m[, z1 := ifelse(ea == a1, z, -z)]          # z for the plink A1 allele, the allele the LD signs refer to
say("Variants matched to the LD panel (position + alleles): %d of %d (min P among matched %.2e, %s)",
    nrow(m), nrow(g), min(m$p), m$rsid[which.min(m$p)])
ids <- file.path(OUT, "ld_ids.txt"); writeLines(m$ldid, ids)
sub <- file.path(OUT, "ld_sub")
system(sprintf("%s --bfile %s --extract %s --make-bed --out %s --silent", PLINK, LDP, ids, sub))
system(sprintf("%s --bfile %s --r square spaces --out %s --silent", PLINK, sub, sub))
sb <- fread(paste0(sub, ".bim"), header = FALSE)$V2
R  <- as.matrix(fread(paste0(sub, ".ld"), header = FALSE)); dimnames(R) <- list(sb, sb)
m  <- m[match(sb, ldid)]; stopifnot(identical(m$ldid, sb))
R[!is.finite(R)] <- 0; diag(R) <- 1

s_mis <- estimate_s_rss(m$z1, R, n = NEFF)
say("LD mismatch estimate s = %.3f (0 = consistent; > 0.2 suggests substantial mismatch between GWAS and LD panel)", s_mis)
kr <- kriging_rss(m$z1, R, n = NEFF)$conditional_dist
flag <- which(kr$logLR > 2 & abs(kr$z) > 2)
say("Possible allele flips / LD outliers (kriging logLR > 2 and |z| > 2): %d", length(flag))
if (length(flag)) { say("  removed: %s", paste(m$rsid[flag], collapse = ", ")); keep <- setdiff(seq_len(nrow(m)), flag); m <- m[keep]; R <- R[keep, keep] }

fit <- susie_rss(z = m$z1, R = R, n = NEFF, L = 10, coverage = 0.95, min_abs_corr = 0.5,
                 estimate_residual_variance = FALSE, max_iter = 500)
cs <- susie_get_cs(fit, Xcorr = R, coverage = 0.95, min_abs_corr = 0.5)
m[, PIP := susie_get_pip(fit)]
fwrite(m[order(-PIP), .(ldid, rsid, pos, ea, oa, p, z = z1, PIP)], file.path(OUT, "susie_pip.tsv"), sep = "\t")
say("SuSiE-RSS converged: %s; credible sets (95%%, purity >= 0.5): %d; max PIP %.3f (%s)",
    fit$converged, length(cs$cs), max(m$PIP), m$rsid[which.max(m$PIP)])
summ <- if (length(cs$cs)) rbindlist(lapply(seq_along(cs$cs), function(k) {
  v <- cs$cs[[k]]; data.table(cs = names(cs$cs)[k], n_variants = length(v), min_abs_corr = cs$purity[k, "min.abs.corr"],
    from_Mb = min(m$pos[v]) / 1e6, to_Mb = max(m$pos[v]) / 1e6, top_rsid = m$rsid[v][which.max(m$PIP[v])],
    top_PIP = max(m$PIP[v]), variants = paste(m$rsid[v], collapse = ";"))
})) else data.table(cs = "none", n_variants = 0L)
fwrite(summ, file.path(OUT, "susie_summary.tsv"), sep = "\t")
print(summ[, setdiff(names(summ), "variants"), with = FALSE])
saveRDS(list(fit = fit, m = m), file.path(OUT, "susie_fit.rds"))
close(diag)
