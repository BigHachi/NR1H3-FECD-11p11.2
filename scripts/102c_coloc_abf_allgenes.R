#!/usr/bin/env Rscript
# ============================================================
# Step 102c: as step 102 but gene list = data/tophits21/genes_coloc.tsv
#   (European S-MultiXcan ranks 1-23, down to NR1H3, + four genes from earlier versions of S6; ranks written by
#   step 19f). Was:
# Step 102: coloc.abf on FULL cis summary statistics for every FECD
#           EUR S-MultiXcan gene down to NR1H3 (ranks 1-21), in the
#           3 TWAS tissues.
#   eQTL   : eQTL Catalogue r6 = GTEx v8, ALL variant-gene pairs (step 101)
#            -> replaces the GTEx-API (significant pairs only, ~1.6%
#               coverage) and FUSION-internal (weight SNPs only) colocs,
#               both shown INVALID in the Sep 30 2026 audit.
#   GWAS   : FECD EUR GCST90475865 (primary) + multi-ancestry GCST90480062.
#   Region : exactly the cis window eQTL Catalogue tested for the gene.
#   Match  : chr + pos + alleles (either orientation; coloc.abf uses Z^2
#            so sign alignment is not required, but allele match stops
#            two different variants at one position being merged).
#   Priors : coloc defaults p1=p2=1e-4, p12=1e-5; PP4 at p12=5e-6
#            (Wallace 2020) also reported.
# ============================================================
suppressPackageStartupMessages({ library(data.table); library(coloc) })
BASE <- Sys.getenv("FECD_BASE", path.expand("~/FECD_GWAS"))
D    <- file.path(BASE, "data/tophits21")
OUT  <- file.path(BASE, "results/coloc_allpairs_v2"); dir.create(OUT, FALSE, TRUE)
TISSUES <- c("fibroblast", "skin_not_sun_exposed", "esophagus_mucosa")

genes <- fread(file.path(D, "genes_coloc.tsv"))
chrs  <- unique(as.character(genes$chr))
# stream-filter the GWAS to the 21 cis windows (TSS +/- 1.1 Mb) before fread:
# avoids decompressing whole genomes into tmp (disk is near full)
win <- genes[, .(chr = as.character(chr), s = tss - 1.1e6, e = tss + 1.1e6)]
win_awk <- function(chrcol, bpcol) sprintf("awk -F'\\t' 'NR==1 || %s'",
  paste(sprintf("($%d==\"%s\" && $%d>=%d && $%d<=%d)", chrcol, win$chr, bpcol,
                as.integer(pmax(win$s, 1)), bpcol, as.integer(win$e)), collapse = " || "))

# ---- GWAS ----
cat("Loading FECD EUR GWAS...\n")
eur <- fread(cmd = paste("gzip -cd", file.path(BASE, "data/gwas/fecd_eur_GCST90475865.tsv.gz"), "|", win_awk(1, 2)),
             select = c("chromosome","base_pair_location","effect_allele","other_allele",
                        "odds_ratio","p_value","rsid","ci_upper","ci_lower"),
             colClasses = list(character = "chromosome"))
setnames(eur, c("CHR","BP","A1","A2","OR","P","SNP","CIU","CIL"))
eur <- eur[CHR %in% chrs]
eur[, `:=`(BETA = log(as.numeric(OR)),
           SE   = (log(as.numeric(CIU)) - log(as.numeric(CIL))) / (2*qnorm(0.975)))]
eur <- eur[is.finite(BETA) & is.finite(SE) & SE > 0, .(CHR, BP=as.integer(BP), A1, A2, BETA, SE, P=as.numeric(P), SNP)]

cat("Loading FECD multi-ancestry GWAS...\n")
mul <- fread(cmd = paste("gzip -cd", file.path(BASE, "data/gwas/fecd_gwas_full.tsv.gz"), "|", win_awk(2, 3)),
             select = c("SNP","CHR","BP","A1","A2","BETA","SE","P"),
             colClasses = list(character = "CHR"))
mul <- mul[CHR %in% chrs & is.finite(BETA) & is.finite(SE) & SE > 0]
mul[, BP := as.integer(BP)]

GW <- list(EUR           = list(d = eur, N = 448712, cases = 3473),
           multiancestry = list(d = mul, N = 570321, cases = 3928))

read_eqtl <- function(gene, tissue) {
  f <- file.path(D, "allpairs", sprintf("%s_%s.tsv", gene, tissue))
  if (!file.exists(f)) return(NULL)
  e <- fread(f, colClasses = list(character = c("chromosome","ref","alt")))
  if (!nrow(e)) return(NULL)
  e <- e[is.finite(beta) & is.finite(se) & se > 0 & maf > 0 & maf < 1]
  e <- e[!duplicated(variant)]
  e[, N := an / 2]
  e
}

run_one <- function(g, tissue, anc) {
  e <- read_eqtl(g$gene, tissue)
  base <- data.table(rank = g$rank, gene = g$gene, ensg = g$ensg, tissue = tissue, ancestry = anc)
  if (is.null(e)) return(cbind(base, status = "no eQTL Catalogue data for gene/tissue"))
  gw <- GW[[anc]]$d[CHR == as.character(g$chr) & BP >= min(e$position) & BP <= max(e$position)]
  m  <- merge(e[, .(BP = position, ref, alt, b_e = beta, se_e = se, p_e = pvalue, maf, N, rsid)],
              gw, by = "BP", allow.cartesian = TRUE)
  m  <- m[(A1 == alt & A2 == ref) | (A1 == ref & A2 == alt)]
  m[, id := paste(BP, ref, alt, sep = ":")]
  m  <- m[!duplicated(id)]
  if (nrow(m) < 50) return(cbind(base, n_eqtl = nrow(e), n_snps = nrow(m), status = "<50 matched variants"))

  d1 <- list(beta = m$BETA, varbeta = m$SE^2, snp = m$id, type = "cc",
             s = GW[[anc]]$cases / GW[[anc]]$N, N = GW[[anc]]$N)
  d2 <- list(beta = m$b_e, varbeta = m$se_e^2, snp = m$id, type = "quant",
             N = round(median(m$N)), MAF = m$maf)
  r  <- suppressMessages(coloc.abf(d1, d2, p1 = 1e-4, p2 = 1e-4, p12 = 1e-5))
  pp <- r$summary
  r5 <- suppressMessages(coloc.abf(d1, d2, p1 = 1e-4, p2 = 1e-4, p12 = 5e-6))$summary
  post <- as.data.table(r$results)
  top  <- post[which.max(SNP.PP.H4)]
  fwrite(merge(post[, .(id = snp, SNP.PP.H4)], m[, .(id, rsid, b_e, p_e, BETA, P)], by = "id")[order(-SNP.PP.H4)],
         file.path(OUT, sprintf("%s_%s_%s_snp_pp.tsv", g$gene, tissue, anc)), sep = "\t")
  cbind(base, n_eqtl = nrow(e), n_snps = nrow(m),
        min_p_gwas = min(m$P), min_p_eqtl = min(m$p_e),
        PP0 = pp[["PP.H0.abf"]], PP1 = pp[["PP.H1.abf"]], PP2 = pp[["PP.H2.abf"]],
        PP3 = pp[["PP.H3.abf"]], PP4 = pp[["PP.H4.abf"]],
        PP4_p12_5e6 = r5[["PP.H4.abf"]],
        PP4_given_H3H4 = pp[["PP.H4.abf"]] / (pp[["PP.H3.abf"]] + pp[["PP.H4.abf"]]),
        top_snp = m[id == top$snp, rsid][1], top_snp_PP4 = top$SNP.PP.H4,
        status = "OK")
}

res <- list()
for (anc in names(GW)) for (i in seq_len(nrow(genes))) for (t in TISSUES) {
  g <- genes[i]
  r <- tryCatch(run_one(g, t, anc), error = function(e)
         data.table(rank = g$rank, gene = g$gene, ensg = g$ensg, tissue = t, ancestry = anc,
                    status = paste("ERROR:", conditionMessage(e))))
  cat(sprintf("%-9s %2d %-10s %-22s %s\n", anc, g$rank, g$gene, t,
              if (!is.null(r$PP4)) sprintf("n=%5d  minP_gwas=%.1e minP_eqtl=%.1e  PP3=%.3f PP4=%.3f",
                                           r$n_snps, r$min_p_gwas, r$min_p_eqtl, r$PP3, r$PP4) else r$status))
  res[[length(res) + 1]] <- r
}
res <- rbindlist(res, fill = TRUE)
setorder(res, ancestry, rank, tissue)
fwrite(res, file.path(BASE, "analysis/fecd_coloc_allpairs_v2.tsv"), sep = "\t")
cat("\nWrote analysis/fecd_coloc_allpairs_v2.tsv\n")
