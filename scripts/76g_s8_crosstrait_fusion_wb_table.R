#!/usr/bin/env Rscript
# ============================================================
# Step 76g (Oct 19 2026): Supplementary Table 8 source — every chr11 gene tested in the
#   cross-trait FUSION whole-blood TWAS (step 76d) with its conditional/joint status.
#   Replaces the hand-built Aug 2026 table so that S8 is regenerated with step 76d
#   (needed after the step-55 refractive-error merge fix).
#   Symbols: the HGNC symbols published in S8 (data/derived/s8_chr11_wb_gene_symbols.tsv),
#   falling back to GENCODE v26 names from the GTEx v8 MASHR 'extra' tables.
#   Output: analysis/suppl_crosstrait_fusion_wb_full.tsv
# ============================================================
suppressPackageStartupMessages({ library(data.table); library(DBI); library(RSQLite) })
BASE <- Sys.getenv("FECD_BASE", path.expand("~/FECD_GWAS")); setwd(BASE)
RES  <- "results/twas_wb_conditional"

dbs <- c(Sys.glob("data/predictdb_all/eqtl/mashr/mashr_*.db"), Sys.glob("data/predictdb/mashr_*.db"))
sym <- unique(rbindlist(lapply(dbs, function(f) {
  con <- dbConnect(SQLite(), f); on.exit(dbDisconnect(con))
  as.data.table(dbGetQuery(con, "SELECT gene AS ENSG, genename AS symbol FROM extra"))
})))[!duplicated(ENSG)]
pub <- fread("data/derived/s8_chr11_wb_gene_symbols.tsv", colClasses = "character")
sym <- rbind(pub, sym[!ENSG %in% pub$ENSG])

out <- rbindlist(lapply(c("glaucoma", "myopia", "iop"), function(tr) {
  d   <- fread(file.path(RES, sprintf("%s_wb_chr11.dat", tr)))[, .(ENSG = ID, TWAS.Z, TWAS.P)]
  inc <- fread(file.path(RES, sprintf("%s_wb_chr11.joint.joint_included.dat", tr)))[
           , .(ENSG = ID, cond_status = "conditionally_independent", joint_or_cond_Z = JOINT.Z, joint_or_cond_P = JOINT.P)]
  drp <- fread(file.path(RES, sprintf("%s_wb_chr11.joint.joint_dropped.dat", tr)))[
           , .(ENSG = ID, cond_status = "conditioned_away", joint_or_cond_Z = COND.Z, joint_or_cond_P = COND.P)]
  x <- merge(d, rbind(inc, drp), by = "ENSG", all.x = TRUE)
  x <- merge(x, sym, by = "ENSG", all.x = TRUE)
  x[, trait := tr][order(TWAS.P)]
}))
setcolorder(out, c("trait", "ENSG", "symbol", "TWAS.Z", "TWAS.P", "cond_status", "joint_or_cond_Z", "joint_or_cond_P"))
out[is.na(symbol), symbol := ""]
out[is.na(cond_status), cond_status := "not_in_joint_model"]   # in neither FUSION joint file
fwrite(out, "analysis/suppl_crosstrait_fusion_wb_full.tsv", sep = "\t", quote = FALSE)
cat(sprintf("S8 rows: %s\n", paste(out[, .N, by = trait][, sprintf("%s %d", trait, N)], collapse = ", ")))
print(out[ENSG == "ENSG00000025434.18"])
