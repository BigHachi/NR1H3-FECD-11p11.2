#!/usr/bin/env Rscript
# ============================================================
# Step 110 (review point 7): NR1H3 (and ZNF408, TCF4; TSPAN10 from release v14.9) expression in FECD versus
#   control corneal endothelium, published bulk RNA-seq (Wieben et al., Mayo):
#   GSE112201 (24 FECD: 18 CTG18.1 RE+, 6 RE-; 4 controls; batches 1-3)
#   + the 2 controls of GSE101872 (whose 5 FECD samples are the GSE112201
#   batch-1 donors and are therefore not added again).
#   DESeq2, design ~ batch + status (FECD vs control). The CTG18.1 expansion status
#   (RE+ / RE-) is recorded in the sample sheet (Supplementary Table S14) only; no subgroup
#   model is fitted (the ~ batch + group comparisons were removed on 30 Oct 2026).
# ============================================================
suppressPackageStartupMessages({ library(data.table); library(DESeq2) })
D <- file.path(Sys.getenv("FECD_BASE", path.expand("~/FECD_GWAS")), "data/endothelium_rnaseq")
a <- fread(file.path(D, "GSE112201_GeneCounts.tsv.gz"))
b <- fread(file.path(D, "GSE101872_GeneCounts.tsv.gz"))
stopifnot(identical(a$GeneID, b$GeneID))
# check donor overlap: GSE101872 2011-024 vs GSE112201 RNA16 (title "2011-024 RNA16")
cat(sprintf("overlap check 2011-024 vs RNA16: identical=%s, cor=%.4f\n",
            identical(b[["2011-024"]], a[["RNA16"]]), cor(b[["2011-024"]], a[["RNA16"]])))
# sample sheet (from series matrices)
meta <- data.table(
  sample = c("RNA16","RNA15","RNA19","RNA10","RNA20","RNA79",
             "RNA90","RNA100","RNA102","RNA111","RNA112","RNA120","RNA27","RNA28","RNA124","RNA121","RNA114","RNA87",
             "RNA183","RNA193","RNA141","RNA142","RNA166","RNA167","RNA115","RNA136","RNA84","RNA89"),
  group  = c(rep("REpos",4),"REneg","REneg", rep("REpos",11),"REneg", "REpos","REpos","REpos","REneg","REneg","REneg", rep("Control",4)),
  batch  = c(rep("1",6), rep("2",12), rep("3",6), rep("2",4)))
cnt <- as.matrix(a[, meta$sample, with = FALSE]); rownames(cnt) <- a$GeneID
cnt <- cbind(cnt, GSE101872_Control_1 = b$Control_1, GSE101872_Control_2 = b$Control_2)
meta <- rbind(meta, data.table(sample = c("GSE101872_Control_1","GSE101872_Control_2"), group = "Control", batch = "1"))
meta[, status := factor(ifelse(group == "Control", "Control", "FECD"), levels = c("Control","FECD"))]
meta[, batch := factor(batch)]
# sample sheet used by every model (Supplementary Table S14)
fwrite(meta[, .(sample_in_analysis = sample,
                group = c(REpos = "FECD, CTG18.1 expansion-positive", REneg = "FECD, expansion-negative",
                          Control = "control")[as.character(group)],
                batch_in_model = as.character(batch))],
       file.path(Sys.getenv("FECD_BASE", path.expand("~/FECD_GWAS")), "analysis/nr1h3_endothelium_samples.tsv"), sep = "\t")
cat(sprintf("samples: %d FECD (%d RE+, %d RE-), %d controls\n",
            sum(meta$status=="FECD"), sum(meta$group=="REpos"), sum(meta$group=="REneg"), sum(meta$status=="Control")))
keep <- rowSums(cnt >= 10) >= 3; cnt <- cnt[keep, ]
dds <- DESeqDataSetFromMatrix(cnt, as.data.frame(meta), ~ batch + status)
dds <- DESeq(dds, quiet = TRUE)
norm <- counts(dds, normalized = TRUE)
out <- list()
for (g in c("NR1H3","ZNF408","TCF4","TSPAN10")) {   # TSPAN10 added 13 Nov 2026 (release v14.9; manuscript Results)
  if (!g %in% rownames(dds)) { cat(g, "filtered (low counts)\n"); next }
  r <- results(dds, name = "status_FECD_vs_Control")[g, ]
  out[[g]] <- data.table(gene = g,
    median_norm_FECD = median(norm[g, meta$status=="FECD"]), median_norm_control = median(norm[g, meta$status=="Control"]),
    log2FC_FECD_vs_control = r$log2FoldChange, lfcSE = r$lfcSE, P = r$pvalue, padj_genomewide = r$padj)
}
res <- rbindlist(out); print(res[, lapply(.SD, function(x) if (is.numeric(x)) signif(x, 3) else x)])
fwrite(res, file.path(Sys.getenv("FECD_BASE", path.expand("~/FECD_GWAS")), "analysis/nr1h3_endothelium_fecd_vs_control.tsv"), sep = "\t")
