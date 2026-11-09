#!/usr/bin/env Rscript
# ============================================================
# Step 00 (R): install the R packages used by the pipeline, at the versions
#   used for the manuscript (R 4.5.0, Bioconductor 3.21), then report any
#   package whose installed version differs. Already-installed packages at the
#   recorded version are left alone.
# ============================================================
repos <- "https://cloud.r-project.org"
for (p in c("remotes", "BiocManager"))
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p, repos = repos)

cran <- c(data.table = "1.18.2.1", ggplot2 = "4.0.2", ggrepel = "0.9.8", patchwork = "1.3.2",
          ggtext = "0.2.0", svglite = "2.2.2", ragg = "1.5.2", ggrastr = "1.0.2",   # figure style and export (fig_common.R, 10 Nov 2026)
          jsonlite = "2.0.0", httr = "1.4.8", coloc = "5.2.3", optparse = "1.7.5",
          RColorBrewer = "1.1.3", glmnet = "4.1.10", here = "1.0.2", susieR = "0.14.2")
bioc <- c(org.Hs.eg.db = "3.21.0", AnnotationDbi = "1.70.0", DESeq2 = "1.48.2", biomaRt = "2.64.0")
github <- c(TwoSampleMR = "MRCIEU/TwoSampleMR@b34cc5e",     # 0.7.0
            ieugwasr    = "MRCIEU/ieugwasr@8d34780",        # 1.1.0.9000
            plink2R     = "gabraham/plink2R/plink2R@d74be01")  # 1.1, required by FUSION
have <- function(p, v) requireNamespace(p, quietly = TRUE) && as.character(packageVersion(p)) == v

for (p in names(cran)) if (!have(p, cran[[p]])) remotes::install_version(p, version = cran[[p]], repos = repos, upgrade = "never")
BiocManager::install(version = "3.21", ask = FALSE, update = FALSE)
for (p in names(bioc)) if (!have(p, bioc[[p]])) BiocManager::install(p, ask = FALSE, update = FALSE)
gh_ver <- c(TwoSampleMR = "0.7.0", ieugwasr = "1.1.0.9000", plink2R = "1.1")
for (p in names(github)) if (!have(p, gh_ver[[p]])) remotes::install_github(github[[p]], upgrade = "never")

want <- c(cran, bioc, gh_ver)
got  <- sapply(names(want), function(p) tryCatch(as.character(packageVersion(p)), error = function(e) "missing"))
bad  <- got != want
cat(sprintf("%-14s recorded %-11s installed %s\n", names(want), want, got)[bad], sep = "")
cat(if (any(bad)) "Some versions differ from the manuscript analysis (see above).\n" else "All R packages match the recorded versions.\n")
