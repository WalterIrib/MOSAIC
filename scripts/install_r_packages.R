#!/usr/bin/env Rscript
# Alternative to conda for the R stack: install the exact Bioconductor/CRAN
# package versions used for the manuscript into a system/user R 4.5.1 library.
# Use this only if you run MOSAIC against a system R instead of the conda r-base.
#
#   Rscript scripts/install_r_packages.R
#
# Versions verified from the environment that produced the canonical results.

want <- list(
  # CRAN
  optparse = "1.7.5",
  dplyr    = "1.2.0",
  ashr     = "2.2.63",
  WGCNA    = "1.74",
  # Bioconductor
  DESeq2          = "1.50.2",
  clusterProfiler = "4.18.4",
  STRINGdb        = "2.22.0"
)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = "https://cloud.r-project.org")

for (pkg in names(want)) {
  cur <- tryCatch(as.character(packageVersion(pkg)), error = function(e) NA)
  if (!is.na(cur) && cur == want[[pkg]]) {
    cat(sprintf("[OK]   %-16s %s already installed\n", pkg, cur)); next
  }
  cat(sprintf("[INFO] installing %s (want %s, have %s)\n",
              pkg, want[[pkg]], ifelse(is.na(cur), "none", cur)))
  BiocManager::install(pkg, update = FALSE, ask = FALSE)
}

cat("\nDone. Re-run to verify versions.\n")
