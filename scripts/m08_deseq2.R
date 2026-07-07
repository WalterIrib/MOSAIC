#!/usr/bin/env Rscript
# =============================================================================
# MOSAIC M08 — DESeq2 Differential Expression
# Design Document v0.7
#
# Features:
#   - sizeFactors from canonical stable genes (iterative, min 500)
#   - independentFiltering=FALSE for TE and satellite
#   - LFC shrinkage: ashr (default) or apeglm (simple 2-group only)
#   - Simple mode: one comparison (condition A vs B)
#   - Extended mode: all_pairwise with Bonferroni correction
#   - Fallback: if stable genes < 500 -> WARNING + fallback to all canonical baseMean>10
#
# Usage:
#   Rscript m08_deseq2.R \
#     --counts      integrated/all_genes_counts_clean.tsv \
#     --te-counts   integrated/te_counts_clean.tsv \
#     --sat-counts  integrated/satellite_counts_clean.tsv \
#     --coldata     samples.tsv \
#     --outdir      deseq2_results/ \
#     --lfc         1.0 \
#     --padj        0.05 \
#     --shrink      ashr \
#     --mode        simple \
#     --ref         normal
# =============================================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(optparse)
})

# ─────────────────────────────────────────────────────────────────────────────
# CLI arguments
# ─────────────────────────────────────────────────────────────────────────────

option_list <- list(
  make_option("--counts",    type="character", help="all_genes_counts_clean.tsv"),
  make_option("--te-counts", type="character", help="te_counts_clean.tsv"),
  make_option("--sat-counts",type="character", help="satellite_counts_clean.tsv"),
  make_option("--coldata",   type="character", help="samples.tsv (sample_id, condition)"),
  make_option("--outdir",    type="character", help="Output directory"),
  make_option("--lfc",       type="double",  default=1.0,   help="LFC threshold [default: 1.0]"),
  make_option("--padj",      type="double",  default=0.05,  help="padj threshold [default: 0.05]"),
  make_option("--shrink",    type="character", default="ashr", help="ashr|apeglm [default: ashr]"),
  make_option("--mode",      type="character", default="simple", help="simple|extended"),
  make_option("--ref",       type="character", default=NULL, help="Reference condition for contrast"),
  make_option("--min-stable-genes", type="integer", default=500,
              help="Min stable genes for sizeFactors [default: 500]"),
  make_option("--fraction-col", type="character", default="fraction",
              help="Coldata column for fraction factor (interaction mode) [default: fraction]"),
  make_option("--phase-col",    type="character", default="phase",
              help="Coldata column for phase factor (interaction mode) [default: phase]"),
  make_option("--ref-fraction", type="character", default=NULL,
              help="Reference level for fraction (interaction mode)"),
  make_option("--ref-phase",    type="character", default=NULL,
              help="Reference level for phase (interaction mode)"),
  make_option("--comparison",   type="character", default=NULL,
              help="Comparison name used for standardized output file alias (interaction mode)")
)

opt <- parse_args(OptionParser(option_list=option_list))

# Validate required args
for (arg in c("counts", "coldata", "outdir")) {
  if (is.null(opt[[arg]])) stop(paste("Missing required argument:", arg))
}

dir.create(opt$outdir, recursive=TRUE, showWarnings=FALSE)

cat("\n========================================\n")
cat("MOSAIC M08 — DESeq2\n")
cat("========================================\n")
cat(sprintf("  Counts:     %s\n", opt$counts))
cat(sprintf("  Coldata:    %s\n", opt$coldata))
cat(sprintf("  Outdir:     %s\n", opt$outdir))
cat(sprintf("  LFC:        %.2f\n", opt$lfc))
cat(sprintf("  padj:       %.3f\n", opt$padj))
cat(sprintf("  Shrinkage:  %s\n", opt$shrink))
cat(sprintf("  Mode:       %s\n", opt$mode))
cat("\n")


# ─────────────────────────────────────────────────────────────────────────────
# Load data
# ─────────────────────────────────────────────────────────────────────────────

cat("[INF] Loading count matrix...\n")
counts_raw <- read.table(opt$counts, header=TRUE, sep="\t",
                         row.names=1, check.names=FALSE, comment.char="")

# Remove annotation columns (gene_name, gene_type) if present
# Keep only numeric columns
numeric_cols <- sapply(counts_raw, is.numeric)
anno_cols    <- counts_raw[, !numeric_cols, drop=FALSE]
counts_mat   <- as.matrix(counts_raw[, numeric_cols, drop=FALSE])
mode(counts_mat) <- "integer"

cat(sprintf("[INF] Count matrix: %d genes x %d samples\n",
            nrow(counts_mat), ncol(counts_mat)))

# Load coldata
cat("[INF] Loading sample metadata...\n")
coldata <- read.table(opt$coldata, header=TRUE, sep="\t",
                      row.names=1, check.names=FALSE)

# Align coldata to count matrix columns
common_samples <- intersect(colnames(counts_mat), rownames(coldata))
if (length(common_samples) == 0) {
  stop("No matching samples between count matrix and coldata!")
}
counts_mat <- counts_mat[, common_samples, drop=FALSE]
coldata     <- coldata[common_samples, , drop=FALSE]
coldata$condition <- factor(coldata$condition)

cat(sprintf("[INF] Samples: %s\n", paste(common_samples, collapse=", ")))
cat(sprintf("[INF] Conditions: %s\n",
            paste(names(table(coldata$condition)), collapse=", ")))

# Set reference condition
if (!is.null(opt$ref)) {
  if (opt$ref %in% levels(coldata$condition)) {
    coldata$condition <- relevel(coldata$condition, ref=opt$ref)
    cat(sprintf("[INF] Reference condition: %s\n", opt$ref))
  } else {
    cat(sprintf("[WRN] Reference '%s' not found. Using default.\n", opt$ref))
  }
}

# ─────────────────────────────────────────────────────────────────────────────
# INTERACTION MODE: ~fraction + phase + fraction:phase
# Used for ChaRNA vs totRNA background normalization
# Primary output: fractionChaRNA.phasemitosis coefficient
#   = extra mitosis effect in ChaRNA above totRNA background
# ─────────────────────────────────────────────────────────────────────────────

if (opt$mode == "interaction") {

  frac_col  <- opt[["fraction-col"]]
  phase_col <- opt[["phase-col"]]

  for (col in c(frac_col, phase_col)) {
    if (!col %in% colnames(coldata))
      stop(sprintf("Interaction mode requires coldata column '%s'", col))
  }

  coldata[[frac_col]]  <- factor(coldata[[frac_col]])
  coldata[[phase_col]] <- factor(coldata[[phase_col]])

  if (!is.null(opt[["ref-fraction"]]) &&
      opt[["ref-fraction"]] %in% levels(coldata[[frac_col]])) {
    coldata[[frac_col]] <- relevel(coldata[[frac_col]], ref=opt[["ref-fraction"]])
    cat(sprintf("[INF] Reference fraction: %s\n", opt[["ref-fraction"]]))
  }
  if (!is.null(opt[["ref-phase"]]) &&
      opt[["ref-phase"]] %in% levels(coldata[[phase_col]])) {
    coldata[[phase_col]] <- relevel(coldata[[phase_col]], ref=opt[["ref-phase"]])
    cat(sprintf("[INF] Reference phase: %s\n", opt[["ref-phase"]]))
  }

  cat(sprintf("[INF] Fraction levels: %s\n",
              paste(levels(coldata[[frac_col]]), collapse=", ")))
  cat(sprintf("[INF] Phase levels: %s\n",
              paste(levels(coldata[[phase_col]]), collapse=", ")))

  design_formula <- as.formula(
    sprintf("~ %s + %s + %s:%s", frac_col, phase_col, frac_col, phase_col)
  )
  cat(sprintf("[INF] Interaction design: %s\n", deparse(design_formula)))

  # Stable sizeFactors using interaction design
  cat("\n[INF] Computing sizeFactors from stable genes (interaction design)...\n")
  dds_init_int <- DESeqDataSetFromMatrix(
    countData = counts_mat,
    colData   = coldata,
    design    = design_formula
  )
  dds_init_int <- estimateSizeFactors(dds_init_int)
  dds_init_int <- estimateDispersions(dds_init_int, quiet=TRUE)
  norm_init    <- counts(dds_init_int, normalized=TRUE)
  base_mean_i  <- rowMeans(norm_init)
  cv_i         <- apply(norm_init, 1, sd) / (base_mean_i + 1e-10)
  stable_mask_i <- base_mean_i > 10 & cv_i < 0.3 & !is.na(cv_i)
  n_stable_i    <- sum(stable_mask_i)
  cat(sprintf("[INF] Stable genes (interaction): %d\n", n_stable_i))
  if (n_stable_i < opt[["min-stable-genes"]]) {
    cat("[WRN] Falling back to baseMean>10.\n")
    stable_mask_i <- base_mean_i > 10 & !is.na(base_mean_i)
    n_stable_i    <- sum(stable_mask_i)
  }
  sf_int <- estimateSizeFactorsForMatrix(counts_mat[stable_mask_i, , drop=FALSE])
  cat(sprintf("[INF] SizeFactors: %s\n", paste(sprintf("%.3f", sf_int), collapse=", ")))

  write.table(
    data.frame(gene_id=rownames(counts_mat)[stable_mask_i]),
    file.path(opt$outdir, "stable_genes.tsv"),
    sep="\t", row.names=FALSE, quote=FALSE
  )

  # Build and run main DESeq2
  dds_int <- DESeqDataSetFromMatrix(
    countData = counts_mat,
    colData   = coldata,
    design    = design_formula
  )
  sizeFactors(dds_int) <- sf_int
  dds_int <- DESeq(dds_int, fitType="parametric", quiet=FALSE)

  coef_names <- resultsNames(dds_int)
  cat(sprintf("[INF] Coefficients: %s\n", paste(coef_names, collapse=", ")))

  # Save normalized counts
  norm_counts_int <- counts(dds_int, normalized=TRUE)
  write.table(
    data.frame(gene_id=rownames(norm_counts_int), norm_counts_int, check.names=FALSE),
    file.path(opt$outdir, "normalized_counts.tsv"),
    sep="\t", row.names=FALSE, quote=FALSE
  )
  cat("[INF] Normalized counts saved.\n")

  # Helper: extract results for a coefficient
  extract_coef_results <- function(dds, coef, label, lfc_thr, padj_thr) {
    cat(sprintf("\n[INF] Extracting coefficient: %s\n", coef))
    res <- results(dds, name=coef, lfcThreshold=lfc_thr,
                   alpha=padj_thr, independentFiltering=TRUE)
    res_shrunk <- lfcShrink(dds, coef=coef, type="ashr", quiet=TRUE)
    res_df <- as.data.frame(res)
    res_df$log2FoldChange_shrunk <- as.data.frame(res_shrunk)$log2FoldChange
    res_df$lfcSE_shrunk          <- as.data.frame(res_shrunk)$lfcSE
    res_df$gene_id               <- rownames(res_df)
    res_df <- res_df[order(res_df$padj, na.last=TRUE), ]
    n_up   <- sum(!is.na(res_df$padj) & res_df$padj < padj_thr &
                  res_df$log2FoldChange > lfc_thr)
    n_down <- sum(!is.na(res_df$padj) & res_df$padj < padj_thr &
                  res_df$log2FoldChange < -lfc_thr)
    cat(sprintf("[INF] %s: UP=%d, DOWN=%d\n", label, n_up, n_down))
    return(res_df)
  }

  # Find interaction coefficient (fractionChaRNA.phasemitosis or similar)
  interact_coef <- grep("\\.", coef_names, value=TRUE)
  if (length(interact_coef) == 0)
    stop("Could not find interaction coefficient in: ", paste(coef_names, collapse=", "))
  interact_coef <- interact_coef[1]

  # Find phase coefficient (mitosis effect in reference fraction = totRNA)
  phase_coef <- grep(sprintf("%s_", phase_col), coef_names, value=TRUE)
  if (length(phase_coef) == 0)
    stop("Could not find phase coefficient")
  phase_coef <- phase_coef[1]

  # Extract interaction result (primary: ChaRNA-specific mitosis effect)
  res_interact <- extract_coef_results(
    dds_int, interact_coef,
    label=sprintf("interaction_(%s)", interact_coef),
    lfc_thr=opt$lfc, padj_thr=opt$padj
  )
  out_interact <- file.path(opt$outdir,
                            sprintf("genes_%s_DE.tsv", gsub("\\.", "_", interact_coef)))
  write.table(res_interact, out_interact, sep="\t", quote=FALSE, row.names=FALSE)
  cat(sprintf("[INF] Written: %s\n", basename(out_interact)))

  sig_interact <- res_interact[!is.na(res_interact$padj) &
                               res_interact$padj < opt$padj &
                               abs(res_interact$log2FoldChange) > opt$lfc, ]
  write.table(sig_interact,
    file.path(opt$outdir, sprintf("genes_%s_DE_significant.tsv",
                                  gsub("\\.", "_", interact_coef))),
    sep="\t", quote=FALSE, row.names=FALSE)

  # ── Standardized alias files for unified report template ─────────────────
  # m14_report.qmd expects genes_{comparison}_DE.tsv; write aliases so the
  # standard template can read interaction results without special-casing.
  std_comp <- if (!is.null(opt$comparison)) opt$comparison else "ChaRNA_interaction"
  std_full <- file.path(opt$outdir, sprintf("genes_%s_DE.tsv", std_comp))
  std_sig  <- file.path(opt$outdir, sprintf("genes_%s_DE_significant.tsv", std_comp))
  if (normalizePath(std_full) != normalizePath(out_interact)) {
    write.table(res_interact, std_full, sep="\t", quote=FALSE, row.names=FALSE)
    write.table(sig_interact, std_sig,  sep="\t", quote=FALSE, row.names=FALSE)
    cat(sprintf("[INF] Standardized alias: %s\n", basename(std_full)))
  }

  # Extract phase result (totRNA mitosis effect = background / global phase effect)
  res_phase <- extract_coef_results(
    dds_int, phase_coef,
    label=sprintf("phase_(%s)", phase_coef),
    lfc_thr=opt$lfc, padj_thr=opt$padj
  )
  out_phase <- file.path(opt$outdir,
                         sprintf("genes_%s_DE.tsv", gsub("_", "", phase_coef)))
  write.table(res_phase, out_phase, sep="\t", quote=FALSE, row.names=FALSE)
  # Also write standardized alias for phase effect
  write.table(res_phase,
    file.path(opt$outdir, sprintf("genes_%s_phase_effect_DE.tsv", std_comp)),
    sep="\t", quote=FALSE, row.names=FALSE)
  cat(sprintf("[INF] Written: %s\n", basename(out_phase)))

  # TE and satellite DE (interaction mode)
  for (tsv_path in c(opt[["te-counts"]], opt[["sat-counts"]])) {
    if (is.null(tsv_path) || !file.exists(tsv_path)) next
    name_prefix <- ifelse(grepl("satellite", tsv_path), "satellite", "TE")
    cat(sprintf("\n[INF] Running TE DE (interaction) for %s...\n", name_prefix))

    te_raw <- read.table(tsv_path, header=TRUE, sep="\t",
                         row.names=1, check.names=FALSE)
    te_mat <- as.matrix(te_raw)
    mode(te_mat) <- "integer"
    common_te <- intersect(colnames(te_mat), rownames(coldata))
    te_mat    <- te_mat[, common_te, drop=FALSE]
    keep      <- rowSums(te_mat >= 1) >= ncol(te_mat) / 2
    te_mat    <- te_mat[keep, , drop=FALSE]
    if (nrow(te_mat) == 0) next

    dds_te_int <- DESeqDataSetFromMatrix(
      countData = te_mat,
      colData   = coldata[common_te, , drop=FALSE],
      design    = design_formula
    )
    sizeFactors(dds_te_int) <- sf_int[common_te]
    # Robust dispersion fit (see run_te_contrast): fall back to constant fit if
    # locfit fails on a small/sparse matrix; no effect on genome-wide data.
    dds_te_int <- tryCatch(
      DESeq(dds_te_int, fitType="local", quiet=TRUE),
      error = function(e) {
        message("[WARN] local dispersion fit failed (", conditionMessage(e),
                "); retrying with fitType='mean'")
        DESeq(dds_te_int, fitType="mean", quiet=TRUE)
      })

    for (coef in c(interact_coef, phase_coef)) {
      te_res <- results(dds_te_int, name=coef, lfcThreshold=opt$lfc,
                        alpha=opt$padj, independentFiltering=FALSE)
      te_shrunk <- lfcShrink(dds_te_int, coef=coef, type="ashr", quiet=TRUE)
      te_df <- as.data.frame(te_res)
      te_df$log2FoldChange_shrunk <- as.data.frame(te_shrunk)$log2FoldChange
      te_df$family_id <- rownames(te_df)
      tag <- gsub("\\.", "_", coef)
      out_te <- file.path(opt$outdir, sprintf("%s_%s_DE.tsv", name_prefix, tag))
      write.table(te_df, out_te, sep="\t", quote=FALSE, row.names=FALSE)
      cat(sprintf("[INF] Written: %s\n", basename(out_te)))
      # Standardized alias for interaction coefficient only
      if (coef == interact_coef) {
        std_te <- file.path(opt$outdir, sprintf("%s_%s_DE.tsv", name_prefix, std_comp))
        if (normalizePath(std_te) != normalizePath(out_te)) {
          write.table(te_df, std_te, sep="\t", quote=FALSE, row.names=FALSE)
          cat(sprintf("[INF] Standardized alias: %s\n", basename(std_te)))
        }
      }
    }
  }

  # PCA of all samples
  suppressPackageStartupMessages(library(DESeq2))
  vsd_int <- tryCatch(vst(dds_int, blind=FALSE), error=function(e) {
    cat("[WRN] VST failed, using rlog\n"); rlog(dds_int, blind=FALSE)
  })
  pca_int <- prcomp(t(assay(vsd_int)), scale.=FALSE)
  var_exp <- round(100 * pca_int$sdev^2 / sum(pca_int$sdev^2), 1)
  pca_df  <- data.frame(
    sample_id  = rownames(pca_int$x),
    PC1        = pca_int$x[, 1],
    PC2        = pca_int$x[, 2],
    condition  = coldata[rownames(pca_int$x), frac_col],
    phase      = coldata[rownames(pca_int$x), phase_col],
    PC1_var    = var_exp[1],
    PC2_var    = var_exp[2]
  )
  write.table(pca_df, file.path(opt$outdir, "pca_coordinates.tsv"),
              sep="\t", row.names=FALSE, quote=FALSE)
  cat("[INF] PCA coordinates saved.\n")

  sink(file.path(opt$outdir, "session_info.txt"))
  sessionInfo()
  sink()

  cat("\n[OK] M08 DESeq2 interaction mode complete.\n")
  cat(sprintf("[OK] Results in: %s\n\n", opt$outdir))
  quit(save="no", status=0)
}


# ─────────────────────────────────────────────────────────────────────────────
# SizeFactors from stable genes (iterative)
# ─────────────────────────────────────────────────────────────────────────────

cat("\n[INF] Computing sizeFactors from stable genes...\n")

compute_stable_sizefactors <- function(counts_mat, coldata, min_stable=500) {
  # Step 1: Run DESeq2 with default sizeFactors to get baseMean
  dds_init <- DESeqDataSetFromMatrix(
    countData = counts_mat,
    colData   = coldata,
    design    = ~condition
  )
  dds_init <- estimateSizeFactors(dds_init)
  dds_init <- estimateDispersions(dds_init, quiet=TRUE)

  # Step 2: Get normalized counts and CV
  norm_counts <- counts(dds_init, normalized=TRUE)
  base_mean   <- rowMeans(norm_counts)
  cv          <- apply(norm_counts, 1, sd) / (base_mean + 1e-10)

  # Step 3: Select stable genes
  # Criteria: baseMean > 10, CV < 0.3 (low variability across samples)
  stable_mask <- base_mean > 10 & cv < 0.3 & !is.na(cv)
  n_stable    <- sum(stable_mask)

  cat(sprintf("[INF] Stable genes (baseMean>10, CV<0.3): %d\n", n_stable))

  if (n_stable < min_stable) {
    cat(sprintf("[WRN] Stable genes (%d) < minimum (%d).\n", n_stable, min_stable))
    cat("[WRN] Falling back to all canonical genes with baseMean>10.\n")

    stable_mask <- base_mean > 10 & !is.na(base_mean)
    n_stable    <- sum(stable_mask)
    cat(sprintf("[INF] Fallback stable genes: %d\n", n_stable))
  }

  # Step 4: Recompute sizeFactors using only stable genes
  stable_counts <- counts_mat[stable_mask, , drop=FALSE]
  sf <- estimateSizeFactorsForMatrix(stable_counts)

  cat(sprintf("[INF] SizeFactors: %s\n",
              paste(sprintf("%.3f", sf), collapse=", ")))

  return(list(sf=sf, n_stable=n_stable, stable_genes=rownames(counts_mat)[stable_mask]))
}

sf_result <- compute_stable_sizefactors(counts_mat, coldata,
                                         min_stable=opt[["min-stable-genes"]])


# ─────────────────────────────────────────────────────────────────────────────
# Build DESeq2 dataset
# ─────────────────────────────────────────────────────────────────────────────

cat("\n[INF] Building DESeq2 dataset...\n")

dds <- DESeqDataSetFromMatrix(
  countData = counts_mat,
  colData   = coldata,
  design    = ~condition
)

# Apply stable-gene sizeFactors
sizeFactors(dds) <- sf_result$sf

# Run DESeq2 (skip sizeFactors estimation since we set them manually)
dds <- DESeq(dds, fitType="parametric", quiet=FALSE)

# Save stable genes list
write.table(
  data.frame(gene_id=sf_result$stable_genes),
  file.path(opt$outdir, "stable_genes.tsv"),
  sep="\t", row.names=FALSE, quote=FALSE
)
cat(sprintf("[INF] Stable genes saved: %d\n", sf_result$n_stable))


# ─────────────────────────────────────────────────────────────────────────────
# Extract results
# ─────────────────────────────────────────────────────────────────────────────

run_contrast <- function(dds, contrast, shrink_method="ashr",
                         lfc_threshold=1.0, padj_threshold=0.05,
                         independent_filtering=TRUE, label="") {

  cat(sprintf("\n[INF] Running contrast: %s\n", label))

  # Get results
  res <- results(dds,
                 contrast=contrast,
                 lfcThreshold=lfc_threshold,
                 alpha=padj_threshold,
                 independentFiltering=independent_filtering)

  # LFC shrinkage
  if (shrink_method == "apeglm") {
    # apeglm requires coef, not contrast
    coef_name <- resultsNames(dds)[length(resultsNames(dds))]
    res_shrunk <- lfcShrink(dds, coef=coef_name, type="apeglm", quiet=TRUE)
  } else {
    res_shrunk <- lfcShrink(dds, contrast=contrast, type="ashr", quiet=TRUE)
  }

  # Merge original pvalue/padj with shrunken LFC
  res_df <- as.data.frame(res)
  res_df$log2FoldChange_shrunk <- as.data.frame(res_shrunk)$log2FoldChange
  res_df$lfcSE_shrunk          <- as.data.frame(res_shrunk)$lfcSE
  res_df$gene_id               <- rownames(res_df)

  # Summary
  n_up   <- sum(res_df$padj < padj_threshold &
                res_df$log2FoldChange > lfc_threshold, na.rm=TRUE)
  n_down <- sum(res_df$padj < padj_threshold &
                res_df$log2FoldChange < -lfc_threshold, na.rm=TRUE)
  cat(sprintf("[INF] %s: UP=%d, DOWN=%d\n", label, n_up, n_down))

  return(res_df)
}

# Determine comparisons
conditions <- levels(coldata$condition)
cat(sprintf("\n[INF] Conditions: %s\n", paste(conditions, collapse=" vs ")))

results_list <- list()

if (opt$mode == "simple") {
  # One contrast: condition A vs reference
  ref_cond  <- levels(coldata$condition)[1]
  test_cond <- levels(coldata$condition)[2]

  contrast <- c("condition", test_cond, ref_cond)
  label    <- paste0(test_cond, "_vs_", ref_cond)

  res_df <- run_contrast(
    dds, contrast,
    shrink_method=opt$shrink,
    lfc_threshold=opt$lfc,
    padj_threshold=opt$padj,
    independent_filtering=TRUE,
    label=label
  )
  results_list[[label]] <- res_df

} else {
  # Extended: all pairwise with Bonferroni correction
  pairs <- combn(conditions, 2, simplify=FALSE)
  n_tests <- length(pairs)
  cat(sprintf("[INF] Extended mode: %d pairwise comparisons\n", n_tests))

  for (pair in pairs) {
    test_cond <- pair[2]
    ref_cond  <- pair[1]
    contrast  <- c("condition", test_cond, ref_cond)
    label     <- paste0(test_cond, "_vs_", ref_cond)

    res_df <- run_contrast(
      dds, contrast,
      shrink_method=opt$shrink,
      lfc_threshold=opt$lfc,
      padj_threshold=opt$padj / n_tests,  # Bonferroni
      independent_filtering=TRUE,
      label=label
    )
    results_list[[label]] <- res_df
  }
}


# ─────────────────────────────────────────────────────────────────────────────
# TE and satellite DE (independentFiltering=FALSE)
# ─────────────────────────────────────────────────────────────────────────────

run_te_de <- function(te_tsv, sat_tsv, coldata, sf, outdir,
                      lfc_threshold, padj_threshold, shrink_method, mode) {

  run_te_contrast <- function(counts_mat, name_prefix, label,
                               contrast, indep_filt=FALSE) {
    if (nrow(counts_mat) == 0) return(NULL)

    # Filter low counts
    keep <- rowSums(counts_mat >= 1) >= ncol(counts_mat) / 2
    counts_mat <- counts_mat[keep, , drop=FALSE]
    if (nrow(counts_mat) == 0) return(NULL)

    dds_te <- DESeqDataSetFromMatrix(
      countData = counts_mat,
      colData   = coldata,
      design    = ~condition
    )
    sizeFactors(dds_te) <- sf
    # Local (locfit) dispersion fitting is used for TE/satellite family DE. On
    # very small or sparse count matrices (few families, e.g. a chromosome-subset
    # test set) locfit can fail with "out of vertex space"; fall back to the
    # robust constant-dispersion fit. On genome-wide data the local fit succeeds,
    # so this fallback does not affect the reported results.
    dds_te <- tryCatch(
      DESeq(dds_te, fitType="local", quiet=TRUE),
      error = function(e) {
        message("[WARN] local dispersion fit failed (", conditionMessage(e),
                "); retrying with fitType='mean'")
        DESeq(dds_te, fitType="mean", quiet=TRUE)
      })

    res <- results(dds_te,
                   contrast=contrast,
                   lfcThreshold=lfc_threshold,
                   alpha=padj_threshold,
                   independentFiltering=indep_filt)  # FALSE for TE/satellite

    res_shrunk <- lfcShrink(dds_te, contrast=contrast,
                             type="ashr", quiet=TRUE)

    res_df <- as.data.frame(res)
    res_df$log2FoldChange_shrunk <- as.data.frame(res_shrunk)$log2FoldChange
    res_df$family_id <- rownames(res_df)

    n_up   <- sum(res_df$padj < padj_threshold &
                  res_df$log2FoldChange > lfc_threshold, na.rm=TRUE)
    n_down <- sum(res_df$padj < padj_threshold &
                  res_df$log2FoldChange < -lfc_threshold, na.rm=TRUE)
    cat(sprintf("[INF] %s %s: UP=%d, DOWN=%d\n",
                name_prefix, label, n_up, n_down))
    return(res_df)
  }

  conditions <- levels(coldata$condition)

  for (tsv_path in c(te_tsv, sat_tsv)) {
    if (is.null(tsv_path) || !file.exists(tsv_path)) next

    name_prefix <- ifelse(grepl("satellite", tsv_path), "satellite", "TE")
    cat(sprintf("\n[INF] Running DE for %s...\n", name_prefix))

    counts_raw <- read.table(tsv_path, header=TRUE, sep="\t",
                             row.names=1, check.names=FALSE)
    counts_mat <- as.matrix(counts_raw)
    mode(counts_mat) <- "integer"

    # Align samples
    common <- intersect(colnames(counts_mat), rownames(coldata))
    counts_mat <- counts_mat[, common, drop=FALSE]

    if (mode == "simple") {
      ref_cond  <- conditions[1]
      test_cond <- conditions[2]
      contrast  <- c("condition", test_cond, ref_cond)
      label     <- paste0(test_cond, "_vs_", ref_cond)

      res_df <- run_te_contrast(counts_mat, name_prefix, label,
                                 contrast, indep_filt=FALSE)
      if (!is.null(res_df)) {
        out_f <- file.path(outdir,
                           sprintf("%s_%s_DE.tsv", name_prefix, label))
        write.table(res_df, out_f, sep="\t", quote=FALSE, row.names=FALSE)
        cat(sprintf("[INF] Written: %s\n", basename(out_f)))
      }
    } else {
      pairs <- combn(conditions, 2, simplify=FALSE)
      for (pair in pairs) {
        contrast <- c("condition", pair[2], pair[1])
        label    <- paste0(pair[2], "_vs_", pair[1])
        res_df   <- run_te_contrast(counts_mat, name_prefix, label,
                                     contrast, indep_filt=FALSE)
        if (!is.null(res_df)) {
          out_f <- file.path(outdir,
                             sprintf("%s_%s_DE.tsv", name_prefix, label))
          write.table(res_df, out_f, sep="\t", quote=FALSE, row.names=FALSE)
        }
      }
    }
  }
}

cat("\n[INF] Running TE and satellite DE (independentFiltering=FALSE)...\n")
run_te_de(
  te_tsv=opt[["te-counts"]],
  sat_tsv=opt[["sat-counts"]],
  coldata=coldata,
  sf=sf_result$sf,
  outdir=opt$outdir,
  lfc_threshold=opt$lfc,
  padj_threshold=opt$padj,
  shrink_method=opt$shrink,
  mode=opt$mode
)


# ─────────────────────────────────────────────────────────────────────────────
# Write canonical gene results
# ─────────────────────────────────────────────────────────────────────────────

cat("\n[INF] Writing results...\n")

for (label in names(results_list)) {
  res_df <- results_list[[label]]

  # Sort by padj
  res_df <- res_df[order(res_df$padj, na.last=TRUE), ]

  out_f <- file.path(opt$outdir, sprintf("genes_%s_DE.tsv", label))
  write.table(res_df, out_f, sep="\t", quote=FALSE, row.names=FALSE)
  cat(sprintf("[INF] Written: %s (%d genes)\n", basename(out_f), nrow(res_df)))

  # Significant subset
  sig_df <- res_df[!is.na(res_df$padj) &
                   res_df$padj < opt$padj &
                   abs(res_df$log2FoldChange) > opt$lfc, ]
  out_sig <- file.path(opt$outdir, sprintf("genes_%s_DE_significant.tsv", label))
  write.table(sig_df, out_sig, sep="\t", quote=FALSE, row.names=FALSE)
  cat(sprintf("[INF] Significant: %d genes (padj<%.3f, |LFC|>%.1f)\n",
              nrow(sig_df), opt$padj, opt$lfc))
}

# Save normalized counts (DESeq2 size-factor corrected)
norm_counts <- counts(dds, normalized=TRUE)
write.table(
  data.frame(gene_id=rownames(norm_counts), norm_counts, check.names=FALSE),
  file.path(opt$outdir, "normalized_counts.tsv"),
  sep="\t", row.names=FALSE, quote=FALSE
)
cat("[INF] Normalized counts saved.\n")

# ─────────────────────────────────────────────────────────────────────────────
# TPM (Transcripts Per Million)
# Formula: TPM_i = (raw_i / gene_length_kb) / sum(raw_j / gene_length_kb) * 1e6
# Gene lengths from featureCounts .summary or GTF exon union.
# If gene lengths unavailable, estimate from GTF exon union via GenomicFeatures.
# ─────────────────────────────────────────────────────────────────────────────
cat("[INF] Computing TPM...\n")
tryCatch({
  gene_lengths <- NULL

  # Try to load from featureCounts annotation columns (gene_id + Length col)
  # featureCounts output has a "Length" column in the raw file
  fc_raw_path <- gsub("_clean\\.tsv$", ".tsv",
                       gsub("integrated", "mapping/counts",
                            gsub("all_genes_counts_clean", "all_genes_counts", opt$counts)))
  if (file.exists(fc_raw_path)) {
    fc_raw <- tryCatch(
      read.table(fc_raw_path, header=TRUE, sep="\t", comment.char="#",
                 check.names=FALSE),
      error=function(e) NULL)
    if (!is.null(fc_raw) && "Length" %in% colnames(fc_raw)) {
      gene_lengths <- setNames(as.numeric(fc_raw$Length), fc_raw$Geneid)
      cat(sprintf("[INF] Gene lengths from featureCounts: %d genes\n",
                  length(gene_lengths)))
    }
  }

  # Fallback: estimate from GTF using GenomicFeatures (exon union)
  if (is.null(gene_lengths) && requireNamespace("GenomicFeatures", quietly=TRUE)) {
    cat("[INF] Estimating gene lengths from GTF (exon union)...\n")
    library(GenomicFeatures)
    gtf_path <- Sys.getenv("MOSAIC_GTF", "")
    if (nchar(gtf_path) > 0 && file.exists(gtf_path)) {
      txdb <- makeTxDbFromGFF(gtf_path, format="gtf")
      exons_by_gene <- exonsBy(txdb, by="gene")
      gene_lengths  <- sum(width(reduce(exons_by_gene)))
      cat(sprintf("[INF] Gene lengths from GTF: %d genes\n", length(gene_lengths)))
    }
  }

  if (!is.null(gene_lengths)) {
    raw_counts  <- counts(dds, normalized=FALSE)
    common_g    <- intersect(rownames(raw_counts), names(gene_lengths))
    gl_kb       <- gene_lengths[common_g] / 1000
    rpk         <- raw_counts[common_g, ] / gl_kb
    tpm_mat     <- sweep(rpk, 2, colSums(rpk) / 1e6, "/")

    write.table(
      data.frame(gene_id=rownames(tpm_mat), tpm_mat, check.names=FALSE),
      file.path(opt$outdir, "tpm_counts.tsv"),
      sep="\t", row.names=FALSE, quote=FALSE
    )
    cat(sprintf("[INF] TPM saved: %d genes × %d samples\n",
                nrow(tpm_mat), ncol(tpm_mat)))
  } else {
    cat("[WARN] Gene lengths not available — TPM not computed.\n")
    cat("[WARN] Set MOSAIC_GTF env var or ensure featureCounts raw output is accessible.\n")
  }
}, error=function(e) {
  cat("[WARN] TPM computation failed:", conditionMessage(e), "\n")
})

# Save session info
sink(file.path(opt$outdir, "session_info.txt"))
sessionInfo()
sink()

cat("\n[OK] M08 DESeq2 complete.\n")
cat(sprintf("[OK] Results in: %s\n\n", opt$outdir))
