#!/usr/bin/env Rscript
# mosaic/scripts/m12_deseq2.R
# M12 — DESeq2 on TE locus-level counts (from m12_te_locus.py locus-count)
#
# Input:  m12/locus_counts.tsv   (loci x samples, all 12 samples)
#         coldata.tsv
# Output: m12/{comparison}_te_de.tsv          (all loci, DE results)
#         m12/{comparison}_te_de_sig.tsv       (significant only, padj<threshold)

suppressPackageStartupMessages({
  library(optparse)
  library(DESeq2)
  library(data.table)
  library(dplyr)
})

option_list <- list(
  make_option("--locus-counts",  type="character", help="m12/locus_counts.tsv"),
  make_option("--coldata",       type="character", help="coldata.tsv"),
  make_option("--outdir",        type="character", help="Output directory (m12/)"),
  make_option("--treatment",     type="character", help="Treatment condition"),
  make_option("--control",       type="character", help="Control condition"),
  make_option("--comparison",    type="character", help="Comparison name (e.g. mito-ChaRNA_vs_ChaRNA)"),
  make_option("--gene-counts",   type="character", default=NULL,
              help="all_genes_counts.tsv from comparison deseq2 dir (for canonical size factors)"),
  make_option("--stable-genes",  type="character", default=NULL,
              help="stable_genes.tsv from comparison deseq2 dir"),
  make_option("--padj",          type="double",    default=0.05),
  make_option("--lfc",           type="double",    default=1.0),
  make_option("--min-counts",    type="integer",   default=10L,
              help="Min total counts across all samples to keep locus [10]"),
  make_option("--threads",       type="integer",   default=4L)
)

opt <- parse_args(OptionParser(option_list=option_list))

BiocParallel::register(BiocParallel::MulticoreParam(opt$threads))

cat("[INF] M12 DESeq2 — locus-level TE differential expression\n")
cat(sprintf("[INF] Comparison: %s vs %s\n", opt$treatment, opt$control))

# ── 1. Load counts ─────────────────────────────────────────────────────────────
cat("[INF] Reading locus_counts.tsv (may take a moment)...\n")
cts_raw <- fread(opt[["locus-counts"]], sep="\t", header=TRUE,
                 data.table=FALSE, check.names=FALSE)

# Strip _Aligned.sortedByCoord.out suffix from sample columns
meta_cols <- c("locus_id","te_name","chrom","start","end","strand")
samp_cols <- setdiff(colnames(cts_raw), meta_cols)
samp_clean <- sub("_Aligned\\.sortedByCoord\\.out$", "", samp_cols)
colnames(cts_raw)[match(samp_cols, colnames(cts_raw))] <- samp_clean

cat(sprintf("[INF] Loci (total): %s | Samples: %d\n",
            format(nrow(cts_raw), big.mark=","), length(samp_clean)))

# ── 2. Load coldata ────────────────────────────────────────────────────────────
coldata <- read.table(opt$coldata, header=TRUE, sep="\t",
                      stringsAsFactors=FALSE, check.names=FALSE)

# Filter to comparison samples only
sids <- coldata$sample_id[coldata$condition %in% c(opt$treatment, opt$control)]
sids <- intersect(sids, samp_clean)

if (length(sids) < 4) {
  cat(sprintf("[ERROR] Need >=4 samples; found %d for %s vs %s\n",
              length(sids), opt$treatment, opt$control))
  quit(status=1)
}

cat(sprintf("[INF] Samples in comparison: %d (%s)\n",
            length(sids), paste(sids, collapse=", ")))

# ── 3. Build count matrix ──────────────────────────────────────────────────────
count_mat <- as.matrix(cts_raw[, sids, drop=FALSE])
rownames(count_mat) <- cts_raw$locus_id

# Filter low-count loci (across ALL 12 samples for stable estimation)
full_mat <- as.matrix(cts_raw[, intersect(samp_clean, colnames(cts_raw)), drop=FALSE])
keep      <- rowSums(full_mat) >= opt[["min-counts"]]
count_mat <- count_mat[keep, , drop=FALSE]
te_meta   <- cts_raw[keep, meta_cols, drop=FALSE]

cat(sprintf("[INF] Loci after filtering (total >= %d reads): %s\n",
            opt[["min-counts"]], format(nrow(count_mat), big.mark=",")))

# ── 4. DESeq2 ──────────────────────────────────────────────────────────────────
dds_cd <- data.frame(
  row.names = sids,
  condition = factor(coldata$condition[match(sids, coldata$sample_id)],
                     levels=c(opt$control, opt$treatment))
)

dds <- DESeqDataSetFromMatrix(
  countData = count_mat,
  colData   = dds_cd,
  design    = ~ condition
)

# ── Size factors from canonical stable genes (same as pipeline DESeq2) ─────────
# Must use gene-derived size factors so TE expression is comparable to canonical
# gene results and not relative to other TE loci.
gene_sf <- NULL

if (!is.null(opt[["gene-counts"]]) && file.exists(opt[["gene-counts"]])) {
  cat("[INF] Loading canonical gene counts for size factor estimation...\n")
  gc_raw <- fread(opt[["gene-counts"]], sep="\t", header=TRUE, data.table=FALSE,
                  check.names=FALSE)
  # Column 1 = Geneid, remaining = samples
  gc_mat <- as.matrix(gc_raw[, sids, drop=FALSE])
  rownames(gc_mat) <- gc_raw[[1]]

  # Use stable genes if provided; otherwise all genes with baseMean > 10
  if (!is.null(opt[["stable-genes"]]) && file.exists(opt[["stable-genes"]])) {
    stable <- read.table(opt[["stable-genes"]], header=TRUE, sep="\t",
                         stringsAsFactors=FALSE)[[1]]
    # Strip version from IDs if needed
    stable_base <- sub("\\..*", "", stable)
    gids_base   <- sub("\\..*", "", rownames(gc_mat))
    use_genes   <- rownames(gc_mat)[gids_base %in% stable_base]
    cat(sprintf("[INF] Using %d stable genes for size factor estimation\n",
                length(use_genes)))
  } else {
    # Fallback: genes with baseMean > 10
    use_genes <- rownames(gc_mat)[rowMeans(gc_mat) > 10]
    cat(sprintf("[INF] No stable_genes.tsv — using %d genes (baseMean>10) for size factors\n",
                length(use_genes)))
  }

  if (length(use_genes) >= 50) {
    # Estimate size factors from canonical genes using geometric mean method
    gc_sub <- gc_mat[use_genes, , drop=FALSE]
    log_gc  <- log(gc_sub)
    log_gc[is.infinite(log_gc)] <- NA
    geom_mean <- rowMeans(log_gc, na.rm=TRUE)
    ratios     <- sweep(log_gc, 1, geom_mean, "-")
    gene_sf    <- exp(apply(ratios, 2, median, na.rm=TRUE))
    cat(sprintf("[INF] Canonical gene size factors: %s\n",
                paste(names(gene_sf), "=", round(gene_sf, 3), collapse="; ")))
  } else {
    cat("[WARN] Too few stable genes — falling back to TE-derived size factors\n")
  }
}

if (!is.null(gene_sf) && all(is.finite(gene_sf)) && all(gene_sf > 0)) {
  sizeFactors(dds) <- gene_sf[sids]
  cat("[INF] Applied canonical gene size factors to TE DESeq2 object\n")
} else {
  cat("[INF] Estimating size factors from TE loci (no gene counts provided)\n")
  dds <- estimateSizeFactors(dds)
}
cat(sprintf("[INF] Size factors: %s\n",
            paste(round(sizeFactors(dds), 3), collapse=", ")))

dds <- DESeq(dds, quiet=TRUE, parallel=FALSE)

# ── 5. Results + shrinkage ────────────────────────────────────────────────────
res <- results(dds,
               contrast=c("condition", opt$treatment, opt$control),
               independentFiltering=FALSE)

res_shrunk <- tryCatch(
  lfcShrink(dds, contrast=c("condition", opt$treatment, opt$control),
            type="ashr", quiet=TRUE),
  error=function(e) {
    cat("[WARN] lfcShrink failed:", conditionMessage(e), "\n")
    NULL
  }
)

# ── 6. Build output table ──────────────────────────────────────────────────────
out_df <- data.frame(
  locus_id   = rownames(res),
  te_name    = te_meta$te_name[match(rownames(res), te_meta$locus_id)],
  baseMean   = round(res$baseMean, 2),
  log2FoldChange = round(res$log2FoldChange, 4),
  log2FoldChange_shrunk = if (!is.null(res_shrunk))
    round(res_shrunk$log2FoldChange, 4) else NA_real_,
  pvalue = res$pvalue,
  padj   = res$padj,
  stringsAsFactors = FALSE
)

# Parse te_name: family|class|subtype or family#class/subtype
out_df$te_family <- sub("\\|.*", "", sub("#.*", "", out_df$te_name))
out_df$te_class  <- {
  # try pipe format: TAR1|Satellite|telo → Satellite
  cl <- sub("^[^|]+\\|([^|]+)\\|.*", "\\1", out_df$te_name)
  # if no pipe, try hash format: L1PA8#LINE/L1 → LINE
  cl[!grepl("\\|", out_df$te_name)] <-
    sub("^[^#]+#([^/]+)/.*", "\\1", out_df$te_name[!grepl("\\|", out_df$te_name)])
  cl
}

# gene_id alias for report compatibility
out_df$gene_id <- out_df$locus_id

out_df <- out_df[order(out_df$padj, na.last=TRUE), ]

dir.create(opt$outdir, recursive=TRUE, showWarnings=FALSE)

comp <- opt$comparison %||% sprintf("%s_vs_%s", opt$treatment, opt$control)
`%||%` <- function(a,b) if (!is.null(a) && !is.na(a)) a else b

# All results
write.table(out_df,
            file.path(opt$outdir, sprintf("%s_te_de.tsv", comp)),
            sep="\t", quote=FALSE, row.names=FALSE)

# Significant
sig_df <- out_df[!is.na(out_df$padj) & out_df$padj < opt$padj &
                 abs(out_df$log2FoldChange_shrunk) > opt$lfc, ]

write.table(sig_df,
            file.path(opt$outdir, sprintf("%s_te_de_sig.tsv", comp)),
            sep="\t", quote=FALSE, row.names=FALSE)

cat(sprintf("[INF] Total loci tested: %s\n", format(nrow(out_df), big.mark=",")))
cat(sprintf("[INF] Significant (padj<%.2f, |LFC|>%.1f): %d\n",
            opt$padj, opt$lfc, nrow(sig_df)))
if (nrow(sig_df) > 0) {
  cat(sprintf("[INF] Top hits:\n"))
  print(head(sig_df[, c("locus_id","te_family","te_class",
                         "log2FoldChange_shrunk","padj")], 10),
        row.names=FALSE)
}
cat(sprintf("[OK] M12 DESeq2 complete. Output: %s\n", opt$outdir))
