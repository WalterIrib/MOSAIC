#!/usr/bin/env Rscript
# mosaic/scripts/m08b_qc_plots.R
# M08b — QC and exploratory visualisation
#   1. PCA of samples (DESeq2 rlog or vst)
#   2. Sample-sample correlation heatmap
#   3. Top DE genes heatmap (z-score)
#
# Run after M08 (DESeq2). Reads normalized_counts.tsv + DE results.

suppressPackageStartupMessages({
  library(optparse)
  library(ggplot2)
  library(dplyr)
})

option_list <- list(
  make_option("--norm-counts", type="character",
              help="normalized_counts.tsv (gene_id + sample cols)"),
  make_option("--de-genes",    type="character",
              help="genes_*_DE.tsv from M08"),
  make_option("--coldata",     type="character", default="/tmp/coldata.tsv"),
  make_option("--gene-map",    type="character", default=NULL,
              help="gene_type_map.tsv for gene symbols"),
  make_option("--outdir",      type="character", default="mosaic_output/deseq2"),
  make_option("--treatment",   type="character", default="tumor"),
  make_option("--control",     type="character", default="normal"),
  make_option("--top-n",       type="integer",   default=50L,
              help="Top DE genes for heatmap [default: 50]"),
  make_option("--padj",        type="double",    default=0.05),
  make_option("--lfc",         type="double",    default=1.0)
)

opt <- parse_args(OptionParser(option_list=option_list))
dir.create(opt$outdir, recursive=TRUE, showWarnings=FALSE)

# ── Helpers ────────────────────────────────────────────────────────────────────
safe_read <- function(path, row.names=NULL, ...) {
  if (!is.null(path) && file.exists(path))
    tryCatch(
      read.table(path, header=TRUE, sep="\t", check.names=FALSE,
                 row.names=row.names, ...),
      error=function(e) { cat("[WARN]", path, ":", conditionMessage(e), "\n"); NULL }
    )
  else { cat("[WARN] Not found:", path, "\n"); NULL }
}

col_cond <- c(tumor="#E74C3C", normal="#2980B9",
              treatment="#E74C3C", control="#2980B9")

# ── Load ───────────────────────────────────────────────────────────────────────
cat("[INF] Loading data...\n")
nc_raw   <- safe_read(opt[["norm-counts"]])
de_genes <- safe_read(opt[["de-genes"]])
coldata  <- safe_read(opt$coldata)
gene_map <- safe_read(opt[["gene-map"]])

if (is.null(nc_raw)) { cat("[ERROR] norm-counts required\n"); quit(status=1) }

# Build count matrix with gene_id as rownames
if ("gene_id" %in% colnames(nc_raw)) {
  nc <- as.matrix(nc_raw[, !colnames(nc_raw) %in% "gene_id"])
  rownames(nc) <- nc_raw$gene_id
} else {
  nc <- as.matrix(nc_raw)
}

# Sample condition vector
if (!is.null(coldata)) {
  cond <- setNames(coldata$condition, coldata$sample_id)
} else {
  cond <- setNames(rep("unknown", ncol(nc)), colnames(nc))
}

# Gene name lookup
gene_names <- if (!is.null(gene_map)) {
  setNames(gene_map$gene_name, gene_map$gene_id)
} else character(0)

label_gene <- function(gids) {
  nm <- gene_names[gids]
  ifelse(!is.na(nm) & nm != "", nm, gids)
}

# log2 transform for PCA/heatmap
nc_log <- log2(nc + 1)

# ═══════════════════════════════════════════════════════════════════════════════
# 1. PCA
# ═══════════════════════════════════════════════════════════════════════════════
cat("[INF] Computing PCA...\n")
tryCatch({
  # Use top 500 most variable genes
  gene_vars <- apply(nc_log, 1, var, na.rm=TRUE)
  top500    <- names(sort(gene_vars, decreasing=TRUE))[1:min(500, length(gene_vars))]
  mat_pca   <- t(nc_log[top500, ])

  pca_res <- prcomp(mat_pca, scale.=TRUE, center=TRUE)
  var_exp <- summary(pca_res)$importance[2, 1:min(6, ncol(pca_res$x))] * 100

  pca_df <- data.frame(
    sample_id = rownames(pca_res$x),
    PC1       = pca_res$x[,1],
    PC2       = pca_res$x[,2],
    condition = cond[rownames(pca_res$x)]
  )

  p_pca <- ggplot(pca_df, aes(x=PC1, y=PC2, color=condition, label=sample_id)) +
    geom_point(size=5, alpha=0.9) +
    ggrepel::geom_text_repel(size=3.5, show.legend=FALSE) +
    scale_color_manual(values=col_cond, na.value="grey60") +
    labs(
      title    = "PCA — Top 500 most variable genes",
      subtitle = sprintf("PC1 explains %.1f%% variance; PC2 explains %.1f%%",
                         var_exp[1], var_exp[2]),
      x = sprintf("PC1 (%.1f%%)", var_exp[1]),
      y = sprintf("PC2 (%.1f%%)", var_exp[2]),
      color = "Condition"
    ) +
    theme_minimal(base_size=13) +
    theme(legend.position="top",
          plot.subtitle=element_text(size=10, color="grey50"))

  ggsave(file.path(opt$outdir, "pca_samples.pdf"),
         p_pca, width=7, height=6)
  ggsave(file.path(opt$outdir, "pca_samples.png"),
         p_pca, width=7, height=6, dpi=150)

  # Scree plot
  scree_df <- data.frame(
    PC       = paste0("PC", seq_along(var_exp)),
    variance = var_exp
  )
  scree_df$PC <- factor(scree_df$PC, levels=scree_df$PC)

  p_scree <- ggplot(scree_df, aes(x=PC, y=variance)) +
    geom_col(fill="#2980B9", alpha=0.8) +
    geom_text(aes(label=sprintf("%.1f%%", variance)), vjust=-0.3, size=3.5) +
    labs(title="PCA scree plot", y="% variance explained", x=NULL) +
    theme_minimal(base_size=12)

  ggsave(file.path(opt$outdir, "pca_scree.pdf"),
         p_scree, width=5, height=4)

  # Save coordinates
  write.table(pca_df, file.path(opt$outdir, "pca_coordinates.tsv"),
              sep="\t", quote=FALSE, row.names=FALSE)

  cat(sprintf("[OK] PCA: PC1=%.1f%%, PC2=%.1f%%\n", var_exp[1], var_exp[2]))
}, error=function(e) cat("[WARN] PCA failed:", conditionMessage(e), "\n"))


# ═══════════════════════════════════════════════════════════════════════════════
# 2. Sample-sample correlation heatmap
# ═══════════════════════════════════════════════════════════════════════════════
cat("[INF] Sample correlation heatmap...\n")
tryCatch({
  if (!requireNamespace("pheatmap", quietly=TRUE)) {
    cat("[WARN] pheatmap not installed\n")
  } else {
    library(pheatmap)
    library(RColorBrewer)

    cor_mat <- cor(nc_log, method="pearson")

    # Annotation
    ann_col <- data.frame(
      Condition = cond[colnames(cor_mat)],
      row.names = colnames(cor_mat)
    )
    ann_colors <- list(
      Condition = col_cond[unique(ann_col$Condition)]
    )
    ann_colors$Condition <- ann_colors$Condition[!is.na(ann_colors$Condition)]

    pheatmap(
      cor_mat,
      color           = colorRampPalette(c("#2C3E50","#85C1E9","white"))(100),
      annotation_col  = ann_col,
      annotation_row  = ann_col,
      annotation_colors= ann_colors,
      display_numbers = TRUE,
      number_format   = "%.3f",
      fontsize_number = 9,
      main            = "Sample-sample Pearson correlation (log2 norm counts)",
      filename        = file.path(opt$outdir, "sample_correlation_heatmap.pdf"),
      width=7, height=6
    )
    cat("[OK] Sample correlation heatmap saved.\n")
  }
}, error=function(e) cat("[WARN] Correlation heatmap failed:", conditionMessage(e), "\n"))


# ═══════════════════════════════════════════════════════════════════════════════
# 3. Top DE genes heatmap (z-score)
# ═══════════════════════════════════════════════════════════════════════════════
cat("[INF] Top DE genes heatmap...\n")
tryCatch({
  if (is.null(de_genes)) {
    cat("[WARN] No DE genes file — skipping heatmap\n")
  } else if (!requireNamespace("pheatmap", quietly=TRUE)) {
    cat("[WARN] pheatmap not installed\n")
  } else {
    library(pheatmap)

    # Select top DE genes by |LFC|, filtered by padj
    top_de <- de_genes |>
      filter(!is.na(padj), padj < opt$padj,
             abs(log2FoldChange) > opt$lfc) |>
      arrange(padj, desc(abs(log2FoldChange))) |>
      head(opt[["top-n"]])

    if (nrow(top_de) < 3) {
      cat("[WARN] Fewer than 3 significant DE genes — skipping heatmap\n")
    } else {
      # gene_id column
      gid_col <- if ("gene_id" %in% colnames(top_de)) "gene_id" else
                 colnames(top_de)[ncol(top_de)]
      gids    <- top_de[[gid_col]]

      # Subset expression matrix
      common_g <- intersect(gids, rownames(nc_log))
      if (length(common_g) < 3) {
        cat("[WARN] Genes not found in count matrix — skipping heatmap\n")
      } else {
        mat_heat <- nc_log[common_g, , drop=FALSE]

        # Z-score across samples
        mat_z <- t(scale(t(mat_heat)))

        # Row labels: gene symbols
        rownames(mat_z) <- label_gene(rownames(mat_z))

        # LFC annotation on rows
        lfc_ann <- top_de |>
          slice(match(common_g, .[[gid_col]])) |>
          select(log2FoldChange) |>
          mutate(Direction = ifelse(log2FoldChange > 0,
                                   opt$treatment, opt$control))
        rownames(lfc_ann) <- label_gene(common_g)
        lfc_ann <- lfc_ann[, "Direction", drop=FALSE]

        # Column annotation
        ann_col <- data.frame(
          Condition = cond[colnames(mat_z)],
          row.names = colnames(mat_z)
        )

        ann_colors <- list(
          Condition = col_cond[unique(ann_col$Condition)],
          Direction = setNames(c("#E74C3C","#2980B9"),
                               c(opt$treatment, opt$control))
        )
        ann_colors <- lapply(ann_colors, function(x) x[!is.na(x)])

        n_genes_label <- nrow(mat_z)
        pheatmap(
          mat_z,
          color            = colorRampPalette(c("#2980B9","white","#E74C3C"))(100),
          annotation_col   = ann_col,
          annotation_row   = lfc_ann,
          annotation_colors= ann_colors,
          cluster_rows     = TRUE,
          cluster_cols     = TRUE,
          show_rownames    = TRUE,
          show_colnames    = TRUE,
          fontsize_row     = max(5, 10 - n_genes_label / 10),
          main             = sprintf("Top %d DE genes (z-score, padj<%.2f, |LFC|>%.1f)",
                                     n_genes_label, opt$padj, opt$lfc),
          filename         = file.path(opt$outdir, "top_de_heatmap.pdf"),
          width=8, height=max(6, n_genes_label * 0.18 + 2)
        )
        cat(sprintf("[OK] Heatmap: %d genes × %d samples\n",
                    nrow(mat_z), ncol(mat_z)))
      }
    }
  }
}, error=function(e) cat("[WARN] DE heatmap failed:", conditionMessage(e), "\n"))

cat("[OK] M08b QC plots complete. Output:", opt$outdir, "\n")
