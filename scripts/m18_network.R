#!/usr/bin/env Rscript
# mosaic/scripts/m18_network.R
# M18 — Network analysis: STRINGdb protein interaction network + WGCNA
#
# Two complementary approaches:
# 1. STRINGdb: PPI network for DE genes (confidence ≥ 400)
#    - Hub gene identification (degree centrality)
#    - Community detection (Louvain)
# 2. WGCNA: co-expression modules from normalised counts
#    - Module-trait correlation
#    - Hub genes per module
#    - Module GO enrichment
#
# Outputs: network TSV, hub genes, module assignments, plots

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

# ── CLI ────────────────────────────────────────────────────────────────────────
option_list <- list(
  make_option("--de-genes",       type="character", help="DESeq2 DE genes TSV"),
  make_option("--norm-counts",    type="character", help="Normalised counts TSV"),
  make_option("--coldata",        type="character", default="/tmp/coldata.tsv"),
  make_option("--gene-map",       type="character", help="gene_type_map.tsv"),
  make_option("--outdir",         type="character", default="mosaic_output/m18"),
  make_option("--organism",       type="character", default="human"),
  make_option("--string-species", type="integer",   default=9606L),
  make_option("--string-score",   type="integer",   default=400L,
              help="STRINGdb min score (0-1000)"),
  make_option("--wgcna-power",    type="integer",   default=0L,
              help="Soft threshold power (0=auto)"),
  make_option("--min-module-size",type="integer",   default=20L),
  make_option("--padj",           type="double",    default=0.05),
  make_option("--lfc",            type="double",    default=1.0),
  make_option("--n-hub-genes",    type="integer",   default=20L,
              help="Number of hub genes to report per analysis")
)

opt <- parse_args(OptionParser(option_list=option_list))
dir.create(opt$outdir, recursive=TRUE, showWarnings=FALSE)

cat("[INF] M18 Network Analysis\n")

# ── Helpers ────────────────────────────────────────────────────────────────────
safe_read <- function(path, ...) {
  if (!is.null(path) && file.exists(path))
    tryCatch(read.table(path, header=TRUE, sep="\t", check.names=FALSE, ...),
             error=function(e) { cat("[WARN]", conditionMessage(e), "\n"); NULL })
  else { cat("[WARN] File not found:", path, "\n"); NULL }
}

write_tsv <- function(df, path) {
  write.table(df, path, sep="\t", quote=FALSE, row.names=FALSE)
  cat("[OK] Saved:", path, "\n")
}

# ── Load data ──────────────────────────────────────────────────────────────────
de_genes    <- safe_read(opt[["de-genes"]])
norm_counts <- safe_read(opt[["norm-counts"]], row.names=1)
coldata     <- safe_read(opt$coldata)
gene_map    <- safe_read(opt[["gene-map"]])

if (is.null(de_genes) || is.null(norm_counts)) {
  cat("[ERROR] de-genes and norm-counts are required\n")
  quit(status=1)
}

# Gene name lookup
gene_names <- if (!is.null(gene_map)) {
  setNames(gene_map$gene_name, gene_map$gene_id)
} else character(0)

# Significant DE genes
de_sig <- de_genes |>
  filter(!is.na(padj), padj < opt$padj,
         abs(log2FoldChange) > opt$lfc)

# Fallback for test data: use top 50 by pvalue if no significant genes
if (nrow(de_sig) == 0) {
  cat("[WARN] No DE genes at threshold — using top 50 by pvalue for testing
")
  de_sig <- de_genes |>
    filter(!is.na(pvalue)) |>
    arrange(pvalue) |>
    head(50)
}

cat(sprintf("[INF] Significant DE genes: %d\n", nrow(de_sig)))

# ═══════════════════════════════════════════════════════════════════════════════
# PART 1: STRINGdb PPI Network
# ═══════════════════════════════════════════════════════════════════════════════
cat("\n[INF] === STRINGdb PPI Network ===\n")

string_results <- tryCatch({
  if (!requireNamespace("STRINGdb", quietly=TRUE)) {
    cat("[WARN] STRINGdb not available — skipping PPI analysis\n")
    NULL
  } else {
    library(STRINGdb)
    library(igraph)

    # Init STRINGdb
    string_db <- STRINGdb$new(
      version      = "12.0",
      species      = opt[["string-species"]],
      score_threshold = opt[["string-score"]],
      input_directory = opt$outdir
    )

    # Map DE genes to STRING IDs
    # Use gene symbols
    de_symbols <- data.frame(
      gene_id   = de_genes$gene_id,
      gene_name = gene_names[de_genes$gene_id],
      log2FC    = de_genes$log2FoldChange,
      padj      = de_genes$padj
    ) |> filter(!is.na(gene_name), gene_name != "")

    mapped <- string_db$map(de_symbols, "gene_name", removeUnmappedRows=TRUE)
    cat(sprintf("[INF] Mapped %d/%d genes to STRING\n",
                nrow(mapped), nrow(de_symbols)))

    if (nrow(mapped) < 3) {
      cat("[WARN] Too few mapped genes for network analysis\n")
      NULL
    } else {
      # Get interactions
      interactions <- string_db$get_interactions(mapped$STRING_id)
      cat(sprintf("[INF] STRING interactions: %d\n", nrow(interactions)))

      # Build igraph
      # Deduplicate mapped vertices and edges
      mapped_uniq <- mapped[!duplicated(mapped[[1]]), ]
      interactions_clean <- interactions[
        interactions$from %in% mapped_uniq[[1]] &
        interactions$to   %in% mapped_uniq[[1]], ]
      interactions_clean <- interactions_clean[
        interactions_clean$from != interactions_clean$to, ]

      g <- graph_from_data_frame(
        d        = interactions_clean[,c("from","to","combined_score")],
        directed = FALSE,
        vertices = mapped_uniq
      )

      # Hub genes by degree
      degrees   <- degree(g)
      hub_genes <- data.frame(
        STRING_id = names(degrees),
        degree    = as.integer(degrees)
      ) |>
        left_join(mapped[,c("STRING_id","gene_name","log2FC","padj")],
                  by="STRING_id") |>
        arrange(desc(degree)) |>
        head(opt[["n-hub-genes"]])

      cat("[INF] Top hub genes (by degree):\n")
      print(hub_genes[,c("gene_name","degree","log2FC","padj")] |> head(10))

      # Community detection (Louvain)
      communities <- cluster_louvain(g)
      community_df <- data.frame(
        STRING_id = V(g)$name,
        module    = membership(communities)
      ) |>
        left_join(mapped[,c("STRING_id","gene_name","log2FC","padj")],
                  by="STRING_id") |>
        arrange(module, desc(abs(log2FC)))

      # Save network edges
      edge_df <- as_data_frame(g, what="edges") |>
        rename(STRING_from=from, STRING_to=to) |>
        left_join(mapped[,c("STRING_id","gene_name")],
                  by=c("STRING_from"="STRING_id")) |>
        rename(gene1=gene_name) |>
        left_join(mapped[,c("STRING_id","gene_name")],
                  by=c("STRING_to"="STRING_id")) |>
        rename(gene2=gene_name)

      write_tsv(edge_df,      file.path(opt$outdir, "string_network_edges.tsv"))
      write_tsv(hub_genes,    file.path(opt$outdir, "string_hub_genes.tsv"))
      write_tsv(community_df, file.path(opt$outdir, "string_communities.tsv"))

      # Network plot (top component)
      pdf(file.path(opt$outdir, "string_network.pdf"), width=12, height=10)
      tryCatch({
        plot(communities, g,
             vertex.size   = log2(degrees + 2) * 3,
             vertex.label  = ifelse(degrees >= quantile(degrees, 0.8),
                                     V(g)$gene_name, NA),
             vertex.label.cex = 0.6,
             vertex.color  = membership(communities),
             edge.width    = E(g)$combined_score / 200,
             main = sprintf("STRINGdb PPI Network — DE genes (n=%d)", vcount(g)))
      }, error=function(e) cat("[WARN] Network plot error:", conditionMessage(e), "\n"))
      dev.off()

      list(hub_genes=hub_genes, communities=community_df,
           n_nodes=vcount(g), n_edges=ecount(g))
    }
  }
}, error=function(e) {
  cat("[WARN] STRINGdb analysis failed:", conditionMessage(e), "\n")
  NULL
})

# ═══════════════════════════════════════════════════════════════════════════════
# PART 2: WGCNA Co-expression Modules
# ═══════════════════════════════════════════════════════════════════════════════
cat("\n[INF] === WGCNA Co-expression Analysis ===\n")

wgcna_results <- tryCatch({
  if (!requireNamespace("WGCNA", quietly=TRUE)) {
    cat("[WARN] WGCNA not available — skipping co-expression analysis\n")
    NULL
  } else {
    library(WGCNA)
    options(stringsAsFactors=FALSE)
    enableWGCNAThreads()

    # Prepare expression matrix
    # Use variance-stabilised or normalised counts
    # Filter: keep top 5000 most variable genes
    nc <- as.matrix(norm_counts)
    nc <- nc[rowSums(nc) > 0, ]

    # Use log2(norm + 1) for WGCNA
    expr_mat <- log2(nc + 1)

    # Select most variable genes
    gene_vars <- apply(expr_mat, 1, var)
    top_genes <- names(sort(gene_vars, decreasing=TRUE))[1:min(5000, nrow(expr_mat))]
    expr_sub  <- t(expr_mat[top_genes, ])

    cat(sprintf("[INF] WGCNA input: %d samples × %d genes\n",
                nrow(expr_sub), ncol(expr_sub)))

    # Goodness-of-fit check
    gsg <- goodSamplesGenes(expr_sub, verbose=0)
    if (!gsg$allOK) {
      expr_sub <- expr_sub[gsg$goodSamples, gsg$goodGenes]
      cat(sprintf("[INF] After QC: %d × %d\n", nrow(expr_sub), ncol(expr_sub)))
    }

    # ── Soft threshold selection ──────────────────────────────────────────────
    powers <- c(1:10, seq(12, 20, 2))
    sft    <- pickSoftThreshold(expr_sub, powerVector=powers,
                                 networkType="unsigned", verbose=0)

    power <- opt[["wgcna-power"]]
    if (power == 0) {
      power <- sft$powerEstimate
      if (is.na(power)) {
        power <- 6L  # fallback
        cat("[WARN] Auto soft threshold failed — using power=6\n")
      }
    }
    cat(sprintf("[INF] WGCNA soft threshold power: %d\n", power))

    # Save scale-free plot
    pdf(file.path(opt$outdir, "wgcna_soft_threshold.pdf"))
    plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
         type="b", xlab="Soft threshold (power)", ylab="Scale-free R²",
         main="Scale-free topology fit")
    abline(h=0.8, col="red")
    dev.off()

    # ── Module detection ──────────────────────────────────────────────────────
    net <- blockwiseModules(
      expr_sub,
      power            = power,
      networkType      = "unsigned",
      TOMType          = "unsigned",
      minModuleSize    = opt[["min-module-size"]],
      mergeCutHeight   = 0.25,
      numericLabels    = FALSE,
      pamRespectsDendro= FALSE,
      saveTOMs         = FALSE,
      verbose          = 0
    )

    module_colors <- net$colors
    n_modules     <- length(unique(module_colors)) - 1  # exclude grey
    cat(sprintf("[INF] WGCNA modules detected: %d (+ grey = unassigned)\n", n_modules))

    # Module membership table
    MEs      <- moduleEigengenes(expr_sub, module_colors)$eigengenes
    kME      <- signedKME(expr_sub, MEs)

    module_df <- data.frame(
      gene_id   = colnames(expr_sub),
      module    = module_colors
    ) |>
      mutate(gene_name = gene_names[gene_id]) |>
      left_join(
        data.frame(gene_id=colnames(expr_sub),
                   kME_inmodule = apply(kME, 1, max)),
        by="gene_id"
      ) |>
      arrange(module, desc(kME_inmodule))

    write_tsv(module_df, file.path(opt$outdir, "wgcna_module_genes.tsv"))

    # ── Module-trait correlation ──────────────────────────────────────────────
    if (!is.null(coldata)) {
      sids  <- intersect(rownames(expr_sub), coldata$sample_id)
      cdata <- coldata[match(sids, coldata$sample_id), ]
      MEs_sub <- MEs[sids, ]

      # Binary trait (1=treatment, 0=control)
      treat <- if (!is.null(opt$treatment)) opt$treatment else "tumor"
      trait <- as.numeric(cdata$condition == treat)
      trait_mat <- matrix(trait, ncol=1,
                          dimnames=list(sids, "condition"))

      module_trait_cor <- cor(MEs_sub, trait_mat, use="pairwise.complete.obs")
      module_trait_p   <- corPvalueStudent(module_trait_cor, nrow(MEs_sub))

      mt_df <- data.frame(
        module      = gsub("ME","", rownames(module_trait_cor)),
        correlation = round(module_trait_cor[,1], 4),
        pvalue      = round(module_trait_p[,1], 6)
      ) |> arrange(pvalue)

      write_tsv(mt_df, file.path(opt$outdir, "wgcna_module_trait_correlation.tsv"))
      cat("[INF] Top module-trait correlations:\n")
      print(head(mt_df, 5))

      # Heatmap
      pdf(file.path(opt$outdir, "wgcna_module_trait_heatmap.pdf"),
          width=6, height=max(4, n_modules/3))
      tryCatch(
        labeledHeatmap(
          Matrix    = module_trait_cor,
          xLabels   = "Condition",
          yLabels   = gsub("ME","", colnames(MEs)),
          ySymbols  = gsub("ME","", colnames(MEs)),
          colorLabels= FALSE,
          colors    = blueWhiteRed(50),
          textMatrix= paste0(round(module_trait_cor,2), "\n(",
                              signif(module_trait_p,1), ")"),
          setStdMargins=FALSE, cex.text=0.8,
          main="Module-Condition correlation"
        ),
        error=function(e) cat("[WARN] Heatmap error:", conditionMessage(e))
      )
      dev.off()
    }

    # ── Hub genes per module ─────────────────────────────────────────────────
    hub_df <- module_df |>
      group_by(module) |>
      slice_max(kME_inmodule, n=10) |>
      ungroup()
    write_tsv(hub_df, file.path(opt$outdir, "wgcna_hub_genes.tsv"))

    # Dendrogram
    pdf(file.path(opt$outdir, "wgcna_dendrogram.pdf"), width=14, height=6)
    plotDendroAndColors(net$dendrograms[[1]], module_colors[net$blockGenes[[1]]],
                        "Module colors", dendroLabels=FALSE, hang=0.03,
                        addGuide=TRUE, guideHang=0.05,
                        main="Gene dendrogram and module colors")
    dev.off()

    list(n_modules=n_modules, module_df=module_df, power=power)
  }
}, error=function(e) {
  cat("[WARN] WGCNA analysis failed:", conditionMessage(e), "\n")
  NULL
})

# ── Summary ────────────────────────────────────────────────────────────────────
cat("\n[INF] M18 Summary:\n")
if (!is.null(string_results))
  cat(sprintf("  STRINGdb: %d nodes, %d edges\n",
              string_results$n_nodes, string_results$n_edges))
if (!is.null(wgcna_results))
  cat(sprintf("  WGCNA: %d modules (power=%d)\n",
              wgcna_results$n_modules, wgcna_results$power))

cat("[OK] M18 complete. Results:", opt$outdir, "\n")
