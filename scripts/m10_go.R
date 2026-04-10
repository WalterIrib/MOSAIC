#!/usr/bin/env Rscript
# =============================================================================
# MOSAIC M10 — GO Enrichment Analysis
# Design Document v0.7
#
# Three analyses:
#   1. Canonical genes (protein_coding) — standard ORA
#   2. lncRNA — separate ORA with lncRNA universe
#   3. TE host genes — genes hosting any expressed TE (custom universe)
#
# GO databases by organism:
#   human   → org.Hs.eg.db / hsa
#   mouse   → org.Mm.eg.db / mmu
#   rat     → org.Rn.eg.db / rno
#   zebrafish → org.Dr.eg.db / dre
#
# Usage:
#   Rscript m10_go.R \
#     --de-genes    deseq2/genes_tumor_vs_normal_DE_significant.tsv \
#     --all-genes   integrated/all_genes_counts_clean.tsv \
#     --lncrna-de   deseq2/lncrna_DE_significant.tsv \        [optional]
#     --te-de       deseq2/TE_tumor_vs_normal_DE.tsv \
#     --gene-type-map integrated/gene_type_map.tsv \
#     --organism    human \
#     --outdir      go_results/ \
#     --padj        0.05 \
#     --min-gs      10 \
#     --max-gs      500
# =============================================================================

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(DOSE)
  library(optparse)
})

# ─────────────────────────────────────────────────────────────────────────────
# CLI
# ─────────────────────────────────────────────────────────────────────────────

option_list <- list(
  make_option("--de-genes",     type="character", help="Significant DE genes TSV (gene_id col)"),
  make_option("--all-genes",    type="character", help="all_genes_counts_clean.tsv (universe)"),
  make_option("--lncrna-de",    type="character", default=NULL, help="lncRNA DE significant TSV"),
  make_option("--te-de",        type="character", default=NULL, help="TE DE results TSV"),
  make_option("--te-host-universe", type="character", default=NULL,
              help="te_host_universe.txt from M11 (one gene_id per line)"),
  make_option("--te-host-genes",   type="character", default=NULL,
              help="te_host_genes.tsv from M11 (for DE gene filtering)"),
  make_option("--te-promoter-intersect", type="character", default=NULL,
              help="te_promoter_intersect.bed from M11b (col4=TE_name, col10=gene_id)"),
  make_option("--gene-type-map",type="character", help="gene_type_map.tsv"),
  make_option("--organism",     type="character", default="human",
              help="human|mouse|rat|zebrafish [default: human]"),
  make_option("--outdir",       type="character", help="Output directory"),
  make_option("--padj",         type="double", default=0.05,
              help="GO padj threshold [default: 0.05]"),
  make_option("--min-gs",       type="integer", default=10,
              help="Min gene set size [default: 10]"),
  make_option("--max-gs",       type="integer", default=500,
              help="Max gene set size [default: 500]"),
  make_option("--kegg-organism", type="character", default="hsa",
              help="KEGG organism code [default: hsa]"),
  make_option("--run-kegg",     action="store_true", default=TRUE,
              help="Run KEGG enrichment [default: TRUE]"),
  make_option("--run-do",       action="store_true", default=TRUE,
              help="Run Disease Ontology enrichment [default: TRUE]")
)

opt <- parse_args(OptionParser(option_list=option_list))

for (arg in c("de-genes", "all-genes", "gene-type-map", "outdir")) {
  if (is.null(opt[[arg]])) stop(paste("Missing required argument:", arg))
}

dir.create(opt$outdir, recursive=TRUE, showWarnings=FALSE)

# ─────────────────────────────────────────────────────────────────────────────
# Organism setup
# ─────────────────────────────────────────────────────────────────────────────

org_map <- list(
  human     = list(db="org.Hs.eg.db", kegg="hsa"),
  mouse     = list(db="org.Mm.eg.db", kegg="mmu"),
  rat       = list(db="org.Rn.eg.db", kegg="rno"),
  zebrafish = list(db="org.Dr.eg.db", kegg="dre")
)

if (!opt$organism %in% names(org_map)) {
  stop(sprintf("Unsupported organism: %s. Use: %s",
               opt$organism, paste(names(org_map), collapse="|")))
}

org_db_name <- org_map[[opt$organism]]$db
suppressPackageStartupMessages(
  library(org_db_name, character.only=TRUE)
)
org_db <- get(org_db_name)

cat("\n========================================\n")
cat("MOSAIC M10 — GO Enrichment\n")
cat("========================================\n")
cat(sprintf("  Organism:  %s (%s)\n", opt$organism, org_db_name))
cat(sprintf("  Outdir:    %s\n", opt$outdir))
cat(sprintf("  padj:      %.3f\n", opt$padj))
cat(sprintf("  Gene sets: %d–%d\n\n", opt[["min-gs"]], opt[["max-gs"]]))


# ─────────────────────────────────────────────────────────────────────────────
# Helper: ENSEMBL → ENTREZ conversion
# ─────────────────────────────────────────────────────────────────────────────

ensembl_to_entrez <- function(gene_ids, org_db) {
  # Strip version suffix (ENSG00000123.4 → ENSG00000123)
  gene_ids_clean <- sub("\\..*", "", gene_ids)

  mapping <- bitr(gene_ids_clean,
                  fromType="ENSEMBL",
                  toType="ENTREZID",
                  OrgDb=org_db,
                  drop=TRUE)
  return(mapping)
}


# ─────────────────────────────────────────────────────────────────────────────
# Helper: run ORA for one gene list
# ─────────────────────────────────────────────────────────────────────────────

run_ora <- function(gene_entrez, universe_entrez, label, outdir,
                    org_db, padj_cutoff=0.05, min_gs=10, max_gs=500) {

  cat(sprintf("\n[INF] GO ORA: %s (%d genes, universe=%d)\n",
              label, length(gene_entrez), length(universe_entrez)))

  if (length(gene_entrez) < 5) {
    cat(sprintf("[WRN] Too few genes (%d) for GO ORA: %s. Skipping.\n",
                length(gene_entrez), label))
    return(NULL)
  }

  results <- list()

  for (ont in c("BP", "MF", "CC")) {
    tryCatch({
      ego <- enrichGO(
        gene          = gene_entrez,
        universe      = universe_entrez,
        OrgDb         = org_db,
        ont           = ont,
        pAdjustMethod = "BH",
        pvalueCutoff  = padj_cutoff,
        qvalueCutoff  = 0.2,
        minGSSize     = min_gs,
        maxGSSize     = max_gs,
        readable      = TRUE
      )

      if (!is.null(ego) && nrow(ego@result) > 0) {
        n_sig <- sum(ego@result$p.adjust < padj_cutoff)
        cat(sprintf("[INF]   %s %s: %d significant terms\n", label, ont, n_sig))

        # Write full results
        out_f <- file.path(outdir, sprintf("%s_GO_%s.tsv", label, ont))
        write.table(ego@result, out_f, sep="\t", quote=FALSE, row.names=FALSE)

        results[[ont]] <- ego
      } else {
        cat(sprintf("[INF]   %s %s: no significant terms\n", label, ont))
      }
    }, error=function(e) {
      cat(sprintf("[WRN]   %s %s: error — %s\n", label, ont, e$message))
    })
  }

  # KEGG pathway enrichment
  tryCatch({
    kegg_org <- getOption("mosaic.kegg_org", "hsa")
    ekegg <- enrichKEGG(
      gene          = gene_entrez,
      organism      = kegg_org,
      pAdjustMethod = "BH",
      pvalueCutoff  = padj_cutoff,
      qvalueCutoff  = 0.2,
      minGSSize     = min_gs,
      maxGSSize     = max_gs
    )
    if (!is.null(ekegg) && nrow(ekegg@result) > 0) {
      n_sig <- sum(ekegg@result$p.adjust < padj_cutoff)
      cat(sprintf("[INF]   %s KEGG: %d significant pathways\n", label, n_sig))
      out_f <- file.path(outdir, sprintf("%s_KEGG.tsv", label))
      write.table(ekegg@result, out_f, sep="\t", quote=FALSE, row.names=FALSE)
      results[["KEGG"]] <- ekegg
    } else {
      cat(sprintf("[INF]   %s KEGG: no significant pathways\n", label))
    }
  }, error=function(e) {
    cat(sprintf("[WRN]   %s KEGG: %s\n", label, e$message))
  })

  # Disease Ontology enrichment
  tryCatch({
    edo <- enrichDO(
      gene          = gene_entrez,
      ont           = "DO",
      pAdjustMethod = "BH",
      pvalueCutoff  = padj_cutoff,
      qvalueCutoff  = 0.2,
      minGSSize     = min_gs,
      maxGSSize     = max_gs,
      readable      = TRUE
    )
    if (!is.null(edo) && nrow(edo@result) > 0) {
      n_sig <- sum(edo@result$p.adjust < padj_cutoff)
      cat(sprintf("[INF]   %s DO: %d significant disease terms\n", label, n_sig))
      out_f <- file.path(outdir, sprintf("%s_DO.tsv", label))
      write.table(edo@result, out_f, sep="\t", quote=FALSE, row.names=FALSE)
      results[["DO"]] <- edo
    } else {
      cat(sprintf("[INF]   %s DO: no significant disease terms\n", label))
    }
  }, error=function(e) {
    cat(sprintf("[WRN]   %s DO: %s\n", label, e$message))
  })

  return(results)
}


# ─────────────────────────────────────────────────────────────────────────────
# Load data
# ─────────────────────────────────────────────────────────────────────────────

cat("[INF] Loading DE genes...\n")
de_genes_df <- read.table(opt[["de-genes"]], header=TRUE, sep="\t",
                           check.names=FALSE)

# Accept both 'gene_id' column or row names
if ("gene_id" %in% colnames(de_genes_df)) {
  de_gene_ids <- de_genes_df$gene_id
} else {
  de_gene_ids <- rownames(de_genes_df)
}
cat(sprintf("[INF] DE genes: %d\n", length(de_gene_ids)))

cat("[INF] Loading universe (all expressed genes)...\n")
all_genes_df <- read.table(opt[["all-genes"]], header=TRUE, sep="\t",
                            row.names=1, check.names=FALSE)
# Keep only genes with at least 1 count in any sample
numeric_cols <- sapply(all_genes_df, is.numeric)
counts_mat   <- as.matrix(all_genes_df[, numeric_cols, drop=FALSE])
expressed    <- rowSums(counts_mat) > 0
universe_ids <- rownames(all_genes_df)[expressed]
cat(sprintf("[INF] Universe (expressed): %d genes\n", length(universe_ids)))

cat("[INF] Loading gene type map...\n")
gene_type_map <- read.table(opt[["gene-type-map"]], header=TRUE, sep="\t",
                             check.names=FALSE)


# ─────────────────────────────────────────────────────────────────────────────
# Convert to ENTREZ
# ─────────────────────────────────────────────────────────────────────────────

cat("[INF] Converting ENSEMBL → ENTREZ...\n")

de_map       <- ensembl_to_entrez(de_gene_ids, org_db)
universe_map <- ensembl_to_entrez(universe_ids, org_db)

de_entrez       <- unique(de_map$ENTREZID)
universe_entrez <- unique(universe_map$ENTREZID)

cat(sprintf("[INF] DE genes mapped: %d / %d\n",
            length(de_entrez), length(de_gene_ids)))
cat(sprintf("[INF] Universe mapped: %d / %d\n",
            length(universe_entrez), length(universe_ids)))


# ─────────────────────────────────────────────────────────────────────────────
# 1. Canonical genes GO ORA
# ─────────────────────────────────────────────────────────────────────────────

cat("\n--- Canonical genes GO ORA ---\n")

# Split DE genes by direction if log2FoldChange column exists
if ("log2FoldChange" %in% colnames(de_genes_df)) {
  up_ids   <- de_gene_ids[!is.na(de_genes_df$log2FoldChange) &
                           de_genes_df$log2FoldChange > 0]
  down_ids <- de_gene_ids[!is.na(de_genes_df$log2FoldChange) &
                           de_genes_df$log2FoldChange < 0]

  up_map   <- ensembl_to_entrez(up_ids, org_db)
  down_map <- ensembl_to_entrez(down_ids, org_db)

  run_ora(unique(up_map$ENTREZID),   universe_entrez, "genes_UP",
          opt$outdir, org_db, opt$padj, opt[["min-gs"]], opt[["max-gs"]])
  run_ora(unique(down_map$ENTREZID), universe_entrez, "genes_DOWN",
          opt$outdir, org_db, opt$padj, opt[["min-gs"]], opt[["max-gs"]])
} else {
  run_ora(de_entrez, universe_entrez, "genes_ALL",
          opt$outdir, org_db, opt$padj, opt[["min-gs"]], opt[["max-gs"]])
}


# ─────────────────────────────────────────────────────────────────────────────
# 2. lncRNA GO ORA (separate universe)
# ─────────────────────────────────────────────────────────────────────────────

if (!is.null(opt[["lncrna-de"]]) && file.exists(opt[["lncrna-de"]])) {
  cat("\n--- lncRNA GO ORA ---\n")

  lncrna_de_df <- read.table(opt[["lncrna-de"]], header=TRUE, sep="\t",
                              check.names=FALSE)
  lncrna_de_ids <- if ("gene_id" %in% colnames(lncrna_de_df))
                     lncrna_de_df$gene_id else rownames(lncrna_de_df)

  # lncRNA universe: all expressed lncRNA genes
  lncrna_universe_ids <- gene_type_map$gene_id[
    gene_type_map$gene_type == "lncRNA" &
    gene_type_map$gene_id %in% universe_ids
  ]

  lncrna_de_map  <- ensembl_to_entrez(lncrna_de_ids, org_db)
  lncrna_uni_map <- ensembl_to_entrez(lncrna_universe_ids, org_db)

  run_ora(unique(lncrna_de_map$ENTREZID),
          unique(lncrna_uni_map$ENTREZID),
          "lncRNA",
          opt$outdir, org_db, opt$padj, opt[["min-gs"]], opt[["max-gs"]])
}


# ─────────────────────────────────────────────────────────────────────────────
# 3. TE host genes GO ORA
# ─────────────────────────────────────────────────────────────────────────────

if (!is.null(opt[["te-de"]]) && file.exists(opt[["te-de"]])) {
  cat("\n--- TE host genes GO ORA ---\n")

  te_de_df <- read.table(opt[["te-de"]], header=TRUE, sep="\t",
                          check.names=FALSE)

  # Significant TE families
  sig_te <- te_de_df[!is.na(te_de_df$padj) & te_de_df$padj < opt$padj, ]
  cat(sprintf("[INF] Significant TE families: %d\n", nrow(sig_te)))

  # ── TE host universe ────────────────────────────────────────────────────────
  # Use M11 te_host_universe.txt if available, otherwise fallback to proxy
  if (!is.null(opt[["te-host-universe"]]) &&
      file.exists(opt[["te-host-universe"]])) {

    cat("[INF] Using M11 TE host universe (rmsk-based).\n")
    te_host_universe_ids <- readLines(opt[["te-host-universe"]])
    # Intersect with expressed genes
    te_host_universe_ids <- te_host_universe_ids[
      te_host_universe_ids %in% universe_ids
    ]
    cat(sprintf("[INF] TE host universe (M11, expressed): %d genes\n",
                length(te_host_universe_ids)))

  } else {
    cat("[WRN] M11 universe not provided — using protein_coding proxy.\n")
    te_host_universe_ids <- gene_type_map$gene_id[
      gene_type_map$gene_type == "protein_coding" &
      gene_type_map$gene_id %in% universe_ids
    ]
    cat(sprintf("[INF] TE host universe (proxy): %d genes\n",
                length(te_host_universe_ids)))
  }

  te_host_uni_map         <- ensembl_to_entrez(te_host_universe_ids, org_db)
  te_host_universe_entrez <- unique(te_host_uni_map$ENTREZID)

  # ── TE host DE genes ─────────────────────────────────────────────────────────
  # Filter de_gene_ids to SIGNIFICANT DE genes only (padj < threshold)
  # This is essential: using all tested genes (including non-significant) as
  # foreground would make it larger than the background, breaking enrichGO.
  sig_mask <- !is.na(de_genes_df$padj) & de_genes_df$padj < opt$padj
  sig_de_ids <- de_gene_ids[sig_mask]
  cat(sprintf("[INF] Significant DE genes (padj<%.2f): %d / %d\n",
              opt$padj, length(sig_de_ids), length(de_gene_ids)))

  # Use M11 te_host_genes.tsv to restrict DE genes to actual TE-hosting genes
  # Foreground = sig DE genes hosting ≥1 SIGNIFICANTLY DE TE family
  if (!is.null(opt[["te-host-genes"]]) &&
      file.exists(opt[["te-host-genes"]])) {

    cat("[INF] Using M11 te_host_genes.tsv for DE gene filtering.\n")
    te_host_df <- read.table(opt[["te-host-genes"]], header=TRUE, sep="\t",
                             check.names=FALSE)

    if (nrow(sig_te) > 0 && "te_families" %in% colnames(te_host_df)) {
      # Filter to genes hosting ≥1 significantly DE TE family
      sig_te_families <- sig_te$family_id
      hosts_de_te <- vapply(strsplit(te_host_df$te_families, ","), function(fams) {
        any(trimws(fams) %in% sig_te_families)
      }, logical(1))
      host_gene_ids <- te_host_df$gene_id[hosts_de_te]
      cat(sprintf("[INF] Genes hosting ≥1 DE TE family: %d / %d (sig TE families: %d)\n",
                  length(host_gene_ids), nrow(te_host_df), length(sig_te_families)))
    } else if (nrow(sig_te) == 0) {
      cat("[INF] No significantly DE TE families — TE host GO foreground will be empty.\n")
      host_gene_ids <- character(0)
    } else {
      # te_families column missing: fall back to any TE host
      host_gene_ids <- te_host_df$gene_id
      cat(sprintf("[WRN] te_families column missing — using all TE hosts: %d genes\n",
                  length(host_gene_ids)))
    }

    te_host_de_ids <- sig_de_ids[sig_de_ids %in% host_gene_ids]
    cat(sprintf("[INF] Significant DE genes hosting ≥1 DE TE family: %d / %d\n",
                length(te_host_de_ids), length(sig_de_ids)))

  } else {
    cat("[WRN] M11 te_host_genes.tsv not provided — using significant DE protein_coding.\n")
    te_host_de_ids <- sig_de_ids[sig_de_ids %in%
                                   gene_type_map$gene_id[
                                     gene_type_map$gene_type == "protein_coding"]]
  }

  if (length(te_host_de_ids) >= 5) {
    te_host_de_map <- ensembl_to_entrez(te_host_de_ids, org_db)

    run_ora(unique(te_host_de_map$ENTREZID),
            te_host_universe_entrez,
            "TE_host",
            opt$outdir, org_db, opt$padj, opt[["min-gs"]], opt[["max-gs"]])
  } else {
    cat(sprintf("[WRN] Too few DE TE-host genes (%d). Skipping TE host GO.\n",
                length(te_host_de_ids)))
  }

  # ── 4. TE carrier genes GO ORA (body + promoter, regardless of gene DE) ────
  # Foreground = ALL expressed genes harboring ≥1 DE TE family in gene body OR
  # promoter (2kb upstream), regardless of whether the gene itself is DE.
  # Rationale: TE differential expression may point to genes and functions that
  # escape traditional DE analysis. This analysis is complementary to the TE host
  # GO (section 4.2 of the report) which requires genes to also be DE.
  cat("\n--- TE carrier genes GO ORA (body+promoter, regardless of gene DE) ---\n")

  carrier_ids      <- character(0)
  prom_carrier_ids <- character(0)

  if (nrow(sig_te) > 0) {
    sig_te_fams_carrier <- sig_te$family_id

    # ── Gene body carriers ────────────────────────────────────────────────────
    if (!is.null(opt[["te-host-genes"]]) && file.exists(opt[["te-host-genes"]])) {
      if (!exists("te_host_df")) {
        te_host_df <- read.table(opt[["te-host-genes"]], header=TRUE, sep="\t",
                                 check.names=FALSE)
      }
      if ("te_families" %in% colnames(te_host_df)) {
        carries_body <- vapply(strsplit(te_host_df$te_families, ","), function(fams) {
          any(trimws(fams) %in% sig_te_fams_carrier)
        }, logical(1))
        carrier_ids <- te_host_df$gene_id[carries_body]
        carrier_ids <- carrier_ids[carrier_ids %in% universe_ids]
        cat(sprintf("[INF] Body carrier genes (≥1 DE TE family in body): %d\n",
                    length(carrier_ids)))
      } else {
        cat("[WRN] te_families column missing — body carrier genes not derived.\n")
      }
    }

    # ── Promoter carriers ─────────────────────────────────────────────────────
    # Derive from te_promoter_intersect.bed: col4=TE_name (family|class|type), col10=gene_id
    if (!is.null(opt[["te-promoter-intersect"]]) &&
        file.exists(opt[["te-promoter-intersect"]])) {
      cat("[INF] Reading te_promoter_intersect.bed for promoter carrier genes...\n")
      prom_bed <- tryCatch(
        read.table(opt[["te-promoter-intersect"]], header=FALSE, sep="\t",
                   check.names=FALSE, stringsAsFactors=FALSE,
                   quote="", fill=TRUE, comment.char=""),
        error=function(e) { cat("[WRN] Failed to read promoter intersect:", e$message, "\n"); NULL }
      )
      if (!is.null(prom_bed) && ncol(prom_bed) >= 10) {
        prom_df <- data.frame(
          te_name = prom_bed[[4]],
          gene_id = prom_bed[[10]],
          stringsAsFactors=FALSE
        )
        # Family name = part before first "|"
        prom_df$family_id <- sub("\\|.*", "", prom_df$te_name)
        prom_carrier_ids <- unique(prom_df$gene_id[prom_df$family_id %in% sig_te_fams_carrier])
        prom_carrier_ids <- prom_carrier_ids[prom_carrier_ids %in% universe_ids]
        cat(sprintf("[INF] Promoter carrier genes (≥1 DE TE family in promoter): %d\n",
                    length(prom_carrier_ids)))
      }
    }

    # ── Body-only ORA ────────────────────────────────────────────────────────
    if (length(carrier_ids) >= 5) {
      body_map <- ensembl_to_entrez(carrier_ids, org_db)
      run_ora(unique(body_map$ENTREZID),
              te_host_universe_entrez,
              "TE_carrier_body",
              opt$outdir, org_db, opt$padj, opt[["min-gs"]], opt[["max-gs"]])
    } else {
      cat(sprintf("[WRN] Too few body carrier genes (%d). Skipping body carrier GO.\n",
                  length(carrier_ids)))
    }

    # ── Promoter-only ORA ────────────────────────────────────────────────────
    if (length(prom_carrier_ids) >= 5) {
      prom_map <- ensembl_to_entrez(prom_carrier_ids, org_db)
      run_ora(unique(prom_map$ENTREZID),
              te_host_universe_entrez,
              "TE_carrier_promoter",
              opt$outdir, org_db, opt$padj, opt[["min-gs"]], opt[["max-gs"]])
    } else {
      cat(sprintf("[WRN] Too few promoter carrier genes (%d). Skipping promoter carrier GO.\n",
                  length(prom_carrier_ids)))
    }

    # ── Combined body + promoter ──────────────────────────────────────────────
    combined_carrier_ids <- unique(c(carrier_ids, prom_carrier_ids))
    cat(sprintf("[INF] Combined carrier genes (body OR promoter): %d\n",
                length(combined_carrier_ids)))

    if (length(combined_carrier_ids) >= 5) {
      carrier_map <- ensembl_to_entrez(combined_carrier_ids, org_db)
      run_ora(unique(carrier_map$ENTREZID),
              te_host_universe_entrez,
              "TE_carrier",
              opt$outdir, org_db, opt$padj, opt[["min-gs"]], opt[["max-gs"]])
    } else {
      cat(sprintf("[WRN] Too few carrier genes (%d body+promoter). Skipping TE carrier GO.\n",
                  length(combined_carrier_ids)))
    }
  } else {
    cat("[INF] No significantly DE TE families — skipping TE carrier GO.\n")
  }
}


# ─────────────────────────────────────────────────────────────────────────────
# Write summary
# ─────────────────────────────────────────────────────────────────────────────

result_files <- list.files(opt$outdir, pattern="_GO_.*\\.tsv$", full.names=TRUE)

summary_df <- do.call(rbind, lapply(result_files, function(f) {
  df <- tryCatch(read.table(f, header=TRUE, sep="\t"), error=function(e) NULL)
  if (is.null(df)) return(NULL)
  sig <- sum(df$p.adjust < opt$padj, na.rm=TRUE)
  data.frame(file=basename(f), n_terms=nrow(df), n_significant=sig)
}))

if (!is.null(summary_df) && nrow(summary_df) > 0) {
  write.table(summary_df,
              file.path(opt$outdir, "go_summary.tsv"),
              sep="\t", quote=FALSE, row.names=FALSE)
  cat("\n[INF] GO Summary:\n")
  print(summary_df)
}

cat(sprintf("\n[OK] M10 GO enrichment complete.\n"))
cat(sprintf("[OK] Results in: %s\n\n", opt$outdir))
