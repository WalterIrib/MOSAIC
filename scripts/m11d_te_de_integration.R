#!/usr/bin/env Rscript
# mosaic/scripts/m11d_te_de_integration.R
# M11d - Integration of DE genes with TE insertions in promoter and gene body
#
# Analysis:
# 1. Intersect DE genes with promoter TE (from M11b) and gene body TE
# 2. Fisher's exact test: are DE genes enriched for promoter/body TE?
# 3. GO enrichment on DE genes WITH promoter TE vs DE genes WITHOUT
# 4. Summary tables for report
#
# Gene body = full gene locus (start to end, strand-aware)
# Promoter  = 2kb upstream of TSS (from M11b)

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

option_list <- list(
  make_option("--de-genes",       type="character", help="genes_*_DE.tsv from M08"),
  make_option("--norm-counts",    type="character", help="normalized_counts.tsv"),
  make_option("--te-prom-genes",  type="character", help="te_promoter_genes.tsv from M11b"),
  make_option("--rmsk-bed",       type="character", help="rmsk_hg38.bed for gene body intersect"),
  make_option("--gtf",            type="character", help="GENCODE GTF for gene body BED"),
  make_option("--gene-map",       type="character", help="gene_type_map.tsv"),
  make_option("--outdir",         type="character", default="mosaic_output/m11d"),
  make_option("--treatment",      type="character", default="tumor"),
  make_option("--control",        type="character", default="normal"),
  make_option("--padj",           type="double",    default=0.5,
              help="padj threshold for DE [default: 0.5 for test data]"),
  make_option("--lfc",            type="double",    default=0.5,
              help="LFC threshold [default: 0.5]"),
  make_option("--top-n",          type="integer",   default=30L,
              help="Top DE genes to highlight"),
  make_option("--organism",       type="character", default="hsa",
              help="KEGG organism code [default: hsa]")
)

opt <- parse_args(OptionParser(option_list=option_list))
dir.create(opt$outdir, recursive=TRUE, showWarnings=FALSE)

cat(sprintf("[INF] M11d - DE-TE Integration (padj<%.2f, |LFC|>%.1f)\n",
            opt$padj, opt$lfc))

# -- Helpers --------------------------------------------------------------------
safe_read <- function(path, ...) {
  if (!is.null(path) && file.exists(path))
    tryCatch(read.table(path, header=TRUE, sep="\t", check.names=FALSE, ...),
             error=function(e) { cat("[WARN]", conditionMessage(e),"\n"); NULL })
  else { cat("[WARN] Not found:", path %||% "NULL", "\n"); NULL }
}
`%||%` <- function(a,b) if(!is.null(a)) a else b
write_tsv <- function(df, path) {
  write.table(df, path, sep="\t", quote=FALSE, row.names=FALSE)
  cat("[OK] Saved:", basename(path), sprintf("(%d rows)\n", nrow(df)))
}

col_macro <- c(
  "SINE"="#3498DB", "LINE"="#E74C3C", "LTR/ERV"="#2ECC71",
  "DNA transposon"="#F39C12", "SVA"="#9B59B6",
  "Satellite"="#1ABC9C", "Other TE"="#95A5A6"
)
col_cond <- c(tumor="#E74C3C", normal="#2980B9")

classify_macro <- function(name) {
  n <- toupper(name)
  dplyr::case_when(
    grepl("^ALU|^MIR|SINE",    n) ~ "SINE",
    grepl("^L1|^L2|LINE|^CR1", n) ~ "LINE",
    grepl("^LTR|^HERV|^ERV|^MLT|^THE|^MER[0-9]", n) ~ "LTR/ERV",
    grepl("^SVA",              n) ~ "SVA",
    grepl("DNA|^TIGGER|^CHARLIE|^HAT|^MARINER|^MER[A-Z]", n) ~ "DNA transposon",
    grepl("SAT|ALPHA|CENTR|ACRO|MICRO", n) ~ "Satellite",
    TRUE ~ "Other TE"
  )
}

# -- Load data ------------------------------------------------------------------
de_raw      <- safe_read(opt[["de-genes"]])
te_prom     <- safe_read(opt[["te-prom-genes"]])
gene_map    <- safe_read(opt[["gene-map"]])
norm_counts <- safe_read(opt[["norm-counts"]], row.names=1)

if (is.null(de_raw)) { cat("[ERROR] --de-genes required\n"); quit(status=1) }

gene_names <- if (!is.null(gene_map))
  setNames(gene_map$gene_name, gene_map$gene_id) else character(0)
label_g <- function(ids) {
  nm <- gene_names[ids]; ifelse(!is.na(nm) & nm!="", nm, ids)
}

# gene_id is last column
gid_col <- colnames(de_raw)[ncol(de_raw)]
de_raw$gene_id <- de_raw[[gid_col]]

# Filter DE genes
de_sig <- de_raw |>
  filter(!is.na(padj), padj < opt$padj, abs(log2FoldChange) > opt$lfc)
cat(sprintf("[INF] DE genes (padj<%.2f, |LFC|>%.1f): %d\n",
            opt$padj, opt$lfc, nrow(de_sig)))

if (nrow(de_sig) == 0) {
  cat("[WARN] No DE genes found - check thresholds\n")
  quit(status=0)
}

# All expressed genes (background)
all_expressed <- unique(de_raw$gene_id[de_raw$baseMean > 5])
cat(sprintf("[INF] Background (baseMean>5): %d genes\n", length(all_expressed)))

# ===============================================================================
# 1. Gene body BED from GTF
# ===============================================================================
gene_body_bed <- file.path(opt$outdir, "gene_body.bed")

if (!file.exists(gene_body_bed) && !is.null(opt$gtf) && file.exists(opt$gtf)) {
  cat("[INF] Building gene body BED from GTF...\n")
  records <- list()
  con <- file(opt$gtf, "r")
  while (TRUE) {
    line <- readLines(con, n=1)
    if (length(line) == 0) break
    if (startsWith(line, "#")) next
    parts <- strsplit(line, "\t")[[1]]
    if (length(parts) < 9 || parts[3] != "gene") next
    chr    <- parts[1]
    start  <- as.integer(parts[4]) - 1L
    end    <- as.integer(parts[5])
    strand <- parts[7]
    gid    <- regmatches(parts[9],
               regexpr('(?<=gene_id ")[^"]+', parts[9], perl=TRUE))
    if (length(gid) == 0) next
    records[[length(records)+1]] <- c(chr, start, end, gid, ".", strand)
  }
  close(con)
  df_bed <- do.call(rbind, records) |> as.data.frame()
  df_bed <- df_bed[order(df_bed[[1]], as.integer(df_bed[[2]])), ]
  write.table(df_bed, gene_body_bed, sep="\t", quote=FALSE,
              row.names=FALSE, col.names=FALSE)
  cat(sprintf("[INF] Gene body BED: %d genes\n", nrow(df_bed)))
} else if (file.exists(gene_body_bed)) {
  cat("[INF] Gene body BED already exists\n")
} else {
  cat("[WARN] GTF not provided - skipping gene body intersect\n")
}

# ===============================================================================
# 2. Gene body TE intersect (bedtools)
# ===============================================================================
te_body_genes <- NULL
gene_body_intersect <- file.path(opt$outdir, "te_body_intersect.bed")

if (file.exists(gene_body_bed) && !is.null(opt[["rmsk-bed"]]) &&
    file.exists(opt[["rmsk-bed"]])) {
  cat("[INF] bedtools intersect: TE in gene body...\n")

  cmd <- paste("bedtools intersect",
               "-a", opt[["rmsk-bed"]],
               "-b", gene_body_bed,
               "-wa -wb -F 0.05",
               ">", gene_body_intersect, "2>/dev/null")
  ret <- system(cmd)
  if (ret != 0) {
    # try without -sorted
    cmd2 <- paste("bedtools intersect",
                  "-a", opt[["rmsk-bed"]],
                  "-b", gene_body_bed,
                  "-wa -wb",
                  ">", gene_body_intersect, "2>/dev/null")
    system(cmd2)
  }

  te_body_genes <- NULL  # initialize before conditional
  if (file.exists(gene_body_intersect) &&
      file.info(gene_body_intersect)$size > 0) {
    bed_body <- tryCatch(
      read.table(gene_body_intersect, header=FALSE, sep="\t",
                 stringsAsFactors=FALSE, fill=TRUE),
      error=function(e) { cat("[WARN] body intersect read failed\n"); NULL }
    )
    if (!is.null(bed_body) && nrow(bed_body) > 0) {
      # col4=TE name, col10=gene_id
      te_body_genes <- bed_body |>
        select(te_name=4, gene_id=10) |>
        mutate(macro=classify_macro(te_name)) |>
        group_by(gene_id) |>
        summarise(
          n_te_body_loci    = n(),
          n_te_body_families= n_distinct(te_name),
          te_body_macros    = paste(sort(unique(macro)), collapse=","),
          .groups="drop"
        )
      cat(sprintf("[INF] Genes with body TE: %d\n", nrow(te_body_genes)))
      write_tsv(te_body_genes,
                file.path(opt$outdir, "te_body_genes.tsv"))
    }
  }
}

# ===============================================================================
# 3. Merge DE + promoter TE + body TE
# ===============================================================================
cat("[INF] Integrating DE with TE annotations...\n")

# Prepare promoter TE summary per gene
prom_summ <- if (!is.null(te_prom)) {
  te_prom |>
    select(gene_id,
           n_te_prom_loci   = n_te_loci,
           n_te_prom_families=n_te_families,
           te_prom_macros   = te_macrocategories)
} else NULL

# Build integrated DE table
de_integrated <- de_sig |>
  select(gene_id, baseMean, log2FoldChange, lfcSE, pvalue, padj) |>
  mutate(
    gene_name  = label_g(gene_id),
    direction  = ifelse(log2FoldChange > 0, opt$treatment, opt$control)
  )

if (!is.null(prom_summ))
  de_integrated <- left_join(de_integrated, prom_summ, by="gene_id")

if (!is.null(te_body_genes))
  de_integrated <- left_join(de_integrated, te_body_genes, by="gene_id")

de_integrated <- de_integrated |>
  mutate(
    has_prom_te  = !is.na(n_te_prom_loci),
    has_body_te  = if ("n_te_body_loci" %in% colnames(de_integrated)) !is.na(n_te_body_loci) else FALSE,
    te_category  = case_when(
      has_prom_te & has_body_te ~ "Promoter + Body TE",
      has_prom_te               ~ "Promoter TE only",
      has_body_te               ~ "Body TE only",
      TRUE                      ~ "No TE"
    )
  ) |>
  arrange(padj, desc(abs(log2FoldChange)))

write_tsv(de_integrated,
          file.path(opt$outdir, "de_te_integrated.tsv"))

cat(sprintf("[INF] DE genes with promoter TE: %d / %d (%.1f%%)\n",
            sum(de_integrated$has_prom_te), nrow(de_integrated),
            100*mean(de_integrated$has_prom_te)))
cat(sprintf("[INF] DE genes with body TE:     %d / %d (%.1f%%)\n",
            sum(de_integrated$has_body_te, na.rm=TRUE), nrow(de_integrated),
            100*mean(de_integrated$has_body_te, na.rm=TRUE)))

# ===============================================================================
# 4. Fisher's exact test - enrichment
# ===============================================================================
cat("[INF] Fisher enrichment tests...\n")

fisher_results <- list()

# Promoter TE enrichment
if (!is.null(prom_summ)) {
  prom_gene_set <- unique(prom_summ$gene_id)

  # 2x2 contingency: DE/notDE x has_prom_TE/not
  de_ids    <- unique(de_sig$gene_id)
  not_de_ids<- setdiff(all_expressed, de_ids)

  a <- sum(de_ids     %in% prom_gene_set)  # DE + prom TE
  b <- sum(not_de_ids %in% prom_gene_set)  # notDE + prom TE
  c <- sum(!de_ids     %in% prom_gene_set) # DE + no prom TE
  d <- sum(!not_de_ids %in% prom_gene_set) # notDE + no prom TE

  ft_prom <- fisher.test(matrix(c(a,b,c,d), nrow=2),
                          alternative="greater")
  fisher_results[["promoter_TE"]] <- data.frame(
    analysis   = "Promoter TE in DE genes",
    n_DE_with  = a,
    n_DE_total = length(de_ids),
    n_bg_with  = b,
    n_bg_total = length(not_de_ids),
    odds_ratio = round(ft_prom$estimate, 3),
    pvalue     = ft_prom$p.value,
    CI_low     = round(ft_prom$conf.int[1], 3),
    CI_high    = round(ft_prom$conf.int[2], 3)
  )
  cat(sprintf("[INF] Promoter TE enrichment: OR=%.2f p=%.4f (%d/%d DE vs %d/%d bg)\n",
              ft_prom$estimate, ft_prom$p.value,
              a, length(de_ids), b, length(not_de_ids)))
}

# Gene body TE enrichment
if (!is.null(te_body_genes)) {
  body_gene_set <- unique(te_body_genes$gene_id)
  de_ids    <- unique(de_sig$gene_id)
  not_de_ids<- setdiff(all_expressed, de_ids)

  a <- sum(de_ids     %in% body_gene_set)
  b <- sum(not_de_ids %in% body_gene_set)
  c <- sum(!de_ids     %in% body_gene_set)
  d <- sum(!not_de_ids %in% body_gene_set)

  ft_body <- fisher.test(matrix(c(a,b,c,d), nrow=2),
                          alternative="greater")
  fisher_results[["body_TE"]] <- data.frame(
    analysis   = "Gene body TE in DE genes",
    n_DE_with  = a,
    n_DE_total = length(de_ids),
    n_bg_with  = b,
    n_bg_total = length(not_de_ids),
    odds_ratio = round(ft_body$estimate, 3),
    pvalue     = ft_body$p.value,
    CI_low     = round(ft_body$conf.int[1], 3),
    CI_high    = round(ft_body$conf.int[2], 3)
  )
  cat(sprintf("[INF] Body TE enrichment:     OR=%.2f p=%.4f\n",
              ft_body$estimate, ft_body$p.value))
}

# Per macrocategory enrichment (promoter)
if (!is.null(prom_summ) && "te_prom_macros" %in% colnames(de_integrated)) {
  de_ids     <- unique(de_sig$gene_id)
  not_de_ids <- setdiff(all_expressed, de_ids)

  macro_fisher <- lapply(
    c("SINE","LINE","LTR/ERV","DNA transposon","SVA","Satellite"),
    function(mac) {
      prom_mac <- te_prom |>
        filter(grepl(mac, te_macrocategories, fixed=TRUE)) |>
        pull(gene_id)
      a <- sum(de_ids     %in% prom_mac)
      b <- sum(not_de_ids %in% prom_mac)
      c <- length(de_ids)     - a
      d <- length(not_de_ids) - b
      if (a + b == 0) return(NULL)
      ft <- fisher.test(matrix(c(a,b,c,d), nrow=2), alternative="greater")
      data.frame(macrocategory=mac, n_DE_with=a, n_bg_with=b,
                 odds_ratio=round(ft$estimate,3), pvalue=ft$p.value)
    }
  )
  macro_fisher_df <- bind_rows(Filter(Negate(is.null), macro_fisher)) |>
    mutate(padj=p.adjust(pvalue, method="BH")) |>
    arrange(pvalue)
  write_tsv(macro_fisher_df,
            file.path(opt$outdir, "promoter_te_macro_enrichment.tsv"))
  cat("[INF] Per-macrocategory promoter enrichment saved.\n")
  fisher_results[["macro"]] <- macro_fisher_df
}

# Save all Fisher results
all_fisher <- bind_rows(
  fisher_results[["promoter_TE"]],
  fisher_results[["body_TE"]]
)
if (nrow(all_fisher) > 0)
  write_tsv(all_fisher, file.path(opt$outdir, "te_enrichment_fisher.tsv"))

# ===============================================================================
# 5. GO enrichment on DE genes WITH promoter TE
# ===============================================================================
cat("[INF] GO enrichment on DE genes with promoter TE...\n")

run_go <- function(gene_ids, universe_ids, label, outdir) {
  if (!requireNamespace("clusterProfiler", quietly=TRUE) ||
      !requireNamespace("org.Hs.eg.db",    quietly=TRUE)) {
    cat("[WARN] clusterProfiler/org.Hs.eg.db not available\n")
    return(NULL)
  }
  library(clusterProfiler); library(org.Hs.eg.db)

  # Strip version suffix (ENSG00000001234.5 -> ENSG00000001234)
  gene_ids_clean    <- sub("\\..*$", "", gene_ids)
  universe_ids_clean<- sub("\\..*$", "", universe_ids)

  # Convert ENSEMBL to ENTREZ
  eg_genes <- bitr(gene_ids_clean,     fromType="ENSEMBL",
                   toType="ENTREZID", OrgDb=org.Hs.eg.db, drop=TRUE)
  eg_univ  <- bitr(universe_ids_clean, fromType="ENSEMBL",
                   toType="ENTREZID", OrgDb=org.Hs.eg.db, drop=TRUE)

  if (nrow(eg_genes) < 3) {
    cat(sprintf("[WARN] Too few mapped genes for GO (%s): %d\n",
                label, nrow(eg_genes)))
    return(NULL)
  }

  results <- list()
  for (ont in c("BP","MF","CC")) {
    tryCatch({
      go <- enrichGO(gene          = eg_genes$ENTREZID,
                     universe      = eg_univ$ENTREZID,
                     OrgDb         = org.Hs.eg.db,
                     ont           = ont,
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.2,
                     qvalueCutoff  = 0.2,
                     readable      = TRUE)
      if (!is.null(go) && nrow(go@result) > 0) {
        res_df <- go@result |> filter(p.adjust < 0.2) |>
          arrange(p.adjust)
        write_tsv(res_df,
                  file.path(outdir, sprintf("go_%s_%s.tsv", label, ont)))
        results[[ont]] <- res_df
        cat(sprintf("[INF] GO %s (%s): %d sig terms\n",
                    ont, label, nrow(res_df)))
      }
    }, error=function(e) cat("[WARN] GO", ont, label, ":", conditionMessage(e), "\n"))
  }
  results
}

de_with_prom_te <- de_integrated |>
  filter(has_prom_te) |>
  pull(gene_id)

de_with_body_te <- de_integrated |>
  filter(has_body_te) |>
  pull(gene_id)

cat(sprintf("[INF] DE+promoter TE: %d genes for GO\n", length(de_with_prom_te)))
cat(sprintf("[INF] DE+body TE:     %d genes for GO\n", length(de_with_body_te)))

go_prom <- run_go(de_with_prom_te, all_expressed,
                  "prom_TE_DE", opt$outdir)
go_body <- run_go(de_with_body_te, all_expressed,
                  "body_TE_DE", opt$outdir)

# ===============================================================================
# 6. Visualisation plots
# ===============================================================================
cat("[INF] Generating plots...\n")

# 6a. TE category breakdown of DE genes
p_cat <- de_integrated |>
  count(direction, te_category) |>
  ggplot(aes(x=direction, y=n, fill=te_category)) +
  geom_col(position="stack", width=0.6) +
  scale_fill_manual(values=c(
    "Promoter + Body TE" = "#8E44AD",
    "Promoter TE only"   = "#3498DB",
    "Body TE only"       = "#E67E22",
    "No TE"              = "#BDC3C7"
  ), name="TE context") +
  labs(title="DE genes by TE insertion context",
       subtitle=sprintf("padj<%.2f, |LFC|>%.1f", opt$padj, opt$lfc),
       x="Direction", y="Number of DE genes") +
  theme_minimal(base_size=13) +
  theme(legend.position="right")

ggsave(file.path(opt$outdir, "de_te_category_barplot.pdf"),
       p_cat, width=7, height=5)
ggsave(file.path(opt$outdir, "de_te_category_barplot.png"),
       p_cat, width=7, height=5, dpi=150)

# 6b. Top DE genes with TE annotation - lollipop
top_de_annot <- de_integrated |>
  arrange(padj, desc(abs(log2FoldChange))) |>
  head(opt[["top-n"]]) |>
  mutate(
    label    = ifelse(gene_name!="" & !is.na(gene_name), gene_name, gene_id),
    te_label = case_when(
      has_prom_te & has_body_te ~ "** Prom+Body",
      has_prom_te               ~ ">> Promoter",
      has_body_te               ~ ">> Body",
      TRUE                      ~ "-- None"
    )
  )

p_lollipop <- ggplot(top_de_annot,
                     aes(x=reorder(label, log2FoldChange),
                         y=log2FoldChange, color=te_label)) +
  geom_segment(aes(xend=reorder(label, log2FoldChange), yend=0),
               linewidth=0.8, alpha=0.6) +
  geom_point(aes(size=-log10(padj+1e-10)), alpha=0.9) +
  geom_hline(yintercept=0, linetype="dashed", color="grey50") +
  coord_flip() +
  scale_color_manual(values=c(
    "** Prom+Body" = "#8E44AD",
    ">> Promoter"  = "#3498DB",
    ">> Body"      = "#E67E22",
    "-- None"      = "#95A5A6"
  ), name="TE context") +
  scale_size_continuous(range=c(2,7), name="-log10(padj)") +
  labs(
    title    = sprintf("Top %d DE genes - TE insertion context", opt[["top-n"]]),
    subtitle = sprintf("%s vs %s | padj<%.2f | |LFC|>%.1f",
                       opt$treatment, opt$control, opt$padj, opt$lfc),
    x        = NULL,
    y        = sprintf("log2FC (%s / %s)", opt$treatment, opt$control)
  ) +
  theme_minimal(base_size=11) +
  theme(legend.position="right")

h <- max(6, opt[["top-n"]] * 0.22 + 1.5)
ggsave(file.path(opt$outdir, "top_de_te_lollipop.pdf"),
       p_lollipop, width=10, height=h)
ggsave(file.path(opt$outdir, "top_de_te_lollipop.png"),
       p_lollipop, width=10, height=h, dpi=150)

# 6c. Macrocategory enrichment barplot
if (!is.null(fisher_results[["macro"]])) {
  p_macro_enrich <- fisher_results[["macro"]] |>
    mutate(sig=ifelse(padj<0.05,"*",""),
           macrocategory=reorder(macrocategory, -log10(pvalue+1e-10))) |>
    ggplot(aes(x=macrocategory, y=-log10(pvalue), fill=macrocategory)) +
    geom_col(show.legend=FALSE) +
    geom_text(aes(label=sig), vjust=-0.3, size=5) +
    geom_hline(yintercept=-log10(0.05), linetype="dashed",
               color="grey40", linewidth=0.6) +
    scale_fill_manual(values=col_macro) +
    labs(title="Promoter TE macrocategory enrichment in DE genes",
         subtitle="Fisher's exact test (one-sided); dashed = p=0.05; * = BH adj p<0.05",
         x=NULL, y="-log10(p-value)") +
    theme_minimal(base_size=12) +
    theme(axis.text.x=element_text(angle=30, hjust=1))

  ggsave(file.path(opt$outdir, "promoter_macro_enrichment.pdf"),
         p_macro_enrich, width=8, height=5)
  ggsave(file.path(opt$outdir, "promoter_macro_enrichment.png"),
         p_macro_enrich, width=8, height=5, dpi=150)
}

# -- Summary --------------------------------------------------------------------
cat("\n[INF] M11d Summary:\n")
cat(sprintf("  DE genes analysed:          %d\n", nrow(de_integrated)))
cat(sprintf("  With promoter TE:           %d (%.1f%%)\n",
            sum(de_integrated$has_prom_te),
            100*mean(de_integrated$has_prom_te)))
cat(sprintf("  With gene body TE:          %d (%.1f%%)\n",
            sum(de_integrated$has_body_te, na.rm=TRUE),
            100*mean(de_integrated$has_body_te, na.rm=TRUE)))
if (!is.null(fisher_results[["promoter_TE"]]))
  cat(sprintf("  Promoter enrichment OR:     %.2f (p=%.4f)\n",
              fisher_results[["promoter_TE"]]$odds_ratio,
              fisher_results[["promoter_TE"]]$pvalue))

cat("[OK] M11d complete. Output:", opt$outdir, "\n")
