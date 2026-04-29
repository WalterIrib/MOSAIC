#!/usr/bin/env Rscript
# mosaic/scripts/m19_editing_plots.R
# M19 — RNA editing visualisation and statistical analysis
#
# Input:  editing_sites_ato_i.tsv (from m19_rna_editing.py)
#         editing_per_gene.tsv
#         editing_summary.tsv
# Output: dEF distribution, scatter normal vs tumor EF,
#         per-gene barplot, GO enrichment on edited genes

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

option_list <- list(
  make_option("--editing-dir", type="character", help="M19 output directory"),
  make_option("--outdir",      type="character", default=NULL),
  make_option("--treatment",   type="character", default="tumor"),
  make_option("--control",     type="character", default="normal"),
  make_option("--min-def",     type="double",    default=0.05,
              help="|dEF| threshold for differential editing [0.05]"),
  make_option("--top-n",       type="integer",   default=30L)
)

opt     <- parse_args(OptionParser(option_list=option_list))
outdir  <- opt$outdir %||% opt[["editing-dir"]]
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

`%||%` <- function(a,b) if(!is.null(a)) a else b

safe_read <- function(path, ...) {
  if (!is.null(path) && file.exists(path))
    tryCatch(read.table(path, header=TRUE, sep="\t",
                        check.names=FALSE, stringsAsFactors=FALSE, ...),
             error=function(e) NULL)
  else NULL
}

write_tsv <- function(df, path) {
  write.table(df, path, sep="\t", quote=FALSE, row.names=FALSE)
  cat("[OK]", basename(path), sprintf("(%d rows)\n", nrow(df)))
}

cat("[INF] M19 editing visualisation\n")

sites     <- safe_read(file.path(opt[["editing-dir"]], "editing_sites_ato_i.tsv"))
per_gene  <- safe_read(file.path(opt[["editing-dir"]], "editing_per_gene.tsv"))
summary_t <- safe_read(file.path(opt[["editing-dir"]], "editing_summary.tsv"))

if (is.null(sites) || nrow(sites) == 0) {
  cat("[WARN] No editing sites found — generating placeholder summary\n")

  # Placeholder summary plot
  p <- ggplot() +
    annotate("text", x=0.5, y=0.5,
             label="No A-to-I editing sites detected\nat current thresholds.\nWith full-depth data, sites will appear here.",
             size=5, hjust=0.5) +
    theme_void() +
    labs(title="RNA editing (A-to-I) — M19")
  ggsave(file.path(outdir, "editing_placeholder.pdf"), p, width=7, height=4)
  ggsave(file.path(outdir, "editing_placeholder.png"), p, width=7, height=4, dpi=150)
  cat("[OK] M19 plots complete (placeholder)\n")
  quit(status=0)
}

cat(sprintf("[INF] Editing sites: %d\n", nrow(sites)))

# ── 1. dEF distribution ────────────────────────────────────────────────────────
p_def <- ggplot(sites, aes(x=dEF)) +
  geom_histogram(bins=60, fill="#3498DB", color="white", alpha=0.8) +
  geom_vline(xintercept=c(-opt[["min-def"]], opt[["min-def"]]),
             linetype="dashed", color="#E74C3C", linewidth=0.8) +
  geom_vline(xintercept=0, linewidth=0.5, color="grey30") +
  labs(
    title="Distribution of differential editing frequency (dEF)",
    subtitle=sprintf("dEF = EF(%s) - EF(%s) | dashed = +/-%.2f threshold | n=%d sites",
                     opt$treatment, opt$control, opt[["min-def"]], nrow(sites)),
    x=sprintf("dEF (%s - %s)", opt$treatment, opt$control),
    y="Number of sites"
  ) +
  theme_minimal(base_size=13)

ggsave(file.path(outdir, "def_distribution.pdf"), p_def, width=9, height=5)
ggsave(file.path(outdir, "def_distribution.png"), p_def, width=9, height=5, dpi=150)

# ── 2. Scatter: EF normal vs EF tumor ─────────────────────────────────────────
sites_plot <- sites |>
  mutate(
    diff_edit = case_when(
      dEF >  opt[["min-def"]] ~ paste0("Increased in ", opt$treatment),
      dEF < -opt[["min-def"]] ~ paste0("Decreased in ", opt$treatment),
      TRUE                     ~ "Similar"
    )
  )

n_inc <- sum(sites_plot$diff_edit == paste0("Increased in ", opt$treatment))
n_dec <- sum(sites_plot$diff_edit == paste0("Decreased in ", opt$treatment))

p_scatter <- ggplot(sites_plot,
                    aes(x=EF1, y=EF2, color=diff_edit)) +
  geom_point(alpha=0.5, size=1.2) +
  geom_abline(slope=1, intercept=0, linetype="dashed", color="grey50") +
  scale_color_manual(
    values=c("#E74C3C","#3498DB","#BDC3C7"),
    name="Editing change"
  ) +
  labs(
    title=sprintf("Editing frequency: %s vs %s", opt$treatment, opt$control),
    subtitle=sprintf("Increased: %d | Decreased: %d | Similar: %d",
                     n_inc, n_dec,
                     nrow(sites_plot) - n_inc - n_dec),
    x=sprintf("EF (%s)", opt$control),
    y=sprintf("EF (%s)", opt$treatment)
  ) +
  theme_minimal(base_size=13)

ggsave(file.path(outdir, "ef_scatter.pdf"),   p_scatter, width=8, height=7)
ggsave(file.path(outdir, "ef_scatter.png"),   p_scatter, width=8, height=7, dpi=150)

# ── 3. Top genes by editing sites ─────────────────────────────────────────────
if (!is.null(per_gene) && nrow(per_gene) > 0) {
  top_genes <- per_gene |>
    arrange(desc(n_sites)) |>
    head(opt[["top-n"]]) |>
    mutate(
      gene_name = reorder(gene_name, n_sites),
      edit_direction = case_when(
        mean_dEF >  opt[["min-def"]] ~ paste0("More in ", opt$treatment),
        mean_dEF < -opt[["min-def"]] ~ paste0("More in ", opt$control),
        TRUE                          ~ "Similar"
      )
    )

  p_genes <- ggplot(top_genes,
                    aes(x=gene_name, y=n_sites, fill=edit_direction)) +
    geom_col() +
    coord_flip() +
    scale_fill_manual(
      values=c("#E74C3C","#3498DB","#BDC3C7"),
      name="Direction"
    ) +
    labs(
      title=sprintf("Top %d genes by A-to-I editing site count", opt[["top-n"]]),
      x=NULL, y="Number of editing sites"
    ) +
    theme_minimal(base_size=11)

  h <- max(5, opt[["top-n"]] * 0.22 + 1.5)
  ggsave(file.path(outdir, "top_genes_editing.pdf"), p_genes, width=9, height=h)
  ggsave(file.path(outdir, "top_genes_editing.png"), p_genes, width=9, height=h, dpi=150)
}

# ── 4. dEF vs score volcano-style ─────────────────────────────────────────────
if ("gene_name" %in% colnames(sites)) {
  top_labels <- sites |>
    filter(!is.na(gene_name)) |>
    arrange(desc(abs(dEF)), desc(score)) |>
    head(20)

  p_vol <- ggplot(sites_plot,
                  aes(x=dEF, y=score, color=diff_edit)) +
    geom_point(alpha=0.5, size=1.5) +
    ggrepel::geom_text_repel(
      data=top_labels,
      aes(x=dEF, y=score, label=gene_name),
      size=2.8, color="black", max.overlaps=15, inherit.aes=FALSE
    ) +
    geom_vline(xintercept=c(-opt[["min-def"]], opt[["min-def"]]),
               linetype="dashed", color="grey50") +
    scale_color_manual(
      values=c("#E74C3C","#3498DB","#BDC3C7"),
      name="Editing change"
    ) +
    labs(
      title="RNA editing sites - dEF vs JACUSA2 score",
      subtitle="Labels = top 20 sites by |dEF| in annotated genes",
      x=sprintf("dEF (%s - %s)", opt$treatment, opt$control),
      y="JACUSA2 score"
    ) +
    theme_minimal(base_size=12)

  ggsave(file.path(outdir, "editing_volcano.pdf"), p_vol, width=9, height=6)
  ggsave(file.path(outdir, "editing_volcano.png"), p_vol, width=9, height=6, dpi=150)
}

# ── 5. GO enrichment on differentially edited genes ───────────────────────────
if ("gene_name" %in% colnames(sites) &&
    requireNamespace("clusterProfiler", quietly=TRUE) &&
    requireNamespace("org.Hs.eg.db", quietly=TRUE)) {

  library(clusterProfiler); library(org.Hs.eg.db)

  diff_genes <- sites |>
    filter(abs(dEF) >= opt[["min-def"]], !is.na(gene_name)) |>
    pull(gene_name) |> unique()

  bg_genes <- sites |>
    filter(!is.na(gene_name)) |>
    pull(gene_name) |> unique()

  cat(sprintf("[INF] GO on differentially edited genes: %d / %d bg\n",
              length(diff_genes), length(bg_genes)))

  # Convert to ENTREZ via gene symbol
  eg_diff <- tryCatch(
    bitr(diff_genes, fromType="SYMBOL", toType="ENTREZID",
         OrgDb=org.Hs.eg.db, drop=TRUE),
    error=function(e) NULL)
  eg_bg <- tryCatch(
    bitr(bg_genes, fromType="SYMBOL", toType="ENTREZID",
         OrgDb=org.Hs.eg.db, drop=TRUE),
    error=function(e) NULL)

  if (!is.null(eg_diff) && nrow(eg_diff) >= 3) {
    for (ont in c("BP","MF")) {
      tryCatch({
        go <- enrichGO(
          gene          = eg_diff$ENTREZID,
          universe      = eg_bg$ENTREZID,
          OrgDb         = org.Hs.eg.db,
          ont           = ont,
          pAdjustMethod = "BH",
          pvalueCutoff  = 0.2,
          readable      = TRUE
        )
        if (!is.null(go) && nrow(go@result) > 0) {
          res <- go@result |> filter(p.adjust < 0.2) |> arrange(p.adjust)
          write_tsv(res, file.path(outdir, sprintf("go_edited_genes_%s.tsv", ont)))
          cat(sprintf("[INF] GO %s: %d terms\n", ont, nrow(res)))
        }
      }, error=function(e) cat("[WARN] GO", ont, ":", conditionMessage(e), "\n"))
    }
  } else {
    cat("[WARN] Too few differentially edited genes for GO\n")
  }
}

# ── 6. Summary table ──────────────────────────────────────────────────────────
if (!is.null(summary_t)) {
  cat("\n[INF] Editing summary:\n")
  print(summary_t, row.names=FALSE)
}

# ── 7. All editing types (C→U/APOBEC, G→A, etc.) ─────────────────────────────
all_edits <- safe_read(file.path(opt[["editing-dir"]], "editing_sites_all.tsv"))

if (!is.null(all_edits) && nrow(all_edits) > 0 && "mutation_type" %in% colnames(all_edits)) {
  cat(sprintf("[INF] All-type editing sites: %d across %d types\n",
              nrow(all_edits),
              length(unique(all_edits$mutation_type))))

  # canonical 6-substitution grouping (pyrimidine base)
  canonical_map <- c(
    "A>G"="A>G (A-to-I)",  "T>C"="A>G (A-to-I)",   # ADAR
    "C>T"="C>T (APOBEC)",  "G>A"="C>T (APOBEC)",    # APOBEC3/AID
    "C>G"="C>G",           "G>C"="C>G",
    "C>A"="C>A",           "G>T"="C>A",
    "T>A"="T>A",           "A>T"="T>A",
    "T>G"="T>G",           "A>C"="T>G"
  )
  type_colors <- c(
    "A>G (A-to-I)"  = "#1F77B4",
    "C>T (APOBEC)"  = "#FF7F0E",
    "C>G"           = "#2CA02C",
    "C>A"           = "#D62728",
    "T>A"           = "#9467BD",
    "T>G"           = "#8C564B"
  )

  all_edits <- all_edits |>
    mutate(canonical = canonical_map[mutation_type],
           canonical = ifelse(is.na(canonical), mutation_type, canonical))

  # 7a. Mutation spectrum barplot
  type_counts <- all_edits |>
    count(mutation_type, canonical) |>
    arrange(desc(n))

  p_spectrum <- ggplot(type_counts,
                       aes(x=reorder(mutation_type, n), y=n, fill=canonical)) +
    geom_col(width=0.7) +
    coord_flip() +
    scale_fill_manual(values=type_colors, name="Substitution class") +
    labs(
      title="RNA editing mutation spectrum (all types)",
      subtitle=sprintf("Total sites: %s | %s vs %s",
                       format(nrow(all_edits), big.mark=","),
                       opt$treatment, opt$control),
      x="Mutation type", y="Number of sites"
    ) +
    theme_minimal(base_size=13) +
    theme(legend.position="right")

  ggsave(file.path(outdir, "editing_spectrum.pdf"), p_spectrum, width=8, height=6)
  ggsave(file.path(outdir, "editing_spectrum.png"), p_spectrum, width=8, height=6, dpi=150)
  cat("[OK] editing_spectrum.png\n")

  # 7b. dEF violin per canonical class (classes with >= 5 sites)
  min_for_violin <- 5L
  cls_ok <- all_edits |>
    count(canonical) |>
    filter(n >= min_for_violin) |>
    pull(canonical)

  if (length(cls_ok) >= 2) {
    violin_df <- all_edits |>
      filter(canonical %in% cls_ok) |>
      mutate(canonical = factor(canonical, levels=names(type_colors)[names(type_colors) %in% cls_ok]))

    p_violin <- ggplot(violin_df, aes(x=canonical, y=dEF, fill=canonical)) +
      geom_violin(trim=TRUE, alpha=0.7, scale="width") +
      geom_boxplot(width=0.12, fill="white", outlier.shape=NA) +
      geom_hline(yintercept=0, linetype="dashed", color="grey40") +
      scale_fill_manual(values=type_colors, name=NULL) +
      labs(
        title="Differential editing frequency (dEF) by substitution class",
        subtitle=sprintf("dEF = EF(%s) − EF(%s) | only classes with ≥%d sites",
                         opt$treatment, opt$control, min_for_violin),
        x=NULL, y=sprintf("dEF (%s − %s)", opt$treatment, opt$control)
      ) +
      theme_minimal(base_size=13) +
      theme(legend.position="none",
            axis.text.x=element_text(angle=25, hjust=1))

    ggsave(file.path(outdir, "editing_dEF_by_class.pdf"), p_violin, width=9, height=6)
    ggsave(file.path(outdir, "editing_dEF_by_class.png"), p_violin, width=9, height=6, dpi=150)
    cat("[OK] editing_dEF_by_class.png\n")
  }

  # 7c. Summary TSV per canonical class
  class_summary <- all_edits |>
    group_by(canonical, mutation_type) |>
    summarise(
      n_sites        = n(),
      mean_EF_ctrl   = round(mean(EF1, na.rm=TRUE), 4),
      mean_EF_treat  = round(mean(EF2, na.rm=TRUE), 4),
      mean_dEF       = round(mean(dEF, na.rm=TRUE), 4),
      n_increased    = sum(dEF >  opt[["min-def"]], na.rm=TRUE),
      n_decreased    = sum(dEF < -opt[["min-def"]], na.rm=TRUE),
      .groups="drop"
    ) |>
    arrange(desc(n_sites))

  write_tsv(class_summary,
            file.path(outdir, "editing_type_summary.tsv"))
  cat("[OK] editing_type_summary.tsv\n")

} else {
  cat("[WARN] editing_sites_all.tsv not found or empty — skipping all-types plots\n")
}

cat("[OK] M19 editing plots complete. Output:", outdir, "\n")
