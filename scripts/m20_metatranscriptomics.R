#!/usr/bin/env Rscript
# mosaic/scripts/m20_metatranscriptomics.R
# M20 - Metatranscriptomics: Kraken2/Bracken taxonomic profiling
#
# Input:  Kraken2 reports (*_kraken2_report.txt)
#         Bracken species tables (*_bracken_species.txt)
#         coldata.tsv
# Output: combined abundance matrix, diversity metrics,
#         differential abundance, plots (barplot, heatmap, PCoA)
#
# NOTE: Homo sapiens reads dominate unmapped fraction.
# Analysis excludes Homo sapiens and focuses on microbial taxa.
# With full data, filter: --min-reads 10, --min-samples 2

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

option_list <- list(
  make_option("--meta-dir",   type="character", help="M20 metatranscriptomics dir"),
  make_option("--coldata",    type="character", help="coldata.tsv"),
  make_option("--outdir",     type="character", default=NULL),
  make_option("--treatment",  type="character", default="tumor"),
  make_option("--control",    type="character", default="normal"),
  make_option("--min-reads",  type="integer",   default=2L,
              help="Min reads to keep taxon [2]"),
  make_option("--min-samples",type="integer",   default=1L,
              help="Min samples with reads to keep taxon [1]"),
  make_option("--top-n",      type="integer",   default=20L,
              help="Top N taxa to show in plots [20]"),
  make_option("--level",      type="character", default="S",
              help="Taxonomic level: S=species, G=genus [S]")
)

opt    <- parse_args(OptionParser(option_list=option_list))
outdir <- if (!is.null(opt$outdir)) { opt$outdir } else { opt[["meta-dir"]] }
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)
`%||%` <- function(a,b) { if (!is.null(a)) a else b }

safe_read <- function(path, ...) {
  if (!is.null(path) && file.exists(path)) {
    tryCatch(read.table(path, header=TRUE, sep="\t",
                        check.names=FALSE, stringsAsFactors=FALSE, ...),
             error=function(e) NULL)
  } else {
    NULL
  }
}
write_tsv <- function(df, path) {
  write.table(df, path, sep="\t", quote=FALSE, row.names=FALSE)
  cat("[OK]", basename(path), sprintf("(%d rows)\n", nrow(df)))
}

cat("[INF] M20 - Metatranscriptomics analysis\n")

coldata <- safe_read(opt$coldata)
if (is.null(coldata)) { cat("[ERROR] --coldata required\n"); quit(status=1) }

# ── 1. Load Bracken species tables ────────────────────────────────────────────
bracken_files <- list.files(opt[["meta-dir"]],
                             pattern="_bracken_species\\.txt$",
                             full.names=TRUE, recursive=TRUE)

if (length(bracken_files) == 0) {
  cat("[ERROR] No Bracken files found in", opt[["meta-dir"]], "\n")
  quit(status=1)
}

cat(sprintf("[INF] Found %d Bracken files\n", length(bracken_files)))

bracken_list <- lapply(bracken_files, function(f) {
  sid <- sub("_bracken_species\\.txt$", "", basename(f))
  df  <- safe_read(f)
  if (is.null(df)) return(NULL)
  # Bracken cols: name, taxonomy_id, taxonomy_lvl, kraken_assigned_reads,
  #               added_reads, new_est_reads, fraction_total_reads
  df$sample_id <- sid
  df
})
bracken_list <- Filter(Negate(is.null), bracken_list)
all_bracken  <- bind_rows(bracken_list)

# ── 2. Load Kraken2 reports for summary stats ─────────────────────────────────
kraken_files <- list.files(opt[["meta-dir"]],
                            pattern="_kraken2_report\\.txt$",
                            full.names=TRUE, recursive=TRUE)

kraken_summary <- bind_rows(lapply(kraken_files, function(f) {
  sid <- sub("_kraken2_report\\.txt$", "", basename(f))
  df  <- tryCatch(read.table(f, header=FALSE, sep="\t",
                              stringsAsFactors=FALSE, fill=TRUE),
                  error=function(e) NULL)
  if (is.null(df)) return(NULL)
  # cols: pct, clade_reads, direct_reads, rank, taxid, name
  total      <- sum(df[[3]], na.rm=TRUE)
  classified <- sum(df[[3]][df[[5]] != 0], na.rm=TRUE)
  unclass    <- sum(df[[3]][df[[5]] == 0], na.rm=TRUE)
  data.frame(
    sample_id       = sid,
    total_reads     = total,
    classified      = classified,
    unclassified    = unclass,
    pct_classified  = round(100 * classified / max(total, 1), 2)
  )
}))

if (!is.null(kraken_summary) && nrow(kraken_summary) > 0) {
  kraken_summary <- left_join(kraken_summary,
                               coldata[,c("sample_id","condition")],
                               by="sample_id")
  write_tsv(kraken_summary,
            file.path(outdir, "kraken2_classification_summary.tsv"))
  cat("[INF] Classification rates:\n")
  print(kraken_summary[,c("sample_id","condition","total_reads","pct_classified")],
        row.names=FALSE)
}

# ── 3. Build abundance matrix (exclude Homo sapiens) ─────────────────────────
# Use new_est_reads (Bracken-corrected counts)
abund_long <- all_bracken |>
  filter(
    !grepl("Homo sapiens|homo sapiens", name, ignore.case=TRUE),
    !grepl("^name$", name),           # skip header artefacts
    new_est_reads >= opt[["min-reads"]]
  ) |>
  select(name, sample_id, new_est_reads)

# Wide matrix
abund_wide <- abund_long |>
  pivot_wider(names_from=sample_id,
              values_from=new_est_reads,
              values_fill=0)

# Filter: taxon must appear in min-samples
sample_cols <- coldata$sample_id[coldata$sample_id %in% colnames(abund_wide)]
n_present   <- rowSums(abund_wide[, sample_cols, drop=FALSE] > 0)
abund_wide  <- abund_wide[n_present >= opt[["min-samples"]], ]

cat(sprintf("[INF] Taxa after filtering: %d\n", nrow(abund_wide)))
write_tsv(abund_wide,
          file.path(outdir, "microbial_abundance_matrix.tsv"))

if (nrow(abund_wide) == 0) {
  cat("[WARN] No microbial taxa after filtering\n")
  cat("[INF] Generating placeholder output\n")

  p <- ggplot() +
    annotate("text", x=0.5, y=0.5,
             label=paste0("No microbial taxa detected above threshold.\n",
                          "With full-depth data, microbial signals\n",
                          "will appear here (viruses, bacteria, fungi)."),
             size=4, hjust=0.5) +
    theme_void() +
    labs(title="Metatranscriptomics (M20) — placeholder")
  ggsave(file.path(outdir, "metatranscriptomics_placeholder.pdf"),
         p, width=7, height=4)
  ggsave(file.path(outdir, "metatranscriptomics_placeholder.png"),
         p, width=7, height=4, dpi=150)
  quit(status=0)
}

# ── 4. Relative abundance per sample ─────────────────────────────────────────
abund_mat <- as.matrix(abund_wide[, sample_cols])
rownames(abund_mat) <- abund_wide$name

# Normalise to relative abundance within microbial reads
rel_abund <- sweep(abund_mat, 2, colSums(abund_mat), "/") * 100

rel_df <- as.data.frame(rel_abund) |>
  tibble::rownames_to_column("taxon") |>
  pivot_longer(-taxon, names_to="sample_id", values_to="rel_abund") |>
  left_join(coldata[,c("sample_id","condition")], by="sample_id")

# ── 5. Stacked barplot — top N taxa ──────────────────────────────────────────
top_taxa <- rel_df |>
  group_by(taxon) |>
  summarise(mean_abund=mean(rel_abund), .groups="drop") |>
  arrange(desc(mean_abund)) |>
  head(opt[["top-n"]]) |>
  pull(taxon)

# colour palette
n_taxa <- length(top_taxa)
tax_cols <- if (n_taxa <= 12) {
  RColorBrewer::brewer.pal(max(3, n_taxa), "Paired")[1:n_taxa]
} else {
  colorRampPalette(RColorBrewer::brewer.pal(12,"Paired"))(n_taxa)
}
names(tax_cols) <- top_taxa

plot_df <- rel_df |>
  mutate(taxon=ifelse(taxon %in% top_taxa, taxon, "Other")) |>
  group_by(sample_id, condition, taxon) |>
  summarise(rel_abund=sum(rel_abund), .groups="drop")

p_bar <- ggplot(plot_df,
                aes(x=sample_id, y=rel_abund, fill=taxon)) +
  geom_col(position="stack") +
  facet_wrap(~condition, scales="free_x", nrow=1) +
  scale_fill_manual(
    values=c(tax_cols, "Other"="#D3D3D3"),
    name="Taxon"
  ) +
  labs(
    title=sprintf("Microbial composition — Top %d species", opt[["top-n"]]),
    subtitle="Relative abundance (% of microbial reads, Homo sapiens excluded)",
    x=NULL, y="Relative abundance (%)"
  ) +
  theme_minimal(base_size=12) +
  theme(axis.text.x=element_text(angle=30, hjust=1),
        legend.text=element_text(size=8))

ggsave(file.path(outdir, "microbial_barplot.pdf"),   p_bar, width=10, height=6)
ggsave(file.path(outdir, "microbial_barplot.png"),   p_bar, width=10, height=6, dpi=150)

# ── 6. Differential abundance (simple fold-change, no stats with n=2) ─────────
cond_mean <- rel_df |>
  group_by(taxon, condition) |>
  summarise(mean_ra=mean(rel_abund), .groups="drop") |>
  pivot_wider(names_from=condition,
              values_from=mean_ra, values_fill=0)

trt <- opt$treatment
ctl <- opt$control

if (trt %in% colnames(cond_mean) && ctl %in% colnames(cond_mean)) {
  cond_mean <- cond_mean |>
    mutate(
      log2FC = log2((!!sym(trt) + 0.01) / (!!sym(ctl) + 0.01)),
      mean_abund = (!!sym(trt) + !!sym(ctl)) / 2
    ) |>
    arrange(desc(abs(log2FC)))

  write_tsv(cond_mean,
            file.path(outdir, "microbial_differential_abundance.tsv"))

  # Lollipop
  top_diff <- cond_mean |>
    filter(mean_abund > 0) |>
    arrange(desc(abs(log2FC))) |>
    head(opt[["top-n"]])

  if (nrow(top_diff) > 0) {
    p_lollipop <- ggplot(top_diff,
                         aes(x=reorder(taxon, log2FC),
                             y=log2FC,
                             color=log2FC > 0)) +
      geom_segment(aes(xend=reorder(taxon, log2FC), yend=0),
                   linewidth=0.8, alpha=0.7) +
      geom_point(aes(size=mean_abund), alpha=0.9) +
      geom_hline(yintercept=0, linetype="dashed", color="grey50") +
      coord_flip() +
      scale_color_manual(
        values=c("TRUE"="#E74C3C","FALSE"="#3498DB"),
        labels=c("TRUE"=paste0("Higher in ",trt),
                 "FALSE"=paste0("Higher in ",ctl)),
        name=NULL
      ) +
      scale_size_continuous(range=c(2,8), name="Mean RA (%)") +
      labs(
        title="Differential microbial abundance",
        subtitle=sprintf("%s vs %s | log2FC (note: n=2 per group, no FDR)",
                         trt, ctl),
        x=NULL,
        y=sprintf("log2FC (%s / %s)", trt, ctl)
      ) +
      theme_minimal(base_size=11)

    h <- max(5, nrow(top_diff) * 0.22 + 1.5)
    ggsave(file.path(outdir, "microbial_differential_lollipop.pdf"),
           p_lollipop, width=10, height=h)
    ggsave(file.path(outdir, "microbial_differential_lollipop.png"),
           p_lollipop, width=10, height=h, dpi=150)
  }
}

# ── 7. Alpha diversity ─────────────────────────────────────────────────────────
shannon <- function(x) {
  p <- x[x > 0] / sum(x)
  -sum(p * log(p))
}
richness <- function(x) sum(x > 0)

div_df <- data.frame(
  sample_id = sample_cols,
  shannon   = apply(abund_mat[, sample_cols, drop=FALSE], 2, shannon),
  richness  = apply(abund_mat[, sample_cols, drop=FALSE], 2, richness)
) |>
  left_join(coldata[,c("sample_id","condition")], by="sample_id")

write_tsv(div_df, file.path(outdir, "alpha_diversity.tsv"))

p_div <- div_df |>
  pivot_longer(c(shannon, richness), names_to="metric", values_to="value") |>
  ggplot(aes(x=condition, y=value, color=condition, label=sample_id)) +
  geom_point(size=4, alpha=0.8) +
  ggrepel::geom_text_repel(size=3, show.legend=FALSE) +
  facet_wrap(~metric, scales="free_y",
             labeller=labeller(metric=c(shannon="Shannon index",
                                        richness="Species richness"))) +
  scale_color_manual(values=c(setNames(c("#2ECC71","#E74C3C"),
                                        c(ctl, trt))), name=NULL) +
  labs(title="Microbial alpha diversity",
       subtitle="Shannon index and species richness per sample",
       x=NULL, y=NULL) +
  theme_minimal(base_size=12) +
  theme(legend.position="top")

ggsave(file.path(outdir, "alpha_diversity.pdf"), p_div, width=8, height=5)
ggsave(file.path(outdir, "alpha_diversity.png"), p_div, width=8, height=5, dpi=150)

# ── 8. Summary ────────────────────────────────────────────────────────────────
cat(sprintf("\n[INF] M20 Summary:\n"))
cat(sprintf("  Microbial taxa detected:  %d\n", nrow(abund_wide)))
cat(sprintf("  Samples:                  %d\n", length(sample_cols)))
cat(sprintf("  Mean Shannon diversity:   %.3f\n", mean(div_df$shannon)))
if (exists("cond_mean"))
  cat(sprintf("  Top enriched in %s:    %s\n", trt,
              paste(head(cond_mean$taxon[cond_mean$log2FC > 0], 3),
                    collapse=", ")))

cat("[OK] M20 complete. Output:", outdir, "\n")
