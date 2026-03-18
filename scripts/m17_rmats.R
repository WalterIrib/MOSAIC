#!/usr/bin/env Rscript
# mosaic/scripts/m17_rmats.R
# M17 - rMATS differential alternative splicing analysis
#
# Input:  rMATS output directory (*.MATS.JC.txt files)
# Output: summary tables, volcano plots, sashimi-ready tables
#
# Event types:
#   SE  = Skipped Exon
#   A3SS= Alternative 3' Splice Site
#   A5SS= Alternative 5' Splice Site
#   MXE = Mutually Exclusive Exons
#   RI  = Retained Intron
#
# NOTE: rMATS was run with --statoff (2 replicates per condition).
# With full data (>=3 replicates), re-run without --statoff for
# proper FDR-controlled statistics.

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

option_list <- list(
  make_option("--rmats-dir",  type="character", help="rMATS output directory"),
  make_option("--outdir",     type="character", default="mosaic_output/m17_rmats"),
  make_option("--treatment",  type="character", default="tumor"),
  make_option("--control",    type="character", default="normal"),
  make_option("--fdr",        type="double",    default=0.05,
              help="FDR threshold [default 0.05]"),
  make_option("--dpsi",       type="double",    default=0.1,
              help="|deltaPSI| threshold [default 0.1]"),
  make_option("--min-reads",  type="integer",   default=10L,
              help="Min junction reads to keep event [default 10]"),
  make_option("--top-n",      type="integer",   default=30L,
              help="Top N events per type to show"),
  make_option("--statoff",    action="store_true", default=FALSE,
              help="rMATS was run with --statoff (no FDR available)")
)

opt <- parse_args(OptionParser(option_list=option_list))
dir.create(opt$outdir, recursive=TRUE, showWarnings=FALSE)

cat("[INF] M17 - rMATS differential splicing\n")
cat(sprintf("[INF] Thresholds: FDR<%.2f, |dPSI|>%.2f, min_reads>=%d\n",
            opt$fdr, opt$dpsi, opt[["min-reads"]]))

write_tsv <- function(df, path) {
  write.table(df, path, sep="\t", quote=FALSE, row.names=FALSE)
  cat("[OK] Saved:", basename(path), sprintf("(%d rows)\n", nrow(df)))
}

event_colors <- c(
  SE   = "#E74C3C",
  A3SS = "#3498DB",
  A5SS = "#2ECC71",
  MXE  = "#F39C12",
  RI   = "#9B59B6"
)

event_labels <- c(
  SE   = "Skipped Exon",
  A3SS = "Alt 3' Splice Site",
  A5SS = "Alt 5' Splice Site",
  MXE  = "Mutually Exclusive Exons",
  RI   = "Retained Intron"
)

# ── 1. Load all event types ────────────────────────────────────────────────────
load_rmats <- function(rmats_dir, event_type, count_type="JC") {
  fname <- file.path(rmats_dir,
                     sprintf("%s.MATS.%s.txt", event_type, count_type))
  if (!file.exists(fname)) {
    cat(sprintf("[WARN] Not found: %s\n", basename(fname)))
    return(NULL)
  }
  df <- tryCatch(
    read.table(fname, header=TRUE, sep="\t", stringsAsFactors=FALSE,
               fill=TRUE, quote=""),
    error=function(e) { cat("[WARN]", conditionMessage(e), "\n"); NULL }
  )
  if (is.null(df)) return(NULL)

  # Standardise column names
  df$event_type <- event_type
  df$geneSymbol <- gsub('"','', df$geneSymbol)
  df$GeneID     <- gsub('"','', df$GeneID)

  # Compute deltaPSI = IncLevel1_mean - IncLevel2_mean
  # b1=control(normal), b2=treatment(tumor) in rMATS
  # so dPSI = PSI_tumor - PSI_normal
  parse_levels <- function(x) {
    sapply(strsplit(x, ","), function(v) mean(as.numeric(v), na.rm=TRUE))
  }

  if ("IncLevel1" %in% colnames(df) && "IncLevel2" %in% colnames(df)) {
    df$psi_control   <- parse_levels(df$IncLevel1)
    df$psi_treatment <- parse_levels(df$IncLevel2)
    df$deltaPSI      <- df$psi_treatment - df$psi_control
  }

  # Junction read counts
  if ("IJC_SAMPLE_1" %in% colnames(df)) {
    parse_counts <- function(x) {
      sapply(strsplit(x, ","), function(v) sum(as.integer(v), na.rm=TRUE))
    }
    df$jc_reads_ctrl <- parse_counts(df$IJC_SAMPLE_1) +
                        parse_counts(df$SJC_SAMPLE_1)
    df$jc_reads_trt  <- parse_counts(df$IJC_SAMPLE_2) +
                        parse_counts(df$SJC_SAMPLE_2)
    df$jc_reads_total<- df$jc_reads_ctrl + df$jc_reads_trt
  }

  # FDR column
  if ("FDR" %in% colnames(df)) {
    df$fdr <- df$FDR
  } else if ("PValue" %in% colnames(df)) {
    df$fdr <- p.adjust(df$PValue, method="BH")
  } else {
    df$fdr <- NA_real_  # --statoff: no statistics
  }

  cat(sprintf("[INF] Loaded %s: %d events\n", event_type, nrow(df)))
  df
}

event_types <- c("SE","A3SS","A5SS","MXE","RI")
rmats_list  <- lapply(event_types, function(et)
  load_rmats(opt[["rmats-dir"]], et, "JC"))
names(rmats_list) <- event_types
rmats_list <- Filter(Negate(is.null), rmats_list)

# ── 2. Filter and annotate ────────────────────────────────────────────────────
filter_events <- function(df, fdr_thr, dpsi_thr, min_reads) {
  df <- df |> filter(jc_reads_total >= min_reads)

  if (opt$statoff || all(is.na(df$fdr))) {
    # No stats available — filter by dPSI only
    df$sig <- abs(df$deltaPSI) >= dpsi_thr & !is.na(df$deltaPSI)
    df$filter_note <- "dPSI only (--statoff)"
  } else {
    df$sig <- !is.na(df$fdr) & df$fdr < fdr_thr &
              !is.na(df$deltaPSI) & abs(df$deltaPSI) >= dpsi_thr
    df$filter_note <- sprintf("FDR<%.2f & |dPSI|>%.2f", fdr_thr, dpsi_thr)
  }
  df
}

rmats_filtered <- lapply(rmats_list, filter_events,
                         fdr_thr=opt$fdr, dpsi_thr=opt$dpsi,
                         min_reads=opt[["min-reads"]])

# ── 3. Summary table ──────────────────────────────────────────────────────────
summary_df <- bind_rows(lapply(names(rmats_filtered), function(et) {
  df <- rmats_filtered[[et]]
  data.frame(
    event_type  = et,
    description = event_labels[et],
    total_events= nrow(df),
    sig_events  = sum(df$sig, na.rm=TRUE),
    sig_enhanced= sum(df$sig & df$deltaPSI > 0, na.rm=TRUE),
    sig_skipped = sum(df$sig & df$deltaPSI < 0, na.rm=TRUE),
    mean_dpsi   = round(mean(abs(df$deltaPSI[df$sig]), na.rm=TRUE), 3)
  )
}))

cat("\n[INF] Splicing event summary:\n")
print(summary_df, row.names=FALSE)
write_tsv(summary_df, file.path(opt$outdir, "splicing_summary.tsv"))

# ── 4. Combined significant events table ──────────────────────────────────────
all_sig <- bind_rows(lapply(names(rmats_filtered), function(et) {
  df <- rmats_filtered[[et]]
  df[df$sig & !is.na(df$sig), ] |>
    select(event_type, GeneID, geneSymbol, chr, strand,
           deltaPSI, psi_control, psi_treatment,
           jc_reads_total, fdr,
           any_of(c("exonStart_0base","exonEnd",
                    "upstreamES","upstreamEE",
                    "downstreamES","downstreamEE"))) |>
    arrange(desc(abs(deltaPSI)))
}))

if (nrow(all_sig) > 0) {
  write_tsv(all_sig, file.path(opt$outdir, "significant_splicing_events.tsv"))
  cat(sprintf("[INF] Total significant events: %d\n", nrow(all_sig)))
} else {
  cat("[INF] No significant events at current thresholds\n")
}

# ── 5. Per-event-type top tables ──────────────────────────────────────────────
for (et in names(rmats_filtered)) {
  df <- rmats_filtered[[et]]
  top <- df |>
    filter(!is.na(deltaPSI)) |>
    arrange(desc(abs(deltaPSI))) |>
    head(opt[["top-n"]]) |>
    select(GeneID, geneSymbol, chr, strand,
           deltaPSI, psi_control, psi_treatment,
           jc_reads_total, fdr, sig)
  write_tsv(top, file.path(opt$outdir, sprintf("top_%s_events.tsv", et)))
}

# ── 6. Plots ──────────────────────────────────────────────────────────────────

# 6a. Summary barplot — events per type
p_summary <- summary_df |>
  select(event_type, sig_enhanced, sig_skipped) |>
  pivot_longer(-event_type, names_to="direction", values_to="n") |>
  mutate(
    n = ifelse(direction=="sig_skipped", -n, n),
    direction = ifelse(direction=="sig_enhanced",
                       paste0("Increased in ", opt$treatment),
                       paste0("Decreased in ", opt$treatment))
  ) |>
  ggplot(aes(x=event_type, y=n, fill=direction)) +
  geom_col() +
  geom_hline(yintercept=0, linewidth=0.5) +
  scale_fill_manual(values=c("#E74C3C","#3498DB"), name=NULL) +
  labs(
    title="Differential alternative splicing events",
    subtitle=sprintf("%s vs %s | |dPSI|>%.2f%s",
                     opt$treatment, opt$control, opt$dpsi,
                     if (opt$statoff) " (dPSI filter only, --statoff)"
                     else sprintf(" & FDR<%.2f", opt$fdr)),
    x="Event type", y="Number of significant events"
  ) +
  theme_minimal(base_size=13)

ggsave(file.path(opt$outdir, "splicing_summary_barplot.pdf"),
       p_summary, width=8, height=5)
ggsave(file.path(opt$outdir, "splicing_summary_barplot.png"),
       p_summary, width=8, height=5, dpi=150)

# 6b. dPSI distribution per event type
all_events <- bind_rows(lapply(names(rmats_filtered), function(et) {
  rmats_filtered[[et]] |>
    filter(!is.na(deltaPSI), jc_reads_total >= opt[["min-reads"]]) |>
    select(event_type, deltaPSI, sig)
}))

if (nrow(all_events) > 0) {
  p_dpsi <- ggplot(all_events,
                   aes(x=deltaPSI, fill=event_type, color=event_type)) +
    geom_density(alpha=0.3, linewidth=0.8) +
    geom_vline(xintercept=c(-opt$dpsi, opt$dpsi),
               linetype="dashed", color="grey40") +
    scale_fill_manual(values=event_colors, name="Event type",
                      labels=event_labels) +
    scale_color_manual(values=event_colors, name="Event type",
                       labels=event_labels) +
    facet_wrap(~event_type, scales="free_y", nrow=2,
               labeller=labeller(event_type=event_labels)) +
    labs(title="deltaPSI distribution by event type",
         subtitle=sprintf("dPSI = PSI(%s) - PSI(%s) | dashed = threshold",
                          opt$treatment, opt$control),
         x="deltaPSI", y="Density") +
    theme_minimal(base_size=11) +
    theme(legend.position="none")

  ggsave(file.path(opt$outdir, "dpsi_distribution.pdf"),
         p_dpsi, width=10, height=6)
  ggsave(file.path(opt$outdir, "dpsi_distribution.png"),
         p_dpsi, width=10, height=6, dpi=150)
}

# 6c. Top SE events lollipop (most common event type)
se_df <- rmats_filtered[["SE"]]
if (!is.null(se_df) && nrow(se_df) > 0) {
  top_se <- se_df |>
    filter(!is.na(deltaPSI), jc_reads_total >= opt[["min-reads"]]) |>
    arrange(desc(abs(deltaPSI))) |>
    head(opt[["top-n"]]) |>
    mutate(
      label = ifelse(geneSymbol != "" & !is.na(geneSymbol),
                     geneSymbol, GeneID),
      direction = ifelse(deltaPSI > 0,
                         paste0("Exon included in ", opt$treatment),
                         paste0("Exon skipped in ", opt$treatment))
    )

  p_se <- ggplot(top_se,
                 aes(x=reorder(label, deltaPSI),
                     y=deltaPSI, color=direction)) +
    geom_segment(aes(xend=reorder(label, deltaPSI), yend=0),
                 linewidth=0.8, alpha=0.7) +
    geom_point(aes(size=jc_reads_total), alpha=0.9) +
    geom_hline(yintercept=0, linetype="dashed", color="grey50") +
    coord_flip() +
    scale_color_manual(values=c("#E74C3C","#3498DB"), name=NULL) +
    scale_size_continuous(range=c(2,7), name="Junction reads") +
    labs(
      title=sprintf("Top %d Skipped Exon events by |dPSI|", opt[["top-n"]]),
      subtitle=sprintf("%s vs %s", opt$treatment, opt$control),
      x=NULL, y=sprintf("dPSI (%s - %s)", opt$treatment, opt$control)
    ) +
    theme_minimal(base_size=11)

  h <- max(6, opt[["top-n"]] * 0.22 + 1.5)
  ggsave(file.path(opt$outdir, "top_SE_lollipop.pdf"),
         p_se, width=10, height=h)
  ggsave(file.path(opt$outdir, "top_SE_lollipop.png"),
         p_se, width=10, height=h, dpi=150)
}

# 6d. Volcano-style: dPSI vs -log10(FDR) for SE
if (!is.null(se_df) && !all(is.na(se_df$fdr))) {
  p_vol_se <- se_df |>
    filter(!is.na(deltaPSI), !is.na(fdr),
           jc_reads_total >= opt[["min-reads"]]) |>
    mutate(
      col = case_when(
        sig & deltaPSI > 0  ~ "Included",
        sig & deltaPSI < 0  ~ "Skipped",
        TRUE                 ~ "NS"
      ),
      lab = ifelse(sig & abs(deltaPSI) > opt$dpsi * 1.5,
                   geneSymbol, NA)
    ) |>
    ggplot(aes(x=deltaPSI, y=-log10(fdr + 1e-10), color=col)) +
    geom_point(alpha=0.6, size=1.5) +
    ggrepel::geom_text_repel(aes(label=lab), size=2.5,
                              max.overlaps=15, show.legend=FALSE) +
    scale_color_manual(values=c(Included="#E74C3C",Skipped="#3498DB",NS="#BDC3C7"),
                       name=NULL) +
    geom_vline(xintercept=c(-opt$dpsi, opt$dpsi),
               linetype="dashed", color="grey50") +
    geom_hline(yintercept=-log10(opt$fdr),
               linetype="dashed", color="grey50") +
    labs(title="Skipped Exon events - volcano plot",
         x=sprintf("dPSI (%s - %s)", opt$treatment, opt$control),
         y="-log10(FDR)") +
    theme_minimal(base_size=12)

  ggsave(file.path(opt$outdir, "SE_volcano.pdf"), p_vol_se, width=8, height=6)
  ggsave(file.path(opt$outdir, "SE_volcano.png"), p_vol_se, width=8, height=6, dpi=150)
}

# ── 7. Gene overlap with DE genes ─────────────────────────────────────────────
# Check if any significantly spliced genes are also DE
de_path <- file.path(dirname(opt[["rmats-dir"]]),
                     "deseq2/genes_tumor_vs_normal_DE.tsv")
if (file.exists(de_path)) {
  de_df   <- read.table(de_path, header=TRUE, sep="\t",
                         stringsAsFactors=FALSE)
  gid_col <- colnames(de_df)[ncol(de_df)]
  de_sig_ids <- de_df[[gid_col]][!is.na(de_df$padj) & de_df$padj < 0.05]

  if (nrow(all_sig) > 0) {
    overlap <- all_sig |>
      filter(GeneID %in% de_sig_ids) |>
      select(event_type, GeneID, geneSymbol, deltaPSI, fdr)
    cat(sprintf("[INF] Spliced genes also DE: %d\n", nrow(overlap)))
    if (nrow(overlap) > 0)
      write_tsv(overlap, file.path(opt$outdir, "spliced_and_de_overlap.tsv"))
  }
}

cat("\n[OK] M17 complete. Output:", opt$outdir, "\n")
cat(sprintf("[INF] Significant events: %d total across %d event types\n",
            sum(summary_df$sig_events), nrow(summary_df)))
