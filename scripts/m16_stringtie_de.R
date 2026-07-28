#!/usr/bin/env Rscript
# mosaic/scripts/m16_stringtie_de.R
# M16 - StringTie de novo transcript assembly and DE analysis
#
# Input:  t_data.ctab per sample (from stringtie -e -b)
#         novel_transcripts.tmap (from gffcompare, class codes u/i/j/x/s/e/p)
#         gffcmp.merged.gtf.tmap (all class codes)
#         featurecounts_transcripts.txt (raw integer counts from featureCounts)
# Output: novel_transcript_fpkm.tsv   (FPKM from StringTie, for display only)
#         novel_transcript_DE.tsv     (DESeq2 on featureCounts raw counts)
#         stringtie_summary.tsv
#         plots: class_code_barplot, novel_de_volcano

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

option_list <- list(
  make_option("--st-dir",         type="character", help="StringTie output dir (m16_stringtie)"),
  make_option("--tmap",           type="character", help="gffcmp.merged.gtf.tmap"),
  make_option("--novel-tmap",     type="character", help="novel_transcripts.tmap"),
  make_option("--feature-counts", type="character", default=NULL,
              help="featureCounts output (raw integer counts per transcript)"),
  make_option("--coldata",        type="character", help="coldata.tsv (sample_id, condition)"),
  make_option("--outdir",         type="character", default="mosaic_output/m16_stringtie"),
  make_option("--treatment",      type="character", default="tumor"),
  make_option("--control",        type="character", default="normal"),
  make_option("--padj",           type="double",    default=0.05),
  make_option("--lfc",            type="double",    default=1.0),
  make_option("--min-cov",        type="double",    default=1.0,
              help="Min mean coverage to keep transcript"),
  make_option("--min-tpm",       type="double",    default=1.0,
              help="Min mean TPM across all samples (from t_data.ctab; applied before DESeq2)"),
  make_option("--min-samples",   type="integer",   default=2L,
              help="Min samples per group with count>0 (featureCounts; applied before DESeq2)")
)

opt <- parse_args(OptionParser(option_list=option_list))
dir.create(opt$outdir, recursive=TRUE, showWarnings=FALSE)

cat("[INF] M16 - StringTie DE analysis\n")

safe_read <- function(path, ...) {
  if (!is.null(path) && file.exists(path))
    tryCatch(read.table(path, header=TRUE, sep="\t",
                        check.names=FALSE, stringsAsFactors=FALSE, ...),
             error=function(e) { cat("[WARN]", conditionMessage(e), "\n"); NULL })
  else { cat("[WARN] Not found:", path %||% "NULL", "\n"); NULL }
}
`%||%` <- function(a,b) if(!is.null(a)) a else b
write_tsv <- function(df, path) {
  write.table(df, path, sep="\t", quote=FALSE, row.names=FALSE)
  cat("[OK] Saved:", basename(path), sprintf("(%d rows)\n", nrow(df)))
}

# DESeq2 requires raw integer counts. This function rejects FPKM/TPM inputs
# before they silently corrupt dispersion estimation and p-values.
validate_count_matrix <- function(mat, label = "count matrix") {
  if (any(mat < 0, na.rm = TRUE))
    stop(sprintf(
      "[ERROR] %s contains negative values. DESeq2 requires non-negative integer counts.",
      label))
  non_int <- mat[!is.na(mat) & mat != round(mat)]
  if (length(non_int) > 0)
    stop(sprintf(paste0(
      "[ERROR] %s contains non-integer values (e.g. %.6g). ",
      "DESeq2 requires raw counts, not FPKM or TPM. ",
      "Ensure featureCounts output is used, not StringTie FPKM."),
      label, non_int[1]))
  invisible(TRUE)
}

# Class code descriptions
class_desc <- c(
  "=" = "Exact match to reference",
  "j" = "Novel junction combination",
  "k" = "Containment of reference",
  "c" = "Contained in reference",
  "m" = "Retained intron (all introns)",
  "n" = "Retained intron (partial)",
  "o" = "Overlapping (different strand)",
  "x" = "Antisense to reference",
  "i" = "Novel intronic transcript",
  "u" = "Novel intergenic transcript",
  "e" = "Single exon overlapping",
  "s" = "Intron-matching on opposite strand",
  "p" = "Possible polymerase run-on"
)

# ── 1. Load tmap and summarise class codes ─────────────────────────────────────
tmap <- safe_read(opt$tmap)
# novel_transcripts.tmap has no header (extracted with awk)
novel_tmap_raw <- if (!is.null(opt[["novel-tmap"]]) && file.exists(opt[["novel-tmap"]]))
  tryCatch(
    read.table(opt[["novel-tmap"]], header=FALSE, sep="	",
               stringsAsFactors=FALSE, fill=TRUE,
               col.names=c("ref_gene_id","ref_id","class_code","qry_gene_id",
                           "qry_id","num_exons","FPKM","TPM","cov","len",
                           "major_iso_id","ref_match_len")),
    error=function(e) NULL
  ) else NULL
novel_tmap <- novel_tmap_raw

if (!is.null(tmap)) {
  cat("[INF] Class code distribution:\n")
  cc_summary <- tmap |>
    count(class_code) |>
    mutate(
      description = class_desc[class_code] %||% "Unknown",
      novel = class_code %in% c("u","i","j","x","n","o")
    ) |>
    arrange(desc(n))
  print(cc_summary)
  write_tsv(cc_summary, file.path(opt$outdir, "class_code_summary.tsv"))

  # Plot
  p_cc <- cc_summary |>
    mutate(label=paste0(class_code, "\n", sub(" .*","",description)),
           label=reorder(label, n)) |>
    ggplot(aes(x=label, y=n, fill=novel)) +
    geom_col() +
    scale_fill_manual(values=c("TRUE"="#E74C3C","FALSE"="#85C1E9"),
                      labels=c("TRUE"="Novel","FALSE"="Annotated"),
                      name=NULL) +
    scale_y_continuous(labels=scales::comma) +
    labs(title="StringTie transcript classification (gffcompare)",
         subtitle="Red = novel/partially novel; blue = matching reference",
         x="Class code", y="Number of transcripts") +
    theme_minimal(base_size=12)
  ggsave(file.path(opt$outdir, "class_code_barplot.pdf"), p_cc, width=9, height=5)
  ggsave(file.path(opt$outdir, "class_code_barplot.png"), p_cc, width=9, height=5, dpi=150)
}

# ── 1b. Load featureCounts raw integer counts ─────────────────────────────────
# featurecounts_transcripts.txt is generated by pipeline.py after gffcompare.
# Column names are full BAM paths; we extract the sample ID from the filename.
fc_counts <- NULL
if (!is.null(opt[["feature-counts"]]) && file.exists(opt[["feature-counts"]])) {
  fc_raw  <- read.table(opt[["feature-counts"]], header=TRUE, sep="\t",
                        comment.char="#", check.names=FALSE, stringsAsFactors=FALSE)
  bam_cols <- grep("\\.bam$", colnames(fc_raw), value=TRUE)
  # Extract sample ID: strip directory prefix and _Aligned.sortedByCoord.out.bam suffix
  sid_names <- sub("^.*/([^/]+)_Aligned\\.sortedByCoord\\.out\\.bam$", "\\1", bam_cols)
  # Fall back to basename without .bam if pattern doesn't match
  unmatched <- sid_names == bam_cols
  if (any(unmatched))
    sid_names[unmatched] <- sub("\\.bam$", "", basename(bam_cols[unmatched]))
  fc_mat <- as.matrix(fc_raw[, bam_cols, drop=FALSE])
  rownames(fc_mat) <- fc_raw$Geneid
  colnames(fc_mat) <- sid_names
  storage.mode(fc_mat) <- "integer"
  cat(sprintf("[INF] featureCounts matrix: %d transcripts × %d samples\n",
              nrow(fc_mat), ncol(fc_mat)))
  fc_counts <- fc_mat
} else {
  cat("[WARN] --feature-counts not provided; novel transcript DE will be skipped.\n")
  cat("[WARN] Re-run pipeline to generate featurecounts_transcripts.txt.\n")
}

# ── 2. Build count matrix from t_data.ctab ────────────────────────────────────
cat("[INF] Building transcript count matrix...\n")
coldata <- safe_read(opt$coldata)
if (is.null(coldata)) {
  cat("[ERROR] --coldata required\n"); quit(status=1)
}

# Find sample directories
sample_dirs <- list.dirs(opt[["st-dir"]], recursive=FALSE, full.names=TRUE)
quant_dirs  <- sample_dirs[grepl("_quant$", sample_dirs)]

if (length(quant_dirs) == 0) {
  cat("[ERROR] No *_quant directories found in", opt[["st-dir"]], "\n")
  quit(status=1)
}

cat(sprintf("[INF] Found %d quantification directories\n", length(quant_dirs)))

# Load t_data.ctab for each sample
ctab_list <- lapply(quant_dirs, function(d) {
  sid  <- sub("_quant$", "", basename(d))
  path <- file.path(d, "t_data.ctab")
  if (!file.exists(path)) { cat("[WARN] Missing:", path, "\n"); return(NULL) }
  df <- read.table(path, header=TRUE, sep="\t", stringsAsFactors=FALSE)
  df$sample_id <- sid
  df
})
ctab_list <- Filter(Negate(is.null), ctab_list)

if (length(ctab_list) == 0) {
  cat("[ERROR] No t_data.ctab files found\n"); quit(status=1)
}

all_ctab <- bind_rows(ctab_list)

# Wide matrix: rows = transcripts, cols = samples (FPKM)
fpkm_wide <- all_ctab |>
  select(t_name, sample_id, FPKM) |>
  pivot_wider(names_from=sample_id, values_from=FPKM, values_fill=0)

# Coverage matrix for filtering
cov_wide <- all_ctab |>
  select(t_name, sample_id, cov) |>
  pivot_wider(names_from=sample_id, values_from=cov, values_fill=0)

# Transcript metadata
tx_meta <- all_ctab |>
  select(t_id, t_name, chr, strand, start, end, num_exons, length,
         gene_id, gene_name) |>
  distinct()

cat(sprintf("[INF] Total transcripts: %d | Samples: %d\n",
            nrow(fpkm_wide), length(quant_dirs)))

# ── 3. Filter to novel transcripts ────────────────────────────────────────────
# novel_tmap may be NULL or 0-row (empty after structural filters).
# Guard against type errors in joins when data frame is empty.
novel_ids <- if (!is.null(novel_tmap) && nrow(novel_tmap) > 0) {
  unique(as.character(novel_tmap$qry_id))
} else {
  character(0)
}

cat(sprintf("[INF] Novel transcript IDs from tmap (structural filter): %d\n",
            length(novel_ids)))

# Filter matrix to novel
novel_fpkm <- fpkm_wide |> filter(t_name %in% novel_ids)
novel_cov   <- cov_wide  |> filter(t_name %in% novel_ids)

cat(sprintf("[INF] Novel transcripts in matrix: %d\n", nrow(novel_fpkm)))

# min_tpm filter: mean TPM across all samples from t_data.ctab (pass-4 re-quant).
# stringtie --merge does not write TPM to merged.gtf, so TMAP TPM is always 0;
# per-sample TPM from -e re-quantification is the correct source.
min_tpm <- opt[["min-tpm"]]
if (nrow(novel_fpkm) > 0 && min_tpm > 0) {
  tpm_mat  <- as.matrix(novel_fpkm[, -1, drop=FALSE])
  mean_tpm <- rowMeans(tpm_mat, na.rm=TRUE)
  keep     <- mean_tpm >= min_tpm
  cat(sprintf("[INF] min_tpm >= %.2f: %d / %d transcripts pass\n",
              min_tpm, sum(keep), nrow(novel_fpkm)))
  novel_fpkm <- novel_fpkm[keep, ]
  novel_ids  <- novel_fpkm$t_name
  novel_cov  <- novel_cov |> filter(t_name %in% novel_ids)
}

# Filter by minimum coverage
if (nrow(novel_cov) > 0) {
  cov_mat    <- as.matrix(novel_cov[,-1])
  mean_cov   <- rowMeans(cov_mat, na.rm=TRUE)
  keep       <- mean_cov >= opt[["min-cov"]]
  novel_fpkm <- novel_fpkm[keep, ]
  cat(sprintf("[INF] After min coverage filter (>= %.1f): %d transcripts\n",
              opt[["min-cov"]], nrow(novel_fpkm)))
}

# Add class code annotation — guard against empty novel_tmap
if (!is.null(novel_tmap) && nrow(novel_tmap) > 0) {
  anno <- novel_tmap |>
    select(t_name=qry_id, class_code, ref_gene_id, num_exons, len) |>
    mutate(t_name=as.character(t_name),
           class_desc=class_desc[class_code])
  novel_fpkm <- novel_fpkm |> left_join(anno, by="t_name")
}

write_tsv(novel_fpkm, file.path(opt$outdir, "novel_transcript_fpkm.tsv"))

# ── 4. DE analysis on novel transcripts (featureCounts raw counts → DESeq2) ───
# DESeq2 dispersion estimation fits a mean-dispersion curve across all features.
# Running on only the ~67 pre-filtered novel transcripts collapses the empirical
# Bayes shrinkage → underestimated dispersions → anti-conservative p-values.
# FIX: build DESeqDataSet on the FULL featureCounts matrix (~507K transcripts),
# estimate dispersions globally, then extract results only for novel transcript IDs.
cat("[INF] DE analysis on novel transcripts...\n")

if (!requireNamespace("DESeq2", quietly=TRUE)) {
  cat("[WARN] DESeq2 not available — skipping DE\n")
} else if (length(novel_ids) < 3) {
  cat("[WARN] Too few novel transcripts for DE\n")
} else if (is.null(fc_counts)) {
  cat("[WARN] featureCounts matrix not available — skipping DESeq2.\n")
  cat("[WARN] FPKM proxy counts are statistically invalid for DESeq2 (mean-variance\n")
  cat("[WARN] relationship is destroyed by length/depth normalisation). Re-run the\n")
  cat("[WARN] pipeline to generate featurecounts_transcripts.txt.\n")
} else {
  library(DESeq2)

  # Subset to only treatment and control samples; coldata may contain
  # additional conditions that would produce NA in the DESeq2 design.
  sample_cols <- coldata$sample_id[
    coldata$sample_id %in% colnames(fc_counts) &
    coldata$condition %in% c(opt$treatment, opt$control)
  ]

  # Full matrix: all transcripts (known + novel) × comparison samples
  full_mat <- fc_counts[, sample_cols, drop=FALSE]
  full_mat <- full_mat[rowSums(full_mat) > 0, ]  # remove all-zero rows
  cat(sprintf("[INF] Full matrix for dispersion estimation: %d transcripts × %d samples\n",
              nrow(full_mat), ncol(full_mat)))

  tryCatch(
    validate_count_matrix(full_mat, "M16 featureCounts (full matrix)"),
    error = function(e) { cat(conditionMessage(e), "\n"); stop(e) }
  )

  if (nrow(full_mat) < 100 || ncol(full_mat) < 4) {
    cat("[WARN] Insufficient transcripts/samples for DESeq2\n")
  } else {
    cd_full <- coldata[match(sample_cols, coldata$sample_id), , drop=FALSE]
    cd_full$condition <- factor(cd_full$condition,
                                levels=c(opt$control, opt$treatment))

    dds_full <- DESeqDataSetFromMatrix(full_mat, cd_full, ~condition)
    dds_full <- DESeq(dds_full, quiet=TRUE)

    # Extract results only for novel transcript IDs (those present in matrix)
    keep_ids <- novel_ids[novel_ids %in% rownames(dds_full)]
    cat(sprintf("[INF] Novel transcripts in full matrix: %d\n", length(keep_ids)))

    if (length(keep_ids) >= 3) {
      # Unshrunken results: padj and pvalue used for significance.
      # padj is NOT changed by lfcShrink — it comes from the Wald test here.
      res <- results(dds_full,
                     contrast=c("condition", opt$treatment, opt$control),
                     alpha=opt$padj)

      # LFC shrinkage (ashr) — same approach as M08.
      # Regularises log2FC estimates for low-count transcripts.
      # Adds log2FoldChange_shrunk, lfcSE_shrunk, svalue columns.
      # Note: extreme padj values (e.g. near-zero within-group variance at n=3)
      # reflect genuine signal but are precision artefacts; padj is unchanged.
      shrunk_df <- NULL
      if (requireNamespace("ashr", quietly=TRUE)) {
        res_shrunk <- lfcShrink(dds_full,
                                contrast=c("condition", opt$treatment, opt$control),
                                type="ashr", quiet=TRUE)
        shrunk_raw <- as.data.frame(res_shrunk) |>
          tibble::rownames_to_column("t_name")
        # svalue is stored in metadata by some DESeq2 versions; extract defensively
        sval_vec <- tryCatch(res_shrunk$svalue, error=function(e) NULL)
        if (!is.null(sval_vec)) {
          shrunk_raw$svalue <- sval_vec
        }
        shrunk_cols <- c("t_name", "log2FoldChange", "lfcSE")
        if ("svalue" %in% colnames(shrunk_raw)) shrunk_cols <- c(shrunk_cols, "svalue")
        shrunk_df <- shrunk_raw |>
          select(all_of(shrunk_cols)) |>
          rename(log2FoldChange_shrunk=log2FoldChange, lfcSE_shrunk=lfcSE)
        cat(sprintf("[INF] lfcShrink (ashr) applied — columns: %s\n",
                    paste(colnames(shrunk_df)[-1], collapse=", ")))
      } else {
        cat("[WARN] ashr not available — reporting unshrunken LFC only\n")
      }

      res_df <- as.data.frame(res[keep_ids, ]) |>
        tibble::rownames_to_column("t_name") |>
        { df <- .
          if (!is.null(shrunk_df)) {
            left_join(df,
                      shrunk_df |> filter(t_name %in% keep_ids),
                      by="t_name")
          } else df
        }() |>
        left_join(novel_fpkm |>
                    select(t_name, any_of(c("class_code","class_desc",
                                            "ref_gene_id","num_exons","len"))),
                  by="t_name") |>
        left_join(tx_meta |> select(t_name, chr, strand, start, end,
                                    gene_id, gene_name),
                  by="t_name") |>
        arrange(padj, desc(abs(log2FoldChange)))

      write_tsv(res_df, file.path(opt$outdir, "novel_transcript_DE.tsv"))

      n_sig <- sum(!is.na(res_df$padj) & res_df$padj < opt$padj &
                     abs(res_df$log2FoldChange) > opt$lfc)
      cat(sprintf("[INF] Significant novel DE transcripts: %d / %d\n",
                  n_sig, nrow(res_df)))

      # Volcano: use shrunken LFC on x-axis if available (less extreme for noisy transcripts)
      lfc_col <- if ("log2FoldChange_shrunk" %in% names(res_df)) "log2FoldChange_shrunk" else "log2FoldChange"

      p_vol <- res_df |>
        filter(!is.na(padj), !is.na(.data[[lfc_col]])) |>
        mutate(
          sig = padj < opt$padj & abs(.data[[lfc_col]]) > opt$lfc,
          lab = ifelse(sig, t_name, NA),
          col = case_when(
            sig & .data[[lfc_col]] > 0 ~ "UP",
            sig & .data[[lfc_col]] < 0 ~ "DOWN",
            TRUE ~ "NS"
          )
        ) |>
        ggplot(aes(x=.data[[lfc_col]], y=-log10(padj), color=col)) +
        geom_point(alpha=0.7, size=1.8) +
        ggrepel::geom_text_repel(aes(label=lab), size=2.5,
                                  max.overlaps=15, show.legend=FALSE) +
        scale_color_manual(values=c(UP="#E74C3C",DOWN="#3498DB",NS="#BDC3C7"),
                           name=NULL) +
        geom_vline(xintercept=c(-opt$lfc, opt$lfc),
                   linetype="dashed", color="grey50") +
        geom_hline(yintercept=-log10(opt$padj),
                   linetype="dashed", color="grey50") +
        labs(title="Novel transcript DE - volcano plot",
             subtitle=sprintf("%s vs %s | padj<%.2f | |LFC|>%.1f | %d sig",
                              opt$treatment, opt$control,
                              opt$padj, opt$lfc, n_sig),
             x=sprintf("log2FC shrunken (%s / %s)", opt$treatment, opt$control),
             y="-log10(padj)") +
        theme_minimal(base_size=12)

      ggsave(file.path(opt$outdir, "novel_de_volcano.pdf"), p_vol, width=8, height=6)
      ggsave(file.path(opt$outdir, "novel_de_volcano.png"), p_vol, width=8, height=6, dpi=150)
    } else {
      cat("[WARN] Too few novel IDs in featureCounts matrix — skipping DE\n")
    }
  }
}

# ── 5. Summary table ──────────────────────────────────────────────────────────
summary_df <- data.frame(
  metric = c(
    "Total merged transcripts",
    "Exact match (=)",
    "Novel junction combo (j)",
    "Novel intronic (i)",
    "Novel intergenic (u)",
    "Antisense (x)",
    "Other novel",
    "Novel candidates (filtered)",
    "Novel after coverage filter"
  ),
  value = c(
    if (!is.null(tmap)) nrow(tmap) else NA,
    if (!is.null(tmap)) sum(tmap$class_code=="=") else NA,
    if (!is.null(tmap)) sum(tmap$class_code=="j") else NA,
    if (!is.null(tmap)) sum(tmap$class_code=="i") else NA,
    if (!is.null(tmap)) sum(tmap$class_code=="u") else NA,
    if (!is.null(tmap)) sum(tmap$class_code=="x") else NA,
    if (!is.null(tmap)) sum(tmap$class_code %in% c("n","o","e","s","p")) else NA,
    length(novel_ids),
    nrow(novel_fpkm)
  )
)
write_tsv(summary_df, file.path(opt$outdir, "stringtie_summary.tsv"))

cat("\n[INF] M16 Summary:\n")
print(summary_df, row.names=FALSE)
cat("[OK] M16 complete. Output:", opt$outdir, "\n")
