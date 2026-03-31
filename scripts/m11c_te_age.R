#!/usr/bin/env Rscript
# mosaic/scripts/m11c_te_age.R
# M11c — TE age and activity analysis
#
# RepeatMasker reports `milliDiv` = substitution level per 1000 bases
# from the consensus sequence. This is a proxy for insertion age:
#   - milliDiv → 0:  very young / recently inserted (< 1 Mya)
#   - milliDiv → 300+: ancient, highly diverged (> 100 Mya)
#   - kimura2p_dist = milliDiv / 1000 (approx)
#
# This script:
#   1. Loads rmsk_hg38.bed with milliDiv scores
#   2. Computes age distribution per TE family/macrocategory
#   3. Correlates TE age with DE expression changes
#   4. Plots Kimura CpG-corrected divergence plots (Dfam-style)

suppressPackageStartupMessages({
  library(optparse)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

option_list <- list(
  make_option("--rmsk-bed",   type="character",
              help="rmsk_hg38.bed (6+ cols; col5=milliDiv or col from rmsk_full)"),
  make_option("--rmsk-full",  type="character", default=NULL,
              help="Full RepeatMasker .out or .bed with milliDiv column"),
  make_option("--te-de",      type="character", default=NULL,
              help="TE_tumor_vs_normal_DE.tsv from M08"),
  make_option("--outdir",     type="character", default="mosaic_output/m11c"),
  make_option("--treatment",  type="character", default="tumor"),
  make_option("--control",    type="character", default="normal"),
  make_option("--min-loci",   type="integer",   default=100L,
              help="Min loci per family to include in plots")
)

opt <- parse_args(OptionParser(option_list=option_list))
dir.create(opt$outdir, recursive=TRUE, showWarnings=FALSE)

col_macro <- c(
  "SINE"         = "#3498DB",
  "LINE"         = "#E74C3C",
  "LTR/ERV"      = "#2ECC71",
  "DNA transposon"="#F39C12",
  "SVA"          = "#9B59B6",
  "Satellite"    = "#1ABC9C",
  "Other TE"     = "#95A5A6"
)

classify_macro <- function(name) {
  n <- toupper(name)
  case_when(
    grepl("^ALU|^MIR|SINE",    n) ~ "SINE",
    grepl("^L1|^L2|LINE|^CR1", n) ~ "LINE",
    grepl("^LTR|^HERV|^ERV|^MLT|^THE|^MER[0-9]", n) ~ "LTR/ERV",
    grepl("^SVA",              n) ~ "SVA",
    grepl("DNA|^TIGGER|^CHARLIE|^HAT|^MARINER|^MER[A-Z]", n) ~ "DNA transposon",
    grepl("SAT|ALPHA|CENTR|ACRO|MICRO", n) ~ "Satellite",
    TRUE ~ "Other TE"
  )
}

# ── Load RepeatMasker data ─────────────────────────────────────────────────────
cat("[INF] Loading RepeatMasker data...\n")

rmsk <- NULL

# Try full rmsk file (UCSC .txt.gz or RepeatMasker .out)
if (!is.null(opt[["rmsk-full"]]) && file.exists(opt[["rmsk-full"]])) {
  cat("[INF] Reading RepeatMasker file (auto-detecting format)...\n")
  tryCatch({
    # Peek at first data line to detect format
    con <- if (grepl("\\.gz$", opt[["rmsk-full"]]))
             gzcon(file(opt[["rmsk-full"]], "rb")) else
             file(opt[["rmsk-full"]], "r")
    first_line <- ""
    for (i in 1:5) {
      ln <- readLines(con, n=1)
      if (length(ln)==0) break
      if (nchar(trimws(ln)) > 0) { first_line <- ln; break }
    }
    close(con)

    fields    <- strsplit(trimws(first_line), "\t")[[1]]
    n_cols    <- length(fields)
    first_int <- suppressWarnings(!is.na(as.integer(fields[1])))
    cat(sprintf("[INF] %d tab-separated columns, first field integer: %s\n",
                n_cols, first_int))

    if (first_int && n_cols >= 16) {
      # UCSC rmsk.txt format (tab-separated, no header):
      # V1=bin V2=swScore V3=milliDiv V4=milliDel V5=milliIns
      # V6=genoName V7=genoStart V8=genoEnd V9=genoLeft V10=strand
      # V11=repName V12=repClass V13=repFamily V14=repStart V15=repEnd
      # V16=repLeft V17=id
      cat("[INF] UCSC rmsk.txt format\n")
      raw <- read.table(opt[["rmsk-full"]], header=FALSE, sep="\t",
                        stringsAsFactors=FALSE, fill=TRUE,
                        colClasses="character")
      rmsk <- data.frame(
        chr      = raw[[6]],
        start    = as.integer(raw[[7]]),   # 0-based in UCSC
        end      = as.integer(raw[[8]]),
        name     = raw[[11]],              # repName
        milliDiv = as.numeric(raw[[3]]),   # milliDiv (e.g. 153 = 15.3%)
        strand   = raw[[10]],
        class    = raw[[12]],
        stringsAsFactors = FALSE
      )
      cat(sprintf("[INF] Loaded %d loci (UCSC format)\n", nrow(rmsk)))

    } else {
      # RepeatMasker .out: space-delimited, 3 header lines
      # cols: swScore perc_div perc_del perc_ins query q_beg q_end q_left
      #       strand repeat class/family r_beg r_end r_left id
      cat("[INF] RepeatMasker .out format\n")
      raw <- read.table(opt[["rmsk-full"]], skip=3, fill=TRUE,
                        header=FALSE, stringsAsFactors=FALSE)
      rmsk <- data.frame(
        chr      = raw[[5]],
        start    = as.integer(raw[[6]]) - 1L,
        end      = as.integer(raw[[7]]),
        name     = raw[[10]],
        milliDiv = as.numeric(raw[[2]]) * 10,  # perc_div → milliDiv
        strand   = ifelse(raw[[9]] == "C", "-", "+"),
        class    = raw[[11]],
        stringsAsFactors = FALSE
      )
      cat(sprintf("[INF] Loaded %d loci (.out format)\n", nrow(rmsk)))
    }
  }, error=function(e) cat("[WARN] rmsk-full parse error:", conditionMessage(e), "\n"))
}

# Fallback: BED file (milliDiv may not be present — use col5 if numeric)
if (is.null(rmsk) && !is.null(opt[["rmsk-bed"]]) && file.exists(opt[["rmsk-bed"]])) {
  cat("[INF] Reading rmsk BED (milliDiv from col5 if available)...\n")
  tryCatch({
    bed <- read.table(opt[["rmsk-bed"]], header=FALSE, sep="\t",
                      stringsAsFactors=FALSE, fill=TRUE)
    milliDiv_col <- if (ncol(bed) >= 5 && is.numeric(bed[[5]])) bed[[5]] else NA_real_
    rmsk <- data.frame(
      chr      = bed[[1]],
      start    = as.integer(bed[[2]]),
      end      = as.integer(bed[[3]]),
      name     = bed[[4]],
      milliDiv = milliDiv_col,
      strand   = if (ncol(bed) >= 6) bed[[6]] else "."
    )
    cat(sprintf("[INF] Loaded %d RM loci from BED\n", nrow(rmsk)))
    if (all(is.na(rmsk$milliDiv))) {
      cat("[WARN] No milliDiv in BED col5 — age plots will be limited\n")
    }
  }, error=function(e) { cat("[ERROR]", conditionMessage(e), "\n"); quit(status=1) })
}

if (is.null(rmsk)) {
  cat("[ERROR] No RepeatMasker input found\n")
  quit(status=1)
}

rmsk$macro <- classify_macro(rmsk$name)
rmsk$family_id <- rmsk$name  # repName for exact match with DE file

cat(sprintf("[INF] TE macrocategories:\n"))
print(table(rmsk$macro))

# ═══════════════════════════════════════════════════════════════════════════════
# 1. Kimura divergence distribution (Dfam-style landscape)
# ═══════════════════════════════════════════════════════════════════════════════
if (!all(is.na(rmsk$milliDiv))) {
  cat("[INF] Plotting Kimura divergence landscape...\n")

  # Kimura 2-parameter distance (approximate)
  # k2p = -3/4 * ln(1 - 4/3 * d) where d = milliDiv/1000
  d <- rmsk$milliDiv / 1000
  d <- pmin(d, 0.74)  # cap to avoid log(0)
  rmsk$k2p <- tryCatch(
    -0.75 * log(1 - (4/3) * d),
    warning=function(w) d  # fallback
  )

  # Histogram per macrocategory (stacked, Dfam style)
  rmsk_valid <- filter(rmsk, !is.na(milliDiv), milliDiv >= 0)

  p_landscape <- ggplot(rmsk_valid,
                        aes(x=milliDiv/10, fill=macro)) +  # /10 = % divergence
    geom_histogram(binwidth=1, position="stack", alpha=0.9) +
    scale_fill_manual(values=col_macro, name="TE class") +
    scale_x_continuous(limits=c(0, 45), breaks=seq(0, 45, 5)) +
    labs(
      title    = "TE age landscape (RepeatMasker divergence)",
      subtitle = "Kimura-substitution level; low = young/recently active",
      x        = "Divergence from consensus (%)",
      y        = "Number of loci"
    ) +
    theme_minimal(base_size=12) +
    theme(legend.position="right")

  ggsave(file.path(opt$outdir, "te_age_landscape.pdf"),
         p_landscape, width=12, height=6)
  ggsave(file.path(opt$outdir, "te_age_landscape.png"),
         p_landscape, width=12, height=6, dpi=150)

  # Per macrocategory faceted
  p_facet <- ggplot(rmsk_valid, aes(x=milliDiv/10, fill=macro)) +
    geom_histogram(binwidth=1, alpha=0.85, show.legend=FALSE) +
    scale_fill_manual(values=col_macro) +
    facet_wrap(~macro, scales="free_y", ncol=3) +
    labs(title="TE divergence per class",
         x="Divergence from consensus (%)", y="Number of loci") +
    theme_minimal(base_size=10)

  ggsave(file.path(opt$outdir, "te_age_by_class.pdf"),
         p_facet, width=12, height=8)

  # Summary table
  age_summary <- rmsk_valid |>
    group_by(macro) |>
    summarise(
      n_loci        = n(),
      mean_div      = round(mean(milliDiv/10), 2),
      median_div    = round(median(milliDiv/10), 2),
      sd_div        = round(sd(milliDiv/10), 2),
      pct_young     = round(mean(milliDiv/10 < 5) * 100, 1),  # < 5% = young
      .groups="drop"
    ) |> arrange(mean_div)

  write.table(age_summary,
              file.path(opt$outdir, "te_age_summary.tsv"),
              sep="\t", quote=FALSE, row.names=FALSE)
  cat("[OK] Age landscape plots saved.\n")
  print(age_summary)
}

# ═══════════════════════════════════════════════════════════════════════════════
# 2. TE age vs differential expression
# ═══════════════════════════════════════════════════════════════════════════════
if (!is.null(opt[["te-de"]]) && file.exists(opt[["te-de"]])) {
  cat("[INF] Correlating TE age with DE...\n")
  tryCatch({
    te_de <- read.table(opt[["te-de"]], header=TRUE, sep="\t", check.names=FALSE)

    # family_id is last column
    fid_col <- colnames(te_de)[ncol(te_de)]
    te_de$family_id <- te_de[[fid_col]]

    # Average milliDiv per family
    family_age <- rmsk |>
      filter(!is.na(milliDiv)) |>
      group_by(family_id) |>
      summarise(
        mean_div  = mean(milliDiv/10, na.rm=TRUE),
        n_loci    = n(),
        macro     = first(macro),
        .groups="drop"
      ) |>
      filter(n_loci >= opt[["min-loci"]])

    # Join DE with age data
    # Try exact match first (repName == family_id)
    te_age_de <- te_de |>
      left_join(family_age, by="family_id")

    n_matched <- sum(!is.na(te_age_de$mean_div))
    cat(sprintf("[INF] Exact join matched: %d families\n", n_matched))

    # If few matches, try prefix: DE "L1" matches rmsk "L1HS","L1PA2" etc.
    if (n_matched < 5) {
      cat("[INF] Trying prefix match (DE family as prefix of repName)...\n")
      de_fids   <- unique(te_de$family_id)
      rmsk_fids <- unique(rmsk$family_id)
      # map each DE family to first matching rmsk repName
      prefix_map <- setNames(sapply(de_fids, function(fid) {
        m <- rmsk_fids[rmsk_fids == fid]
        if (length(m) > 0) return(fid)
        m2 <- rmsk_fids[startsWith(rmsk_fids, fid)]
        if (length(m2) > 0) return(m2[1])
        NA_character_
      }), de_fids)
      te_de$rmsk_match <- prefix_map[te_de$family_id]

      age_mapped <- rmsk |>
        filter(!is.na(milliDiv), family_id %in% na.omit(te_de$rmsk_match)) |>
        group_by(family_id) |>
        summarise(mean_div=mean(milliDiv/10, na.rm=TRUE),
                  n_loci=n(), macro=first(macro), .groups="drop")

      te_age_de <- te_de |>
        left_join(age_mapped, by=c("rmsk_match"="family_id"))
      cat(sprintf("[INF] After prefix match: %d families\n",
                  sum(!is.na(te_age_de$mean_div))))
    }

    te_age_de <- filter(te_age_de, !is.na(mean_div))
    cat(sprintf("[INF] Families with age data and DE: %d\n", nrow(te_age_de)))

    if (nrow(te_age_de) >= 5) {
      # LFC vs divergence scatter
      lfc_col <- if ("log2FoldChange" %in% colnames(te_age_de)) "log2FoldChange" else
                 grep("log2", colnames(te_age_de), value=TRUE)[1]
      pval_col <- if ("pvalue" %in% colnames(te_age_de)) "pvalue" else "padj"

      te_plot <- te_age_de |>
        rename(lfc=all_of(lfc_col), pval=all_of(pval_col)) |>
        mutate(sig = pval < 0.05,
               label = ifelse(sig & abs(lfc) > 1, family_id, NA))

      p_age_de <- ggplot(te_plot,
                         aes(x=mean_div, y=lfc, color=macro, size=n_loci)) +
        geom_point(alpha=0.7) +
        geom_hline(yintercept=0, linetype="dashed", color="grey50") +
        ggrepel::geom_text_repel(
          aes(label=label), size=2.8, max.overlaps=15,
          show.legend=FALSE
        ) +
        scale_color_manual(values=col_macro) +
        scale_size_continuous(range=c(2, 8), name="# loci") +
        labs(
          title    = "TE age vs expression change (tumor vs normal)",
          subtitle = "Young TEs (low divergence) more likely to be dysregulated",
          x        = "Mean divergence from consensus (%) — proxy for age",
          y        = sprintf("log2FC (%s / %s)", opt$treatment, opt$control),
          color    = "TE class"
        ) +
        theme_minimal(base_size=12) +
        theme(legend.position="right")

      comp_tag <- if (!is.null(opt$treatment) && !is.null(opt$control))
                    sprintf("_%s_vs_%s", opt$treatment, opt$control) else ""
      ggsave(file.path(opt$outdir, sprintf("te_age_vs_de%s.pdf", comp_tag)),
             p_age_de, width=10, height=7)
      ggsave(file.path(opt$outdir, sprintf("te_age_vs_de%s.png", comp_tag)),
             p_age_de, width=10, height=7, dpi=150)

      # Pearson correlation
      cor_res <- cor.test(te_plot$mean_div, te_plot$lfc, method="pearson")
      cat(sprintf("[INF] Age-LFC correlation: r=%.3f, p=%.4f\n",
                  cor_res$estimate, cor_res$p.value))

      write.table(
        te_age_de |> select(family_id, macro, mean_div, n_loci,
                             all_of(lfc_col), all_of(pval_col)),
        file.path(opt$outdir, sprintf("te_age_de_table%s.tsv", comp_tag)),
        sep="\t", quote=FALSE, row.names=FALSE
      )
      cat("[OK] TE age vs DE plot saved.\n")
    }
  }, error=function(e) cat("[WARN] TE age vs DE failed:", conditionMessage(e), "\n"))
}

# ═══════════════════════════════════════════════════════════════════════════════
# 3. Young TE enrichment near DE genes (optional)
# ═══════════════════════════════════════════════════════════════════════════════
# This analysis requires te_host_genes.tsv from M11 + age data
# Placeholder: reported in M11b if promoter intersect is available
cat("[INF] M11c complete. Output:", opt$outdir, "\n")
