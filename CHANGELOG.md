# Changelog

All notable changes to MOSAIC are documented here.

## [1.1.0] — 2026-07-06 — Reproducibility & robustness release

This release is a cleaned, debugged, and reproducible version of the pipeline used
for the manuscript. It fixes **reproducibility and robustness** issues only; the
declared **scientific behaviour is unchanged** (direct-overlap TE counting, stable
canonical-gene size factors, |log₂FC|≥1 & padj<0.05 thresholds, strand-aware
editing, DFAM pseudochromosomes, unconditional execution of all core modules).

### Fixed — reproducibility

- **TE-locus significance was thresholded on a rounded fold change** (`m12_deseq2.R`).
  `log2FoldChange_shrunk` was rounded to 4 decimals *before* the `|LFC| > 1` test,
  creating a boundary artefact at exactly `|LFC| == 1.0000`: loci with a true
  fold change just above 1 were dropped. The significance mask now uses the
  **full-precision** shrunk fold change with `|LFC| >= 1` (matching Methods §4.4);
  the written table keeps the rounded column for display only.
  *Effect (locus level only; genes/families/satellites/splicing unaffected):* a
  small upward correction of the significant-locus counts — up to +29
  (SEN_L vs EP), +24 (SEN_E vs EP), +23 (SEN_L vs SEN_E), +8 (mito-ChaRNA vs
  ChaRNA). Exact post-fix totals require a genome-wide rerun. See README /
  `check_canonical.py`.

- **Hardcoded absolute paths removed.** All paths now come from `config.yaml`.
  Removed `/home/walter/...` defaults from `m18_network.R` (STRING DB search;
  now `--outdir`-derived plus optional `MOSAIC_STRING_DIR`) and from the
  `m14_report.qmd` `out_dir` parameter default (now the neutral `mosaic_output`).

- **Portable config paths.** `config.py` now expands `${ENV}` variables and a
  leading `~` in all config strings, so the same config runs on any machine after
  exporting e.g. `MOSAIC_REF`.

- **Per-comparison figure directories in the HTML report.** `m14_report.qmd` sets
  a comparison-specific `knitr` `fig.path` (`m14_report_files_<comparison>/…`).
  This eliminates a cross-render contamination bug in which sequential Quarto
  renders sharing one figure directory could reuse a previous comparison's
  cached plot (e.g. an S1 volcano showing S5 data). Verified by the smoke test.

- **Deterministic seeds.** The only stochastic step, WGCNA (`m18_network.R`), is
  seeded (`set.seed(42)`, from `project.random_seed`). DESeq2, rMATS, STAR and
  bowtie2 are deterministic for a fixed input and thread count.

- **Removed dead code.** Deleted the unused `m14_interaction_report.qmd`
  (the pipeline renders only `m14_report.qmd`, with `comparison_type: interaction`
  for the factorial model), which also carried a hardcoded path.

### Fixed — robustness

- **A single module failure no longer aborts the whole pipeline.** Only the
  backbone modules (`m02_trim`, `m04_star`, `m05_te_count`, `m08_deseq2`) are
  fatal on error. Any other module (heavy/optional/screening — `m18` WGCNA,
  `m20`/`m20b` metatranscriptomics, `m23` insertions, `m16`, `m06`, `m15`) that
  fails, e.g. out-of-memory, is logged, recorded in the checkpoint state and the
  run continues with the remaining core results. Set `continue_on_error: false`
  in `config.yaml` to make all failures fatal. Reruns resume from the
  `.mosaic_state.json` checkpoint and retry only failed modules.

- **Metatranscriptomics off by default.** `m20_meta`/`m20b_humann` now default to
  `false` in the config template and example configs, matching the manuscript
  ("only metatranscriptomics is optional and off by default").

- **Robust dispersion fitting for TE/satellite DE** (`m08_deseq2.R`). The
  family-level DE uses DESeq2 local (locfit) dispersion fitting, which can fail
  with `newsplit: out of vertex space` on very small/sparse count matrices (e.g.
  a chromosome-subset test set with few satellite families). Both TE and
  satellite DE calls now fall back to the constant `fitType="mean"` fit on error.
  On genome-wide data the local fit succeeds, so reported results are unchanged;
  this only rescues degenerate small inputs (found via the smoke test).

### Fixed — reproducibility gap (M12 locus-level DE was not wired in)

- **`m12_deseq2.R` (locus-level TE differential expression) was never called by
  the pipeline.** `mosaic run` produced only the locus count matrix
  (`m12/locus_counts.tsv`); the per-comparison DESeq2 step that yields
  `comparison_<c>/m12/<c>_te_de_sig.tsv` (the manuscript's TE-locus results, e.g.
  376,164 loci for SEN_L vs EP) had to be run by hand. It is now invoked
  automatically per comparison at the end of `m12_tetranscripts`, using the
  canonical stable-gene size factors from m08 (`all_genes_counts.tsv` +
  `stable_genes.tsv`). `m12`'s dependencies were corrected to include
  `m11_te_host` (provides the rmsk BED) and `m08_deseq2` (provides size factors).
  This is the module that also carries the rounding-artefact fix above.

### Investigated — NOT changed (behaviour is intentional)

- **Splicing junction-read filter (`reads ≥ 10`).** rMATS events must be supported
  by ≥10 junction reads summed across replicates (`m17_rmats.R`, configurable via
  `thresholds.min_reads_splicing`). This filter is **required** to reproduce the
  manuscript splicing counts (mito-ChaRNA vs ChaRNA: 566 events, RI = 164;
  SEN_L vs EP: 12,123; SEN_E vs EP: 13,449). It was initially suspected to be the
  cause of a 565/RI163 discrepancy, but that discrepancy came from a *separate
  figure-plotting script* (`make_fig4_composite.py`, which used a strict `>0.1`
  cutoff on the raw rMATS tables), not from the pipeline. The filter is retained
  and now documented in code.
  **Manuscript action:** Methods §4.4 should state the ≥10-read criterion
  explicitly (currently lists only FDR<0.05 & |ΔΨ|≥0.1).

### Added

- `environment.yml` (main env) and `environment_telescope.yml` (benchmark) with
  versions verified from the live environment.
- `configs/config_charna.yaml` example (portable `${MOSAIC_REF}`/`${MOSAIC_DATA}`).
- `test_data/` generator and a CI-ready functional smoke test (`tests/`).
- `check_canonical.py` for manual genome-wide numerical regression against the
  paper's reported values.
- `scripts/install_r_packages.R` (BiocManager alternative to conda for the R stack).
