# Changelog

All notable changes to MOSAIC are documented here.

## [Unreleased] — 2026-07-13/14 — Post-release fixes

### Fixed — GO foreground (M10)

- **`pipeline.py`: M10 GO enrichment was using all tested genes as foreground
  instead of only DE-significant genes.** `m10_go.R` was called with
  `genes_{comp}_DE.tsv` (~78,691 genes) as `--de-genes`. The correct input is
  `genes_{comp}_DE_significant.tsv` (padj<0.05, |LFC|>1 only). Because
  `m10_go.R` uses the file content directly as foreground with no internal
  filter, this inflated the foreground 15–30× and produced spurious GO term
  associations. Fixed by changing the `de_file` path in `_run_m10_go()`.
  Impact: supplementary GO figures (Fig. S3) and any Results/caption text
  referencing GO term counts are affected; corrected HTML reports regenerated
  for all METACAR and De Cecco 2019 comparisons.

### Fixed — M14 report: M12 path mismatch and interaction fallback messages

- **`m14_report.qmd`: M12 locus files were looked up in `out/m12/` (global
  output root) instead of `comp_dir/m12/` (per-comparison subdirectory).**
  The pipeline writes `{comp}_te_de_sig.tsv` and `{comp}_te_de.tsv` under
  `comparison_{comp}/m12/`, not under the global `mosaic_output/m12/`. This
  caused a false-negative "not yet generated" for M12 locus data in all
  comparisons even when the files existed. Fixed by changing `file.path(out,
  "m12", ...)` to `file.path(comp_dir, "m12", ...)` for both `tet_sig` and
  `tet_all`.

- **`m14_report.qmd`: M17/M18/M19 fallback messages said "run M18 first" or
  "not available" on interaction comparisons**, implying a missing execution
  rather than an inapplicable analysis. Added `is_interaction` guard to show
  "Not applicable to interaction comparisons (M17 rMATS / M18 / M19 JACUSA2
  requires a two-group contrast)" for rMATS (two instances), WGCNA hub genes,
  WGCNA modules, STRING PPI, and all three RNA-editing subsections.

### Fixed — M15 Arriba not configured for METACAR

- **`configs/config_metacar_gw.yaml` was missing `arriba_bin` and `arriba_db`**
  in the `databases` section, so M15 fusion detection was never executed on the
  METACAR dataset despite the Arriba binary and database being present and
  fully functional (confirmed on De Cecco 2019). Added the two entries pointing
  to the installed Arriba v2.4.0. Note: `config_template.yaml` and
  `configs/config_charna.yaml` already had these keys with correct
  `/path/to/...` / `${MOSAIC_REF}/` placeholders.

---

## [1.1.0] — 2026-07-06 through 2026-07-09 — Reproducibility & robustness release

This release is a cleaned, debugged, and reproducible version of the pipeline used
for the manuscript. It fixes **reproducibility and robustness** issues only; the
declared **scientific behaviour is unchanged** (direct-overlap TE counting, stable
canonical-gene size factors, DFAM pseudochromosomes, unconditional execution of
all core modules).

The work documented here proceeded in chronological order: m12 locus DE wired
automatically → initial genome-wide run (revealed gene_type_map.tsv/m10_go
divergences) → strandedness diagnosis and fix → pipeline.py bug
(self.out/self.runner.out, found during strandedness re-run) → canonical number
re-verification (splicing root-cause, TE families tolerance) → functional smoke
test (GAPs 1–8) → robustness fixes on edge cases.

---

### Fixed — strandedness detection for Total RNA ribo-depleted libraries

- **Auto-detection of strandedness (`m07_strandedness`) produced incorrect results
  for Total RNA ribo-depleted (rRNA-depletion) libraries**, including the ChaRNA and
  De Cecco 2019 datasets. RSeQC `infer_experiment.py` reported RF≈0.35 / FR≈0.39,
  below the 0.7 majority-vote threshold, and the pipeline called `unstranded`.
  The correct strandedness for dUTP-based prep is `FR` in MOSAIC config notation
  (`strandedness: FR` → featureCounts `-s 2`, reversely stranded). Note: MOSAIC's
  `FR` label corresponds to fr-firststrand/dUTP in featureCounts and RSeQC
  terminology; it should not be confused with the `RF` label (which maps to `-s 1`,
  fr-secondstrand, and is incorrect for dUTP libraries).

  Root cause: the 0.7 threshold is calibrated for poly-A selected libraries, which
  give strong strand signal because all reads come from mature mRNA. Ribo-depleted
  Total RNA includes abundant intronic reads, which are unspliced and map in both
  orientations, diluting the strand-specific signal below 0.7 even though the
  library is directional.

  Fix: explicit `strandedness: FR` in the ChaRNA and De Cecco configs (overrides
  auto-detection). The pipeline auto-detection logic and threshold are **not**
  changed — the 0.7 value is appropriate for poly-A libraries (the majority use
  case); ribo-depleted experiments should set strandedness explicitly.

  Effect: featureCounts gene counts (m08, the official gene DE input) changed from
  unstranded (`-s 0`) to dUTP-specific (`-s 2`), producing correct strand-specific
  read assignments for the ribo-depleted library.

---

### Fixed — canonical result verification and canonical counts

- **Splicing count for mito-ChaRNA vs ChaRNA corrected: 566 → 2032 (RI 164 → 735).**
  The value 566 was produced at the canonical commit (818c8484) because
  `--variable-read-length` was absent from the pipeline's rMATS invocation at
  that commit — it was not added to the `mosaic run` code path until this release.
  Reads trimmed by fastp to slightly unequal lengths were rejected by rMATS
  post-alignment read-length filtering. With `--variable-read-length` included in
  the rMATS call, the correct count is 2032 significant splicing events (RI=735).
  All other comparisons (SEN_L vs EP: 12,123; SEN_E vs EP: 13,449) are unaffected.
  The 566 count reflects the default pipeline behaviour at commit 818c8484, not a
  deviation from it.

  The value 2032 is now the canonical anchor in `check_canonical.py`.

- **TE family DE counts: tolerance ±5 families added to `check_canonical.py`
  for mito-ChaRNA vs ChaRNA (canonical 444 families).** DESeq2 local
  dispersion fitting (locfit) is sensitive to floating-point differences in the
  NB MLE across R patch versions. A handful of families near the padj = 0.05
  significance boundary can shift in or out of significance depending on whether
  the environment uses R 4.5.1 vs 4.5.2. The ±5 tolerance prevents spurious CI failures on these
  boundary cases without masking genuine regressions (the tolerance is tight
  relative to the 444 total).

  The canonical anchor (444) is retained as the published value; `check_canonical.py`
  reports the delta and pass/fail against the tolerance band.

---

### Fixed — reproducibility

- **TE-locus significance threshold uses full-precision shrunk LFC**
  (`m12_deseq2.R`). The significance mask tests `|log2FoldChange_shrunk| >= lfc`
  on the full-precision value; the rounded-to-4-decimals copy written to the
  output table is for display only. Using the rounded value for thresholding
  would create a boundary artefact: loci with a true |LFC_shrunk| in the interval
  (1.0, 1.00005) would round to 1.0000, which fails a `> 1.0` test and would
  incorrectly be excluded.

  The first committed version of `m12_deseq2.R` (added to the predecessor repo
  in April 2026, after canonical commit 818c8484) filtered on the already-rounded
  `log2FoldChange_shrunk` column. The v1.1 script separates a full-precision
  `lfc_shrunk_raw` variable for thresholding from the rounded output column, fixing
  this artefact. At commit 818c8484 itself no `m12_deseq2.R` existed — the m12
  DESeq2 step was run by hand using an ad-hoc script not in version control — so it
  cannot be established whether the published canonical loci counts used rounded or
  full-precision LFC for thresholding.

  `check_canonical.py` carries a defensive tolerance band (`loci_boundary`): up to
  +29 (SEN_L vs EP), +24 (SEN_E vs EP), +23 (SEN_L vs SEN_E), +8
  (mito-ChaRNA vs ChaRNA), representing the maximum number of loci that could
  plausibly shift near the threshold boundary under the rounding artefact. The
  genome-wide run reproduced the canonical published loci counts exactly (delta=+0
  for all verified comparisons), confirming that no loci in these datasets fall in
  the critical boundary zone in practice.

- **Hardcoded absolute paths removed.** All paths now come from `config.yaml`.
  Removed `/home/walter/...` defaults from `m18_network.R` (STRING DB search;
  now `--outdir`-derived plus optional `MOSAIC_STRING_DIR`) and from
  `m14_report.qmd` `out_dir` parameter default (now the neutral `mosaic_output`).

- **Portable config paths.** `config.py` now expands `${ENV}` variables and a
  leading `~` in all config strings, so the same config runs on any machine
  after exporting e.g. `MOSAIC_REF`.

- **Per-comparison figure directories in the HTML report.** `m14_report.qmd`
  sets a comparison-specific `knitr` `fig.path`
  (`m14_report_files_<comparison>/…`). This eliminates a cross-render
  contamination bug in which sequential Quarto renders sharing one figure
  directory could reuse a previous comparison's cached plot (e.g. an S1
  volcano showing S5 data). Verified by the smoke test.

- **Deterministic seeds.** The only stochastic step, WGCNA (`m18_network.R`),
  is seeded (`set.seed(42)`, from `project.random_seed`). DESeq2, rMATS,
  STAR and bowtie2 are deterministic for a fixed input and thread count.

- **Removed dead code.** Deleted the unused `m14_interaction_report.qmd`
  (the pipeline renders only `m14_report.qmd`, with
  `comparison_type: interaction` for the factorial model), which also carried
  a hardcoded path.

---

### Fixed — reproducibility gap (M12 locus-level DE was not wired in)

- **`m12_deseq2.R` (locus-level TE differential expression) was never called by
  the pipeline.** `mosaic run` produced only the locus count matrix
  (`m12/locus_counts.tsv`); the per-comparison DESeq2 step that yields
  `comparison_<c>/m12/<c>_te_de_sig.tsv` (the manuscript's TE-locus results,
  e.g. 376,164 loci for SEN_L vs EP) had to be run by hand. It is now invoked
  automatically per comparison at the end of `m12_tetranscripts`, using the
  canonical stable-gene size factors from m08 (`all_genes_counts.tsv` +
  `stable_genes.tsv`). `m12`'s dependencies were corrected to include
  `m11_te_host` (provides the rmsk BED) and `m08_deseq2` (provides size
  factors).

---

### Fixed — pipeline.py internal reference bugs

- **`self.out` vs `self.runner.out` reference error in `Pipeline`.** The outer
  `Pipeline` class stores the output root as `self.runner.out` (via
  `ModuleRunner`). One code path used `self.out` which did not exist on the
  `Pipeline` object; fixed to `self.runner.out`.

- **`gene_type_map.tsv` generated proactively.** This biotype map (built from
  the GTF, consumed by m08b QC and downstream modules) was generated lazily
  inside the m04/m08 step. If any module that needed it ran before m08
  completed (e.g. in a `--modules` partial run), the map was absent. It is now
  generated at pipeline start, before the per-comparison loop.

---

### Fixed — config loading

- **`load_coldata()` crashed with `IsADirectoryError` when `sequencing.coldata`
  was unset** (`mosaic/config.py`). `Path("")` resolves to `"."` (the current
  directory), which exists but is not a file, so the `not p.exists()` guard did
  not catch it and `open(".")` raised `IsADirectoryError`. Fixed to
  `not p.is_file()`, which correctly treats an unset or empty `coldata` as
  "no sample table provided" and skips loading.

---

### Removed

- **`mosaic-refs list` (and all `mosaic-refs` subcommands)** — the legacy
  `mosaic-refs` CLI group is kept as a stub for backward compatibility but no
  longer has subcommands registered. Use `mosaic db list` (and
  `mosaic db check` / `mosaic db download`) instead.

---

### Fixed — robustness

- **A single module failure no longer aborts the whole pipeline.** Only the
  backbone modules (`m02_trim`, `m04_star`, `m05_te_count`, `m08_deseq2`) are
  fatal on error. Any other module (heavy/optional/screening — `m18` WGCNA,
  `m20`/`m20b` metatranscriptomics, `m23` insertions, `m16`, `m06`, `m15`)
  that fails, e.g. out-of-memory, is logged, recorded in the checkpoint state
  and the run continues. Set `continue_on_error: false` in `config.yaml` to
  make all failures fatal. Reruns resume from the `.mosaic_state.json`
  checkpoint and retry only failed modules.

- **Metatranscriptomics off by default.** `m20_meta`/`m20b_humann` now default
  to `false` in the config template and example configs.

- **Robust dispersion fitting for TE/satellite family DE** (`m08_deseq2.R`).
  TE and satellite family-level DE now has a three-level fallback chain for
  DESeq2 dispersion estimation:
  1. `fitType="local"` (locfit) — always used on genome-wide data; succeeds
     in all real runs.
  2. `fitType="mean"` (constant dispersion) — fallback if locfit fails (e.g.
     `newsplit: out of vertex space` on very small/sparse matrices in tests).
  3. Gene-wise dispersion estimates set directly (`estimateDispersionsGeneEst`
     + `dispersions(dds) <- mcols(dds)$dispGeneEst` + `nbinomWaldTest`) —
     second fallback if even the constant fit fails. DESeq2 itself recommends
     this path in its error message. Produces statistically valid (if
     conservatively estimated) p-values; does not fabricate values.
  On genome-wide data level 1 always succeeds; levels 2–3 are unreachable.
  Level 3 was added to handle the degenerate case found in the synthetic smoke
  test (satellite counts have near-zero variance across 4 wgsim samples).

- **Gene DE table written before TE/satellite DE** (`m08_deseq2.R`). The
  canonical gene results (`genes_*_DE.tsv`, computed in memory before this
  point) are now written to disk before `run_te_de()` is called. Previously
  an unhandled error in satellite DE could exit the script before the gene
  table was saved. The reordering is pure I/O; no statistical values change.

- **`splicing_summary.tsv` schema preserved at zero events** (`m17_rmats.R`).
  When rMATS produces no events in any category (header-only output files),
  `bind_rows` on an empty list returns a 0×0 tibble with no column names.
  `write_tsv` on this produced a 1-byte file (bare newline) with no header,
  causing downstream schema checks to fail. Now an explicit zero-row
  `data.frame` with the correct seven columns is substituted before writing.
  On real data (events > 0 in at least one category) the `if` branch is never
  entered and behaviour is unchanged.

- **`editing_summary.tsv` written on the "no A-to-I sites" exit path**
  (`m19_rna_editing.py`). When JACUSA2 finds candidate sites but all are
  filtered out (score / EF thresholds), the Python module exited with
  `sys.exit(0)` before reaching the `compute_summary()` call that writes
  `editing_summary.tsv`. Added a zero-valued summary write (7 rows, same
  schema as the real summary) before the early exit. On real data `df_ai` is
  never empty (thousands of A-to-I sites pass filters), so this path is
  unreachable and behaviour is unchanged.

---

### Added — functional smoke test: self-contained, no external data required

The smoke test now uses **synthetic paired-end reads** (generated by `wgsim`)
on a 500 kb region of chr22 (22–22.5 Mb from GRCh38, renamed `chrSmoke`) that
is **shipped with the repo** under `test_data/synth_ref/`. No reference genome
download, no raw FASTQs from ENA, no `MOSAIC_REF` required. All 12 tests pass
(`12/12 PASS`, exit 0) with the conda env and nothing else.

Eight gaps (GAPs 1–8) identified during dry-run testing were resolved:

| GAP | Severity | Fix |
|-----|----------|-----|
| **GAP 1** | **BLOCKER** | Self-contained synthetic dataset: `synth_ref/` (genome.fa, annotation.gtf, TE.gtf, rmsk.txt.gz, dfam/ FASTAs, 2.6 MB total) committed to repo; `generate_test_data.sh` rewritten to build STAR/bowtie2/bwa indexes and generate FASTQs via wgsim with no external inputs |
| GAP 2 | CRITICAL | `smoke_test.py` derived `WD` from a hardcoded path; now reads `project.wd` from `config_test.yaml` via `_get_wd()` — test is valid at any install path |
| GAP 3 | BLOCKER | Missing R/Bioconductor packages added to `environment.yml` after full audit: `bioconductor-org.hs.eg.db` (m10 GO), `bioconductor-genomicfeatures` (m08 TPM fallback), `r-pheatmap` (m14 heatmap), `r-kableextra` (m14 HTML tables), `r-dt` (m14 interactive tables), `bwa=0.7.19` (m23 TE insertions) |
| **GAP 4** | **MAJOR** | `r-base` updated `4.5.1 → 4.5.2`; all pinned Bioconductor/CRAN packages verified to have `r45` build tags on bioconda/conda-forge at this version |
| GAP 5 | MINOR | `__version__` corrected `0.1.0 → 1.1.0` in `mosaic/__init__.py` |
| GAP 6 | MINOR | `pytest` added to `environment.yml` (was missing; `pytest tests/` failed immediately) |
| GAP 7 | MINOR | `generate_test_data.sh` hardcoded env name `mosaic` replaced with `${CONDA_DEFAULT_ENV:-mosaic}` |
| GAP 8 | INFO | `blast` documented as a known transitive conda dependency (not listed explicitly; arrives via bioconductor packages) |

GAPs 1 and 4 are the most substantial: together they resolve the full
self-containment requirement ("anyone who clones from GitHub can run the smoke
test with only the conda env").

---

### Added

- `environment.yml` (main env) and `environment_telescope.yml` (benchmark)
  with versions verified from the live environment.
- `configs/config_charna.yaml` example (portable `${MOSAIC_REF}`/`${MOSAIC_DATA}`).
- `test_data/synth_ref/` — 2.6 MB synthetic reference (chrSmoke 500 kb,
  GENCODE v49 subset, RepeatMasker TE.gtf subset, DFAM consensus FASTAs).
  Enables a fully self-contained smoke test with no external downloads.
- `tests/smoke_test.py` — 12-test functional smoke test; runs as
  `pytest tests/ -v` or standalone `python tests/smoke_test.py`.
- `check_canonical.py` — manual genome-wide numerical regression against the
  paper's reported values (genes, TE families, satellites, loci, splicing,
  editing dEF per comparison); supports `±families_tol` tolerance band.
- `scripts/install_r_packages.R` (BiocManager alternative to conda for the
  R stack).
- `tests/test_config.py` — 23-test unit suite for `mosaic/config.py` (brought
  forward from the predecessor repo and updated for the v1.1 API).

---

### Known issues / limitations

- **M20b (HUMAnN3 functional profiling) is not supported in the v1.1.0 conda
  environment.** HUMAnN3 conflicts with the pinned R 4.5.x / Bioconductor stack
  in `environment.yml` and must be installed in a separate `humann3` env. The
  pipeline handles this gracefully: if M20b is enabled but the `humann3` env is
  absent or `humann` is not importable, a warning is logged and M20b is skipped;
  all other modules continue normally. M20 (Kraken2/Bracken taxonomic
  classification) is fully supported and validated. See README for setup
  instructions.
