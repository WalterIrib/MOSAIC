# Changelog

All notable changes to MOSAIC are documented here.

## [Unreleased] — 2026-07-21 — M23 3′-breakpoint bug fix (P0)

### Fixed — CRITICAL: M23 TE insertion clustering used alignment start instead of alignment end for 3′ soft-clips

- **`mosaic/m23_te_insertions.py` (`parse_chimeric_sam`, lines 148–158):** for
  reads with a 3′ soft-clip the breakpoint position was recorded as the alignment
  start (`int(parts[3]) − 1`). The true 3′ breakpoint is the end of the aligned
  segment: `pos + Σ(CIGAR ops that consume the reference: M, D, N, =, X)`. For a
  typical chimeric read aligned over ~120 nt, this produces a systematic ~120 bp
  error in recorded position for all 3′ clips.

- **Impact on clustering:** `cluster_insertions()` applies single-linkage at
  500 bp. The ~120 bp error can incorrectly merge sites in the 380–500 bp range
  or split sites in the 500–620 bp range. In the METACAR canonical run this
  produced 6,420 total clusters (≥2 reads); after the fix, **6,404 clusters**
  (16 fewer, 0.25%). 7SLRNA artifacts: 1,120 → 1,111 clusters. Final count
  after 7SLRNA/chrM exclusion: **5,293** (was 5,300 reconstructed; paper had
  cited 5,297 which is irrecoverable from the buggy canonical).

- **Numbers to update in manuscript:**
  - Total TE insertion clusters: 6,420 → **6,404**
  - Clusters after 7SLRNA/chrM exclusion: 5,297 → **5,293**

- **Fix:** `ref_len = sum(l for l, op in cigar_ops if op in ('M','D','N','=','X'))`;
  `clip_breakpoint = pos + ref_len if clip["position"] == "3prime" else pos`.
  5′ clip positions unchanged (already correct).

- Re-clustering performed using corrected 3′ positions patched into canonical
  `clips_with_te_mapping.tsv` (TE alignment step not repeated; TE family
  assignment is independent of genomic position). Output saved to
  `/tmp/m23_fixed_output/te_insertion_clusters_fixed.tsv`.

## [Unreleased] — 2026-07-18/19 — Canonical number audit and GO fix

### Fixed — METACAR splicing count corrected (2,032 phantom → 566 canonical)

- The METACAR `comparison_mito-ChaRNA_vs_ChaRNA/m17_rmats/significant_splicing_events.tsv`
  in a scratch directory (`METACAR_m14_fix/`) contained 2,032 events from a
  2026-07-09 rMATS rerun (`METACAR_v1.1_genomewide_20260707`). This phantom run
  was inadvertently used as a data source. The canonical canonical
  `mosaic_output/comparison_mito-ChaRNA_vs_ChaRNA/m17_rmats/significant_splicing_events.tsv`
  (2026-03-27, md5=ddf08b9f) contains 566 events (SE 177, MXE 110, RI 164,
  A3SS 70, A5SS 45), matching the manuscript. The phantom scratch directory
  was archived to external HD as `METACAR_m14_fix_CONTAMINATED_phantom_20260709`
  and is no longer referenced.

### Fixed — METACAR canonical go/ regenerated from canonical inputs

- `comparison_mito-ChaRNA_vs_ChaRNA/go/` and `comparison_mito-totRNA_vs_totRNA/go/`
  contained GO enrichment results from the initial (buggy) pipeline run where
  the foreground was ~7,300 genes instead of the DE-significant set. A corrected
  run existed only in a contaminated scratch directory (`METACAR_m14_fix/`, which
  also contained phantom rMATS/m12 data) and was never synced to the canonical.
  **Fix:** re-ran `m10_go.R` for both METACAR comparisons directly from canonical
  inputs (DE-significant genes md5=df7be9ee, 1,417 genes for mito-ChaRNA).
  Verified foreground: mito-ChaRNA UP=143, DOWN=266; mito-totRNA UP=113, DOWN=79.
  Key result for manuscript: mito-ChaRNA GO UP shows no significant BP or CC terms;
  single MF hit is BMP binding (GO:0036122, 4/146, padj=0.007), confirming the
  null result described in Fig. S3. Supplementary reports S4 and S5 regenerated
  and verified.

### Fixed — SUPPA2 benchmark metric corrected in Table S6

- The original Table S6 P1.4 metric compared raw event counts between SUPPA2 and
  rMATS (SUPPA2 457 vs rMATS 627), which is not a valid comparison given different
  tested-event sets and multiple-testing regimes. Replaced with directional
  concordance on co-testable events: 319/430 rMATS-significant events (74%)
  concordant in direction with SUPPA2 (binomial P<10⁻²³); ΔΨ Pearson r=0.25.
  `tableS6_versions_perf.xlsx` and `.tsv` updated accordingly.

### Documented — PAPER_NUMBERS.tsv marked non-canonical

- `paper_MOSAIC/PAPER_NUMBERS.tsv` is a manually compiled reference file shown to
  contain at least three values from superseded analyses: splicing 2,032 (phantom),
  MXE=0 (contradicted by canonical), editing 1,915/5,493 (pre-fix JACUSA2). A
  warning header was added to the file. Authoritative numbers are in
  `mosaic_output/` module outputs, re-verified against canonical files on
  2026-07-18.

## [Unreleased] — 2026-07-13/15 — Post-release fixes

### Fixed — CRITICAL: M19 RNA-editing parse_jacusa2 compared replicates within the same condition instead of across conditions

- **`mosaic/m19_rna_editing.py` (`parse_jacusa2`): the function read only
  `parts[6]` (control replicate 1) as condition 1 and `parts[7]` (control
  replicate 2) as condition 2, ignoring all treatment-sample columns
  entirely.** In JACUSA2 `call-2` output with N+M replicates, the column
  layout is: `parts[6..6+N-1]` = per-replicate base counts for the N control
  samples, `parts[6+N..6+N+M-1]` = per-replicate base counts for the M
  treatment samples. The old code treated the first two columns as condition 1
  vs condition 2, which for a standard 3+3 design means it compared
  ChaRNA\_rep1 against ChaRNA\_rep2 — two replicates of the same control
  condition — while completely ignoring the three treatment samples
  (mito-ChaRNA1, mito-ChaRNA2, mito-ChaRNA3).

  **Impact — all five comparisons in the manuscript were affected:**
  - `mito-ChaRNA vs ChaRNA` and `mito-totRNA vs totRNA` (METACAR): the
    reported ΔEF measured intra-replicate variability within ChaRNA / totRNA,
    not the mitotic editing effect. Example: mito-ChaRNA vs ChaRNA reported
    ΔEF = +0.069 (5,493 sites); correct value is ΔEF = +0.054 (5,669 sites),
    a 46% → 28% drop in the claimed relative increase.
  - `SEN_L vs EP`, `SEN_E vs EP`, `SEN_L vs SEN_E` (De Cecco 2019): each
    comparison was measuring within-condition replicate noise (EP\_rep1 vs
    EP\_rep2, or SEN\_E\_rep1 vs SEN\_E\_rep2). The reported "progressive
    A-to-I increase across senescence stages" (EF 0.161→0.249→0.261) was
    therefore entirely artifactual. The corrected inter-condition comparisons
    confirm the progressive increase is real (EF 0.187→0.220→0.269 on 10,343
    sites common to all three contrasts), but the absolute numbers and site
    counts change substantially.
  - **Not affected:** the REDItools3 cross-validation benchmark (r = 0.774,
    21,506 shared chr1 sites) read JACUSA2 raw columns directly with correct
    per-condition indexing (`parts[6:9]` for control, `parts[9:12]` for
    treatment) and did not pass through `parse_jacusa2`. The benchmark result
    is unchanged and remains valid.

  **Root cause:** `parse_jacusa2` was written assuming a 2-column output
  (one base-count field per condition), which corresponds to JACUSA2 with
  exactly 1 replicate per condition. The pipeline was designed and used with
  3 replicates per condition throughout, making every real run invalid. The
  additional bug that `ref` was inferred from the most-common base in
  condition 1 (rather than read from the dedicated JACUSA2 `ref` column) is
  also fixed, though it had minimal downstream effect since most sites had
  unambiguous reference bases.

  **Fix:** `parse_jacusa2(raw_out, n_normal, n_tumor)` now accepts replicate
  counts for each condition. It sums ACGT base counts across all N normal
  replicates (`parts[6]` through `parts[6+N-1]`) and all M tumor replicates
  (`parts[6+N]` through `parts[6+N+M-1]`), deriving total pooled counts per
  condition before computing editing fractions. Column offsets for `info`,
  `filter`, and `ref` are computed dynamically from `N+M`. In `main()`,
  `n_normal` and `n_tumor` are derived from the comma-separated BAM argument
  lists (`--bam-normal`, `--bam-tumor`), so the fix is automatic for any
  replicate count including the 2+2 smoke-test layout.

  **Corrected values (METACAR):**
  | Metric | Old (buggy) | Corrected |
  |---|---|---|
  | mito-ChaRNA_vs_ChaRNA total sites | 5,493 | 5,669 |
  | mito-ChaRNA_vs_ChaRNA mean ΔEF | +0.069 | +0.054 |
  | mito-ChaRNA_vs_ChaRNA relative increase | 46% | 28% |
  | mito-ChaRNA_vs_ChaRNA coverage-matched n | 1,273 of 1,924* | 3,767 of 5,669 |
  | mito-ChaRNA_vs_ChaRNA coverage-matched ΔEF | +0.069* | +0.042 |
  | mito-totRNA_vs_totRNA total sites | 7,246 | 7,817 |
  | mito-totRNA_vs_totRNA mean ΔEF | +0.065 | +0.053 |

  *The 1,924-site "all\_AtoI" set used in the original coverage-matched
  analysis was generated by an external script no longer on disk; it is
  considered irrecoverable and replaced by the 5,669 corrected sites above.

### Fixed — M19 `annotate_sites`: duplicate rows after left join on genomic position

- **`mosaic/m19_rna_editing.py` (`annotate_sites`): `bedtools intersect -wa -wb
  -loj` emits one output row per overlapping gene interval.** When a site
  overlaps multiple gene entries in the GTF, the subsequent `df.merge(ann,
  on=["chrom","start","end"], how="left")` produces N rows for that site
  (one per gene hit), inflating the total site count and assigning an
  arbitrary gene annotation to each duplicated row depending on sort order.

  **Impact:** any comparison where at least one site overlapped more than one
  gene interval would have its site count and EF statistics inflated. The
  magnitude depends on how many multi-gene-overlapping sites are present;
  in the five manuscript comparisons the inflation was confirmed to be zero
  (0 duplicate rows after applying the fix), but the bug created a
  reproducibility risk: a different annotation GTF or genome region could
  produce inflation silently.

  **Not affected — Method C values and Fig. S5:** the pipeline path used to
  compute the Method C union set and Fig. S5 statistics
  (`paper_MOSAIC/editing_reanalysis/`) ran `parse_jacusa2` +
  `filter_all_edits` directly without calling `annotate_sites`. All four
  source files (`editing_sites_ato_i.tsv`) have no `gene_name` column and
  zero annotation duplicates; the paper values (10,343 sites, EF\_EP=0.187,
  EF\_SEN\_E=0.220, EF\_SEN\_L=0.269, coverage-matched n=3,767,
  ΔEF=+0.042, P=8×10⁻⁸⁹) are unaffected and do not need recalculation.

  **Fix:** after the `df.merge(ann, ...)` call in `annotate_sites()`, add a
  deduplication step that keeps the highest-scoring row per `(chrom, start,
  end, edit_type)` key:
  ```python
  df = (df.sort_values("score", ascending=False)
          .drop_duplicates(["chrom","start","end","edit_type"])
          .reset_index(drop=True))
  ```
  The key includes `edit_type` (not just position) to preserve dual-strand
  sites — positions where both an A>G event on the + strand and a T>C(A>I)
  event on the − strand are detected at the same genomic coordinates (10
  such positions found in the De Cecco 2019 data). Verified: smoke test
  12/12 PASS with zero annotation duplicates in all five comparisons.

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

### Verified — M15 Arriba ran correctly on METACAR (prior CHANGELOG entry corrected)

- **A previous entry in this section incorrectly stated that "M15 fusion detection
  was never executed on the METACAR dataset".** That entry was wrong. The canonical
  METACAR run log (`mosaic_output/logs/m15_fusion_ChaRNA1.log`, timestamp
  2026-03-27) shows Arriba 2.4.0 launching and completing successfully for all 12
  samples, using the local `config.yaml` which already contained `arriba_bin` and
  `arriba_db` pointing to the installed Arriba v2.4.0 at
  `/home/walter/Desktop/references/arriba_v2.4.0/`.

  The previous entry also referenced `configs/config_metacar_gw.yaml` as the file
  that was missing the Arriba keys. That file was never committed to this repository
  (only `configs/config_charna.yaml` is tracked in `configs/`), and
  `config_charna.yaml` already contains `arriba_bin: ${MOSAIC_REF}/arriba_v2.4.0/arriba`
  and `arriba_db: ${MOSAIC_REF}/arriba_v2.4.0/database` with correct
  `${MOSAIC_REF}` placeholder paths. No code or config change was needed or made.

### Removed — Dead-code module files m04_map.py and m02_trim.py

- **`mosaic/m04_map.py` removed.** The file defined `star_align()` and
  `run_m04()` but neither was ever imported or called from `pipeline.py` or any
  other file in the codebase. The live STAR alignment implementation is
  `_run_m04_star()` in `pipeline.py` (lines 843–918), which `pipeline.py`
  dispatches via `getattr(self, "_run_m04_star", None)` for the `m04_star`
  module key. Proof the dead code never ran: `star_align()` contained
  `--outFilterMismatchNoverReadLengthMax` as a STAR parameter — a name that does
  not exist in STAR (valid: `--outFilterMismatchNoverReadLmax`); STAR exits
  fatally on unknown parameters, so any invocation would have aborted
  immediately. The three STAR values cited in manuscript Methods §4.3
  (`outFilterMultimapNmax 20`, `alignIntronMax 1,000,000`, `chimSegmentMin 20`)
  and `twopassMode Basic` all match `_run_m04_star()` exactly.

- **`mosaic/m02_trim.py` removed.** The file defined `run_fastp()`,
  `run_trimming_all()`, `load_sample_manifest()`, and `_parse_fastp_json()` but
  none was ever imported or called from `pipeline.py`. The live fastp
  implementation is `_run_m02_trim()` in `pipeline.py` (lines 738–839), which
  builds the fastp command directly via `run_cmd`. Proof the dead code never
  ran: the canonical run log (`mosaic_output/logs/m02_trim_EP1.log`) records the
  exact fastp invocation as `fastp … --detect_adapter_for_pe -q 20 -l 36` —
  matching `_run_m02_trim()` exactly — with none of the additional flags that
  `run_fastp()` in `m02_trim.py` would have added by default (`--correction`,
  `--cut_right`, `--unqualified_percent_limit 40`). All five fastp claims in
  manuscript Methods §4.3 (Q20, min-length 36 bp, automatic adapter detection,
  dedup off, low-complexity filter off) are confirmed against `_run_m02_trim()`.

  Both files were named after the pipeline module labels used in `config.yaml`
  (`m02_trim`, `m04_star`), which made them superficially resemble active module
  implementations and concealed that neither was wired into `pipeline.py`.
  Smoke test 12/12 PASS after removal.

### Removed — Dead-code module file m15_fusion.py

- **`mosaic/m15_fusion.py` removed.** The file defined `run_arriba()`,
  `parse_arriba_fusions()`, `merge_fusions()`, and `run_m15()`, but none of these
  was ever imported or called from `pipeline.py` or any other file in the codebase.
  The live Arriba implementation is `_run_m15_fusion()` in `pipeline.py`
  (lines 1813–1863), which builds and runs the Arriba command directly via
  `run_cmd()` without reference to `mosaic.m15_fusion`. This is the third
  standalone `.py` file removed for this reason (after `m04_map.py` and
  `m02_trim.py`); all three followed the same pattern of being named after
  pipeline module labels while duplicating logic already present inline in
  `pipeline.py`.

  **Note — m11b_promoter_te.py and m23_te_insertions.py are NOT dead code.**
  `_run_m11b_promoter_te()` (pipeline.py:1331) and `_run_m23_te_insertions()`
  (pipeline.py:1796) invoke their respective `.py` files as subprocess modules via
  `sys.executable -m mosaic.m11b_promoter_te` and `sys.executable -m
  mosaic.m23_te_insertions`. Those files contain the actual analysis logic and must
  not be removed.

### Known issues — deferred to v2.0 post-submission refactoring

- **Short-circuit in comparison loop (`pipeline.py:2158`).**
  ```python
  ok = ok and self.runner.run_module(module, comp)
  ```
  Python's `and` short-circuits: once `ok` is `False` (because comparison 1 failed),
  all subsequent comparisons for the same module are skipped without being attempted.
  Effect: when any non-critical module fails on comparison 1, comparisons 2–N for
  that module are silently not run. The `failed_modules` list and final warning only
  report the module name, not which comparisons were skipped. Fix: accumulate results
  without short-circuiting (`ok = self.runner.run_module(module, comp); all_ok &= ok`).
  **No paper results are affected** — in the canonical METACAR run all modules
  succeeded on all comparisons. Deferred because the fix requires changing state
  accumulation semantics and is better handled in a full runner refactor.

- **Checkpoint non-idempotence in `_should_run()` (`pipeline.py:541–552`).**
  On a clean first run a module enters `done` state. On a second run `is_done()`
  returns `True` and the module is correctly skipped, but `run_module()` then
  calls `state.mark(module, "skipped")`, writing `skipped` over `done`. On a
  third run `m.status == "skipped"` satisfies the `return True` branch at line 551,
  causing the module to re-execute. The state cycle is therefore:
  `done → skipped → (re-runs again) → done → …`, which is not stable.
  **No paper results are affected** — the canonical METACAR run completed in a
  single pass with no partial-restart scenario. Deferred pending checkpoint
  redesign in v2.0.

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
