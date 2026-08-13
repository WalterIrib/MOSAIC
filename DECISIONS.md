# MOSAIC — Technical Decisions

## Deferred fixes (targeted for v1.1.2)

### D1: M12 dead code — `stranded` variable never passed to m12_te_locus.py

**Location:** `mosaic/pipeline.py:1445–1446`

**What it does:** Computes `stranded = "no" if strand == "unstranded" else ("reverse" if strand == "RF" else "yes")`.
For FR libraries, this yields `"yes"` — which is incorrect (FR/dUTP should map to `"reverse"` in TEtranscripts convention).

**Why deferred:** The variable `stranded` is never passed to `m12_te_locus.py` (see `run_cmd` call at line 1453–1459; no `--stranded` argument). The script therefore ran with its own default strandedness behavior for all configs, producing identical results regardless of config strandedness. Zero divergence between buggy and fixed paper runs.

**Fix for v1.1.2:** Pass `--stranded stranded` to `m12_te_locus.py` and correct the FR→"reverse" mapping.

---

### D2: M19 latent lib-type — hardcoded `RF-FIRSTSTRAND` not passed from config

**Location:** `mosaic/m19_rna_editing.py:52` (default), `mosaic/pipeline.py:1814–1825` (call site)

**What it does:** `m19_rna_editing.py` defaults `--lib-type` to `RF-FIRSTSTRAND`. `pipeline.py` never passes `--lib-type` from config strandedness.

**Why deferred:** Both paper datasets (METACAR, DECECCA2019) are FR/dUTP. In JACUSA2 convention, `RF-FIRSTSTRAND` is the correct lib-type for dUTP/FR libraries (R2 on sense strand). The hardcoded default was therefore correct for all paper results.

**Latent failure mode:** An RF-strandedness dataset (non-dUTP stranded, R1 on sense) with `m19_editing: true` would receive `RF-FIRSTSTRAND` instead of the correct `FR-SECONDSTRAND`, producing erroneous editing site calls.

**Fix for v1.1.2:** Add `--lib-type` passthrough in `_run_m19_editing()` mapping MOSAIC strandedness → JACUSA2 lib-type (FR→RF-FIRSTSTRAND, RF→FR-SECONDSTRAND, unstranded→UNSTRANDED).

---

### D3: M14 module-state inferred from filesystem, not from config or pipeline state

**Location:** `scripts/m14_report.qmd:3704`

**What it does:** `.m20_enabled <- dir.exists(file.path(out, "m20_metatranscriptomics"))`. The report decides
whether to render the M20 section by checking whether the output directory exists on disk, not by reading
`m20_meta: true/false` from the config or from `.mosaic_state.json`.

**Observed failure mode:** When M14 was re-rendered for DECECCA2019 on 2026-08-07 (16:10–16:46), the
`m20_metatranscriptomics/` directory did not yet exist (the M20 De Cecco run completed later the same day).
The section therefore rendered the "disabled" fallback message despite the module being fully supported. The
HTMLs were promoted to canonical before the inconsistency was noticed. A subsequent re-render after
`m20_metatranscriptomics/` was populated corrected the output.

**Same family as:** D1/D2 — pipeline produces output silently inconsistent with what a reader of the config
would expect. D1/D2 affect quantitative results; D3 affects only report rendering completeness.

**Fix for v2.0 (scope model):** M14 should receive an explicit module-state map (derived from
`.mosaic_state.json` or from `Pipeline.state`) rather than probing the filesystem directly. This is the
broader "scope model" refactor: each module declares its outputs and status explicitly; downstream consumers
(including the report) read that declaration, not the disk.

---

### D4: Paper figures contain hardcoded numeric values not derived from canonical outputs at render time

**Location:** `DECECCA2019/make_validation_figures.py`

**What it does:** Several numeric values displayed in figures are written as Python literals in the plotting
script rather than read from the files that contain the authoritative source data:

- `n_sig = [16, 263, 0]` (Fig 2C) — TE family sig counts for SEN_E_vs_EP and SEN_L_vs_EP; source is
  `figure_data/b_te_macro_counts_*.tsv`, the same files Fig 3B reads at runtime. The Fig 2C hardcoded
  values were correct as of 2026-08-10 (16 = 13+2+1; 263 = 162+72+28+1), but will silently diverge if
  `b_te_macro_counts_*.tsv` is regenerated.

- `ef_vals = {'EP': 0.187, 'SEN_E': 0.221, 'SEN_L': 0.269}` and `"10,343 sites"` (Fig 3D) — mean editing
  frequencies and site count; source is `mosaic_output/ef_trajectory_3cond.tsv`. The title count was
  wrong (10,343 vs actual 10,351). Fixed 2026-08-10: Fig 3D now reads ef_trajectory_3cond.tsv at runtime.

**Observed failure mode:** Fig 3D title said "10,343 sites" while the source file had 10,351 rows. The SEN_E
EF label printed "0.220" instead of "0.221" (truncation of 0.2205… before the fix in the prior session).
Both errors were invisible to diff-based review because the values were plausible and close to correct.

**Same family as:** D3 — state (or data) inferred/hardcoded rather than received from an authoritative source.

**Fix for v2.0:** Figures should be generated entirely from canonical outputs with no hardcoded numeric
constants derived from those outputs. Where a value appears in multiple places (figure label, title, paper
text), it should flow from a single source file so that all representations stay in sync automatically.

---

### D5: M14 locus-level (M12) results located by directory convention, not by config

**Location:** `scripts/m14_report.qmd:185` (path construction), `mosaic/pipeline.py` (output layout)

**What it does:** The report QMD resolves M12 locus-level results via:

```r
f_new <- file.path(comp_dir, "m12", sprintf("%s_te_de_sig.tsv", comp))
```

where `comp_dir = file.path(out_dir, paste0("comparison_", comp))`. This assumes M12 writes its output
into a per-comparison `comparison_*/m12/` subdirectory. In practice (METACAR), M12 writes to a flat
`mosaic_output/m12/` directory shared across all comparisons, with file names disambiguated by
comparison prefix.

**Observed failure mode:** DataS4–S6 (METACAR) §6.1 showed "No significant TE loci found" and
"M12 locus-level results not yet generated" despite 25,489 (DataS4), 12,597 (DataS5), and 25,015
(DataS6) significant loci existing in `mosaic_output/m12/`. The per-comparison `m12/` subdirectory
did not exist. Fix applied 2026-08-13: created relative symlinks
`comparison_*/m12 → ../m12` in all three METACAR comparison directories.

**Same family as:** D3 — pipeline module state (or output location) inferred from filesystem layout
rather than from an explicit declaration in config or `.mosaic_state.json`.

**Fix for v2.0 (scope model):** M12 module should declare its output path in `.mosaic_state.json`;
M14 should read that declaration instead of constructing a path by convention. Alternatively,
M12 should write per-comparison subdirectories to match M14's expectation.
