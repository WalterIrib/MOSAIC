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
