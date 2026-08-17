## v1.1.7 — Fix: single-source strandedness resolution; explicit validation

**Tag:** v1.1.7 — commit `63bbe55`
**Date:** 2026-08-17
**Supersedes:** v1.1.6 (372ef2f)

### Bugs fixed

#### B4 — `_get_consensus_strandedness()`: silent fallback to "unstranded" on resumed run

**Trigger:** pipeline resumed with `strandedness: auto` in config and m07 state=done
(skipped). `_strandedness_config = "auto"` failed the `!= "auto"` guard;
`cfg._strandedness = {}` (never set by m07) failed the falsy guard; function
fell through to `return "unstranded"`.

**Impact:**
- **M08 featureCounts**: `-s 0` instead of `-s 2` for the 2 new comparisons
  (WT_vs_parental, AA_vs_parental) in the Lucibello v1.1.7 rerun.
- **M16 featureCounts**: would have gotten `-s 0` for all 6 comparisons (file
  is created once globally; new directory → featurecounts_transcripts.txt did
  not exist → would re-run with wrong flag).
- **M17 rMATS**: separate systematic bug — read `_strandedness_config` at
  parse-time instead of the resolved strandedness. `"auto"` fell through
  to `"fr-unstranded"` for **all** comparisons, including old ones.
  (Canonical METACAR/DECECCA2019 runs unaffected: their config declares
  `strandedness: FR` explicitly, not `auto`.)

#### Fix B — TSV as third resolution path

`_get_consensus_strandedness()` now reads `mapping/strandedness.tsv` (written
by m07) as a third resolution path, before any final fallback. A resumed run
finds the correct strandedness without re-running m07. Same approach already
used by M16 StringTie assembly (pipeline.py:1581-1590 since v1.1.0).

Priority order:
1. Explicit config value (`strandedness` ≠ `"auto"`) — validated.
2. In-memory map set by m07 (`cfg._strandedness`) — validated.
3. TSV on disk (`mapping/strandedness.tsv`) — validated.
4. `RuntimeError` — no silent fallback to any default.

#### Fix A — Explicit validation via `_validate_strand_label()`

New module-level function `_validate_strand_label(label, source)` raises
`RuntimeError` on any label not in `{"FR", "RF", "unstranded"}` — covering
`"undetermined"`, `"auto"`, or any future typo. Called from all three
resolution paths (config, strand_map, TSV) and from `strandedness_to_featurecounts`.

New module-level `_STRAND_CONSENSUS_THRESHOLD = 2/3`: if the majority label
covers < 67% of samples, `_consensus_from_map()` raises rather than silently
proceeding.

Non-unanimous consensus logs a warning with per-sample breakdown (AI2).

#### Fix A (M17 source) — M17 reads resolved strandedness, not parse-time config

`_run_m17_rmats` now calls `_strand_to_rmats_libtype(_get_consensus_strandedness())`
instead of reading `cfg.library._strandedness_config` directly. The old code
was a systematic bug independent of m07: `"auto"` → `.get("auto","fr-unstranded")`
= `"fr-unstranded"` for all comparisons regardless of run history.

New module-level `_strand_to_rmats_libtype(strand)` maps FR→fr-firststrand,
RF→fr-secondstrand, unstranded→fr-unstranded, and raises on any other label.

#### Also: dead code M12 removed; `strandedness_to_featurecounts` explicit

- M12: `stranded` variable (computed but never passed to `m12_te_locus.py locus-count`)
  removed.
- `strandedness_to_featurecounts`: replaced `else: return "2"` catch-all with
  explicit dict lookup via `_validate_strand_label`, so "undetermined" no longer
  silently returns "2" (FR).

### Tests added

`TestStrandednessResolution` in `tests/test_strandedness_regression.py`:

| Test | What it covers |
|------|---------------|
| `test_validate_raises_on_undetermined` | AI1: no silent pass-through |
| `test_validate_raises_on_auto` | AI1: "auto" is not a resolved label |
| `test_strand_to_rmats_*` | _strand_to_rmats_libtype unit tests |
| `test_consensus_warns_on_non_unanimous` | AI2: non-unanimous warns |
| `test_consensus_raises_on_weak_majority` | AI2: < 67% majority raises |
| `test_consensus_reads_tsv_when_m07_skipped` | Core regression (Fix B) |
| `test_consensus_raises_without_tsv_and_no_map` | Fix A: no silent unstranded |
| `test_m17_delegates_to_strand_to_rmats_libtype` | AI3 source regression check |

65 tests pass; 12 pre-existing smoke_test errors unchanged.

### Does this change pipeline results?

Yes, intentionally. Runs with `strandedness: auto` and m07 already done
will now correctly use the strandedness detected by m07 (via TSV) instead
of silently defaulting to "unstranded". Runs with explicit strandedness in
config are unaffected (path 1 resolves immediately, same as before).

### Tag and citation

The manuscript cites **v1.1.2** (commit `e022cf6`), unchanged.
v1.1.7 is not cited in the manuscript.

### Rerun scope for Lucibello (after applying this fix)

State entries to clear before resuming (see AI5 in session notes):
- `m08_deseq2/WT_vs_parental`, `m08_deseq2/AA_vs_parental` (wrong `-s 0`)
- `m08b_qc`, `m10_go`, `m11d_te_de` for the same 2 comparisons
- `m17_rmats` for all 4 old comparisons (`fr-unstranded`)

`m16_stringtie/parental_vs_pBabe` is state=running → will retry automatically.
`featurecounts_transcripts.txt` does not exist → featureCounts will re-run
with correct `-s 2`.

Do NOT delete: `merged.gtf`, `gffcmp.*`, per-sample quant GTFs — all valid.
