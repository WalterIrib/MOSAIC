## v1.1.3 — Bug fix: numeric sample IDs coerced to character

### Bug fixed (B1 — numeric sample_id)

When `coldata$sample_id` is read by R as integer (e.g. `1881793`), using it
directly as a matrix/data.frame subscript causes either a crash
("undefined columns selected") or an OOB positional-indexing failure.
Added `as.character()` at every affected site.

**Files changed:**

- `scripts/m14_report.qmd` — 7 sites: biotype-wilcox, sat-wilcox-per-family-run,
  sat-aggregate-macro-data, and four M20-section vectors
  (`sids_meta`, `sample_cols`, `sample_cols_m`, `sample_cols_b`)
- `scripts/m20_metatranscriptomics.R` — 1 site: `sample_cols` used to filter
  the abundance matrix before rowSums and downstream plots

**Not changed:** `scripts/m16_stringtie_de.R` (fix already in v1.1.2 at
lines 298–299) and `scripts/m12_deseq2.R` (safe: `intersect(integer,
character)` returns character in R via `c(x[match(...)>0L], y[0L])`).

### Regression tests added

Four new test classes in `tests/test_strandedness_regression.py`:

- `TestB1_M14_SatWilcox` — column OOB crash pattern (treat_sids)
- `TestB1_M14_BiotypeWilcox` — OOB crash pattern (treat_s)
- `TestB1_M20_SampleCols` — column OOB crash pattern (sample_cols)
- `TestB1_M12_IntersectSafe` — confirms m12 correct without as.character()

All 22 tests pass.

### Documentation

- `DECISIONS.md`: added D4 (intersect safety proof) and D5 (m16 fix
  already in v1.1.2)
- `scripts/m14_report.qmd`: extended sat-wilcox callout note to clarify
  that the significantly DE satellite families in the main text come from
  the aggregate family-level DESeq2 (§3.4), not from the per-family test

---

### Does this change the paper results?

**No.** The canonical datasets (De Cecco 2019 and METACAR) use character
sample IDs (e.g. `SEN_L1`). The `as.character()` wrapping is a no-op for
character inputs. All results and figures in the submitted manuscript were
produced with v1.1.2 and are unaffected.

### Tag and citation

The submitted manuscript cites tag v1.1.2. That tag is unchanged. **v1.1.3 is not cited in the
manuscript** and is provided solely for reproducibility of pipelines with
numeric sample IDs (e.g. ENA/SRA accession numbers used directly as sample
identifiers).

Nothing has been submitted to Zenodo.
