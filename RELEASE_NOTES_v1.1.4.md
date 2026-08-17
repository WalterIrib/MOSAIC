## v1.1.4 — Test fix: use data.frame for biotype-wilcox OOB regression test

**Tag:** v1.1.4 — commit `41dc43d`
**Date:** 2026-08-17
**Supersedes:** v1.1.3 (e39731f)

### What changed

Three tests introduced in v1.1.3 used a `matrix` to represent `pct_df`,
but the actual pipeline object is a `data.frame`. The distinction matters:

| Object | Row OOB (large integer index) | Behaviour |
|--------|-------------------------------|-----------|
| `matrix` | `matrix[999L, col]` | ERROR: subscript out of bounds |
| `data.frame` | `data.frame[999L, col]` | `NA_real_` silently |

The tests were asserting the wrong (matrix) behaviour; the as.character()
fix in v1.1.3 remained correct. Tests rewritten with a real `data.frame`
built via `sapply(...) |> as.data.frame()`, matching the pipeline code.
Also fixed an off-by-one in the matrix dimension used in the intersect
test (`nrow=3` → `nrow=2`).

**Files changed:**
- `tests/test_strandedness_regression.py`:
  - `TestB1_M14_BiotypeWilcox`: three tests rewritten with `data.frame`
  - `TestB1_M12_IntersectSafe`: matrix dimension corrected to `nrow=2`

### Does this change pipeline results?

No. The production code in `m14_report.qmd` was correct since v1.1.3;
only the tests were wrong. All 22 tests pass on v1.1.4.

### Tag and citation

The manuscript cites **v1.1.2** (commit `e022cf6`), unchanged.
v1.1.4 is a test-suite maintenance release.
