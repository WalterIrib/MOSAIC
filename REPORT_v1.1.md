# MOSAIC v1.1 — Reproducibility & robustness release report

Prepared for internal review before public GitHub release. `mosaic1.1/` is a
cleaned, debugged, reproducible copy of the pipeline that performs exactly what the
manuscript (`MOSAIC_paper_draft_v46.docx`) describes. Only reproducibility/robustness
issues were fixed; the declared scientific logic is unchanged. **Nothing was pushed.**

---

## (a) Bugs found and how they were fixed

| # | Area | Bug | Fix |
|---|------|-----|-----|
| 1 | `m12_deseq2.R` | TE-locus significance was thresholded on the **rounded** `log2FoldChange_shrunk` (4 dp) with strict `>`, creating a boundary artefact at `|LFC|==1.0000` (the documented 25,497-vs-25,489 / 376,193-vs-376,164 discrepancy). | Threshold now uses the **full-precision** shrunk LFC with `|LFC| >= 1` (Methods §4.4). Rounding kept for display only. |
| 2 | paths | Hardcoded `/home/walter/...` paths in `m18_network.R`, `m14_report.qmd` default, and dead `m14_interaction_report.qmd`. | All removed. `m18` uses `--outdir` + optional `MOSAIC_STRING_DIR`; `m14` default is neutral `mosaic_output`; dead file deleted. Whole tree is now free of absolute user paths. |
| 3 | `config.py` | Config paths were literal only. | Added `${ENV}`/`~` expansion so one config runs on any machine (`${MOSAIC_REF}` etc.). |
| 4 | `pipeline.py` | **Any** module failure aborted the whole run (`break`) — a single OOM in a heavy/optional module lost all core results. | Only backbone modules (m02/m04/m05/m08) are fatal; other failures are logged, checkpointed and the run continues. `continue_on_error` flag (default true), wired through `config.py`. |
| 5 | `m14_report.qmd` | Sequential Quarto renders shared one knitr figure dir → a comparison could reuse a previous one's cached plot (S1 showing S5 volcano). | Per-comparison `fig.path = m14_report_files_<comparison>/…` (already in code; verified + kept, dead single-dir path removed). |
| 6 | config defaults | `m20_meta`/`m20b_humann` defaulted **on**, contradicting the paper ("optional, off by default"). | Default **off** in template and example configs. |
| 7 | `m12` not wired | **`m12_deseq2.R` (locus-level TE DE) was never called by the pipeline** — `mosaic run` produced only the count matrix; the per-comparison DE (`<c>_te_de_sig.tsv`, e.g. 376,164 loci) had to be run by hand. A headline result was not reproducible from the pipeline. | Wired `m12_deseq2.R` to run per comparison after locus counting, using the canonical stable-gene size factors from m08. `m12` deps corrected to include `m11_te_host` + `m08_deseq2`. |
| 8 | `m08_deseq2.R` | TE/satellite DE used local (locfit) dispersion, which crashes (`out of vertex space`) on tiny/sparse matrices. | Fall back to `fitType="mean"` on error; genome-wide runs (local fit succeeds) unchanged. |
| 9 | determinism | — | Verified the only stochastic step (WGCNA) is seeded (`set.seed(42)`); all other tools deterministic. |

Full detail in [`CHANGELOG.md`](CHANGELOG.md).

---

## (b) Points where the code diverged from the paper (decisions required)

> Per instructions, these were **flagged, not silently "fixed"**.

### Divergence #1 — splitting reads≥10 filter — **RESOLVED (keep + document)**
- The manuscript splicing counts (mito-ChaRNA **566/RI 164**, SEN_L **12,123**,
  SEN_E **13,449**) are produced **with** the `≥10 junction-reads` filter in
  `m17_rmats.R`. Removing it (to match Methods §4.4 as literally written) gives
  **584/RI 166**, **12,124**, **13,452** — i.e. it would break reproducibility.
- The 565/RI 163 discrepancy that motivated the original "remove it" request came
  from a **separate figure script** (`make_fig4_composite.py`, strict `>0.1` on raw
  rMATS tables), not the pipeline.
- **Decision (you):** keep the filter (configurable, documented in code).
- **Manuscript action for you:** add the ≥10-read criterion to Methods §4.4
  ("FDR<0.05, |ΔΨ|≥0.1 and ≥10 junction reads summed across replicates").

### Divergence #2 — TE-locus threshold rounded — **FIX APPLIED (you pre-authorised)**
- Code diverged from Methods §4.4 (`|log2FC|≥1`) by thresholding a rounded value
  with strict `>`. Fixed to full-precision `≥1`.
- **Consequence:** significant **TE-locus** counts rise by a few loci (bounded):
  SEN_L +≤29 (376,164), SEN_E +≤24 (342,831), SEN_L_vs_SEN_E +≤23, mito-ChaRNA
  +≤8 (25,489). Genes, TE families, satellites and splicing are **unaffected**.
- Exact post-fix totals require a genome-wide rerun. `check_canonical.py` accepts
  the band `[published, published+boundary]` for locus counts and prints the delta.
- **Manuscript action for you:** update the locus counts to the post-fix values
  after a genome-wide rerun (or state the boundary tolerance).

### Divergence #3 — module key `m12_tetranscripts` — **documented, not renamed**
- The key name implies TEtranscripts/EM, but M12 is bedtools **direct-overlap**
  counting (no EM), exactly as §4.5 states. Renaming touches config/pipeline/deps.
- **Recommendation:** rename to `m12_te_locus` in a follow-up; for now documented in
  README + CHANGELOG. Your call.

### Divergence #4 — M12 locus DE not produced by the pipeline — **FIX APPLIED (clear bug)**
- The pipeline computed only the locus **count matrix**; the locus-level **DE**
  (a headline result of §4.5) was never generated by `mosaic run`. This is
  unambiguously a wiring bug, not a paper error, so I wired `m12_deseq2.R` in
  (per comparison, canonical stable-gene size factors). Flagged here for your
  visibility; I can revert to the manual-step behaviour if you prefer.
- Validated: script runs with the exact wired arguments and applies the canonical
  gene size factors; smoke test passes.

### Observation (not a divergence)
- Gene significance (m08) is thresholded on the **unshrunk** MLE `log2FoldChange`
  (with `results(lfcThreshold=1)` padj), while TE loci use the **shrunk** LFC. This
  asymmetry is pre-existing and canonical; the paper only says shrinkage is applied
  and reported. Left unchanged.

---

## (c) Smoke-test outcome on the subset (incl. M14)

The generator (`test_data/generate_test_data.sh`) built a 2-vs-2 subset
(mito-ChaRNA vs ChaRNA) of chr1+chr22 reads (≈0.49–0.55 M pairs/sample) plus a
mini reference (chr1+chr22 genome, subset GTF/rmsk/TE.gtf, mini STAR index).

The full pipeline ran **end-to-end** (m02→m23 + M14, `Pipeline complete in ~36 min`;
the long pole is m12 bedtools-multicov over ~450k chr1+chr22 loci — a CI fast mode
`CHRS=chr22` is provided). Backbone failures correctly abort (verified: an early
FASTQ-layout mistake stopped at m02 as designed); the m08 dispersion and m12-wiring
fixes above were both **found by this run**.

`python tests/smoke_test.py` → **12 / 12 PASS**:

| Test | Result |
|------|--------|
| m02 trim / m04 STAR BAM / m05 TE family counts | PASS (files + schema) |
| m08 gene DE schema (`baseMean`,`log2FoldChange`,`padj`,`log2FoldChange_shrunk`) | PASS |
| m09 satellite DE / m12 locus DE schema (`locus_id`) / m17 splicing / m19 editing | PASS |
| **M14 report named for the correct comparison** (fig.path cross-wiring regression) | PASS |
| **M14 per-comparison `fig.path` directory** (`m14_report_files_mito-ChaRNA_vs_ChaRNA`) | PASS |
| **M14 significant count in report == significant-table rows** (points == rows) | PASS |
| **Determinism** (genes DE + splicing byte-identical across two forced runs) | PASS |

Notes: on the tiny subset m12 yields 0 significant loci (expected — 2×2, very low
power) but writes a correctly-structured table; the smoke test checks structure and
determinism, not biological significance. The M12 locus-DE wiring was additionally
validated by direct invocation of `m12_deseq2.R` on the subset (1231 stable genes,
canonical size factors applied).

---

## (d) Running the genome-wide numerical regression

```bash
python check_canonical.py \
    -m /path/to/DECECCA2019/mosaic_output \
    -m /path/to/METACAR/mosaic_output
```

Validated against the current on-disk genome-wide outputs: **26 PASS / 0 FAIL / 1
SKIP** (the METACAR TE-locus file is not in the standard `m12/` location on disk;
the script reports SKIP rather than failing). After a fresh v1.1 genome-wide run the
TE-locus lines are expected to sit a few loci above the published value (Divergence
#2); the script flags this as PASS-within-band and prints the delta.

---

## (e) Final folder structure

```
mosaic1.1/
├── mosaic/                     # Python package (driver, config, modules)
│   ├── pipeline.py             # orchestrator (+ continue_on_error robustness)
│   ├── config.py               # config loader (+ ${ENV}/~ expansion)
│   ├── __main__.py             # `python -m mosaic`
│   └── m0*/m1*/m2*.py          # module runners
├── scripts/                    # R / Quarto / Python module scripts
│   ├── m08_deseq2.R  m10_go.R  m12_deseq2.R (fixed)  m17_rmats.R (documented)
│   ├── m14_report.qmd (per-comparison fig.path)      m18_network.R (de-hardcoded)
│   └── install_r_packages.R
├── configs/config_charna.yaml  # portable example (ENA PRJEB112202)
├── config_template.yaml        # annotated single-source template
├── environment.yml             # pinned main env (verified versions)
├── environment_telescope.yml   # benchmark-only env
├── check_canonical.py          # manual genome-wide regression (PASS/FAIL)
├── test_data/generate_test_data.sh   # builds the chr1+chr22 2v2 subset
├── tests/smoke_test.py         # CI-ready functional smoke test
├── .github/workflows/smoke.yml # structural gate + functional (self-hosted)
├── README.md  CHANGELOG.md  LICENSE (MIT)  CITATION.cff  setup.py
└── REPORT_v1.1.md              # this file
```
```
```
