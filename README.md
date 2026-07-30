# MOSAIC — Multi-Omic Sequencing Analysis and Integrated Characterization

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python 3.11](https://img.shields.io/badge/Python-3.11-3776AB?logo=python&logoColor=white)](https://www.python.org/)
[![R 4.5](https://img.shields.io/badge/R-4.5-276DC3?logo=r&logoColor=white)](https://www.r-project.org/)

MOSAIC is a modular, dependency-aware RNA-seq pipeline that turns paired-end reads
and a single YAML file into differential gene expression, transposable-element and
satellite biology (family **and** locus level), alternative splicing, A-to-I RNA
editing, non-canonical RNA, optional metatranscriptomics, and a self-contained
interactive HTML report per comparison. All core modules run unconditionally.

This is **v1.1**, a reproducibility- and robustness-hardened release of the pipeline
used for the manuscript. Scientific behaviour is unchanged; see
[`CHANGELOG.md`](CHANGELOG.md) for every fix.

---

## Modules and where they appear in the manuscript

All module methods are described in **§4.2** of the paper. The table maps each module
to its tool and to the Results/figures it feeds.

| Module | Function | Tool | Manuscript |
|--------|----------|------|-----------|
| M02 | Read trimming & QC (dedup/low-complexity **off** to keep repeat reads) | fastp 1.1.0 | §4.2, §4.3 |
| M04 | Two-pass genome alignment | STAR 2.7.10b | §4.2; Fig 1 |
| M05 | TE **&** satellite family quantification (DFAM consensus pseudochromosomes) | bowtie2 2.5.4 | §4.2, §4.5 |
| M06 | circRNA detection (screening) | CIRI2 2.0.6 | §4.2 |
| M07 | Strandedness inference | RSeQC 5.0.4 | §4.2 |
| M08 / M08b | Gene DE + factorial interaction; satellite family DE (+ Wilcoxon); QC/biotype | featureCounts 2.1.1 + DESeq2 1.50.2 (ashr) | §4.2, §4.4; Tab 2, Fig 2–3 |
| M10 | GO / KEGG enrichment | clusterProfiler 4.18.4 | §4.2, §4.4 |
| M11 / M11b–d | TE→host mapping, promoter/body TE, TE age, Fisher enrichment | bedtools 2.31.1 / R | §4.2 |
| M12 | **TE locus-level DE by direct overlap counting (no EM)** | bedtools 2.31.1 + DESeq2 | §4.2, §4.5 |
| M14 | Interactive HTML report (per comparison) | Quarto 1.6 | §4.2; Suppl. Data S1–S6 |
| M15 | Fusion transcripts (screening) | Arriba 2.4.0 | §4.2 |
| M16 | Novel transcript assembly | StringTie 2.1.7 + gffcompare | §4.2 |
| M17 | Alternative splicing (SE/MXE/A3SS/A5SS/RI) | rMATS 4.3.0 | §4.2, §4.4; Fig 4C |
| M18 | WGCNA co-expression + STRING v12 PPI | WGCNA 1.74 + STRINGdb 2.22.0 | §4.2 |
| M19 | A-to-I RNA editing (strand-aware) | JACUSA2 2.1.16 | §4.2, §4.4; Fig 3D |
| M20 / M20b | Metatranscriptomics (**optional, off by default**) | Kraken2 2.17.1 / HUMAnN3 | §4.2 |
| M23 | Non-reference TE-junction screening | bowtie2 / custom | §4.2 |

> **Note on M12 naming.** The internal module key is `m12_tetranscripts` for backward
> compatibility, but the method is **direct overlap counting with `bedtools`**, *not*
> expectation-maximization (TEtranscripts/Telescope). Multi-mapping reads may be
> counted at each alignment site, as stated in §4.5. See `CHANGELOG.md`.

---

## Installation

```bash
git clone https://github.com/WalterIrib/MOSAIC.git
cd MOSAIC
conda env create -f environment.yml      # main env "mosaic" (pinned versions)
conda activate mosaic
pip install -e .
mosaic --version
```

Optional separate environments (not needed for the core pipeline):

```bash
conda env create -f environment_telescope.yml   # Telescope benchmark only (§4.8)
# HUMAnN3 (optional M20b) → conda env "humann3"
```

Binary/script tools set via `config.yaml`: Arriba 2.4.0 (M15), CIRI2 2.0.6 (M06).

### Known limitations (v1.1.0)

- **M20b (HUMAnN3 functional profiling) is experimental and not tested in this
  release.** It requires a separate `humann3` conda environment (HUMAnN3 is not
  compatible with the pinned R 4.5.x / Bioconductor stack in `environment.yml`).
  M20 (Kraken2/Bracken taxonomic classification of unmapped reads) is fully
  supported and validated. To use M20b, create the `humann3` env following the
  [HUMAnN3 installation guide](https://huttenhower.sph.harvard.edu/humann/) and
  set `environments.humann: humann3` in `config.yaml`. If M20b is enabled but
  the env is unavailable, the pipeline logs a warning and continues; all other
  modules are unaffected.

---

## Configuration — a single YAML

Everything is driven by one config file; there are **no hardcoded paths**. Paths may
use `${ENV}` variables and `~` for portability:

```bash
export MOSAIC_REF=/data/references        # reference databases
export MOSAIC_DATA=/data/charna_fastq      # raw FASTQs
cp config_template.yaml config.yaml        # then edit, OR use a ready example:
mosaic check   configs/config_charna.yaml
mosaic run     configs/config_charna.yaml
```

Ready-made example: [`configs/config_charna.yaml`](configs/config_charna.yaml)
(the manuscript's chromatin-RNA experiment, ENA PRJEB112202).

Key knobs:

- `thresholds:` — `padj 0.05`, `lfc 1.0`, `min_reads_splicing 10` (≥10 junction
  reads for rMATS events; see §4.4 / CHANGELOG), editing coverage/score/EF.
- `continue_on_error: true` — a non-critical module failure (OOM in m18/m20/m23…)
  is logged and the run continues; backbone modules (m02/m04/m05/m08) are always
  fatal. Reruns resume from the `.mosaic_state.json` checkpoint.
- `modules:` — enable/disable each module; M20/M20b are **off by default**.

---

## Running

```bash
mosaic run    config.yaml                    # full pipeline
mosaic run    config.yaml --modules m17_rmats,m19_editing
mosaic status config.yaml                     # checkpoint status
```

Outputs per comparison: `deseq2/` (gene/TE/satellite DE TSVs), `m12/` (locus DE),
`m17_rmats/` (splicing), `m19_editing/` (editing), `go/`, `m18/`, and the
`m14_report.html`.

---

## Reproducibility & testing

### Functional smoke test — self-contained, no external data needed

The smoke test uses **synthetic PE reads** (generated by `wgsim`) on a tiny 500 kb
region of chr22 (22–22.5 Mb from GRCh38) that is **shipped with the repo** under
`test_data/synth_ref/`. No reference genome download, no raw FASTQs needed. It
verifies that every module produces non-empty, correctly-schema'd output, that
results are deterministic, and that the M14 report correctly names the comparison.

```bash
# One-time setup: build STAR/bowtie2 indexes + generate synthetic FASTQs
# (takes ~1 min, only uses tools already in the mosaic conda env)
bash test_data/generate_test_data.sh

# Run the smoke test (~20-30 min end-to-end)
conda run -n mosaic python tests/smoke_test.py
```

> **Important:** always invoke the smoke test via `conda run -n mosaic` (or from
> inside an activated `mosaic` environment). Running `python tests/smoke_test.py`
> directly — e.g. with the system Python or a different env — makes `featureCounts`
> unavailable in the subprocess PATH, causing `test_determinism_analysis_modules` to
> fail spuriously with "Module m08_deseq2 failed" even though the code is correct.
> pytest users: `conda run -n mosaic pytest tests/ -v`.

The synthetic data intentionally produces few significant DE genes and zero
splicing events — the goal is code-path verification, not biological validation.
See `.github/workflows/smoke.yml` for CI integration.

### Paper reproduction — requires full reference data and raw FASTQs

To reproduce the manuscript numbers on the complete METACAR dataset:

**1. Download reference files** (~50 GB):

| File | Source |
|------|--------|
| `GRCh38.primary_assembly.genome.fa` | [GENCODE GRCh38](https://www.gencodegenes.org/human/) |
| `gencode.v49.annotation.gtf` | [GENCODE v49](https://www.gencodegenes.org/human/) |
| DFAM TE/satellite pseudochromosomes | Generated by `scripts/build_dfam_pseudochromosomes.py` (see §4.2) |
| `rmsk/TE.gtf`, `rmsk/rmsk.txt.gz` | UCSC RepeatMasker track for hg38 |

Set paths in your config or via environment variable:
```bash
export MOSAIC_REF=/path/to/references
```

**2. Raw FASTQs** — the METACAR chromatin-RNA libraries are deposited at
ENA under accession **PRJEB112202** (see `configs/config_charna.yaml`). Download
with `mosaic run configs/config_charna.yaml` (module M00, or use `fastq-dump`/`fasterq-dump`).

**3. Numerical regression check** — after the full run, compare against manuscript values:
```bash
python check_canonical.py --mosaic-output /path/to/mosaic_output
```
Prints PASS/FAIL with expected-vs-observed for gene / TE-family / satellite /
locus / splicing / editing counts per comparison. TE-locus counts are expected to be
a few loci above published values (boundary-artefact fix; see CHANGELOG).

---

## Citation

> Sciaraffa N, Russo R, Zito F, Cusimano A, Deidda I, Arancio W. *MOSAIC: an automated
> platform for unbiased multi-omic RNA-seq analysis…* **Genomics** (2026), in press.

A `CITATION.cff` is provided. Source: <https://github.com/WalterIrib/MOSAIC> (MIT).

## Author

**Walter Arancio** — IRIB, CNR, Palermo, Italy — walter.arancio@irib.cnr.it
