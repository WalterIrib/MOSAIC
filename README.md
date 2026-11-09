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
| M06 | circRNA detection (screening) | DCC 0.5.0 | §4.2 |
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

## Installation and computational requirements

```bash
git clone https://github.com/WalterIrib/MOSAIC.git
cd MOSAIC
CONDA_CHANNEL_PRIORITY=flexible conda env create -f environment.yml   # main env "mosaic" (pinned versions)
conda activate mosaic
pip install -e .
mosaic --version
```

(or set it once with `conda config --set channel_priority flexible`; with strict
priority libmamba excludes the bioconda build of DCC)

Optional separate environments (not needed for the core pipeline):

```bash
conda env create -f environment_telescope.yml   # Telescope benchmark only (§4.8)
# HUMAnN3 (optional M20b) → conda env "humann3"
```

Binary/script tools set via `config.yaml`: Arriba 2.4.0 (M15).

**Software.** `environment.yml` pins the command-line tools and Python libraries
of the environment that produced the manuscript results (full export:
`versions.txt`). All published runs used **R 4.5.1** (system R; package
versions and build R in `r_versions.txt`); `environment.yml` pins `r-base=4.5.1`
and the same R package versions. Use R 4.5.1: under R 4.5.2 a few dozen TE loci
at |log2FC| = 1 change significance (see "Reproducibility of the manuscript
numbers").

**Hardware used for the published runs.** One desktop workstation: Intel Core
i5-13500 (20 threads), 32 GB RAM, no GPU. No cluster is required.

**Disk.** Reference data (`references/` used for the manuscript: genome, GTF,
STAR / HISAT2 / BWA indexes, DFAM, RepeatMasker, Kraken2 Standard DB, Arriba):
about 67 GB. The optional HUMAnN3 databases (ChocoPhlAn, MetaPhlAn, UniRef;
M20b only) add about 111 GB.

**Runtime** (chromatin-RNA dataset, 12 PE libraries, 20 threads; wall-clock from
the pipeline logs):

| Step | Wall-clock |
|------|-----------|
| M02 fastp, 12 samples | 0 h 25 min |
| M04 STAR, 12 samples | 7 h 30 min |
| M08 featureCounts, 36 runs (3 comparisons × 12 samples) | 0 h 36 min |
| M17 rMATS, one comparison (6 samples) | 0 h 03 min |
| M19 JACUSA2, one comparison (6 samples) | 1 h 01 min |

Steps without timestamps in their logs (DESeq2, bowtie2 TE/satellite counting,
`bedtools multicov`) are not listed; peak memory was not recorded for these modules.

### Reference data

All paths are set under `databases:` in the YAML (see `config_template.yaml`).
The set used for the manuscript (human GRCh38, GENCODE release 49) occupies
**about 67 GB**; the optional HUMAnN3 databases (M20b only) add about 111 GB.

| Component | Needed by | Where it comes from |
|-----------|-----------|--------------------|
| Genome FASTA `GRCh38.primary_assembly.genome.fa` and annotation `gencode.v49.annotation.gtf` | alignment, quantification and annotation modules | [GENCODE release 49](https://www.gencodegenes.org/human/release_49.html) (also fetched by `mosaic-refs`) |
| STAR index (~29 GB) | M04 | built locally from the FASTA and GTF (`STAR --runMode genomeGenerate`) |
| RepeatMasker table `rmsk/rmsk.txt.gz` (+ `rmsk/TE.gtf` derived from it) | M11, M12 | UCSC, hg38 `rmsk` track: <https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz> |
| Dfam library and TE / satellite pseudochromosomes with their bowtie2 indexes | M05, M09 | [Dfam](https://www.dfam.org) 3.9; pseudochromosomes built with `scripts/build_pseudochromosomes.py` |
| Kraken2 Standard database (~11 GB; the manuscript used the build of 2024-09-13) | M20 | prebuilt Kraken2/Bracken indexes: <https://benlangmead.github.io/aws-indexes/k2> |
| Arriba 2.4.0 database (`database/` of the Arriba release) | M15 | [Arriba v2.4.0 release](https://github.com/suhrig/arriba/releases/tag/v2.4.0) |
| BWA index | M23 | built locally (`mosaic-refs`) |
| HISAT2 index | not used by the pipeline (HISAT2 comparison in the manuscript only) | built locally (`mosaic-refs`) |
| HUMAnN3 databases (ChocoPhlAn, MetaPhlAn, UniRef) | M20b (optional) | `humann_databases --download …` (HUMAnN3 documentation) |

Assembly constraints in v1.1:

- **Other GENCODE releases of GRCh38** can be used by pointing `genome_fa` and
  `gtf` to them and rebuilding the STAR index.
- **The Arriba database is assembly-specific.** M15 passes the GRCh38 files of
  Arriba 2.4.0 by name (`blacklist_hg38_GRCh38_v2.4.0.tsv.gz`,
  `known_fusions_hg38_GRCh38_v2.4.0.tsv.gz`,
  `protein_domains_hg38_GRCh38_v2.4.0.gff3`), so fusion screening runs on GRCh38 only.
- **M12 needs the RepeatMasker table of the same assembly** as the genome FASTA
  and GTF: locus coordinates and TE classes are taken from it.
- **T2T-CHM13 is not supported in v1.1**: Arriba 2.4.0 ships no CHM13 database,
  and M12 has not been run or tested with a CHM13 RepeatMasker annotation.

### Known limitations (v1.1)

- **M20b (HUMAnN3 functional profiling) is experimental and not tested in this
  release.** It requires a separate `humann3` conda environment (HUMAnN3 is not
  compatible with the pinned R 4.5.1 / Bioconductor stack in `environment.yml`).
  M20 (Kraken2/Bracken taxonomic classification of unmapped reads) is fully
  supported and validated. To use M20b, create the `humann3` env following the
  [HUMAnN3 installation guide](https://huttenhower.sph.harvard.edu/humann/) and
  set `environments.humann: humann3` in `config.yaml`. If M20b is enabled but
  the env is unavailable, the pipeline logs a warning and continues; all other
  modules are unaffected.

- **`scripts/b1_sync_dcc.sh` contains hardcoded paths to the METACAR and
  DECECCA2019 canonical directories** (`/home/walter/Desktop/METACAR/mosaic_output`,
  `/home/walter/Desktop/DECECCA2019/mosaic_output`). The script is a one-off
  utility for syncing DCC circRNA results back to canonical trees after an
  external rerun; it is **not invoked by the pipeline runner** and is not called
  from any module. The canonical directories are `chmod a-w`, so an accidental
  invocation would fail harmlessly. Noted here so future maintainers are not
  surprised by the grep output.

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
- `trimming:` — fastp `-q min_phred` (default 20) and `-l min_length` (default 36);
  `cut_right` (default `false`) adds `--cut_right --cut_window_size --cut_mean_quality`;
  `pe_correction` (default `false`) adds `-c`; `disable_adapter_trimming` (default
  `false`) replaces `--detect_adapter_for_pe` with `-A`.
- `circrna.min_bsj_reads` — minimum back-splice junction reads for DCC (default 2).
- `project.random_seed` — seed for bowtie2 (M05), the M18 network layout and the
  `mosaic-refs` index builds (default 42).

Fixed in the code, not configurable (keys with these names were accepted but
ignored by earlier versions and have been removed from the templates):

- JACUSA2 `call-2` runs with its default statistical model (no `-u DirMult:PHRED`;
  formerly `thresholds.editing_model`).
- STAR coordinate sorting uses `--limitBAMsortRAM 27000000000` (~27 GB; formerly
  `memory_gb`).
- The RSeQC BED12 is reused from `mosaic_output/mapping/annotation/genes.bed12` or
  built from `databases.gtf` (formerly `databases.gtf_bed12`).
- circRNA detection uses DCC only (`databases.ciri2_path` was never used).
- The annotation version is whatever GTF `databases.gtf` points to
  (`project.annotation_version` was a label only).
- The M14 report title and format are fixed in `scripts/m14_report.qmd` (`report.*`).

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

# Run the smoke test (about one minute on the shipped synthetic data)
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

**1. Download reference files** (~67 GB with all indexes and Kraken2 DB; see "Disk" above):

| File | Source |
|------|--------|
| `GRCh38.primary_assembly.genome.fa` | [GENCODE GRCh38](https://www.gencodegenes.org/human/) |
| `gencode.v49.annotation.gtf` | [GENCODE v49](https://www.gencodegenes.org/human/) |
| DFAM TE/satellite pseudochromosomes | Generated by `scripts/build_pseudochromosomes.py` (see §4.2) |
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
locus / splicing / editing counts per comparison.

**Reproducibility of the manuscript numbers.** The manuscript counts were
produced before the first tagged release with the |log2FC| >= 1 threshold
evaluated as a strict inequality on values rounded to four decimals, whereas all
tagged releases evaluate the inclusive threshold on the full-precision value:
8, 29, 24 and 23 significant loci (adjusted P < 0.05) round to ±1.0000 (9, 66,
211 and 26 in total) in mito-ChaRNA vs ChaRNA, SEN_L vs EP, SEN_E vs EP and
SEN_L vs SEN_E. With release v1.1.8 in the environment created from
`environment.yml` (R 4.5.1), the significant sets of all gene-, TE family- and
satellite-level results are identical to the reported ones and the locus-level
sets differ only by the threshold rule (+4, +16, +13, +2 significant loci); log2
fold changes differ from the reported values by at most 0.0006 and adjusted P
values by at most 2% (relative), differences attributable to the numerical
libraries of the conda build. Under the system R 4.5.2 the gene-, family-
and satellite-level tables are byte-identical, whereas ashr-shrunk locus-level
fold changes differ by up to 0.04 log2 units (net +2, +16, +46, +2; fewer than
0.02% of tested loci, at most about 0.1% of the significant set), with identical
adjusted P values. `check_canonical.py` accepts [published, published + boundary
count] and passes under R 4.5.1; under R 4.5.2 SEN_E vs EP exceeds the band by
22 loci.

### Determinism and fixed parameters

- Seeds: `project.random_seed` (default 42) is passed to bowtie2 (M05), to the
  M18 network layout and to the `mosaic-refs` index builders.
- M18 community detection (`cluster_louvain` in `scripts/m18_network.R`) runs
  without a seed: community labels and membership can differ between runs.
- `mosaic-refs` (`mosaic/refs.py`): only the human reference is current
  (GENCODE v49 / Ensembl 115, as used for the manuscript). The other organisms
  still point to earlier releases: mouse GRCm39 GENCODE vM34 / Ensembl 111,
  rat mRatBN7.2 (rn7) Ensembl 111, zebrafish GRCz11 Ensembl 111.

---

## Archived versions

| Version | Archive |
|---------|---------|
| v1.1.2 (the version cited in the manuscript) | Zenodo, DOI [10.5281/zenodo.23979119](https://doi.org/10.5281/zenodo.23979119) |
| v1.1.8 | Zenodo, DOI [10.5281/zenodo.23983678](https://doi.org/10.5281/zenodo.23983678) |

The v1.1.8 row (and the DOI identifiers in `CITATION.cff`) were added on master
after the tag, once Zenodo had assigned the DOI: the README and `CITATION.cff`
inside the v1.1.8 archive list only the v1.1.2 DOI. That commit (41279a5)
changes only README.md and CITATION.cff; the pipeline code is that of the tag.

## Citation

> Sciaraffa N, Russo R, Zito F, Cusimano A, Arancio W. *MOSAIC: an automated RNA-seq
> platform integrating transposable element, RNA editing and splicing dynamics across
> cellular senescence and mitosis.* **Methods** (submitted, under review).

A `CITATION.cff` is provided. Source: <https://github.com/WalterIrib/MOSAIC> (MIT).

Software author: Walter Arancio (sole author of the code archived on Zenodo); the
manuscript authors are listed in the preferred-citation of CITATION.cff.

## Author

**Walter Arancio** — IRIB, CNR, Palermo, Italy — walter.arancio@irib.cnr.it
