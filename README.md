# MOSAIC — Multi-Omic Sequencing Analysis and Integrated Characterization

[![Python ≥3.10](https://img.shields.io/badge/Python-%E2%89%A53.10-3776AB?logo=python&logoColor=white)](https://www.python.org/)
[![R ≥4.2](https://img.shields.io/badge/R-%E2%89%A54.2-276DC3?logo=r&logoColor=white)](https://www.r-project.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

MOSAIC is a modular, automated pipeline for comprehensive multi-omic RNA-seq analysis.
Starting from raw paired-end FASTQ files and a single YAML configuration, it executes
more than 20 analysis modules (M02–M23) covering differential gene expression, transposable
element biology, alternative splicing, RNA editing, circular RNA detection, fusion genes,
co-expression networks, and metatranscriptomics. Results are delivered as interactive
HTML reports and structured TSV/CSV output files.

---

## Features

| Module | Name | Tools |
|--------|------|-------|
| M02 | Read trimming & QC | fastp |
| M04 | Two-pass alignment | STAR |
| M05 | TE quantification | featureCounts (pseudochromosome) |
| M06 | circRNA detection | CIRI2 |
| M07 | Strandedness inference | RSeQC |
| M08 | Differential gene expression | DESeq2 |
| M08b | QC plots & PCA | R/ggplot2 |
| M09 | Satellite repeat DE | DESeq2 |
| M10 | GO / KEGG / DO enrichment | clusterProfiler |
| M11 | TE–host gene integration | bedtools |
| M11b | TE in promoter regions | bedtools |
| M11c | TE evolutionary age analysis | R |
| M11d | DE genes × TE annotation | R |
| M12 | TE locus-level DE | TEtranscripts / DESeq2 |
| M14 | Interactive HTML report | Quarto |
| M15 | Fusion gene detection | Arriba |
| M16 | De novo transcript assembly | StringTie |
| M17 | Differential alternative splicing | rMATS |
| M18 | Co-expression network & PPI | WGCNA + STRINGdb |
| M19 | RNA editing (A-to-I and other) | JACUSA2 |
| M20 | Metatranscriptomics | Kraken2 / Bracken |
| M20b | Functional profiling | HUMAnN3 |
| M23 | TE somatic insertion detection | custom / bowtie2 |

Modules can be enabled or disabled individually in the YAML config file.
The pipeline supports single-comparison and multi-comparison (interaction model) designs,
and can download FASTQ files automatically from SRA/ENA.

---

## Requirements

### System tools
- [STAR](https://github.com/alexdobin/STAR) ≥ 2.7
- [bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/) ≥ 2.4
- [fastp](https://github.com/OpenGene/fastp) ≥ 0.23
- [rMATS](https://github.com/Xinglab/rmats-turbo) ≥ 4.1
- [JACUSA2](https://github.com/dieterich-lab/JACUSA2) ≥ 2.0
- [CIRI2](https://sourceforge.net/projects/ciri/) ≥ 2.0
- [Arriba](https://github.com/suhrig/arriba) ≥ 2.4
- [StringTie](https://ccb.jhu.edu/software/stringtie/) ≥ 2.2
- [Kraken2](https://ccb.jhu.edu/software/kraken2/) ≥ 2.1
- [HUMAnN3](https://github.com/biobakery/humann) ≥ 3.8
- [bedtools](https://bedtools.readthedocs.io/) ≥ 2.30
- [Quarto](https://quarto.org/) ≥ 1.4
- [featureCounts](https://subread.sourceforge.net/) (Subread ≥ 2.0)

### R packages (≥ 4.2)
DESeq2, WGCNA, clusterProfiler, STRINGdb, ggplot2, pheatmap, ggrepel, DT, RColorBrewer

### Python packages (≥ 3.10)
click, PyYAML, requests

---

## Installation

```bash
git clone https://github.com/WalterIrib/MOSAIC.git
cd MOSAIC
pip install -e .
```

Verify installation:

```bash
mosaic --version
mosaic check config.yaml
```

Reference databases (genome, GTF, RMSK, Kraken2 DB, etc.) can be downloaded automatically:

```bash
mosaic preflight config.yaml
```

---

## Quick start

**1. Initialise a new project**

```bash
mosaic init my_project
cd my_project
```

**2. Edit `config.yaml`**

```yaml
project:
  name: my_project
  wd: /path/to/output
  organism: human
  genome_build: hg38
  annotation_version: gencode_v49

library:
  layout: paired
  strandedness: auto

samples:
  - id: control_1
    condition: control
    R1: /data/control_1_R1.fastq.gz
    R2: /data/control_1_R2.fastq.gz
  - id: treated_1
    condition: treated
    R1: /data/treated_1_R1.fastq.gz
    R2: /data/treated_1_R2.fastq.gz

comparisons:
  - name: treated_vs_control
    treatment: treated
    control: control

modules:
  m02_trim:    true
  m04_star:    true
  m08_deseq2:  true
  m10_go:      true
  m14_report:  true
```

**3. Run the pipeline**

```bash
mosaic run config.yaml
```

Run specific modules only:

```bash
mosaic run config.yaml --modules m17_rmats,m19_editing
```

Check pipeline status:

```bash
mosaic status config.yaml
```

---

## Output

Each comparison generates a self-contained output directory containing:

- **`m14_report.html`** — interactive HTML report with volcano plots, PCA, TE heatmaps,
  splicing summaries, RNA editing landscape, network visualisations, and satellite plots
- **`deseq2/`** — normalised counts, DE results (TSV), MA/volcano plots
- **`m17_rmats/`** — rMATS output tables and summary plots
- **`m19_editing/`** — JACUSA2 editing sites and differential editing tables
- **`m20_metatranscriptomics/`** — Kraken2/Bracken species tables and diversity metrics
- **`integrated/`** — cross-module integrated matrices (TE × gene co-occurrence, etc.)

All tabular results are in TSV format, compatible with R and Python directly.

---

## Citation

If you use MOSAIC in your research, please cite:

> Sciaraffa N, Russo R, Zito F, Cusimano A, Deidda I, Arancio W. *MOSAIC: an automated
> platform for unbiased multi-omic RNA-seq analysis reveals mitotic regulation of splicing,
> transposable elements, and RNA editing in chromatin-associated RNA.*
> **Genomics** (2026). DOI: [10.1016/j.ygeno.2026.XXXXXX](https://doi.org/10.1016/j.ygeno.2026.XXXXXX) *(in press)*

A `CITATION.cff` file is provided for automated citation tools (Zenodo, GitHub, Zotero).

---

## Authors

**Walter Arancio**
Institute of Biomedical Research and Innovation (IRIB), CNR — Palermo, Italy
[walter.arancio@irib.cnr.it](mailto:walter.arancio@irib.cnr.it)

---

## License

This project is released under the [MIT License](LICENSE).
