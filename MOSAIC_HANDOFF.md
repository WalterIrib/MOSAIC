# MOSAIC — Project Handoff Document
# Generated: 2026-03-19
# For: Claude Code session continuation

## WHAT IS MOSAIC

MOSAIC (Multi-Omic Analysis of Small-RNA, AI-assisted Integration and Characterization)
is a complete RNA-seq analysis pipeline focused on:
- Canonical gene differential expression (DESeq2)
- Repetitive element (TE/satellite) quantification and DE
- circRNA detection and DE
- Metatranscriptomics (Kraken2/Bracken + HUMAnN3)
- RNA editing (JACUSA2)
- Fusion genes (Arriba)
- Alternative splicing (rMATS)
- De novo transcripts (StringTie)
- Network analysis (STRINGdb + WGCNA)
- TE somatic insertions (M23)
- HTML report (Quarto)

Key biological focus: chrNA (chromatin-associated RNA) metafase vs interfase, TE biology preservation
(NEVER use low-complexity filter or poly-X trimming — removes TE/satellite reads)

## SYSTEM

- OS: Ubuntu 25.10
- Machine: ~/Desktop/mosaic1/ (project code)
- Output: ~/Desktop/mosaic/mosaic_output/
- conda env: mosaic (Python 3.11, R 4.5.1)
- conda env: humann3 (Python 3.9, HUMAnN3 + MetaPhlAn)
- References: ~/Desktop/references/

## PROJECT STRUCTURE

```
~/Desktop/mosaic1/
├── mosaic/              # Python package
│   ├── cli.py           # Main CLI: mosaic init/check/preflight/run/compare/download/db
│   ├── pipeline.py      # Pipeline orchestrator with all module runners
│   ├── config.py        # Config schema + load_config()
│   ├── db_manager.py    # Database download/verify manager
│   ├── m02_trim.py      # fastp trimming
│   ├── m04_map.py       # STAR alignment
│   ├── m05_te_count.py  # TE quantification (bowtie2 pseudochromosomes)
│   ├── m06_circrna.py   # CIRI2 circRNA detection
│   ├── m11b_promoter_te.py
│   ├── m15_fusion.py    # Arriba fusion detection
│   ├── m19_rna_editing.py # JACUSA2
│   └── m23_te_insertions.py
├── scripts/             # R scripts
│   ├── m08_deseq2.R
│   ├── m08b_qc_plots.R  # PCA, heatmap, strandedness QC
│   ├── m10_go.R
│   ├── m11_te_host_genes.py
│   ├── m11c_te_age.R
│   ├── m11d_te_de_integration.R
│   ├── m12_te_locus.py
│   ├── m14_report.qmd   # Main Quarto report
│   ├── m16_stringtie_de.R
│   ├── m17_rmats.R
│   ├── m18_network.R
│   ├── m19_editing_plots.R
│   └── m20_metatranscriptomics.R
├── config.yaml          # Active project config
├── setup.py
└── tests/

~/Desktop/mosaic/mosaic_output/
├── trimmed/             # fastp output + trim_stats.tsv
├── mapping/bam/         # STAR BAM files (4 samples)
│   ├── SRR19388325/     # normal
│   ├── SRR19388327/     # normal
│   ├── SRR19388330/     # tumor
│   └── SRR19388331/     # tumor
├── integrated/          # te_counts_clean.tsv, gene_type_map.tsv
├── comparison_tumor_vs_normal/
│   ├── deseq2/          # DE results, normalized_counts, PCA
│   ├── go/              # GO/KEGG enrichment
│   ├── m11d/            # TE-DE integration
│   ├── m17_rmats/       # Alternative splicing
│   ├── m18/             # STRINGdb + WGCNA
│   └── report/          # Final HTML report
├── m11/                 # te_host_genes.tsv, te_host_universe.txt
├── m11b/                # te_promoter_genes.tsv
├── m11c/                # te_age_summary.tsv
├── m12/                 # TEtranscripts locus count
├── m15_fusions/         # Arriba fusion genes
├── m16_stringtie/       # StringTie de novo
├── m19_editing/         # JACUSA2 RNA editing
├── m20_metatranscriptomics/ # Kraken2/Bracken + HUMAnN3
├── m23_te_insertions/   # TE somatic insertions
└── .mosaic_state.json   # Pipeline state tracking

~/Desktop/references/
├── GRCh38.primary_assembly.genome.fa
├── gencode.v49.annotation.gtf
├── star_index/
├── dfam/pseudochromosomes/
│   ├── consensus_TE.fa + .bt2 indices
│   └── consensus_satellite.fa + .bt2 indices
├── rmsk/rmsk.txt.gz + TE.gtf
├── kraken2_db/
├── humann3_db/ (metaphlan_db, chocophlan, uniref)
└── arriba_v2.4.0/ (binary + database)
```

## CLI USAGE

```bash
conda activate mosaic
cd ~/Desktop/mosaic1

# Full pipeline
mosaic run config.yaml --comparison tumor_vs_normal

# Specific modules
mosaic run config.yaml --modules m08_deseq2,m10_go,m14_report --comparison tumor_vs_normal

# New comparison without re-aligning
mosaic compare config.yaml --treatment normal --control tumor

# Resume from module
mosaic run config.yaml --from-module m08_deseq2 --comparison tumor_vs_normal

# Force re-run
mosaic run config.yaml --modules m14_report --comparison tumor_vs_normal --force

# Check databases
mosaic db check config.yaml

# Download FASTQs
mosaic download config.yaml --source ena --accession SRR19388325

# Status
mosaic status config.yaml
```

## SAMPLES (test dataset)

| sample_id | condition |
|-----------|-----------|
| SRR19388325 | normal |
| SRR19388327 | normal |
| SRR19388330 | tumor |
| SRR19388331 | tumor |

Test data: ~49M reads per sample (downloaded from ENA) — usato solo per validare il pipeline.
Coldata: ~/Desktop/mosaic/mosaic_output/coldata.tsv

## DATASET REALE (articolo)

Studio: chrNA (chromatin-associated RNA) — ciclo cellulare metafase vs interfase

| Confronto | Nome | Descrizione |
|-----------|------|-------------|
| 1 | chrNA_meta_vs_inter | chrNA metafase vs chrNA interfase |
| 2 | total_meta_vs_inter | RNA totale metafase vs RNA totale interfase |
| 3 | chrNA_meta_vs_inter_bg | chrNA meta vs inter con RNA totale come background |

Il confronto 3 richiede design statistico speciale (background correction).
I nomi "tumor_vs_normal" nel codice/output sono provvisori dal test dataset — non rispecchiano la biologia reale.

## CURRENT STATUS (2026-03-27)

Pipeline METACAR operativa. Analisi in corso.

Ultimi fix committati:
- commit 016d392: strand_fc mapping RF→1/FR→2, apeglm→ashr, n_te_prom_loci guard
- commit 7027602: _run_m06_circrna aggiunto, _get_consensus_strandedness legge
  library._strandedness_config (era ignorato → featureCounts usava -s 0 su lib FR),
  fallback default cambiato da "RF" a "unstranded"

Dataset METACAR (reale, per articolo): /home/walter/Desktop/METACAR/
- 12 sample: ChaRNA/mito-ChaRNA/totRNA/mito-totRNA × 3 repliche
- strandedness: FR (sequencing.strandedness: FR in config)
- featureCounts: -s 2 (reversely stranded) su tutti e 3 i confronti ✓

Stato confronti METACAR:
- mito-ChaRNA_vs_ChaRNA: m08✓ m08b? m10? m11d? m14?
- mito-totRNA_vs_totRNA: m08✓ m08b✓ m10✓ m11d✓ m14?
- ChaRNA_interaction: m08✓ (-s 2, ashr) m08b? m10? m11d? m14?
  → 48 geni con interazione significativa (padj<0.05, LFC>1)

In esecuzione ora: m08b+m10+m11d+m14 per tutti i confronti mancanti.

## NOTE OPERATIVE

- coldata.tsv era in /tmp/ (si cancella al riavvio) → spostato in
  ~/Desktop/mosaic/mosaic_output/coldata.tsv e config.yaml aggiornato
- Se la pipeline non trova i sample, verificare che coldata esista

## KEY DESIGN DECISIONS

1. **TE quantification**: pseudochromosome approach (Bowtie2 --very-sensitive-local)
   NOT TEtranscripts (used only for locus-level in M12)
   NEVER filter low-complexity or poly-X (removes satellite reads TTAGGG etc.)

2. **Size factors**: from stable genes (CV<0.3, baseMean>10)
   NOT from all genes (TE expression skews normalization)

3. **strandedness**: auto-detect con RSeQC per sample; questo dataset = unstranded

4. **Multi-comparison**: each comparison in own subdirectory
   State tracked per-module per-comparison in .mosaic_state.json

5. **Multi-env**: humann3 env for HUMAnN3/MetaPhlAn
   Arriba binary direct (not conda — htslib conflict)

## PRODUCTION CHECKLIST (key items)

| Module | Test param | Production param |
|--------|-----------|-----------------|
| DESeq2 | padj=0.5, lfc=0.5 | padj=0.05, lfc=1.0 |
| rMATS | --statoff | remove --statoff |
| Editing | min-cov=3 | min-cov=5 |
| Kraken2 | k2_minusb | k2_standard |
| HUMAnN3 | DEMO db | full ChocoPhlAn+UniRef90 |

(Pipeline runner already uses production values via thresholds in config.yaml)

## REPORT STRUCTURE (m14_report.qmd)

Sections: 1 QC → 2 DE genes → 3 DE RS (TE) → 4 GO → 5 circRNA →
6 TE locus → 7 StringTie → 8 rMATS → 9 Fusion genes →
10 Network → 11 RNA Editing → 12 Metatranscriptomics →
13 TE Somatic → 14 Methods → 15 Session Info

Parameters passed via Quarto -P flags:
- out_dir: ~/Desktop/mosaic/mosaic_output
- comparison: tumor_vs_normal
- treatment: tumor
- control: normal
- padj_threshold: 0.05
- lfc_threshold: 1.0

## GIT STATUS

```
~/Desktop/mosaic1/ è un git repo (branch master)
Ultimo commit: "pipeline: fix strandedness→featureCounts mapping and shrinkage estimator" (016d392)
```
