
## M20 Metatranscriptomics (Kraken2/Bracken)
### Parametri test → produzione
- min-reads:   1  → 10   (reads minime per taxon)
- min-samples: 1  → 2    (campioni minimi con reads)
- Database: k2_minusb (8GB) → k2_standard (70GB) per maggiore sensibilità
  wget https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20240904.tar.gz
- Aggiungere --outReadsUnmapped Fastx al comando STAR principale
- Nota: Homo sapiens domina unmapped fraction — è normale
- Con dati reali cercare: virus oncogenici (EBV/HPV/HBV taxid specifici),
  microbioma tumorale, funghi opportunisti

## M23 TE Somatic Insertions
### Parametri test → produzione
- min-clip:    25  → 30   (lunghezza minima soft-clip)
- min-support: 1   → 2    (reads minime per cluster)
- min-mapq:    0   → 0    (mantenere 0 - chimeric reads STAR hanno MAPQ basso)
- Aggiungere --outReadsUnmapped Fastx al STAR principale (già nel checklist)
- Fix minore: colonna samples/strand incollate - da correggere in m23_te_insertions.py
- Nota: 0 somatic candidates con dati test = atteso (bassa profondità)
  Con dati reali attesi decine-centinaia di candidati somatic

## Report - Aggiornamenti post-produzione
- Sezione QC: aggiornare % mapping dopo riallineamento con --outReadsUnmapped Fastx
- Sezione 4.8: rimuovere nota "test data padj<0.5"
- Sezione M16: rimuovere nota "test data"
- Sezione M17: rimuovere nota "test data --statoff"
- Sezione M19: aggiornare con mean_dEF realistico (atteso 0.05-0.15)

## M20b HUMAnN3 Functional Profiling
### Setup (già fatto, da replicare con dati reali)
- env: humann3 (Python 3.9 separato da mosaic)
- humann v3.9 installato via pip
- metaphlan v4.2.4 installato via pip
- Database MetaPhlAn: mpa_vJan25_CHOCOPhlAnSGB_202503 (~32GB)
- Database ChocoPhlAn: DEMO → produzione: full_chocophlan (~15GB)
  humann_databases --download chocophlan full ~/Desktop/references/humann3_db
- Database UniRef: DEMO_diamond → produzione: uniref90_diamond (~20GB)
  humann_databases --download uniref uniref90_diamond ~/Desktop/references/humann3_db

### Fix applicati (da mantenere)
- metaphlan wrapper in ~/miniconda3/envs/humann3/bin/metaphlan
  traduce --bowtie2out -> --mapout e aggiunge --db_dir automaticamente
- humann config.py: major=4, line=0, metaphlan_v4_db_version="vJan25"

### Fix da fare in produzione
- humann_join_tables produce colonne duplicate se si copiano i file
  Soluzione: usare solo i file copiati, non mescolare con gli originali
  O meglio: usare --search-subdirectories se disponibile nella versione

### Comando produzione per campione
conda activate humann3
humann --input unmapped_R1.fastq --output outdir/SAMPLE \
  --threads 20 \
  --nucleotide-database ~/Desktop/references/humann3_db/chocophlan \
  --protein-database ~/Desktop/references/humann3_db/uniref
