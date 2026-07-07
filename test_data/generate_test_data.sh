#!/usr/bin/env bash
# =============================================================================
# generate_test_data.sh — build the MOSAIC functional-smoke-test subset.
#
# Produces a tiny, self-contained 2-vs-2 dataset (mito-ChaRNA vs ChaRNA) whose
# reads are restricted to chr1 + chr22, so the full pipeline (m02–m23 + M14)
# runs end-to-end in a few minutes. chr1 gives dense gene/TE (SINE/LINE)
# content; chr22 adds a pericentromeric satellite region.
#
# What it builds under test_data/:
#   reference/   mini genome (chr1+chr22), subset GTF / rmsk BED / TE.gtf,
#                mini STAR index; DFAM pseudochromosomes are symlinked (family-
#                level M05/M09 are chromosome-independent).
#   fastq/       4 small paired FASTQs: reads from the raw ChaRNA libraries that
#                map to chr1/chr22 (deterministic: first N pairs, fixed).
#   coldata.tsv  sample sheet
#   config_test.yaml   ready-to-run config pointing at the mini reference
#
# Inputs (override via env):
#   MOSAIC_REF        reference root   (default /home/walter/Desktop/references)
#   CHARNA_FASTQ_DIR  raw ChaRNA FASTQ (default USB METACAR_fastq_raw)
#   N_PAIRS           input read pairs to scan per sample (default 4000000)
#   THREADS           (default 8)
#
# Requires the "mosaic" conda env (STAR, samtools) on PATH.
# =============================================================================
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MOSAIC_REF="${MOSAIC_REF:-/home/walter/Desktop/references}"
CHARNA_FASTQ_DIR="${CHARNA_FASTQ_DIR:-/media/walter/EXTERNAL_USB/METACAR_fastq_raw}"
N_PAIRS="${N_PAIRS:-4000000}"
THREADS="${THREADS:-8}"
# Chromosomes for the subset. Default chr1+chr22 (dense genes + a satellite
# region). For a faster CI run set CHRS="chr22" (~10x fewer TE loci, so the
# m12 bedtools-multicov step finishes in ~1-2 min instead of ~25-40 min).
CHRS="${CHRS:-chr1 chr22}"

REF="$HERE/reference"
FQ="$HERE/fastq"
WORK="$HERE/.work"
mkdir -p "$REF" "$FQ" "$WORK"

GENOME="$MOSAIC_REF/GRCh38.primary_assembly.genome.fa"
GTF="$MOSAIC_REF/gencode.v49.annotation.gtf"
TE_GTF="$MOSAIC_REF/rmsk/TE.gtf"
RMSK="$MOSAIC_REF/rmsk/rmsk.txt.gz"          # optional
DFAM_DIR="$MOSAIC_REF/dfam/pseudochromosomes"

echo "=== [1/5] mini genome ($CHRS) ==="
if [[ ! -s "$REF/genome.fa" ]]; then
  [[ -s "$GENOME.fai" ]] || samtools faidx "$GENOME"
  samtools faidx "$GENOME" $CHRS > "$REF/genome.fa"
  samtools faidx "$REF/genome.fa"
fi

echo "=== [2/5] subset annotations to $CHRS ==="
# GTF (keep header + chr1/chr22 records)
if [[ ! -s "$REF/annotation.gtf" ]]; then
  awk 'BEGIN{OFS="\t"} /^#/{print; next} $1=="chr1"||$1=="chr22"{print}' "$GTF" > "$REF/annotation.gtf"
fi
# TE.gtf
if [[ -s "$TE_GTF" && ! -s "$REF/TE.gtf" ]]; then
  awk 'BEGIN{OFS="\t"} /^#/{print; next} $1=="chr1"||$1=="chr22"{print}' "$TE_GTF" > "$REF/TE.gtf"
fi
# rmsk BED (gzip table -> keep chr1/chr22). UCSC rmsk.txt.gz columns: genoName=col6.
if [[ -s "$RMSK" && ! -s "$REF/rmsk.txt.gz" ]]; then
  zcat "$RMSK" | awk 'BEGIN{OFS="\t"} $6=="chr1"||$6=="chr22"{print}' | gzip > "$REF/rmsk.txt.gz"
fi
# DFAM pseudochromosomes: family-level, chromosome-independent -> reuse as-is.
if [[ -d "$DFAM_DIR" && ! -e "$REF/dfam" ]]; then
  ln -s "$DFAM_DIR" "$REF/dfam"
fi

echo "=== [3/5] mini STAR index ==="
if [[ ! -s "$REF/star_index/SAindex" ]]; then
  mkdir -p "$REF/star_index"
  # genomeSAindexNbases must be reduced for a small genome: min(14, log2(L)/2 - 1)
  GLEN=$(awk '/^[^>]/{l+=length($0)} END{print l}' "$REF/genome.fa")
  SAIDX=$(python3 -c "import math;print(min(14,int(math.log2($GLEN)/2-1)))")
  STAR --runMode genomeGenerate --runThreadN "$THREADS" \
       --genomeDir "$REF/star_index" \
       --genomeFastaFiles "$REF/genome.fa" \
       --sjdbGTFfile "$REF/annotation.gtf" \
       --sjdbOverhang 149 \
       --genomeSAindexNbases "$SAIDX" > "$WORK/star_index.log" 2>&1
fi

echo "=== [4/5] extract chr1/chr22 reads into small FASTQs (N_PAIRS=$N_PAIRS) ==="
# Sample sheet: control=ChaRNA{1,2}, treatment=mito-ChaRNA{1,2}
declare -A R1 R2
resolve() { ls "$CHARNA_FASTQ_DIR/$1"_R"$2".fastq*.gz 2>/dev/null | head -1; }
for s in ChaRNA1 ChaRNA2 mito-ChaRNA1 mito-ChaRNA2; do
  R1[$s]=$(resolve "$s" 1); R2[$s]=$(resolve "$s" 2)
  [[ -s "${R1[$s]}" && -s "${R2[$s]}" ]] || { echo "FASTQ missing for $s"; exit 1; }
done

extract() {
  local s="$1"
  if [[ -s "$FQ/${s}_R1.fastq.gz" ]]; then echo "  $s: already done"; return; fi
  echo "  $s: head -${N_PAIRS} + align + keep mapped"
  local nlines=$(( N_PAIRS * 4 ))
  # head closes the pipe early -> zcat gets SIGPIPE (141); isolate from pipefail.
  set +o pipefail
  zcat "${R1[$s]}" | head -n "$nlines" | gzip > "$WORK/${s}_sub_R1.fq.gz"
  zcat "${R2[$s]}" | head -n "$nlines" | gzip > "$WORK/${s}_sub_R2.fq.gz"
  set -o pipefail
  STAR --runThreadN "$THREADS" --genomeDir "$REF/star_index" \
       --readFilesIn "$WORK/${s}_sub_R1.fq.gz" "$WORK/${s}_sub_R2.fq.gz" \
       --readFilesCommand zcat \
       --outSAMtype BAM Unsorted \
       --outSAMunmapped None \
       --outFileNamePrefix "$WORK/${s}_" > "$WORK/${s}_star.log" 2>&1
  # keep only properly-paired mapped reads, write back to paired FASTQ
  samtools collate -@ "$THREADS" -O -u "$WORK/${s}_Aligned.out.bam" "$WORK/${s}_col" \
    | samtools fastq -@ "$THREADS" -f 0x2 -F 0x900 \
        -1 "$FQ/${s}_R1.fastq.gz" -2 "$FQ/${s}_R2.fastq.gz" \
        -0 /dev/null -s /dev/null -n - 2> "$WORK/${s}_fastq.log"
  local n=$(( $(zcat "$FQ/${s}_R1.fastq.gz" | wc -l) / 4 ))
  echo "    -> $n read pairs on $CHRS"
}
for s in ChaRNA1 ChaRNA2 mito-ChaRNA1 mito-ChaRNA2; do extract "$s"; done

echo "=== [5/5] coldata + test config ==="
# coldata carries explicit R1/R2 paths (flat FASTQ dir layout).
{
  printf 'sample_id\tcondition\tfraction\tphase\tR1\tR2\n'
  for s in ChaRNA1 ChaRNA2 mito-ChaRNA1 mito-ChaRNA2; do
    case "$s" in ChaRNA*) cond=ChaRNA; ph=interphase;; esac
    case "$s" in mito-*) cond=mito-ChaRNA; ph=mitosis;; esac
    printf '%s\t%s\tchromatin\t%s\t%s\t%s\n' \
      "$s" "$cond" "$ph" "$FQ/${s}_R1.fastq.gz" "$FQ/${s}_R2.fastq.gz"
  done
} > "$HERE/coldata.tsv"

cat > "$HERE/config_test.yaml" <<YAML
# Auto-generated by generate_test_data.sh — MOSAIC functional smoke test.
project:
  name: mosaic_smoke
  wd: ${HERE}          # pipeline appends /mosaic_output
  organism: human
  genome_build: hg38
  annotation_version: gencode_v49
  random_seed: 42
sequencing:
  type: PE
  read_length: 150
  strandedness: auto
library:
  prep: rRNA_depletion
  layout: PE
samples:
  fastq_dir: ${FQ}
  coldata:   ${HERE}/coldata.tsv
comparisons:
  - treatment: mito-ChaRNA
    control:   ChaRNA
databases:
  db_root:      ${REF}
  genome_fa:    ${REF}/genome.fa
  gtf:          ${REF}/annotation.gtf
  star_index:   ${REF}/star_index
  te_consensus: ${REF}/dfam/consensus_TE.fa
  te_satellite: ${REF}/dfam/consensus_satellite.fa
  rmsk_bed:     ${REF}/rmsk.txt.gz
  te_gtf:       ${REF}/TE.gtf
  auto_download: false
environments:
  main:   mosaic
  humann: humann3
  arriba: null
thresholds:
  padj: 0.05
  lfc: 1.0
  min_cov_editing: 5
  min_score_editing: 1.5
  min_ef_editing: 0.05
  editing_model: "DirMult:PHRED"
  min_reads_splicing: 10
  min_clip_te: 30
  min_support_te: 2
  min_reads_meta: 10
  min_samples_meta: 2
continue_on_error: true
modules:
  m00_download: false
  m02_trim: true
  m04_star: true
  m05_te_count: true
  m06_circrna: false        # CIRI2 external Perl script; off in smoke test
  m07_strandedness: true
  m08_deseq2: true
  m09_satellite: true
  m10_go: true
  m11_te_host: true
  m11b_promoter_te: true
  m11c_te_age: true
  m11d_te_de: true
  m12_tetranscripts: true
  m15_fusion: false         # Arriba external binary; off in smoke test
  m16_stringtie: true
  m17_rmats: true
  m18_network: true
  m19_editing: true
  m20_meta: false
  m20b_humann: false
  m23_te_insertions: true
  m14_report: true
threads: ${THREADS}
memory_gb: 12
trimming:
  min_phred: 20
  min_length: 36
  low_complexity_filter: false
  trim_poly_x: false
report:
  title: "MOSAIC smoke test"
  remove_test_notes: false
  produce_pdf: false
YAML

echo ""
echo "DONE. Test config: $HERE/config_test.yaml"
echo "Run the smoke test:  pytest tests/ -v   (or python tests/smoke_test.py)"
