#!/usr/bin/env python3
"""
scripts/m12_te_locus.py

M12 — TE locus-specific analysis.

Two components:
  1. rmsk_to_gtf(): converts rmsk_hg38.bed → TE.gtf for TEtranscripts
  2. LocusCounter: pysam-based locus-specific TE fragment counter

Design Document v0.7:
  - TEtranscripts handles multi-mapping via EM (family-level, complements M05)
  - Locus counter uses M05 BAMs + rmsk coordinates for locus-specific counts
  - Output: locus_counts.tsv (one row per TE locus, columns = samples)

Usage:
  # Step 1: build GTF (one-time)
  python3 scripts/m12_te_locus.py build-gtf \
      --rmsk-bed  mosaic_output/m11/rmsk_hg38.bed \
      --outdir    references/rmsk

  # Step 2: run TEtranscripts (called via subprocess)
  python3 scripts/m12_te_locus.py tetranscripts \
      --bam-dir   mosaic_output/mapping/bam \
      --coldata   samples.tsv \
      --gene-gtf  references/gencode.v49.annotation.gtf \
      --te-gtf    references/rmsk/TE.gtf \
      --outdir    mosaic_output/m12

  # Step 3: locus-specific counts via pysam
  python3 scripts/m12_te_locus.py locus-count \
      --bam-dir   mosaic_output/te_counts/bam \
      --rmsk-bed  mosaic_output/m11/rmsk_hg38.bed \
      --outdir    mosaic_output/m12 \
      --min-mapq  1
"""
from __future__ import annotations

import argparse
import csv
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

import pysam


# ─────────────────────────────────────────────────────────────────────────────
# 1. rmsk BED → TEtranscripts GTF
# ─────────────────────────────────────────────────────────────────────────────

def rmsk_to_gtf(rmsk_bed: Path, out_gtf: Path) -> int:
    """
    Convert rmsk_hg38.bed (BED6, name=repName|repClass|repFamily)
    to TEtranscripts-compatible GTF.

    TEtranscripts GTF format:
      chrom  RMSK  exon  start  end  .  strand  .
      gene_id "family_id"; transcript_id "locus_id"; family_id "family"; class_id "class";
    """
    out_gtf.parent.mkdir(parents=True, exist_ok=True)
    n = 0

    print(f"Converting rmsk BED → TEtranscripts GTF...", file=sys.stderr)

    with open(rmsk_bed) as fh, open(out_gtf, "w") as out:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6:
                continue

            chrom  = parts[0]
            start  = int(parts[1]) + 1  # BED 0-based → GTF 1-based
            end    = parts[2]
            name   = parts[3]           # repName|repClass|repFamily
            strand = parts[5]

            nparts     = name.split("|")
            rep_name   = nparts[0]
            rep_class  = nparts[1] if len(nparts) > 1 else "Unknown"
            rep_family = nparts[2] if len(nparts) > 2 else rep_class

            # Unique locus ID: repName_chrom_start_end_strand
            locus_id = f"{rep_name}_{chrom}_{parts[1]}_{end}_{strand}"
            # Family ID used as gene_id
            family_id = f"{rep_name}:{rep_family}"

            attr = (
                f'gene_id "{family_id}"; '
                f'transcript_id "{locus_id}"; '
                f'family_id "{rep_family}"; '
                f'class_id "{rep_class}";'
            )

            out.write(
                f"{chrom}\tRMSK\texon\t{start}\t{end}\t.\t{strand}\t.\t{attr}\n"
            )
            n += 1

    print(f"  GTF entries: {n:,}", file=sys.stderr)
    print(f"  Written: {out_gtf}", file=sys.stderr)
    return n


# ─────────────────────────────────────────────────────────────────────────────
# 2. TEtranscripts runner
# ─────────────────────────────────────────────────────────────────────────────

def run_tetranscripts(
    bam_dir: Path,
    coldata: Path,
    gene_gtf: Path,
    te_gtf: Path,
    out_dir: Path,
    threads: int = 4,
    stranded: str = "no",
) -> None:
    """
    Run TEtranscripts for DE analysis.
    Reads coldata (sample_id, condition) to split into treatment/control groups.
    """
    out_dir.mkdir(parents=True, exist_ok=True)

    # Parse coldata
    samples: dict[str, list] = defaultdict(list)
    with open(coldata) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            sid  = row.get("sample_id", "").strip()
            cond = row.get("condition", "").strip()
            if sid and cond:
                samples[cond].append(sid)

    conditions = list(samples.keys())
    if len(conditions) != 2:
        raise ValueError(
            f"TEtranscripts requires exactly 2 conditions, got: {conditions}"
        )

    # Identify treatment vs control (non-reference = treatment)
    # Convention: first condition alphabetically = control
    conditions_sorted = sorted(conditions)
    control   = conditions_sorted[0]
    treatment = conditions_sorted[1]

    print(f"  Control:   {control} ({len(samples[control])} samples)",
          file=sys.stderr)
    print(f"  Treatment: {treatment} ({len(samples[treatment])} samples)",
          file=sys.stderr)

    def bam_paths(sids):
        paths = []
        for sid in sids:
            bam = bam_dir / f"{sid}.bam"
            if not bam.exists():
                # Try subdirectory
                bam = bam_dir / sid / f"{sid}.bam"
            if bam.exists():
                paths.append(str(bam))
            else:
                print(f"[WRN] BAM not found for {sid}: {bam}", file=sys.stderr)
        return paths

    control_bams   = bam_paths(samples[control])
    treatment_bams = bam_paths(samples[treatment])

    if not control_bams or not treatment_bams:
        raise RuntimeError("No BAM files found for TEtranscripts")

    cmd = [
        "TEtranscripts",
        "--treatment", *treatment_bams,
        "--control",   *control_bams,
        "--GTF",       str(gene_gtf),
        "--TE",        str(te_gtf),
        "--mode",      "multi",           # EM multi-mapping
        "--stranded",  stranded,
        "--project",   str(out_dir / f"{treatment}_vs_{control}"),
        "--outdir",    str(out_dir),
        "--format",    "BAM",
        "--sortByPos",
    ]

    print(f"\nRunning TEtranscripts...", file=sys.stderr)
    print(f"CMD: {' '.join(cmd)}", file=sys.stderr)

    result = subprocess.run(cmd, stderr=subprocess.PIPE)
    if result.returncode != 0:
        print(f"[ERR] TEtranscripts failed:", file=sys.stderr)
        print(result.stderr.decode()[-1000:], file=sys.stderr)
        raise RuntimeError("TEtranscripts failed")

    print(f"[OK] TEtranscripts complete. Output in {out_dir}", file=sys.stderr)


# ─────────────────────────────────────────────────────────────────────────────
# 3. Locus-specific counter via pysam
# ─────────────────────────────────────────────────────────────────────────────

def locus_count(
    bam_dir: Path,
    rmsk_bed: Path,
    out_dir: Path,
    min_mapq: int = 1,
    max_loci: int | None = None,
) -> None:
    """
    Count reads mapping to each TE locus using pysam.
    Uses the M05 BAMs (already aligned to pseudochromosomes — NOT genome BAMs).
    
    NOTE: For true locus-specific counting we need genome-aligned BAMs.
    This function uses genome BAMs from M04 and intersects with rmsk loci.
    """
    out_dir.mkdir(parents=True, exist_ok=True)

    # Find all genome BAMs
    bam_files = sorted(bam_dir.glob("*.bam"))
    if not bam_files:
        # Try subdirectories
        bam_files = sorted(bam_dir.glob("*/*.bam"))
        bam_files = [b for b in bam_files if not b.name.endswith(".bai")]

    if not bam_files:
        raise RuntimeError(f"No BAM files found in {bam_dir}")

    sample_ids = [b.stem for b in bam_files]
    print(f"  Samples: {sample_ids}", file=sys.stderr)

    # Load rmsk loci into memory (grouped by chrom for efficiency)
    print(f"  Loading rmsk loci...", file=sys.stderr)
    loci_by_chrom: dict[str, list] = defaultdict(list)
    n_loci = 0

    import gzip as _gzip
    _opener = _gzip.open if str(rmsk_bed).endswith(".gz") else open
    with _opener(rmsk_bed, "rt") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6:
                continue
            chrom  = parts[0]
            start  = int(parts[1])
            end    = int(parts[2])
            name   = parts[3]
            strand = parts[5]
            locus_id = f"{name}_{chrom}_{start}_{end}_{strand}"
            loci_by_chrom[chrom].append((start, end, locus_id, name))
            n_loci += 1
            if max_loci and n_loci >= max_loci:
                break

    print(f"  Loci loaded: {n_loci:,}", file=sys.stderr)

    # Count reads per locus per sample
    # For efficiency: iterate BAM pileup per chrom, intersect with loci
    locus_counts: dict[str, dict[str, int]] = defaultdict(
        lambda: {sid: 0 for sid in sample_ids}
    )

    for bam_path, sample_id in zip(bam_files, sample_ids):
        print(f"  Counting {sample_id}...", file=sys.stderr)
        bam = pysam.AlignmentFile(str(bam_path), "rb")

        for chrom, loci in loci_by_chrom.items():
            try:
                for read in bam.fetch(chrom):
                    if read.is_unmapped:
                        continue
                    if read.mapping_quality < min_mapq:
                        continue
                    if read.is_duplicate:
                        continue

                    rstart = read.reference_start
                    rend   = read.reference_end or rstart + 1

                    # Check overlap with each locus on this chrom
                    for lstart, lend, locus_id, _ in loci:
                        if rstart < lend and rend > lstart:
                            locus_counts[locus_id][sample_id] += 1
                            break  # count each read once (first overlapping locus)
            except (ValueError, KeyError):
                pass  # chrom not in BAM

        bam.close()

    # Write output
    out_tsv = out_dir / "locus_counts.tsv"
    with open(out_tsv, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["locus_id", "te_name", "chrom", "start", "end",
                         "strand"] + sample_ids)

        for chrom, loci in sorted(loci_by_chrom.items()):
            for start, end, locus_id, name in loci:
                nparts = name.split("|")
                te_name = nparts[0]
                strand  = locus_id.split("_")[-1]
                counts  = [locus_counts[locus_id][sid] for sid in sample_ids]
                if sum(counts) > 0:  # only write expressed loci
                    writer.writerow([locus_id, te_name, chrom, start, end,
                                     strand] + counts)

    n_expressed = sum(1 for c in locus_counts.values() if sum(c.values()) > 0)
    print(f"  Expressed loci: {n_expressed:,}", file=sys.stderr)
    print(f"  Written: {out_tsv}", file=sys.stderr)


# ─────────────────────────────────────────────────────────────────────────────
# CLI
# ─────────────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="MOSAIC M12 — TE locus-specific analysis"
    )
    sub = parser.add_subparsers(dest="cmd", required=True)

    # build-gtf
    p1 = sub.add_parser("build-gtf", help="Convert rmsk BED to TEtranscripts GTF")
    p1.add_argument("--rmsk-bed", required=True)
    p1.add_argument("--outdir",   required=True)

    # tetranscripts
    p2 = sub.add_parser("tetranscripts", help="Run TEtranscripts DE")
    p2.add_argument("--bam-dir",   required=True)
    p2.add_argument("--coldata",   required=True)
    p2.add_argument("--gene-gtf",  required=True)
    p2.add_argument("--te-gtf",    required=True)
    p2.add_argument("--outdir",    required=True)
    p2.add_argument("--threads",   type=int, default=4)
    p2.add_argument("--stranded",  default="no")

    # locus-count
    p3 = sub.add_parser("locus-count", help="Locus-specific TE counts via pysam")
    p3.add_argument("--bam-dir",   required=True)
    p3.add_argument("--rmsk-bed",  required=True)
    p3.add_argument("--outdir",    required=True)
    p3.add_argument("--min-mapq",  type=int, default=1)
    p3.add_argument("--max-loci",  type=int, default=None,
                    help="Limit loci for testing")

    args = parser.parse_args()

    print(f"\nMOSAIC M12 — TE locus analysis [{args.cmd}]", file=sys.stderr)

    if args.cmd == "build-gtf":
        rmsk_to_gtf(Path(args.rmsk_bed),
                    Path(args.outdir) / "TE.gtf")

    elif args.cmd == "tetranscripts":
        run_tetranscripts(
            bam_dir   = Path(args.bam_dir),
            coldata   = Path(args.coldata),
            gene_gtf  = Path(args.gene_gtf),
            te_gtf    = Path(args.te_gtf),
            out_dir   = Path(args.outdir),
            threads   = args.threads,
            stranded  = args.stranded,
        )

    elif args.cmd == "locus-count":
        locus_count(
            bam_dir   = Path(args.bam_dir),
            rmsk_bed  = Path(args.rmsk_bed),
            out_dir   = Path(args.outdir),
            min_mapq  = args.min_mapq,
            max_loci  = args.max_loci,
        )


if __name__ == "__main__":
    main()
