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

    # Count reads per locus per sample using bedtools multicov.
    # This is a sweep-line algorithm → orders of magnitude faster than pysam loops.
    import subprocess as _sp
    import tempfile as _tmp

    print(f"  Running bedtools multicov on {len(bam_files)} BAMs...", file=sys.stderr)

    # Build a minimal BED3+name+score+strand from the loaded loci (preserve locus_id)
    bed_tmp = _tmp.NamedTemporaryFile(mode="w", suffix=".bed", delete=False)
    bed_tmp_path = Path(bed_tmp.name)
    for chrom, loci in sorted(loci_by_chrom.items()):
        for start, end, locus_id, name in loci:
            strand = locus_id.split("_")[-1]
            bed_tmp.write(f"{chrom}\t{start}\t{end}\t{locus_id}\t0\t{strand}\n")
    bed_tmp.close()

    multicov_cmd = [
        "bedtools", "multicov",
        "-bams", *[str(b) for b in bam_files],
        "-bed", str(bed_tmp_path),
        "-q", str(min_mapq),
        "-D",   # count duplicate reads
        "-s",   # strand-specific (match strand)
    ]
    result = _sp.run(multicov_cmd, capture_output=True, text=True)
    bed_tmp_path.unlink(missing_ok=True)

    if result.returncode != 0:
        # Fallback: retry without strand specificity
        print("  bedtools multicov failed with -s, retrying without strand filter...",
              file=sys.stderr)
        multicov_cmd_ns = [c for c in multicov_cmd if c != "-s"]
        result = _sp.run(multicov_cmd_ns, capture_output=True, text=True)

    if result.returncode != 0:
        raise RuntimeError(f"bedtools multicov failed:\n{result.stderr}")

    # Parse output: chrom, start, end, locus_id, 0, strand, count1, count2, ...
    out_tsv = out_dir / "locus_counts.tsv"
    n_expressed = 0
    with open(out_tsv, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["locus_id", "te_name", "chrom", "start", "end",
                         "strand"] + sample_ids)
        for line in result.stdout.splitlines():
            parts = line.split("\t")
            if len(parts) < 6 + len(bam_files):
                continue
            chrom_out = parts[0]
            start_out = parts[1]
            end_out   = parts[2]
            locus_id  = parts[3]
            strand_out = parts[5]
            counts    = [int(x) for x in parts[6:6 + len(bam_files)]]
            if sum(counts) == 0:
                continue
            # Extract te_name from locus_id (format: name_chrom_start_end_strand)
            te_name = locus_id.rsplit("_", 4)[0]
            writer.writerow([locus_id, te_name, chrom_out, start_out, end_out,
                             strand_out] + counts)
            n_expressed += 1

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
