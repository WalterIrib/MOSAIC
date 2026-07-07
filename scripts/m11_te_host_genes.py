#!/usr/bin/env python3
"""
scripts/m11_te_host_genes.py

M11 — TE host genes identification.

Uses RepeatMasker coordinates (rmsk_hg38.bed) to identify genes
that host TE/satellite elements. This provides the correct GO universe
for TE host gene enrichment analysis (M10).

Design Document v0.7:
  - rmsk_hg38.bed: genomic coordinates from RepeatMasker (NOT DFAM consensus)
  - Intersection with gene bodies (GTF) via bedtools
  - Output: te_host_genes.tsv (gene_id, gene_name, gene_type, te_families hosted)
  - GO universe: genes hosting any expressed TE family

Input:
  --rmsk      rmsk.txt.gz (UCSC RepeatMasker table)
  --gtf       gencode annotation GTF
  --te-counts te_counts_clean.tsv (to filter for expressed TE families)
  --gene-type-map gene_type_map.tsv
  --outdir    output directory

Output:
  rmsk_hg38.bed              RepeatMasker BED6 (TE + satellite only)
  genes.bed                  Gene body BED from GTF
  te_host_genes.tsv          gene_id | gene_name | gene_type | n_te_families | te_families
  te_host_universe.txt       gene_ids hosting expressed TE (one per line)
"""
from __future__ import annotations

import argparse
import csv
import gzip
import subprocess
import sys
from collections import defaultdict
from pathlib import Path


# ─────────────────────────────────────────────────────────────────────────────
# UCSC rmsk.txt → BED6
# ─────────────────────────────────────────────────────────────────────────────

# rmsk.txt columns (0-based):
# 0:bin 1:swScore 2:milliDiv 3:milliDel 4:milliIns
# 5:genoName 6:genoStart 7:genoEnd 8:genoLeft
# 9:strand 10:repName 11:repClass 12:repFamily
# 13:repStart 14:repEnd 15:repLeft 16:id

# Classes to KEEP (TE and satellite only)
KEEP_CLASSES = {
    "LINE", "SINE", "LTR", "DNA", "RC",
    "Retroposon", "Satellite",
    "Unknown",  # some TE families classified as Unknown
}

# Classes to EXCLUDE
EXCLUDE_CLASSES = {
    "Simple_repeat", "Low_complexity", "rRNA", "scRNA", "snRNA",
    "srpRNA", "tRNA", "RNA", "Other",
}


def rmsk_to_bed(rmsk_gz: Path, out_bed: Path, verbose: bool = True) -> int:
    """
    Convert UCSC rmsk.txt.gz to BED6 keeping only TE and satellite entries.
    BED6 columns: chrom, start, end, name, score, strand
    name = repName|repClass|repFamily
    """
    out_bed.parent.mkdir(parents=True, exist_ok=True)

    n_total = n_kept = 0
    opener = gzip.open if str(rmsk_gz).endswith(".gz") else open

    print(f"Converting rmsk → BED6: {rmsk_gz.name}...", file=sys.stderr)

    with opener(rmsk_gz, "rt") as fh, open(out_bed, "w") as out:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 13:
                continue
            n_total += 1

            chrom      = parts[5]
            start      = parts[6]
            end        = parts[7]
            strand     = parts[9]
            rep_name   = parts[10]
            rep_class  = parts[11]
            rep_family = parts[12]

            # Skip non-standard chromosomes for speed
            if "_" in chrom or chrom == "chrM":
                continue

            # Filter by class
            if rep_class in EXCLUDE_CLASSES:
                continue
            if rep_class not in KEEP_CLASSES:
                continue

            n_kept += 1
            name = f"{rep_name}|{rep_class}|{rep_family}"
            out.write(f"{chrom}\t{start}\t{end}\t{name}\t0\t{strand}\n")

    print(f"  rmsk: {n_total:,} total → {n_kept:,} TE/satellite kept",
          file=sys.stderr)
    return n_kept


# ─────────────────────────────────────────────────────────────────────────────
# GTF → gene body BED
# ─────────────────────────────────────────────────────────────────────────────

def gtf_to_gene_bed(gtf: Path, out_bed: Path,
                    gene_type_map: dict | None = None) -> int:
    """
    Extract gene body coordinates from GTF.
    One line per gene: chrom, start, end, gene_id, gene_name|gene_type, strand
    """
    out_bed.parent.mkdir(parents=True, exist_ok=True)
    opener = gzip.open if str(gtf).endswith(".gz") else open

    print(f"Extracting gene bodies from GTF...", file=sys.stderr)
    genes: dict[str, dict] = {}

    with opener(gtf, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "gene":
                continue

            attrs    = _parse_gtf_attrs(parts[8])
            gene_id  = attrs.get("gene_id", "")
            if not gene_id:
                continue

            chrom      = parts[0]
            start      = int(parts[3]) - 1  # 0-based
            end        = int(parts[4])
            strand     = parts[6]
            gene_name  = attrs.get("gene_name", gene_id)
            gene_type  = attrs.get("gene_type", attrs.get("gene_biotype", ""))

            if "_" in chrom or chrom == "chrM":
                continue

            genes[gene_id] = {
                "chrom":     chrom,
                "start":     start,
                "end":       end,
                "strand":    strand,
                "gene_id":   gene_id,
                "gene_name": gene_name,
                "gene_type": gene_type,
            }

    n_genes = len(genes)
    print(f"  Gene bodies: {n_genes:,}", file=sys.stderr)

    with open(out_bed, "w") as fh:
        for g in sorted(genes.values(),
                        key=lambda x: (x["chrom"], x["start"])):
            name = f"{g['gene_id']}|{g['gene_name']}|{g['gene_type']}"
            fh.write(f"{g['chrom']}\t{g['start']}\t{g['end']}\t"
                     f"{name}\t0\t{g['strand']}\n")

    return n_genes


def _parse_gtf_attrs(attr_str: str) -> dict[str, str]:
    d: dict[str, str] = {}
    for part in attr_str.strip().rstrip(";").split(";"):
        part = part.strip()
        if not part:
            continue
        sp = part.split(None, 1)
        if len(sp) == 2:
            d[sp[0]] = sp[1].strip('"')
    return d


# ─────────────────────────────────────────────────────────────────────────────
# bedtools intersect
# ─────────────────────────────────────────────────────────────────────────────

def intersect_te_genes(rmsk_bed: Path, genes_bed: Path,
                       out_intersect: Path) -> int:
    """
    Run bedtools intersect: TE loci × gene bodies.
    -wa -wb: write both A (TE) and B (gene) entries.
    -f 0.25: TE must overlap at least 25% of its length with gene.
    """
    out_intersect.parent.mkdir(parents=True, exist_ok=True)

    cmd = [
        "bedtools", "intersect",
        "-a", str(rmsk_bed),
        "-b", str(genes_bed),
        "-wa", "-wb",
        "-f", "0.10",   # 10% overlap threshold
    ]

    print(f"Running bedtools intersect...", file=sys.stderr)
    print(f"CMD: {' '.join(cmd)}", file=sys.stderr)

    with open(out_intersect, "w") as out:
        result = subprocess.run(cmd, stdout=out, stderr=subprocess.PIPE)

    if result.returncode != 0:
        print(f"ERROR: bedtools intersect failed: {result.stderr.decode()}",
              file=sys.stderr)
        raise RuntimeError("bedtools intersect failed")

    # Count lines
    with open(out_intersect) as fh:
        n_lines = sum(1 for _ in fh)

    print(f"  Intersections: {n_lines:,}", file=sys.stderr)
    return n_lines


# ─────────────────────────────────────────────────────────────────────────────
# Parse intersections → TE host genes table
# ─────────────────────────────────────────────────────────────────────────────

def parse_intersections(
    intersect_bed: Path,
    expressed_te_families: set[str] | None,
    out_tsv: Path,
    out_universe: Path,
) -> dict:
    """
    Parse bedtools intersect output.
    Columns: [TE BED6 cols 0-5] [gene BED6 cols 6-11]
    Returns summary dict.
    """
    # gene_id → set of TE family names
    gene_te_map: dict[str, set] = defaultdict(set)
    # gene_id → gene metadata
    gene_meta: dict[str, dict] = {}

    with open(intersect_bed) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 12:
                continue

            # TE info (cols 0-5)
            te_name = parts[3]  # repName|repClass|repFamily
            rep_name = te_name.split("|")[0]

            # Gene info (cols 6-11)
            gene_name_field = parts[9]  # gene_id|gene_name|gene_type
            gparts = gene_name_field.split("|")
            if len(gparts) < 3:
                continue

            gene_id   = gparts[0]
            gene_name = gparts[1]
            gene_type = gparts[2]

            # Filter for expressed TE families only
            if expressed_te_families and rep_name not in expressed_te_families:
                continue

            gene_te_map[gene_id].add(rep_name)
            gene_meta[gene_id] = {
                "gene_id":   gene_id,
                "gene_name": gene_name,
                "gene_type": gene_type,
            }

    # Write TE host genes table
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    with open(out_tsv, "w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=["gene_id", "gene_name", "gene_type",
                        "n_te_families", "te_families"],
            delimiter="\t"
        )
        writer.writeheader()
        for gene_id, te_set in sorted(gene_te_map.items()):
            rec = gene_meta[gene_id].copy()
            rec["n_te_families"] = len(te_set)
            rec["te_families"]   = ",".join(sorted(te_set))
            writer.writerow(rec)

    # Write universe file (one gene_id per line)
    with open(out_universe, "w") as fh:
        for gene_id in sorted(gene_te_map.keys()):
            fh.write(gene_id + "\n")

    n_host_genes = len(gene_te_map)
    print(f"  TE host genes: {n_host_genes:,}", file=sys.stderr)
    print(f"  Written: {out_tsv.name}", file=sys.stderr)
    print(f"  Universe: {out_universe.name}", file=sys.stderr)

    return {
        "n_host_genes": n_host_genes,
        "gene_te_map":  dict(gene_te_map),
    }


# ─────────────────────────────────────────────────────────────────────────────
# Load expressed TE families
# ─────────────────────────────────────────────────────────────────────────────

def load_expressed_te_families(te_counts_tsv: Path,
                                min_total_count: int = 1) -> set[str]:
    """
    Load TE families with total count > min_total_count across all samples.
    """
    expressed = set()
    with open(te_counts_tsv) as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader)
        for row in reader:
            if not row:
                continue
            family_id = row[0]
            total = sum(int(x) for x in row[1:] if x.isdigit())
            if total >= min_total_count:
                expressed.add(family_id)

    print(f"  Expressed TE families: {len(expressed)}", file=sys.stderr)
    return expressed


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="MOSAIC M11 — TE host genes identification"
    )
    parser.add_argument("--rmsk",         required=True,
                        help="rmsk.txt.gz (UCSC RepeatMasker)")
    parser.add_argument("--gtf",          required=True,
                        help="Genome annotation GTF")
    parser.add_argument("--te-counts",    required=False, default=None,
                        help="te_counts_clean.tsv (filter expressed families)")
    parser.add_argument("--outdir",       required=True,
                        help="Output directory")
    parser.add_argument("--min-count",    type=int, default=1,
                        help="Min total TE count to consider expressed [default: 1]")
    args = parser.parse_args()

    out_dir = Path(args.outdir)
    out_dir.mkdir(parents=True, exist_ok=True)

    print("\nMOSAIC M11 — TE Host Genes", file=sys.stderr)
    print(f"  rmsk:    {args.rmsk}", file=sys.stderr)
    print(f"  GTF:     {args.gtf}", file=sys.stderr)
    print(f"  outdir:  {out_dir}", file=sys.stderr)
    print("", file=sys.stderr)

    # Step 1: rmsk → BED
    rmsk_bed = out_dir / "rmsk_hg38.bed"
    if rmsk_bed.exists() and rmsk_bed.stat().st_size > 0:
        print(f"rmsk BED already exists: {rmsk_bed.name}", file=sys.stderr)
    else:
        rmsk_to_bed(Path(args.rmsk), rmsk_bed)

    # Step 2: GTF → gene BED
    genes_bed = out_dir / "genes.bed"
    if genes_bed.exists() and genes_bed.stat().st_size > 0:
        print(f"Gene BED already exists: {genes_bed.name}", file=sys.stderr)
    else:
        gtf_to_gene_bed(Path(args.gtf), genes_bed)

    # Step 3: bedtools intersect
    intersect_bed = out_dir / "te_gene_intersect.bed"
    if intersect_bed.exists() and intersect_bed.stat().st_size > 0:
        print(f"Intersect BED already exists: {intersect_bed.name}",
              file=sys.stderr)
    else:
        intersect_te_genes(rmsk_bed, genes_bed, intersect_bed)

    # Step 4: Load expressed TE families
    expressed_te = None
    if args.te_counts and Path(args.te_counts).exists():
        print("\nFiltering for expressed TE families...", file=sys.stderr)
        expressed_te = load_expressed_te_families(
            Path(args.te_counts), args.min_count
        )

    # Step 5: Parse intersections
    print("\nBuilding TE host gene table...", file=sys.stderr)
    result = parse_intersections(
        intersect_bed,
        expressed_te,
        out_dir / "te_host_genes.tsv",
        out_dir / "te_host_universe.txt",
    )

    print(f"\nDone.", file=sys.stderr)
    print(f"  TE host genes:    {result['n_host_genes']:,}", file=sys.stderr)
    print(f"  Output: {out_dir}", file=sys.stderr)


if __name__ == "__main__":
    main()
