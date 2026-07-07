#!/usr/bin/env python3
"""
m12_gene_hits.py — Gene hit ranking from m12 locus-level DESeq2 results.

For each significant DE TE locus, finds the overlapping protein-coding / lncRNA
gene body (bedtools intersect) and counts loci per gene (n_up / n_dn).

Input:
    m12/{comp}_te_de_sig.tsv   — significant DE loci from m12_deseq2.R
    integrated/gene_type_map.tsv — gene coordinates and types

Output:
    m12/{comp}_gene_hits.tsv   — columns: gene_id, gene_name, gene_type,
                                           n_loci, n_up, n_dn
                                 sorted by n_loci desc
"""

import re
import sys
import csv
import subprocess
import tempfile
import os
from pathlib import Path
from collections import defaultdict


def parse_locus_coords(locus_id: str):
    """Extract chrom, start, end from locus_id: {te_name}_{chrom}_{start}_{end}_{strand}"""
    m = re.search(r'_([A-Za-z0-9]+)_(\d+)_(\d+)_([+-])$', locus_id)
    if m:
        return m.group(1), int(m.group(2)), int(m.group(3))
    return None, None, None


def run(sig_tsv: Path, gene_map_tsv: Path, out_tsv: Path) -> int:
    # ── 1. Read significant loci ───────────────────────────────────────────────
    loci = []
    with open(sig_tsv) as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        for row in reader:
            chrom, start, end = parse_locus_coords(row['locus_id'])
            if chrom is None:
                continue
            lfc = float(row.get('log2FoldChange_shrunk') or row.get('log2FoldChange') or 0)
            direction = 'up' if lfc >= 0 else 'dn'
            loci.append((chrom, start, end, row['locus_id'], direction))

    if not loci:
        print(f"[WARN] No loci parsed from {sig_tsv.name}", file=sys.stderr)
        out_tsv.parent.mkdir(parents=True, exist_ok=True)
        out_tsv.write_text("gene_id\tgene_name\tgene_type\tn_loci\tn_up\tn_dn\n")
        return 0

    print(f"[INF] Parsed {len(loci):,} significant loci", file=sys.stderr)

    # ── 2. Read gene bodies (protein_coding + lncRNA only) ────────────────────
    GENE_TYPES = {'protein_coding', 'lncRNA'}
    genes = []
    gene_meta = {}
    with open(gene_map_tsv) as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        for row in reader:
            if row['gene_type'] not in GENE_TYPES:
                continue
            try:
                start = int(row['start'])
                end   = int(row['end'])
            except (ValueError, KeyError):
                continue
            gid = row['gene_id']
            genes.append((row['chr'], start, end, gid))
            gene_meta[gid] = (row['gene_name'], row['gene_type'])

    print(f"[INF] {len(genes):,} protein_coding/lncRNA gene bodies loaded", file=sys.stderr)

    # ── 3. bedtools intersect ─────────────────────────────────────────────────
    with tempfile.TemporaryDirectory() as tmp:
        loci_bed   = Path(tmp) / "loci.bed"
        genes_bed  = Path(tmp) / "genes.bed"
        isect_out  = Path(tmp) / "intersect.tsv"

        # BED: chrom start end name lfc_direction
        with open(loci_bed, 'w') as fh:
            for chrom, start, end, lid, direction in sorted(loci):
                fh.write(f"{chrom}\t{start}\t{end}\t{lid}\t{direction}\n")

        with open(genes_bed, 'w') as fh:
            for chrom, start, end, gid in sorted(genes):
                fh.write(f"{chrom}\t{start}\t{end}\t{gid}\n")

        result = subprocess.run(
            ["bedtools", "intersect",
             "-a", str(loci_bed),
             "-b", str(genes_bed),
             "-wa", "-wb"],
            capture_output=True, text=True
        )
        if result.returncode != 0:
            print(f"[ERR] bedtools intersect failed: {result.stderr[:300]}", file=sys.stderr)
            return 1

        # ── 4. Count per gene ─────────────────────────────────────────────────
        gene_counts = defaultdict(lambda: {'n_loci': 0, 'n_up': 0, 'n_dn': 0})
        seen = set()  # avoid double-counting a locus in overlapping gene isoforms

        for line in result.stdout.splitlines():
            parts = line.split('\t')
            if len(parts) < 9:
                continue
            locus_id  = parts[3]
            direction = parts[4]
            gene_id   = parts[8]
            key = (locus_id, gene_id)
            if key in seen:
                continue
            seen.add(key)
            gene_counts[gene_id]['n_loci'] += 1
            gene_counts[gene_id][f'n_{direction}'] += 1

    # ── 5. Write output ───────────────────────────────────────────────────────
    rows = []
    for gid, counts in gene_counts.items():
        name, gtype = gene_meta.get(gid, (gid, 'unknown'))
        rows.append({
            'gene_id':   gid,
            'gene_name': name,
            'gene_type': gtype,
            'n_loci':    counts['n_loci'],
            'n_up':      counts['n_up'],
            'n_dn':      counts['n_dn'],
        })

    rows.sort(key=lambda r: r['n_loci'], reverse=True)

    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    with open(out_tsv, 'w', newline='') as fh:
        writer = csv.DictWriter(fh, fieldnames=['gene_id','gene_name','gene_type',
                                                 'n_loci','n_up','n_dn'],
                                delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)

    print(f"[INF] {len(rows):,} genes with ≥1 significant DE locus → {out_tsv.name}",
          file=sys.stderr)
    return 0


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(description="Gene hit ranking from m12 DESeq2 loci")
    ap.add_argument("--sig",      required=True, help="m12 *_te_de_sig.tsv")
    ap.add_argument("--gene-map", required=True, help="integrated/gene_type_map.tsv")
    ap.add_argument("--out",      required=True, help="output *_gene_hits.tsv")
    args = ap.parse_args()
    sys.exit(run(Path(args.sig), Path(getattr(args, 'gene_map')), Path(args.out)))
