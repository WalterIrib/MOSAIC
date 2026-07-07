"""
mosaic/build_integrated_matrix.py

Fase 4 — Build integrated count matrix.

Merges:
  1. featureCounts output (canonical genes + lncRNA)
  2. TE fragment counts (te_counts.tsv)
  3. Satellite fragment counts (satellite_counts.tsv)

Output:
  all_genes_counts_clean.tsv   — canonical genes only (no TE/satellite)
  lncrna_counts.tsv            — lncRNA subset (for GO M10)
  te_counts_clean.tsv          — TE counts aligned to sample order
  satellite_counts_clean.tsv   — satellite counts aligned to sample order
  integrated_matrix.tsv        — ALL features merged (genes + TE + satellite)
  gene_type_map.tsv            — gene_id | gene_name | gene_type | chr | start | end | strand

Design Document v0.7:
  - all_genes_counts.tsv contains ALL genes including lncRNA (no double counting)
  - gene_type_map.tsv assigns category from GTF gene_type field
  - lncrna_counts.tsv used ONLY for GO (M10)
  - TE and satellite are SEPARATE from canonical gene counts
"""
from __future__ import annotations

import csv
import gzip
from pathlib import Path
from typing import Optional

from .log import info, warn, error, success, section


# ─────────────────────────────────────────────────────────────────────────────
# Parse featureCounts output
# ─────────────────────────────────────────────────────────────────────────────

def parse_featurecounts(fc_tsv: Path) -> tuple[list[str], dict[str, dict]]:
    """
    Parse featureCounts output TSV.
    Returns:
      sample_ids: list of sample names (stripped of path and .bam)
      genes: dict {gene_id: {gene_name, gene_type, chr, start, end, strand, length,
                              counts: {sample_id: count}}}
    """
    genes: dict[str, dict] = {}
    sample_ids: list[str] = []

    with open(fc_tsv) as fh:
        for line in fh:
            line = line.rstrip("\n")
            # Skip comment lines
            if line.startswith("#"):
                continue
            parts = line.split("\t")

            # Header line
            if parts[0] == "Geneid":
                # Columns: Geneid Chr Start End Strand Length [extraAttrs...] sample1 sample2...
                # Find where BAM paths start (contain .bam)
                bam_cols = [i for i, p in enumerate(parts) if p.endswith(".bam")]
                extra_cols = [i for i, p in enumerate(parts)
                              if i > 5 and not p.endswith(".bam")]

                # Extract sample IDs from BAM paths
                sample_ids = [Path(parts[i]).stem for i in bam_cols]
                bam_indices = bam_cols
                extra_indices = extra_cols
                header_parts = parts
                continue

            if not sample_ids:
                continue

            gene_id = parts[0]
            chrom   = parts[1].split(";")[0]  # take first chr for multi-exon
            start   = parts[3].split(";")[0]
            end     = parts[4].split(";")[0]
            strand  = parts[5].split(";")[0] if len(parts) > 5 else "."
            length  = parts[6] if len(parts) > 6 else "0"

            # Extra attributes (gene_name, gene_type)
            gene_name = ""
            gene_type = ""
            if len(extra_indices) >= 1 and len(parts) > extra_indices[0]:
                # First extra col after length — depends on --extraAttributes order
                # We specified: gene_name,gene_type
                idx = extra_indices[0]
                if idx < len(parts):
                    gene_name = parts[idx]
                if idx + 1 < len(parts) and len(extra_indices) >= 2:
                    gene_type = parts[extra_indices[1]] if extra_indices[1] < len(parts) else ""

            # Count values
            counts = {}
            for i, sid in zip(bam_indices, sample_ids):
                counts[sid] = int(parts[i]) if i < len(parts) else 0

            genes[gene_id] = {
                "gene_id":   gene_id,
                "gene_name": gene_name,
                "gene_type": gene_type,
                "chr":       chrom,
                "start":     start,
                "end":       end,
                "strand":    strand,
                "length":    length,
                "counts":    counts,
            }

    info(f"Parsed featureCounts: {len(genes)} genes, {len(sample_ids)} samples")
    return sample_ids, genes


# ─────────────────────────────────────────────────────────────────────────────
# Parse TE / satellite count matrices
# ─────────────────────────────────────────────────────────────────────────────

def parse_count_matrix(tsv: Path) -> tuple[list[str], dict[str, dict[str, int]]]:
    """
    Parse te_counts.tsv or satellite_counts.tsv.
    Returns: sample_ids, {family_id: {sample_id: count}}
    """
    sample_ids: list[str] = []
    counts: dict[str, dict[str, int]] = {}

    with open(tsv) as fh:
        reader = csv.reader(fh, delimiter="\t")
        for i, row in enumerate(reader):
            if i == 0:
                sample_ids = row[1:]
                continue
            family_id = row[0]
            counts[family_id] = {
                sid: int(row[j + 1]) for j, sid in enumerate(sample_ids)
            }

    info(f"Parsed {tsv.name}: {len(counts)} features, {len(sample_ids)} samples")
    return sample_ids, counts


# ─────────────────────────────────────────────────────────────────────────────
# Write clean count matrices
# ─────────────────────────────────────────────────────────────────────────────

def write_count_matrix(
    features: dict[str, dict],
    sample_ids: list[str],
    out_tsv: Path,
    id_field: str = "gene_id",
    extra_fields: list[str] = None,
) -> None:
    """Write a count matrix TSV with optional extra annotation columns."""
    extra_fields = extra_fields or []
    out_tsv.parent.mkdir(parents=True, exist_ok=True)

    with open(out_tsv, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow([id_field] + extra_fields + sample_ids)
        for fid, rec in features.items():
            extras = [rec.get(f, "") for f in extra_fields]
            counts_row = [rec["counts"].get(sid, 0) for sid in sample_ids]
            writer.writerow([fid] + extras + counts_row)

    info(f"Written: {out_tsv.name} "
         f"({len(features)} features × {len(sample_ids)} samples)")


def write_simple_matrix(
    counts: dict[str, dict[str, int]],
    sample_ids: list[str],
    out_tsv: Path,
) -> None:
    """Write a simple family_id × sample count matrix."""
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    with open(out_tsv, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["family_id"] + sample_ids)
        for fid, sample_counts in sorted(counts.items()):
            writer.writerow([fid] + [sample_counts.get(sid, 0)
                                     for sid in sample_ids])
    info(f"Written: {out_tsv.name} "
         f"({len(counts)} features × {len(sample_ids)} samples)")


# ─────────────────────────────────────────────────────────────────────────────
# Gene type map
# ─────────────────────────────────────────────────────────────────────────────

def write_gene_type_map(genes: dict[str, dict], out_tsv: Path) -> None:
    """
    Write gene_type_map.tsv: gene_id | gene_name | gene_type | chr | start | end | strand
    """
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    with open(out_tsv, "w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=["gene_id", "gene_name", "gene_type",
                        "chr", "start", "end", "strand"],
            delimiter="\t",
            extrasaction="ignore",
        )
        writer.writeheader()
        for gene_id, rec in sorted(genes.items()):
            writer.writerow(rec)
    info(f"Gene type map: {out_tsv.name} ({len(genes)} entries)")


# ─────────────────────────────────────────────────────────────────────────────
# Integrated matrix builder
# ─────────────────────────────────────────────────────────────────────────────

def build_integrated_matrix(
    fc_tsv: str | Path,
    te_tsv: str | Path,
    satellite_tsv: str | Path,
    out_dir: str | Path,
    sample_order: Optional[list[str]] = None,
) -> dict:
    """
    Build integrated count matrix merging canonical genes + TE + satellite.

    Returns dict with paths to all output files.
    """
    section("Building Integrated Count Matrix")
    out_dir  = Path(out_dir)
    fc_tsv   = Path(fc_tsv)
    te_tsv   = Path(te_tsv)
    sat_tsv  = Path(satellite_tsv)

    # ── Parse inputs ─────────────────────────────────────────────────────────
    fc_samples,  genes      = parse_featurecounts(fc_tsv)
    te_samples,  te_counts  = parse_count_matrix(te_tsv)
    sat_samples, sat_counts = parse_count_matrix(sat_tsv)

    # Reconcile sample order
    if sample_order:
        samples = sample_order
    else:
        # Use featureCounts order as canonical
        samples = fc_samples
        # Verify TE and satellite have same samples
        for sid in samples:
            if sid not in te_samples:
                warn(f"Sample '{sid}' in featureCounts but not in TE counts")

    info(f"Sample order: {samples}")

    # ── Canonical genes ───────────────────────────────────────────────────────
    # Re-align counts to sample order
    for gene_id, rec in genes.items():
        rec["counts"] = {sid: rec["counts"].get(sid, 0) for sid in samples}

    # Separate lncRNA
    lncrna_genes = {gid: rec for gid, rec in genes.items()
                    if "lncRNA" in rec.get("gene_type", "").lower() or
                       "lncrna" in rec.get("gene_type", "").lower()}

    protein_coding = {gid: rec for gid, rec in genes.items()
                      if rec.get("gene_type", "") == "protein_coding"}

    info(f"Gene types: {len(genes)} total, "
         f"{len(protein_coding)} protein_coding, "
         f"{len(lncrna_genes)} lncRNA")

    # ── Write canonical matrix (all genes including lncRNA) ───────────────────
    all_genes_clean = out_dir / "all_genes_counts_clean.tsv"
    write_count_matrix(
        genes, samples, all_genes_clean,
        id_field="gene_id",
        extra_fields=["gene_name", "gene_type"],
    )

    # ── Write lncRNA matrix (for GO M10 only) ─────────────────────────────────
    lncrna_path = out_dir / "lncrna_counts.tsv"
    write_count_matrix(
        lncrna_genes, samples, lncrna_path,
        id_field="gene_id",
        extra_fields=["gene_name", "gene_type"],
    )

    # ── Write gene type map ────────────────────────────────────────────────────
    gene_type_map_path = out_dir / "gene_type_map.tsv"
    write_gene_type_map(genes, gene_type_map_path)

    # ── Align TE/satellite to sample order ────────────────────────────────────
    # Exclude non-TE entries that appear in RepeatMasker but are RNA genes or
    # mapping artefacts: 7SLRNA (SRP RNA), chrM (mitochondrial mapping noise).
    NON_TE_EXCLUDE = {"7SLRNA"}
    te_aligned  = {fid: {sid: counts.get(sid, 0) for sid in samples}
                   for fid, counts in te_counts.items()
                   if fid not in NON_TE_EXCLUDE}
    sat_aligned = {fid: {sid: counts.get(sid, 0) for sid in samples}
                   for fid, counts in sat_counts.items()}

    excluded_non_te = set(te_counts.keys()) & NON_TE_EXCLUDE
    if excluded_non_te:
        info(f"Excluded non-TE entries from te_counts_clean: {sorted(excluded_non_te)}")

    te_clean_path  = out_dir / "te_counts_clean.tsv"
    sat_clean_path = out_dir / "satellite_counts_clean.tsv"
    write_simple_matrix(te_aligned,  samples, te_clean_path)
    write_simple_matrix(sat_aligned, samples, sat_clean_path)

    # ── Integrated matrix: genes + TE + satellite ─────────────────────────────
    integrated_path = out_dir / "integrated_matrix.tsv"
    _write_integrated_matrix(
        genes, te_aligned, sat_aligned, samples, integrated_path
    )

    # ── Summary ───────────────────────────────────────────────────────────────
    section("Integrated Matrix Summary")
    info(f"  Canonical genes:  {len(genes):>8,}")
    info(f"    protein_coding: {len(protein_coding):>8,}")
    info(f"    lncRNA:         {len(lncrna_genes):>8,}")
    info(f"  TE families:      {len(te_counts):>8,}")
    info(f"  Satellite:        {len(sat_counts):>8,}")
    info(f"  TOTAL features:   {len(genes)+len(te_counts)+len(sat_counts):>8,}")
    info(f"  Samples:          {len(samples):>8,}")
    info(f"")
    info(f"  Output directory: {out_dir}")

    return {
        "all_genes_clean":    str(all_genes_clean),
        "lncrna_counts":      str(lncrna_path),
        "gene_type_map":      str(gene_type_map_path),
        "te_counts_clean":    str(te_clean_path),
        "satellite_counts":   str(sat_clean_path),
        "integrated_matrix":  str(integrated_path),
        "samples":            samples,
        "n_genes":            len(genes),
        "n_te":               len(te_counts),
        "n_satellite":        len(sat_counts),
    }


def _write_integrated_matrix(
    genes: dict[str, dict],
    te_counts: dict[str, dict],
    sat_counts: dict[str, dict],
    samples: list[str],
    out_tsv: Path,
) -> None:
    """Write integrated matrix with feature_type column."""
    out_tsv.parent.mkdir(parents=True, exist_ok=True)

    with open(out_tsv, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["feature_id", "feature_type",
                         "gene_name", "gene_type"] + samples)

        # Canonical genes
        for gene_id, rec in sorted(genes.items()):
            gt = rec.get("gene_type", "")
            feat_type = "lncRNA" if "lnc" in gt.lower() else gt or "gene"
            counts_row = [rec["counts"].get(sid, 0) for sid in samples]
            writer.writerow([gene_id, feat_type,
                             rec.get("gene_name", ""),
                             gt] + counts_row)

        # TE families
        for fid, sample_counts in sorted(te_counts.items()):
            counts_row = [sample_counts.get(sid, 0) for sid in samples]
            writer.writerow([fid, "TE", fid, "TE"] + counts_row)

        # Satellite families
        for fid, sample_counts in sorted(sat_counts.items()):
            counts_row = [sample_counts.get(sid, 0) for sid in samples]
            writer.writerow([fid, "satellite", fid, "satellite"] + counts_row)

    total = len(genes) + len(te_counts) + len(sat_counts)
    success(f"Integrated matrix: {out_tsv.name} "
            f"({total:,} features × {len(samples)} samples)")
