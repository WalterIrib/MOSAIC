#!/usr/bin/env python3
"""
scripts/build_pseudochromosomes.py

Fase 3 — DFAM filtering, classification, pseudochromosome construction.

Input:  Dfam-RepeatMasker.lib  (FASTA, all organisms)
Output: consensus_TE.fa
        consensus_satellite.fa
        pseudochromosome_TE.fa
        pseudochromosome_satellite.fa
        classification.tsv
        pseudochromosome_params.yaml

Design Document v0.7:
  - Spacer: 200 N between consecutive consensus sequences
  - Order:  lexicographic by family name (deterministic)
  - Separate files for TE and satellite
  - classification.tsv: family_id | category | repeat_type | repeat_subtype
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import sys
from datetime import datetime
from pathlib import Path

import yaml


# ─────────────────────────────────────────────────────────────────────────────
# Clade sets — which clades are present in hg38
# ─────────────────────────────────────────────────────────────────────────────

# Clades whose TE families are annotated in hg38 (RepeatMasker)
# Ordered from most specific to most general
HUMAN_CLADES = {
    # Human-specific
    "Homo_sapiens", "Homo sapiens",
    # Great apes
    "Homininae", "Hominidae", "Hominoidea",
    # Broader primates
    "Catarrhini", "Simiiformes", "Haplorrhini", "Primates",
    # Broader placental mammals
    "Euarchontoglires", "Boreoeutheria", "Eutheria",
    # All mammals
    "Theria", "Mammalia",
    # Amniotes and broader vertebrates
    "Amniota", "Tetrapoda", "Sarcopterygii", "Euteleostomi",
    "Teleostomi", "Vertebrata", "Craniata", "Chordata",
    # Very broad (present in all metazoans including human)
    "Deuterostomia", "Bilateria", "Eumetazoa", "Metazoa",
    "Opisthokonta", "Eukaryota",
    # Generic vertebrate tags used in RepeatMasker lib
    "Vertebrata_vertebrates",
}

# ─────────────────────────────────────────────────────────────────────────────
# Classification rules
# ─────────────────────────────────────────────────────────────────────────────

# Map repeat type (from FASTA header #TYPE/SUBTYPE) to MOSAIC category
def classify_family(name: str, repeat_class: str) -> str | None:
    """
    Returns category string or None if family should be excluded.
    Categories: TE_LINE | TE_SINE | TE_LTR | TE_DNA | satellite | other
    """
    rc = repeat_class.upper()

    # Satellites — always separate
    if "SATELLITE" in rc or "SAT" in rc:
        return "satellite"
    if "SIMPLE_REPEAT" in rc or "SIMPLE" in rc:
        return None   # exclude simple repeats
    if "LOW_COMPLEXITY" in rc:
        return None   # exclude low complexity
    if "RRNA" in rc or "SCRNA" in rc or "SNRNA" in rc:
        return None   # exclude structural RNAs

    # TE categories
    if "LINE" in rc:
        return "TE_LINE"
    if "SINE" in rc:
        return "TE_SINE"
    if "LTR" in rc:
        return "TE_LTR"
    if rc.startswith("DNA"):
        return "TE_DNA"
    if "RETROTRANSPOSON" in rc or "RETROPOSON" in rc:
        return "TE_LTR"   # generic retrotransposon -> LTR bucket
    if "RC" in rc:        # rolling circle
        return "TE_DNA"
    if rc in ("UNKNOWN", "OTHER", "TRANSPOSABLE_ELEMENT"):
        return "TE_DNA"   # catch-all -> DNA bucket

    return None   # rRNA, tRNA, etc. — exclude


# ─────────────────────────────────────────────────────────────────────────────
# FASTA parser
# ─────────────────────────────────────────────────────────────────────────────

def parse_repeatmasker_fasta(fasta_path: Path, target_clades: set[str],
                             verbose: bool = True):
    """
    Parse Dfam-RepeatMasker.lib FASTA.
    Header format: >NAME#TYPE/SUBTYPE @CLADE [S:score,...]
                or >NAME#TYPE/SUBTYPE @CLADE1,CLADE2,...

    Yields dicts for families present in target_clades.
    """
    current_header = None
    current_seq: list[str] = []
    total = kept = 0

    def _emit(header, seq):
        nonlocal total, kept
        total += 1
        rec = _parse_header(header, seq)
        if rec is None:
            return None
        # Check clade
        if not _is_target_clade(rec["clades"], target_clades):
            return None
        # Classify
        cat = classify_family(rec["name"], rec["repeat_class"])
        if cat is None:
            return None
        rec["category"] = cat
        kept += 1
        return rec

    with open(fasta_path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if current_header is not None:
                    rec = _emit(current_header, "".join(current_seq))
                    if rec:
                        yield rec
                current_header = line[1:]
                current_seq = []
            else:
                current_seq.append(line)

    if current_header is not None:
        rec = _emit(current_header, "".join(current_seq))
        if rec:
            yield rec

    if verbose:
        print(f"  Parsed {total} families, kept {kept} for target clades",
              file=sys.stderr)


def _parse_header(header: str, seq: str) -> dict | None:
    """
    Parse RepeatMasker FASTA header.
    Returns dict or None on parse failure.
    """
    # Split on spaces: first token is NAME#CLASS, rest is metadata
    parts = header.split()
    if not parts:
        return None

    name_class = parts[0]
    if "#" in name_class:
        name, repeat_class = name_class.split("#", 1)
    else:
        name = name_class
        repeat_class = "Unknown"

    # Parse clade from @CLADE token
    clades: list[str] = []
    for p in parts[1:]:
        if p.startswith("@"):
            clade_str = p[1:].rstrip(",").strip()
            # May be comma-separated
            clades = [c.strip() for c in clade_str.split(",") if c.strip()]
            break

    if not seq:
        return None

    return {
        "name":          name,
        "repeat_class":  repeat_class,
        "clades":        clades,
        "sequence":      seq.upper().replace(" ", ""),
        "length":        len(seq.replace(" ", "").replace("\n", "")),
    }


def _is_target_clade(family_clades: list[str], target_clades: set[str]) -> bool:
    """Check if any of the family's clades is in our target set."""
    for c in family_clades:
        # Exact match
        if c in target_clades:
            return True
        # Partial match for multi-word clades
        c_clean = c.replace("_", " ")
        if c_clean in target_clades:
            return True
        # Check if any target is a substring (e.g. "Vertebrata" in "Vertebrata_vertebrates")
        for t in target_clades:
            if t in c or c in t:
                return True
    return False


# ─────────────────────────────────────────────────────────────────────────────
# Pseudochromosome builder
# ─────────────────────────────────────────────────────────────────────────────

SPACER_N = 200

def build_pseudochromosome(families: list[dict], out_fasta: Path,
                           spacer_n: int = SPACER_N) -> dict:
    """
    Build pseudochromosome FASTA from list of family dicts.
    Families are sorted lexicographically by name (deterministic).
    Spacer: spacer_n N characters between consecutive sequences.

    Returns stats dict.
    """
    spacer = "N" * spacer_n
    families_sorted = sorted(families, key=lambda r: r["name"].lower())

    out_fasta.parent.mkdir(parents=True, exist_ok=True)
    total_len = 0
    n_families = 0

    with open(out_fasta, "w") as fh:
        for i, rec in enumerate(families_sorted):
            # Write FASTA header
            fh.write(f">{rec['name']}#{rec['repeat_class']} "
                     f"[category={rec['category']}]\n")
            # Write sequence in 60-char lines
            seq = rec["sequence"]
            for j in range(0, len(seq), 60):
                fh.write(seq[j:j+60] + "\n")
            # Write spacer (except after last)
            if i < len(families_sorted) - 1:
                fh.write(spacer + "\n")
            total_len += len(seq) + (spacer_n if i < len(families_sorted) - 1 else 0)
            n_families += 1

    md5 = _md5(out_fasta)
    print(f"  {out_fasta.name}: {n_families} families, "
          f"{total_len:,} bp, md5={md5[:8]}...", file=sys.stderr)

    return {
        "n_families":  n_families,
        "total_length": total_len,
        "md5":         md5,
        "path":        str(out_fasta),
    }


def _md5(path: Path) -> str:
    h = hashlib.md5()
    with open(path, "rb") as fh:
        while True:
            b = fh.read(1 << 20)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


# ─────────────────────────────────────────────────────────────────────────────
# Classification TSV
# ─────────────────────────────────────────────────────────────────────────────

def write_classification(families: list[dict], out_tsv: Path) -> None:
    """Write classification.tsv: family_id | category | repeat_type | repeat_subtype"""
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    with open(out_tsv, "w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=["family_id", "category", "repeat_class",
                        "clades", "length"],
            delimiter="\t"
        )
        writer.writeheader()
        for rec in sorted(families, key=lambda r: r["name"].lower()):
            writer.writerow({
                "family_id":    rec["name"],
                "category":     rec["category"],
                "repeat_class": rec["repeat_class"],
                "clades":       ",".join(rec["clades"]),
                "length":       rec["length"],
            })
    print(f"  Classification: {out_tsv.name} ({len(families)} entries)",
          file=sys.stderr)


# ─────────────────────────────────────────────────────────────────────────────
# Params YAML
# ─────────────────────────────────────────────────────────────────────────────

def write_params(out_yaml: Path, te_stats: dict, sat_stats: dict,
                 dfam_source: str, organism: str) -> None:
    params = {
        "build_date":        datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        "dfam_source":       dfam_source,
        "organism":          organism,
        "spacer_N":          SPACER_N,
        "family_order":      "lexicographic_by_name",
        "TE": {
            "n_families":    te_stats["n_families"],
            "total_length":  te_stats["total_length"],
            "md5":           te_stats["md5"],
            "path":          te_stats["path"],
        },
        "satellite": {
            "n_families":    sat_stats["n_families"],
            "total_length":  sat_stats["total_length"],
            "md5":           sat_stats["md5"],
            "path":          sat_stats["path"],
        },
    }
    out_yaml.parent.mkdir(parents=True, exist_ok=True)
    with open(out_yaml, "w") as fh:
        yaml.dump(params, fh, default_flow_style=False, sort_keys=False)
    print(f"  Params: {out_yaml.name}", file=sys.stderr)


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Build DFAM pseudochromosomes for MOSAIC"
    )
    parser.add_argument("--dfam-lib", required=True,
                        help="Path to Dfam-RepeatMasker.lib")
    parser.add_argument("--outdir", required=True,
                        help="Output directory")
    parser.add_argument("--organism", default="human",
                        help="Organism (human|mouse|...) — affects clade filter")
    parser.add_argument("--spacer-n", type=int, default=SPACER_N,
                        help=f"N spacer length (default: {SPACER_N})")
    args = parser.parse_args()

    dfam_lib = Path(args.dfam_lib)
    out_dir  = Path(args.outdir)
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"\nDFAM Pseudochromosome Builder — MOSAIC v0.7", file=sys.stderr)
    print(f"  Input:    {dfam_lib}", file=sys.stderr)
    print(f"  Outdir:   {out_dir}", file=sys.stderr)
    print(f"  Organism: {args.organism}", file=sys.stderr)
    print(f"  Spacer:   {args.spacer_n} N", file=sys.stderr)
    print(f"  Order:    lexicographic by family name", file=sys.stderr)
    print("", file=sys.stderr)

    # Parse and filter families
    print("Parsing DFAM library...", file=sys.stderr)
    all_families = list(parse_repeatmasker_fasta(
        dfam_lib, HUMAN_CLADES, verbose=True
    ))

    te_families  = [f for f in all_families if f["category"].startswith("TE_")]
    sat_families = [f for f in all_families if f["category"] == "satellite"]

    print(f"\nFiltered families:", file=sys.stderr)
    print(f"  TE:        {len(te_families)}", file=sys.stderr)
    print(f"  Satellite: {len(sat_families)}", file=sys.stderr)

    # Category breakdown
    from collections import Counter
    cats = Counter(f["category"] for f in all_families)
    for cat, n in sorted(cats.items()):
        print(f"    {cat:20} {n}", file=sys.stderr)

    if len(te_families) < 100:
        print(f"\nWARNING: Only {len(te_families)} TE families found. "
              "Expected >500 for human.", file=sys.stderr)

    # Write consensus FASTAs
    print("\nWriting consensus FASTAs...", file=sys.stderr)
    consensus_te  = out_dir / "consensus_TE.fa"
    consensus_sat = out_dir / "consensus_satellite.fa"

    with open(consensus_te, "w") as fh:
        for rec in sorted(te_families, key=lambda r: r["name"].lower()):
            fh.write(f">{rec['name']}#{rec['repeat_class']} "
                     f"[category={rec['category']}]\n")
            seq = rec["sequence"]
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i+60] + "\n")

    with open(consensus_sat, "w") as fh:
        for rec in sorted(sat_families, key=lambda r: r["name"].lower()):
            fh.write(f">{rec['name']}#{rec['repeat_class']} "
                     f"[category={rec['category']}]\n")
            seq = rec["sequence"]
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i+60] + "\n")

    # Build pseudochromosomes
    print("\nBuilding pseudochromosomes...", file=sys.stderr)
    pseudo_te  = out_dir / "pseudochromosome_TE.fa"
    pseudo_sat = out_dir / "pseudochromosome_satellite.fa"

    te_stats  = build_pseudochromosome(te_families,  pseudo_te,  args.spacer_n)
    sat_stats = build_pseudochromosome(sat_families, pseudo_sat, args.spacer_n)

    # Classification TSV
    write_classification(all_families, out_dir / "classification.tsv")

    # Params YAML
    write_params(
        out_dir / "pseudochromosome_params.yaml",
        te_stats, sat_stats,
        dfam_source=str(dfam_lib),
        organism=args.organism,
    )

    print("\nDone.", file=sys.stderr)
    print(f"  TE pseudochromosome:        {pseudo_te}", file=sys.stderr)
    print(f"  Satellite pseudochromosome: {pseudo_sat}", file=sys.stderr)
    print(f"  Classification:             {out_dir / 'classification.tsv'}", file=sys.stderr)
    print(f"  Params:                     {out_dir / 'pseudochromosome_params.yaml'}", file=sys.stderr)


if __name__ == "__main__":
    main()
