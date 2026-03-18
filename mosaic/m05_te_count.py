"""
mosaic/m05_te_count.py

M05 — TE and Satellite fragment counting via Bowtie2 SE (R1 + R2 separately).

Design Document v0.7:
  - Two separate Bowtie2 SE runs per sample: R1 and R2
  - te_fragment_count.py deduplicates by read name → 1 count per fragment
  - --very-sensitive-local alignment
  - Separate counting for TE families and satellite families
  - Output: te_counts.tsv, satellite_counts.tsv (one column per sample)
"""
from __future__ import annotations

import csv
import subprocess
from collections import defaultdict
from pathlib import Path

from .log import info, warn, error, success, section


# ─────────────────────────────────────────────────────────────────────────────
# Bowtie2 SE mapping
# ─────────────────────────────────────────────────────────────────────────────

def run_bowtie2_se(
    sample_id: str,
    fastq: Path,
    index: str,
    out_bam: Path,
    mate: str,       # "R1" or "R2"
    threads: int = 4,
) -> dict:
    """
    Run Bowtie2 SE on one FASTQ file, sort and index BAM.
    Uses --very-sensitive-local.
    Returns mapping stats dict.
    """
    out_bam.parent.mkdir(parents=True, exist_ok=True)

    cmd_bowtie2 = [
        "bowtie2",
        "--very-sensitive-local",
        "--no-unal",          # discard unmapped reads
        "-p", str(threads),
        "--seed", "42",
        "-x", index,
        "-U", str(fastq),
    ]

    info(f"[{sample_id}] Bowtie2 SE {mate}...")

    # Use shell pipe — more reliable than Python Popen chaining for large streams
    bt2_log = out_bam.parent / f"{out_bam.stem}.bt2.log"
    cmd_list = cmd_bowtie2 + ["2>" + str(bt2_log), "|",
                               "samtools", "sort", "-@", str(threads), "-o", str(out_bam)]
    shell_cmd = (
        " ".join(cmd_bowtie2) +
        f" 2>{bt2_log}" +
        f" | samtools sort -@ {threads} -o {out_bam}"
    )
    result = subprocess.run(shell_cmd, shell=True)

    if result.returncode != 0:
        error(f"[{sample_id}] Bowtie2 {mate} failed")
        if bt2_log.exists():
            error(bt2_log.read_text()[-500:])
        raise RuntimeError(f"Bowtie2 SE failed for {sample_id} {mate}")

    bt2_stderr_text = bt2_log.read_text() if bt2_log.exists() else ""

    subprocess.run(["samtools", "index", str(out_bam)], check=True)

    # Parse alignment rate
    rate = _parse_bowtie2_rate(bt2_stderr_text)
    info(f"[{sample_id}] {mate} alignment rate: {rate:.1f}%")

    return {"sample_id": sample_id, "mate": mate,
            "alignment_rate": rate, "bam": str(out_bam)}


def _parse_bowtie2_rate(stderr: str) -> float:
    import re
    m = re.search(r"([\d.]+)% overall alignment rate", stderr)
    return float(m.group(1)) if m else 0.0


# ─────────────────────────────────────────────────────────────────────────────
# Fragment deduplication counter
# ─────────────────────────────────────────────────────────────────────────────

def count_te_fragments(
    bam_r1: Path,
    bam_r2: Path,
    out_tsv: Path,
    sample_id: str,
) -> dict[str, int]:
    """
    Merge R1 and R2 BAMs and count fragments per TE family.
    Deduplication: for each read name (= fragment), count only once
    per reference (family). If R1 and R2 both map to same family → 1 count.

    Returns dict: {family_name: count}
    """
    # Read R1 hits: read_name -> set of references
    r1_hits: dict[str, set] = defaultdict(set)
    r2_hits: dict[str, set] = defaultdict(set)

    def _read_bam(bam: Path, hits: dict):
        result = subprocess.run(
            ["samtools", "view", str(bam)],
            capture_output=True, text=True
        )
        for line in result.stdout.splitlines():
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            read_name = parts[0]
            ref_name  = parts[2]
            # Strip family name from "FAMILY#CLASS [category=...]"
            family = ref_name.split("#")[0]
            hits[read_name].add(family)

    _read_bam(bam_r1, r1_hits)
    _read_bam(bam_r2, r2_hits)

    # Merge: union of all read names
    all_reads = set(r1_hits.keys()) | set(r2_hits.keys())

    # Count: for each read, count once per family
    family_counts: dict[str, int] = defaultdict(int)
    for read in all_reads:
        families = r1_hits.get(read, set()) | r2_hits.get(read, set())
        for fam in families:
            family_counts[fam] += 1

    # Write TSV
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    with open(out_tsv, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["family_id", sample_id])
        for fam, cnt in sorted(family_counts.items()):
            writer.writerow([fam, cnt])

    info(f"[{sample_id}] Fragment counts: {len(family_counts)} families with >0 counts")
    return dict(family_counts)


# ─────────────────────────────────────────────────────────────────────────────
# Merge per-sample counts into matrix
# ─────────────────────────────────────────────────────────────────────────────

def merge_count_matrices(
    per_sample_counts: list[tuple[str, dict[str, int]]],
    out_tsv: Path,
) -> None:
    """
    Merge per-sample count dicts into a single TSV matrix.
    Rows = families, Columns = samples.
    """
    # Collect all families
    all_families: set[str] = set()
    for _, counts in per_sample_counts:
        all_families.update(counts.keys())

    sample_ids = [sid for sid, _ in per_sample_counts]

    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    with open(out_tsv, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["family_id"] + sample_ids)
        for fam in sorted(all_families):
            row = [fam] + [per_sample_counts[i][1].get(fam, 0)
                           for i in range(len(per_sample_counts))]
            writer.writerow(row)

    success(f"Count matrix: {out_tsv.name} "
            f"({len(all_families)} families × {len(sample_ids)} samples)")


# ─────────────────────────────────────────────────────────────────────────────
# Full M05 pipeline
# ─────────────────────────────────────────────────────────────────────────────

def run_m05(
    trim_stats: list[dict],
    te_index: str,
    satellite_index: str,
    classification_tsv: str | Path,
    out_dir: str | Path,
    threads: int = 4,
) -> dict:
    """
    Full M05:
      For each sample:
        1. Bowtie2 SE R1 → TE pseudochromosome
        2. Bowtie2 SE R2 → TE pseudochromosome
        3. Fragment deduplication → te_fragments per sample
        4. Same for satellite
      Merge into count matrices.
    """
    section("M05 — TE/Satellite Fragment Counting (Bowtie2 SE)")
    out_dir = Path(out_dir)

    te_counts_all:  list[tuple[str, dict]] = []
    sat_counts_all: list[tuple[str, dict]] = []

    for s in trim_stats:
        sid = s["sample_id"]
        r1  = Path(s["r1_trimmed"])
        r2  = Path(s["r2_trimmed"])
        bam_dir = out_dir / "bam" / sid

        # ── TE mapping ──────────────────────────────────────────────────────
        te_bam_r1 = bam_dir / f"{sid}_TE_R1.bam"
        te_bam_r2 = bam_dir / f"{sid}_TE_R2.bam"

        run_bowtie2_se(sid, r1, te_index, te_bam_r1, "R1", threads)
        run_bowtie2_se(sid, r2, te_index, te_bam_r2, "R2", threads)

        te_counts = count_te_fragments(
            te_bam_r1, te_bam_r2,
            out_dir / "per_sample" / f"{sid}_te_counts.tsv",
            sid,
        )
        te_counts_all.append((sid, te_counts))

        # ── Satellite mapping ────────────────────────────────────────────────
        sat_bam_r1 = bam_dir / f"{sid}_satellite_R1.bam"
        sat_bam_r2 = bam_dir / f"{sid}_satellite_R2.bam"

        run_bowtie2_se(sid, r1, satellite_index, sat_bam_r1, "R1", threads)
        run_bowtie2_se(sid, r2, satellite_index, sat_bam_r2, "R2", threads)

        sat_counts = count_te_fragments(
            sat_bam_r1, sat_bam_r2,
            out_dir / "per_sample" / f"{sid}_satellite_counts.tsv",
            sid,
        )
        sat_counts_all.append((sid, sat_counts))

        success(f"[{sid}] M05 done — "
                f"TE families: {len(te_counts)}, "
                f"Satellite families: {len(sat_counts)}")

    # Merge matrices
    te_matrix_path  = out_dir / "te_counts.tsv"
    sat_matrix_path = out_dir / "satellite_counts.tsv"

    merge_count_matrices(te_counts_all,  te_matrix_path)
    merge_count_matrices(sat_counts_all, sat_matrix_path)

    # Summary
    section("M05 — Summary")
    total_te_fams  = len(set().union(*[c for _, c in te_counts_all]))
    total_sat_fams = len(set().union(*[c for _, c in sat_counts_all]))
    info(f"  TE families with counts:        {total_te_fams}")
    info(f"  Satellite families with counts: {total_sat_fams}")
    info(f"  TE matrix:        {te_matrix_path}")
    info(f"  Satellite matrix: {sat_matrix_path}")

    return {
        "te_matrix":        str(te_matrix_path),
        "satellite_matrix": str(sat_matrix_path),
        "te_counts":        te_counts_all,
        "sat_counts":       sat_counts_all,
    }
