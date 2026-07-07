"""
mosaic/m02_trim.py
M02 — Quality Control e Trimming con fastp.
Parametri fissi da Design Document v0.7.
"""
from __future__ import annotations

import csv
import json
import subprocess
import sys
from pathlib import Path
from typing import Optional

from .log import info, warn, error, success, section


# ─────────────────────────────────────────────────────────────────────────────
# Sample manifest
# ─────────────────────────────────────────────────────────────────────────────

def load_sample_manifest(tsv_path: str | Path) -> list[dict]:
    """
    Load sample manifest TSV.
    Required columns: sample_id, condition, R1, R2
    Returns list of dicts.
    """
    p = Path(tsv_path)
    if not p.exists():
        raise FileNotFoundError(f"Sample manifest not found: {p}")

    samples = []
    with open(p) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            # Validate required columns
            for col in ("sample_id", "R1", "R2"):
                if col not in row or not row[col].strip():
                    raise ValueError(f"Missing column '{col}' in manifest row: {row}")
            # Validate files exist
            for fq in (row["R1"], row["R2"]):
                if not Path(fq.strip()).exists():
                    raise FileNotFoundError(f"FASTQ not found: {fq}")
            samples.append({k: v.strip() for k, v in row.items()})

    if not samples:
        raise ValueError("Sample manifest is empty.")

    info(f"Loaded {len(samples)} samples from manifest.")
    return samples


# ─────────────────────────────────────────────────────────────────────────────
# fastp runner
# ─────────────────────────────────────────────────────────────────────────────

def run_fastp(
    sample_id: str,
    r1_in: str | Path,
    r2_in: str | Path,
    out_dir: str | Path,
    threads: int = 4,
    min_phred: int = 20,
    min_length: int = 36,
    min_length_warning: int = 75,
    pe_correction: bool = True,
    cut_right: bool = True,
    cut_window_size: int = 4,
    cut_mean_quality: int = 20,
    disable_adapter_trimming: bool = False,
) -> dict:
    """
    Run fastp on one PE sample with MOSAIC-fixed parameters.

    CRITICAL parameters (non-negotiable per design doc v0.7):
      --low_complexity_filter: always OFF  (TE/satellite safety)
      --trim_poly_x:           always OFF  (satellite: TTAGGG, CATTC)

    Returns dict with trimming stats.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    r1_out  = out_dir / f"{sample_id}_R1_trimmed.fastq.gz"
    r2_out  = out_dir / f"{sample_id}_R2_trimmed.fastq.gz"
    json_out = out_dir / f"{sample_id}.fastp.json"
    html_out = out_dir / f"{sample_id}.fastp.html"

    cmd = [
        "fastp",
        "--in1",  str(r1_in),
        "--in2",  str(r2_in),
        "--out1", str(r1_out),
        "--out2", str(r2_out),
        "--json", str(json_out),
        "--html", str(html_out),
        # Quality filtering
        "--qualified_quality_phred", str(min_phred),
        "--unqualified_percent_limit", "40",
        "--length_required", str(min_length),
        # CRITICAL: always OFF
        # trim_poly_x is OFF by default in fastp 0.22 (no flag needed)
        # PE correction
        "--thread", str(threads),
        "--detect_adapter_for_pe",
    ]

    if pe_correction:
        cmd.append("--correction")

    if cut_right:
        cmd += [
            "--cut_right",
            "--cut_window_size", str(cut_window_size),
            "--cut_mean_quality", str(cut_mean_quality),
        ]

    if disable_adapter_trimming:
        cmd.append("--disable_adapter_trimming")

    info(f"[{sample_id}] Running fastp...")
    info("CMD: " + " ".join(cmd))

    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        error(f"[{sample_id}] fastp failed:")
        error(result.stderr[-2000:])
        raise RuntimeError(f"fastp failed for {sample_id}")

    # Parse JSON stats
    stats = _parse_fastp_json(json_out, sample_id, min_length_warning)
    stats["r1_trimmed"] = str(r1_out)
    stats["r2_trimmed"] = str(r2_out)
    stats["json"] = str(json_out)
    stats["html"] = str(html_out)

    success(f"[{sample_id}] fastp OK — "
            f"{stats['reads_after']:,} reads kept "
            f"({stats['pass_rate']:.1f}%)")

    return stats


def _parse_fastp_json(json_path: Path, sample_id: str,
                      min_length_warning: int = 75) -> dict:
    """Parse fastp JSON output and check QC thresholds."""
    with open(json_path) as fh:
        data = json.load(fh)

    summary = data.get("summary", {})
    before  = summary.get("before_filtering", {})
    after   = summary.get("after_filtering", {})

    reads_before = before.get("total_reads", 0)
    reads_after  = after.get("total_reads", 0)
    pass_rate    = (reads_after / reads_before * 100) if reads_before else 0
    mean_len     = after.get("read1_mean_length", 0)
    q30_rate     = after.get("q30_rate", 0) * 100

    # Check for short reads after trimming
    if mean_len < min_length_warning:
        warn(f"[{sample_id}] Mean read length after trimming = {mean_len}bp "
             f"(< {min_length_warning}bp warning threshold). "
             "Possible RNA degradation or short insert size.")

    if pass_rate < 50:
        warn(f"[{sample_id}] Low pass rate: {pass_rate:.1f}%. "
             "Check library quality.")

    return {
        "sample_id":     sample_id,
        "reads_before":  reads_before,
        "reads_after":   reads_after,
        "pass_rate":     pass_rate,
        "mean_len_r1":   mean_len,
        "q30_rate":      q30_rate,
        "duplication":   data.get("duplication", {}).get("rate", 0) * 100,
    }


# ─────────────────────────────────────────────────────────────────────────────
# Run all samples
# ─────────────────────────────────────────────────────────────────────────────

def run_trimming_all(
    manifest_path: str | Path,
    out_dir: str | Path,
    threads: int = 4,
    min_phred: int = 20,
    min_length: int = 36,
    min_length_warning: int = 75,
) -> list[dict]:
    """
    Run fastp on all samples in manifest.
    Returns list of per-sample stats dicts.
    """
    section("M02 — Trimming (fastp)")
    samples = load_sample_manifest(manifest_path)
    all_stats = []

    for s in samples:
        stats = run_fastp(
            sample_id=s["sample_id"],
            r1_in=s["R1"],
            r2_in=s["R2"],
            out_dir=Path(out_dir) / s["sample_id"],
            threads=threads,
            min_phred=min_phred,
            min_length=min_length,
            min_length_warning=min_length_warning,
        )
        stats["condition"] = s.get("condition", "unknown")
        all_stats.append(stats)

    # Write summary TSV
    summary_path = Path(out_dir) / "trim_stats.tsv"
    _write_trim_summary(all_stats, summary_path)

    section("M02 — Trimming summary")
    for s in all_stats:
        info(f"  {s['sample_id']:20} {s['condition']:10} "
             f"reads: {s['reads_before']:>10,} -> {s['reads_after']:>10,} "
             f"({s['pass_rate']:.1f}%)  "
             f"Q30: {s['q30_rate']:.1f}%")

    return all_stats


def _write_trim_summary(stats_list: list[dict], out_path: Path) -> None:
    if not stats_list:
        return
    fields = ["sample_id", "condition", "reads_before", "reads_after",
              "pass_rate", "mean_len_r1", "q30_rate", "duplication",
              "r1_trimmed", "r2_trimmed"]
    with open(out_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields,
                                delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(stats_list)
    info(f"Trim summary written: {out_path}")


# ─────────────────────────────────────────────────────────────────────────────
# CLI entry point (stand-alone use)
# ─────────────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="MOSAIC M02 — fastp trimming")
    parser.add_argument("--manifest", required=True, help="samples.tsv")
    parser.add_argument("--outdir",   required=True, help="Output directory")
    parser.add_argument("--threads",  type=int, default=4)
    args = parser.parse_args()

    from .log import setup_logging
    setup_logging()
    run_trimming_all(args.manifest, args.outdir, threads=args.threads)
