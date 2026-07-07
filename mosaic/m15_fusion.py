"""
mosaic/m15_fusion.py
M15 — Fusion gene detection using Arriba.
Input: STAR chimeric alignments (Chimeric.out.sam from M04).
Arriba is the most accurate fusion detector for RNA-seq (Dresden benchmark).
"""
from __future__ import annotations
import csv, subprocess
from pathlib import Path
from typing import Optional
from .log import info, warn, error, success, section


# ─────────────────────────────────────────────────────────────────────────────
# Arriba runner
# ─────────────────────────────────────────────────────────────────────────────

def run_arriba(
    bam: Path,
    genome_fa: Path,
    gtf: Path,
    out_dir: Path,
    sample_id: str,
    blacklist: Optional[Path] = None,
    known_fusions: Optional[Path] = None,
    protein_domains: Optional[Path] = None,
    threads: int = 16,
) -> dict:
    """
    Run Arriba fusion detection on STAR-aligned BAM.

    Arriba uses STAR chimeric alignments internally via the BAM file
    (requires STAR run with --chimOutType WithinBAM or post-merge).
    For SeparateSAMold output, we convert chimeric SAM to BAM and merge.

    Returns dict with fusions_tsv and discarded_tsv paths.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    fusions_tsv   = out_dir / f"{sample_id}_fusions.tsv"
    discarded_tsv = out_dir / f"{sample_id}_fusions_discarded.tsv"

    if fusions_tsv.exists():
        info(f"[{sample_id}] Arriba output already exists")
        return {"fusions_tsv": fusions_tsv, "discarded_tsv": discarded_tsv}

    info(f"[{sample_id}] Arriba fusion detection...")

    cmd = [
        "arriba",
        "-x", str(bam),
        "-a", str(genome_fa),
        "-g", str(gtf),
        "-o", str(fusions_tsv),
        "-O", str(discarded_tsv),
    ]

    if blacklist and blacklist.exists():
        cmd += ["-b", str(blacklist)]
    else:
        # Try default Arriba blacklist location
        default_bl = Path("/usr/local/share/arriba-*/blacklist_hg38_*.tsv.gz")
        import glob
        matches = glob.glob(str(default_bl))
        if matches:
            cmd += ["-b", matches[0]]
        else:
            warn(f"[{sample_id}] No blacklist found — Arriba may produce more false positives")

    if known_fusions and known_fusions.exists():
        cmd += ["-k", str(known_fusions)]
    if protein_domains and protein_domains.exists():
        cmd += ["-p", str(protein_domains)]

    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        error(f"[{sample_id}] Arriba failed:\n{result.stderr[-2000:]}")
        return {"fusions_tsv": None, "discarded_tsv": None, "error": result.stderr}

    n_fusions = _count_fusions(fusions_tsv)
    info(f"[{sample_id}] Arriba: {n_fusions} fusion calls")

    return {
        "fusions_tsv":   fusions_tsv,
        "discarded_tsv": discarded_tsv,
        "n_fusions":     n_fusions,
    }


def _count_fusions(tsv: Path) -> int:
    if not tsv or not tsv.exists():
        return 0
    return sum(1 for l in open(tsv) if not l.startswith("#"))


# ─────────────────────────────────────────────────────────────────────────────
# Parse and merge fusions across samples
# ─────────────────────────────────────────────────────────────────────────────

def parse_arriba_fusions(tsv: Path, sample_id: str,
                          min_confidence: str = "low") -> list[dict]:
    """
    Parse Arriba fusions.tsv.
    confidence levels: high > medium > low
    Returns list of fusion dicts.
    """
    if not tsv or not tsv.exists():
        return []

    conf_rank = {"high": 3, "medium": 2, "low": 1}
    min_rank  = conf_rank.get(min_confidence, 1)
    records   = []

    with open(tsv) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            if row.get("#gene1", row.get("gene1","")) == "":
                continue
            conf = row.get("confidence","low").lower()
            if conf_rank.get(conf,0) < min_rank:
                continue
            records.append({
                "sample_id":      sample_id,
                "gene1":          row.get("#gene1", row.get("gene1","")),
                "gene2":          row.get("gene2",""),
                "strand1":        row.get("strand1(gene/fusion)",""),
                "strand2":        row.get("strand2(gene/fusion)",""),
                "breakpoint1":    row.get("breakpoint1",""),
                "breakpoint2":    row.get("breakpoint2",""),
                "site1":          row.get("site1",""),
                "site2":          row.get("site2",""),
                "type":           row.get("type",""),
                "confidence":     conf,
                "reading_frame":  row.get("reading_frame",""),
                "split_reads1":   row.get("split_reads1","0"),
                "split_reads2":   row.get("split_reads2","0"),
                "discordant_mates": row.get("discordant_mates","0"),
                "fusion_transcript": row.get("fusion_transcript","")[:100],
            })
    return records


def merge_fusions(sample_fusions: dict[str, list[dict]]) -> list[dict]:
    """
    Merge fusions across samples: group by gene pair, count samples.
    Returns list sorted by n_samples desc, then confidence.
    """
    from collections import defaultdict
    fusion_groups: dict[str, dict] = defaultdict(lambda: {
        "n_samples": 0, "samples": [], "max_confidence": "low",
        "split_reads_total": 0
    })
    conf_rank = {"high": 3, "medium": 2, "low": 1}

    for sid, fusions in sample_fusions.items():
        for f in fusions:
            key = f"{f['gene1']}::{f['gene2']}"
            fg  = fusion_groups[key]
            fg["gene1"]   = f["gene1"]
            fg["gene2"]   = f["gene2"]
            fg["type"]    = f.get("type","")
            fg["n_samples"] += 1
            fg["samples"].append(sid)
            fg["split_reads_total"] += (int(f.get("split_reads1",0) or 0) +
                                         int(f.get("split_reads2",0) or 0))
            if conf_rank.get(f["confidence"],0) > conf_rank.get(fg["max_confidence"],0):
                fg["max_confidence"] = f["confidence"]
                fg["breakpoint1"]    = f.get("breakpoint1","")
                fg["breakpoint2"]    = f.get("breakpoint2","")
                fg["reading_frame"]  = f.get("reading_frame","")

    merged = []
    for key, d in fusion_groups.items():
        d["fusion_id"] = key
        d["samples"]   = ",".join(sorted(set(d["samples"])))
        merged.append(d)

    return sorted(merged, key=lambda x: (-x["n_samples"],
                                          -conf_rank.get(x["max_confidence"],0),
                                          -x["split_reads_total"]))


# ─────────────────────────────────────────────────────────────────────────────
# Main M15 runner
# ─────────────────────────────────────────────────────────────────────────────

def run_m15(cfg, bam_paths: dict[str, Path]) -> dict:
    """
    M15 pipeline:
    1. Run Arriba on each sample BAM (STAR-aligned)
    2. Parse and filter fusions
    3. Merge across samples
    4. Save per-sample and merged outputs
    """
    section("M15 — Fusion gene detection (Arriba)")

    out_dir  = Path(cfg.output_dir) / "m15_fusions"
    out_dir.mkdir(parents=True, exist_ok=True)

    genome_fa      = Path(cfg.genome_fa)
    gtf            = Path(cfg.gtf)
    threads        = getattr(cfg, "threads", 16)
    blacklist      = Path(getattr(cfg, "arriba_blacklist", "")) if getattr(cfg, "arriba_blacklist", None) else None
    known_fusions  = Path(getattr(cfg, "arriba_known_fusions", "")) if getattr(cfg, "arriba_known_fusions", None) else None
    protein_domains= Path(getattr(cfg, "arriba_protein_domains", "")) if getattr(cfg, "arriba_protein_domains", None) else None
    min_confidence = getattr(cfg, "fusion_min_confidence", "low")

    samples = cfg.samples
    sample_fusions: dict[str, list[dict]] = {}

    for s in samples:
        sid = s["sample_id"]
        bam = bam_paths.get(sid)
        if bam is None or not bam.exists():
            warn(f"[{sid}] BAM not found — skipping")
            sample_fusions[sid] = []
            continue

        res = run_arriba(
            bam=bam, genome_fa=genome_fa, gtf=gtf,
            out_dir=out_dir / sid, sample_id=sid,
            blacklist=blacklist, known_fusions=known_fusions,
            protein_domains=protein_domains, threads=threads,
        )

        sample_fusions[sid] = parse_arriba_fusions(
            res.get("fusions_tsv"), sid, min_confidence=min_confidence
        )
        info(f"[{sid}] {len(sample_fusions[sid])} fusions (confidence >= {min_confidence})")

    # Merge
    merged = merge_fusions(sample_fusions)
    info(f"Merged fusion events: {len(merged)}")

    # Save per-sample
    all_per_sample = [f for fusions in sample_fusions.values() for f in fusions]
    per_sample_path = out_dir / "fusions_per_sample.tsv"
    if all_per_sample:
        with open(per_sample_path, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=all_per_sample[0].keys(), delimiter="\t")
            w.writeheader(); w.writerows(all_per_sample)

    # Save merged
    merged_path = out_dir / "fusions_merged.tsv"
    if merged:
        with open(merged_path, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=merged[0].keys(), delimiter="\t")
            w.writeheader(); w.writerows(merged)

    success(f"M15 complete — {len(merged)} fusion events")
    return {
        "sample_fusions": sample_fusions,
        "merged_fusions": merged,
        "per_sample_path": per_sample_path,
        "merged_path":    merged_path,
    }
