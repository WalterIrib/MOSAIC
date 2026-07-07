"""
mosaic/m06_circrna.py
M06 — circRNA detection (CIRI2) using STAR chimeric alignments.
STAR chimeric SAMs from M04 replace BWA-MEM alignments.
Eliminates BWA dependency; chimeric detection is more sensitive with STAR.
"""
from __future__ import annotations
import csv, subprocess
from pathlib import Path
from typing import Optional
from .log import info, warn, error, success, section


# ─────────────────────────────────────────────────────────────────────────────
# CIRI2 runner (input: STAR Chimeric.out.sam)
# ─────────────────────────────────────────────────────────────────────────────

def run_ciri2(
    chimeric_sam: Path,
    genome_fa: Path,
    gtf: Path,
    out_dir: Path,
    sample_id: str,
    ciri2_pl: str = "CIRI2.pl",
    threads: int = 16,
    min_junction_reads: int = 2,
) -> Optional[Path]:
    """
    Run CIRI2 on STAR Chimeric.out.sam.

    STAR chimeric SAM is compatible with CIRI2 as BWA-MEM SAM replacement.
    Key: STAR must be run with --chimOutType SeparateSAMold and
         --chimSegmentMin >=20 for optimal junction detection.

    Returns path to .ciri output file, or None on failure.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    out_file = out_dir / f"{sample_id}.ciri"
    log_file = out_dir / f"{sample_id}_ciri2.log"

    if not chimeric_sam.exists() or chimeric_sam.stat().st_size == 0:
        warn(f"[{sample_id}] No chimeric SAM — skipping CIRI2")
        return None

    info(f"[{sample_id}] CIRI2 on {chimeric_sam.name}...")

    cmd = [
        "perl", ciri2_pl,
        "-I", str(chimeric_sam),
        "-O", str(out_file),
        "-F", str(genome_fa),
        "-A", str(gtf),
        "-T", str(threads),
        "-0",   # report all junction reads (min=0; we filter later)
    ]

    with open(log_file, "w") as lf:
        result = subprocess.run(cmd, stdout=lf, stderr=subprocess.STDOUT, text=True)

    if result.returncode != 0:
        error(f"[{sample_id}] CIRI2 failed (check {log_file})")
        return None

    if not out_file.exists() or out_file.stat().st_size == 0:
        warn(f"[{sample_id}] CIRI2 produced empty output")
        return None

    n = _count_ciri_entries(out_file, min_junction_reads)
    info(f"[{sample_id}] CIRI2: {n} circRNA (>= {min_junction_reads} junction reads)")
    return out_file


def _count_ciri_entries(ciri_file: Path, min_junc: int = 2) -> int:
    count = 0
    with open(ciri_file) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) >= 5:
                try:
                    if int(parts[4]) >= min_junc:
                        count += 1
                except ValueError:
                    pass
    return count


# ─────────────────────────────────────────────────────────────────────────────
# Parse CIRI2 output
# ─────────────────────────────────────────────────────────────────────────────

def parse_ciri2(ciri_file: Path, sample_id: str,
                min_junction_reads: int = 2) -> list[dict]:
    """
    Parse CIRI2 .ciri output.
    Returns list of dicts with circRNA locus info.

    CIRI2 columns (v2.0.6):
    0: circRNA_ID, 1: chr, 2: circRNA_start, 3: circRNA_end,
    4: #junction_reads, 5: SM_MS_SMS, 6: #non_junction_reads,
    7: junction_reads_ratio, 8: circRNA_type, 9: gene_id,
    10: strand, 11: junction_reads_ID
    """
    records = []
    if ciri_file is None or not ciri_file.exists():
        return records

    with open(ciri_file) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 10:
                continue
            try:
                junc_reads = int(parts[4])
            except ValueError:
                continue
            if junc_reads < min_junction_reads:
                continue
            records.append({
                "sample_id":      sample_id,
                "circRNA_id":     parts[0],
                "chr":            parts[1],
                "start":          int(parts[2]),
                "end":            int(parts[3]),
                "junction_reads": junc_reads,
                "circ_type":      parts[8] if len(parts) > 8 else "unknown",
                "gene_id":        parts[9] if len(parts) > 9 else "",
                "strand":         parts[10] if len(parts) > 10 else ".",
            })
    return records


# ─────────────────────────────────────────────────────────────────────────────
# Merge circRNA across samples
# ─────────────────────────────────────────────────────────────────────────────

def merge_circrna(
    sample_records: dict[str, list[dict]],
    coord_tolerance: int = 10,
) -> tuple[dict, dict]:
    """
    Merge circRNA loci across samples using coordinate tolerance.
    Returns (count_matrix, annotation).
    count_matrix: {locus_id: {sample_id: n_junc_reads}}
    annotation:   {locus_id: {chr, start, end, circ_type, gene_id, strand, n_samples}}
    """
    # Collect all loci
    all_loci = {}  # locus_id -> canonical coords

    for sid, records in sample_records.items():
        for r in records:
            matched = None
            for lid, locus in all_loci.items():
                if (r["chr"] == locus["chr"] and
                    abs(r["start"] - locus["start"]) <= coord_tolerance and
                    abs(r["end"]   - locus["end"])   <= coord_tolerance):
                    matched = lid
                    break
            if matched is None:
                lid = f"{r['chr']}:{r['start']}-{r['end']}"
                all_loci[lid] = {
                    "chr":       r["chr"],
                    "start":     r["start"],
                    "end":       r["end"],
                    "circ_type": r["circ_type"],
                    "gene_id":   r["gene_id"],
                    "strand":    r["strand"],
                }

    # Build count matrix
    all_samples = list(sample_records.keys())
    count_matrix = {lid: {s: 0 for s in all_samples} for lid in all_loci}

    for sid, records in sample_records.items():
        for r in records:
            for lid, locus in all_loci.items():
                if (r["chr"] == locus["chr"] and
                    abs(r["start"] - locus["start"]) <= coord_tolerance and
                    abs(r["end"]   - locus["end"])   <= coord_tolerance):
                    count_matrix[lid][sid] += r["junction_reads"]
                    break

    # Annotation with n_samples_detected
    annotation = {}
    for lid, locus in all_loci.items():
        n_det = sum(1 for s in all_samples if count_matrix[lid][s] > 0)
        annotation[lid] = {**locus, "n_samples_detected": n_det}

    return count_matrix, annotation


# ─────────────────────────────────────────────────────────────────────────────
# Differential expression (Wilcoxon fallback for n < 3)
# ─────────────────────────────────────────────────────────────────────────────

def circrna_de(
    count_matrix: dict,
    annotation: dict,
    conditions: dict[str, str],
    treatment: str,
    control: str,
    padj_threshold: float = 0.05,
) -> list[dict]:
    """
    DE for circRNA using Wilcoxon rank-sum (scipy) — robust for small n.
    Falls back to fold-change only if n < 2 per group.
    """
    try:
        from scipy.stats import mannwhitneyu
        import numpy as np
    except ImportError:
        warn("scipy not available — skipping circRNA DE")
        return []

    treat_samples = [s for s, c in conditions.items() if c == treatment]
    ctrl_samples  = [s for s, c in conditions.items() if c == control]

    results = []
    loci_ids = list(count_matrix.keys())
    pvalues  = []

    for lid in loci_ids:
        t_counts = [count_matrix[lid].get(s, 0) for s in treat_samples]
        c_counts = [count_matrix[lid].get(s, 0) for s in ctrl_samples]

        mean_t = np.mean(t_counts) if t_counts else 0
        mean_c = np.mean(c_counts) if c_counts else 0
        lfc = np.log2((mean_t + 1) / (mean_c + 1))

        if len(treat_samples) >= 2 and len(ctrl_samples) >= 2:
            _, pval = mannwhitneyu(t_counts, c_counts, alternative="two-sided")
        else:
            pval = 1.0

        pvalues.append(pval)
        results.append({
            "circRNA_id":       lid,
            "chr":              annotation[lid]["chr"],
            "start":            annotation[lid]["start"],
            "end":              annotation[lid]["end"],
            "circ_type":        annotation[lid]["circ_type"],
            "gene_id":          annotation[lid]["gene_id"],
            "n_samples_detected": annotation[lid]["n_samples_detected"],
            "mean_treatment":   round(mean_t, 3),
            "mean_control":     round(mean_c, 3),
            "log2FoldChange":   round(lfc, 4),
            "pvalue":           pval,
        })

    # BH correction
    try:
        from statsmodels.stats.multitest import multipletests
        import numpy as np
        _, padjs, _, _ = multipletests(np.array(pvalues), method="fdr_bh")
    except ImportError:
        # manual BH
        import numpy as np
        n = len(pvalues)
        ranked = np.argsort(pvalues)
        padjs  = np.empty(n)
        for i, idx in enumerate(ranked):
            padjs[idx] = min(pvalues[idx] * n / (i + 1), 1.0)

    for r, padj in zip(results, padjs):
        r["padj"] = round(float(padj), 6)

    results.sort(key=lambda x: x["pvalue"])
    return results


# ─────────────────────────────────────────────────────────────────────────────
# Main M06 runner
# ─────────────────────────────────────────────────────────────────────────────

def run_m06(cfg, chimeric_sams: dict[str, Path]) -> dict:
    """
    M06 pipeline:
    1. Run CIRI2 on each STAR chimeric SAM
    2. Merge circRNA loci across samples
    3. Build count matrix and annotation
    4. Differential expression (Wilcoxon)
    """
    section("M06 — circRNA detection (CIRI2 + STAR chimeric)")

    out_dir        = Path(cfg.output_dir) / "circrna"
    merged_dir     = Path(cfg.output_dir) / "circrna_merged"
    out_dir.mkdir(parents=True, exist_ok=True)
    merged_dir.mkdir(parents=True, exist_ok=True)

    genome_fa = Path(cfg.genome_fa)
    gtf       = Path(cfg.gtf)
    ciri2_pl  = getattr(cfg, "ciri2_path", "CIRI2.pl")
    threads   = getattr(cfg, "threads", 16)
    min_junc  = getattr(cfg, "min_junction_reads", 2)
    treatment = getattr(cfg, "treatment", "tumor")
    control   = getattr(cfg, "control",   "normal")
    samples   = cfg.samples

    conditions = {s["sample_id"]: s["condition"] for s in samples}

    # 1. CIRI2 per sample
    sample_records: dict[str, list[dict]] = {}
    for s in samples:
        sid = s["sample_id"]
        chim = chimeric_sams.get(sid)
        if chim is None:
            warn(f"[{sid}] No chimeric SAM from M04 — skipping")
            sample_records[sid] = []
            continue

        ciri_out = run_ciri2(
            chimeric_sam=chim,
            genome_fa=genome_fa,
            gtf=gtf,
            out_dir=out_dir,
            sample_id=sid,
            ciri2_pl=ciri2_pl,
            threads=threads,
            min_junction_reads=min_junc,
        )
        sample_records[sid] = parse_ciri2(ciri_out, sid, min_junc) if ciri_out else []

    total = sum(len(v) for v in sample_records.values())
    info(f"CIRI2 total records across samples: {total}")

    # 2. Merge
    count_matrix, annotation = merge_circrna(sample_records, coord_tolerance=10)
    info(f"Merged circRNA loci: {len(count_matrix)}")

    # 3. Save count matrix
    all_samples = [s["sample_id"] for s in samples]
    counts_path = merged_dir / "circrna_counts.tsv"
    with open(counts_path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["circRNA_id"] + all_samples)
        for lid in sorted(count_matrix):
            w.writerow([lid] + [count_matrix[lid].get(s, 0) for s in all_samples])

    # 4. Save annotation
    annot_path = merged_dir / "circrna_annotation.tsv"
    with open(annot_path, "w", newline="") as f:
        if annotation:
            w = csv.DictWriter(f, fieldnames=list(next(iter(annotation.values())).keys()) +
                               ["circRNA_id"], delimiter="\t")
            w.writeheader()
            for lid, ann in sorted(annotation.items()):
                w.writerow({**ann, "circRNA_id": lid})

    # 5. DE
    de_results = circrna_de(count_matrix, annotation, conditions,
                            treatment, control)
    de_path = merged_dir / "circrna_DE_wilcoxon.tsv"
    if de_results:
        with open(de_path, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=de_results[0].keys(), delimiter="\t")
            w.writeheader(); w.writerows(de_results)

    success(f"M06 complete — {len(count_matrix)} circRNA loci")
    return {
        "count_matrix": count_matrix,
        "annotation":   annotation,
        "de_results":   de_results,
        "counts_path":  counts_path,
        "annot_path":   annot_path,
        "de_path":      de_path,
    }
