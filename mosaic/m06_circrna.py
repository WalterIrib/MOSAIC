"""
mosaic/m06_circrna.py
M06 — circRNA detection (DCC) using STAR chimeric junction files.

M04 must produce Chimeric.out.junction via:
  --chimOutType SeparateSAMold Junctions
Chimeric.out.sam is unchanged (SeparateSAMold remains); Junctions adds the
junction file required by DCC.  Zero impact on SAM content (Opzione C,
validated: 411,561 alignment lines, identical to canonical).

DCC ref: Zeng et al., Nucleic Acids Res. 2017; tool: bioconda dcc=0.5.0.
"""
from __future__ import annotations
import csv, subprocess
from pathlib import Path
from typing import Optional
from .log import info, warn, error, success, section

# Maps M07 strandedness labels to DCC flags.
# Convention in this codebase (pipeline.py strandedness_to_featurecounts):
#   RF = "1++,1--" = R1 same as gene = featureCounts -s 1
#   FR = "1+-,1-+" = R2 same as gene = featureCounts -s 2
# DCC -s  = sense  (R1 on same strand as circRNA) → matches RF
# DCC -ss = anti-sense (R1 opposite to circRNA)   → matches FR
# Validated: De Cecco M07-FR samples ran correctly with DCC -ss in batch script.
_STRAND_MAP = {
    "unstranded": "-N",
    "auto":       "-N",   # safe fallback when M07 not yet run
    "rf":         "-s",   # R1 sense → DCC -s
    "fr":         "-ss",  # R2 sense, R1 antisense → DCC -ss
}


def _strand_flag(strandedness: str) -> str:
    return _STRAND_MAP.get(strandedness.lower().strip(), "-N")


# ─────────────────────────────────────────────────────────────────────────────
# DCC runner
# ─────────────────────────────────────────────────────────────────────────────

def run_dcc(
    junction_file: Path,
    gtf: Path,
    out_dir: Path,
    sample_id: str,
    dcc_bin: str = "DCC",
    strandedness: str = "unstranded",
    threads: int = 8,
    min_junction_reads: int = 2,
) -> Optional[Path]:
    """
    Run DCC on STAR Chimeric.out.junction.
    Returns path to CircRNACount, or None on failure.

    DCC exits RC=1 at the optional CircSkipJunction step (GENCODE exon_id
    trailing-quote bug in Circ_nonCirc_Exon_Match.py line 222) — this occurs
    AFTER CircRNACount is written and is non-fatal for our use case.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    count_file = out_dir / "CircRNACount"
    log_file   = out_dir / "dcc.log"

    if not junction_file.exists() or junction_file.stat().st_size == 0:
        warn(f"[{sample_id}] No junction file — skipping DCC")
        return None

    sflag = _strand_flag(strandedness)
    info(f"[{sample_id}] DCC on {junction_file.name} ({sflag})...")

    cmd = [
        dcc_bin,
        str(junction_file),
        "-D",
        sflag,
        "-F",
        "-Nr", str(min_junction_reads), "1",
        "-an", str(gtf),
        "-T", str(threads),
        "-O", str(out_dir),
    ]

    with open(log_file, "w") as lf:
        subprocess.run(cmd, stdout=lf, stderr=subprocess.STDOUT, text=True)

    if not count_file.exists() or count_file.stat().st_size == 0:
        error(f"[{sample_id}] DCC failed — CircRNACount missing (check {log_file})")
        return None

    n = _count_entries(count_file, min_junction_reads)
    info(f"[{sample_id}] DCC: {n} circRNA (>= {min_junction_reads} junction reads)")
    return count_file


def _count_entries(count_file: Path, min_junc: int = 2) -> int:
    count = 0
    with open(count_file) as f:
        for i, line in enumerate(f):
            if i == 0:
                continue
            parts = line.strip().split("\t")
            if len(parts) >= 4:
                try:
                    if int(parts[3]) >= min_junc:
                        count += 1
                except ValueError:
                    pass
    return count


# ─────────────────────────────────────────────────────────────────────────────
# Parse DCC output
# ─────────────────────────────────────────────────────────────────────────────

def parse_dcc(
    out_dir: Path,
    sample_id: str,
    min_junction_reads: int = 2,
) -> list[dict]:
    """
    Parse DCC CircRNACount + CircCoordinates from out_dir.

    CircRNACount columns: Chr, Start, End, <sample>_Chimeric.out.junction
    CircCoordinates cols: Chr, Start, End, Gene, JunctionType, Strand, ...
    """
    count_file = out_dir / "CircRNACount"
    coord_file = out_dir / "CircCoordinates"
    records: list[dict] = []

    if not count_file.exists():
        return records

    # Load gene annotation from CircCoordinates (keyed by chr:start:end)
    coord_map: dict[tuple, dict] = {}
    if coord_file.exists():
        with open(coord_file) as f:
            for i, line in enumerate(f):
                if i == 0:
                    continue
                parts = line.strip().split("\t")
                if len(parts) >= 6:
                    key = (parts[0], parts[1], parts[2])
                    coord_map[key] = {
                        "gene_id": parts[3],
                        "strand":  parts[5],
                    }

    with open(count_file) as f:
        for i, line in enumerate(f):
            if i == 0:
                continue
            parts = line.strip().split("\t")
            if len(parts) < 4:
                continue
            try:
                junc_reads = int(parts[3])
            except ValueError:
                continue
            if junc_reads < min_junction_reads:
                continue
            key = (parts[0], parts[1], parts[2])
            ann = coord_map.get(key, {})
            records.append({
                "sample_id":      sample_id,
                "circRNA_id":     f"{parts[0]}:{parts[1]}-{parts[2]}",
                "chr":            parts[0],
                "start":          int(parts[1]),
                "end":            int(parts[2]),
                "junction_reads": junc_reads,
                "gene_id":        ann.get("gene_id", "not_annotated"),
                "strand":         ann.get("strand", "."),
                "circ_type":      "exon",
            })
    return records


# ─────────────────────────────────────────────────────────────────────────────
# Merge across samples
# ─────────────────────────────────────────────────────────────────────────────

def merge_circrna(
    sample_records: dict[str, list[dict]],
) -> tuple[dict, dict]:
    """
    Merge circRNA loci across samples using exact coordinates.
    DCC uses the same junction file format for all samples, so coordinates
    are reproducible and exact matching is correct (no tolerance needed).

    Returns (count_matrix, annotation):
      count_matrix: {locus_id: {sample_id: junction_reads}}
      annotation:   {locus_id: {chr, start, end, gene_id, strand, n_samples_detected}}
    """
    all_loci: dict[str, dict] = {}

    for sid, records in sample_records.items():
        for r in records:
            lid = r["circRNA_id"]
            if lid not in all_loci:
                all_loci[lid] = {
                    "chr":       r["chr"],
                    "start":     r["start"],
                    "end":       r["end"],
                    "gene_id":   r["gene_id"],
                    "strand":    r["strand"],
                    "circ_type": r["circ_type"],
                }

    all_samples = list(sample_records.keys())
    count_matrix = {lid: {s: 0 for s in all_samples} for lid in all_loci}

    for sid, records in sample_records.items():
        for r in records:
            lid = r["circRNA_id"]
            if lid in count_matrix:
                count_matrix[lid][sid] = r["junction_reads"]

    annotation = {}
    for lid, locus in all_loci.items():
        n_det = sum(1 for s in all_samples if count_matrix[lid][s] > 0)
        annotation[lid] = {**locus, "n_samples_detected": n_det}

    return count_matrix, annotation


# ─────────────────────────────────────────────────────────────────────────────
# Differential expression (Wilcoxon fallback for small n)
# ─────────────────────────────────────────────────────────────────────────────

def circrna_de(
    count_matrix: dict,
    annotation: dict,
    conditions: dict[str, str],
    treatment: str,
    control: str,
) -> list[dict]:
    """
    DE for circRNA using Wilcoxon rank-sum (scipy) — robust for small n.
    Falls back to fold-change only if n < 2 per group.
    BH correction applied when statsmodels is available.
    """
    try:
        from scipy.stats import mannwhitneyu
        import numpy as np
    except ImportError:
        warn("scipy not available — skipping circRNA DE")
        return []

    treat_samples = [s for s, c in conditions.items() if c == treatment]
    ctrl_samples  = [s for s, c in conditions.items() if c == control]

    results: list[dict] = []
    pvalues: list[float] = []

    for lid in count_matrix:
        t_counts = [count_matrix[lid].get(s, 0) for s in treat_samples]
        c_counts = [count_matrix[lid].get(s, 0) for s in ctrl_samples]

        mean_t = float(np.mean(t_counts)) if t_counts else 0.0
        mean_c = float(np.mean(c_counts)) if c_counts else 0.0
        lfc    = float(np.log2((mean_t + 1) / (mean_c + 1)))

        if len(treat_samples) >= 2 and len(ctrl_samples) >= 2:
            _, pval = mannwhitneyu(t_counts, c_counts, alternative="two-sided")
        else:
            pval = 1.0

        pvalues.append(pval)
        results.append({
            "circRNA_id":         lid,
            "chr":                annotation[lid]["chr"],
            "start":              annotation[lid]["start"],
            "end":                annotation[lid]["end"],
            "gene_id":            annotation[lid]["gene_id"],
            "strand":             annotation[lid]["strand"],
            "n_samples_detected": annotation[lid]["n_samples_detected"],
            "mean_treatment":     round(mean_t, 3),
            "mean_control":       round(mean_c, 3),
            "log2FoldChange":     round(lfc, 4),
            "pvalue":             pval,
        })

    try:
        from statsmodels.stats.multitest import multipletests
        _, padjs, _, _ = multipletests(np.array(pvalues), method="fdr_bh")
    except ImportError:
        n = len(pvalues)
        ranked = np.argsort(pvalues)
        padjs = np.empty(n)
        for i, idx in enumerate(ranked):
            padjs[idx] = min(pvalues[idx] * n / (i + 1), 1.0)

    for r, padj in zip(results, padjs):
        r["padj"] = round(float(padj), 6)

    results.sort(key=lambda x: x["pvalue"])
    return results


# ─────────────────────────────────────────────────────────────────────────────
# Main M06 entry point
# ─────────────────────────────────────────────────────────────────────────────

def run_m06(cfg, chimeric_junctions: dict[str, Path]) -> dict:
    """
    M06 pipeline:
    1. Run DCC on each Chimeric.out.junction
    2. Merge circRNA loci across samples (exact coordinates)
    3. Build count matrix + annotation TSVs
    4. Differential expression (Wilcoxon rank-sum + BH)
    """
    section("M06 — circRNA detection (DCC + STAR chimeric junctions)")

    out_dir    = Path(cfg.output_dir) / "circrna"
    merged_dir = Path(cfg.output_dir) / "circrna_merged"
    out_dir.mkdir(parents=True, exist_ok=True)
    merged_dir.mkdir(parents=True, exist_ok=True)

    gtf          = Path(cfg.gtf)
    dcc_bin      = getattr(cfg, "dcc_path",           "DCC")
    threads      = getattr(cfg, "threads",             8)
    min_junc     = getattr(cfg, "min_junction_reads",  2)
    treatment    = getattr(cfg, "treatment",           "tumor")
    control      = getattr(cfg, "control",             "normal")
    strandedness = getattr(cfg, "strandedness",        "unstranded")
    samples      = cfg.samples

    conditions = {s["sample_id"]: s["condition"] for s in samples}

    # 1. DCC per sample
    sample_records: dict[str, list[dict]] = {}
    for s in samples:
        sid = s["sample_id"]
        jxn = chimeric_junctions.get(sid)
        if jxn is None:
            warn(f"[{sid}] No junction file from M04 — skipping")
            sample_records[sid] = []
            continue

        sample_out = out_dir / sid
        # Per-sample strandedness override from coldata if available
        s_strand = s.get("strandedness", strandedness)
        count_file = run_dcc(
            junction_file=jxn,
            gtf=gtf,
            out_dir=sample_out,
            sample_id=sid,
            dcc_bin=dcc_bin,
            strandedness=s_strand,
            threads=threads,
            min_junction_reads=min_junc,
        )
        sample_records[sid] = parse_dcc(sample_out, sid, min_junc) if count_file else []

    total = sum(len(v) for v in sample_records.values())
    info(f"DCC total records across samples: {total}")

    # 2. Merge loci
    count_matrix, annotation = merge_circrna(sample_records)
    info(f"Merged circRNA loci: {len(count_matrix)}")

    # 3. Count matrix TSV
    all_samples = [s["sample_id"] for s in samples]
    counts_path = merged_dir / "circrna_counts.tsv"
    with open(counts_path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["circRNA_id"] + all_samples)
        for lid in sorted(count_matrix):
            w.writerow([lid] + [count_matrix[lid].get(s, 0) for s in all_samples])

    # 4. Annotation TSV
    annot_path = merged_dir / "circrna_annotation.tsv"
    with open(annot_path, "w", newline="") as f:
        if annotation:
            first = next(iter(annotation.values()))
            w = csv.DictWriter(
                f,
                fieldnames=["circRNA_id"] + list(first.keys()),
                delimiter="\t",
            )
            w.writeheader()
            for lid, ann in sorted(annotation.items()):
                w.writerow({"circRNA_id": lid, **ann})

    # 5. DE
    de_results = circrna_de(count_matrix, annotation, conditions, treatment, control)
    de_path = merged_dir / "circrna_DE_wilcoxon.tsv"
    if de_results:
        with open(de_path, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=de_results[0].keys(), delimiter="\t")
            w.writeheader()
            w.writerows(de_results)

    success(f"M06 complete — {len(count_matrix)} circRNA loci")
    return {
        "count_matrix": count_matrix,
        "annotation":   annotation,
        "de_results":   de_results,
        "counts_path":  counts_path,
        "annot_path":   annot_path,
        "de_path":      de_path,
    }
