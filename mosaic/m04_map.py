"""
mosaic/m04_map.py
M04 — Mapping canonico (STAR two-pass) + strandedness (RSeQC) + featureCounts.
Sostituisce HISAT2+BWA con STAR two-pass mode.
STAR genera anche chimeric alignments per CIRI2 (M06), eliminando BWA.
"""
from __future__ import annotations
import csv, os, re, shutil, subprocess
from pathlib import Path
from typing import Optional
from .log import info, warn, error, success, section


# ─────────────────────────────────────────────────────────────────────────────
# STAR genome index
# ─────────────────────────────────────────────────────────────────────────────

def build_star_index(
    genome_fa: Path,
    gtf: Path,
    index_dir: Path,
    threads: int = 16,
    ram_gb: int = 28,
    overhang: int = 149,
) -> Path:
    """
    Build STAR genome index.
    overhang = readLength - 1 (default 149 for 150bp reads).
    hg38 requires ~30GB RAM.
    """
    index_dir.mkdir(parents=True, exist_ok=True)
    if (index_dir / "SA").exists():
        info(f"STAR index already exists: {index_dir}")
        return index_dir

    section("Building STAR genome index")
    cmd = [
        "STAR", "--runMode", "genomeGenerate",
        "--genomeDir",              str(index_dir),
        "--genomeFastaFiles",       str(genome_fa),
        "--sjdbGTFfile",            str(gtf),
        "--sjdbOverhang",           str(overhang),
        "--runThreadN",             str(threads),
        "--limitGenomeGenerateRAM", str(ram_gb * 1_073_741_824),
        "--genomeSAindexNbases",    "14",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        error(f"STAR genomeGenerate failed:\n{result.stderr[-2000:]}")
        raise RuntimeError("STAR index build failed")
    success(f"STAR index built: {index_dir}")
    return index_dir


# ─────────────────────────────────────────────────────────────────────────────
# STAR two-pass alignment
# ─────────────────────────────────────────────────────────────────────────────

def star_align(
    r1: Path,
    r2: Optional[Path],
    index_dir: Path,
    out_dir: Path,
    sample_id: str,
    threads: int = 16,
    ram_gb: int = 28,
    chimeric: bool = True,
    chim_segment_min: int = 20,
    chim_score_drop_max: int = 20,
    two_pass: bool = True,
) -> dict:
    """
    STAR alignment with optional two-pass mode and chimeric output.

    Two-pass mode: discovers novel junctions in pass1, re-aligns in pass2.
    Chimeric output: SeparateSAMold produces Chimeric.out.sam for CIRI2.
    Returns dict: bam, chimeric_sam, log_final, mapping_rate.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    prefix = str(out_dir / f"{sample_id}_")

    reads_arg = [str(r1), str(r2)] if r2 else [str(r1)]
    gz = str(r1).endswith(".gz")

    cmd = [
        "STAR",
        "--runMode",            "alignReads",
        "--genomeDir",          str(index_dir),
        "--readFilesIn",        *reads_arg,
        "--outFileNamePrefix",  prefix,
        "--runThreadN",         str(threads),
        "--limitBAMsortRAM",    str(ram_gb * 1_073_741_824),
        # Output format
        "--outSAMtype",         "BAM", "SortedByCoordinate",
        "--outSAMattributes",   "NH", "HI", "AS", "NM", "MD",
        "--outSAMstrandField",  "intronMotif",
        "--outFilterIntronMotifs", "RemoveNoncanonical",
        # Mapping params (ENCODE standard)
        "--twopassMode",        "Basic" if two_pass else "None",
        "--outFilterMultimapNmax", "20",
        "--alignSJoverhangMin",  "8",
        "--alignSJDBoverhangMin", "1",
        "--outFilterMismatchNmax", "999",
        "--outFilterMismatchNoverReadLengthMax", "0.04",
        "--alignIntronMin",     "20",
        "--alignIntronMax",     "1000000",
        "--alignMatesGapMax",   "1000000",
        "--readFilesCommand",   "zcat" if gz else "-",
        # Gene counts (produces ReadsPerGene.out.tab)
        "--quantMode",          "GeneCounts",
    ]

    if chimeric:
        cmd += [
            "--chimOutType",               "SeparateSAMold",
            "--chimSegmentMin",            str(chim_segment_min),
            "--chimScoreDropMax",          str(chim_score_drop_max),
            "--chimJunctionOverhangMin",   "15",
            "--chimScoreMin",              "1",
            "--chimScoreSeparation",       "10",
            "--chimMultimapScoreRange",    "1",
            "--chimMultimapNmax",          "20",
            "--chimNonchimScoreDropMin",   "10",
        ]

    info(f"[{sample_id}] STAR two-pass (chimeric={chimeric})...")
    result = subprocess.run(cmd, capture_output=True, text=True)

    bam_path  = Path(prefix + "Aligned.sortedByCoord.out.bam")
    chim_sam  = Path(prefix + "Chimeric.out.sam")
    log_final = Path(prefix + "Log.final.out")

    if result.returncode != 0:
        error(f"[{sample_id}] STAR failed:\n{result.stderr[-2000:]}")
        raise RuntimeError(f"STAR alignment failed for {sample_id}")

    subprocess.run(["samtools", "index", str(bam_path)], check=True)
    rate = _parse_star_log(log_final)
    info(f"[{sample_id}] Uniquely mapped: {rate:.1f}%")

    return {
        "bam":          bam_path,
        "chimeric_sam": chim_sam if chimeric and chim_sam.exists() else None,
        "log_final":    log_final,
        "mapping_rate": rate,
    }


def _parse_star_log(log_path: Path) -> float:
    if not log_path.exists():
        return 0.0
    for line in log_path.read_text().splitlines():
        if "Uniquely mapped reads %" in line:
            try:
                return float(line.split("|")[1].strip().rstrip("%"))
            except Exception:
                pass
    return 0.0


# ─────────────────────────────────────────────────────────────────────────────
# Strandedness (RSeQC)
# ─────────────────────────────────────────────────────────────────────────────

def infer_strandedness(
    bam: Path,
    bed12: Path,
    sample_id: str,
    n_reads: int = 200_000,
) -> dict:
    cmd = ["infer_experiment.py", "-r", str(bed12), "-i", str(bam), "-s", str(n_reads)]
    info(f"[{sample_id}] Inferring strandedness...")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        warn(f"[{sample_id}] infer_experiment.py failed — assuming unstranded")
        return {"sample_id": sample_id, "strand": "unstranded",
                "fc_strand": 0, "frac_fr": 0.0, "frac_rf": 0.0, "frac_failed": 0.0}
    return _parse_infer_experiment(result.stdout, sample_id)


def _parse_infer_experiment(output: str, sample_id: str) -> dict:
    frac_fr = frac_rf = frac_failed = 0.0
    for line in output.splitlines():
        line = line.strip()
        if "failed to determine" in line.lower():
            m = re.search(r"([\d.]+)%?$", line)
            if m:
                frac_failed = float(m.group(1))
                if frac_failed > 1: frac_failed /= 100
        elif "++" in line and "+-" not in line:
            m = re.search(r"([\d.]+)$", line)
            if m: frac_fr = float(m.group(1))
        elif "+-" in line and "++" not in line:
            m = re.search(r"([\d.]+)$", line)
            if m: frac_rf = float(m.group(1))

    strand = "FR" if frac_fr >= 0.70 else "RF" if frac_rf >= 0.70 else "unstranded"
    fc_strand = {"FR": 1, "RF": 2, "unstranded": 0}[strand]
    info(f"[{sample_id}] Strandedness: {strand} (FR={frac_fr:.2f} RF={frac_rf:.2f})")
    return {"sample_id": sample_id, "strand": strand, "fc_strand": fc_strand,
            "frac_fr": frac_fr, "frac_rf": frac_rf, "frac_failed": frac_failed}


def consensus_strand(strand_results: list[dict]) -> int:
    votes = [r["fc_strand"] for r in strand_results]
    c = max(set(votes), key=votes.count)
    label = {0: "unstranded", 1: "FR", 2: "RF"}[c]
    info(f"Strandedness consensus: {label} (featureCounts -s {c})")
    return c


# ─────────────────────────────────────────────────────────────────────────────
# featureCounts
# ─────────────────────────────────────────────────────────────────────────────

def run_featurecounts(
    bam_files: list[Path],
    gtf: Path,
    out_dir: Path,
    strand: int = 0,
    threads: int = 16,
    paired: bool = True,
    prefix: str = "all_genes",
) -> Path:
    out_dir.mkdir(parents=True, exist_ok=True)
    out_file = out_dir / f"{prefix}_counts.tsv"

    cmd = [
        "featureCounts",
        "-T", str(threads),
        "-a", str(gtf),
        "-o", str(out_file),
        "-t", "exon",
        "-g", "gene_id",
        "-s", str(strand),
        "--fracOverlap", "0.1",
        "--extraAttributes", "gene_name,gene_type",
    ]
    if paired:
        cmd += ["-p", "--countReadPairs"]
    cmd += [str(b) for b in bam_files]

    info(f"featureCounts: {len(bam_files)} BAMs → {out_file.name}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        error(f"featureCounts failed:\n{result.stderr[-1000:]}")
        raise RuntimeError("featureCounts failed")
    success(f"featureCounts: {out_file}")
    return out_file


# ─────────────────────────────────────────────────────────────────────────────
# Main M04 runner
# ─────────────────────────────────────────────────────────────────────────────

def run_m04(cfg) -> dict:
    """
    M04 pipeline:
    1. STAR index (build if missing)
    2. STAR two-pass alignment per sample + chimeric SAM for CIRI2
    3. RSeQC strandedness
    4. featureCounts
    """
    section("M04 — STAR two-pass + featureCounts")

    out_map = Path(cfg.output_dir) / "mapping"
    bam_dir = out_map / "bam"
    cnt_dir = out_map / "counts"
    bam_dir.mkdir(parents=True, exist_ok=True)
    cnt_dir.mkdir(parents=True, exist_ok=True)

    star_idx = Path(cfg.star_index)
    gtf      = Path(cfg.gtf)
    bed12    = Path(cfg.bed12) if getattr(cfg, "bed12", None) else None
    threads  = getattr(cfg, "threads", 16)
    ram_gb   = getattr(cfg, "star_ram_gb", 28)
    samples  = cfg.samples

    # 1. Build index if needed
    if not (star_idx / "SA").exists():
        build_star_index(Path(cfg.genome_fa), gtf, star_idx,
                         threads=threads, ram_gb=ram_gb)

    # 2. Align
    bam_paths, chimeric_sams = [], {}
    strand_results, mapping_stats = [], []
    paired = bool(samples[0].get("r2"))

    for s in samples:
        sid = s["sample_id"]
        r1  = Path(s["r1"])
        r2  = Path(s["r2"]) if s.get("r2") else None
        sout = bam_dir / sid

        res = star_align(r1=r1, r2=r2, index_dir=star_idx,
                         out_dir=sout, sample_id=sid,
                         threads=threads, ram_gb=ram_gb, chimeric=True)

        final_bam = bam_dir / f"{sid}.bam"
        shutil.copy2(str(res["bam"]), str(final_bam))
        subprocess.run(["samtools", "index", str(final_bam)], check=True)
        bam_paths.append(final_bam)

        if res["chimeric_sam"]:
            chimeric_sams[sid] = res["chimeric_sam"]

        mapping_stats.append({"sample_id": sid, "condition": s.get("condition",""),
                               "mapping_rate": res["mapping_rate"]})

        if bed12 and bed12.exists():
            strand_results.append(infer_strandedness(final_bam, bed12, sid))
        else:
            warn(f"[{sid}] No BED12 — assuming unstranded")
            strand_results.append({"sample_id": sid, "strand": "unstranded",
                                   "fc_strand": 0, "frac_fr": 0, "frac_rf": 0,
                                   "frac_failed": 0})

    # Save stats
    for path, rows in [(out_map/"strandedness.tsv", strand_results),
                       (out_map/"mapping_stats.tsv", mapping_stats)]:
        with open(path, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=rows[0].keys(), delimiter="\t")
            w.writeheader(); w.writerows(rows)

    # 3. featureCounts
    fc_strand   = consensus_strand(strand_results)
    counts_file = run_featurecounts(bam_files=bam_paths, gtf=gtf,
                                    out_dir=cnt_dir, strand=fc_strand,
                                    threads=threads, paired=paired)

    success("M04 complete")
    return {"bam_paths": bam_paths, "chimeric_sams": chimeric_sams,
            "counts_file": counts_file, "strandedness": strand_results,
            "fc_strand": fc_strand, "mapping_stats": mapping_stats}
