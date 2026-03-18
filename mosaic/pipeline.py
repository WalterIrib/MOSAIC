"""
mosaic/pipeline.py
MOSAIC pipeline orchestrator.

Runs modules in dependency order, supports:
- Multi-comparison (multiple treatment/control pairs)
- Incremental execution (skip completed modules)
- PE/SE and stranded/unstranded handling
- Multi-conda-env tool execution
- FASTQ download from ENA/SRA/SSH/URL
"""
from __future__ import annotations

import json
import logging
import os
import subprocess
import sys
import time
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import yaml

from .log import section, info, warn, error, success

log = logging.getLogger(__name__)


# ─────────────────────────────────────────────────────────────────────────────
# Module dependency graph
# ─────────────────────────────────────────────────────────────────────────────

MODULE_DEPS: Dict[str, List[str]] = {
    "m00_download":    [],
    "m02_trim":        ["m00_download"],
    "m04_star":        ["m02_trim"],
    "m05_te_count":    ["m02_trim"],
    "m06_circrna":     ["m04_star"],
    "m07_strandedness":["m04_star"],
    "m08_deseq2":      ["m04_star", "m05_te_count"],
    "m09_satellite":   ["m02_trim"],
    "m10_go":          ["m08_deseq2"],
    "m11_te_host":     ["m04_star"],
    "m11b_promoter_te":["m11_te_host"],
    "m11c_te_age":     ["m11_te_host", "m08_deseq2"],
    "m11d_te_de":      ["m08_deseq2", "m11b_promoter_te"],
    "m12_tetranscripts":["m04_star"],
    "m15_fusion":      ["m04_star"],
    "m16_stringtie":   ["m04_star"],
    "m17_rmats":       ["m04_star"],
    "m18_network":     ["m08_deseq2"],
    "m19_editing":     ["m04_star"],
    "m20_meta":        ["m04_star"],
    "m20b_humann":     ["m20_meta"],
    "m23_te_insertions":["m04_star"],
    "m14_report":      [
        "m08_deseq2", "m10_go", "m11b_promoter_te", "m11c_te_age",
        "m11d_te_de", "m12_tetranscripts", "m15_fusion", "m16_stringtie",
        "m17_rmats", "m18_network", "m19_editing", "m20_meta",
        "m20b_humann", "m23_te_insertions"
    ],
}

# Tools that require a specific conda env
TOOL_ENVS: Dict[str, str] = {
    "m20b_humann": "humann3",
    "m15_fusion":  "arriba_env",  # if not in main env
}


# ─────────────────────────────────────────────────────────────────────────────
# State tracking
# ─────────────────────────────────────────────────────────────────────────────

@dataclass
class ModuleState:
    name: str
    status: str = "pending"   # pending | running | done | failed | skipped
    start_time: Optional[float] = None
    end_time: Optional[float] = None
    comparison: Optional[str] = None
    error_msg: Optional[str] = None

    @property
    def elapsed(self) -> Optional[float]:
        if self.start_time and self.end_time:
            return self.end_time - self.start_time
        return None

    def to_dict(self) -> dict:
        return {
            "name": self.name,
            "status": self.status,
            "comparison": self.comparison,
            "elapsed_s": round(self.elapsed, 1) if self.elapsed else None,
            "error": self.error_msg,
        }


class PipelineState:
    """Tracks and persists pipeline execution state."""

    def __init__(self, state_file: Path):
        self.state_file = state_file
        self.modules: Dict[str, ModuleState] = {}
        self._load()

    def _load(self):
        if self.state_file.exists():
            try:
                data = json.loads(self.state_file.read_text())
                for name, d in data.items():
                    m = ModuleState(name=name)
                    m.status     = d.get("status", "pending")
                    m.comparison = d.get("comparison")
                    m.error_msg  = d.get("error")
                    self.modules[name] = m
            except Exception:
                warn(f"Could not load state from {self.state_file}")

    def save(self):
        self.state_file.parent.mkdir(parents=True, exist_ok=True)
        self.state_file.write_text(
            json.dumps({k: v.to_dict() for k, v in self.modules.items()},
                       indent=2))

    def is_done(self, module: str, comparison: Optional[str] = None) -> bool:
        key = f"{module}/{comparison}" if comparison else module
        m = self.modules.get(key)
        return m is not None and m.status == "done"

    def mark(self, module: str, status: str,
             comparison: Optional[str] = None,
             error: Optional[str] = None):
        key = f"{module}/{comparison}" if comparison else module
        m = self.modules.setdefault(key, ModuleState(name=key))
        m.status = status
        m.comparison = comparison
        m.error_msg = error
        if status == "running":
            m.start_time = time.time()
        elif status in ("done", "failed"):
            m.end_time = time.time()
        self.save()


# ─────────────────────────────────────────────────────────────────────────────
# Subprocess runner with conda env support
# ─────────────────────────────────────────────────────────────────────────────

def run_cmd(cmd: List[str], env_name: Optional[str] = None,
            cwd: Optional[Path] = None,
            log_file: Optional[Path] = None) -> int:
    """
    Run a command, optionally in a specific conda env.
    Returns exit code.
    """
    if env_name:
        cmd = ["conda", "run", "-n", env_name, "--no-capture-output"] + cmd

    cmd_str = " ".join(str(c) for c in cmd)
    info(f"$ {cmd_str}")

    kwargs: dict = dict(cwd=cwd)
    if log_file:
        log_file.parent.mkdir(parents=True, exist_ok=True)
        fh = open(log_file, "w")
        kwargs["stdout"] = fh
        kwargs["stderr"] = subprocess.STDOUT
    else:
        kwargs["stdout"] = None
        kwargs["stderr"] = None

    proc = subprocess.run(cmd, **kwargs)

    if log_file:
        fh.close()

    return proc.returncode


def run_r_script(script: Path, args: List[str],
                 log_file: Optional[Path] = None,
                 cwd: Optional[Path] = None) -> int:
    """Run an R script with Rscript."""
    return run_cmd(["Rscript", str(script)] + args,
                   log_file=log_file, cwd=cwd)


# ─────────────────────────────────────────────────────────────────────────────
# FASTQ downloader
# ─────────────────────────────────────────────────────────────────────────────

def download_fastq(source: str, accession: str,
                   outdir: Path, threads: int = 8) -> bool:
    """
    Download FASTQ files from various sources.

    source: ena | sra | url | ssh | local
    """
    outdir.mkdir(parents=True, exist_ok=True)

    if source == "ena":
        return _download_ena(accession, outdir, threads)
    elif source == "sra":
        return _download_sra(accession, outdir, threads)
    elif source.startswith("http") or source.startswith("ftp"):
        return _download_url(source, outdir)
    elif source == "ssh":
        return _download_ssh(accession, outdir)
    else:
        error(f"Unknown download source: {source}")
        return False


def _download_ena(accession: str, outdir: Path, threads: int) -> bool:
    """Download from ENA using wget with FTP URLs."""
    # ENA FTP URL pattern
    prefix = accession[:6]
    if len(accession) == 10:
        subdir = f"00{accession[-1]}"
    elif len(accession) == 11:
        subdir = f"0{accession[-2:]}"
    else:
        subdir = ""

    base_url = f"ftp://ftp.sra.ebi.ac.uk/vol1/fastq/{prefix}/{subdir}/{accession}"

    info(f"Downloading from ENA: {accession}")
    ret = run_cmd([
        "wget", "-r", "-nd", "-np",
        "-A", "*.fastq.gz",
        "--cut-dirs=7",
        "-P", str(outdir),
        base_url + "/"
    ])

    if ret != 0:
        # Try aspera or alternative
        warn("wget failed, trying direct file download")
        for suffix in ["_1.fastq.gz", "_2.fastq.gz", ".fastq.gz"]:
            run_cmd(["wget", "-P", str(outdir),
                     f"{base_url}/{accession}{suffix}"])

    return (outdir / f"{accession}_1.fastq.gz").exists() or \
           (outdir / f"{accession}.fastq.gz").exists()


def _download_sra(accession: str, outdir: Path, threads: int) -> bool:
    """Download from SRA using fasterq-dump."""
    info(f"Downloading from SRA: {accession}")
    ret = run_cmd([
        "fasterq-dump", accession,
        "--outdir", str(outdir),
        "--threads", str(threads),
        "--split-files",
        "--progress"
    ])
    if ret == 0:
        # Compress
        run_cmd(["pigz", "-p", str(threads)] +
                [str(f) for f in outdir.glob(f"{accession}*.fastq")])
    return ret == 0


def _download_url(url: str, outdir: Path) -> bool:
    """Download from HTTP/FTP URL."""
    info(f"Downloading: {url}")
    return run_cmd(["wget", "-P", str(outdir), "-c", url]) == 0


def _download_ssh(remote_path: str, outdir: Path) -> bool:
    """Download from SSH remote path (user@host:/path)."""
    info(f"Downloading via SSH: {remote_path}")
    return run_cmd(["rsync", "-avz", "--progress",
                    remote_path, str(outdir)]) == 0


# ─────────────────────────────────────────────────────────────────────────────
# Strandedness detection
# ─────────────────────────────────────────────────────────────────────────────

def detect_strandedness(bam: Path, bed12: Path) -> str:
    """
    Run RSeQC infer_experiment.py and return strandedness.
    Returns: "RF" | "FR" | "unstranded"
    """
    result = subprocess.run(
        ["infer_experiment.py", "-i", str(bam), "-r", str(bed12), "-s", "200000"],
        capture_output=True, text=True
    )
    out = result.stdout + result.stderr

    # Parse fractions
    rf_frac = 0.0
    fr_frac = 0.0
    for line in out.split("\n"):
        if "1++,1--,2+-,2-+" in line or "Fraction of reads failed" not in line:
            try:
                if "1++,1--" in line:
                    rf_frac = float(line.strip().split(":")[-1])
                elif "1+-,1-+" in line:
                    fr_frac = float(line.strip().split(":")[-1])
            except ValueError:
                pass

    if rf_frac > 0.7:
        return "RF"
    elif fr_frac > 0.7:
        return "FR"
    else:
        return "unstranded"


def strandedness_to_featurecounts(strand: str, layout: str) -> str:
    """Convert strandedness to featureCounts -s parameter."""
    if strand == "unstranded":
        return "0"
    elif strand == "RF":
        return "2"  # reverse stranded
    else:
        return "1"  # forward stranded


def strandedness_to_star(strand: str) -> List[str]:
    """Return STAR strandedness-related flags."""
    if strand == "unstranded":
        return []
    else:
        return ["--outSAMstrandField", "intronMotif"]


def strandedness_to_rmats(strand: str, layout: str) -> List[str]:
    """Return rMATS -t flag."""
    t = "paired" if layout == "PE" else "single"
    return ["-t", t]


# ─────────────────────────────────────────────────────────────────────────────
# Core module runners
# ─────────────────────────────────────────────────────────────────────────────

class ModuleRunner:
    """
    Runs individual MOSAIC modules with proper parameters.
    All methods receive the unified config and comparison context.
    """

    def __init__(self, cfg: Any, state: PipelineState,
                 scripts_dir: Path, force: bool = False):
        self.cfg = cfg
        self.state = state
        self.scripts = scripts_dir
        self.force = force
        self.out = Path(cfg.project.wd) / "mosaic_output"

    def _log(self, module: str, comparison: Optional[str] = None) -> Path:
        tag = f"{module}_{comparison}" if comparison else module
        return self.out / "logs" / f"{tag}.log"

    def _comparison_dir(self, comparison: Optional[str]) -> Path:
        if comparison:
            return self.out / f"comparison_{comparison}"
        return self.out

    def _should_run(self, module: str, comparison: Optional[str] = None) -> bool:
        if self.force:
            return True
        return not self.state.is_done(module, comparison)

    def run_module(self, module: str,
                   comparison: Optional[str] = None) -> bool:
        """Run a module, tracking state. Returns True on success."""
        if not self._should_run(module, comparison):
            info(f"[SKIP] {module} already complete")
            self.state.mark(module, "skipped", comparison)
            return True

        section(f"Running {module}" + (f" [{comparison}]" if comparison else ""))
        self.state.mark(module, "running", comparison)

        try:
            method = getattr(self, f"_run_{module}", None)
            if method is None:
                warn(f"No runner for {module} — skipping")
                self.state.mark(module, "skipped", comparison)
                return True

            ok = method(comparison=comparison)
            status = "done" if ok else "failed"
            self.state.mark(module, status, comparison,
                            error=None if ok else f"{module} returned error")
            return ok

        except Exception as e:
            self.state.mark(module, "failed", comparison, error=str(e))
            error(f"{module} failed: {e}")
            return False

    # ── M00: FASTQ download ──────────────────────────────────────────────────

    def _run_m00_download(self, comparison=None) -> bool:
        cfg = self.cfg
        if not hasattr(cfg, "samples") or not hasattr(cfg.samples, "download"):
            info("No download config — assuming FASTQs are local")
            return True

        fastq_dir = Path(cfg.samples.fastq_dir)
        all_ok = True
        for sample in cfg.samples.download:
            sid    = sample["sample_id"]
            source = sample.get("source", "ena")
            acc    = sample.get("accession", sid)
            outdir = fastq_dir / sid
            ok = download_fastq(source, acc, outdir, cfg.threads)
            if not ok:
                warn(f"Download failed for {sid}")
                all_ok = False
        return all_ok

    # ── M02: Trimming ────────────────────────────────────────────────────────

    def _run_m02_trim(self, comparison=None) -> bool:
        cfg = self.cfg
        fastq_dir = Path(cfg.samples.fastq_dir)
        out_trim  = self.out / "trimmed"
        out_trim.mkdir(parents=True, exist_ok=True)
        layout    = cfg.library.layout  # PE | SE

        all_ok = True
        for sid in cfg.samples.sample_ids:
            out_sid = out_trim / sid
            out_sid.mkdir(exist_ok=True)

            if layout == "PE":
                r1 = fastq_dir / sid / f"{sid}_R1.fastq.gz"
                r2 = fastq_dir / sid / f"{sid}_R2.fastq.gz"
                # Try alternative naming
                if not r1.exists():
                    r1 = fastq_dir / sid / f"{sid}_1.fastq.gz"
                    r2 = fastq_dir / sid / f"{sid}_2.fastq.gz"
                cmd = [
                    "fastp",
                    "-i", str(r1), "-I", str(r2),
                    "-o", str(out_sid / f"{sid}_R1_trimmed.fastq.gz"),
                    "-O", str(out_sid / f"{sid}_R2_trimmed.fastq.gz"),
                    "-j", str(out_sid / f"{sid}.fastp.json"),
                    "-h", str(out_sid / f"{sid}.fastp.html"),
                    "--thread", str(min(cfg.threads, 16)),
                    "--detect_adapter_for_pe",
                    "-q", "20", "-l", "36"
                ]
            else:
                r1 = fastq_dir / sid / f"{sid}.fastq.gz"
                cmd = [
                    "fastp",
                    "-i", str(r1),
                    "-o", str(out_sid / f"{sid}_R1_trimmed.fastq.gz"),
                    "-j", str(out_sid / f"{sid}.fastp.json"),
                    "-h", str(out_sid / f"{sid}.fastp.html"),
                    "--thread", str(min(cfg.threads, 16)),
                    "-q", "20", "-l", "36"
                ]

            ret = run_cmd(cmd, log_file=self._log("m02_trim", sid))
            if ret != 0:
                warn(f"Trimming failed for {sid}")
                all_ok = False

        return all_ok

    # ── M04: STAR alignment ──────────────────────────────────────────────────

    def _run_m04_star(self, comparison=None) -> bool:
        cfg = self.cfg
        star_idx  = Path(cfg.databases.star_index)
        trim_dir  = self.out / "trimmed"
        bam_dir   = self.out / "mapping" / "bam"
        layout    = cfg.library.layout
        bam_dir.mkdir(parents=True, exist_ok=True)

        all_ok = True
        for sid in cfg.samples.sample_ids:
            sid_bam = bam_dir / sid
            sid_bam.mkdir(exist_ok=True)

            r1 = trim_dir / sid / f"{sid}_R1_trimmed.fastq.gz"
            reads = [str(r1)]
            if layout == "PE":
                r2 = trim_dir / sid / f"{sid}_R2_trimmed.fastq.gz"
                reads.append(str(r2))

            # Strandedness flags
            strand_flags = strandedness_to_star(
                getattr(cfg, "_strandedness", {}).get(sid, "RF"))

            cmd = [
                "bash", "-c",
                f"ulimit -n 65536 && STAR "
                f"--runMode alignReads "
                f"--genomeDir {star_idx} "
                f"--readFilesIn {' '.join(reads)} "
                f"--outFileNamePrefix {sid_bam}/{sid}_ "
                f"--runThreadN {cfg.threads} "
                f"--limitBAMsortRAM 27000000000 "
                f"--outSAMtype BAM SortedByCoordinate "
                f"--outSAMattributes NH HI AS NM MD "
                f"--outFilterIntronMotifs RemoveNoncanonical "
                f"--twopassMode Basic "
                f"--outFilterMultimapNmax 20 "
                f"--outFilterMismatchNmax 999 "
                f"--outFilterMismatchNoverLmax 0.04 "
                f"--alignIntronMin 20 "
                f"--alignIntronMax 1000000 "
                f"--alignMatesGapMax 1000000 "
                f"--readFilesCommand zcat "
                f"--quantMode GeneCounts "
                f"--chimOutType SeparateSAMold "
                f"--chimSegmentMin 20 "
                f"--chimScoreDropMax 20 "
                f"--chimJunctionOverhangMin 15 "
                f"--chimScoreMin 1 "
                f"--chimScoreSeparation 10 "
                f"--outReadsUnmapped Fastx "
                f"{'--outSAMstrandField intronMotif' if strand_flags else ''} "
                f"2>&1 | tee {self._log('m04_star', sid)}"
            ]

            ret = subprocess.run(" ".join(cmd), shell=True).returncode
            if ret == 0:
                # Index BAM
                run_cmd(["samtools", "index",
                         str(sid_bam / f"{sid}_Aligned.sortedByCoord.out.bam")])
            else:
                warn(f"STAR failed for {sid}")
                all_ok = False

        return all_ok

    # ── M07: Strandedness detection ──────────────────────────────────────────

    def _run_m07_strandedness(self, comparison=None) -> bool:
        cfg     = self.cfg
        bam_dir = self.out / "mapping" / "bam"
        bed12   = Path(cfg.databases.get("gtf_bed12", ""))

        strand_map = {}
        for sid in cfg.samples.sample_ids:
            bam = bam_dir / sid / f"{sid}_Aligned.sortedByCoord.out.bam"
            if bam.exists() and bed12.exists():
                s = detect_strandedness(bam, bed12)
                info(f"  {sid}: {s}")
            else:
                s = cfg.library.strandedness or "RF"
                info(f"  {sid}: {s} (from config, no auto-detect)")
            strand_map[sid] = s

        # Save strandedness
        out_f = self.out / "mapping" / "strandedness.tsv"
        out_f.parent.mkdir(exist_ok=True)
        with open(out_f, "w") as fh:
            fh.write("sample_id\tstrandedness\n")
            for sid, s in strand_map.items():
                fh.write(f"{sid}\t{s}\n")

        # Attach to config for downstream use
        self.cfg._strandedness = strand_map
        return True

    # ── M08: DESeq2 (per comparison) ────────────────────────────────────────

    def _run_m08_deseq2(self, comparison: str = None) -> bool:
        cfg = self.cfg
        comp = self._parse_comparison(comparison)
        if not comp:
            return False

        treatment, control = comp
        out_de = self._comparison_dir(comparison) / "deseq2"
        out_de.mkdir(parents=True, exist_ok=True)

        # Write coldata for this comparison
        coldata_path = self._write_coldata(comparison)

        # Strandedness → featureCounts -s
        strand = self._get_consensus_strandedness()
        s_flag = strandedness_to_featurecounts(strand, cfg.library.layout)

        # featureCounts
        bam_files = self._get_bam_files()
        counts_out = out_de / "all_genes_counts.tsv"

        ret = run_cmd([
            "featureCounts",
            "-T", str(cfg.threads),
            "-a", str(Path(cfg.databases.gtf)),
            "-o", str(counts_out),
            "-t", "exon",
            "-g", "gene_id",
            "-s", s_flag,
        ] + (["-p", "--countReadPairs"] if cfg.library.layout == "PE" else []) + [
            "--fracOverlap", "0.1",
        ] + bam_files,
            log_file=self._log("m08_featurecounts", comparison)
        )
        if ret != 0:
            return False

        # Post-process featureCounts output:
        # 1. Remove comment line (starts with #)
        # 2. Rename BAM path columns to sample_ids
        # 3. Keep only Geneid + count columns (drop Chr/Start/End/Strand/Length)
        self._clean_featurecounts(counts_out)

        # DESeq2
        thresholds = cfg.thresholds
        ret = run_r_script(
            self.scripts / "m08_deseq2.R",
            [
                "--counts",     str(counts_out),
                "--te-counts",  str(self.out / "integrated" / "te_counts_clean.tsv"),
                "--sat-counts", str(self.out / "integrated" / "satellite_counts_clean.tsv"),
                "--coldata",    str(coldata_path),
                "--outdir",     str(out_de),
                "--lfc",        str(thresholds.lfc),
                "--padj",       str(thresholds.padj),
                "--shrink",     "ashr",
                "--mode",       "simple",
                "--ref",        control,
            ],
            log_file=self._log("m08_deseq2", comparison)
        )
        return ret == 0

    # ── M17: rMATS (per comparison) ─────────────────────────────────────────

    def _run_m17_rmats(self, comparison: str = None) -> bool:
        cfg  = self.cfg
        comp = self._parse_comparison(comparison)
        if not comp:
            return False

        treatment, control = comp
        out_dir = self._comparison_dir(comparison) / "m17_rmats"
        out_dir.mkdir(parents=True, exist_ok=True)
        tmp_dir = out_dir / "tmp"
        tmp_dir.mkdir(exist_ok=True)

        # Get BAMs per condition
        normal_bams = self._get_condition_bams(control)
        tumor_bams  = self._get_condition_bams(treatment)

        b1_file = out_dir / "b1_normal.txt"
        b2_file = out_dir / "b2_tumor.txt"
        b1_file.write_text(",".join(normal_bams))
        b2_file.write_text(",".join(tumor_bams))

        # Determine if we have enough replicates for statistics
        use_statoff = len(normal_bams) < 3 or len(tumor_bams) < 3
        if use_statoff:
            warn(f"< 3 replicates per condition — running rMATS with --statoff")

        cmd = [
            "rmats.py",
            "--b1", str(b1_file),
            "--b2", str(b2_file),
            "--gtf", str(Path(cfg.databases.gtf)),
            "--od", str(out_dir),
            "--tmp", str(tmp_dir),
            "-t", "paired" if cfg.library.layout == "PE" else "single",
            "--readLength", str(cfg.sequencing.read_length),
            "--nthread", str(cfg.threads),
            "--novelSS",
        ]
        if use_statoff:
            cmd.append("--statoff")

        ret = run_cmd(cmd, log_file=self._log("m17_rmats", comparison))
        if ret != 0:
            return False

        # Run R analysis
        thresholds = cfg.thresholds
        r_args = [
            "--rmats-dir", str(out_dir),
            "--outdir",    str(out_dir),
            "--treatment", treatment,
            "--control",   control,
            "--dpsi",      "0.1",
            "--min-reads", str(thresholds.min_reads_splicing),
        ]
        if not use_statoff:
            r_args += ["--fdr", str(thresholds.padj)]
        else:
            r_args.append("--statoff")

        ret = run_r_script(
            self.scripts / "m17_rmats.R",
            r_args,
            log_file=self._log("m17_rmats_r", comparison)
        )
        return ret == 0

    # ── M14: Report (per comparison) ────────────────────────────────────────

    def _run_m14_report(self, comparison: str = None) -> bool:
        cfg  = self.cfg
        comp = self._parse_comparison(comparison)
        if not comp:
            return False

        treatment, control = comp
        out_dir = self._comparison_dir(comparison) / "report"
        out_dir.mkdir(parents=True, exist_ok=True)

        # Determine out_dir for report (comparison-specific)
        comp_out = self._comparison_dir(comparison)

        thresholds = cfg.thresholds
        # out_dir = main output dir where all module results live
        # report output goes into comparison subdir
        main_out = self.out
        ret = run_cmd([
            "quarto", "render",
            str(self.scripts / "m14_report.qmd"),
            "-P", f"out_dir:{main_out}",
            "-P", f"comparison:{comparison}",
            "-P", f"treatment:{treatment}",
            "-P", f"control:{control}",
            "-P", f"padj_threshold:{thresholds.padj}",
            "-P", f"lfc_threshold:{thresholds.lfc}",
            "--output-dir", str(out_dir),
        ], log_file=self._log("m14_report", comparison))
        return ret == 0

    # ── Helpers ──────────────────────────────────────────────────────────────

    def _clean_featurecounts(self, counts_file: Path) -> None:
        """
        Clean featureCounts output:
        - Remove leading comment line (# Program:...)
        - Replace full BAM paths with sample_ids
        - Drop annotation columns (Chr, Start, End, Strand, Length)
        """
        lines = counts_file.read_text().splitlines()
        # Skip comment lines
        data_lines = [l for l in lines if not l.startswith("#")]
        if not data_lines:
            return

        # Parse header — replace BAM paths with sample_ids
        header = data_lines[0].split("\t")
        # header = [Geneid, Chr, Start, End, Strand, Length, /path/bam1, ...]
        # Find count column indices (BAM paths)
        annot_cols = {"Geneid", "Chr", "Start", "End", "Strand", "Length"}
        new_header = []
        keep_idx = []
        for i, col in enumerate(header):
            if col in annot_cols:
                if col == "Geneid":
                    new_header.append(col)
                    keep_idx.append(i)
                # skip other annotation cols
            else:
                # BAM path -> extract sample_id from path
                bam_path = Path(col)
                sid = bam_path.stem.replace("_Aligned.sortedByCoord.out", "")
                new_header.append(sid)
                keep_idx.append(i)

        # Rebuild file with cleaned header and only kept columns
        out_lines = ["\t".join(new_header)]
        for line in data_lines[1:]:
            parts = line.split("\t")
            out_lines.append("\t".join(parts[i] for i in keep_idx))

        counts_file.write_text("\n".join(out_lines) + "\n")


    def _parse_comparison(self, comparison: Optional[str]) -> Optional[Tuple[str, str]]:
        if not comparison:
            # Use first comparison from config
            if self.cfg.comparisons:
                c = self.cfg.comparisons[0]
                return c["treatment"], c["control"]
            return None
        parts = comparison.split("_vs_")
        if len(parts) == 2:
            return parts[0], parts[1]
        return None

    def _get_bam_files(self) -> List[str]:
        bam_dir = self.out / "mapping" / "bam"
        bams = []
        for sid in self.cfg.samples.sample_ids:
            bam = bam_dir / sid / f"{sid}_Aligned.sortedByCoord.out.bam"
            if bam.exists():
                bams.append(str(bam))
        return bams

    def _get_condition_bams(self, condition: str) -> List[str]:
        """Get BAM files for a specific condition."""
        coldata = self.cfg.samples.coldata_df or []
        sids = [r["sample_id"] for r in coldata
                if r.get("condition") == condition]
        bam_dir = self.out / "mapping" / "bam"
        return [str(bam_dir / sid / f"{sid}_Aligned.sortedByCoord.out.bam")
                for sid in sids]

    def _get_consensus_strandedness(self) -> str:
        strand_map = getattr(self.cfg, "_strandedness", {})
        if strand_map:
            from collections import Counter
            return Counter(strand_map.values()).most_common(1)[0][0]
        # Try strandedness from config
        strand = (getattr(self.cfg.library, "strandedness", None) or
                  getattr(self.cfg.library, "_strandedness_config", None) or
                  getattr(self.cfg, "strandedness", None) or "RF")
        if strand == "auto":
            strand = "RF"
        return strand

    def _write_coldata(self, comparison: str) -> Path:
        """Write coldata TSV for a specific comparison."""
        comp = self._parse_comparison(comparison)
        if not comp:
            return Path(self.cfg.samples.coldata)

        treatment, control = comp
        coldata_rows = self.cfg.samples.coldata_df  # list of dicts
        if not coldata_rows:
            # Fallback to original coldata file
            return Path(self.cfg.samples.coldata)

        subset = [r for r in coldata_rows
                  if r.get("condition") in (treatment, control)]

        out_path = self.out / f"coldata_{comparison}.tsv"
        out_path.parent.mkdir(parents=True, exist_ok=True)
        if subset:
            keys = list(subset[0].keys())
            with open(out_path, "w") as fh:
                fh.write("\t".join(keys) + "\n")
                for row in subset:
                    fh.write("\t".join(str(row.get(k,"")) for k in keys) + "\n")
        return out_path


# ─────────────────────────────────────────────────────────────────────────────
# Pipeline orchestrator
# ─────────────────────────────────────────────────────────────────────────────

class Pipeline:
    """
    Main pipeline orchestrator.
    Resolves dependencies, runs modules in order,
    handles multiple comparisons.
    """

    def __init__(self, cfg: Any, scripts_dir: Path,
                 force: bool = False,
                 from_module: Optional[str] = None):
        self.cfg      = cfg
        self.scripts  = scripts_dir
        self.force    = force
        self.from_mod = from_module

        outdir = Path(cfg.project.wd) / "mosaic_output"
        outdir.mkdir(parents=True, exist_ok=True)

        state_file = outdir / ".mosaic_state.json"
        self.state  = PipelineState(state_file)
        self.runner = ModuleRunner(cfg, self.state, scripts_dir, force)

    def _resolve_order(self, modules: List[str]) -> List[str]:
        """Topological sort of modules by dependencies."""
        visited = set()
        order   = []

        def visit(m):
            if m in visited:
                return
            visited.add(m)
            for dep in MODULE_DEPS.get(m, []):
                if dep in modules:
                    visit(dep)
            order.append(m)

        for m in modules:
            visit(m)
        return order

    def run(self, modules: Optional[List[str]] = None) -> bool:
        """Run the full pipeline or specified modules."""
        # Determine which modules to run
        enabled = []
        if modules:
            enabled = modules
        else:
            mods_cfg = self.cfg.modules
            enabled  = [m for m, v in vars(mods_cfg).items() if v]

        # If --from-module specified, mark earlier modules as done
        if self.from_mod and not self.force:
            all_mods = self._resolve_order(list(MODULE_DEPS.keys()))
            idx = all_mods.index(self.from_mod) if self.from_mod in all_mods else 0
            for m in all_mods[:idx]:
                if not self.state.is_done(m):
                    self.state.mark(m, "skipped")

        # Resolve execution order
        order = self._resolve_order(enabled)

        # Get comparisons
        comparisons = [
            f"{c['treatment']}_vs_{c['control']}"
            for c in getattr(self.cfg, "comparisons", [])
        ]
        if not comparisons:
            comparisons = [None]  # single comparison from config

        # Comparison-independent modules
        COMP_INDEPENDENT = {
            "m00_download", "m02_trim", "m04_star",
            "m05_te_count", "m07_strandedness",
            "m09_satellite", "m11_te_host",
            "m11b_promoter_te", "m11c_te_age",
            "m19_editing", "m20_meta", "m20b_humann",
            "m23_te_insertions", "m12_tetranscripts",
            "m15_fusion", "m16_stringtie",
        }

        section("MOSAIC Pipeline")
        info(f"Modules: {', '.join(order)}")
        info(f"Comparisons: {', '.join(str(c) for c in comparisons)}")
        info(f"Force re-run: {self.force}")

        all_ok = True
        start  = datetime.now()

        for module in order:
            if module in COMP_INDEPENDENT:
                ok = self.runner.run_module(module)
            else:
                # Run once per comparison
                ok = True
                for comp in comparisons:
                    ok = ok and self.runner.run_module(module, comp)

            if not ok:
                error(f"Module {module} failed — stopping pipeline")
                all_ok = False
                break

        elapsed = (datetime.now() - start).total_seconds()
        if all_ok:
            success(f"Pipeline complete in {elapsed/60:.1f} min")
        else:
            error(f"Pipeline failed after {elapsed/60:.1f} min")

        return all_ok

    def status(self) -> None:
        """Print current pipeline status."""
        section("Pipeline Status")
        for key, m in sorted(self.state.modules.items()):
            icon = {"done": "✓", "failed": "✗", "running": "▶",
                    "skipped": "→", "pending": "○"}.get(m.status, "?")
            elapsed = f" ({m.elapsed:.0f}s)" if m.elapsed else ""
            print(f"  {icon} {key:<40} {m.status}{elapsed}")
