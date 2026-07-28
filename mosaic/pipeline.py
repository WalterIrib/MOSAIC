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
    "m08b_qc":         ["m08_deseq2"],
    "m09_satellite":   ["m02_trim"],
    "m10_go":          ["m08_deseq2"],
    "m11_te_host":     ["m04_star"],
    "m11b_promoter_te":["m11_te_host"],
    "m11c_te_age":     ["m11_te_host", "m08_deseq2"],
    "m11d_te_de":      ["m08_deseq2", "m11b_promoter_te"],
    "m12_tetranscripts":["m04_star", "m11_te_host", "m08_deseq2"],
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
            log_file: Optional[Path] = None,
            extra_env: Optional[dict] = None) -> int:
    """
    Run a command, optionally in a specific conda env.
    Returns exit code.
    """
    if env_name:
        cmd = ["conda", "run", "-n", env_name, "--no-capture-output"] + cmd

    cmd_str = " ".join(str(c) for c in cmd)
    info(f"$ {cmd_str}")

    import os
    run_env = os.environ.copy()
    if extra_env:
        run_env.update(extra_env)

    kwargs: dict = dict(cwd=cwd, env=run_env)
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
    import shutil
    if not shutil.which("fasterq-dump"):
        warn("fasterq-dump not found. Install SRA-tools:")
        warn("  conda install -c bioconda sra-tools -y")
        warn("  OR: mosaic download --source ena instead")
        return False

    info(f"Downloading from SRA: {accession}")
    ret = run_cmd([
        "fasterq-dump", accession,
        "--outdir", str(outdir),
        "--threads", str(threads),
        "--split-files",
        "--progress"
    ])
    if ret == 0:
        # Compress with pigz or gzip
        compressor = "pigz" if shutil.which("pigz") else "gzip"
        for f in outdir.glob(f"{accession}*.fastq"):
            run_cmd([compressor] + (["-p", str(threads)] if compressor == "pigz" else []) + [str(f)])
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
    Run RSeQC infer_experiment.py and return strandedness call.
    Returns: "RF" | "FR" | "unstranded"
    """
    stats = detect_strandedness_full(bam, bed12)
    return stats["strandedness"]


def detect_strandedness_full(bam: Path, bed12: Path) -> dict:
    """
    Run RSeQC infer_experiment.py and return full stats dict:
      frac_rf, frac_fr, frac_failed, strandedness, strand_fc
    """
    import re as _re
    result = subprocess.run(
        ["infer_experiment.py", "-i", str(bam), "-r", str(bed12), "-s", "200000"],
        capture_output=True, text=True
    )
    out = result.stdout + result.stderr

    frac_rf = frac_fr = frac_fail = 0.0
    for line in out.splitlines():
        m = _re.search(r"([\d.]+)\s*$", line)
        if not m:
            continue
        val = float(m.group(1))
        if "1++,1--,2+-,2-+" in line:
            frac_rf = val
        elif "1+-,1-+,2++,2--" in line:
            frac_fr = val
        elif "failed to determine" in line:
            frac_fail = val

    # frac_rf = "1++,1--,2+-,2-+" = R1 on same strand as gene = standard FR-forward
    # frac_fr = "1+-,1-+,2++,2--" = R2 on same strand as gene = standard RF-dUTP
    # featureCounts: -s 1 = R1 same as feature (fr-secondstrand)
    #                -s 2 = R2 same as feature (fr-firststrand / dUTP / TruSeq Stranded)
    if frac_rf > 0.7:
        strand = "RF"
    elif frac_fr > 0.7:
        strand = "FR"
    else:
        strand = "unstranded"

    return {
        "frac_rf":    frac_rf,
        "frac_fr":    frac_fr,
        "frac_failed": frac_fail,
        "strandedness": strand,
        "strand_fc":  {"unstranded": "0", "RF": "1", "FR": "2"}[strand],
    }


def parse_star_log(log_path: Path) -> dict:
    """Parse STAR Log.final.out and return key mapping stats."""
    import re as _re
    stats = {}
    if not log_path.exists():
        return stats
    txt = log_path.read_text()

    def get(pattern):
        m = _re.search(pattern, txt)
        return m.group(1).strip() if m else "NA"

    stats["n_input_reads"]       = get(r"Number of input reads \|\s+(\d+)")
    stats["pct_uniquely_mapped"] = get(r"Uniquely mapped reads % \|\s+([\d.]+%)")
    stats["pct_multimapped"]     = get(r"% of reads mapped to multiple loci \|\s+([\d.]+%)")
    return stats


def _find_or_build_bed12(cfg: Any, out_dir: Path) -> Optional[Path]:
    """
    Return path to BED12 for RSeQC.
    Priority:
      1. cfg.databases.gtf_bed12 (explicit config key)
      2. known output location from previous run
      3. Build from GTF using gtfToGenePred / bedtools / Python fallback
    """
    import shutil as _shutil

    # 1. Explicit config key
    explicit = cfg.databases.get("gtf_bed12", "")
    if explicit and Path(explicit).exists():
        return Path(explicit)

    # 2. Known output location (generated by previous run)
    candidate = out_dir / "mapping" / "annotation" / "genes.bed12"
    if candidate.exists():
        return candidate

    # 3. Build from GTF
    gtf = Path(cfg.databases.gtf)
    if not gtf.exists():
        warn("GTF not found — cannot build BED12 for RSeQC")
        return None

    annot_dir = out_dir / "mapping" / "annotation"
    annot_dir.mkdir(parents=True, exist_ok=True)
    bed12 = annot_dir / "genes.bed12"

    if _shutil.which("gtfToGenePred") and _shutil.which("genePredToBed"):
        gp = annot_dir / "genes.genePred"
        r1 = subprocess.run(
            ["gtfToGenePred", "-genePredExt", str(gtf), str(gp)],
            capture_output=True)
        r2 = subprocess.run(
            ["genePredToBed", str(gp), str(bed12)],
            capture_output=True)
        if r1.returncode == 0 and r2.returncode == 0 and bed12.exists():
            info(f"Built BED12 via gtfToGenePred: {bed12}")
            return bed12

    # Python fallback: parse GTF → BED12
    info("Building BED12 from GTF (Python parser)...")
    try:
        _gtf_to_bed12(gtf, bed12)
        info(f"Built BED12: {bed12} ({bed12.stat().st_size // 1024}KB)")
        return bed12
    except Exception as e:
        warn(f"BED12 build failed: {e}")
        return None


def _gtf_to_bed12(gtf_path: Path, bed12_path: Path) -> None:
    """Convert GTF to BED12 format (transcript-level, for RSeQC)."""
    from collections import defaultdict

    # transcript -> {chrom, strand, txStart, txEnd, exons: [(start,end)], name}
    transcripts: dict = defaultdict(lambda: {
        "chrom": "", "strand": ".", "txStart": int(1e18),
        "txEnd": 0, "name": "", "exons": []
    })

    with open(gtf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            feat = parts[2]
            if feat not in ("transcript", "exon"):
                continue
            chrom  = parts[0]
            start  = int(parts[3]) - 1   # 0-based
            end    = int(parts[4])
            strand = parts[6]
            attrs  = parts[8]

            # Extract transcript_id
            m = __import__("re").search(r'transcript_id "([^"]+)"', attrs)
            if not m:
                continue
            tid = m.group(1)

            t = transcripts[tid]
            t["chrom"]  = chrom
            t["strand"] = strand
            if feat == "transcript":
                m2 = __import__("re").search(r'gene_name "([^"]+)"', attrs)
                t["name"] = m2.group(1) if m2 else tid
                t["txStart"] = min(t["txStart"], start)
                t["txEnd"]   = max(t["txEnd"], end)
            elif feat == "exon":
                t["exons"].append((start, end))
                t["txStart"] = min(t["txStart"], start)
                t["txEnd"]   = max(t["txEnd"], end)

    with open(bed12_path, "w") as out:
        for tid, t in transcripts.items():
            if not t["exons"] or not t["chrom"]:
                continue
            exons = sorted(t["exons"])
            tx_start = t["txStart"]
            tx_end   = t["txEnd"]
            block_count  = len(exons)
            block_sizes  = ",".join(str(e - s) for s, e in exons) + ","
            block_starts = ",".join(str(s - tx_start) for s, e in exons) + ","
            out.write(
                f"{t['chrom']}\t{tx_start}\t{tx_end}\t{t['name']}\t0\t{t['strand']}\t"
                f"{tx_start}\t{tx_end}\t0\t{block_count}\t{block_sizes}\t{block_starts}\n"
            )


def strandedness_to_featurecounts(strand: str, layout: str) -> str:
    """Convert strandedness to featureCounts -s parameter.
    RF = "1++,1--,2+-,2-+" = R1 same as gene = fr-secondstrand = -s 1
    FR = "1+-,1-+,2++,2--" = R2 same as gene = fr-firststrand/dUTP/TruSeq = -s 2
    """
    if strand == "unstranded":
        return "0"
    elif strand == "RF":
        return "1"  # R1 on sense strand (fr-secondstrand)
    else:  # FR
        return "2"  # R2 on sense strand (fr-firststrand / dUTP)


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
        key = f"{module}/{comparison}" if comparison else module
        m = self.state.modules.get(key) or self.state.modules.get(module)
        # --force overrides everything except modules skipped by --from-module
        # (those are marked with a special flag)
        if m and m.status == "skipped" and getattr(m, "_from_module_skip", False):
            return False  # explicitly skipped by --from-module, never re-run
        if self.force:
            return True
        if m and m.status == "skipped":
            return True   # skipped previously but not by from_module → run now
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
            # Interaction comparisons: modules that can't run (no treatment/control)
            # return False but should be treated as skip, not failure
            if not ok and self._is_interaction_comparison(comparison):
                warn(f"[SKIP] {module} not applicable to interaction comparison")
                self.state.mark(module, "skipped", comparison)
                return True
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

    # ── M05: TE/Satellite fragment counting ──────────────────────────────────

    def _run_m05_te_count(self, comparison=None) -> bool:
        from .m05_te_count import run_m05
        cfg      = self.cfg
        trim_dir = self.out / "trimmed"
        out_dir  = self.out / "te_counts"
        out_dir.mkdir(parents=True, exist_ok=True)

        te_fa  = Path(cfg.databases.get("te_consensus",
                      str(Path(cfg.databases.te_satellite).parent / "consensus_TE.fa")))
        sat_fa = Path(cfg.databases.te_satellite)

        # Build bowtie2 index names (stem without extension)
        te_idx  = str(te_fa.parent / te_fa.stem)
        sat_idx = str(sat_fa.parent / sat_fa.stem)

        # Build TE index if missing
        import shutil as _shutil
        if not Path(f"{te_idx}.1.bt2").exists() and not Path(f"{te_idx}.1.bt2l").exists():
            info("Building bowtie2 index for TE consensus...")
            run_cmd(["bowtie2-build", str(te_fa), te_idx,
                     "--threads", str(min(cfg.threads, 8))])

        # Read trim_stats to get per-sample FASTQ paths
        trim_stats_f = trim_dir / "trim_stats.tsv"
        if not trim_stats_f.exists():
            warn("trim_stats.tsv not found — run m02_trim first")
            return False
        import csv as _csv
        with open(trim_stats_f) as fh:
            trim_stats = list(_csv.DictReader(fh, delimiter="\t"))

        # Filter to samples that have trimmed FASTQs
        valid = [s for s in trim_stats if Path(s.get("r1_trimmed","")).exists()]
        if not valid:
            warn("No trimmed FASTQs found — run m02_trim first")
            return False

        class_tsv = Path(cfg.databases.get("te_classification",
                         str(te_fa.parent / "classification.tsv")))

        run_m05(
            trim_stats   = valid,
            te_index     = te_idx,
            satellite_index = sat_idx,
            classification_tsv = str(class_tsv) if class_tsv.exists() else "",
            out_dir      = out_dir,
            threads      = cfg.threads,
        )

        # Copy matrices to integrated/, filtering non-TE artefacts from te_counts.
        # 7SLRNA is an SRP RNA that DFAM labels SINE/Alu — it is NOT a TE and
        # produces massive spurious read counts that dominate TE expression tables.
        _NON_TE_NAMES = {"7SLRNA"}

        integrated = self.out / "integrated"
        integrated.mkdir(exist_ok=True)
        import shutil as _sh

        # Satellite counts — copy as-is
        sat_src = out_dir / "satellite_counts.tsv"
        if sat_src.exists():
            _sh.copy(sat_src, integrated / "satellite_counts_clean.tsv")

        # TE counts — filter out non-TE entries by row name (first column)
        te_src = out_dir / "te_counts.tsv"
        te_dst = integrated / "te_counts_clean.tsv"
        if te_src.exists():
            excluded = []
            with open(te_src) as fh_in, open(te_dst, "w") as fh_out:
                for line in fh_in:
                    name = line.split("\t", 1)[0]
                    if name in _NON_TE_NAMES:
                        excluded.append(name)
                    else:
                        fh_out.write(line)
            if excluded:
                info(f"Excluded non-TE entries from te_counts_clean: {sorted(excluded)}")

        return True

    # ── M06: circRNA detection (DCC) ────────────────────────────────────────

    def _run_m06_circrna(self, comparison=None) -> bool:
        from .m06_circrna import run_m06
        cfg     = self.cfg
        bam_dir = self.out / "mapping" / "bam"

        dcc_bin = cfg.databases.dcc_path or "DCC"
        import shutil as _sh
        if not Path(dcc_bin).exists() and not _sh.which(dcc_bin):
            warn("DCC not found — skipping m06_circrna")
            return True   # non-fatal skip

        chimeric_junctions = {}
        for sid in cfg.samples.sample_ids:
            jxn = bam_dir / sid / f"{sid}_Chimeric.out.junction"
            if jxn.exists():
                chimeric_junctions[sid] = jxn

        if not chimeric_junctions:
            warn("No Chimeric.out.junction files found — run m04_star first")
            return False

        comp = self._parse_comparison(comparison)
        treatment, control = comp if comp else (None, None)

        # Global fallback strandedness (used only when M07 hasn't run for a sample)
        global_strand = "unstranded"
        if hasattr(cfg, "sequencing") and hasattr(cfg.sequencing, "strandedness"):
            global_strand = cfg.sequencing.strandedness or "unstranded"

        # Per-sample strandedness from M07 strandedness.tsv (preferred)
        import csv as _csv
        strand_m07: dict[str, str] = {}
        strand_tsv = self.out / "mapping" / "strandedness.tsv"
        if strand_tsv.exists():
            with open(strand_tsv) as fh:
                for row in _csv.DictReader(fh, delimiter="\t"):
                    strand_m07[row["sample_id"]] = row["strandedness"].lower()
            info(f"M06: loaded per-sample strandedness from {strand_tsv.name} "
                 f"({len(strand_m07)} samples)")
        else:
            warn(f"M06: M07 strandedness.tsv not found — using global '{global_strand}' for all samples")

        # Patch cfg so run_m06 can access what it needs
        cfg_m06                    = type("C", (), {})()
        cfg_m06.output_dir         = str(self.out)
        cfg_m06.gtf                = str(Path(cfg.databases.gtf))
        cfg_m06.dcc_path           = dcc_bin
        cfg_m06.threads            = cfg.threads
        cfg_m06.min_junction_reads = 2
        cfg_m06.treatment          = treatment or "tumor"
        cfg_m06.control            = control  or "normal"
        cfg_m06.strandedness       = global_strand  # fallback only

        coldata_f = Path(cfg.samples.coldata)
        with open(coldata_f) as fh:
            samples = list(_csv.DictReader(fh, delimiter="\t"))

        # Inject per-sample strandedness from M07 (overrides global fallback in run_m06)
        for s in samples:
            if s["sample_id"] in strand_m07:
                s["strandedness"] = strand_m07[s["sample_id"]]
        cfg_m06.samples = samples

        result = run_m06(cfg_m06, chimeric_junctions)
        return bool(result)

    # ── M02: Trimming ────────────────────────────────────────────────────────

    def _run_m02_trim(self, comparison=None) -> bool:
        cfg = self.cfg
        fastq_dir = Path(cfg.samples.fastq_dir)
        out_trim  = self.out / "trimmed"
        out_trim.mkdir(parents=True, exist_ok=True)
        layout    = cfg.library.layout  # PE | SE

        # Build R1/R2 lookup from coldata when available (flat FASTQ dirs)
        coldata_fq = {}
        for row in (cfg.samples.coldata_df or []):
            sid = row.get("sample_id", "")
            if row.get("R1") and row.get("R2"):
                coldata_fq[sid] = (Path(row["R1"]), Path(row["R2"]))

        all_ok = True
        for sid in cfg.samples.sample_ids:
            out_sid = out_trim / sid
            out_sid.mkdir(exist_ok=True)

            if layout == "PE":
                if sid in coldata_fq and coldata_fq[sid][0].exists():
                    r1, r2 = coldata_fq[sid]
                else:
                    r1 = fastq_dir / sid / f"{sid}_R1.fastq.gz"
                    r2 = fastq_dir / sid / f"{sid}_R2.fastq.gz"
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
                if sid in coldata_fq and coldata_fq[sid][0].exists():
                    r1 = coldata_fq[sid][0]
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

        # Write trim_stats.tsv from fastp JSON outputs
        self._write_trim_stats(out_trim, cfg)
        return all_ok

    def _write_trim_stats(self, trim_dir: Path, cfg) -> None:
        """Parse fastp JSON files and write trim_stats.tsv."""
        import json as _json
        cond_map = {r["sample_id"]: r.get("condition", "")
                    for r in (cfg.samples.coldata_df or [])}
        rows = []
        cols = ["sample_id", "condition", "reads_before", "reads_after",
                "pass_rate", "mean_len_r1", "q30_rate", "duplication",
                "r1_trimmed", "r2_trimmed"]
        for sid in cfg.samples.sample_ids:
            jf = trim_dir / sid / f"{sid}.fastp.json"
            if not jf.exists():
                continue
            try:
                d = _json.loads(jf.read_text())
                s = d["summary"]
                rb = s["before_filtering"]["total_reads"]
                ra = s["after_filtering"]["total_reads"]
                rows.append({
                    "sample_id":    sid,
                    "condition":    cond_map.get(sid, ""),
                    "reads_before": rb,
                    "reads_after":  ra,
                    "pass_rate":    round(ra / rb * 100, 3) if rb else 0,
                    "mean_len_r1":  s["after_filtering"]["read1_mean_length"],
                    "q30_rate":     round(s["after_filtering"]["q30_rate"] * 100, 4),
                    "duplication":  round(d["duplication"]["rate"] * 100, 4),
                    "r1_trimmed":   str(trim_dir / sid / f"{sid}_R1_trimmed.fastq.gz"),
                    "r2_trimmed":   str(trim_dir / sid / f"{sid}_R2_trimmed.fastq.gz"),
                })
            except Exception as e:
                warn(f"Could not parse fastp JSON for {sid}: {e}")
        if rows:
            out = trim_dir / "trim_stats.tsv"
            with open(out, "w") as fh:
                fh.write("\t".join(cols) + "\n")
                for r in rows:
                    fh.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")
            info(f"trim_stats.tsv written: {len(rows)} samples")

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
                f"--chimOutType SeparateSAMold Junctions "
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

        # Build gene_type_map.tsv from GTF (needed by downstream modules)
        gene_map = self.out / "integrated" / "gene_type_map.tsv"
        if not gene_map.exists():
            try:
                from .refs import extract_gene_type_map
                (self.out / "integrated").mkdir(parents=True, exist_ok=True)
                extract_gene_type_map(Path(cfg.databases.gtf), gene_map)
                info(f"gene_type_map.tsv written: {gene_map}")
            except Exception as e:
                warn(f"Could not build gene_type_map.tsv: {e}")

        return all_ok

    # ── M07: Strandedness detection ──────────────────────────────────────────

    def _run_m07_strandedness(self, comparison=None) -> bool:
        cfg     = self.cfg
        bam_dir = self.out / "mapping" / "bam"

        # Find or build BED12 for RSeQC
        bed12 = _find_or_build_bed12(cfg, self.out)

        cols = ["sample_id", "n_input_reads", "pct_uniquely_mapped",
                "pct_multimapped", "frac_rf", "frac_fr", "frac_failed",
                "strandedness", "strand_fc"]
        rows = []
        strand_map = {}

        for sid in cfg.samples.sample_ids:
            bam = bam_dir / sid / f"{sid}_Aligned.sortedByCoord.out.bam"
            log = bam_dir / sid / f"{sid}_Log.final.out"

            star = parse_star_log(log)

            if bam.exists() and bed12 is not None:
                st = detect_strandedness_full(bam, bed12)
                info(f"  {sid}: {st['strandedness']} "
                     f"(RF={st['frac_rf']:.3f} FR={st['frac_fr']:.3f})")
            else:
                cfg_strand = getattr(cfg.library, "strandedness", None) or "unstranded"
                if cfg_strand == "auto":
                    cfg_strand = "unstranded"
                st = {
                    "frac_rf": 0.0, "frac_fr": 0.0, "frac_failed": 0.0,
                    "strandedness": cfg_strand,
                    "strand_fc": {"unstranded": "0", "RF": "1", "FR": "2"}.get(cfg_strand, "0"),
                }
                info(f"  {sid}: {cfg_strand} (from config, no auto-detect)")

            strand_map[sid] = st["strandedness"]
            rows.append({
                "sample_id":           sid,
                "n_input_reads":       star.get("n_input_reads", "NA"),
                "pct_uniquely_mapped": star.get("pct_uniquely_mapped", "NA"),
                "pct_multimapped":     star.get("pct_multimapped", "NA"),
                "frac_rf":             f"{st['frac_rf']:.4f}",
                "frac_fr":             f"{st['frac_fr']:.4f}",
                "frac_failed":         f"{st['frac_failed']:.4f}",
                "strandedness":        st["strandedness"],
                "strand_fc":           st["strand_fc"],
            })

        # Write comprehensive strandedness.tsv
        out_f = self.out / "mapping" / "strandedness.tsv"
        out_f.parent.mkdir(exist_ok=True)
        with open(out_f, "w") as fh:
            fh.write("\t".join(cols) + "\n")
            for r in rows:
                fh.write("\t".join(r[c] for c in cols) + "\n")

        # Attach to config for downstream use
        self.cfg._strandedness = strand_map
        return True

    # ── M08: DESeq2 (per comparison) ────────────────────────────────────────

    def _run_m08_deseq2(self, comparison: str = None) -> bool:
        cfg = self.cfg
        comp_cfg = self._get_comparison_config(comparison)
        is_interaction = comp_cfg and comp_cfg.get("type") == "interaction"

        # For non-interaction: require treatment/control
        if not is_interaction:
            comp = self._parse_comparison(comparison)
            if not comp:
                return False
            treatment, control = comp

        out_de = self._comparison_dir(comparison) / "deseq2"
        out_de.mkdir(parents=True, exist_ok=True)

        coldata_path = self._write_coldata(comparison)

        # Strandedness → featureCounts -s
        strand = self._get_consensus_strandedness()
        s_flag = strandedness_to_featurecounts(strand, cfg.library.layout)

        # featureCounts (all samples for both simple and interaction)
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

        self._clean_featurecounts(counts_out)

        thresholds = cfg.thresholds

        if is_interaction:
            deseq2_args = [
                "--counts",        str(counts_out),
                "--te-counts",     str(self.out / "integrated" / "te_counts_clean.tsv"),
                "--sat-counts",    str(self.out / "integrated" / "satellite_counts_clean.tsv"),
                "--coldata",       str(coldata_path),
                "--outdir",        str(out_de),
                "--lfc",           str(thresholds.lfc),
                "--padj",          str(thresholds.padj),
                "--shrink",        "ashr",
                "--mode",          "interaction",
                "--fraction-col",  "fraction",
                "--phase-col",     "phase",
                "--ref-fraction",  comp_cfg.get("ref_fraction", "totRNA"),
                "--ref-phase",     comp_cfg.get("ref_phase", "interphase"),
                "--comparison",    comparison,
            ]
        else:
            deseq2_args = [
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
            ]

        ret = run_r_script(
            self.scripts / "m08_deseq2.R",
            deseq2_args,
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
        # Always wipe tmp: reusing prep-step .rmats files from a previous run with
        # different --readLength causes "post step read length mismatch" and 0 events.
        import shutil as _sh_rmats
        if tmp_dir.exists():
            _sh_rmats.rmtree(tmp_dir)
        tmp_dir.mkdir()

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

        # Map MOSAIC strandedness to rMATS --libType
        # FR = dUTP/TruSeq Stranded (R2 same strand as gene) = fr-firststrand
        # RF = forward-stranded (R1 same strand as gene)     = fr-secondstrand
        _strand_cfg = getattr(cfg.library, "_strandedness_config", "unstranded")
        _rmats_libtype = {
            "FR": "fr-firststrand",
            "RF": "fr-secondstrand",
        }.get(_strand_cfg, "fr-unstranded")

        cmd = [
            "rmats.py",
            "--b1", str(b1_file),
            "--b2", str(b2_file),
            "--gtf", str(Path(cfg.databases.gtf)),
            "--od", str(out_dir),
            "--tmp", str(tmp_dir),
            "-t", "paired" if cfg.library.layout == "PE" else "single",
            "--readLength", str(cfg.library.read_length),
            "--variable-read-length",
            "--libType", _rmats_libtype,
            "--nthread", str(cfg.threads),
            "--novelSS",
        ]
        if use_statoff:
            cmd.append("--statoff")

        ret = run_cmd(cmd, log_file=self._log("m17_rmats", comparison))
        # rMATS may return non-zero even on partial success — check output exists
        rmats_out = out_dir / "SE.MATS.JCEC.txt"
        if ret != 0 and not rmats_out.exists():
            warn(f"rMATS failed (exit {ret}) and no output found")
            return False
        if ret != 0:
            warn(f"rMATS returned exit {ret} but output exists — continuing")

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
        thresholds = cfg.thresholds
        out_dir = self._comparison_dir(comparison) / "report"
        out_dir.mkdir(parents=True, exist_ok=True)

        # Interaction comparison → unified template with comparison_type=interaction
        if self._is_interaction_comparison(comparison):
            comp_cfg = next(
                (c for c in cfg.comparisons if c.get("name") == comparison), {}
            ) if hasattr(cfg, "comparisons") else {}
            # Derive meaningful treatment/control labels for axis labels
            # treatment = chromatin fraction in mitosis; control = reference
            frac_of_interest = comp_cfg.get("fraction_of_interest", "ChaRNA")
            treatment_label  = f"mito-{frac_of_interest}"
            control_label    = frac_of_interest
            phase_file = str(
                self._comparison_dir(comparison) / "deseq2" /
                f"genes_{comparison}_phase_effect_DE.tsv"
            )
            ret = run_cmd([
                "quarto", "render",
                "m14_report.qmd",
                "-P", f"out_dir:{self.out}",
                "-P", f"comparison:{comparison}",
                "-P", f"treatment:{treatment_label}",
                "-P", f"control:{control_label}",
                "-P", f"padj_threshold:{thresholds.padj}",
                "-P", f"lfc_threshold:{thresholds.lfc}",
                "-P", f"coldata_path:{cfg.samples.coldata}",
                "-P", "comparison_type:interaction",
                "-P", f"phase_effect_file:{phase_file}",
                "--output-dir", str(out_dir),
            ], cwd=self.scripts, log_file=self._log("m14_report", comparison))
            return ret == 0

        comp = self._parse_comparison(comparison)
        if not comp:
            return False

        treatment, control = comp
        main_out = self.out
        ret = run_cmd([
            "quarto", "render",
            "m14_report.qmd",   # relative path — CWD is self.scripts
            "-P", f"out_dir:{main_out}",
            "-P", f"comparison:{comparison}",
            "-P", f"treatment:{treatment}",
            "-P", f"control:{control}",
            "-P", f"padj_threshold:{thresholds.padj}",
            "-P", f"lfc_threshold:{thresholds.lfc}",
            "-P", f"coldata_path:{cfg.samples.coldata}",
            "--output-dir", str(out_dir),
        ], cwd=self.scripts, log_file=self._log("m14_report", comparison))
        return ret == 0



    # ── M08b: QC plots (PCA, heatmap, strandedness) ─────────────────────────

    def _run_m08b_qc(self, comparison=None) -> bool:
        cfg      = self.cfg
        comp     = self._parse_comparison(comparison)
        if not comp:
            return False
        treatment, control = comp
        comp_dir = self._comparison_dir(comparison)
        outdir   = comp_dir / "deseq2"

        return run_r_script(self.scripts / "m08b_qc_plots.R", [
            "--norm-counts", str(comp_dir / "deseq2" / "normalized_counts.tsv"),
            "--de-genes",    str(comp_dir / "deseq2" / f"genes_{comparison}_DE.tsv"),
            "--coldata",     str(self._write_coldata(comparison)),
            "--gene-map",    str(self.out / "integrated" / "gene_type_map.tsv"),
            "--outdir",      str(outdir),
            "--treatment",   treatment,
            "--control",     control,
        ], log_file=self._log("m08b_qc", comparison)) == 0

    # ── M09: Satellite quantification ────────────────────────────────────────

    def _run_m09_satellite(self, comparison=None) -> bool:
        cfg     = self.cfg
        trim_dir = self.out / "trimmed"
        out_dir  = self.out / "satellite_counts"
        out_dir.mkdir(parents=True, exist_ok=True)
        sat_fa   = Path(cfg.databases.te_satellite)
        layout   = cfg.library.layout

        if not sat_fa.exists():
            warn(f"Satellite consensus not found: {sat_fa}")
            return True  # non-blocking

        # Build bowtie2 index if missing
        bt2_idx = sat_fa.parent / sat_fa.stem
        if not (sat_fa.parent / f"{sat_fa.stem}.1.bt2").exists() and            not (sat_fa.parent / f"{sat_fa.stem}.1.bt2l").exists():
            info(f"Building bowtie2 index for satellite consensus...")
            ret = run_cmd(["bowtie2-build", str(sat_fa), str(bt2_idx)],
                          log_file=self._log("m09_bt2_build"))
            if ret != 0:
                warn("bowtie2-build failed for satellite — skipping")
                return True  # non-blocking

        all_ok = True
        processed = 0
        fastq_dir = Path(cfg.samples.fastq_dir)
        for sid in cfg.samples.sample_ids:
            r1 = trim_dir / sid / f"{sid}_R1_trimmed.fastq.gz"
            # Fallback to original FASTQ if trimmed not available
            if not r1.exists():
                r1 = fastq_dir / sid / f"{sid}_R1.fastq.gz"
            if not r1.exists():
                r1 = fastq_dir / sid / f"{sid}_1.fastq.gz"
            if not r1.exists():
                warn(f"No reads found for {sid} — skipping satellite count")
                continue
            processed += 1
            sid_out = out_dir / sid
            sid_out.mkdir(exist_ok=True)
            reads = [str(r1)]
            if layout == "PE":
                r2 = trim_dir / sid / f"{sid}_R2_trimmed.fastq.gz"
                if r2.exists():
                    reads.append(str(r2))

            ret = run_cmd([
                "bowtie2",
                "--very-sensitive-local",
                "-x", str(sat_fa.parent / sat_fa.stem),
                "-1" if layout == "PE" else "-U", reads[0],
            ] + (["-2", reads[1]] if len(reads) > 1 else []) + [
                "--no-unal", "-p", str(cfg.threads),
                "-S", str(sid_out / f"{sid}_satellite.sam"),
            ], log_file=self._log("m09_satellite", sid))
            if ret != 0:
                all_ok = False
        if processed == 0:
            warn("No reads found for satellite counting — skipping")
            return True  # non-blocking
        return all_ok

    # ── M11: TE host genes ────────────────────────────────────────────────────

    def _run_m11_te_host(self, comparison=None) -> bool:
        cfg     = self.cfg
        out_dir = self.out / "m11"
        out_dir.mkdir(parents=True, exist_ok=True)
        te_counts = self.out / "integrated" / "te_counts_clean.tsv"

        return run_cmd([
            sys.executable, str(self.scripts / "m11_te_host_genes.py"),
            "--rmsk",      str(Path(cfg.databases.rmsk_bed)),
            "--gtf",       str(Path(cfg.databases.gtf)),
            "--outdir",    str(out_dir),
        ] + (["--te-counts", str(te_counts)] if te_counts.exists() else []),
            log_file=self._log("m11_te_host")) == 0

    # ── M11b: Promoter TE ────────────────────────────────────────────────────

    def _run_m11b_promoter_te(self, comparison=None) -> bool:
        cfg     = self.cfg
        out_dir = self.out / "m11b"
        out_dir.mkdir(parents=True, exist_ok=True)

        rmsk_hg38_bed = self.out / "m11" / "rmsk_hg38.bed"
        if not rmsk_hg38_bed.exists():
            warn(f"rmsk_hg38.bed not found at {rmsk_hg38_bed} — run M11 first")
            return False

        genome_fai = Path(cfg.databases.genome_fa).with_suffix(
            Path(cfg.databases.genome_fa).suffix + ".fai"
        )
        cmd = [
            sys.executable, "-m", "mosaic.m11b_promoter_te",
            "--rmsk-bed",  str(rmsk_hg38_bed),
            "--gtf",       str(Path(cfg.databases.gtf)),
            "--outdir",    str(out_dir),
        ]
        if genome_fai.exists():
            cmd += ["--genome-fai", str(genome_fai)]
        return run_cmd(cmd, cwd=self.scripts.parent,
                       log_file=self._log("m11b_promoter_te")) == 0

    # ── M11c: TE age ─────────────────────────────────────────────────────────

    def _run_m11c_te_age(self, comparison=None) -> bool:
        cfg = self.cfg
        comp = self._parse_comparison(comparison)
        if not comp:
            return False
        treatment, control = comp

        # Resolve comparison name for file paths (None → first non-interaction)
        if comparison is None:
            comp_name = None
            for c in (cfg.comparisons or []):
                if c.get("type") != "interaction":
                    comp_name = f"{c['treatment']}_vs_{c['control']}"
                    break
            if comp_name is None:
                return False
        else:
            comp_name = comparison

        comp_dir = self._comparison_dir(comp_name)
        te_de    = comp_dir / "deseq2" / f"TE_{comp_name}_DE.tsv"
        out_dir  = self.out / "m11c"
        out_dir.mkdir(parents=True, exist_ok=True)

        return run_r_script(self.scripts / "m11c_te_age.R", [
            "--te-de",     str(te_de),
            "--rmsk-full", str(Path(cfg.databases.rmsk_bed)),  # UCSC gz for milliDiv
            "--outdir",    str(out_dir),
            "--treatment", treatment,
            "--control",   control,
        ], log_file=self._log("m11c_te_age", comp_name)) == 0

    # ── M11d: TE DE integration ───────────────────────────────────────────────

    def _run_m11d_te_de(self, comparison=None) -> bool:
        cfg      = self.cfg
        comp     = self._parse_comparison(comparison)
        # Interaction comparisons: use fraction_of_interest as treatment, its base as control
        if not comp:
            if not self._is_interaction_comparison(comparison):
                return False
            comp_cfg = self._get_comparison_config(comparison)
            frac = comp_cfg.get("fraction_of_interest", "ChaRNA") if comp_cfg else "ChaRNA"
            treatment, control = f"mito-{frac}", frac
        else:
            treatment, control = comp
        comp_dir = self._comparison_dir(comparison)
        out_dir  = comp_dir / "m11d"
        out_dir.mkdir(parents=True, exist_ok=True)

        rmsk_hg38_bed = self.out / "m11" / "rmsk_hg38.bed"
        if not rmsk_hg38_bed.exists():
            warn(f"rmsk_hg38.bed not found at {rmsk_hg38_bed} — run M11 first")
            return False

        return run_r_script(self.scripts / "m11d_te_de_integration.R", [
            "--de-genes",      str(comp_dir / "deseq2" / f"genes_{comparison}_DE.tsv"),
            "--norm-counts",   str(comp_dir / "deseq2" / "normalized_counts.tsv"),
            "--te-prom-genes", str(self.out / "m11b" / "te_promoter_genes.tsv"),
            "--rmsk-bed",      str(rmsk_hg38_bed),
            "--gtf",           str(Path(cfg.databases.gtf)),
            "--gene-map",      str(self.out / "integrated" / "gene_type_map.tsv"),
            "--outdir",        str(out_dir),
            "--treatment",     treatment,
            "--control",       control,
            "--padj",          str(cfg.thresholds.padj),
            "--lfc",           str(cfg.thresholds.lfc),
        ], log_file=self._log("m11d_te_de", comparison)) == 0

    # ── M12: TEtranscripts locus-level ────────────────────────────────────────

    def _run_m12_tetranscripts(self, comparison=None) -> bool:
        cfg     = self.cfg
        bam_dir = self.out / "mapping" / "bam"
        out_dir = self.out / "m12"
        out_dir.mkdir(parents=True, exist_ok=True)
        strand  = self._get_consensus_strandedness()
        stranded = "no" if strand == "unstranded" else (
            "reverse" if strand == "RF" else "yes")

        # Use the processed 6-col BED (m11/rmsk_hg38.bed), NOT the raw UCSC gz
        rmsk_hg38_bed = self.out / "m11" / "rmsk_hg38.bed"
        if not rmsk_hg38_bed.exists():
            warn(f"rmsk_hg38.bed not found at {rmsk_hg38_bed} — run M11 first")
            return False
        rc = run_cmd([
            sys.executable, str(self.scripts / "m12_te_locus.py"),
            "locus-count",
            "--bam-dir",  str(bam_dir),
            "--rmsk-bed", str(rmsk_hg38_bed),
            "--outdir",   str(out_dir),
        ], log_file=self._log("m12_tetranscripts"))
        if rc != 0:
            return False

        locus_counts = out_dir / "locus_counts.tsv"
        if not locus_counts.exists():
            warn("m12: locus_counts.tsv not produced — skipping locus DE")
            return False

        # Per-comparison locus-level DESeq2 (m12_deseq2.R). Produces
        # comparison_<comp>/m12/<comp>_te_de_sig.tsv, using the canonical
        # stable-gene size factors from m08. Failures are per-comparison and
        # non-fatal; the count matrix (module core product) is already written.
        for c in (cfg.comparisons or []):
            if c.get("type") == "interaction":
                continue  # locus DE needs a 2-group contrast
            treat, ctrl = c["treatment"], c["control"]
            comp = f"{treat}_vs_{ctrl}"
            comp_m12 = self._comparison_dir(comp) / "m12"
            comp_m12.mkdir(parents=True, exist_ok=True)
            deseq2_dir = self._comparison_dir(comp) / "deseq2"
            coldata = self.out / f"coldata_{comp}.tsv"
            if not coldata.exists():
                coldata = Path(cfg.samples.coldata)
            args = [
                "--locus-counts", str(locus_counts),
                "--coldata",      str(coldata),
                "--outdir",       str(comp_m12),
                "--treatment",    treat,
                "--control",      ctrl,
                "--comparison",   comp,
                "--padj",         str(cfg.thresholds.padj),
                "--lfc",          str(cfg.thresholds.lfc),
                "--threads",      str(min(cfg.threads, 4)),
            ]
            gc = deseq2_dir / "all_genes_counts.tsv"
            sg = deseq2_dir / "stable_genes.tsv"
            if gc.exists():
                args += ["--gene-counts", str(gc)]
            if sg.exists():
                args += ["--stable-genes", str(sg)]
            if run_r_script(self.scripts / "m12_deseq2.R", args,
                            log_file=self._log("m12_deseq2", comp)) != 0:
                warn(f"m12 locus-level DE failed for {comp} (non-fatal)")
        return True


    # ── M10: GO/KEGG enrichment ──────────────────────────────────────────────

    def _run_m10_go(self, comparison=None) -> bool:
        cfg  = self.cfg
        comp = self._parse_comparison(comparison)
        # Interaction comparisons don't have treatment/control; still run GO
        if not comp and not self._is_interaction_comparison(comparison):
            return False
        comp_dir = self._comparison_dir(comparison)
        de_file  = comp_dir / "deseq2" / f"genes_{comparison}_DE_significant.tsv"
        if not de_file.exists():
            warn(f"DE significant file not found: {de_file} — run m08_deseq2 first")
            return False
        out_dir = comp_dir / "go"
        out_dir.mkdir(parents=True, exist_ok=True)
        comp_dir = self._comparison_dir(comparison)
        all_genes = comp_dir / "deseq2" / "all_genes_counts.tsv"
        te_de     = comp_dir / "deseq2" / f"TE_{comparison}_DE.tsv"
        te_host       = self.out / "m11" / "te_host_universe.txt"
        te_host_g     = self.out / "m11" / "te_host_genes.tsv"
        te_prom_int   = self.out / "m11b" / "te_promoter_intersect.bed"
        gene_map      = self.out / "integrated" / "gene_type_map.tsv"

        args = [
            "--de-genes",      str(de_file),
            "--all-genes",     str(all_genes),
            "--gene-type-map", str(gene_map),
            "--outdir",        str(out_dir),
            "--organism",      cfg.project.organism,
            "--padj",          str(cfg.thresholds.padj),
        ]
        if te_de.exists():
            args += ["--te-de", str(te_de)]
        if te_host.exists():
            args += ["--te-host-universe", str(te_host)]
        if te_host_g.exists():
            args += ["--te-host-genes", str(te_host_g)]
        if te_prom_int.exists():
            args += ["--te-promoter-intersect", str(te_prom_int)]

        return run_r_script(self.scripts / "m10_go.R", args,
                            log_file=self._log("m10_go", comparison)) == 0

    # ── M16: StringTie de novo ───────────────────────────────────────────────

    def _run_m16_stringtie(self, comparison=None) -> bool:
        cfg     = self.cfg
        bam_dir = self.out / "mapping" / "bam"
        out_dir = self.out / "m16_stringtie"        # global: assembly + merge
        out_dir.mkdir(parents=True, exist_ok=True)
        comp    = self._parse_comparison(comparison)
        if not comp:
            return False
        treatment, control = comp
        de_dir  = self._comparison_dir(comparison) / "m16_stringtie"  # per-comparison DE
        de_dir.mkdir(parents=True, exist_ok=True)

        # StringTie per-sample assembly
        gtf = Path(cfg.databases.gtf)
        for sid in cfg.samples.sample_ids:
            bam = bam_dir / sid / f"{sid}_Aligned.sortedByCoord.out.bam"
            if not bam.exists():
                continue
            sid_out = out_dir / f"{sid}_quant"
            sid_out.mkdir(exist_ok=True)
            ret = run_cmd([
                "stringtie", str(bam),
                "-G", str(gtf),
                "-o", str(sid_out / f"{sid}.gtf"),
                "-e", "-b", str(sid_out),
                "-p", str(cfg.threads),
            ], log_file=self._log("m16_stringtie_asm", sid))
            if ret != 0:
                warn(f"StringTie assembly failed for {sid}")

        # StringTie merge with gffcmp
        merged_gtf = out_dir / "merged.gtf"
        gtf_list   = out_dir / "gtf_list.txt"
        sample_gtfs = list(out_dir.glob("*/*.gtf"))

        if sample_gtfs:
            gtf_list.write_text("\n".join(str(g) for g in sample_gtfs))
            run_cmd(["stringtie", "--merge",
                     "-G", str(gtf),
                     "-o", str(merged_gtf),
                     str(gtf_list)],
                    log_file=self._log("m16_stringtie_merge"))
            run_cmd(["gffcompare",
                     "-r", str(gtf),
                     "-o", str(out_dir / "gffcmp"),
                     str(merged_gtf)],
                    log_file=self._log("m16_gffcmp"))

        # ── novel_transcripts.tmap: filter tmap to novel class codes ──────────
        tmap       = out_dir / "gffcmp.merged.gtf.tmap"
        novel_tmap = out_dir / "novel_transcripts.tmap"
        novel_codes = {"u", "i", "j", "x", "s", "e", "p"}
        if tmap.exists() and not novel_tmap.exists():
            novel_lines = []
            with open(tmap) as fh:
                fh.readline()  # skip header — R reads without header
                for line in fh:
                    parts = line.split("\t")
                    if len(parts) >= 3 and parts[2] in novel_codes:
                        novel_lines.append(line)
            novel_tmap.write_text("".join(novel_lines))
            info(f"M16: {len(novel_lines)} novel transcripts → novel_transcripts.tmap")

        # ── featureCounts: raw integer counts from all BAMs vs merged GTF ────
        # Runs once (idempotent); output is reused across comparisons.
        fc_output = out_dir / "featurecounts_transcripts.txt"
        merged_gtf = out_dir / "merged.gtf"
        if merged_gtf.exists() and not fc_output.exists():
            all_bams = [
                bam_dir / sid / f"{sid}_Aligned.sortedByCoord.out.bam"
                for sid in cfg.samples.sample_ids
            ]
            all_bams = [b for b in all_bams if b.exists()]
            if all_bams:
                fc_cmd = [
                    "featureCounts",
                    "-T", str(min(cfg.threads, 4)),
                    "-g", "transcript_id",
                    "-f",          # exon-level features, grouped by transcript_id
                    "-a", str(merged_gtf),
                    "-o", str(fc_output),
                ]
                if cfg.library.layout == "PE":
                    fc_cmd.append("-p")
                fc_cmd += [str(b) for b in all_bams]
                run_cmd(fc_cmd, log_file=self._log("m16_featurecounts"))
                if fc_output.exists():
                    info(f"M16: featureCounts complete → {fc_output.name}")
                else:
                    warn("M16: featureCounts did not produce output — novel DE will be skipped")

        # DE analysis — results go to per-comparison dir; assembly files from global dir
        args = [
            "--st-dir",    str(out_dir),
            "--coldata",   str(self._write_coldata(comparison)),
            "--outdir",    str(de_dir),
            "--treatment", treatment,
            "--control",   control,
            "--padj",      str(cfg.thresholds.padj),
            "--lfc",       str(cfg.thresholds.lfc),
        ]
        if tmap.exists():
            args += ["--tmap", str(tmap)]
        if novel_tmap.exists():
            args += ["--novel-tmap", str(novel_tmap)]
        if fc_output.exists():
            args += ["--feature-counts", str(fc_output)]
        return run_r_script(self.scripts / "m16_stringtie_de.R", args,
                            log_file=self._log("m16_stringtie_r", comparison)) == 0

    # ── M17: rMATS (already defined above, keep stub) ───────────────────────

    # ── M18: Network (STRINGdb + WGCNA) ────────────────────────────────────

    def _run_m18_network(self, comparison=None) -> bool:
        cfg      = self.cfg
        comp     = self._parse_comparison(comparison)
        if not comp:
            return False
        treatment, control = comp
        comp_dir = self._comparison_dir(comparison)
        de_file  = comp_dir / "deseq2" / f"genes_{comparison}_DE.tsv"
        out_dir  = comp_dir / "m18"
        out_dir.mkdir(parents=True, exist_ok=True)
        return run_r_script(self.scripts / "m18_network.R", [
            "--de-genes",    str(de_file),
            "--norm-counts", str(comp_dir / "deseq2" / "normalized_counts.tsv"),
            "--gene-map",    str(self.out / "integrated" / "gene_type_map.tsv"),
            "--coldata",     str(self._write_coldata(comparison)),
            "--outdir",      str(out_dir),
            "--padj",        str(cfg.thresholds.padj),
            "--lfc",         str(cfg.thresholds.lfc),
        ], log_file=self._log("m18_network", comparison)) == 0

    # ── M19: RNA editing ────────────────────────────────────────────────────

    def _run_m19_editing(self, comparison=None) -> bool:
        cfg     = self.cfg
        bam_dir = self.out / "mapping" / "bam"
        comp    = self._parse_comparison(comparison)
        if not comp:
            return False
        treatment, control = comp
        out_dir = self._comparison_dir(comparison) / "m19_editing"
        out_dir.mkdir(parents=True, exist_ok=True)

        normal_bams = self._get_condition_bams(control)
        tumor_bams  = self._get_condition_bams(treatment)

        # m19 takes comma-separated BAM lists
        ret = run_cmd([
            sys.executable, "-m", "mosaic.m19_rna_editing",
            "--bam-normal",  ",".join(normal_bams),
            "--bam-tumor",  ",".join(tumor_bams),
            "--ref",        str(Path(cfg.databases.genome_fa)),
            "--gtf",        str(Path(cfg.databases.gtf)),
            "--outdir",     str(out_dir),
            "--threads",    str(cfg.threads),
            "--min-cov",    str(cfg.thresholds.min_cov_editing),
            "--min-score",  str(cfg.thresholds.min_score_editing),
            "--min-ef",     str(cfg.thresholds.min_ef_editing),
        ], cwd=self.scripts.parent, log_file=self._log("m19_editing", comparison))
        if ret != 0:
            return False

        # R plots
        return run_r_script(self.scripts / "m19_editing_plots.R", [
            "--editing-dir", str(out_dir),
            "--outdir",      str(out_dir),
            "--treatment",   treatment,
            "--control",     control,
        ], log_file=self._log("m19_editing_r", comparison)) == 0

    # ── M20: Metatranscriptomics (Kraken2/Bracken) ──────────────────────────

    def _run_m20_meta(self, comparison=None) -> bool:
        cfg     = self.cfg
        out_dir = self.out / "m20_metatranscriptomics"
        out_dir.mkdir(parents=True, exist_ok=True)
        bam_dir = self.out / "mapping" / "bam"

        for sid in cfg.samples.sample_ids:
            unmapped_r1 = bam_dir / sid / f"{sid}_Unmapped.out.mate1"
            if not unmapped_r1.exists():
                warn(f"Unmapped reads not found for {sid} — was STAR run with --outReadsUnmapped Fastx?")
                continue
            sid_out = out_dir / sid
            sid_out.mkdir(exist_ok=True)
            # Kraken2 — STAR unmapped files are plain FASTQ (not gzip)
            k2_report = sid_out / f"{sid}_kraken2_report.txt"
            unmapped_r2 = bam_dir / sid / f"{sid}_Unmapped.out.mate2"
            k2_cmd = [
                "kraken2",
                "--db",      str(Path(cfg.databases.kraken2_db)),
                "--threads", str(cfg.threads),
                "--report",  str(k2_report),
                "--output",  str(sid_out / f"{sid}_kraken2_output.txt"),
            ]
            if unmapped_r2.exists():
                k2_cmd += ["--paired", str(unmapped_r1), str(unmapped_r2)]
            else:
                k2_cmd.append(str(unmapped_r1))
            run_cmd(k2_cmd, log_file=self._log("m20_kraken2", sid))
            # Bracken abundance re-estimation at species level
            read_len = getattr(cfg.library, "read_length", 150)
            bracken_out = sid_out / f"{sid}_bracken_species.txt"
            run_cmd([
                "bracken",
                "-d",  str(Path(cfg.databases.kraken2_db)),
                "-i",  str(k2_report),
                "-o",  str(bracken_out),
                "-r",  str(read_len),
                "-l",  "S",
                "-t",  str(cfg.thresholds.min_reads_meta),
            ], log_file=self._log("m20_bracken", sid))

        # R summary + diversity — use full coldata (all samples, not comparison-specific)
        return run_r_script(self.scripts / "m20_metatranscriptomics.R", [
            "--meta-dir",  str(out_dir),
            "--outdir",    str(out_dir),
            "--coldata",   str(Path(cfg.samples.coldata)),
            "--min-reads", str(cfg.thresholds.min_reads_meta),
        ], log_file=self._log("m20_meta")) == 0

    # ── M20b: HUMAnN3 functional profiling ──────────────────────────────────

    def _run_m20b_humann(self, comparison=None) -> bool:
        cfg     = self.cfg
        bam_dir = self.out / "mapping" / "bam"
        out_dir = self.out / "m20_metatranscriptomics" / "humann3"
        out_dir.mkdir(parents=True, exist_ok=True)
        env     = cfg.environments.humann or "humann3"

        # Check humann env
        import subprocess as _sp
        env_check = _sp.run(["conda", "run", "-n", env,
                              "humann", "--version"],
                             capture_output=True, text=True)
        if env_check.returncode != 0:
            warn("HUMAnN3 functional profiling (M20b) is not available in the "
                 f"current environment ('{env}'). M20b is experimental and "
                 "requires a separate conda environment; see README §Optional "
                 "environments. Skipping M20b — the rest of the pipeline continues.")
            return True  # non-blocking

        all_ok = True
        for sid in cfg.samples.sample_ids:
            unmapped = bam_dir / sid / f"{sid}_Unmapped.out.mate1"
            if not unmapped.exists():
                warn(f"Unmapped reads not found for {sid}")
                continue
            sid_out = out_dir / sid
            sid_out.mkdir(exist_ok=True)
            ret = run_cmd([
                "humann",
                "--input",               str(unmapped),
                "--output",              str(sid_out),
                "--threads",             str(cfg.threads),
                "--nucleotide-database", str(Path(cfg.databases.chocophlan_db)),
                "--protein-database",    str(Path(cfg.databases.uniref_db)),
            ], env_name=env, log_file=self._log("m20b_humann", sid))
            if ret != 0:
                all_ok = False

        if not all_ok:
            return False

        # Merge tables
        for table_type in ["pathabundance", "genefamilies", "pathcoverage"]:
            for sid in cfg.samples.sample_ids:
                src = out_dir / sid
                for f in src.glob(f"*_{table_type}.tsv"):
                    import shutil
                    shutil.copy(f, out_dir / f"{sid}_{table_type}.tsv")
            run_cmd([
                "humann_join_tables",
                "--input",   str(out_dir),
                "--output",  str(out_dir / f"all_{table_type}.tsv"),
                "--file_name", table_type,
            ], env_name=env)

        # Normalize
        run_cmd([
            "humann_renorm_table",
            "--input",  str(out_dir / "all_pathabundance.tsv"),
            "--output", str(out_dir / "all_pathabundance_cpm.tsv"),
            "--units",  "cpm",
        ], env_name=env)
        return True

    # ── M23: TE somatic insertions ───────────────────────────────────────────

    def _run_m23_te_insertions(self, comparison=None) -> bool:
        cfg     = self.cfg
        bam_dir = self.out / "mapping" / "bam"
        comp    = self._parse_comparison(comparison)
        if not comp:
            return False
        treatment, control = comp

        # Per-comparison output directory to avoid overwrites across comparisons
        comp_name = comparison or f"{treatment}_vs_{control}"
        out_dir   = self.out / "m23_te_insertions" / comp_name
        out_dir.mkdir(parents=True, exist_ok=True)

        normal_chim = ",".join(
            str(bam_dir / sid / f"{sid}_Chimeric.out.sam")
            for sid in cfg.samples.get_condition_samples(control))
        tumor_chim  = ",".join(
            str(bam_dir / sid / f"{sid}_Chimeric.out.sam")
            for sid in cfg.samples.get_condition_samples(treatment))

        return run_cmd([
            sys.executable, "-m", "mosaic.m23_te_insertions",
            "--chimeric-normal",       normal_chim,
            "--chimeric-tumor",        tumor_chim,
            "--te-consensus",          str(Path(cfg.databases.te_consensus)),
            "--genome-fa",             str(Path(cfg.databases.genome_fa)),
            "--gtf",                   str(Path(cfg.databases.gtf)),
            "--outdir",                str(out_dir),
            "--threads",               str(cfg.threads),
            "--min-clip",              str(cfg.thresholds.min_clip_te),
            "--min-support",           str(cfg.thresholds.min_support_te),
            "--min-mapq",              "0",
            "--condition-control",     control,
            "--condition-treatment",   treatment,
        ], cwd=self.scripts.parent, log_file=self._log("m23_te_insertions", comp_name)) == 0

    # ── M15: Fusion genes (Arriba) ───────────────────────────────────────────

    def _run_m15_fusion(self, comparison=None) -> bool:
        cfg     = self.cfg
        bam_dir = self.out / "mapping" / "bam"
        out_dir = self.out / "m15_fusions"
        out_dir.mkdir(parents=True, exist_ok=True)
        arriba  = Path(cfg.databases.arriba_bin)
        db      = Path(cfg.databases.arriba_db)

        all_ok = True
        for sid in cfg.samples.sample_ids:
            bam  = bam_dir / sid / f"{sid}_Aligned.sortedByCoord.out.bam"
            chim = bam_dir / sid / f"{sid}_Chimeric.out.sam"
            if not bam.exists():
                continue
            cmd = [
                str(arriba),
                "-x", str(bam),
                "-o", str(out_dir / f"{sid}_fusions.tsv"),
                "-O", str(out_dir / f"{sid}_fusions_discarded.tsv"),
                "-a", str(Path(cfg.databases.genome_fa)),
                "-g", str(Path(cfg.databases.gtf)),
                "-b", str(db / "blacklist_hg38_GRCh38_v2.4.0.tsv.gz"),
                "-k", str(db / "known_fusions_hg38_GRCh38_v2.4.0.tsv.gz"),
                "-t", str(db / "known_fusions_hg38_GRCh38_v2.4.0.tsv.gz"),
                "-p", str(db / "protein_domains_hg38_GRCh38_v2.4.0.gff3"),
            ]
            if chim.exists():
                cmd += ["-c", str(chim)]
            ret = run_cmd(cmd, log_file=self._log("m15_fusion", sid))
            if ret != 0:
                all_ok = False

        if not all_ok:
            return False

        # Merge all fusions — prepend sample_id column to header and data
        merged = out_dir / "fusions_merged.tsv"
        header_written = False
        with open(merged, "w") as fout:
            for sid in cfg.samples.sample_ids:
                f = out_dir / f"{sid}_fusions.tsv"
                if not f.exists():
                    continue
                with open(f) as fin:
                    for i, line in enumerate(fin):
                        if i == 0:
                            if not header_written:
                                # Remove leading # and prepend sample_id column
                                fout.write("sample_id\t" + line.lstrip("#"))
                                header_written = True
                        else:
                            if line.strip():
                                fout.write(f"{sid}\t{line}")
        return True


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
            # Use first non-interaction comparison from config
            for c in (self.cfg.comparisons or []):
                if c.get("type") != "interaction":
                    return c["treatment"], c["control"]
            return None
        # Check if this is a named interaction comparison
        for c in (self.cfg.comparisons or []):
            if c.get("type") == "interaction" and c.get("name") == comparison:
                return None  # Not a simple 2-group comparison
        parts = comparison.split("_vs_")
        if len(parts) == 2:
            return parts[0], parts[1]
        return None

    def _is_interaction_comparison(self, comparison: Optional[str]) -> bool:
        """Return True if the given comparison is an interaction-type comparison."""
        c = self._get_comparison_config(comparison)
        return bool(c and c.get("type") == "interaction")

    def _get_comparison_config(self, comparison: Optional[str]) -> Optional[dict]:
        """Return the full comparison config dict for a given comparison name."""
        if not comparison:
            return (self.cfg.comparisons or [None])[0]
        for c in (self.cfg.comparisons or []):
            name = c.get("name", f"{c.get('treatment','')}_vs_{c.get('control','')}")
            if name == comparison:
                return c
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
        # If config has an explicit value (not "auto"), it always wins
        cfg_strand = (getattr(self.cfg.library, "_strandedness_config", None) or
                      getattr(self.cfg.library, "strandedness", None) or "")
        if cfg_strand and cfg_strand != "auto":
            return cfg_strand

        # Otherwise use auto-detection result from m07
        strand_map = getattr(self.cfg, "_strandedness", {})
        if strand_map:
            from collections import Counter
            return Counter(strand_map.values()).most_common(1)[0][0]

        # Final fallback: unstranded (safe default)
        return "unstranded"

    def _write_coldata(self, comparison: str) -> Path:
        """Write coldata TSV for a specific comparison."""
        comp_cfg = self._get_comparison_config(comparison)
        coldata_rows = self.cfg.samples.coldata_df  # list of dicts

        # Interaction comparison: include all samples (no subsetting)
        if comp_cfg and comp_cfg.get("type") == "interaction":
            if not coldata_rows:
                return Path(self.cfg.samples.coldata)
            out_path = self.out / f"coldata_{comparison}.tsv"
            out_path.parent.mkdir(parents=True, exist_ok=True)
            keys = list(coldata_rows[0].keys())
            with open(out_path, "w") as fh:
                fh.write("\t".join(keys) + "\n")
                for row in coldata_rows:
                    fh.write("\t".join(str(row.get(k, "")) for k in keys) + "\n")
            return out_path

        comp = self._parse_comparison(comparison)
        if not comp:
            return Path(self.cfg.samples.coldata)

        treatment, control = comp
        if not coldata_rows:
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
            # Use v2 modules config if available
            mods_v2 = getattr(self.cfg, "modules_v2", None)
            if mods_v2 and hasattr(mods_v2, "enabled"):
                enabled = mods_v2.enabled()
            else:
                # Fallback to old modules config
                mods_cfg = self.cfg.modules
                enabled = [m for m, v in vars(mods_cfg).items()
                           if isinstance(v, bool) and v]

        # If --from-module specified, mark earlier modules as skipped
        # (regardless of --force — we never re-run modules before from_module)
        if self.from_mod:
            all_mods = self._resolve_order(list(MODULE_DEPS.keys()))
            idx = all_mods.index(self.from_mod) if self.from_mod in all_mods else 0
            for mod_name in all_mods[:idx]:
                self.state.mark(mod_name, "skipped")
                m = self.state.modules.get(mod_name)
                if m:
                    m._from_module_skip = True
            self.state.save()

        # Resolve execution order
        order = self._resolve_order(enabled)

        # Get comparisons
        comparisons = []
        for c in getattr(self.cfg, "comparisons", []):
            if c.get("type") == "interaction":
                comparisons.append(c.get("name", "interaction"))
            else:
                comparisons.append(f"{c['treatment']}_vs_{c['control']}")
        if not comparisons:
            comparisons = [None]  # single comparison from config

        # Comparison-independent modules
        COMP_INDEPENDENT = {
            "m00_download", "m02_trim", "m04_star",
            "m05_te_count", "m07_strandedness",
            "m09_satellite", "m11_te_host",
            "m11b_promoter_te",
            "m20_meta", "m20b_humann",
            "m12_tetranscripts",
            "m15_fusion", "m16_stringtie",
        }
        # These run per-comparison (need comparison-specific DE data)
        # m11c_te_age, m11d_te_de, m17_rmats, m18_network,
        # m19_editing, m23_te_insertions → handled in else branch

        section("MOSAIC Pipeline")
        info(f"Modules: {', '.join(order)}")
        info(f"Comparisons: {', '.join(str(c) for c in comparisons)}")
        info(f"Force re-run: {self.force}")

        all_ok = True
        start  = datetime.now()

        # gene_type_map.tsv is normally generated inside _run_m04_star, but any
        # run that skips m04 via --from-module or explicit --modules must also have
        # it (needed by m10_go, m11_te_host, m18_network). Generate it here if absent.
        gene_map = self.runner.out / "integrated" / "gene_type_map.tsv"
        if not gene_map.exists():
            try:
                from .refs import extract_gene_type_map
                (self.runner.out / "integrated").mkdir(parents=True, exist_ok=True)
                extract_gene_type_map(Path(self.cfg.databases.gtf), gene_map)
                info(f"gene_type_map.tsv written: {gene_map}")
            except Exception as e:
                warn(f"Could not build gene_type_map.tsv: {e}")

        # Robustness policy: only the backbone modules are critical — without them
        # nothing downstream can run, so their failure aborts the pipeline. Every
        # other module (heavy/optional/screening: m18 WGCNA, m20 metatranscriptomics,
        # m23 insertions, m16 StringTie, m06 circRNA, m15 fusion, …) is non-critical:
        # if it fails (e.g. out-of-memory) we log it, mark it failed in the checkpoint
        # state, and continue so the remaining core results are still produced. A
        # rerun resumes from the checkpoint and re-attempts only the failed modules.
        # Set continue_on_error: false in config to make ALL failures fatal.
        CRITICAL_MODULES = {"m02_trim", "m04_star", "m05_te_count", "m08_deseq2"}
        continue_on_error = bool(getattr(self.cfg, "continue_on_error", True))
        failed_modules: List[str] = []

        for module in order:
            if module in COMP_INDEPENDENT:
                ok = self.runner.run_module(module)
            else:
                # Run once per comparison
                ok = True
                for comp in comparisons:
                    ok = ok and self.runner.run_module(module, comp)

            if not ok:
                all_ok = False
                if module in CRITICAL_MODULES or not continue_on_error:
                    error(f"Module {module} failed — stopping pipeline"
                          + (" (critical backbone module)"
                             if module in CRITICAL_MODULES else ""))
                    break
                # Non-critical failure: keep the checkpoint, warn, carry on.
                failed_modules.append(module)
                warn(f"Module {module} failed (non-critical) — continuing; "
                     f"downstream modules needing its output will skip or fail "
                     f"gracefully. Rerun to retry from checkpoint.")

        if failed_modules:
            warn("Completed with non-fatal module failures: "
                 + ", ".join(failed_modules))

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
