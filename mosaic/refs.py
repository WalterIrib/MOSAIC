"""
mosaic/refs.py
Reference DB management system.
Implements the full mosaic-refs logic from Design Document v0.7 §2.
"""
from __future__ import annotations

import csv
import hashlib
import json
import os
import shutil
import subprocess
import sys
import time
import urllib.request
import urllib.error
from dataclasses import dataclass, field, asdict
from datetime import datetime
from pathlib import Path
from typing import Optional

import yaml

from .log import get_logger, info, warn, error, success, section

# ─────────────────────────────────────────────────────────────────────────────
# Known organisms catalogue
# ─────────────────────────────────────────────────────────────────────────────

@dataclass
class OrganismSpec:
    organism:     str
    genome_build: str
    ensembl_release: str        # e.g. "111"
    ensembl_fasta_url: str
    ensembl_gtf_url:  str
    annotation_source: str      # GENCODE | Ensembl
    gencode_release: Optional[str] = None
    gencode_fasta_url: Optional[str] = None
    gencode_gtf_url:  Optional[str] = None
    dfam_assembly: str = ""     # DFAM assembly name for this build


# URLs are templates — filled at runtime with release numbers
_ORGANISM_CATALOGUE: dict[str, OrganismSpec] = {
    "human": OrganismSpec(
        organism="human",
        genome_build="hg38",
        ensembl_release="111",
        ensembl_fasta_url=(
            "https://ftp.ensembl.org/pub/release-111/fasta/"
            "homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz"
        ),
        ensembl_gtf_url=(
            "https://ftp.ensembl.org/pub/release-111/gtf/"
            "homo_sapiens/Homo_sapiens.GRCh38.111.gtf.gz"
        ),
        annotation_source="GENCODE",
        gencode_release="v44",
        gencode_fasta_url=(
            "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/"
            "release_44/GRCh38.primary_assembly.genome.fa.gz"
        ),
        gencode_gtf_url=(
            "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/"
            "release_44/gencode.v44.primary_assembly.annotation.gtf.gz"
        ),
        dfam_assembly="hg38",
    ),
    "mouse": OrganismSpec(
        organism="mouse",
        genome_build="mm39",
        ensembl_release="111",
        ensembl_fasta_url=(
            "https://ftp.ensembl.org/pub/release-111/fasta/"
            "mus_musculus/dna/Mus_musculus.GRCm39.dna_sm.primary_assembly.fa.gz"
        ),
        ensembl_gtf_url=(
            "https://ftp.ensembl.org/pub/release-111/gtf/"
            "mus_musculus/Mus_musculus.GRCm39.111.gtf.gz"
        ),
        annotation_source="GENCODE",
        gencode_release="vM34",
        gencode_fasta_url=(
            "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/"
            "release_M34/GRCm39.primary_assembly.genome.fa.gz"
        ),
        gencode_gtf_url=(
            "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/"
            "release_M34/gencode.vM34.primary_assembly.annotation.gtf.gz"
        ),
        dfam_assembly="mm39",
    ),
    "rat": OrganismSpec(
        organism="rat",
        genome_build="rn7",
        ensembl_release="111",
        ensembl_fasta_url=(
            "https://ftp.ensembl.org/pub/release-111/fasta/"
            "rattus_norvegicus/dna/Rattus_norvegicus.mRatBN7.2.dna_sm.primary_assembly.fa.gz"
        ),
        ensembl_gtf_url=(
            "https://ftp.ensembl.org/pub/release-111/gtf/"
            "rattus_norvegicus/Rattus_norvegicus.mRatBN7.2.111.gtf.gz"
        ),
        annotation_source="Ensembl",
        dfam_assembly="rn7",
    ),
    "zebrafish": OrganismSpec(
        organism="zebrafish",
        genome_build="GRCz11",
        ensembl_release="111",
        ensembl_fasta_url=(
            "https://ftp.ensembl.org/pub/release-111/fasta/"
            "danio_rerio/dna/Danio_rerio.GRCz11.dna_sm.primary_assembly.fa.gz"
        ),
        ensembl_gtf_url=(
            "https://ftp.ensembl.org/pub/release-111/gtf/"
            "danio_rerio/Danio_rerio.GRCz11.111.gtf.gz"
        ),
        annotation_source="Ensembl",
        dfam_assembly="GRCz11",
    ),
}

# ─────────────────────────────────────────────────────────────────────────────
# DB path helpers
# ─────────────────────────────────────────────────────────────────────────────

def org_path(db_root: Path, organism: str, build: str) -> Path:
    return db_root / organism / build


def manifest_path(db_root: Path, organism: str, build: str) -> Path:
    return org_path(db_root, organism, build) / "manifest.yaml"


def update_log_path(db_root: Path, organism: str, build: str) -> Path:
    return org_path(db_root, organism, build) / "update_log.tsv"


# ─────────────────────────────────────────────────────────────────────────────
# Manifest read / write
# ─────────────────────────────────────────────────────────────────────────────

def read_manifest(db_root: Path, organism: str, build: str) -> dict:
    p = manifest_path(db_root, organism, build)
    if not p.exists():
        return {}
    with open(p) as fh:
        return yaml.safe_load(fh) or {}


def write_manifest(db_root: Path, organism: str, build: str,
                   data: dict) -> None:
    p = manifest_path(db_root, organism, build)
    p.parent.mkdir(parents=True, exist_ok=True)
    with open(p, "w") as fh:
        yaml.dump(data, fh, default_flow_style=False, sort_keys=False)


# ─────────────────────────────────────────────────────────────────────────────
# Update log
# ─────────────────────────────────────────────────────────────────────────────

def append_update_log(db_root: Path, organism: str, build: str,
                      component: str, old_ver: str, new_ver: str,
                      action: str, note: str = "") -> None:
    p = update_log_path(db_root, organism, build)
    p.parent.mkdir(parents=True, exist_ok=True)
    header = "timestamp\toperator\tcomponent\told_version\tnew_version\taction\tnote\n"
    write_header = not p.exists()

    ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    operator = os.environ.get("USER", os.environ.get("USERNAME", "unknown"))
    row = f"{ts}\t{operator}\t{component}\t{old_ver}\t{new_ver}\t{action}\t{note}\n"

    with open(p, "a") as fh:
        if write_header:
            fh.write(header)
        fh.write(row)


# ─────────────────────────────────────────────────────────────────────────────
# NFS-safe lock
# ─────────────────────────────────────────────────────────────────────────────

class RefLock:
    """Atomic lock using mkdir (NFS-safe). Cleans stale locks (mtime > 2h)."""

    STALE_SECONDS = 7200  # 2 hours

    def __init__(self, db_root: Path, organism: str, build: str):
        self._lock_dir = org_path(db_root, organism, build) / ".lock"

    def __enter__(self):
        self._acquire()
        return self

    def __exit__(self, *_):
        self._release()

    def _acquire(self):
        ld = self._lock_dir
        # Check for stale lock
        if ld.exists():
            age = time.time() - ld.stat().st_mtime
            if age > self.STALE_SECONDS:
                warn(f"Removing stale lock ({age/3600:.1f}h old): {ld}")
                shutil.rmtree(ld, ignore_errors=True)
            else:
                raise RuntimeError(
                    f"Reference DB is locked by another process: {ld}\n"
                    f"Lock age: {age:.0f}s. If no process is running, "
                    f"delete manually: rm -rf {ld}"
                )
        ld.mkdir(parents=True)

    def _release(self):
        shutil.rmtree(self._lock_dir, ignore_errors=True)


# ─────────────────────────────────────────────────────────────────────────────
# File utilities
# ─────────────────────────────────────────────────────────────────────────────

def md5sum(path: Path, chunk: int = 1 << 20) -> str:
    h = hashlib.md5()
    with open(path, "rb") as fh:
        while True:
            block = fh.read(chunk)
            if not block:
                break
            h.update(block)
    return h.hexdigest()


def download_file(url: str, dest: Path, show_progress: bool = True) -> None:
    """Download a file, showing progress if on a tty."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".tmp")
    try:
        info(f"Downloading: {url}")
        info(f"         -> {dest}")
        with urllib.request.urlopen(url, timeout=60) as resp:
            total = int(resp.headers.get("Content-Length", 0))
            downloaded = 0
            with open(tmp, "wb") as fh:
                while True:
                    chunk = resp.read(1 << 16)
                    if not chunk:
                        break
                    fh.write(chunk)
                    downloaded += len(chunk)
                    if show_progress and total and sys.stderr.isatty():
                        pct = downloaded / total * 100
                        bar = int(pct / 2)
                        sys.stderr.write(
                            f"\r  [{'#'*bar}{' '*(50-bar)}] {pct:5.1f}%  "
                            f"{downloaded//1024//1024}MB/{total//1024//1024}MB"
                        )
                        sys.stderr.flush()
        if sys.stderr.isatty():
            sys.stderr.write("\n")
        tmp.rename(dest)
    except Exception:
        tmp.unlink(missing_ok=True)
        raise


def run_cmd(cmd: list[str], check: bool = True, **kwargs) -> subprocess.CompletedProcess:
    """Run a shell command, logging it, raising on failure."""
    info("CMD: " + " ".join(str(c) for c in cmd))
    result = subprocess.run(cmd, **kwargs)
    if check and result.returncode != 0:
        raise RuntimeError(
            f"Command failed (exit {result.returncode}): {' '.join(str(c) for c in cmd)}"
        )
    return result


def decompress_gz(src: Path, dst: Path) -> None:
    """Decompress .gz file."""
    import gzip
    info(f"Decompressing {src.name} ...")
    with gzip.open(src, "rb") as fi, open(dst, "wb") as fo:
        shutil.copyfileobj(fi, fo)


# ─────────────────────────────────────────────────────────────────────────────
# Index builders
# ─────────────────────────────────────────────────────────────────────────────

def _tool_available(name: str) -> bool:
    return shutil.which(name) is not None


def build_hisat2_index(genome_fa: Path, index_dir: Path,
                       gtf: Path, with_dta: bool = False) -> None:
    index_dir.mkdir(parents=True, exist_ok=True)
    index_base = index_dir / "genome"

    # Extract splice sites and exons from GTF
    ss_file   = index_dir / "splice_sites.txt"
    exon_file = index_dir / "exons.txt"

    if _tool_available("hisat2_extract_splice_sites.py"):
        run_cmd(["hisat2_extract_splice_sites.py", str(gtf)],
                stdout=open(ss_file, "w"))
        run_cmd(["hisat2_extract_exons.py", str(gtf)],
                stdout=open(exon_file, "w"))
        extra = ["--known-splicesite-infile", str(ss_file),
                 "--exon-infile", str(exon_file)]
    else:
        warn("hisat2_extract_splice_sites.py not found — building index without splice sites.")
        extra = []

    cmd = ["hisat2-build", "-p", "8",
           "--seed", "42"] + extra + [str(genome_fa), str(index_base)]
    run_cmd(cmd)

    # Log build parameters
    params = {"dta": with_dta, "seed": 42}
    with open(index_dir / "build_params.yaml", "w") as fh:
        yaml.dump(params, fh)

    success(f"HISAT2 index built: {index_base}")


def build_bwa_index(genome_fa: Path, index_dir: Path) -> None:
    index_dir.mkdir(parents=True, exist_ok=True)
    # BWA mem wants the FASTA in the index dir (symlink is fine)
    dest_fa = index_dir / genome_fa.name
    if not dest_fa.exists():
        dest_fa.symlink_to(genome_fa.resolve())
    run_cmd(["bwa", "index", str(dest_fa)])
    success(f"BWA index built: {dest_fa}")


def build_bowtie2_index(fasta: Path, index_dir: Path, tag: str) -> None:
    index_dir.mkdir(parents=True, exist_ok=True)
    index_base = index_dir / tag
    run_cmd(["bowtie2-build", "--seed", "42", "-p", "8",
             str(fasta), str(index_base)])
    success(f"Bowtie2 index built: {index_base}")


def build_samtools_fai(fasta: Path) -> None:
    run_cmd(["samtools", "faidx", str(fasta)])


# ─────────────────────────────────────────────────────────────────────────────
# gene_type_map extractor
# ─────────────────────────────────────────────────────────────────────────────

def extract_gene_type_map(gtf: Path, out: Path) -> None:
    """
    Parse GTF and write gene_id -> gene_type mapping.
    Used by build_integrated_matrix.py to assign category (canonical vs lncRNA).
    """
    info(f"Extracting gene_type map from {gtf.name} ...")
    seen: dict[str, dict] = {}
    opener = open

    # Handle gzipped GTF
    import gzip as _gzip
    if gtf.suffix == ".gz":
        opener = _gzip.open

    with opener(gtf, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "gene":
                continue
            attrs = _parse_gtf_attrs(parts[8])
            gid   = attrs.get("gene_id", "")
            gtype = attrs.get("gene_type") or attrs.get("gene_biotype", "")
            gname = attrs.get("gene_name", "")
            if gid and gid not in seen:
                seen[gid] = {
                    "gene_id":   gid,
                    "gene_name": gname,
                    "gene_type": gtype,
                    "chr":       parts[0],
                    "start":     parts[3],
                    "end":       parts[4],
                    "strand":    parts[6],
                }

    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=["gene_id","gene_name","gene_type","chr","start","end","strand"],
            delimiter="\t"
        )
        writer.writeheader()
        for row in seen.values():
            writer.writerow(row)

    info(f"  -> {len(seen)} genes written to {out.name}")


def _parse_gtf_attrs(attr_str: str) -> dict[str, str]:
    """Parse GTF attribute string into dict."""
    d: dict[str, str] = {}
    for part in attr_str.strip().rstrip(";").split(";"):
        part = part.strip()
        if not part:
            continue
        sp = part.split(None, 1)
        if len(sp) == 2:
            d[sp[0]] = sp[1].strip('"')
    return d


# ─────────────────────────────────────────────────────────────────────────────
# Reference presence check
# ─────────────────────────────────────────────────────────────────────────────

def is_reference_present(db_root: Path, organism: str, build: str) -> bool:
    """
    Returns True if the reference for this organism/build appears complete
    (manifest.yaml exists and genome FASTA is present).
    """
    m = read_manifest(db_root, organism, build)
    if not m:
        return False
    p = org_path(db_root, organism, build)
    genome_dir = p / "genome"
    if not any(genome_dir.glob("*.fa")):
        return False
    return True


# ─────────────────────────────────────────────────────────────────────────────
# Core download / build flow
# ─────────────────────────────────────────────────────────────────────────────

class ReferenceBuilder:
    """
    Orchestrates download and index building for one organism/build.
    Used both by mosaic-refs CLI and by the pipeline pre-flight.
    """

    def __init__(self, db_root: Path, organism: str, build: str):
        self.db_root  = Path(db_root)
        self.organism = organism
        self.build    = build
        self.base     = org_path(self.db_root, organism, build)
        self.spec     = _ORGANISM_CATALOGUE.get(organism)
        self.log      = get_logger()

    def _dirs(self) -> dict[str, Path]:
        b = self.base
        return {
            "genome":           b / "genome",
            "annotation":       b / "annotation",
            "hisat2_index":     b / "hisat2_index",
            "bwa_index":        b / "bwa_index",
            "dfam":             b / "dfam",
            "bowtie2_TE":       b / "dfam" / "bowtie2_index_TE",
            "bowtie2_sat":      b / "dfam" / "bowtie2_index_satellite",
        }

    def download_genome(self) -> Path:
        dirs = self._dirs()
        dirs["genome"].mkdir(parents=True, exist_ok=True)
        fa_gz = dirs["genome"] / f"{self.organism}_{self.build}.fa.gz"
        fa    = dirs["genome"] / f"{self.organism}_{self.build}.fa"

        if fa.exists():
            info(f"Genome FASTA already present: {fa.name}")
            return fa

        if self.spec is None:
            raise RuntimeError(
                f"Unknown organism {self.organism!r}. "
                "Use 'mosaic-refs add' to register a custom organism."
            )

        # Prefer GENCODE for human/mouse, fall back to Ensembl
        url = (self.spec.gencode_fasta_url or self.spec.ensembl_fasta_url)
        download_file(url, fa_gz)
        decompress_gz(fa_gz, fa)
        fa_gz.unlink()

        build_samtools_fai(fa)

        m5 = md5sum(fa)
        append_update_log(self.db_root, self.organism, self.build,
                          "genome", "--", self.spec.ensembl_release,
                          "download", f"md5={m5}")
        return fa

    def download_annotation(self) -> Path:
        dirs = self._dirs()
        dirs["annotation"].mkdir(parents=True, exist_ok=True)

        if self.spec is None:
            raise RuntimeError(f"Unknown organism {self.organism!r}.")

        src = self.spec.gencode_release or self.spec.ensembl_release
        url = (self.spec.gencode_gtf_url or self.spec.ensembl_gtf_url)
        gtf_gz = dirs["annotation"] / f"{self.organism}_{self.build}_{src}.gtf.gz"
        gtf    = dirs["annotation"] / f"{self.organism}_{self.build}_{src}.gtf"
        gmap   = dirs["annotation"] / "gene_type_map.tsv"

        if gtf.exists() and gmap.exists():
            info(f"GTF annotation already present: {gtf.name}")
            return gtf

        download_file(url, gtf_gz)
        decompress_gz(gtf_gz, gtf)
        gtf_gz.unlink()

        extract_gene_type_map(gtf, gmap)

        m5 = md5sum(gtf)
        append_update_log(self.db_root, self.organism, self.build,
                          "annotation", "--", src,
                          "download", f"md5={m5}")
        return gtf

    def build_indices(self, genome_fa: Path, gtf: Path,
                      with_dta: bool = False) -> None:
        dirs = self._dirs()

        # HISAT2
        if not any((dirs["hisat2_index"]).glob("genome.*.ht2")):
            build_hisat2_index(genome_fa, dirs["hisat2_index"], gtf, with_dta)
            append_update_log(self.db_root, self.organism, self.build,
                              "hisat2_index", "--", "built",
                              "build_index", f"dta={with_dta}")
        else:
            info("HISAT2 index already present.")

        # BWA
        if not (dirs["bwa_index"] / f"{self.organism}_{self.build}.fa.bwt").exists():
            # symlink genome FA into bwa_index dir
            dirs["bwa_index"].mkdir(parents=True, exist_ok=True)
            link = dirs["bwa_index"] / genome_fa.name
            if not link.exists():
                link.symlink_to(genome_fa.resolve())
            build_bwa_index(genome_fa, dirs["bwa_index"])
            append_update_log(self.db_root, self.organism, self.build,
                              "bwa_index", "--", "built", "build_index", "")
        else:
            info("BWA index already present.")

    def write_manifest(self, genome_fa: Path, gtf: Path) -> None:
        if self.spec is None:
            src = "custom"
            annot_ver = "unknown"
        else:
            src = self.spec.annotation_source
            annot_ver = self.spec.gencode_release or self.spec.ensembl_release

        data = {
            "organism":     self.organism,
            "genome_build": self.build,
            "last_updated": datetime.now().strftime("%Y-%m-%d"),
            "genome": {
                "source":          getattr(self.spec, "annotation_source", "Ensembl"),
                "release":         getattr(self.spec, "ensembl_release", "unknown"),
                "file":            str(genome_fa.relative_to(self.base)),
                "md5":             md5sum(genome_fa),
                "date_downloaded": datetime.now().strftime("%Y-%m-%d"),
            },
            "annotation": {
                "source":          src,
                "release":         annot_ver,
                "file":            str(gtf.relative_to(self.base)),
                "md5":             md5sum(gtf),
                "date_downloaded": datetime.now().strftime("%Y-%m-%d"),
            },
            "dfam": {
                "version":     "pending",
                "n_te_families":        0,
                "n_satellite_families": 0,
            },
            "indices": {
                "hisat2_dta":    False,
                "hisat2_version": _tool_version("hisat2"),
                "bowtie2_version": _tool_version("bowtie2"),
                "bwa_version":   _tool_version("bwa"),
            },
        }
        write_manifest(self.db_root, self.organism, self.build, data)
        success(f"Manifest written: {manifest_path(self.db_root, self.organism, self.build)}")

    def full_download_and_index(self, with_dta: bool = False) -> None:
        """Full pipeline: genome + annotation + indices + manifest."""
        section(f"Building reference: {self.organism} / {self.build}")
        with RefLock(self.db_root, self.organism, self.build):
            genome_fa = self.download_genome()
            gtf       = self.download_annotation()
            self.build_indices(genome_fa, gtf, with_dta)
            self.write_manifest(genome_fa, gtf)
        success(f"Reference ready: {self.base}")


def _tool_version(name: str) -> str:
    if not shutil.which(name):
        return "not_found"
    try:
        r = subprocess.run([name, "--version"], capture_output=True, text=True, timeout=5)
        first_line = (r.stdout or r.stderr or "").splitlines()
        return first_line[0].strip() if first_line else "unknown"
    except Exception:
        return "unknown"


# ─────────────────────────────────────────────────────────────────────────────
# Pre-flight check (called from M00 launcher)
# ─────────────────────────────────────────────────────────────────────────────

def preflight_reference_check(db_root: Path, organism: str, build: str,
                               auto_update: bool = False,
                               project_wd: Optional[Path] = None) -> None:
    """
    M00 pre-flight: check reference, auto-download if missing,
    optionally prompt for update, copy snapshot to project WD.
    """
    section(f"Reference DB check: {organism} / {build}")
    db_root = Path(db_root)
    present = is_reference_present(db_root, organism, build)

    if not present:
        info(f"Reference for '{organism}/{build}' not found in {db_root}.")
        info("Starting automatic download and index build...")
        builder = ReferenceBuilder(db_root, organism, build)
        builder.full_download_and_index()
    else:
        manifest = read_manifest(db_root, organism, build)
        info(f"Reference present: {organism}/{build}")
        info(f"  Genome:     {manifest.get('genome',{}).get('source','')} "
             f"release {manifest.get('genome',{}).get('release','?')}")
        info(f"  Annotation: {manifest.get('annotation',{}).get('source','')} "
             f"release {manifest.get('annotation',{}).get('release','?')}")
        info(f"  Last updated: {manifest.get('last_updated','?')}")

        if auto_update:
            info("auto_update=true: updating reference...")
            builder = ReferenceBuilder(db_root, organism, build)
            builder.full_download_and_index()

        elif sys.stdout.isatty():
            # Interactive terminal — ask user
            try:
                ans = input("\nUpdate reference DB? [y/N]: ").strip().lower()
            except (EOFError, KeyboardInterrupt):
                ans = "n"
            if ans in ("y", "yes"):
                builder = ReferenceBuilder(db_root, organism, build)
                builder.full_download_and_index()
            else:
                info("Using existing reference (no update).")
        else:
            # HPC / non-interactive — use existing silently
            info("HPC mode: using existing reference without update.")

    # Copy snapshot to project WD
    if project_wd is not None:
        project_wd = Path(project_wd)
        project_wd.mkdir(parents=True, exist_ok=True)
        snap_src = manifest_path(db_root, organism, build)
        snap_dst = project_wd / "reference_snapshot.yaml"
        if snap_src.exists() and not snap_dst.exists():
            shutil.copy2(snap_src, snap_dst)
            info(f"Reference snapshot saved: {snap_dst}")
        elif snap_dst.exists():
            info(f"Reference snapshot already exists (immutable): {snap_dst}")


# ─────────────────────────────────────────────────────────────────────────────
# mosaic-refs CLI helpers (called from cli.py)
# ─────────────────────────────────────────────────────────────────────────────

def list_references(db_root: Path) -> list[dict]:
    """Return list of available organisms with version info."""
    db_root = Path(db_root)
    results = []
    if not db_root.exists():
        return results
    for org_dir in sorted(db_root.iterdir()):
        if not org_dir.is_dir():
            continue
        for build_dir in sorted(org_dir.iterdir()):
            if not build_dir.is_dir():
                continue
            m = {}
            mp = build_dir / "manifest.yaml"
            if mp.exists():
                with open(mp) as fh:
                    m = yaml.safe_load(fh) or {}
            results.append({
                "organism":     org_dir.name,
                "build":        build_dir.name,
                "last_updated": m.get("last_updated", "?"),
                "genome_release": m.get("genome", {}).get("release", "?"),
                "annotation_release": m.get("annotation", {}).get("release", "?"),
                "dfam_version": m.get("dfam", {}).get("version", "?"),
                "path":         str(build_dir),
            })
    return results


def check_integrity(db_root: Path, organism: str, build: str) -> bool:
    """
    Verify MD5 checksums for all files listed in manifest.
    Returns True if all pass.
    """
    m = read_manifest(db_root, organism, build)
    if not m:
        error(f"No manifest found for {organism}/{build}.")
        return False

    base = org_path(db_root, organism, build)
    all_ok = True

    for component in ("genome", "annotation"):
        comp = m.get(component, {})
        rel_path = comp.get("file", "")
        expected_md5 = comp.get("md5", "")
        if not rel_path or not expected_md5:
            continue
        fp = base / rel_path
        if not fp.exists():
            error(f"  MISSING: {fp}")
            all_ok = False
            continue
        actual = md5sum(fp)
        if actual == expected_md5:
            success(f"  OK  {component}: {fp.name}")
        else:
            error(f"  FAIL {component}: {fp.name} (expected {expected_md5}, got {actual})")
            all_ok = False

    return all_ok


def read_update_log(db_root: Path, organism: str, build: str,
                    since: Optional[str] = None) -> list[dict]:
    """Read update log entries, optionally filtered by date."""
    p = update_log_path(db_root, organism, build)
    if not p.exists():
        return []
    rows = []
    with open(p) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            if since and row.get("timestamp", "") < since:
                continue
            rows.append(row)
    return rows
