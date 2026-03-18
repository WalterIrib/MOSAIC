"""
mosaic/db_manager.py
MOSAIC Database Manager.

Handles download, verification, and version management of:
- Genome FASTA + GTF (GENCODE/Ensembl)
- STAR index
- DFAM TE consensus / pseudochromosomes
- RepeatMasker (UCSC rmsk)
- Kraken2 database
- MetaPhlAn database
- HUMAnN3 databases (ChocoPhlAn, UniRef)
- Arriba databases
- STRINGdb (via R package, auto)

All downloads are tracked in a manifest.json for reproducibility.
"""
from __future__ import annotations

import hashlib
import json
import os
import subprocess
import sys
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional

from .log import info, warn, error, success, section


# ─────────────────────────────────────────────────────────────────────────────
# Database catalogue
# ─────────────────────────────────────────────────────────────────────────────

DB_CATALOGUE: Dict[str, Dict] = {
    # Genome + annotation
    "genome_hg38": {
        "description": "GRCh38 primary assembly (GENCODE v49)",
        "url": "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/GRCh38.primary_assembly.genome.fa.gz",
        "type": "fasta",
        "organism": "human",
    },
    "gtf_hg38_v49": {
        "description": "GENCODE v49 annotation GTF",
        "url": "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.annotation.gtf.gz",
        "type": "gtf",
        "organism": "human",
    },
    "genome_mm39": {
        "description": "GRCm39 primary assembly (GENCODE M36)",
        "url": "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M36/GRCm39.primary_assembly.genome.fa.gz",
        "type": "fasta",
        "organism": "mouse",
    },
    "gtf_mm39_M36": {
        "description": "GENCODE M36 mouse annotation",
        "url": "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M36/gencode.vM36.annotation.gtf.gz",
        "type": "gtf",
        "organism": "mouse",
    },
    # Kraken2
    "kraken2_standard": {
        "description": "Kraken2 standard database (bacteria+virus+archaea+fungi)",
        "url": "https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20240904.tar.gz",
        "type": "kraken2",
        "size_gb": 70,
    },
    "kraken2_standard_8gb": {
        "description": "Kraken2 standard 8GB (reduced)",
        "url": "https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08gb_20240904.tar.gz",
        "type": "kraken2",
        "size_gb": 8,
    },
    "kraken2_minusb": {
        "description": "Kraken2 minus bacteria (viruses+archaea+fungi)",
        "url": "https://genome-idx.s3.amazonaws.com/kraken/k2_minusb_20240904.tar.gz",
        "type": "kraken2",
        "size_gb": 8,
    },
    # Arriba
    "arriba_v2.4.0": {
        "description": "Arriba fusion detection v2.4.0 (binary + databases)",
        "url": "https://github.com/suhrig/arriba/releases/download/v2.4.0/arriba_v2.4.0.tar.gz",
        "type": "arriba",
        "size_gb": 0.5,
    },
}

# Per-organism genome build defaults
GENOME_DEFAULTS: Dict[str, Dict] = {
    "human": {
        "build":    "hg38",
        "gencode":  "49",
        "genome":   "genome_hg38",
        "gtf":      "gtf_hg38_v49",
        "kegg":     "hsa",
        "string_taxon": 9606,
    },
    "mouse": {
        "build":    "mm39",
        "gencode":  "M36",
        "genome":   "genome_mm39",
        "gtf":      "gtf_mm39_M36",
        "kegg":     "mmu",
        "string_taxon": 10090,
    },
}


# ─────────────────────────────────────────────────────────────────────────────
# Utilities
# ─────────────────────────────────────────────────────────────────────────────

def md5sum(path: Path) -> str:
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            h.update(chunk)
    return h.hexdigest()


def wget(url: str, outpath: Path, show_progress: bool = True) -> bool:
    cmd = ["wget", "-c", "-O", str(outpath), url]
    if not show_progress:
        cmd.insert(1, "-q")
    ret = subprocess.run(cmd).returncode
    return ret == 0


def decompress(src: Path, dst: Optional[Path] = None) -> Path:
    """Decompress .gz or .tar.gz file."""
    if src.name.endswith(".tar.gz") or src.name.endswith(".tgz"):
        dst = dst or src.parent
        subprocess.run(["tar", "-xzf", str(src), "-C", str(dst)], check=True)
        return dst
    elif src.name.endswith(".gz"):
        dst = dst or src.with_suffix("")
        subprocess.run(["gunzip", "-c", str(src)], stdout=open(dst, "wb"),
                       check=True)
        return dst
    return src


def run_r(code: str) -> bool:
    ret = subprocess.run(["Rscript", "-e", code], capture_output=True)
    return ret.returncode == 0


# ─────────────────────────────────────────────────────────────────────────────
# Manifest (version tracking)
# ─────────────────────────────────────────────────────────────────────────────

class Manifest:
    """Tracks downloaded database versions for reproducibility."""

    def __init__(self, db_root: Path):
        self.path = db_root / "manifest.json"
        self.data: Dict = {}
        self._load()

    def _load(self):
        if self.path.exists():
            try:
                self.data = json.loads(self.path.read_text())
            except Exception:
                self.data = {}

    def save(self):
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self.path.write_text(json.dumps(self.data, indent=2))

    def record(self, key: str, version: str, url: str,
               path: str, checksum: Optional[str] = None):
        self.data[key] = {
            "version":       version,
            "url":           url,
            "local_path":    path,
            "checksum_md5":  checksum,
            "date_downloaded": datetime.now().strftime("%Y-%m-%d"),
        }
        self.save()

    def get(self, key: str) -> Optional[Dict]:
        return self.data.get(key)

    def is_present(self, key: str) -> bool:
        entry = self.data.get(key)
        if not entry:
            return False
        p = Path(entry["local_path"])
        return p.exists()


# ─────────────────────────────────────────────────────────────────────────────
# Database Manager
# ─────────────────────────────────────────────────────────────────────────────

class DatabaseManager:
    """
    Manages all MOSAIC reference databases.
    Handles download, verification, and version tracking.
    """

    def __init__(self, cfg: Any, organism: str):
        self.cfg      = cfg
        self.organism = organism
        self.db_root  = Path(cfg.databases.db_root
                              if hasattr(cfg.databases, "db_root")
                              else str(Path(cfg.project.wd) / "references"))
        self.db_root.mkdir(parents=True, exist_ok=True)
        self.manifest = Manifest(self.db_root)
        self.org_info = GENOME_DEFAULTS.get(organism, {})

    def check_all(self) -> Dict[str, bool]:
        """Check presence of all required databases."""
        db = self.cfg.databases
        # Only check path-type fields (skip booleans like auto_download)
        path_fields = {
            "genome_fa", "gtf", "star_index", "te_consensus",
            "te_satellite", "rmsk_bed", "te_gtf", "kraken2_db",
            "metaphlan_db", "chocophlan_db", "uniref_db",
            "arriba_bin", "arriba_db",
        }
        required = {
            k: Path(getattr(db, k, "") or "")
            for k in path_fields
        }

        status = {}
        for name, path in required.items():
            p_str = str(path)
            if not p_str or p_str in (".", "/path/to/references"):
                status[name] = None  # not configured = skip
                continue
            if path.exists():
                status[name] = True
                info(f"  ✓ {name}: {path}")
            else:
                status[name] = False
                warn(f"  ✗ {name}: MISSING ({path})")

        return status

    def download(self, component: str = "all", version: str = "latest"):
        """Download specified component(s)."""
        section(f"Downloading: {component} (version: {version})")

        if component == "all":
            self.download_genome(version)
            self.download_star_index()
            self.download_kraken2(version)
            self.download_arriba(version)
            # HUMAnN3 handled separately (different env)
            info("Note: HUMAnN3 databases (MetaPhlAn/ChocoPhlAn/UniRef) require")
            info("      the humann3 conda env. Run:")
            info("      conda activate humann3 && mosaic db download --component metaphlan")
        elif component == "genome":
            self.download_genome(version)
        elif component == "star_index":
            self.download_star_index()
        elif component == "kraken2":
            self.download_kraken2(version)
        elif component == "metaphlan":
            self.download_metaphlan(version)
        elif component == "chocophlan":
            self.download_chocophlan(version)
        elif component == "uniref":
            self.download_uniref(version)
        elif component == "arriba":
            self.download_arriba(version)
        elif component == "dfam":
            self.download_dfam(version)
        else:
            error(f"Unknown component: {component}")

    def download_genome(self, version: str = "latest"):
        """Download genome FASTA and GTF."""
        org_key = self.org_info.get("genome", "genome_hg38")
        gtf_key = self.org_info.get("gtf", "gtf_hg38_v49")

        for key in [org_key, gtf_key]:
            entry = DB_CATALOGUE.get(key, {})
            if not entry:
                warn(f"No catalogue entry for {key}")
                continue

            url  = entry["url"]
            name = key.replace("/", "_")
            ext  = ".fa.gz" if entry["type"] == "fasta" else ".gtf.gz"
            dest = self.db_root / "genome" / f"{name}{ext}"
            dest.parent.mkdir(exist_ok=True)

            if dest.exists() or dest.with_suffix("").exists():
                info(f"Already exists: {dest.stem}")
                continue

            info(f"Downloading {key}...")
            if wget(url, dest):
                out = decompress(dest)
                # Index FASTA
                if entry["type"] == "fasta" and out.suffix == ".fa":
                    subprocess.run(["samtools", "faidx", str(out)])
                self.manifest.record(key, version, url, str(out))
                success(f"Downloaded: {out}")
            else:
                error(f"Failed to download {key}")

    def download_star_index(self):
        """Build STAR index from genome + GTF."""
        star_dir = self.db_root / "star_index"
        if star_dir.exists() and list(star_dir.glob("*.bin")):
            info("STAR index already exists")
            return

        genome_fa = self._find_genome_fa()
        gtf       = self._find_gtf()
        if not genome_fa or not gtf:
            error("Genome FASTA or GTF not found. Download genome first.")
            return

        star_dir.mkdir(exist_ok=True)
        read_len = getattr(self.cfg.sequencing, "read_length", 150)
        overhang = read_len - 1

        info("Building STAR index (requires ~30GB RAM, ~1 hour)...")
        cmd = [
            "bash", "-c",
            f"ulimit -n 65536 && STAR "
            f"--runMode genomeGenerate "
            f"--genomeDir {star_dir} "
            f"--genomeFastaFiles {genome_fa} "
            f"--sjdbGTFfile {gtf} "
            f"--sjdbOverhang {overhang} "
            f"--runThreadN {getattr(self.cfg, 'threads', 20)} "
            f"--limitGenomeGenerateRAM 27000000000"
        ]
        ret = subprocess.run(" ".join(cmd), shell=True).returncode
        if ret == 0:
            self.manifest.record("star_index", "built", "local", str(star_dir))
            success(f"STAR index built: {star_dir}")
        else:
            error("STAR index build failed")

    def download_kraken2(self, version: str = "standard"):
        """Download Kraken2 database."""
        db_key = f"kraken2_{version}" if version in DB_CATALOGUE else "kraken2_standard"
        entry  = DB_CATALOGUE.get(db_key, DB_CATALOGUE["kraken2_standard"])

        kr_dir = self.db_root / "kraken2_db"
        if kr_dir.exists() and (kr_dir / "hash.k2d").exists():
            info(f"Kraken2 database already exists: {kr_dir}")
            return

        kr_dir.mkdir(exist_ok=True)
        url     = entry["url"]
        tarball = kr_dir / "kraken2.tar.gz"

        size_gb = entry.get("size_gb", 70)
        info(f"Downloading Kraken2 database (~{size_gb}GB)...")
        if wget(url, tarball):
            decompress(tarball, kr_dir)
            tarball.unlink()
            self.manifest.record("kraken2", version, url, str(kr_dir))
            success(f"Kraken2 database ready: {kr_dir}")
        else:
            error("Kraken2 download failed")

    def download_metaphlan(self, version: str = "latest"):
        """Download MetaPhlAn database (in humann3 env)."""
        meta_dir = self.db_root / "humann3_db" / "metaphlan_db"
        meta_dir.mkdir(parents=True, exist_ok=True)

        # Check if already installed
        existing = list(meta_dir.glob("*.bt2l"))
        if existing:
            info(f"MetaPhlAn database already exists: {meta_dir}")
            return

        info("Installing MetaPhlAn database...")
        ret = subprocess.run([
            "conda", "run", "-n", "humann3", "--no-capture-output",
            "metaphlan", "--install",
            "--db_dir", str(meta_dir),
            "--index", version if version != "latest" else "latest",
        ]).returncode

        if ret == 0:
            self.manifest.record("metaphlan", version, "metaphlan_install",
                                  str(meta_dir))
            success(f"MetaPhlAn database installed: {meta_dir}")
            # Apply humann3 config fix
            self._fix_humann3_config()
        else:
            error("MetaPhlAn database installation failed")

    def download_chocophlan(self, version: str = "full"):
        """Download HUMAnN3 ChocoPhlAn nucleotide database."""
        cho_dir = self.db_root / "humann3_db" / "chocophlan"
        if cho_dir.exists() and list(cho_dir.glob("*.ffn.gz")):
            info("ChocoPhlAn database already exists")
            return

        info(f"Downloading ChocoPhlAn ({version}, ~15GB)...")
        ret = subprocess.run([
            "conda", "run", "-n", "humann3", "--no-capture-output",
            "humann_databases",
            "--download", "chocophlan", version,
            str(self.db_root / "humann3_db"),
        ]).returncode

        if ret == 0:
            self.manifest.record("chocophlan", version,
                                  "humann_databases", str(cho_dir))
            success(f"ChocoPhlAn ready: {cho_dir}")
        else:
            error("ChocoPhlAn download failed")

    def download_uniref(self, version: str = "uniref90_diamond"):
        """Download HUMAnN3 UniRef protein database."""
        uni_dir = self.db_root / "humann3_db" / "uniref"
        if uni_dir.exists() and list(uni_dir.glob("*.dmnd")):
            info("UniRef database already exists")
            return

        info(f"Downloading UniRef90 diamond (~20GB)...")
        ret = subprocess.run([
            "conda", "run", "-n", "humann3", "--no-capture-output",
            "humann_databases",
            "--download", "uniref", version,
            str(self.db_root / "humann3_db"),
        ]).returncode

        if ret == 0:
            self.manifest.record("uniref", version,
                                  "humann_databases", str(uni_dir))
            success(f"UniRef ready: {uni_dir}")
        else:
            error("UniRef download failed")

    def download_arriba(self, version: str = "2.4.0"):
        """Download Arriba fusion detection binary and databases."""
        arriba_dir = self.db_root / f"arriba_v{version}"
        if arriba_dir.exists() and (arriba_dir / "arriba").exists():
            info(f"Arriba already exists: {arriba_dir}")
            return

        url     = DB_CATALOGUE.get(f"arriba_v{version}", {}).get(
            "url",
            f"https://github.com/suhrig/arriba/releases/download/"
            f"v{version}/arriba_v{version}.tar.gz"
        )
        tarball = self.db_root / f"arriba_v{version}.tar.gz"

        info(f"Downloading Arriba v{version}...")
        if wget(url, tarball):
            decompress(tarball, self.db_root)
            tarball.unlink()
            self.manifest.record("arriba", version, url, str(arriba_dir))
            success(f"Arriba ready: {arriba_dir}")
        else:
            error("Arriba download failed")

    def download_dfam(self, version: str = "3.9"):
        """Download DFAM TE consensus sequences."""
        dfam_dir = self.db_root / "dfam"
        dfam_dir.mkdir(exist_ok=True)

        info(f"Downloading DFAM {version} TE families...")
        # DFAM requires building pseudochromosomes from Dfam.h5
        # This is handled by build_pseudochromosomes.py
        url = (f"https://www.dfam.org/releases/Dfam_{version}/"
               f"families/Dfam_curatedonly.h5.gz")
        dest = dfam_dir / f"Dfam_{version}_curatedonly.h5.gz"

        if wget(url, dest):
            out = decompress(dest)
            self.manifest.record("dfam", version, url, str(out))
            info(f"DFAM downloaded: {out}")
            info("Run build_pseudochromosomes.py to generate TE consensus FASTA")
        else:
            error("DFAM download failed")

    def list_versions(self):
        """List all available database versions from manifest."""
        section("Database Manifest")
        if not self.manifest.data:
            warn("No databases recorded in manifest.")
            return

        fmt = "{:<20} {:<12} {:<25} {:<12}"
        print(fmt.format("COMPONENT", "VERSION", "PATH", "DATE"))
        print("─" * 72)
        for key, entry in sorted(self.manifest.data.items()):
            p = Path(entry["local_path"])
            exists = "✓" if p.exists() else "✗"
            print(fmt.format(
                key[:20],
                entry.get("version", "?")[:12],
                f"{exists} {str(p)[:24]}",
                entry.get("date_downloaded", "?")
            ))

    # ── Private helpers ───────────────────────────────────────────────────────

    def _find_genome_fa(self) -> Optional[Path]:
        p = getattr(self.cfg.databases, "genome_fa", None)
        if p and Path(p).exists():
            return Path(p)
        # Search db_root
        for fa in (self.db_root / "genome").glob("*.fa"):
            return fa
        return None

    def _find_gtf(self) -> Optional[Path]:
        p = getattr(self.cfg.databases, "gtf", None)
        if p and Path(p).exists():
            return Path(p)
        for gtf in (self.db_root / "genome").glob("*.gtf"):
            return gtf
        return None

    def _fix_humann3_config(self):
        """
        Apply required fixes to HUMAnN3 config for MetaPhlAn v4 compatibility.
        - Set metaphlan major version to 4
        - Set metaphlan line parser to 0
        - Update metaphlan_v4_db_version to match installed db
        - Create metaphlan wrapper script to pass --db_dir automatically
        """
        import glob
        info("Applying HUMAnN3/MetaPhlAn v4 compatibility fixes...")

        # Find humann config
        result = subprocess.run(
            ["conda", "run", "-n", "humann3",
             "python", "-c",
             "import humann; import os; "
             "print(os.path.join(os.path.dirname(humann.__file__), 'config.py'))"],
            capture_output=True, text=True
        )
        if result.returncode != 0:
            warn("Could not locate humann config.py")
            return

        cfg_path = result.stdout.strip()
        if not os.path.exists(cfg_path):
            warn(f"humann config.py not found: {cfg_path}")
            return

        # Read and patch
        content = open(cfg_path).read()
        content = content.replace('"major" : 3', '"major" : 4')
        content = content.replace('"line" : -1', '"line" : 0')

        # Detect installed MetaPhlAn db version
        meta_dir = self.db_root / "humann3_db" / "metaphlan_db"
        mpa_files = list(meta_dir.glob("mpa_v*.pkl"))
        if mpa_files:
            db_ver = mpa_files[0].stem.split("_CHOCOPhl")[0].replace(
                "mpa_", "")  # e.g. "vJan25"
            content = content.replace(
                'metaphlan_v4_db_version="vJun23"',
                f'metaphlan_v4_db_version="{db_ver}"'
            )
            info(f"  MetaPhlAn db version set to: {db_ver}")

        open(cfg_path, "w").write(content)
        info("  humann config.py patched")

        # Create metaphlan wrapper
        result2 = subprocess.run(
            ["conda", "run", "-n", "humann3",
             "python", "-c", "import sys; print(sys.prefix)"],
            capture_output=True, text=True
        )
        if result2.returncode == 0:
            env_bin = Path(result2.stdout.strip()) / "bin"
            orig    = env_bin / "metaphlan_orig"
            wrapper = env_bin / "metaphlan"

            if not orig.exists() and wrapper.exists():
                # Backup original
                import shutil
                shutil.copy(wrapper, orig)

            # Write wrapper
            wrapper_content = f"""#!/bin/bash
# MOSAIC metaphlan wrapper - translates MetaPhlAn 3 args to v4
ARGS=()
for arg in "$@"; do
    case "$arg" in
        --bowtie2out) ARGS+=("--mapout") ;;
        --bowtie2db)  ARGS+=("--db_dir") ;;
        *)            ARGS+=("$arg") ;;
    esac
done
exec {orig} \\
  --db_dir {meta_dir} \\
  "${{ARGS[@]}}"
"""
            wrapper.write_text(wrapper_content)
            wrapper.chmod(0o755)
            info(f"  metaphlan wrapper created: {wrapper}")

        success("HUMAnN3/MetaPhlAn compatibility fixes applied")
