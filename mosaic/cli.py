"""
mosaic/cli.py  (v2)
CLI entrypoints for MOSAIC.

Commands:
  mosaic init        — create new project
  mosaic check       — validate config
  mosaic preflight   — check + download missing references
  mosaic run         — run pipeline (all or specific modules)
  mosaic status      — show pipeline state
  mosaic download    — download FASTQs
  mosaic db          — database management subgroup
  mosaic compare     — run a new comparison on existing alignments
"""
from __future__ import annotations

import sys
from pathlib import Path
from typing import Optional

import click

from .config import load_config, generate_default_config, ConfigError
from .log import setup_logging, section, info, warn, error, success
from . import __version__


# ─────────────────────────────────────────────────────────────────────────────
# Main group
# ─────────────────────────────────────────────────────────────────────────────

@click.group(name="mosaic")
@click.version_option(__version__, prog_name="mosaic")
def cli() -> None:
    """MOSAIC — Multi-Omic Analysis of Small-RNA, AI-assisted Integration
    and Characterization.

    A modular RNA-seq analysis pipeline with support for TE biology,
    metatranscriptomics, RNA editing, fusion genes, and more.
    """
    pass


# ─────────────────────────────────────────────────────────────────────────────
# mosaic init
# ─────────────────────────────────────────────────────────────────────────────

@cli.command("init")
@click.argument("project_dir", type=click.Path())
@click.option("--organism", "-o", default="human",
              help="Organism (human|mouse|rat|zebrafish|custom)")
@click.option("--context", "-c", default="tumor",
              type=click.Choice(["tumor", "aging", "congenital", "other"]),
              help="Biological context")
@click.option("--layout", "-l", default="PE",
              type=click.Choice(["PE", "SE"]),
              help="Library layout: PE (paired-end) or SE (single-end)")
def cmd_init(project_dir: str, organism: str,
             context: str, layout: str) -> None:
    """Initialise a new MOSAIC project directory with a template config."""
    setup_logging()
    p = Path(project_dir)
    p.mkdir(parents=True, exist_ok=True)
    cfg_path = p / "config.yaml"
    if cfg_path.exists():
        warn(f"Config already exists: {cfg_path}. Not overwriting.")
    else:
        cfg_path.write_text(generate_default_config(organism, context, layout))
        success(f"Config written: {cfg_path}")
    info(f"Next steps:")
    info(f"  1. Edit {cfg_path}")
    info(f"  2. mosaic preflight {cfg_path}")
    info(f"  3. mosaic run {cfg_path}")


# ─────────────────────────────────────────────────────────────────────────────
# mosaic check
# ─────────────────────────────────────────────────────────────────────────────

@cli.command("check")
@click.argument("config_file", type=click.Path(exists=True))
def cmd_check(config_file: str) -> None:
    """Validate a MOSAIC config file."""
    setup_logging()
    try:
        cfg = load_config(config_file)
        success("Config is valid.")
        info(f"  Organism:      {cfg.project.organism}")
        info(f"  Genome build:  {cfg.genome_build_resolved}")
        info(f"  Context:       {cfg.project.biological_context}")
        info(f"  Library:       {cfg.library.prep} / {cfg.library.layout}")
        comparisons = getattr(cfg, "comparisons", [])
        info(f"  Comparisons:   {len(comparisons)}")
        for c in comparisons:
            info(f"    {c.get('treatment','?')} vs {c.get('control','?')}")
    except ConfigError as e:
        error(f"Config validation failed:\n  {e}")
        sys.exit(1)


# ─────────────────────────────────────────────────────────────────────────────
# mosaic preflight
# ─────────────────────────────────────────────────────────────────────────────

@cli.command("preflight")
@click.argument("config_file", type=click.Path(exists=True))
@click.option("--auto-download", is_flag=True, default=False,
              help="Automatically download missing references")
def cmd_preflight(config_file: str, auto_download: bool) -> None:
    """Check config, verify all databases, optionally download missing ones."""
    setup_logging()
    try:
        cfg = load_config(config_file)
    except ConfigError as e:
        error(f"Config error: {e}")
        sys.exit(1)

    section("MOSAIC Pre-flight Check")
    info(f"Project:    {cfg.project.name}")
    info(f"Organism:   {cfg.project.organism} / {cfg.genome_build_resolved}")
    info(f"Context:    {cfg.project.biological_context}")
    info(f"Library:    {cfg.library.prep} / {cfg.library.layout}")
    info(f"Comparisons:")
    for c in getattr(cfg, "comparisons", []):
        if c.get("type") == "interaction":
            info(f"  [{c.get('name','interaction')}] interaction model")
        else:
            info(f"  {c.get('treatment')} vs {c.get('control')}")

    # Check databases
    from .db_manager import DatabaseManager
    db_mgr = DatabaseManager(cfg, cfg.project.organism)
    status = db_mgr.check_all()

    missing = [k for k, v in status.items() if v is False]

    if missing and auto_download:
        info(f"\nDownloading {len(missing)} missing databases...")
        for component in missing:
            db_mgr.download(component)
        # Re-check
        status = db_mgr.check_all()
        missing = [k for k, v in status.items() if v is False]

    if missing:
        warn(f"\n{len(missing)} databases missing: {', '.join(missing)}")
        warn("Run: mosaic preflight --auto-download config.yaml")
        sys.exit(1)
    else:
        success("\nAll databases present. Ready to run pipeline.")
        info("Next: mosaic run config.yaml")


# ─────────────────────────────────────────────────────────────────────────────
# mosaic run
# ─────────────────────────────────────────────────────────────────────────────

@cli.command("run")
@click.argument("config_file", type=click.Path(exists=True))
@click.option("--modules", "-m", default=None,
              help="Comma-separated list of modules to run (default: all enabled)")
@click.option("--from-module", default=None,
              help="Start from this module (skip earlier ones)")
@click.option("--force", "-f", is_flag=True, default=False,
              help="Force re-run of completed modules")
@click.option("--comparison", default=None,
              help="Run only this comparison (e.g. tumor_vs_normal)")
@click.option("--threads", "-t", default=None, type=int,
              help="Override thread count from config")
def cmd_run(config_file: str, modules: Optional[str],
            from_module: Optional[str], force: bool,
            comparison: Optional[str], threads: Optional[int]) -> None:
    """Run the MOSAIC pipeline."""
    setup_logging()
    try:
        cfg = load_config(config_file)
    except ConfigError as e:
        error(f"Config error: {e}")
        sys.exit(1)

    if threads:
        cfg.threads = threads

    from .pipeline import Pipeline
    scripts_dir = Path(__file__).parent.parent / "scripts"

    # Provenance: re-init logging with log_dir so this appears in session.log
    import datetime as _dt
    _log_dir = Path(cfg.project.wd) / "mosaic_output" / "logs"
    setup_logging(log_dir=_log_dir)
    info(f"[PROVENANCE] {_dt.datetime.now().isoformat()} "
         f"mosaic.__file__={Path(__file__).resolve()} "
         f"scripts_dir={scripts_dir.resolve()}")

    pipeline = Pipeline(
        cfg=cfg,
        scripts_dir=scripts_dir,
        force=force,
        from_module=from_module,
    )

    mods_list = modules.split(",") if modules else None

    # Filter comparisons if specified
    if comparison:
        parts = comparison.split("_vs_")
        if len(parts) == 2:
            cfg.comparisons = [{"treatment": parts[0], "control": parts[1]}]
        else:
            # Named comparison (e.g. interaction): find by name in config
            matches = [c for c in (cfg.comparisons or [])
                       if c.get("name") == comparison]
            if matches:
                cfg.comparisons = matches

    ok = pipeline.run(mods_list)
    sys.exit(0 if ok else 1)


# ─────────────────────────────────────────────────────────────────────────────
# mosaic status
# ─────────────────────────────────────────────────────────────────────────────

@cli.command("status")
@click.argument("config_file", type=click.Path(exists=True))
def cmd_status(config_file: str) -> None:
    """Show current pipeline execution status."""
    setup_logging()
    try:
        cfg = load_config(config_file)
    except ConfigError as e:
        error(f"Config error: {e}")
        sys.exit(1)

    from .pipeline import Pipeline
    scripts_dir = Path(__file__).parent.parent / "scripts"
    pipeline = Pipeline(cfg=cfg, scripts_dir=scripts_dir)
    pipeline.status()


# ─────────────────────────────────────────────────────────────────────────────
# mosaic compare  — run new comparison on existing alignments
# ─────────────────────────────────────────────────────────────────────────────

@cli.command("compare")
@click.argument("config_file", type=click.Path(exists=True))
@click.option("--treatment", "-t", required=True,
              help="Treatment condition name")
@click.option("--control", "-c", required=True,
              help="Control/reference condition name")
@click.option("--modules", "-m", default=None,
              help="Comma-separated modules (default: DE+GO+report)")
def cmd_compare(config_file: str, treatment: str,
                control: str, modules: Optional[str]) -> None:
    """
    Run a new treatment vs control comparison on existing alignments.
    Skips STAR alignment and TE quantification — reuses existing BAMs.
    """
    setup_logging()
    try:
        cfg = load_config(config_file)
    except ConfigError as e:
        error(f"Config error: {e}")
        sys.exit(1)

    # Override comparisons
    cfg.comparisons = [{"treatment": treatment, "control": control}]

    # Default modules for a new comparison
    default_mods = [
        "m08_deseq2", "m10_go", "m11d_te_de",
        "m17_rmats", "m18_network", "m14_report"
    ]
    mods_list = modules.split(",") if modules else default_mods

    info(f"New comparison: {treatment} vs {control}")
    info(f"Modules: {', '.join(mods_list)}")

    from .pipeline import Pipeline
    scripts_dir = Path(__file__).parent.parent / "scripts"
    pipeline = Pipeline(cfg=cfg, scripts_dir=scripts_dir, force=False)
    ok = pipeline.run(mods_list)
    sys.exit(0 if ok else 1)


# ─────────────────────────────────────────────────────────────────────────────
# mosaic download  — FASTQ download
# ─────────────────────────────────────────────────────────────────────────────

@cli.command("download")
@click.argument("config_file", type=click.Path(exists=True))
@click.option("--source", "-s", default="ena",
              type=click.Choice(["ena", "sra", "url", "ssh", "local"]),
              help="Download source")
@click.option("--accession", "-a", default=None,
              help="Single accession to download (e.g. SRR12345678)")
@click.option("--all", "download_all", is_flag=True, default=False,
              help="Download all samples defined in config")
@click.option("--threads", "-t", default=8, type=int,
              help="Download threads")
def cmd_download(config_file: str, source: str, accession: Optional[str],
                 download_all: bool, threads: int) -> None:
    """Download FASTQ files from ENA, SRA, URL, SSH, or local."""
    setup_logging()
    try:
        cfg = load_config(config_file)
    except ConfigError as e:
        error(f"Config error: {e}")
        sys.exit(1)

    from .pipeline import download_fastq

    fastq_dir = Path(cfg.samples.fastq_dir)

    if accession:
        ok = download_fastq(source, accession,
                            fastq_dir / accession, threads)
        sys.exit(0 if ok else 1)

    elif download_all:
        if not hasattr(cfg.samples, "download"):
            error("No 'samples.download' list in config.")
            sys.exit(1)
        all_ok = True
        for sample in cfg.samples.download:
            sid = sample["sample_id"]
            src = sample.get("source", source)
            acc = sample.get("accession", sid)
            ok  = download_fastq(src, acc, fastq_dir / sid, threads)
            if not ok:
                warn(f"Failed: {sid}")
                all_ok = False
        sys.exit(0 if all_ok else 1)

    else:
        error("Specify --accession or --all")
        sys.exit(1)


# ─────────────────────────────────────────────────────────────────────────────
# mosaic db  — database management
# ─────────────────────────────────────────────────────────────────────────────

@cli.group("db")
def cmd_db() -> None:
    """Manage MOSAIC reference databases."""
    pass


@cmd_db.command("check")
@click.argument("config_file", type=click.Path(exists=True))
def db_check(config_file: str) -> None:
    """Check which databases are present and up-to-date."""
    setup_logging()
    cfg = load_config(config_file)
    from .db_manager import DatabaseManager
    db_mgr = DatabaseManager(cfg, cfg.project.organism)
    db_mgr.check_all()


@cmd_db.command("download")
@click.argument("config_file", type=click.Path(exists=True))
@click.option("--component", "-c", default="all",
              type=click.Choice([
                  "all", "genome", "star_index", "kraken2",
                  "metaphlan", "chocophlan", "uniref", "arriba", "dfam"
              ]),
              help="Which component to download")
@click.option("--version", "-v", default="latest",
              help="Version to download (default: latest)")
@click.option("--organism", "-o", default=None,
              help="Override organism from config")
def db_download(config_file: str, component: str,
                version: str, organism: Optional[str]) -> None:
    """Download missing or outdated reference databases."""
    setup_logging()
    cfg = load_config(config_file)

    from .db_manager import DatabaseManager
    org = organism or cfg.project.organism
    db_mgr = DatabaseManager(cfg, org)
    db_mgr.download(component, version)


@cmd_db.command("list")
@click.argument("config_file", type=click.Path(exists=True))
def db_list(config_file: str) -> None:
    """List all available database versions."""
    setup_logging()
    cfg = load_config(config_file)
    from .db_manager import DatabaseManager
    db_mgr = DatabaseManager(cfg, cfg.project.organism)
    db_mgr.list_versions()


# ─────────────────────────────────────────────────────────────────────────────
# Keep backward-compatible entrypoints
# ─────────────────────────────────────────────────────────────────────────────

# mosaic-refs group (legacy compatibility)
@click.group(name="mosaic-refs")
@click.version_option(__version__, prog_name="mosaic-refs")
@click.option("--db", default=None, type=click.Path())
@click.pass_context
def refs_cli(ctx: click.Context, db: str | None) -> None:
    """MOSAIC reference database manager (legacy interface)."""
    setup_logging()
    ctx.ensure_object(dict)
    ctx.obj["db"] = Path(db) if db else Path("/shared/reference_db")


# mosaic-run group (legacy compatibility)
@click.group(name="mosaic-run")
@click.version_option(__version__, prog_name="mosaic-run")
def run_cli() -> None:
    """MOSAIC pipeline launcher (legacy interface)."""
    pass


@run_cli.command("init")
@click.argument("project_dir", type=click.Path())
@click.option("--organism", "-o", default="human")
@click.option("--context", "-c", default="tumor",
              type=click.Choice(["tumor", "aging", "congenital", "other"]))
def run_init(project_dir: str, organism: str, context: str) -> None:
    """Initialise a new MOSAIC project (legacy)."""
    setup_logging()
    p = Path(project_dir)
    p.mkdir(parents=True, exist_ok=True)
    cfg_path = p / "config.yaml"
    if not cfg_path.exists():
        cfg_path.write_text(generate_default_config(organism, context))
        success(f"Config written: {cfg_path}")
    info(f"Edit {cfg_path} then run: mosaic run {cfg_path}")


@run_cli.command("check")
@click.argument("config_file", type=click.Path(exists=True))
def run_check(config_file: str) -> None:
    """Validate config (legacy)."""
    setup_logging()
    try:
        cfg = load_config(config_file)
        success("Config valid.")
    except ConfigError as e:
        error(str(e))
        sys.exit(1)
