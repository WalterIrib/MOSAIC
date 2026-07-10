"""
tests/test_config.py
Unit tests for MOSAIC config loading and validation.
"""
import sys
import textwrap
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from mosaic.config import (
    load_config, generate_default_config, normalise_organism,
    ConfigError, MosaicConfig, ProjectConfig, LibraryConfig,
    _GENOME_BUILD_DEFAULT, _GO_DB,
)

# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────

def write_yaml(content: str) -> Path:
    """Write YAML content to a temp file, return its path."""
    f = tempfile.NamedTemporaryFile(mode="w", suffix=".yaml", delete=False)
    f.write(textwrap.dedent(content))
    f.close()
    return Path(f.name)


def minimal_config(**overrides) -> str:
    """Return a minimal valid YAML config string with optional overrides."""
    base = {
        "project.name": "test_project",
        "project.wd":   "/tmp/mosaic_test",
        "project.organism": "human",
        "project.genome_build": "hg38",
        "library.prep": "rRNA_depletion",
    }
    base.update(overrides)

    lines = [
        "project:",
        f"  name: {base['project.name']}",
        f"  wd: {base['project.wd']}",
        f"  organism: {base['project.organism']}",
        f"  genome_build: {base['project.genome_build']}",
        "library:",
        f"  prep: {base['library.prep']}",
        "  layout: PE",
    ]
    return "\n".join(lines) + "\n"


# ─────────────────────────────────────────────────────────────────────────────
# Test: organism normalisation
# ─────────────────────────────────────────────────────────────────────────────

def test_organism_normalisation():
    cases = [
        ("human",        "human"),
        ("Human",        "human"),
        ("Homo sapiens", "human"),
        ("homo sapiens", "human"),
        ("9606",         "human"),
        ("mouse",        "mouse"),
        ("Mus musculus", "mouse"),
        ("10090",        "mouse"),
        ("rat",          "rat"),
        ("10116",        "rat"),
        ("zebrafish",    "zebrafish"),
        ("7955",         "zebrafish"),
    ]
    for raw, expected in cases:
        result = normalise_organism(raw)
        assert result == expected, f"normalise_organism({raw!r}) = {result!r}, expected {expected!r}"
    print("  PASS test_organism_normalisation")


# ─────────────────────────────────────────────────────────────────────────────
# Test: valid minimal config loads
# ─────────────────────────────────────────────────────────────────────────────

def test_valid_minimal_config():
    p = write_yaml(minimal_config())
    cfg = load_config(p)
    assert cfg.project.organism == "human"
    assert cfg.genome_build_resolved == "hg38"
    assert cfg.library.layout == "PE"
    assert cfg.go_r_db == "org.Hs.eg.db"
    assert cfg.go_kegg == "hsa"
    print("  PASS test_valid_minimal_config")


def test_genome_build_auto_resolution():
    """genome_build: auto should resolve to default for organism."""
    p = write_yaml(minimal_config(**{"project.genome_build": "auto"}))
    cfg = load_config(p)
    assert cfg.genome_build_resolved == "hg38"

    p2 = write_yaml(minimal_config(**{
        "project.organism": "mouse",
        "project.genome_build": "auto",
    }))
    cfg2 = load_config(p2)
    assert cfg2.genome_build_resolved == "mm39"
    print("  PASS test_genome_build_auto_resolution")


# ─────────────────────────────────────────────────────────────────────────────
# Test: GO database auto-mapping
# ─────────────────────────────────────────────────────────────────────────────

def test_go_db_mapping():
    expected = {
        "human":     ("org.Hs.eg.db", "hsa"),
        "mouse":     ("org.Mm.eg.db", "mmu"),
        "rat":       ("org.Rn.eg.db", "rno"),
        "zebrafish": ("org.Dr.eg.db", "dre"),
    }
    for org, (rdb, kegg) in expected.items():
        build = _GENOME_BUILD_DEFAULT[org]
        p = write_yaml(minimal_config(**{
            "project.organism": org,
            "project.genome_build": build,
        }))
        cfg = load_config(p)
        assert cfg.go_r_db == rdb, f"{org}: expected {rdb}, got {cfg.go_r_db}"
        assert cfg.go_kegg == kegg
    print("  PASS test_go_db_mapping")


def test_unknown_organism_go_warning():
    """Unknown organism gets a warning placeholder, not an error."""
    p = write_yaml(minimal_config(**{
        "project.organism": "pig",
        "project.genome_build": "Sscrofa11.1",
    }))
    cfg = load_config(p)
    assert "UNKNOWN" in cfg.go_r_db
    print("  PASS test_unknown_organism_go_warning")


# ─────────────────────────────────────────────────────────────────────────────
# Test: CRITICAL safety checks
# ─────────────────────────────────────────────────────────────────────────────

def test_polya_blocked():
    """polyA selection must raise ConfigError."""
    yaml_str = minimal_config() + "  prep: polyA\n"
    p = write_yaml(yaml_str.replace("rRNA_depletion", "polyA"))
    try:
        load_config(p)
        assert False, "Expected ConfigError for polyA"
    except ConfigError as e:
        assert "polyA" in str(e) or "rRNA" in str(e)
    print("  PASS test_polya_blocked")


def test_low_complexity_filter_blocked():
    """low_complexity_filter: true must raise ConfigError."""
    extra = "\ntrimming:\n  low_complexity_filter: true\n"
    p = write_yaml(minimal_config() + extra)
    try:
        load_config(p)
        assert False, "Expected ConfigError for low_complexity_filter"
    except ConfigError as e:
        assert "low_complexity" in str(e).lower()
    print("  PASS test_low_complexity_filter_blocked")


def test_trim_poly_x_blocked():
    """trim_poly_x: true must raise ConfigError."""
    extra = "\ntrimming:\n  trim_poly_x: true\n"
    p = write_yaml(minimal_config() + extra)
    try:
        load_config(p)
        assert False, "Expected ConfigError for trim_poly_x"
    except ConfigError as e:
        assert "poly_x" in str(e).lower() or "satellite" in str(e).lower()
    print("  PASS test_trim_poly_x_blocked")


def test_apeglm_blocked_for_multipairwise():
    """apeglm + all_pairwise must raise ConfigError."""
    extra = textwrap.dedent("""
    de:
      mode: extended
      lfc_shrink_method: apeglm
      extended:
        all_pairwise: true
    """)
    p = write_yaml(minimal_config() + extra)
    try:
        load_config(p)
        assert False, "Expected ConfigError for apeglm + all_pairwise"
    except ConfigError as e:
        assert "apeglm" in str(e).lower()
    print("  PASS test_apeglm_blocked_for_multipairwise")


# ─────────────────────────────────────────────────────────────────────────────
# Test: defaults
# ─────────────────────────────────────────────────────────────────────────────

def test_defaults():
    """Verify all safety-critical defaults."""
    p = write_yaml(minimal_config())
    cfg = load_config(p)

    assert cfg.trimming.low_complexity_filter is False
    assert cfg.trimming.trim_poly_x is False
    assert cfg.trimming.cut_right is True
    assert cfg.trimming.min_length == 36
    assert cfg.mapping.bowtie2_score_min == "L,0,1.6"
    assert cfg.mapping.hisat2_allow_discordant is True
    assert cfg.quality.strandedness_concordance_threshold == 0.90
    assert cfg.de.lfc_shrink_method == "ashr"
    assert cfg.de.independent_filtering_te is False
    assert cfg.de.sizefactors_min_stable_genes == 500
    assert cfg.modules.microbial is False   # OPTIONAL by default
    assert cfg.modules.de_novo is False
    assert cfg.hisat2_dta is False          # auto resolved to False when de_novo=False
    print("  PASS test_defaults")


def test_hisat2_dta_auto_true_when_denovo():
    """hisat2_dta: auto with de_novo: true should resolve to True."""
    extra = "\nmodules:\n  de_novo: true\n"
    p = write_yaml(minimal_config() + extra)
    cfg = load_config(p)
    assert cfg.hisat2_dta is True
    print("  PASS test_hisat2_dta_auto_true_when_denovo")


def test_tumor_context_flag():
    import textwrap as _tw
    yaml_str = _tw.dedent("""
        project:
          name: test_project
          wd: /tmp/mosaic_test
          organism: human
          genome_build: hg38
          biological_context: tumor
        library:
          prep: rRNA_depletion
          layout: PE
    """)
    p = write_yaml(yaml_str)
    cfg = load_config(p)
    assert cfg.is_tumor is True
    print("  PASS test_tumor_context_flag")

def test_generate_default_config_roundtrip():
    """Generated default config must load and validate without errors."""
    for org in ("human", "mouse"):
        yaml_str = generate_default_config(organism=org, context="tumor")
        # Replace placeholder wd with real temp dir
        yaml_str = yaml_str.replace("/path/to/workdir", "/tmp/mosaic_test")
        p = write_yaml(yaml_str)
        cfg = load_config(p)
        assert cfg.project.organism == org
        assert "UNKNOWN" not in cfg.go_r_db
    print("  PASS test_generate_default_config_roundtrip")


# ─────────────────────────────────────────────────────────────────────────────
# Test: reference DB helpers
# ─────────────────────────────────────────────────────────────────────────────

def test_reference_presence_check():
    """is_reference_present returns False for empty/missing directory."""
    from mosaic.refs import is_reference_present
    import tempfile, os
    with tempfile.TemporaryDirectory() as td:
        db = Path(td)
        assert is_reference_present(db, "human", "hg38") is False
    print("  PASS test_reference_presence_check")


def test_manifest_roundtrip():
    """write_manifest + read_manifest roundtrip."""
    from mosaic.refs import write_manifest, read_manifest
    import tempfile
    with tempfile.TemporaryDirectory() as td:
        db = Path(td)
        data = {"organism": "human", "test_key": "test_value"}
        write_manifest(db, "human", "hg38", data)
        result = read_manifest(db, "human", "hg38")
        assert result["organism"] == "human"
        assert result["test_key"] == "test_value"
    print("  PASS test_manifest_roundtrip")


def test_update_log_append():
    """append_update_log creates and appends entries correctly."""
    from mosaic.refs import append_update_log, read_update_log
    import tempfile
    with tempfile.TemporaryDirectory() as td:
        db = Path(td)
        append_update_log(db, "human", "hg38", "genome", "--", "111",
                          "download", "test note")
        append_update_log(db, "human", "hg38", "dfam", "--", "3.8",
                          "download", "")
        rows = read_update_log(db, "human", "hg38")
        assert len(rows) == 2
        assert rows[0]["component"] == "genome"
        assert rows[1]["component"] == "dfam"
        assert rows[0]["new_version"] == "111"
    print("  PASS test_update_log_append")


def test_update_log_since_filter():
    """since filter works correctly."""
    from mosaic.refs import append_update_log, read_update_log
    import tempfile
    with tempfile.TemporaryDirectory() as td:
        db = Path(td)
        append_update_log(db, "human", "hg38", "genome", "--", "111",
                          "download", "old")
        rows_all  = read_update_log(db, "human", "hg38")
        rows_fut  = read_update_log(db, "human", "hg38", since="2099-01-01")
        assert len(rows_all) == 1
        assert len(rows_fut) == 0
    print("  PASS test_update_log_since_filter")


def test_nfs_lock():
    """RefLock acquires and releases correctly. Double-acquire fails."""
    from mosaic.refs import RefLock
    import tempfile
    with tempfile.TemporaryDirectory() as td:
        db = Path(td)
        with RefLock(db, "human", "hg38"):
            # lock is held
            lock_dir = db / "human" / "hg38" / ".lock"
            assert lock_dir.exists()
            # second acquire while locked should fail
            try:
                with RefLock(db, "human", "hg38"):
                    assert False, "Expected RuntimeError on double lock"
            except RuntimeError:
                pass
        # lock released
        assert not lock_dir.exists()
    print("  PASS test_nfs_lock")


def test_gene_type_map_extraction():
    """extract_gene_type_map parses a minimal GTF correctly."""
    from mosaic.refs import extract_gene_type_map
    import tempfile

    gtf_content = textwrap.dedent("""\
        #!genome-build GRCh38
        1\tHAVANA\tgene\t11869\t14409\t.\t+\t.\tgene_id "ENSG00000223972"; gene_name "DDX11L1"; gene_type "transcribed_unprocessed_pseudogene";
        1\tHAVANA\tgene\t14404\t29570\t.\t-\t.\tgene_id "ENSG00000227232"; gene_name "WASH7P"; gene_biotype "unprocessed_pseudogene";
        1\tENSEMBL\tgene\t65419\t71585\t.\t+\t.\tgene_id "ENSG00000186092"; gene_name "OR4F5"; gene_type "protein_coding";
    """)

    with tempfile.TemporaryDirectory() as td:
        gtf_path = Path(td) / "test.gtf"
        out_path  = Path(td) / "gene_type_map.tsv"
        gtf_path.write_text(gtf_content)
        extract_gene_type_map(gtf_path, out_path)

        import csv
        rows = list(csv.DictReader(open(out_path), delimiter="\t"))
        assert len(rows) == 3
        ids = {r["gene_id"] for r in rows}
        assert "ENSG00000186092" in ids
        # gene_biotype fallback
        wash = next(r for r in rows if r["gene_id"] == "ENSG00000227232")
        assert wash["gene_type"] == "unprocessed_pseudogene"

    print("  PASS test_gene_type_map_extraction")


def test_md5sum():
    """md5sum produces correct hash."""
    import tempfile, hashlib
    from mosaic.refs import md5sum
    with tempfile.NamedTemporaryFile(delete=False, mode="wb") as f:
        f.write(b"hello mosaic\n")
        fname = Path(f.name)
    expected = hashlib.md5(b"hello mosaic\n").hexdigest()
    assert md5sum(fname) == expected
    print("  PASS test_md5sum")


# ─────────────────────────────────────────────────────────────────────────────
# Test: CLI smoke tests (no network required)
# ─────────────────────────────────────────────────────────────────────────────

def test_cli_init_and_check():
    """mosaic-run init creates a config; mosaic-run check validates it."""
    import subprocess, tempfile
    with tempfile.TemporaryDirectory() as td:
        proj = Path(td) / "myproject"
        # init
        r = subprocess.run(
            ["mosaic-run", "init", str(proj), "--organism", "human", "--context", "tumor"],
            capture_output=True, text=True
        )
        assert r.returncode == 0, r.stderr
        cfg_path = proj / "config.yaml"
        assert cfg_path.exists()

        # Patch wd to a real path so validation passes
        content = cfg_path.read_text()
        content = content.replace("/path/to/workdir", str(proj))
        cfg_path.write_text(content)

        # check
        r2 = subprocess.run(
            ["mosaic-run", "check", str(cfg_path)],
            capture_output=True, text=True
        )
        assert r2.returncode == 0, r2.stderr
        assert "valid" in r2.stderr.lower() or "valid" in r2.stdout.lower()
    print("  PASS test_cli_init_and_check")


def test_cli_refs_list_empty():
    """mosaic-refs list on empty dir exits 0 and prints nothing."""
    import subprocess, tempfile
    with tempfile.TemporaryDirectory() as td:
        r = subprocess.run(
            ["mosaic-refs", "--db", td, "list"],
            capture_output=True, text=True
        )
        assert r.returncode == 0
    print("  PASS test_cli_refs_list_empty")


def test_cli_refs_check_missing():
    """mosaic-refs check on missing organism exits non-zero."""
    import subprocess, tempfile
    with tempfile.TemporaryDirectory() as td:
        r = subprocess.run(
            ["mosaic-refs", "--db", td, "check", "--organism", "human"],
            capture_output=True, text=True
        )
        assert r.returncode != 0
    print("  PASS test_cli_refs_check_missing")


# ─────────────────────────────────────────────────────────────────────────────
# Runner
# ─────────────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    tests = [
        test_organism_normalisation,
        test_valid_minimal_config,
        test_genome_build_auto_resolution,
        test_go_db_mapping,
        test_unknown_organism_go_warning,
        test_polya_blocked,
        test_low_complexity_filter_blocked,
        test_trim_poly_x_blocked,
        test_apeglm_blocked_for_multipairwise,
        test_defaults,
        test_hisat2_dta_auto_true_when_denovo,
        test_tumor_context_flag,
        test_generate_default_config_roundtrip,
        test_reference_presence_check,
        test_manifest_roundtrip,
        test_update_log_append,
        test_update_log_since_filter,
        test_nfs_lock,
        test_gene_type_map_extraction,
        test_md5sum,
        test_cli_init_and_check,
        test_cli_refs_list_empty,
        test_cli_refs_check_missing,
    ]

    passed = failed = 0
    print(f"\nRunning {len(tests)} tests for MOSAIC Fase 1...\n")
    for t in tests:
        try:
            t()
            passed += 1
        except Exception as e:
            print(f"  FAIL {t.__name__}: {e}")
            import traceback; traceback.print_exc()
            failed += 1

    print(f"\n{'─'*50}")
    print(f"Results: {passed} passed, {failed} failed")
    sys.exit(0 if failed == 0 else 1)
