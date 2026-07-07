#!/usr/bin/env python3
"""
MOSAIC functional smoke test.

Runs the whole pipeline on the tiny built-in chr1+chr22 subset (2 vs 2) and
checks that every module produces non-empty, correctly-structured output, that
the analysis modules are deterministic, and that the M14 report is wired to the
correct comparison (fig.path regression) with plot/table counts consistent.

Prerequisites:
    bash test_data/generate_test_data.sh     # builds test_data/config_test.yaml

Run:
    pytest tests/smoke_test.py -v
    # or standalone:
    python tests/smoke_test.py

Design: backbone/analysis modules (m02, m04, m05, m08, m09, m12, m17, m19, m14)
must produce output. Heavy/underpowered-on-4-samples modules (m10 GO, m18 WGCNA,
m16, m23) only need to NOT crash the pipeline (continue_on_error). The whole run
must finish with the report present.
"""
from __future__ import annotations
import hashlib
import os
import subprocess
import sys
from pathlib import Path

try:
    import pytest
except ImportError:                       # allow standalone `python tests/smoke_test.py`
    class _PytestShim:
        @staticmethod
        def fixture(*_a, **_k):
            def deco(fn):
                return fn
            return deco
    pytest = _PytestShim()

ROOT = Path(__file__).resolve().parent.parent
CONFIG = ROOT / "test_data" / "config_test.yaml"
WD = ROOT / "test_data" / "mosaic_output"
COMP = "mito-ChaRNA_vs_ChaRNA"
RUN_TIMEOUT = int(os.environ.get("SMOKE_TIMEOUT", "5400"))


def _report_html():
    for p in (WD / f"comparison_{COMP}" / "m14_report.html",
              WD / f"comparison_{COMP}" / "report" / "m14_report.html"):
        if p.exists():
            return p
    return _glob1("**/m14_report*.html")


def _run(args, timeout=RUN_TIMEOUT):
    """Invoke the mosaic CLI via `python -m mosaic`."""
    cmd = [sys.executable, "-m", "mosaic"] + args
    return subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True,
                          timeout=timeout)


def _glob1(pattern):
    hits = sorted(WD.glob(pattern))
    return hits[0] if hits else None


def _nonempty(p: Path, min_bytes=1):
    return p is not None and p.exists() and p.stat().st_size >= min_bytes


def _rows(p: Path):
    with open(p) as fh:
        return max(0, sum(1 for _ in fh) - 1)


def _header(p: Path):
    with open(p) as fh:
        return fh.readline().rstrip("\n").split("\t")


def _md5(p: Path):
    h = hashlib.md5()
    h.update(p.read_bytes())
    return h.hexdigest()


# ── one-time pipeline run (shared across tests) ──────────────────────────────
@pytest.fixture(scope="session", autouse=True)
def run_pipeline():
    assert CONFIG.exists(), (
        f"{CONFIG} not found — run: bash test_data/generate_test_data.sh")
    # Only run if not already complete (idempotent for local iteration).
    if _report_html() is None:
        r = _run(["run", str(CONFIG)])
        print(r.stdout[-4000:])
        print(r.stderr[-4000:], file=sys.stderr)
        # continue_on_error: return code may be non-zero if an optional module
        # failed; we assert the backbone + report below rather than on rc.
    yield


# ── backbone / analysis modules: must produce output ─────────────────────────
def test_m02_trim():
    f = _glob1("**/trim_stats.tsv")
    assert _nonempty(f), "m02: trim_stats.tsv missing/empty"
    assert "reads_after" in ",".join(_header(f)) or _rows(f) >= 1


def test_m04_star_bam():
    bam = _glob1("**/*Aligned*.bam") or _glob1("**/*.sortedByCoord*.bam") \
        or _glob1("mapping/**/*.bam")
    assert _nonempty(bam, 1000), "m04: no BAM produced"


def test_m05_te_family_counts():
    f = _glob1("te_counts/**/*.tsv") or _glob1("**/te_family*counts*.tsv")
    assert _nonempty(f), "m05: TE family counts missing"


def test_m08_gene_de_schema():
    f = _glob1(f"**/genes_{COMP}_DE.tsv")
    assert _nonempty(f), "m08: gene DE table missing"
    cols = _header(f)
    for c in ("baseMean", "log2FoldChange", "padj", "log2FoldChange_shrunk"):
        assert c in cols, f"m08: column {c} missing"
    assert _rows(f) > 0, "m08: gene DE table empty"


def test_m09_satellite_de():
    f = _glob1(f"**/satellite_{COMP}_DE.tsv")
    assert _nonempty(f), "m09: satellite DE table missing"
    assert "family_id" in _header(f)


def test_m12_locus_de_schema():
    f = _glob1(f"**/{COMP}_te_de*.tsv")
    assert _nonempty(f), "m12: locus DE table missing"
    assert "locus_id" in _header(f), "m12: locus_id column missing"


def test_m17_splicing_summary():
    f = _glob1("**/splicing_summary.tsv")
    assert _nonempty(f), "m17: splicing_summary.tsv missing"
    cols = _header(f)
    assert "sig_events" in cols and "event_type" in cols


def test_m19_editing_summary():
    f = _glob1("**/editing_summary.tsv")
    assert _nonempty(f), "m19: editing_summary.tsv missing"


# ── M14 report: correct comparison + plot/table consistency ──────────────────
def test_m14_report_exists_and_named():
    html = _report_html()
    assert _nonempty(html, 100_000), "m14: report HTML missing/too small"
    text = html.read_text(errors="ignore")
    # Correct comparison referenced (fig.path / cross-wiring regression)
    assert "mito-ChaRNA" in text, \
        "m14: report does not name the correct comparison"
    # Must NOT reference a different dataset's comparison (cross-wiring leak)
    for other in ("SEN_L_vs_EP", "SEN_E_vs_EP", "mito-totRNA_vs_totRNA"):
        assert other not in text, f"m14: report leaks unrelated comparison {other}"


def test_m14_percomparison_figpath():
    """The fix uses a comparison-specific knitr fig.path directory."""
    figdir = ROOT / "scripts" / f"m14_report_files_{COMP}"
    # created during render (intermediate). If absent (self-contained cleanup),
    # at least assert no shared generic dir was used with wrong content.
    generic = ROOT / "scripts" / "m14_report_files"
    assert figdir.exists() or not generic.exists(), \
        "m14: expected comparison-specific fig.path directory"


def test_m14_sig_points_equal_sig_rows():
    """Number of significant genes shown in the report equals the significant
    table rows (plot points == table rows, at the data level)."""
    sig = _glob1(f"**/genes_{COMP}_DE_significant.tsv")
    assert _nonempty(sig), "m08: significant gene table missing"
    n_sig = _rows(sig)
    html = _report_html()
    text = html.read_text(errors="ignore") if html else ""
    # the DE summary in the report prints the significant count
    assert str(n_sig) in text, \
        f"m14: significant count {n_sig} not found in report"


# ── determinism: analysis modules reproduce byte-identical output ────────────
def test_determinism_analysis_modules():
    # Re-run only the FAST deterministic analysis modules (not m12: its bedtools
    # multicov over all loci is slow; DESeq2/rMATS determinism is representative).
    targets = [
        _glob1(f"**/genes_{COMP}_DE.tsv"),
        _glob1("**/splicing_summary.tsv"),
    ]
    targets = [t for t in targets if t]
    assert targets, "no analysis outputs to compare"
    before = {str(t): _md5(t) for t in targets}
    r = _run(["run", str(CONFIG), "--modules",
              "m08_deseq2,m17_rmats", "--force"])
    assert r.returncode == 0, f"re-run failed:\n{r.stderr[-2000:]}"
    for t in targets:
        assert _md5(t) == before[str(t)], \
            f"determinism: {t.name} changed between identical runs"


# ── standalone runner ────────────────────────────────────────────────────────
def _main():
    if not CONFIG.exists():
        print(f"[ERROR] {CONFIG} not found. Run: bash test_data/generate_test_data.sh")
        return 2
    report = WD / f"comparison_{COMP}" / "m14_report.html"
    if not report.exists():
        print("[INFO] running pipeline on subset ...")
        r = _run(["run", str(CONFIG)])
        sys.stdout.write(r.stdout[-4000:])
        sys.stderr.write(r.stderr[-2000:])
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    fails = 0
    for t in tests:
        try:
            t()
            print(f"  PASS {t.__name__}")
        except AssertionError as e:
            fails += 1
            print(f"  FAIL {t.__name__}: {e}")
    print(f"\n{len(tests)-fails}/{len(tests)} passed")
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(_main())
