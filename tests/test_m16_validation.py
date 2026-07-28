"""
tests/test_m16_validation.py
Unit tests for M16 fixes:
  1. novel_transcripts.tmap generation (Python, no external tools)
  2. validate_count_matrix in m16_stringtie_de.R rejects FPKM/TPM input
"""
import subprocess
import sys
import textwrap
import tempfile
from pathlib import Path

REPO = Path(__file__).parent.parent
R_SCRIPT = REPO / "scripts" / "m16_stringtie_de.R"

# ── helpers ───────────────────────────────────────────────────────────────────

NOVEL_CODES = {"u", "i", "j", "x", "s", "e", "p"}

def _filter_novel_tmap(tmap_lines: list[str]) -> list[str]:
    """Mirror of the Python logic added to pipeline.py _run_m16_stringtie."""
    result = []
    for line in tmap_lines:
        parts = line.split("\t")
        if len(parts) >= 3 and parts[2] in NOVEL_CODES:
            result.append(line)
    return result


def _run_r_validate(matrix_r_code: str) -> subprocess.CompletedProcess:
    """Run validate_count_matrix from m16_stringtie_de.R in a subprocess."""
    r_code = textwrap.dedent(f"""
        # Source only the validate_count_matrix function from the script
        lines <- readLines("{R_SCRIPT}")
        # Extract lines that define the function
        start <- grep("validate_count_matrix <- function", lines)
        end   <- NA
        if (length(start) > 0) {{
          depth <- 0
          for (i in seq(start[1], length(lines))) {{
            depth <- depth + (nchar(gsub("[^{{]","",lines[i])) -
                               nchar(gsub("[^}}]","",lines[i])))
            if (i > start[1] && depth <= 0) {{ end <- i; break }}
          }}
        }}
        if (!is.na(end)) {{
          eval(parse(text=paste(lines[start[1]:end], collapse="\\n")))
        }} else {{
          stop("Could not extract validate_count_matrix from script")
        }}
        {matrix_r_code}
    """)
    return subprocess.run(
        ["conda", "run", "-n", "mosaic", "Rscript", "-e", r_code],
        capture_output=True, text=True
    )


# ── tests: novel_tmap filter ─────────────────────────────────────────────────

def test_novel_tmap_keeps_novel_codes():
    tmap_lines = [
        "ENSG1\tENST1\t=\tMSTRG.1\tMSTRG.1.1\t3\t10.0\t8.0\t5.0\t500\tENST1\t500\n",
        "ENSG2\tENST2\tu\tMSTRG.2\tMSTRG.2.1\t2\t5.0\t4.0\t3.0\t300\tENST2\t0\n",
        "ENSG3\tENST3\ti\tMSTRG.3\tMSTRG.3.1\t1\t2.0\t1.5\t1.0\t200\t.\t0\n",
        "ENSG4\tENST4\tk\tMSTRG.4\tMSTRG.4.1\t4\t20.0\t18.0\t15.0\t800\tENST4\t800\n",
        "ENSG5\tENST5\tj\tMSTRG.5\tMSTRG.5.1\t3\t8.0\t7.0\t4.0\t400\t.\t0\n",
    ]
    result = _filter_novel_tmap(tmap_lines)
    assert len(result) == 3
    assert all(line.split("\t")[2] in NOVEL_CODES for line in result)


def test_novel_tmap_empty_when_all_annotated():
    tmap_lines = [
        "ENSG1\tENST1\t=\tMSTRG.1\tMSTRG.1.1\t3\t10.0\t8.0\t5.0\t500\tENST1\t500\n",
        "ENSG2\tENST2\tk\tMSTRG.2\tMSTRG.2.1\t4\t20.0\t18.0\t15.0\t800\tENST2\t800\n",
        "ENSG3\tENST3\tc\tMSTRG.3\tMSTRG.3.1\t2\t5.0\t4.0\t3.0\t300\tENST3\t300\n",
    ]
    result = _filter_novel_tmap(tmap_lines)
    assert result == []


def test_novel_tmap_all_novel_codes_included():
    expected_novel = ["u", "i", "j", "x", "s", "e", "p"]
    tmap_lines = [
        f"ENSG{n}\tENST{n}\t{code}\tMSTRG.{n}\tMSTRG.{n}.1\t1\t1.0\t1.0\t1.0\t100\t.\t0\n"
        for n, code in enumerate(expected_novel + ["=", "k", "c"])
    ]
    result = _filter_novel_tmap(tmap_lines)
    result_codes = {line.split("\t")[2] for line in result}
    assert result_codes == set(expected_novel)


# ── tests: R validate_count_matrix ───────────────────────────────────────────

def test_r_validate_accepts_integer_counts():
    proc = _run_r_validate("""
        mat <- matrix(c(10L, 20L, 0L, 5L, 100L, 200L, 0L, 50L), nrow=2)
        validate_count_matrix(mat, "test")
        cat("PASS\\n")
    """)
    assert proc.returncode == 0, proc.stderr
    assert "PASS" in proc.stdout


def test_r_validate_rejects_fpkm():
    """FPKM values (non-integer floats) must be rejected."""
    proc = _run_r_validate("""
        mat <- matrix(c(12.34, 5.67, 0.89, 100.1), nrow=2)
        tryCatch(
            validate_count_matrix(mat, "FPKM test"),
            error = function(e) cat("REJECTED:", conditionMessage(e), "\\n")
        )
    """)
    assert proc.returncode == 0, proc.stderr
    assert "REJECTED" in proc.stdout
    assert "non-integer" in proc.stdout.lower() or "fpkm" in proc.stdout.lower()


def test_r_validate_rejects_negative():
    proc = _run_r_validate("""
        mat <- matrix(c(10L, -1L, 5L, 20L), nrow=2)
        tryCatch(
            validate_count_matrix(mat, "neg test"),
            error = function(e) cat("REJECTED:", conditionMessage(e), "\\n")
        )
    """)
    assert proc.returncode == 0, proc.stderr
    assert "REJECTED" in proc.stdout
    assert "negative" in proc.stdout.lower()


def test_r_validate_accepts_zero_counts():
    """Zero is a valid count (unexpressed transcript)."""
    proc = _run_r_validate("""
        mat <- matrix(c(0L, 0L, 100L, 200L), nrow=2)
        validate_count_matrix(mat, "zero test")
        cat("PASS\\n")
    """)
    assert proc.returncode == 0, proc.stderr
    assert "PASS" in proc.stdout


# ── runner ────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    tests = [
        test_novel_tmap_keeps_novel_codes,
        test_novel_tmap_empty_when_all_annotated,
        test_novel_tmap_all_novel_codes_included,
        test_r_validate_accepts_integer_counts,
        test_r_validate_rejects_fpkm,
        test_r_validate_rejects_negative,
        test_r_validate_accepts_zero_counts,
    ]
    passed = failed = 0
    print(f"\nRunning {len(tests)} M16 validation tests\n")
    for t in tests:
        try:
            t()
            print(f"  PASS  {t.__name__}")
            passed += 1
        except Exception as e:
            print(f"  FAIL  {t.__name__}: {e}")
            failed += 1
    print(f"\n{passed} passed, {failed} failed")
    sys.exit(0 if failed == 0 else 1)
