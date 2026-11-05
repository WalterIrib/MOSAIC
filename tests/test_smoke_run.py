#!/usr/bin/env python3
"""
Startup test for `mosaic run` (cmd_run).

Runs `mosaic run <config>` in a clean temporary directory and passes as soon as
the pipeline enters M02 (or the run completes). It fails if the process exits
before reaching M02 — e.g. an exception raised by cmd_run at startup, as in the
retired v1.1.5 (AttributeError on cfg.wd) — or if M02 is not reached within
the timeout.

Config: test_data/config_test.yaml when present (built by
test_data/generate_test_data.sh); otherwise a minimal config on the shipped
test_data/synth_ref with tiny FASTQs written to the temporary directory, so the
test also runs on a bare CI runner. project.wd always points to the temporary
directory: nothing is written into the repository.

Run:
    pytest tests/test_smoke_run.py -v
"""
from __future__ import annotations
import gzip
import os
import signal
import subprocess
import sys
import tempfile
import threading
from pathlib import Path

import yaml

REPO = Path(__file__).resolve().parent.parent
TEST_CONFIG = REPO / "test_data" / "config_test.yaml"
SYNTH = REPO / "test_data" / "synth_ref"
TIMEOUT_S = 300
M02_MARKER = "Running m02_trim"

# Command prefix for the CLI; equivalent to the `mosaic` entry point, but always
# runs the code of this checkout (PYTHONPATH below).
LAUNCHER = [sys.executable, "-m", "mosaic"]


def _write_fastq(path: Path, n: int = 4) -> None:
    with gzip.open(path, "wt") as fh:
        for i in range(n):
            fh.write(f"@r{i}\n{'ACGT' * 25}\n+\n{'I' * 100}\n")


def _minimal_config(tmp: Path) -> dict:
    fq = tmp / "fastq"
    fq.mkdir()
    rows = ["sample_id\tcondition\tR1\tR2"]
    for sid, cond in (("ctl1", "ctl"), ("ctl2", "ctl"), ("trt1", "trt"), ("trt2", "trt")):
        r1, r2 = fq / f"{sid}_R1.fastq.gz", fq / f"{sid}_R2.fastq.gz"
        _write_fastq(r1)
        _write_fastq(r2)
        rows.append(f"{sid}\t{cond}\t{r1}\t{r2}")
    (tmp / "coldata.tsv").write_text("\n".join(rows) + "\n")
    return {
        "project": {"name": "startup_test", "organism": "human", "genome_build": "hg38"},
        "sequencing": {"type": "PE", "read_length": 100, "strandedness": "auto"},
        "library": {"prep": "rRNA_depletion", "layout": "PE"},
        "samples": {"fastq_dir": str(fq), "coldata": str(tmp / "coldata.tsv")},
        "comparisons": [{"treatment": "trt", "control": "ctl"}],
        "databases": {
            "db_root": str(SYNTH),
            "genome_fa": str(SYNTH / "genome.fa"),
            "gtf": str(SYNTH / "annotation.gtf"),
            "star_index": str(tmp / "star_index"),
            "te_consensus": str(SYNTH / "dfam" / "consensus_TE.fa"),
            "te_satellite": str(SYNTH / "dfam" / "consensus_satellite.fa"),
            "rmsk_bed": str(SYNTH / "rmsk.txt.gz"),
            "te_gtf": str(SYNTH / "TE.gtf"),
            "auto_download": False,
        },
        "threads": 2,
    }


def _killpg(proc: subprocess.Popen, sig: int) -> None:
    try:
        os.killpg(proc.pid, sig)
    except ProcessLookupError:
        pass


def run_until_m02(launcher: list[str], timeout: int = TIMEOUT_S) -> tuple[bool, str]:
    """Start `<launcher> run <config>` in a clean temp dir.

    Returns (reached_m02_or_completed, captured_output). The process group is
    killed once M02 starts, so no tool keeps running after the test.
    """
    with tempfile.TemporaryDirectory(prefix="mosaic_startup_") as td:
        tmp = Path(td)
        if TEST_CONFIG.exists():
            cfg = yaml.safe_load(TEST_CONFIG.read_text())
        else:
            cfg = _minimal_config(tmp)
        cfg.setdefault("project", {})["wd"] = str(tmp / "run")
        cfg.setdefault("modules", {})["m00_download"] = False
        cfg_path = tmp / "config.yaml"
        cfg_path.write_text(yaml.safe_dump(cfg, sort_keys=False))

        env = dict(os.environ, PYTHONPATH=str(REPO))
        proc = subprocess.Popen(
            launcher + ["run", str(cfg_path)],
            cwd=tmp, env=env, text=True,
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            start_new_session=True,
        )
        # Watchdog: a silent hang is killed at the timeout and counts as a failure
        watchdog = threading.Timer(timeout, lambda: _killpg(proc, signal.SIGKILL))
        watchdog.start()
        out: list[str] = []
        reached = False
        try:
            for line in proc.stdout:
                out.append(line)
                if M02_MARKER in line:
                    reached = True
                    break
        finally:
            watchdog.cancel()
            if proc.poll() is None:
                _killpg(proc, signal.SIGTERM)
                try:
                    proc.wait(timeout=30)
                except subprocess.TimeoutExpired:
                    _killpg(proc, signal.SIGKILL)
                    proc.wait()
        if not reached and proc.returncode == 0:
            reached = True      # completed without printing the M02 header
        return reached, "".join(out)


def test_mosaic_run_starts():
    reached, output = run_until_m02(LAUNCHER)
    assert reached, ("`mosaic run` did not reach M02:\n" + output[-3000:])


if __name__ == "__main__":
    test_mosaic_run_starts()
    print("PASS test_mosaic_run_starts")
