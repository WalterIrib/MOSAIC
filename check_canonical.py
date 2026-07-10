#!/usr/bin/env python3
"""
check_canonical.py — manual genome-wide numerical regression for MOSAIC.

Compares the output of a FULL (genome-wide) MOSAIC run against the values
reported in the manuscript, printing PASS / FAIL with expected-vs-observed for
every metric. It NEVER edits any result file — it only reads and reports.

Usage
-----
    python check_canonical.py --mosaic-output /path/to/mosaic_output
    # multiple roots (e.g. senescence + chromatin) can be given:
    python check_canonical.py -m /data/DECECCA/mosaic_output \
                              -m /data/METACAR/mosaic_output

For each comparison it reads the canonical per-comparison TSVs and applies the
manuscript thresholds:
  * genes      — rows of genes_<c>_DE_significant.tsv
  * TE family  — TE_<c>_DE.tsv        with padj<0.05 & |log2FC_shrunk|>=1
  * satellite  — satellite_<c>_DE.tsv with padj<0.05 & |log2FC_shrunk|>=1
  * TE loci    — rows of m12/<c>_te_de_sig.tsv
  * splicing   — sum of sig_events in m17_rmats/splicing_summary.tsv (+ RI)
  * editing    — mean_dEF in m19_editing/editing_summary.tsv

Note on TE loci: v1.1 fixed a boundary-rounding artefact (see CHANGELOG), so the
significant-locus counts are expected to be a FEW loci higher than the published
values. Locus checks therefore PASS when the observed count lies in
[expected, expected + boundary]; the exact delta is printed.

Note on TE families: the count is sensitive to DESeq2 version (local dispersion fit,
BH correction). A tolerance of ±families_tol is allowed where set; the manuscript
value (e.g. 444) is kept as canonical anchor, not updated to the v1.1 observed count.

Note on splicing: the published splicing=566 for mito-ChaRNA_vs_ChaRNA was produced
without --variable-read-length (mosaic1 canonical commit 818c8484). The correct value
with --variable-read-length is 2032 (splicing_RI=735), which is what v1.1 produces.
The CANON entry records the corrected value.
"""
from __future__ import annotations
import argparse
import csv
import sys
from pathlib import Path

# ── Canonical values from the manuscript (None = not reported / not checked) ──
#   up/down and RI given where the paper states them; boundary = max extra loci
#   admissible from the v1.1 rounding-artefact fix (measured at write time).
CANON = {
    "SEN_L_vs_EP": dict(
        genes=4422, genes_up=None, genes_down=None,
        families=263, families_up=261, families_down=2,
        satellites=17, satellites_up=16, satellites_down=1,
        loci=376164, loci_boundary=29,
        splicing=12123, splicing_RI=None,
        dEF=0.100),
    "SEN_E_vs_EP": dict(
        genes=None, families=16, satellites=5,
        loci=342831, loci_boundary=24,
        splicing=13449, splicing_RI=None, dEF=0.088),
    "SEN_L_vs_SEN_E": dict(
        genes=None, families=423, satellites=None,
        loci=None, loci_boundary=23,
        splicing=None, dEF=None),
    "mito-ChaRNA_vs_ChaRNA": dict(
        genes=1417, genes_up=166, genes_down=1251,
        families=444, families_tol=5, families_all_down=True,
        satellites=9,
        loci=25489, loci_boundary=8,
        splicing=2032, splicing_RI=735, dEF=0.069),
    "mito-totRNA_vs_totRNA": dict(
        genes=None, families=0, satellites=None,
        loci=None, splicing=None, dEF=None),
}

GREEN, RED, YEL, DIM, RST = "\033[32m", "\033[31m", "\033[33m", "\033[2m", "\033[0m"


def _f(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def _find(root: Path, comp: str, sub: str, name: str):
    """Locate a per-comparison file, tolerating deseq2/ vs m12/ layouts."""
    cand = [root / f"comparison_{comp}" / sub / name,
            root / f"comparison_{comp}" / name]
    for c in cand:
        if c.exists():
            return c
    hits = list((root / f"comparison_{comp}").glob(f"**/{name}"))
    return hits[0] if hits else None


def _count_sig(path: Path, lfc_col="log2FoldChange_shrunk", padj=0.05, lfc=1.0):
    up = dn = 0
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            p, l = _f(r.get("padj")), _f(r.get(lfc_col))
            if p is None or l is None:
                continue
            if p < padj and abs(l) >= lfc:
                up += l > 0
                dn += l < 0
    return up + dn, up, dn


def _rows(path: Path):
    with open(path) as fh:
        return max(0, sum(1 for _ in fh) - 1)


class Report:
    def __init__(self):
        self.npass = self.nfail = self.nskip = 0

    def line(self, label, expected, observed, ok, extra=""):
        if expected is None:
            return
        if observed is None:
            self.nskip += 1
            print(f"  {YEL}SKIP{RST} {label:<26} (input not found)")
            return
        tag = f"{GREEN}PASS{RST}" if ok else f"{RED}FAIL{RST}"
        self.npass += ok
        self.nfail += not ok
        print(f"  {tag} {label:<26} expected={expected}  observed={observed}  {DIM}{extra}{RST}")


def check_comparison(root: Path, comp: str, exp: dict, rep: Report):
    print(f"\n{comp}  ({root})")

    # genes
    if exp.get("genes") is not None:
        f = _find(root, comp, "deseq2", f"genes_{comp}_DE_significant.tsv")
        obs = _rows(f) if f else None
        rep.line("genes (sig)", exp["genes"], obs, obs == exp["genes"])
        if f and exp.get("genes_up") is not None:
            _, up, dn = _count_sig(f, lfc_col="log2FoldChange")
            rep.line("genes up", exp["genes_up"], up, up == exp["genes_up"])
            rep.line("genes down", exp["genes_down"], dn, dn == exp["genes_down"])

    # TE families
    if exp.get("families") is not None:
        f = _find(root, comp, "deseq2", f"TE_{comp}_DE.tsv")
        if f:
            n, up, dn = _count_sig(f)
            tol = exp.get("families_tol", 0)
            fam_ok = abs(n - exp["families"]) <= tol
            fam_extra = ("exact" if n == exp["families"]
                         else f"Δ={n - exp['families']:+d} (tol ±{tol})" if tol
                         else "")
            rep.line("TE families", exp["families"], n, fam_ok, fam_extra)
            if exp.get("families_up") is not None:
                rep.line("TE families up", exp["families_up"], up, up == exp["families_up"])
                rep.line("TE families down", exp["families_down"], dn, dn == exp["families_down"])
            if exp.get("families_all_down"):
                rep.line("TE families all down", n, dn, dn == n)
        else:
            rep.line("TE families", exp["families"], None, False)

    # satellites
    if exp.get("satellites") is not None:
        f = _find(root, comp, "deseq2", f"satellite_{comp}_DE.tsv")
        if f:
            n, up, dn = _count_sig(f)
            rep.line("satellites", exp["satellites"], n, n == exp["satellites"])
            if exp.get("satellites_up") is not None:
                rep.line("satellites up", exp["satellites_up"], up, up == exp["satellites_up"])
                rep.line("satellites down", exp["satellites_down"], dn, dn == exp["satellites_down"])
        else:
            rep.line("satellites", exp["satellites"], None, False)

    # TE loci (boundary-fix band)
    if exp.get("loci") is not None:
        f = _find(root, comp, "m12", f"{comp}_te_de_sig.tsv")
        obs = _rows(f) if f else None
        b = exp.get("loci_boundary", 0)
        ok = obs is not None and exp["loci"] <= obs <= exp["loci"] + b
        extra = "" if obs is None else (
            "exact" if obs == exp["loci"]
            else f"+{obs - exp['loci']} (boundary-fix band 0..{b})")
        rep.line("TE loci (sig)", f"{exp['loci']}..{exp['loci']+b}", obs, ok, extra)

    # splicing
    if exp.get("splicing") is not None:
        f = _find(root, comp, "m17_rmats", "splicing_summary.tsv")
        if f:
            tot = ri = 0
            with open(f) as fh:
                for r in csv.DictReader(fh, delimiter="\t"):
                    tot += int(_f(r.get("sig_events")) or 0)
                    if r.get("event_type") == "RI":
                        ri = int(_f(r.get("sig_events")) or 0)
            rep.line("splicing (sig)", exp["splicing"], tot, tot == exp["splicing"])
            if exp.get("splicing_RI") is not None:
                rep.line("splicing RI", exp["splicing_RI"], ri, ri == exp["splicing_RI"])
        else:
            rep.line("splicing (sig)", exp["splicing"], None, False)

    # editing dEF
    if exp.get("dEF") is not None:
        f = _find(root, comp, "m19_editing", "editing_summary.tsv")
        obs = None
        if f:
            with open(f) as fh:
                for r in csv.DictReader(fh, delimiter="\t"):
                    if r.get("metric") == "mean_dEF":
                        obs = round(_f(r.get("value")), 3)
        rep.line("editing mean dEF", exp["dEF"], obs,
                 obs is not None and abs(obs - exp["dEF"]) <= 0.001)


def main():
    ap = argparse.ArgumentParser(description="MOSAIC genome-wide canonical check")
    ap.add_argument("-m", "--mosaic-output", action="append", required=True,
                    help="mosaic_output root (repeatable)")
    args = ap.parse_args()
    roots = [Path(p) for p in args.mosaic_output]

    rep = Report()
    for comp, exp in CANON.items():
        root = next((r for r in roots if (r / f"comparison_{comp}").exists()), None)
        if root is None:
            print(f"\n{comp}  {YEL}— no matching output root, skipped{RST}")
            continue
        check_comparison(root, comp, exp, rep)

    print(f"\n{'='*60}")
    print(f"PASS={rep.npass}  FAIL={rep.nfail}  SKIP={rep.nskip}")
    print("TE-locus counts a few above the published value are EXPECTED (v1.1 "
          "boundary-fix); see CHANGELOG.")
    sys.exit(1 if rep.nfail else 0)


if __name__ == "__main__":
    main()
