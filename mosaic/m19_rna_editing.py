#!/usr/bin/env python3
"""
mosaic/m19_rna_editing.py
M19 — RNA editing analysis using JACUSA2

Workflow:
1. Run JACUSA2 call-2 (normal vs tumor)
2. Parse output → filter A-to-I (A>G on + strand, T>C on - strand)
3. Annotate: gene, region (exon/intron/UTR), known ADAR sites (REDIportal)
4. Statistical filtering: score threshold, coverage, strand consistency
5. Output: editing sites TSV, summary, plots

Requirements: JACUSA2, samtools, bedtools, pysam, pandas
"""

import os
import sys
import subprocess
import argparse
import logging
from pathlib import Path

import pandas as pd
import numpy as np

logging.basicConfig(level=logging.INFO,
                    format="%(levelname)s\t%(asctime)s  %(message)s",
                    datefmt="%H:%M:%S")
log = logging.getLogger("M19")


# ── Argument parsing ───────────────────────────────────────────────────────────
def parse_args():
    p = argparse.ArgumentParser(
        description="M19 — RNA editing (JACUSA2 call-2)")
    p.add_argument("--bam-normal", required=True,
                   help="Comma-separated BAM files for normal condition")
    p.add_argument("--bam-tumor",  required=True,
                   help="Comma-separated BAM files for tumor condition")
    p.add_argument("--ref",        required=True,
                   help="Reference genome FASTA (indexed)")
    p.add_argument("--gtf",        required=True,
                   help="GENCODE GTF for annotation")
    p.add_argument("--outdir",     default="mosaic_output/m19_editing")
    p.add_argument("--threads",    type=int, default=20)
    p.add_argument("--min-cov",    type=int, default=5,
                   help="Min coverage per condition [5]")
    p.add_argument("--min-score",  type=float, default=1.5,
                   help="Min JACUSA2 score [1.5]")
    p.add_argument("--min-ef",     type=float, default=0.05,
                   help="Min editing frequency [0.05]")
    p.add_argument("--lib-type",   default="RF-FIRSTSTRAND",
                   choices=["RF-FIRSTSTRAND","FR-SECONDSTRAND","UNSTRANDED"])
    p.add_argument("--chromosomes", default=None,
                   help="Comma-separated chromosomes to analyse [all]")
    p.add_argument("--known-sites", default=None,
                   help="REDIportal BED/TSV of known A-to-I sites (optional)")
    p.add_argument("--jacusa2",    default="JACUSA2",
                   help="Path to JACUSA2 executable [JACUSA2]")
    return p.parse_args()


# ── JACUSA2 runner ─────────────────────────────────────────────────────────────
def run_jacusa2(args, outdir):
    raw_out = os.path.join(outdir, "jacusa2_raw.tsv")

    if os.path.exists(raw_out) and os.path.getsize(raw_out) > 1000:
        log.info("JACUSA2 output exists — skipping re-run")
        return raw_out

    # Expand ~ in BAM paths (subprocess doesn't expand shell tildes)
    def expand_bams(bam_str):
        return ",".join(os.path.expanduser(b.strip())
                        for b in bam_str.split(","))

    bam1 = expand_bams(args.bam_normal)
    bam2 = expand_bams(args.bam_tumor)

    cmd = [
        args.jacusa2, "call-2",
        "-r", raw_out,
        "-p", str(args.threads),
        "-R", os.path.expanduser(args.ref),
        "-P", args.lib_type,
        "-c", str(args.min_cov),
        "-q", "20",
        "-m", "20",
        "-f", "B",               # BED6-generic output (ref base in col 11)
        "-filterNH", "1",        # uniquely mapped only (NH=1)
        "-filterNM", "6",        # max 6 mismatches
        bam1,
        bam2
    ]
    # Note: DirMult:PHRED statistical model omitted — requires high coverage.
    # With full data (>30x per site) add: "-u", "DirMult:PHRED"

    if args.chromosomes:
        # JACUSA2 doesn't support multi-region directly — use BED
        chroms = args.chromosomes.split(",")
        region_bed = os.path.join(outdir, "regions.bed")
        # Get chromosome sizes from fai
        fai = args.ref + ".fai"
        if os.path.exists(fai):
            sizes = {}
            with open(fai) as f:
                for line in f:
                    parts = line.strip().split("\t")
                    sizes[parts[0]] = int(parts[1])
            with open(region_bed, "w") as f:
                for c in chroms:
                    if c in sizes:
                        f.write(f"{c}\t0\t{sizes[c]}\n")
            cmd += ["-b", region_bed]

    log.info("Running JACUSA2 call-2...")
    log.info("Command: " + " ".join(cmd))

    ret = subprocess.run(cmd, capture_output=False)
    if ret.returncode != 0:
        log.error(f"JACUSA2 failed (exit {ret.returncode})")
        sys.exit(1)

    log.info(f"JACUSA2 done: {raw_out}")
    return raw_out


# ── Parse JACUSA2 output ───────────────────────────────────────────────────────
def parse_jacusa2(raw_out, n_normal, n_tumor):
    """Parse JACUSA2 BED6-generic output into DataFrame.

    JACUSA2 call-2 with N normal + M tumor replicates outputs per-replicate
    base counts in separate columns:
      parts[6 .. 6+N-1]       = ACGT counts for each normal replicate
      parts[6+N .. 6+N+M-1]   = ACGT counts for each tumor replicate
      parts[6+N+M]             = info field
      parts[6+N+M+1]           = filter field
      parts[6+N+M+2]           = genomic ref base

    Aggregates (sums) base counts across all replicates for each condition
    before returning, so downstream EF computation uses pooled condition counts.
    """

    def parse_bases(b):
        v = list(map(int, b.split(",")))
        if len(v) == 4:
            return {"A": v[0], "C": v[1], "G": v[2], "T": v[3]}
        return {"A": 0, "C": 0, "G": 0, "T": 0}

    def aggregate_bases(parts, col_start, n_reps):
        agg = {"A": 0, "C": 0, "G": 0, "T": 0}
        for i in range(n_reps):
            b = parse_bases(parts[col_start + i])
            for base in agg:
                agg[base] += b[base]
        return agg

    n_rep_cols = n_normal + n_tumor
    col_info   = 6 + n_rep_cols          # info column
    col_filt   = 6 + n_rep_cols + 1      # filter column
    col_ref    = 6 + n_rep_cols + 2      # genomic reference base
    min_cols   = col_ref + 1

    rows = []
    with open(raw_out) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < min_cols:
                continue

            try:
                chrom  = parts[0]
                start  = int(parts[1])
                end    = int(parts[2])
                score  = float(parts[4])
                strand = parts[5]
                info   = parts[col_info]
                filt   = parts[col_filt]
                ref    = parts[col_ref]

                b1 = aggregate_bases(parts, 6,            n_normal)
                b2 = aggregate_bases(parts, 6 + n_normal, n_tumor)

                rows.append({
                    "chrom":  chrom,
                    "start":  start,
                    "end":    end,
                    "score":  score,
                    "strand": strand,
                    "ref":    ref,
                    "A1": b1["A"], "C1": b1["C"], "G1": b1["G"], "T1": b1["T"],
                    "A2": b2["A"], "C2": b2["C"], "G2": b2["G"], "T2": b2["T"],
                    "info":   info,
                    "filter": filt,
                })
            except (ValueError, IndexError):
                continue

    df = pd.DataFrame(rows)
    log.info(f"Parsed {len(df)} JACUSA2 sites "
             f"({n_normal} normal + {n_tumor} tumor replicates aggregated)")
    return df


# ── Filter for A-to-I editing sites ───────────────────────────────────────────
def filter_ato_i(df, min_score, min_ef):
    """
    A-to-I editing appears as:
    - A>G on + strand (ref=A, alt=G in RNA)
    - T>C on - strand (ref=T on genomic +, but RNA reads show C due to editing)

    We focus on sites where:
    1. ref = A (+ strand) or ref = T (- strand, complement)
    2. Edited base (G on +, C on -) present in tumor (condition 2)
    3. Editing frequency (EF) = edited_reads / total_reads >= min_ef
    4. JACUSA2 score >= min_score
    """
    if df.empty:
        return df

    # Compute total coverage and editing frequency per condition
    df["total1"] = df["A1"] + df["C1"] + df["G1"] + df["T1"]
    df["total2"] = df["A2"] + df["C2"] + df["G2"] + df["T2"]

    # A-to-I: ref=A, edited=G
    mask_plus = (df["ref"] == "A") & (df["strand"] == "+")
    df.loc[mask_plus, "edit_reads1"] = df.loc[mask_plus, "G1"]
    df.loc[mask_plus, "edit_reads2"] = df.loc[mask_plus, "G2"]
    df.loc[mask_plus, "edit_type"]   = "A>G"

    # A-to-I on - strand: ref=T (genomic), editing appears as T>C in reads
    mask_minus = (df["ref"] == "T") & (df["strand"] == "-")
    df.loc[mask_minus, "edit_reads1"] = df.loc[mask_minus, "C1"]
    df.loc[mask_minus, "edit_reads2"] = df.loc[mask_minus, "C2"]
    df.loc[mask_minus, "edit_type"]   = "T>C(A>I)"

    # Keep only A-to-I candidates
    df_ai = df[df["edit_type"].notna()].copy()

    # Editing frequency
    df_ai["EF1"] = df_ai["edit_reads1"] / df_ai["total1"].clip(lower=1)
    df_ai["EF2"] = df_ai["edit_reads2"] / df_ai["total2"].clip(lower=1)
    df_ai["dEF"] = df_ai["EF2"] - df_ai["EF1"]  # tumor - normal

    # Apply filters
    df_filt = df_ai[
        (df_ai["score"]  >= min_score) &
        (df_ai["EF2"]    >= min_ef)    &  # editing present in tumor
        (df_ai["total1"] >= 5)         &
        (df_ai["total2"] >= 5)
    ].copy()

    log.info(f"A-to-I candidates: {len(df_ai)} | After filters: {len(df_filt)}")
    return df_filt


# ── Filter all mismatch types ─────────────────────────────────────────────────
def filter_all_edits(df, min_score, min_ef, min_cov=3):
    """
    Extract all base-substitution mismatch types from JACUSA2 output.
    For each site, identifies the predominant non-ref base in condition 2 (treatment).
    Returns one row per (site, mutation_type) with editing frequency per condition.
    Filters: JACUSA2 score >= min_score, EF2 >= min_ef, coverage >= min_cov.
    """
    if df.empty:
        return pd.DataFrame()

    bases  = ["A", "C", "G", "T"]
    rows   = []

    df["total1"] = df["A1"] + df["C1"] + df["G1"] + df["T1"]
    df["total2"] = df["A2"] + df["C2"] + df["G2"] + df["T2"]

    for _, s in df.iterrows():
        ref = s["ref"]
        if ref not in bases:
            continue
        if s["total1"] < min_cov or s["total2"] < min_cov:
            continue
        if s["score"] < min_score:
            continue
        for alt in bases:
            if alt == ref:
                continue
            col1 = f"{alt}1"
            col2 = f"{alt}2"
            if col1 not in s.index or col2 not in s.index:
                continue
            ef2 = s[col2] / max(s["total2"], 1)
            if ef2 < min_ef:
                continue
            ef1 = s[col1] / max(s["total1"], 1)
            rows.append({
                "chrom":          s["chrom"],
                "start":          s["start"],
                "end":            s["end"],
                "strand":         s.get("strand", "."),
                "ref":            ref,
                "alt":            alt,
                "mutation_type":  f"{ref}>{alt}",
                "score":          s["score"],
                "EF1":            round(ef1, 4),
                "EF2":            round(ef2, 4),
                "dEF":            round(ef2 - ef1, 4),
                "cov1":           int(s["total1"]),
                "cov2":           int(s["total2"]),
            })

    result = pd.DataFrame(rows)
    log.info(f"All-type editing sites: {len(result)} (across all mismatch types)")
    return result


# ── Annotate with gene/region ──────────────────────────────────────────────────
def annotate_sites(df, gtf, outdir):
    """Annotate editing sites with gene name and genomic region using bedtools."""
    if df.empty:
        return df

    # Write sites as BED
    sites_bed = os.path.join(outdir, "editing_sites.bed")
    df[["chrom","start","end","edit_type","score","strand"]].to_csv(
        sites_bed, sep="\t", header=False, index=False)

    # Build gene BED from GTF
    gene_bed = os.path.join(outdir, "genes.bed")
    if not os.path.exists(gene_bed):
        log.info("Building gene BED from GTF...")
        records = []
        with open(gtf) as f:
            for line in f:
                if line.startswith("#"): continue
                parts = line.strip().split("\t")
                if len(parts) < 9 or parts[2] != "gene": continue
                import re
                gid  = re.search(r'gene_id "([^"]+)"', parts[8])
                gname= re.search(r'gene_name "([^"]+)"', parts[8])
                if gid:
                    records.append([
                        parts[0],
                        str(int(parts[3])-1),
                        parts[4],
                        gid.group(1),
                        gname.group(1) if gname else gid.group(1),
                        parts[6]
                    ])
        with open(gene_bed, "w") as f:
            for r in records:
                f.write("\t".join(r) + "\n")

    # bedtools intersect
    annot_bed = os.path.join(outdir, "editing_annotated.bed")
    cmd = (f"bedtools intersect -a {sites_bed} -b {gene_bed} "
           f"-wa -wb -loj > {annot_bed}")
    subprocess.run(cmd, shell=True, check=True)

    # Parse annotation
    try:
        ann = pd.read_csv(annot_bed, sep="\t", header=None,
                          names=["chrom","start","end","edit_type",
                                 "score","strand",
                                 "g_chrom","g_start","g_end",
                                 "gene_id","gene_name","gene_strand"])
        ann = ann.drop_duplicates(["chrom","start","end","edit_type"])

        df = df.merge(
            ann[["chrom","start","end","gene_id","gene_name"]],
            on=["chrom","start","end"], how="left"
        )
        log.info(f"Annotated: {df['gene_name'].notna().sum()} / {len(df)} sites in genes")
    except Exception as e:
        log.warning(f"Annotation failed: {e}")

    return df


# ── Known sites annotation (REDIportal) ───────────────────────────────────────
def annotate_known(df, known_sites_path, outdir):
    """Flag known A-to-I sites from REDIportal or similar database."""
    if known_sites_path is None or not os.path.exists(known_sites_path):
        df["known_site"] = False
        return df

    log.info("Annotating known editing sites...")
    sites_bed = os.path.join(outdir, "editing_sites.bed")
    known_bed = os.path.join(outdir, "known_intersect.bed")

    cmd = (f"bedtools intersect -a {sites_bed} -b {known_sites_path} "
           f"-wa -u > {known_bed}")
    subprocess.run(cmd, shell=True)

    if os.path.exists(known_bed) and os.path.getsize(known_bed) > 0:
        known_df = pd.read_csv(known_bed, sep="\t", header=None,
                               usecols=[0,1,2],
                               names=["chrom","start","end"])
        known_set = set(zip(known_df.chrom, known_df.start, known_df.end))
        df["known_site"] = df.apply(
            lambda r: (r.chrom, r.start, r.end) in known_set, axis=1)
        log.info(f"Known sites: {df['known_site'].sum()} / {len(df)}")
    else:
        df["known_site"] = False

    return df


# ── Summary statistics ────────────────────────────────────────────────────────
def compute_summary(df):
    """Compute per-gene and global editing statistics."""
    summary = {
        "total_sites":       len(df),
        "sites_in_genes":    int(df["gene_name"].notna().sum()) if "gene_name" in df else 0,
        "mean_EF_normal":    round(float(df["EF1"].mean()), 4),
        "mean_EF_tumor":     round(float(df["EF2"].mean()), 4),
        "mean_dEF":          round(float(df["dEF"].mean()), 4),
        "sites_increased":   int((df["dEF"] > 0.05).sum()),
        "sites_decreased":   int((df["dEF"] < -0.05).sum()),
    }
    if "known_site" in df:
        summary["known_sites"] = int(df["known_site"].sum())

    return summary


# ── Main ───────────────────────────────────────────────────────────────────────
def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    log.info("M19 — RNA editing analysis")
    log.info(f"Normal BAMs: {args.bam_normal}")
    log.info(f"Tumor BAMs:  {args.bam_tumor}")

    # 1. Run JACUSA2
    raw_out = run_jacusa2(args, args.outdir)

    # 2. Parse — count replicates from BAM lists so parse_jacusa2 knows column layout
    n_normal = len([b for b in args.bam_normal.split(",") if b.strip()])
    n_tumor  = len([b for b in args.bam_tumor.split(",")  if b.strip()])
    log.info(f"Replicates: {n_normal} normal, {n_tumor} tumor")
    df_raw = parse_jacusa2(raw_out, n_normal, n_tumor)
    if df_raw.empty:
        log.warning("No sites found — check coverage and parameters")
        sys.exit(0)

    # 3. Filter A-to-I
    df_ai = filter_ato_i(df_raw, args.min_score, args.min_ef)

    # 3b. Extract all mismatch types
    df_all_edits = filter_all_edits(df_raw, args.min_score, args.min_ef)
    if not df_all_edits.empty:
        out_all = os.path.join(args.outdir, "editing_sites_all.tsv")
        df_all_edits.to_csv(out_all, sep="\t", index=False)
        log.info(f"Saved all-type editing sites: {out_all} ({len(df_all_edits)} rows)")

    if df_ai.empty:
        log.warning("No A-to-I sites after filtering")
        pd.DataFrame().to_csv(
            os.path.join(args.outdir, "editing_sites_ato_i.tsv"),
            sep="\t", index=False)
        # Write zero-valued summary so downstream checks and report always find the file
        pd.DataFrame([{
            "metric": k, "value": v
        } for k, v in {
            "total_sites": 0, "sites_in_genes": 0,
            "mean_EF_normal": 0, "mean_EF_tumor": 0, "mean_dEF": 0,
            "sites_increased": 0, "sites_decreased": 0,
        }.items()]).to_csv(
            os.path.join(args.outdir, "editing_summary.tsv"),
            sep="\t", index=False)
        sys.exit(0)

    # 4. Annotate genes
    df_ai = annotate_sites(df_ai, args.gtf, args.outdir)

    # 5. Known sites
    df_ai = annotate_known(df_ai, args.known_sites, args.outdir)

    # 6. Save full table
    out_full = os.path.join(args.outdir, "editing_sites_ato_i.tsv")
    df_ai.to_csv(out_full, sep="\t", index=False)
    log.info(f"Saved: {out_full} ({len(df_ai)} sites)")

    # 7. Top sites (highest dEF)
    top_sites = df_ai.nlargest(200, "score")
    top_sites.to_csv(
        os.path.join(args.outdir, "editing_top200.tsv"),
        sep="\t", index=False)

    # 8. Per-gene summary
    if "gene_name" in df_ai.columns:
        gene_summary = df_ai.groupby("gene_name").agg(
            n_sites       =("chrom","count"),
            mean_EF_normal=("EF1","mean"),
            mean_EF_tumor =("EF2","mean"),
            mean_dEF      =("dEF","mean"),
            max_score     =("score","max")
        ).round(4).sort_values("n_sites", ascending=False).reset_index()
        gene_summary.to_csv(
            os.path.join(args.outdir, "editing_per_gene.tsv"),
            sep="\t", index=False)
        log.info(f"Genes with editing sites: {len(gene_summary)}")

    # 9. Summary
    summary = compute_summary(df_ai)
    summary_df = pd.DataFrame(
        list(summary.items()), columns=["metric","value"])
    summary_df.to_csv(
        os.path.join(args.outdir, "editing_summary.tsv"),
        sep="\t", index=False)

    log.info("\nM19 Summary:")
    for k, v in summary.items():
        log.info(f"  {k}: {v}")
    log.info(f"[OK] M19 complete. Output: {args.outdir}")


if __name__ == "__main__":
    main()
