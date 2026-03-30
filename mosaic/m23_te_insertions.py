#!/usr/bin/env python3
"""
mosaic/m23_te_insertions.py
M23 - Somatic TE insertion detection from STAR chimeric reads

Workflow:
1. Parse chimeric SAM files -> extract soft-clipped sequences
2. Align clipped sequences to TE consensus (pseudochromosomes) via minimap2/bwa
3. Cluster insertion sites by genomic position
4. Compare tumor vs normal -> candidate somatic insertions
5. Annotate: TE family, insertion site gene, TSD (target site duplication)

Requirements: pysam, pandas, bwa or minimap2
"""

import os
import sys
import re
import subprocess
import argparse
import logging
from pathlib import Path
from collections import defaultdict

import pysam
import pandas as pd
import numpy as np

logging.basicConfig(level=logging.INFO,
                    format="%(levelname)s\t%(asctime)s  %(message)s",
                    datefmt="%H:%M:%S")
log = logging.getLogger("M23")


def parse_args():
    p = argparse.ArgumentParser(
        description="M23 - TE somatic insertion detection")
    p.add_argument("--chimeric-normal", required=True,
                   help="Comma-separated Chimeric.out.sam files for normal")
    p.add_argument("--chimeric-tumor",  required=True,
                   help="Comma-separated Chimeric.out.sam files for tumor")
    p.add_argument("--te-consensus",    required=True,
                   help="TE consensus FASTA (pseudochromosomes)")
    p.add_argument("--genome-fa",       required=True,
                   help="Reference genome FASTA (indexed)")
    p.add_argument("--gtf",             required=True,
                   help="GENCODE GTF for annotation")
    p.add_argument("--outdir",          default="mosaic_output/m23_te_insertions")
    p.add_argument("--threads",         type=int, default=20)
    p.add_argument("--min-clip",        type=int, default=30,
                   help="Min soft-clip length to consider [30]")
    p.add_argument("--min-mapq",        type=int, default=10,
                   help="Min MAPQ for clipped read [10]")
    p.add_argument("--cluster-dist",    type=int, default=500,
                   help="Max distance to cluster insertion sites [500]")
    p.add_argument("--min-support",     type=int, default=2,
                   help="Min reads supporting insertion [2]")
    p.add_argument("--aligner",         default="bwa",
                   choices=["bwa","minimap2"],
                   help="Aligner for clipped sequences [bwa]")
    p.add_argument("--condition-control",   default="control",
                   help="Name of the control condition (replaces 'normal' labels)")
    p.add_argument("--condition-treatment", default="treatment",
                   help="Name of the treatment condition (replaces 'tumor' labels)")
    return p.parse_args()


# ── CIGAR parsing ──────────────────────────────────────────────────────────────
def parse_cigar(cigar):
    """Parse CIGAR string into list of (length, op) tuples."""
    return [(int(n), op) for n, op in re.findall(r'(\d+)([MIDNSHP=X])', cigar)]


def get_soft_clips(cigar_ops, seq, qual):
    """Extract soft-clipped sequences and their positions (5' and 3')."""
    clips = []
    if not cigar_ops:
        return clips

    # 5' clip
    if cigar_ops[0][1] == 'S' and cigar_ops[0][0] >= 1:
        clip_len = cigar_ops[0][0]
        clips.append({
            "seq":      seq[:clip_len],
            "qual":     qual[:clip_len] if qual else "",
            "position": "5prime",
            "length":   clip_len
        })

    # 3' clip
    if cigar_ops[-1][1] == 'S' and cigar_ops[-1][0] >= 1:
        clip_len = cigar_ops[-1][0]
        clips.append({
            "seq":      seq[-clip_len:],
            "qual":     qual[-clip_len:] if qual else "",
            "position": "3prime",
            "length":   clip_len
        })
    return clips


# ── Parse chimeric SAM ─────────────────────────────────────────────────────────
def parse_chimeric_sam(sam_files, condition, min_clip, min_mapq):
    """
    Parse STAR chimeric SAM files and extract soft-clipped sequences.
    Returns DataFrame with clipped sequences and their genomic positions.
    """
    all_clips = []
    total_reads = 0

    for sam_path in sam_files:
        sam_path = os.path.expanduser(sam_path.strip())
        if not os.path.exists(sam_path):
            log.warning(f"Not found: {sam_path}")
            continue

        log.info(f"Parsing: {os.path.basename(sam_path)}")
        sample_id = os.path.basename(sam_path).split("_Chimeric")[0]

        try:
            # Parse SAM directly (avoids pysam index requirement)
            with open(sam_path, "r") as fh:
                for line in fh:
                    if line.startswith("@"):
                        continue
                    parts = line.rstrip().split("	")
                    if len(parts) < 10:
                        continue

                    read_id = parts[0]
                    flag    = int(parts[1])
                    chrom   = parts[2]
                    pos     = int(parts[3]) - 1  # 0-based
                    mapq    = int(parts[4])
                    cigar   = parts[5]
                    seq     = parts[9]

                    # Skip unmapped
                    if flag & 4 or chrom == "*" or cigar == "*":
                        continue
                    if mapq < min_mapq:
                        continue

                    total_reads += 1
                    strand    = "-" if (flag & 16) else "+"
                    cigar_ops = parse_cigar(cigar)
                    clips     = get_soft_clips(cigar_ops, seq, "")

                    for clip in clips:
                        if clip["length"] < min_clip:
                            continue
                        all_clips.append({
                            "read_id":   read_id,
                            "sample_id": sample_id,
                            "condition": condition,
                            "chrom":     chrom,
                            "pos":       pos,
                            "strand":    strand,
                            "mapq":      mapq,
                            "clip_pos":  clip["position"],
                            "clip_seq":  clip["seq"],
                            "clip_len":  clip["length"]
                        })
        except Exception as e:
            log.warning(f"Error parsing {sam_path}: {e}")
            continue

    log.info(f"{condition}: {total_reads} chimeric reads, "
             f"{len(all_clips)} soft-clips >= {min_clip}bp")
    return pd.DataFrame(all_clips)


# ── Index and align clipped sequences ─────────────────────────────────────────
def index_te_consensus(te_fa, outdir, aligner="bwa"):
    """Index TE consensus FASTA for alignment."""
    te_fa = os.path.expanduser(te_fa)
    idx_done = os.path.join(outdir, ".te_index_done")
    if os.path.exists(idx_done):
        log.info("TE index exists — skipping")
        return te_fa

    log.info(f"Indexing TE consensus with {aligner}...")
    if aligner == "bwa":
        ret = subprocess.run(["bwa", "index", te_fa], capture_output=True)
        if ret.returncode != 0:
            log.error("BWA index failed")
            log.error(ret.stderr.decode())
            return None
    else:
        # minimap2 doesn't need pre-indexing
        pass

    Path(idx_done).touch()
    return te_fa


def align_clips_to_te(clips_df, te_fa, outdir, threads, aligner="bwa"):
    """
    Align soft-clipped sequences to TE consensus.
    Returns DataFrame with TE mapping info for each clip.
    """
    if clips_df.empty:
        return pd.DataFrame()

    # Write clipped sequences as FASTA
    fasta_path = os.path.join(outdir, "clipped_sequences.fa")
    with open(fasta_path, "w") as f:
        for i, row in clips_df.iterrows():
            f.write(f">{row['read_id']}_{row['condition']}_{row['sample_id']}_"
                    f"{row['chrom']}_{row['pos']}_{row['clip_pos']}\n"
                    f"{row['clip_seq']}\n")

    log.info(f"Aligning {len(clips_df)} clipped sequences to TE consensus...")

    sam_out = os.path.join(outdir, "clips_vs_te.sam")
    te_fa_exp = os.path.expanduser(te_fa)

    if aligner == "bwa":
        cmd = (f"bwa mem -t {threads} -k 15 -T 20 "
               f"{te_fa_exp} {fasta_path} > {sam_out} 2>/dev/null")
    else:
        cmd = (f"minimap2 -a -x sr -t {threads} "
               f"{te_fa_exp} {fasta_path} > {sam_out} 2>/dev/null")

    ret = subprocess.run(cmd, shell=True)
    if ret.returncode != 0:
        log.error(f"Alignment failed: {cmd}")
        return pd.DataFrame()

    # Parse alignment results
    te_hits = []
    try:
        sam = pysam.AlignmentFile(sam_out, "r", check_sq=False)
        for read in sam.fetch(until_eof=True):
            if read.is_unmapped:
                continue
            if read.mapping_quality < 5:
                continue

            # Parse read name to get original genomic info
            parts = read.query_name.split("_")
            # Format: readid_condition_sampleid_chrom_pos_clippos
            te_hits.append({
                "clip_id":    read.query_name,
                "te_chrom":   read.reference_name,   # TE family name
                "te_start":   read.reference_start,
                "te_mapq":    read.mapping_quality,
                "te_strand":  "-" if read.is_reverse else "+"
            })
        sam.close()
    except Exception as e:
        log.warning(f"Error parsing TE alignment: {e}")

    log.info(f"Clips mapping to TE consensus: {len(te_hits)}")
    return pd.DataFrame(te_hits)


# ── Extract TE family from pseudochromosome name ───────────────────────────────
def parse_te_family(te_chrom):
    """Extract TE family/class from pseudochromosome coordinate."""
    # Pseudochromosome format: family_name|class|... or just family_name
    if te_chrom is None:
        return "Unknown", "Unknown"
    parts = str(te_chrom).split("|")
    name  = parts[0] if parts else te_chrom
    cls   = parts[1] if len(parts) > 1 else "Unknown"

    # Classify
    n = name.upper()
    if any(x in n for x in ["ALU","MIR","SINE"]):
        macro = "SINE"
    elif any(x in n for x in ["L1","L2","LINE","CR1"]):
        macro = "LINE"
    elif any(x in n for x in ["LTR","HERV","ERV","MLT","THE"]):
        macro = "LTR/ERV"
    elif "SVA" in n:
        macro = "SVA"
    elif any(x in n for x in ["DNA","TIGGER","CHARLIE","HAT","MARINER"]):
        macro = "DNA transposon"
    elif any(x in n for x in ["SAT","ALPHA","CENTR"]):
        macro = "Satellite"
    else:
        macro = "Other TE"

    return name, macro


# ── Cluster insertion sites ────────────────────────────────────────────────────
def cluster_insertions(merged_df, cluster_dist,
                       cond_control="control", cond_treatment="treatment"):
    """
    Cluster nearby insertion sites on same chromosome.
    Returns DataFrame with one row per cluster.
    """
    if merged_df.empty:
        return pd.DataFrame()

    clusters = []
    for chrom, grp in merged_df.groupby("chrom"):
        grp_sorted = grp.sort_values("pos")
        current_cluster = []

        for _, row in grp_sorted.iterrows():
            if not current_cluster:
                current_cluster.append(row)
            elif row["pos"] - current_cluster[-1]["pos"] <= cluster_dist:
                current_cluster.append(row)
            else:
                clusters.append(current_cluster)
                current_cluster = [row]
        if current_cluster:
            clusters.append(current_cluster)

    cluster_rows = []
    for clust in clusters:
        df_c = pd.DataFrame(clust)
        te_fam, te_macro = parse_te_family(
            df_c["te_chrom"].mode()[0] if "te_chrom" in df_c else None)
        cluster_rows.append({
            "chrom":          df_c["chrom"].iloc[0],
            "start":          int(df_c["pos"].min()),
            "end":            int(df_c["pos"].max() + 150),
            "n_reads_total":  len(df_c),
            "n_control":      int((df_c["condition"]==cond_control).sum()),
            "n_treatment":    int((df_c["condition"]==cond_treatment).sum()),
            "te_family":      te_fam,
            "te_macro":       te_macro,
            "samples":        ",".join(sorted(df_c["sample_id"].unique())),
            "strand":         df_c["strand"].mode()[0] if "strand" in df_c else "."
        })

    df_clusters = pd.DataFrame(cluster_rows)
    if not df_clusters.empty:
        df_clusters["somatic_candidate"] = (
            (df_clusters["n_treatment"] >= 1) &
            (df_clusters["n_control"] == 0)
        )
    return df_clusters


# ── Annotate with nearby genes ────────────────────────────────────────────────
def annotate_insertions(df, gtf, outdir):
    """Annotate insertion clusters with nearest gene using bedtools."""
    if df.empty:
        return df

    ins_bed = os.path.join(outdir, "insertion_clusters.bed")
    df[["chrom","start","end","te_family","n_reads_total","strand"]].to_csv(
        ins_bed, sep="\t", header=False, index=False)

    gene_bed = os.path.join(outdir, "genes_annot.bed")
    if not os.path.exists(gene_bed):
        records = []
        gtf_exp = os.path.expanduser(gtf)
        with open(gtf_exp) as f:
            for line in f:
                if line.startswith("#"): continue
                parts = line.strip().split("	")
                if len(parts) < 9 or parts[2] != "gene": continue
                gid  = re.search(r'gene_id "([^"]+)"', parts[8])
                gname= re.search(r'gene_name "([^"]+)"', parts[8])
                if gid:
                    records.append([parts[0], str(int(parts[3])-1), parts[4],
                                    gid.group(1),
                                    gname.group(1) if gname else gid.group(1),
                                    parts[6]])
        with open(gene_bed, "w") as f:
            for r in records:
                f.write("\t".join(r) + "\n")

    annot_out = os.path.join(outdir, "insertions_annotated.bed")
    cmd = (f"bedtools closest -a {ins_bed} -b {gene_bed} "
           f"-D ref -t first > {annot_out} 2>/dev/null")
    subprocess.run(cmd, shell=True)

    try:
        ann = pd.read_csv(annot_out, sep="\t", header=None,
                          names=["chrom","start","end","te_family",
                                 "n_reads","strand",
                                 "g_chrom","g_start","g_end",
                                 "gene_id","gene_name","gene_strand","dist"])
        ann = ann[["chrom","start","end","gene_id","gene_name","dist"]]
        df  = df.merge(ann, on=["chrom","start","end"], how="left")
        log.info(f"Annotated {df['gene_name'].notna().sum()} / {len(df)} clusters")
    except Exception as e:
        log.warning(f"Annotation failed: {e}")

    return df


# ── Main ───────────────────────────────────────────────────────────────────────
def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    log.info("M23 - TE somatic insertion detection")

    cond_control   = args.condition_control
    cond_treatment = args.condition_treatment
    normal_files   = args.chimeric_normal.split(",")
    tumor_files    = args.chimeric_tumor.split(",")

    # 1. Parse chimeric SAMs
    df_normal = parse_chimeric_sam(normal_files, cond_control,
                                    args.min_clip, args.min_mapq)
    df_tumor  = parse_chimeric_sam(tumor_files,  cond_treatment,
                                    args.min_clip, args.min_mapq)

    df_all = pd.concat([df_normal, df_tumor], ignore_index=True)
    log.info(f"Total soft-clips: {len(df_all)} "
             f"({cond_control}: {len(df_normal)}, {cond_treatment}: {len(df_tumor)})")

    if df_all.empty:
        log.warning("No soft-clips found — check --min-clip parameter")
        pd.DataFrame().to_csv(
            os.path.join(args.outdir, "te_insertion_clusters.tsv"),
            sep="\t", index=False)
        sys.exit(0)

    # Save all clips
    df_all.to_csv(os.path.join(args.outdir, "all_softclips.tsv"),
                  sep="\t", index=False)

    # 2. Index TE consensus and align
    te_fa = index_te_consensus(args.te_consensus, args.outdir, args.aligner)
    if te_fa is None:
        log.error("TE indexing failed")
        sys.exit(1)

    te_hits = align_clips_to_te(df_all, te_fa, args.outdir,
                                  args.threads, args.aligner)

    if te_hits.empty:
        log.warning("No clips mapped to TE consensus")
        pd.DataFrame().to_csv(
            os.path.join(args.outdir, "te_insertion_clusters.tsv"),
            sep="\t", index=False)
        sys.exit(0)

    # 3. Merge clip info with TE mapping
    # Build lookup from clip_id to genomic info
    clip_lookup = {}
    for _, row in df_all.iterrows():
        key = (f"{row['read_id']}_{row['condition']}_{row['sample_id']}_"
               f"{row['chrom']}_{row['pos']}_{row['clip_pos']}")
        clip_lookup[key] = row.to_dict()

    merged_rows = []
    for _, hit in te_hits.iterrows():
        info = clip_lookup.get(hit["clip_id"], {})
        if info:
            merged_rows.append({**info,
                                  "te_chrom":  hit["te_chrom"],
                                  "te_start":  hit["te_start"],
                                  "te_mapq":   hit["te_mapq"]})

    merged_df = pd.DataFrame(merged_rows)
    log.info(f"Clips with TE mapping: {len(merged_df)}")
    merged_df.to_csv(os.path.join(args.outdir, "clips_with_te_mapping.tsv"),
                     sep="\t", index=False)

    # 4. Cluster insertion sites
    clusters = cluster_insertions(merged_df, args.cluster_dist,
                                  cond_control, cond_treatment)

    if clusters.empty:
        log.warning("No insertion clusters found")
        sys.exit(0)

    # Filter by min support
    clusters = clusters[clusters["n_reads_total"] >= args.min_support]
    log.info(f"Insertion clusters (>={args.min_support} reads): {len(clusters)}")

    # 5. Annotate
    clusters = annotate_insertions(clusters, args.gtf, args.outdir)

    # 6. Save results
    clusters.to_csv(os.path.join(args.outdir, "te_insertion_clusters.tsv"),
                    sep="\t", index=False)

    somatic = clusters[clusters["somatic_candidate"]]
    somatic.to_csv(os.path.join(args.outdir, "somatic_te_insertions.tsv"),
                   sep="\t", index=False)

    # 7. Summary
    summary = pd.DataFrame([
        {"metric": "Total soft-clips analysed",        "value": len(df_all)},
        {"metric": "Clips mapping to TE consensus",    "value": len(te_hits)},
        {"metric": "Insertion clusters (all)",         "value": len(clusters)},
        {"metric": f"Somatic candidates ({cond_treatment}-only)",  "value": len(somatic)},
        {"metric": "SINE insertions",
         "value": int((clusters["te_macro"]=="SINE").sum())},
        {"metric": "LINE insertions",
         "value": int((clusters["te_macro"]=="LINE").sum())},
        {"metric": "LTR/ERV insertions",
         "value": int((clusters["te_macro"]=="LTR/ERV").sum())},
    ])
    summary.to_csv(os.path.join(args.outdir, "insertion_summary.tsv"),
                   sep="\t", index=False)

    log.info("\nM23 Summary:")
    for _, row in summary.iterrows():
        log.info(f"  {row['metric']}: {row['value']}")

    if len(somatic) > 0:
        log.info(f"\nTop somatic TE insertion candidates:")
        top = somatic.sort_values("n_treatment", ascending=False).head(10)
        for _, row in top.iterrows():
            gene = row.get("gene_name", "intergenic")
            log.info(f"  {row['chrom']}:{row['start']}-{row['end']} | "
                     f"{row['te_family']} ({row['te_macro']}) | "
                     f"gene: {gene} | "
                     f"{cond_treatment} reads: {row['n_treatment']}")

    log.info(f"[OK] M23 complete. Output: {args.outdir}")


if __name__ == "__main__":
    main()
