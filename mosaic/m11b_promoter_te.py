"""
mosaic/m11b_promoter_te.py
M11b — TE insertions in promoter regions (±2kb upstream of TSS).

Strategy:
1. Build promoter BED from GENCODE GTF (2kb upstream of all TSSs)
2. Intersect with RepeatMasker hg38 BED (from M11)
3. Annotate genes with promoter TEs
4. Output: te_promoter_genes.tsv, te_promoter_counts.tsv

Uses ENCODE-style promoter definition: TSS ±2kb (or upstream-only variant).
Falls back to FANTOM5/EPDnew curated promoter database if available.
"""
from __future__ import annotations
import csv, subprocess, tempfile
from pathlib import Path
from typing import Optional
from .log import info, warn, error, success, section


# ─────────────────────────────────────────────────────────────────────────────
# Build promoter BED from GTF
# ─────────────────────────────────────────────────────────────────────────────

def build_promoter_bed(
    gtf: Path,
    out_bed: Path,
    upstream: int   = 2000,
    downstream: int = 200,
    genome_fai: Optional[Path] = None,
) -> Path:
    """
    Extract TSS coordinates from GTF and define promoter windows.
    Promoter = [TSS - upstream, TSS + downstream], strand-aware.
    upstream=2000, downstream=200 is a common conservative definition.

    Writes a 6-column BED: chr, start, end, gene_id, score, strand
    """
    if out_bed.exists():
        info(f"Promoter BED already exists: {out_bed}")
        return out_bed

    info(f"Building promoter BED from {gtf.name} (±{upstream}/{downstream}bp)...")
    out_bed.parent.mkdir(parents=True, exist_ok=True)

    # Read chromosome sizes from .fai if provided (for clamping coords)
    chrom_sizes: dict[str, int] = {}
    if genome_fai and genome_fai.exists():
        with open(genome_fai) as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) >= 2:
                    chrom_sizes[parts[0]] = int(parts[1])

    seen_tss = set()
    records  = []

    with open(gtf) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9 or parts[2] != "gene":
                continue

            chrom  = parts[0]
            strand = parts[6]
            start  = int(parts[3]) - 1  # GTF is 1-based
            end    = int(parts[4])

            # Extract gene_id
            gene_id = ""
            for attr in parts[8].split(";"):
                attr = attr.strip()
                if attr.startswith("gene_id"):
                    gene_id = attr.split('"')[1] if '"' in attr else attr.split()[1]
                    break

            # TSS
            tss = start if strand == "+" else end - 1
            key = (chrom, tss, strand)
            if key in seen_tss:
                continue
            seen_tss.add(key)

            # Promoter window
            if strand == "+":
                prom_start = max(0, tss - upstream)
                prom_end   = tss + downstream
            else:
                prom_start = max(0, tss - downstream)
                prom_end   = tss + upstream

            # Clamp to chromosome size
            if chrom in chrom_sizes:
                prom_end = min(prom_end, chrom_sizes[chrom])

            if prom_start >= prom_end:
                continue

            records.append((chrom, prom_start, prom_end, gene_id, ".", strand))

    # Sort and write
    records.sort(key=lambda x: (x[0], x[1]))
    with open(out_bed, "w") as f:
        for r in records:
            f.write("\t".join(str(x) for x in r) + "\n")

    info(f"Promoter BED: {len(records):,} TSS windows → {out_bed}")
    return out_bed


# ─────────────────────────────────────────────────────────────────────────────
# Intersect TE loci with promoters
# ─────────────────────────────────────────────────────────────────────────────

def intersect_te_promoters(
    rmsk_bed: Path,
    promoter_bed: Path,
    out_intersect: Path,
    min_overlap_frac: float = 0.10,
) -> Path:
    """
    bedtools intersect: find TE loci overlapping promoter regions.
    -F min_overlap_frac: fraction of TE that must overlap promoter.
    Writes: chr_te, start_te, end_te, name_te, score_te, strand_te,
            chr_prom, start_prom, end_prom, gene_id, ., strand_prom
    """
    if out_intersect.exists():
        info(f"TE-promoter intersect already exists: {out_intersect}")
        return out_intersect

    out_intersect.parent.mkdir(parents=True, exist_ok=True)
    info(f"bedtools intersect: TE ∩ promoters (min overlap {min_overlap_frac:.0%})...")

    cmd = [
        "bedtools", "intersect",
        "-a", str(rmsk_bed),
        "-b", str(promoter_bed),
        "-wa", "-wb",
        "-F", str(min_overlap_frac),
        "-sorted",
    ]

    with open(out_intersect, "w") as f:
        result = subprocess.run(cmd, stdout=f, stderr=subprocess.PIPE, text=True)

    if result.returncode != 0:
        # Try without -sorted (files may not be sorted)
        warn("bedtools -sorted failed, retrying without sorting flag...")
        cmd.remove("-sorted")
        with open(out_intersect, "w") as f:
            result = subprocess.run(cmd, stdout=f, stderr=subprocess.PIPE, text=True)
        if result.returncode != 0:
            error(f"bedtools intersect failed: {result.stderr}")
            raise RuntimeError("bedtools intersect failed")

    n_lines = sum(1 for _ in open(out_intersect))
    info(f"TE-promoter intersections: {n_lines:,}")
    return out_intersect


# ─────────────────────────────────────────────────────────────────────────────
# Classify TE macrocategory
# ─────────────────────────────────────────────────────────────────────────────

def classify_te_macro(te_name: str) -> str:
    n = te_name.upper()
    if any(x in n for x in ["ALU", "MIR", "SINE"]):
        return "SINE"
    if any(x in n for x in ["L1", "L2", "LINE", "CR1"]):
        return "LINE"
    if any(x in n for x in ["LTR", "HERV", "ERV", "MLT", "THE"]):
        return "LTR/ERV"
    if any(x in n for x in ["DNA", "TIGGER", "CHARLIE", "HAT", "MARINER"]):
        return "DNA transposon"
    if "SVA" in n:
        return "SVA"
    if any(x in n for x in ["SAT", "ALPHA", "CENTR", "ACRO"]):
        return "Satellite"
    return "Other TE"


# ─────────────────────────────────────────────────────────────────────────────
# Build promoter TE gene table
# ─────────────────────────────────────────────────────────────────────────────

def build_promoter_te_table(
    intersect_bed: Path,
    out_genes: Path,
    out_summary: Path,
) -> tuple[Path, Path]:
    """
    Parse bedtools intersect output and build:
    1. te_promoter_genes.tsv: one row per gene with TE annotations
    2. te_promoter_summary.tsv: count of genes per TE macrocategory
    """
    # gene_id -> {set of TE names, set of macros, count}
    gene_te: dict[str, dict] = {}
    macro_counts: dict[str, set] = {}

    with open(intersect_bed) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 10:
                continue
            te_name  = parts[3]  # TE name from rmsk BED
            gene_id  = parts[9]  # gene_id from promoter BED
            macro    = classify_te_macro(te_name)

            if gene_id not in gene_te:
                gene_te[gene_id] = {
                    "gene_id":       gene_id,
                    "n_te_loci":     0,
                    "te_names":      set(),
                    "te_macros":     set(),
                    "chr":           parts[0],
                    "prom_start":    parts[7],
                    "prom_end":      parts[8],
                    "strand":        parts[11] if len(parts) > 11 else ".",
                }
            gene_te[gene_id]["n_te_loci"]  += 1
            gene_te[gene_id]["te_names"].add(te_name)
            gene_te[gene_id]["te_macros"].add(macro)

            if macro not in macro_counts:
                macro_counts[macro] = set()
            macro_counts[macro].add(gene_id)

    # Write gene table
    out_genes.parent.mkdir(parents=True, exist_ok=True)
    with open(out_genes, "w", newline="") as f:
        fieldnames = ["gene_id", "chr", "prom_start", "prom_end", "strand",
                      "n_te_loci", "n_te_families", "te_macrocategories", "te_names_sample"]
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        for gid, d in sorted(gene_te.items()):
            w.writerow({
                "gene_id":           gid,
                "chr":               d["chr"],
                "prom_start":        d["prom_start"],
                "prom_end":          d["prom_end"],
                "strand":            d["strand"],
                "n_te_loci":         d["n_te_loci"],
                "n_te_families":     len(d["te_names"]),
                "te_macrocategories": ",".join(sorted(d["te_macros"])),
                "te_names_sample":   ",".join(sorted(d["te_names"])[:5]),
            })

    # Write summary
    with open(out_summary, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["te_macrocategory", "n_genes_with_promoter_te"])
        for macro, genes in sorted(macro_counts.items(), key=lambda x: -len(x[1])):
            w.writerow([macro, len(genes)])

    info(f"Promoter TE genes: {len(gene_te):,}")
    info(f"TE macrocategory distribution: {dict((k,len(v)) for k,v in macro_counts.items())}")
    return out_genes, out_summary


# ─────────────────────────────────────────────────────────────────────────────
# Main M11b runner
# ─────────────────────────────────────────────────────────────────────────────

def run_m11b(cfg) -> dict:
    """
    M11b pipeline:
    1. Build promoter BED from GENCODE GTF
    2. Intersect TE loci (rmsk) with promoter regions
    3. Build gene-level promoter TE annotation table
    Returns dict with output paths.
    """
    section("M11b — TE in promoter regions")

    out_dir  = Path(cfg.output_dir) / "m11b"
    m11_dir  = Path(cfg.output_dir) / "m11"
    out_dir.mkdir(parents=True, exist_ok=True)

    gtf       = Path(cfg.gtf)
    rmsk_bed  = m11_dir / "rmsk_hg38.bed"
    genome_fa = Path(cfg.genome_fa)
    genome_fai = Path(str(genome_fa) + ".fai")

    if not rmsk_bed.exists():
        error(f"rmsk_hg38.bed not found at {rmsk_bed} — run M11 first")
        raise FileNotFoundError(f"Missing {rmsk_bed}")

    upstream   = getattr(cfg, "promoter_upstream",   2000)
    downstream = getattr(cfg, "promoter_downstream",  200)
    min_ovlp   = getattr(cfg, "promoter_min_overlap", 0.10)

    # 1. Promoter BED
    promoter_bed = out_dir / "promoters_2kb.bed"
    build_promoter_bed(
        gtf=gtf,
        out_bed=promoter_bed,
        upstream=upstream,
        downstream=downstream,
        genome_fai=genome_fai if genome_fai.exists() else None,
    )

    # 2. Intersect
    intersect_bed = out_dir / "te_promoter_intersect.bed"
    intersect_te_promoters(
        rmsk_bed=rmsk_bed,
        promoter_bed=promoter_bed,
        out_intersect=intersect_bed,
        min_overlap_frac=min_ovlp,
    )

    # 3. Build tables
    te_promoter_genes   = out_dir / "te_promoter_genes.tsv"
    te_promoter_summary = out_dir / "te_promoter_summary.tsv"
    build_promoter_te_table(intersect_bed, te_promoter_genes, te_promoter_summary)

    # 4. Write universe file for GO enrichment
    universe_file = out_dir / "te_promoter_universe.txt"
    gene_ids = set()
    with open(te_promoter_genes) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            gene_ids.add(row["gene_id"])
    universe_file.write_text("\n".join(sorted(gene_ids)))
    info(f"Promoter TE universe: {len(gene_ids):,} genes")

    success(f"M11b complete — {len(gene_ids):,} genes with promoter TE insertions")
    return {
        "promoter_bed":       promoter_bed,
        "intersect_bed":      intersect_bed,
        "te_promoter_genes":  te_promoter_genes,
        "te_promoter_summary": te_promoter_summary,
        "universe_file":      universe_file,
        "n_genes":            len(gene_ids),
    }
