#!/usr/bin/env python3
"""S6: versions + key commands + performance (three-block xlsx + tsv).

--edition s6 (default) reproduces the original Table S6 unchanged.
--edition s10 builds Table S10: adds the rows of tools_extra.tsv and applies
the S10 corrections (see S10_* below).

usage (Table S10 of the revised manuscript, from the repository root):
  python3 scripts/make_tableS10.py --edition s10 --perf docs/table_s10/perf_table.tsv \
      --extra-tools docs/table_s10/tools_extra.tsv --out <output directory>
"""
import pandas as pd, openpyxl, os, argparse
from openpyxl.styles import Font, Alignment, PatternFill, Border, Side
from openpyxl.utils import get_column_letter

ap = argparse.ArgumentParser()
ap.add_argument("--out", required=True,
                help="output directory")
ap.add_argument("--perf", required=True,
                help="performance table for Block 3 (docs/table_s10/perf_table.tsv)")
ap.add_argument("--edition", choices=["s6", "s10"], default="s6")
ap.add_argument("--extra-tools", required=True,
                help="extra Block 1 rows for --edition s10 (docs/table_s10/tools_extra.tsv; "
                     "columns: after, category, tool, version, note)")
args = ap.parse_args()
OUT = args.out

# ── style helpers ─────────────────────────────────────────────────────────────
H1_FILL   = PatternFill("solid", fgColor="08306B")   # dark navy
H2_FILL   = PatternFill("solid", fgColor="2171B5")   # medium blue
H3_FILL   = PatternFill("solid", fgColor="6BAED6")   # light blue
GREY_FILL = PatternFill("solid", fgColor="F2F2F2")
WHT_FONT  = Font(color="FFFFFF", bold=True, name="Calibri", size=10)
BLD_FONT  = Font(bold=True, name="Calibri", size=10)
NRM_FONT  = Font(name="Calibri", size=10)
MON_FONT  = Font(name="Courier New", size=9)

thin = Side(style="thin", color="BBBBBB")
BRD  = Border(top=thin, bottom=thin, left=thin, right=thin)

def c(ws, row, col, val, font=None, fill=None, align="left", wrap=False, mono=False):
    cell = ws.cell(row=row, column=col, value=val)
    cell.font = font or (MON_FONT if mono else NRM_FONT)
    if fill: cell.fill = fill
    cell.alignment = Alignment(horizontal=align, vertical="center", wrap_text=wrap)
    cell.border = BRD
    return cell

def hdr_row(ws, row, labels, fill, n_cols=4):
    for i, lbl in enumerate(labels, 1):
        c(ws, row, i, lbl, font=WHT_FONT, fill=fill, align="center")
    for j in range(len(labels)+1, n_cols+1):
        c(ws, row, j, None, fill=fill)

def section_title(ws, row, text, fill, n_cols=4):
    ws.merge_cells(start_row=row, start_column=1, end_row=row, end_column=n_cols)
    cell = ws.cell(row=row, column=1, value=text)
    cell.font = Font(bold=True, color="FFFFFF", name="Calibri", size=11)
    cell.fill = fill
    cell.alignment = Alignment(horizontal="left", vertical="center")
    cell.border = BRD

# ── Block 1: Tool versions ────────────────────────────────────────────────────
VERSIONS = [
    # (Category, Tool/DB, Version, Note)
    ("Read trimming",   "fastp",              "1.1.0",  "verified: fastp --version"),
    ("Alignment",       "STAR",               "2.7.10b","verified: STAR --version; compiled 2022-11-01"),
    ("Alignment",       "HISAT2",             "2.2.2",  "benchmark comparison; verified: hisat2 --version"),
    ("TE quantification","bedtools multicov", "2.31.1", "verified: bedtools --version"),
    ("Gene quantification","featureCounts",   "2.1.1",  "subread package; verified: featureCounts -v"),
    ("Satellite quantification","bowtie2",    "2.5.4",  "verified: bowtie2 --version"),
    ("Splicing",        "rMATS",              "4.3.0",  "verified: rmats.py --version"),
    ("Splicing",        "SUPPA2",             "2.4",    "verified: suppa.py --version"),
    ("Gene quant bench","salmon",             "2.2.1",  "bioconda package hfa8f182_0; verified: salmon --version"),
    ("Gene quant bench","tximport",           "1.38.2", "Bioconductor; verified: packageVersion('tximport')"),
    ("TE quant bench",  "Telescope",          "1.0.3",  "Python 3.7.12; verified: telescope --version"),
    ("RNA editing",     "JACUSA2",            "2.1.16", "verified: jar manifest; conda env mosaic"),
    ("RNA editing",     "REDItools3",         "3.7",    "verified: importlib.metadata; conda env mosaic"),
    ("Statistics",      "R",                  "4.5.1",  "2025-06-13; verified: R.version.string. "
                                                        "Published results generated with R 4.5.1; "
                                                        "public repository (github.com/WalterIrib/MOSAIC) pins R 4.5.2 "
                                                        "for improved cross-platform Bioconductor compatibility. "
                                                        "A ±5-family tolerance in check_canonical.py accounts for "
                                                        "expected floating-point differences in DESeq2 locfit "
                                                        "dispersion fitting between R patch versions — see CHANGELOG."),
    ("Statistics",      "DESeq2",             "1.50.2", "Bioconductor; verified: packageVersion('DESeq2')"),
    ("Statistics",      "ashr",               "2.2.63", "lfcShrink type='ashr'; verified: packageVersion('ashr')"),
    ("General",         "samtools",           "1.23",   "verified: samtools --version"),
    ("General",         "Python",             "3.11.9", "conda env mosaic; verified: python --version"),
    ("Databases",       "RepeatMasker hg38",  "UCSC rmsk","4,519,330 loci; file date 2026-03-10"),
    ("Databases",       "GENCODE annotation", "v49",    "Homo sapiens; GTF + FASTA"),
    ("Databases",       "DFAM",               "3.9",    "used by RepeatMasker for TE classification"),
]

# ── Block 2: Key commands / parameters ────────────────────────────────────────
COMMANDS = [
    # (Analysis step, Tool, Key command / parameters)
    ("P1.1a — Unique-only vs multi-mapping",
     "bedtools multicov + DESeq2",
     "bash locus_control/join_and_deseq2.sh\n"
     "  bedtools multicov -bams *.bam -bed TE.bed [-q 255 for unique]\n"
     "  DESeq2: design=~condition; lfcShrink(type='ashr'); padj<0.05; |LFC|≥1",
     ""),
    ("P1.1b — Permutation FDR (genes)",
     "DESeq2 (C(6,3)=20 permutations)",
     "conda run -n mosaic Rscript perm_fdr/run_permutation_fdr.R\n"
     "  18 non-trivial permutations; 0/18 with ≥1 significant gene; empirical FDR=0.000",
     ""),
    ("P1.1b — Permutation FDR (TE loci)",
     "DESeq2 (C(6,3)=20 permutations)",
     "conda run -n mosaic Rscript perm_fdr/run_permutation_fdr_te.R\n"
     "  18 non-trivial permutations; 6/18 with ≥1 significant locus; max 999; FDR=1.26%",
     ""),
    ("P1.2 — REDItools3 vs JACUSA2",
     "REDItools3 + JACUSA2 2.1.16",
     "bash /tmp/run_p12_reditools.sh   # chr1 only; 2 threads\n"
     "conda run -n mosaic python bench_editing/compare_reditools_jacusa.py\n"
     "  JACUSA2 dEF: G/(A+G) for + strand; C/(C+T) for - strand",
     ""),
    ("P1.3 — Telescope vs MOSAIC M12",
     "Telescope 1.0.3 (EM)",
     "bash /tmp/run_p13_telescope.sh   # chr1+chr2 TE-filtered BAMs\n"
     "  Direct binary: telescope_env37/bin/telescope --ncpu 1\n"
     "  samtools view -b BAM chr1 chr2 | telescope assign\n"
     "conda run -n mosaic Rscript bench_te/process_telescope_lfc.R\n"
     "  Caveat: chr-based BAM filter excludes multi-mapping secondary alignments on chr3–22",
     ""),
    ("P1.4 — SUPPA2 vs rMATS",
     "SUPPA2 2.4 + rMATS 4.3.0",
     "bash bench_splicing/run_suppa2.sh\n"
     "  7 event types (SE,A5,A3,MX,RI,AF,AL); empirical diffSplice; pval<0.05; |dPSI|>0.1\n"
     "  TPM matrix header: NO 'transcript_id' label; PSI split at data position 0",
     "Metric: Directional concordance on co-testable events (non-NaN SUPPA2 ΔΨ)\n"
     "290/385 rMATS-significant events (75%) concordant in direction with SUPPA2\n"
     "(binomial P<10⁻²³, exact P=2.3×10⁻²⁴); ΔΨ Pearson r=0.27.\n"
     "Event counts not directly comparable (rMATS ~46k vs SUPPA2 ~165k\n"
     "tested events, different multiple-testing regimes).\n"
     "Agreement by event type: RI 86% (P~10⁻¹⁸), A3SS 79%, SE 73%,\n"
     "A5SS 62% (ns), MXE 58% (ns)."),
    ("P1.5 — Salmon vs MOSAIC/STAR",
     "salmon 2.2.1 + tximport 1.38.2",
     "salmon quant -p 4 --gcBias --validateMappings -l A\n"
     "conda run -n mosaic Rscript bench_genes_salmon/run_deseq2_salmon.R\n"
     "  tximport(ignoreTxVersion=FALSE) required for versioned IDs\n"
     "  85.3% concordance with MOSAIC/STAR; LFC r=0.87 on MOSAIC-sig genes",
     ""),
    ("M12 — MOSAIC TE module (main pipeline)",
     "bedtools multicov + DESeq2",
     "bedtools multicov -bams *.bam -bed RepeatMasker_hg38_rmsk.bed\n"
     "  DESeq2: stable-gene size factors; lfcShrink(type='ashr')\n"
     "  Threshold: padj<0.05; |log₂FC_shrunk|≥1; 4,519,330 loci tested",
     ""),
    ("M17 — rMATS (main pipeline)",
     "rMATS 4.3.0",
     "rmats.py --b1 ctrl.txt --b2 mito.txt --gtf gencode.v49.gtf\n"
     "  --od outdir -t paired --readLength 150 --novelSS\n"
     "  FDR<0.05; |dPSI|>0.1",
     ""),
    ("M19 — JACUSA2 (main pipeline)",
     "JACUSA2 2.1.16",
     "java -jar JACUSA_v2.1.16.jar call-2 -p 20 ...\n"
     "  dEF = G/(A+G) [+ strand] or C/(C+T) [- strand]; padj<0.05; |dEF|>0.1",
     ""),
]

# ── Block 3: Performance (from perf_table.tsv) ────────────────────────────────
perf_df = pd.read_csv(args.perf, sep="\t")

LABEL = "Table S6"

# ── S10 edition ───────────────────────────────────────────────────────────────
if args.edition == "s10":
    LABEL = "Table S10"

    # Block 1: per-tool corrections to existing rows
    S10_VERSION_FIXES = {
        "bowtie2":  {"cat": "TE and satellite quantification",
                     "note": "M05; verified: bowtie2 --version"},
        "tximport": {"note": "Bioconductor; verified: packageVersion('tximport'); "
                             "installed in the user R library ~/R/library (built under R 4.5.1)"},
        "R":        {"note": "2025-06-13; verified: R.version.string. "
                             "Published results generated with R 4.5.1; "
                             "the repository environment pins R 4.5.1 (environment.yml, MOSAIC v1.1.8)"},
    }
    VERSIONS = [
        (S10_VERSION_FIXES.get(t, {}).get("cat", cat), t, v,
         S10_VERSION_FIXES.get(t, {}).get("note", n))
        for cat, t, v, n in VERSIONS
    ]
    # bowtie2 (M05) listed next to bedtools, before featureCounts
    order = [r[1] for r in VERSIONS]
    VERSIONS.insert(order.index("featureCounts"), VERSIONS.pop(order.index("bowtie2")))

    # Block 1: extra rows, each inserted after the tool named in column 'after'
    extra = pd.read_csv(args.extra_tools, sep="\t", dtype=str).fillna("")
    for _, e in extra.iterrows():
        pos = [r[1] for r in VERSIONS].index(e["after"]) + 1
        VERSIONS.insert(pos, (e["category"], e["tool"], e["version"], e["note"]))

    # Block 2: commands as implemented in the code
    S10_COMMAND_FIXES = {
        "M12 — MOSAIC TE module (main pipeline)":
            "bedtools multicov -q 1 -bams *.bam -bed RepeatMasker_hg38_rmsk.bed\n"
            "  DESeq2: stable-gene size factors; lfcShrink(type='ashr')\n"
            "  Threshold: padj<0.05; |log₂FC_shrunk|≥1\n"
            "  2,630,539 loci tested (non-NA padj; independentFiltering=FALSE)\n"
            "  Canonical counts predate the v1.1.0 threshold-precision fix (Table S8; Fig. S8 legend)",
        "M19 — JACUSA2 (main pipeline)":
            "java -jar JACUSA_v2.1.16.jar call-2 -p 20 -c 5 ...\n"
            "  Filters: JACUSA2 score≥1.5; EF (treated)≥0.05; coverage≥5\n"
            "  dEF = G/(A+G) [+ strand] or C/(C+T) [- strand]; |dEF|>0.05 = increased/decreased",
    }
    COMMANDS = [(s, t, S10_COMMAND_FIXES.get(s, cmd), m) for s, t, cmd, m in COMMANDS]

    # Block 3: drop the retired M09 module and the report re-render rows (S1–S6)
    perf_df = perf_df[~perf_df["module"].isin(["m09_satellite", "S1", "S2", "S3", "S4", "S5", "S6"])]

# ── Build workbook ─────────────────────────────────────────────────────────────
wb = openpyxl.Workbook()
ws = wb.active
ws.title = LABEL

# Title
ws.merge_cells("A1:D1")
t = ws.cell(row=1, column=1,
            value=f"{LABEL}. Software versions, key commands, and computational performance")
t.font = Font(bold=True, name="Calibri", size=12)
t.alignment = Alignment(horizontal="left", vertical="center")
ws.row_dimensions[1].height = 26

cur = 3  # current row

# ── BLOCK 1: Versions ─────────────────────────────────────────────────────────
section_title(ws, cur, "Block 1 — Software and database versions", H1_FILL, n_cols=4)
cur += 1
hdr_row(ws, cur, ["Category", "Tool / Database", "Version", "Verification / notes"],
        H2_FILL, n_cols=4)
cur += 1
for i, (cat, tool, ver, note) in enumerate(VERSIONS):
    fill = GREY_FILL if i % 2 == 0 else None
    c(ws, cur, 1, cat,  fill=fill, font=BLD_FONT)
    c(ws, cur, 2, tool, fill=fill, font=BLD_FONT)
    c(ws, cur, 3, ver,  fill=fill, align="center", mono=True)
    c(ws, cur, 4, note, fill=fill, wrap=True)
    ws.row_dimensions[cur].height = 16
    cur += 1

cur += 1  # blank spacer

# ── BLOCK 2: Commands ─────────────────────────────────────────────────────────
section_title(ws, cur, "Block 2 — Key analysis commands and parameters", H1_FILL, n_cols=4)
cur += 1
hdr_row(ws, cur, ["Analysis step", "Tool",
                   "Key command / parameters", "Metric(s) / benchmark outcome"],
        H2_FILL, n_cols=4)
cur += 1
for i, (step, tool, cmd, metric) in enumerate(COMMANDS):
    fill = GREY_FILL if i % 2 == 0 else None
    c(ws, cur, 1, step,   fill=fill, font=BLD_FONT, wrap=True)
    c(ws, cur, 2, tool,   fill=fill)
    c(ws, cur, 3, cmd,    fill=fill, wrap=True, mono=True)
    c(ws, cur, 4, metric if metric else None, fill=fill, wrap=True)
    h = max(cmd.count("\n"), metric.count("\n") if metric else 0) * 14 + 14
    ws.row_dimensions[cur].height = max(16, h)
    cur += 1

cur += 1

# ── BLOCK 3: Performance ──────────────────────────────────────────────────────
section_title(ws, cur, "Block 3 — Computational performance", H1_FILL, n_cols=4)
cur += 1
perf_cols = ["module", "tool", "wall_hms_total", "note"]
perf_labels = ["Module", "Tool (version)", "Wall time (total)", "Notes / peak RAM"]
hdr_row(ws, cur, perf_labels, H2_FILL, n_cols=4)
cur += 1

for i, (_, row) in enumerate(perf_df.iterrows()):
    fill = GREY_FILL if i % 2 == 0 else None
    c(ws, cur, 1, str(row["module"]), fill=fill, font=BLD_FONT)
    c(ws, cur, 2, str(row["tool"]),   fill=fill)
    wt = str(row["wall_hms_total"]) if pd.notna(row["wall_hms_total"]) else "N/A"
    c(ws, cur, 3, wt, fill=fill, align="center")
    c(ws, cur, 4, str(row["note"]) if pd.notna(row["note"]) else "", fill=fill, wrap=True)
    ws.row_dimensions[cur].height = 16
    cur += 1

# ── Column widths ─────────────────────────────────────────────────────────────
ws.column_dimensions["A"].width = 36
ws.column_dimensions["B"].width = 30
ws.column_dimensions["C"].width = 45
ws.column_dimensions["D"].width = 48

# ── Save ──────────────────────────────────────────────────────────────────────
BASENAME = "TableS10_tool_versions" if args.edition == "s10" else "tableS6_versions_perf"
xlsx_path = os.path.join(OUT, f"{BASENAME}.xlsx")
wb.save(xlsx_path)
print(f"[saved] {xlsx_path}")

# TSV: three sections concatenated
rows_tsv = []
rows_tsv.append(["BLOCK1_VERSIONS","","",""])
rows_tsv.append(["Category","Tool","Version","Notes"])
for cat, tool, ver, note in VERSIONS:
    rows_tsv.append([cat, tool, ver, note])
rows_tsv.append(["","","",""])
rows_tsv.append(["BLOCK2_COMMANDS","","",""])
rows_tsv.append(["Analysis_step","Tool","Command","Metric/benchmark_outcome"])
for step, tool, cmd, metric in COMMANDS:
    rows_tsv.append([step, tool, cmd.replace("\n"," | "), metric.replace("\n"," | ")])
rows_tsv.append(["","","",""])
rows_tsv.append(["BLOCK3_PERFORMANCE","","",""])
rows_tsv.append(list(perf_df.columns))
for _, row in perf_df.iterrows():
    rows_tsv.append([str(v) for v in row])

tsv_path = os.path.join(OUT, f"{BASENAME}.tsv")
with open(tsv_path, "w") as f:
    for r in rows_tsv:
        f.write("\t".join(str(x) for x in r) + "\n")
print(f"[saved] {tsv_path}")
