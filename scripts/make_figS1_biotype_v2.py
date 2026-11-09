#!/usr/bin/env python3
"""
make_figS1_biotype_v2.py — Supplementary Figure S1 (revision)
4-panel stacked bar (mean ± SD per condition) biotype composition.
Panels A,B = METACAR ChaRNA; C,D = De Cecco senescence.
A,C = full (0-100%); B,D = non-protein-coding zoom (~0-12%).

Values are percentages of the per-sample sum of RAW counts: featureCounts gene
counts (all_genes_counts.tsv, RN7SL* excluded) + M05 bowtie2 TE family counts
(7SLRNA -> other_ncRNA) + M05 satellite counts. Unmapped reads are not counted.

v2 (2026-11-09), from make_figS1_biotype.py (md5 42e09bb8, 2026-06-22):
  - gene_type rRNA and Mt_rRNA -> separate class 'rRNA' (was other_ncRNA);
    rRNA_pseudogene stays in pseudogene, Mt_tRNA stays in other_ncRNA;
  - y-axis label '% of summed raw counts' (was '% normalized counts');
  - --outdir (default ~/Desktop/DECECCA2019/figS1_v2); outputs figS1_v2.{png,pdf}
    and figS1_v2_fractions.tsv (mean, SD per dataset/condition/class).
Everything else is unchanged.

Usage: python make_figS1_biotype_v2.py [--outdir DIR]
"""
import argparse
import csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from pathlib import Path
from collections import defaultdict

# ─── Paths ───────────────────────────────────────────────────────────────────
META   = Path("/home/walter/Desktop/METACAR/mosaic_output")
DECCA  = Path("/home/walter/Desktop/DECECCA2019/mosaic_output")
_ap = argparse.ArgumentParser(description=__doc__.split('\n')[1])
_ap.add_argument('--outdir', default="/home/walter/Desktop/DECECCA2019/figS1_v2")
OUTDIR = Path(_ap.parse_args().outdir)
OUTDIR.mkdir(parents=True, exist_ok=True)

# ─── Biotype grouping ────────────────────────────────────────────────────────
BIOTYPE_MAP = {
    'protein_coding':                       'protein_coding',
    'lncRNA':                               'lncRNA',
    'processed_pseudogene':                 'pseudogene',
    'unprocessed_pseudogene':               'pseudogene',
    'transcribed_unprocessed_pseudogene':   'pseudogene',
    'transcribed_processed_pseudogene':     'pseudogene',
    'transcribed_unitary_pseudogene':       'pseudogene',
    'rRNA_pseudogene':                      'pseudogene',
    'IG_V_pseudogene':                      'pseudogene',
    'TR_V_pseudogene':                      'pseudogene',
    'IG_C_pseudogene':                      'pseudogene',
    'unitary_pseudogene':                   'pseudogene',
    'snRNA':                                'snRNA/snoRNA',
    'snoRNA':                               'snRNA/snoRNA',
    'scaRNA':                               'snRNA/snoRNA',
    'miRNA':                                'miRNA',
    'rRNA':                                 'rRNA',
    'Mt_rRNA':                              'rRNA',
    # misc_RNA → other_ncRNA (dominated by RN7SL* which is excluded separately)
    'misc_RNA':                             'other_ncRNA',
}
def map_biotype(gt):
    if gt in BIOTYPE_MAP:
        return BIOTYPE_MAP[gt]
    return 'other_ncRNA'

# Display order (bottom → top in stacked bar)
CATEGORIES_FULL = [
    'protein_coding', 'TE', 'lncRNA',
    'pseudogene', 'snRNA/snoRNA', 'miRNA', 'satellite', 'rRNA', 'other_ncRNA',
]
CATEGORIES_NC = [c for c in CATEGORIES_FULL if c != 'protein_coding']

# ─── Colors ──────────────────────────────────────────────────────────────────
COLORS = {
    'protein_coding': '#4D82B8',
    'TE':             '#E07B00',
    'lncRNA':         '#2ECC71',
    'pseudogene':     '#9B59B6',
    'snRNA/snoRNA':   '#1ABC9C',
    'miRNA':          '#E74C3C',
    'satellite':      '#8B4513',   # brown — distinct from lncRNA green
    'rRNA':           '#F1C40F',   # yellow — distinct from all other classes
    'other_ncRNA':    '#95A5A6',
}

# Genes to exclude entirely from biotype composition (technical artefacts)
# RN7SL* = 7SL RNA / SRP — ribo-depletion inefficiency, not biological signal
EXCLUDE_GENE_PREFIXES = ('RN7SL',)

# ─── Load helpers ────────────────────────────────────────────────────────────
def load_gene_type_map(path):
    """gene_id → biotype_group"""
    gm = {}
    with open(path) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            gm[row['gene_id']] = map_biotype(row['gene_type'])
    return gm

def load_count_matrix(path):
    """Returns {sample_id: {row_id: count}}"""
    mat = defaultdict(dict)
    with open(path) as f:
        reader = csv.DictReader(f, delimiter='\t')
        id_col = reader.fieldnames[0]   # 'Geneid' or 'gene_id' or 'family_id'
        samples = reader.fieldnames[1:]
        for row in reader:
            rid = row[id_col]
            for s in samples:
                mat[s][rid] = float(row[s] or 0)
    return dict(mat), samples

def aggregate_biotype(gene_counts, gene_type_map, te_counts, sat_counts,
                      gene_name_map=None, exclude_te=None):
    """
    gene_counts:   {gene_id: count}
    te_counts:     {family_id: count}  (7SLRNA → other_ncRNA)
    sat_counts:    {family_id: count}
    gene_name_map: {gene_id: gene_name}  for prefix exclusion
    Returns {biotype_group: total_count}
    """
    agg = defaultdict(float)
    for gid, cnt in gene_counts.items():
        # Exclude RN7SL* genes (SRP RNA — ribo-depletion artefact)
        if gene_name_map:
            gname = gene_name_map.get(gid, '')
            if any(gname.startswith(pfx) for pfx in EXCLUDE_GENE_PREFIXES):
                continue
        bg = gene_type_map.get(gid, 'other_ncRNA')
        agg[bg] += cnt
    for fid, cnt in te_counts.items():
        if fid == '7SLRNA':
            agg['other_ncRNA'] += cnt
        else:
            agg['TE'] += cnt
    for fid, cnt in sat_counts.items():
        agg['satellite'] += cnt
    return dict(agg)

def compute_fractions(agg_by_sample):
    """Returns {sample: {biotype: fraction%}} and total per sample"""
    frac = {}
    for sid, agg in agg_by_sample.items():
        total = sum(agg.values())
        frac[sid] = {bg: (cnt / total * 100) if total > 0 else 0.0
                     for bg, cnt in agg.items()}
    return frac

def group_mean_sd(frac_by_sample, condition_map):
    """
    condition_map: {sample_id: condition}
    Returns {condition: {biotype: (mean, sd)}}
    """
    cond_vals = defaultdict(lambda: defaultdict(list))
    for sid, fracs in frac_by_sample.items():
        cond = condition_map.get(sid)
        if cond is None:
            continue
        for bg, pct in fracs.items():
            cond_vals[cond][bg].append(pct)
    result = {}
    for cond, bg_dict in cond_vals.items():
        result[cond] = {}
        for bg, vals in bg_dict.items():
            result[cond][bg] = (np.mean(vals), np.std(vals, ddof=1) if len(vals) > 1 else 0.0)
    return result

# ─── Load METACAR ────────────────────────────────────────────────────────────
print("[INF] Loading METACAR data...")
meta_gene_map = load_gene_type_map(META / "integrated/gene_type_map.tsv")
meta_gene_name_map = {}
with open(META / "integrated/gene_type_map.tsv") as f:
    for row in csv.DictReader(f, delimiter='\t'):
        meta_gene_name_map[row['gene_id']] = row['gene_name']

# Use ChaRNA comparison for ChaRNA samples, totRNA for the other half
# all_genes_counts has all 12 samples — use ChaRNA comparison (both halves present)
meta_gene_counts, meta_gene_samples = load_count_matrix(
    META / "comparison_mito-ChaRNA_vs_ChaRNA/deseq2/all_genes_counts.tsv")
# totRNA samples only in their comparison
meta_gene_counts2, meta_gene_samples2 = load_count_matrix(
    META / "comparison_mito-totRNA_vs_totRNA/deseq2/all_genes_counts.tsv")
# Merge (both TSVs have all 12 but from different DESeq2 runs — use raw counts, identical)
# Actually all_genes_counts = raw featureCounts, same regardless of comparison
# Take union of samples from both files
meta_all_samples = list(dict.fromkeys(meta_gene_samples + meta_gene_samples2))
meta_gene_merged = {}
for s in meta_all_samples:
    if s in meta_gene_counts:
        meta_gene_merged[s] = meta_gene_counts[s]
    else:
        meta_gene_merged[s] = meta_gene_counts2[s]

meta_te_raw, _ = load_count_matrix(META / "integrated/te_counts_clean.tsv")
meta_sat_raw, _ = load_count_matrix(META / "integrated/satellite_counts_clean.tsv")

meta_cond_map = {
    'ChaRNA1': 'ChaRNA', 'ChaRNA2': 'ChaRNA', 'ChaRNA3': 'ChaRNA',
    'mito-ChaRNA1': 'mito-ChaRNA', 'mito-ChaRNA2': 'mito-ChaRNA', 'mito-ChaRNA3': 'mito-ChaRNA',
    'totRNA1': 'totRNA', 'totRNA2': 'totRNA', 'totRNA3': 'totRNA',
    'mito-totRNA1': 'mito-totRNA', 'mito-totRNA2': 'mito-totRNA', 'mito-totRNA3': 'mito-totRNA',
}
META_COND_ORDER = ['ChaRNA', 'mito-ChaRNA', 'totRNA', 'mito-totRNA']
META_LABELS = ['ChaRNA\n(interphase)', 'mito-ChaRNA\n(mitosis)',
               'totRNA\n(interphase)', 'mito-totRNA\n(mitosis)']

meta_agg_by_sample = {}
for sid in meta_all_samples:
    meta_agg_by_sample[sid] = aggregate_biotype(
        meta_gene_merged.get(sid, {}),
        meta_gene_map,
        meta_te_raw.get(sid, {}),
        meta_sat_raw.get(sid, {}),
        gene_name_map=meta_gene_name_map,
    )

meta_frac = compute_fractions(meta_agg_by_sample)
meta_stats = group_mean_sd(meta_frac, meta_cond_map)

# ─── Load De Cecco ───────────────────────────────────────────────────────────
print("[INF] Loading De Cecco data...")
decca_gene_map = load_gene_type_map(DECCA / "integrated/gene_type_map.tsv")
decca_gene_name_map = {}
with open(DECCA / "integrated/gene_type_map.tsv") as f:
    for row in csv.DictReader(f, delimiter='\t'):
        decca_gene_name_map[row['gene_id']] = row['gene_name']

decca_gene_counts, _ = load_count_matrix(
    DECCA / "comparison_SEN_L_vs_EP/deseq2/all_genes_counts.tsv")
# SEN_E samples from their comparison
decca_gene_counts2, _ = load_count_matrix(
    DECCA / "comparison_SEN_E_vs_EP/deseq2/all_genes_counts.tsv")
decca_gene_merged = {**decca_gene_counts, **decca_gene_counts2}

decca_te_raw, _ = load_count_matrix(DECCA / "integrated/te_counts_clean.tsv")
decca_sat_raw, _ = load_count_matrix(DECCA / "integrated/satellite_counts_clean.tsv")

decca_cond_map = {
    'EP1': 'EP', 'EP2': 'EP', 'EP3': 'EP',
    'SEN_E1': 'SEN_E', 'SEN_E2': 'SEN_E', 'SEN_E3': 'SEN_E',
    'SEN_L1': 'SEN_L', 'SEN_L2': 'SEN_L', 'SEN_L3': 'SEN_L',
}
DECCA_COND_ORDER = ['EP', 'SEN_E', 'SEN_L']
DECCA_LABELS = ['EP\n(proliferating)', 'SEN_E\n(early)', 'SEN_L\n(late)']

decca_agg_by_sample = {}
for sid in list(decca_gene_merged.keys()):
    decca_agg_by_sample[sid] = aggregate_biotype(
        decca_gene_merged.get(sid, {}),
        decca_gene_map,
        decca_te_raw.get(sid, {}),
        decca_sat_raw.get(sid, {}),
        gene_name_map=decca_gene_name_map,
    )

decca_frac = compute_fractions(decca_agg_by_sample)
decca_stats = group_mean_sd(decca_frac, decca_cond_map)

# ─── Print key values ────────────────────────────────────────────────────────
print("\n=== Key values (mean % ± SD) ===")
for label, stats, cond_order in [
    ("METACAR", meta_stats, META_COND_ORDER),
    ("De Cecco", decca_stats, DECCA_COND_ORDER),
]:
    print(f"\n{label}:")
    for cond in cond_order:
        d = stats.get(cond, {})
        te_m, te_s   = d.get('TE', (0,0))
        ln_m, ln_s   = d.get('lncRNA', (0,0))
        pc_m, pc_s   = d.get('protein_coding', (0,0))
        sat_m, sat_s = d.get('satellite', (0,0))
        print(f"  {cond:20s}  protein_coding={pc_m:.2f}±{pc_s:.2f}%  "
              f"TE={te_m:.2f}±{te_s:.2f}%  lncRNA={ln_m:.2f}±{ln_s:.2f}%  "
              f"satellite={sat_m:.3f}±{sat_s:.3f}%")

# ─── Plotting helpers ────────────────────────────────────────────────────────
def stacked_bars_with_errorbars(ax, stats, cond_order, labels, categories,
                                ylim, ylabel, title, bar_width=0.6):
    x = np.arange(len(cond_order))
    bottoms = np.zeros(len(cond_order))

    for cat in categories:
        means = np.array([stats.get(c, {}).get(cat, (0, 0))[0] for c in cond_order])
        sds   = np.array([stats.get(c, {}).get(cat, (0, 0))[1] for c in cond_order])
        bars  = ax.bar(x, means, bar_width, bottom=bottoms,
                       color=COLORS[cat], label=cat, zorder=2)
        # Error bars at top of each segment (cumulative)
        tops = bottoms + means
        ax.errorbar(x, tops, yerr=sds, fmt='none',
                    ecolor='#333333', elinewidth=0.8, capsize=3, zorder=3)
        bottoms += means

    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=9)
    ax.set_ylim(ylim)
    ax.set_ylabel(ylabel, fontsize=9)
    ax.set_title(title, fontsize=10, fontweight='bold', pad=6)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.tick_params(axis='both', labelsize=8)
    ax.yaxis.grid(True, linestyle='--', alpha=0.4, zorder=0)
    ax.set_axisbelow(True)

# ─── Figure ──────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(2, 2, figsize=(14, 11))
fig.patch.set_facecolor('white')

# Panel A — METACAR full
stacked_bars_with_errorbars(
    axes[0, 0], meta_stats, META_COND_ORDER, META_LABELS,
    CATEGORIES_FULL,
    ylim=(0, 102), ylabel='% of summed raw counts',
    title='A   ChaRNA — full biotype composition',
)

# Panel B — METACAR non-coding zoom
stacked_bars_with_errorbars(
    axes[0, 1], meta_stats, META_COND_ORDER, META_LABELS,
    CATEGORIES_NC,
    ylim=(0, 15), ylabel='% of summed raw counts',
    title='B   ChaRNA — non-protein-coding (expanded)',
)

# Panel C — De Cecco full
stacked_bars_with_errorbars(
    axes[1, 0], decca_stats, DECCA_COND_ORDER, DECCA_LABELS,
    CATEGORIES_FULL,
    ylim=(0, 102), ylabel='% of summed raw counts',
    title='C   Senescence — full biotype composition',
)

# Panel D — De Cecco non-coding zoom
stacked_bars_with_errorbars(
    axes[1, 1], decca_stats, DECCA_COND_ORDER, DECCA_LABELS,
    CATEGORIES_NC,
    ylim=(0, 12), ylabel='% of summed raw counts',
    title='D   Senescence — non-protein-coding (expanded)',
)

# ─── Shared legend ───────────────────────────────────────────────────────────
legend_handles = [
    mpatches.Patch(color=COLORS[c], label=c) for c in CATEGORIES_FULL
]
fig.legend(
    handles=legend_handles,
    loc='lower center',
    ncol=5,
    fontsize=9,
    frameon=True,
    framealpha=0.9,
    edgecolor='#CCCCCC',
    bbox_to_anchor=(0.5, 0.01),
)

fig.suptitle(
    'Supplementary Figure S1 — RNA biotype composition',
    fontsize=12, fontweight='bold', y=1.002
)
fig.tight_layout(rect=[0, 0.07, 1, 1])

# ─── Save ────────────────────────────────────────────────────────────────────
for ext, kw in [('png', {'dpi': 300}), ('pdf', {})]:
    out = OUTDIR / f"figS1_v2.{ext}"
    fig.savefig(out, bbox_inches='tight', facecolor='white', **kw)
    print(f"[OK] {out}")

plt.close(fig)
print("\n[DONE]")

TSV_NAME = "figS1_v2_fractions.tsv"

# ─── Fractions TSV (mean, SD per dataset/condition/class) ───────────────────
_tsv = OUTDIR / TSV_NAME
with open(_tsv, 'w') as _f:
    _f.write('dataset\tcondition\tclass\tmean_pct\tsd_pct\tn_samples\n')
    for _ds, _stats, _order, _cmap in [
        ('METACAR', meta_stats, META_COND_ORDER, meta_cond_map),
        ('DeCecco', decca_stats, DECCA_COND_ORDER, decca_cond_map),
    ]:
        for _cond in _order:
            _n = sum(1 for v in _cmap.values() if v == _cond)
            for _cat in CATEGORIES_FULL:
                _m, _s = _stats.get(_cond, {}).get(_cat, (0.0, 0.0))
                _f.write(f'{_ds}\t{_cond}\t{_cat}\t{_m!r}\t{_s!r}\t{_n}\n')
print(f"[OK] {_tsv}")
