#!/usr/bin/env python3
"""
GSE212630 Trigger Stage Analysis v1
====================================
TDP staging = outcome. Define independent "Trigger stage" (Pre-TDP axis).

Key premise: True trigger events may precede TDP-43 pathology.
"""

import pandas as pd
import numpy as np
from pathlib import Path
from scipy.stats import mannwhitneyu, spearmanr
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import gzip
import gc
import warnings
warnings.filterwarnings('ignore')

# ============================================================================
# CONFIGURATION
# ============================================================================
BASE_DIR = Path('/home/akaco/als/motor_cortex_analysis/ids_causal_analysis')
OUT_DIR = BASE_DIR / 'results_gse212630_triggerstage_v1'
OUT_DIR.mkdir(exist_ok=True)
(OUT_DIR / 'tables').mkdir(exist_ok=True)
(OUT_DIR / 'figures').mkdir(exist_ok=True)
(OUT_DIR / 'reports').mkdir(exist_ok=True)

GSE212630_DIR = BASE_DIR / 'GSE_212630_raw_expression_transposed'

# Gene sets (expanded)
GENE_SETS = {
    'DAMP': ['HMGB1', 'HSPA1A', 'HSPA1B', 'HSP90AA1', 'HSP90AB1',
             'S100A8', 'S100A9', 'S100A12', 'IL33', 'LGALS3', 'ANXA1', 'CALR', 'HSPD1'],
    'PRR': ['IFIH1', 'DDX58', 'DHX58', 'MAVS', 'EIF2AK2', 'TLR3', 'TLR7',
            'STING1', 'MB21D1', 'TREX1', 'SAMHD1', 'ADAR'],
    'ISG': ['ISG15', 'OAS1', 'OAS2', 'OAS3', 'MX1', 'MX2', 'IFIT1', 'IFIT2', 'IFIT3',
            'IFI44', 'IFI44L', 'RSAD2', 'IRF7', 'STAT1', 'STAT2', 'USP18', 'IFITM1', 'IFITM3', 'BST2'],
    'DDR': ['ATM', 'ATR', 'TP53BP1', 'H2AFX', 'BRCA1', 'BRCA2', 'PARP1', 'PARP2',
            'RAD51', 'XRCC1', 'XRCC5', 'XRCC6', 'NBN', 'MRE11', 'CHEK1', 'CHEK2', 'MDC1'],
    'RNAproc': ['XRN2', 'EXOSC2', 'DDX3X', 'DDX5', 'UPF1', 'SRSF1', 'SRSF2',
                'HNRNPA1', 'HNRNPK', 'PABPN1', 'FUS', 'TARDBP', 'TIA1'],
    'Mito': ['TFAM', 'POLG', 'NDUFA1', 'NDUFA2', 'SDHA', 'SDHB', 'UQCRB',
             'COX4I1', 'COX5A', 'ATP5F1A', 'ATP5F1B', 'VDAC1', 'SOD2', 'PINK1'],
    'Excitability': ['FOS', 'JUN', 'JUNB', 'EGR1', 'EGR2', 'ATF3', 'ARC', 'NPAS4'],
}

ALL_GENES = set()
for genes in GENE_SETS.values():
    ALL_GENES.update(genes)
ALL_GENES = list(ALL_GENES)

TDP_STAGES = ['Control', 'TDPneg', 'TDPmed', 'TDPhigh']
TDP_ORDER = {s: i for i, s in enumerate(TDP_STAGES)}

# Analysis mode thresholds
MODE_A_CELLS = 500
MODE_A_DONORS = 5
MODE_B_CELLS = 100

def read_genes_only(filepath, genes_of_interest):
    """Read only rows matching genes of interest from a gzipped CSV"""
    gene_set = set(genes_of_interest)
    data = {}
    header = None

    with gzip.open(filepath, 'rt') as f:
        for i, line in enumerate(f):
            parts = line.strip().split(',')
            if i == 0:
                header = parts
                continue
            gene = parts[0].strip('"')
            if gene in gene_set:
                data[gene] = [float(x) if x else 0.0 for x in parts[1:]]

    if not data:
        return None, None

    df = pd.DataFrame(data)
    cell_ids = [h.strip('"') for h in header[1:]]
    df.index = cell_ids
    return df, list(data.keys())

def z_normalize_within_celltype(df, score_col, celltype_col='cell_type'):
    """Z-normalize scores within each cell type"""
    result = pd.Series(index=df.index, dtype=float)
    for ct in df[celltype_col].unique():
        mask = df[celltype_col] == ct
        vals = df.loc[mask, score_col]
        if len(vals) > 1 and vals.std() > 0:
            result.loc[mask] = (vals - vals.mean()) / vals.std()
        else:
            result.loc[mask] = 0
    return result

def cliffs_delta(x, y):
    """Calculate Cliff's delta effect size"""
    n1, n2 = len(x), len(y)
    if n1 == 0 or n2 == 0:
        return np.nan
    count = 0
    for a in x:
        for b in y:
            if a > b:
                count += 1
            elif a < b:
                count -= 1
    return count / (n1 * n2)

# ============================================================================
# STEP 1: Load GSE212630 and calculate all scores
# ============================================================================
print("=" * 70)
print("STEP 1: Processing GSE212630 (all trigger axes)")
print("=" * 70)

gse_expr_files = sorted(GSE212630_DIR.glob('*_expression_transposed.csv.gz'))
print(f"Found {len(gse_expr_files)} cell type files")

all_cells = []
for i, expr_file in enumerate(gse_expr_files):
    ct_name = expr_file.stem.replace('_expression_transposed.csv', '').replace('.gz', '')
    meta_file = GSE212630_DIR / f"{ct_name}_metadata.csv"

    if not meta_file.exists():
        continue

    print(f"  [{i+1}/{len(gse_expr_files)}] {ct_name}...", end=' ', flush=True)

    meta = pd.read_csv(meta_file)
    expr, found_genes = read_genes_only(expr_file, ALL_GENES)
    if expr is None:
        print("no genes!")
        continue

    # Strip suffix from cell IDs
    expr.index = expr.index.str.replace(r'_\d+_(Control|TDPneg|TDPmed|TDPhigh)$', '', regex=True)

    common_cells = [c for c in meta['cell_id'] if c in expr.index]
    if len(common_cells) == 0:
        print("no cells!")
        continue

    expr_aligned = expr.loc[common_cells]
    meta_aligned = meta[meta['cell_id'].isin(common_cells)].set_index('cell_id').loc[common_cells].reset_index()

    # Calculate scores (both mean and detection)
    scores = pd.DataFrame({'cell_id': common_cells})
    for set_name, genes in GENE_SETS.items():
        gene_cols = [g for g in genes if g in expr_aligned.columns]
        if len(gene_cols) > 0:
            scores[f'{set_name}_det'] = (expr_aligned[gene_cols] > 0).sum(axis=1).values / len(gene_cols)
            scores[f'{set_name}_mean'] = expr_aligned[gene_cols].mean(axis=1).values

    # NSA = PRR - ISG
    if 'PRR_det' in scores.columns and 'ISG_det' in scores.columns:
        scores['NSA_det'] = scores['PRR_det'] - scores['ISG_det']
        scores['NSA_mean'] = scores['PRR_mean'] - scores['ISG_mean']

    # QC proxy: number of detected genes (for DAMP bias check)
    scores['n_genes_detected'] = (expr_aligned > 0).sum(axis=1).values

    scores['sample'] = meta_aligned['sample'].values
    scores['condition'] = meta_aligned['condition'].values
    scores['cell_type'] = meta_aligned['cell_type'].values

    all_cells.append(scores)
    print(f"{len(common_cells)} cells, {len(found_genes)} genes")

    del expr, expr_aligned
    gc.collect()

df = pd.concat(all_cells, ignore_index=True)
print(f"\nTotal: {len(df):,} cells")
print(f"TDP stages: {df['condition'].value_counts().to_dict()}")

# Z-normalize within cell type
score_cols = ['DAMP_det', 'PRR_det', 'ISG_det', 'DDR_det', 'RNAproc_det', 'Mito_det', 'NSA_det', 'Excitability_det']
for col in score_cols:
    if col in df.columns:
        df[f'{col}_z'] = z_normalize_within_celltype(df, col)

# ============================================================================
# STEP 2: Create celltype_counts_by_stage.csv (A)
# ============================================================================
print("\n" + "=" * 70)
print("STEP 2: Creating celltype_counts_by_stage.csv")
print("=" * 70)

counts = df.groupby(['cell_type', 'condition']).agg(
    n_cells=('cell_id', 'count'),
    n_donors=('sample', 'nunique')
).reset_index()

# Pivot for easier viewing
counts_pivot = counts.pivot(index='cell_type', columns='condition', values='n_cells').fillna(0).astype(int)
counts_pivot = counts_pivot[TDP_STAGES]
counts_pivot['total'] = counts_pivot.sum(axis=1)

# Add donor counts
donors_pivot = counts.pivot(index='cell_type', columns='condition', values='n_donors').fillna(0).astype(int)
donors_pivot = donors_pivot[TDP_STAGES]
donors_pivot.columns = [f'{c}_donors' for c in TDP_STAGES]
counts_pivot = counts_pivot.join(donors_pivot)

# Determine analysis mode
def get_mode(row):
    total_cells = row['total']
    min_donors = min(row[[f'{s}_donors' for s in TDP_STAGES]].values)
    if total_cells >= MODE_A_CELLS and min_donors >= MODE_A_DONORS:
        return 'A'
    elif total_cells >= MODE_B_CELLS:
        return 'B'
    else:
        return 'C'

counts_pivot['analysis_mode'] = counts_pivot.apply(get_mode, axis=1)
counts_pivot = counts_pivot.sort_values('total', ascending=False)
counts_pivot.to_csv(OUT_DIR / 'tables/celltype_counts_by_stage.csv')

print(f"Saved: celltype_counts_by_stage.csv")
print(f"\nAnalysis modes:")
print(f"  Mode A (main): {(counts_pivot['analysis_mode'] == 'A').sum()} cell types")
print(f"  Mode B (secondary): {(counts_pivot['analysis_mode'] == 'B').sum()} cell types")
print(f"  Mode C (exploratory): {(counts_pivot['analysis_mode'] == 'C').sum()} cell types")

# ============================================================================
# STEP 3: Define Trigger Gate and calculate fractions (B)
# ============================================================================
print("\n" + "=" * 70)
print("STEP 3: Defining Trigger Gate and fractions")
print("=" * 70)

def define_trigger_gate(data, nsa_q=0.75, ddr_q=0.50, damp_q=0.50):
    """
    Trigger Gate: |NSA| >= q, DDR >= q, DAMP <= q (within cell type)
    """
    data = data.copy()
    data['trigger_gate'] = False

    for ct in data['cell_type'].unique():
        mask = data['cell_type'] == ct
        ct_df = data.loc[mask]

        if len(ct_df) < 50:
            continue

        nsa_thresh = ct_df['NSA_det_z'].abs().quantile(nsa_q)
        ddr_thresh = ct_df['DDR_det_z'].quantile(ddr_q)
        damp_thresh = ct_df['DAMP_det_z'].quantile(damp_q)

        gate = (
            (ct_df['NSA_det_z'].abs() >= nsa_thresh) &
            (ct_df['DDR_det_z'] >= ddr_thresh) &
            (ct_df['DAMP_det_z'] <= damp_thresh)
        )
        data.loc[mask, 'trigger_gate'] = gate.values

    return data

# Apply default gate
df = define_trigger_gate(df)
print(f"Trigger Gate cells: {df['trigger_gate'].sum():,} ({100*df['trigger_gate'].mean():.1f}%)")

# Calculate fraction by celltype x stage
gate_frac = df.groupby(['cell_type', 'condition']).agg(
    gate_count=('trigger_gate', 'sum'),
    total_count=('trigger_gate', 'count'),
    gate_fraction=('trigger_gate', 'mean')
).reset_index()

# Pivot
gate_pivot = gate_frac.pivot(index='cell_type', columns='condition', values='gate_fraction')
gate_pivot = gate_pivot[TDP_STAGES]
gate_pivot.columns = [f'{c}_frac' for c in TDP_STAGES]

# Add counts
count_pivot = gate_frac.pivot(index='cell_type', columns='condition', values='total_count').fillna(0).astype(int)
count_pivot = count_pivot[TDP_STAGES]
count_pivot.columns = [f'{c}_n' for c in TDP_STAGES]
gate_pivot = gate_pivot.join(count_pivot)

# Add analysis mode
gate_pivot['analysis_mode'] = counts_pivot['analysis_mode']

# Calculate Control -> TDPneg change
gate_pivot['ctrl_to_neg_diff'] = gate_pivot['TDPneg_frac'] - gate_pivot['Control_frac']
gate_pivot['ctrl_to_neg_pct'] = 100 * (gate_pivot['TDPneg_frac'] - gate_pivot['Control_frac']) / gate_pivot['Control_frac'].replace(0, np.nan)

gate_pivot = gate_pivot.sort_values('ctrl_to_neg_diff', ascending=False)
gate_pivot.to_csv(OUT_DIR / 'tables/trigger_gate_fraction_by_celltype_stage.csv')
print(f"Saved: trigger_gate_fraction_by_celltype_stage.csv")

# ============================================================================
# STEP 4: Create Pre-TDP PCA map (C)
# ============================================================================
print("\n" + "=" * 70)
print("STEP 4: Creating Pre-TDP PCA map (celltype signature)")
print("=" * 70)

# Aggregate by celltype x condition
feature_cols = [c for c in df.columns if c.endswith('_z') and not c.startswith('n_')]
ct_cond_agg = df.groupby(['cell_type', 'condition'])[feature_cols].mean().reset_index()

# Create feature matrix for PCA
X_ct = ct_cond_agg[feature_cols].values
scaler = StandardScaler()
X_ct_scaled = scaler.fit_transform(X_ct)

# PCA
pca = PCA(n_components=min(5, len(feature_cols)))
pca_result = pca.fit_transform(X_ct_scaled)
ct_cond_agg['PC1'] = pca_result[:, 0]
ct_cond_agg['PC2'] = pca_result[:, 1]

# Plot
fig, ax = plt.subplots(figsize=(10, 8))
colors = {'Control': 'blue', 'TDPneg': 'green', 'TDPmed': 'orange', 'TDPhigh': 'red'}
for stage in TDP_STAGES:
    stage_data = ct_cond_agg[ct_cond_agg['condition'] == stage]
    ax.scatter(stage_data['PC1'], stage_data['PC2'], c=colors[stage],
               label=stage, alpha=0.7, s=100)
    # Add cell type labels for key types
    for _, row in stage_data.iterrows():
        if row['cell_type'] in ['Glia.Astro.GFAP.pos', 'Glia.OPC', 'Glia.Micro', 'Glia.Oligo']:
            ax.annotate(row['cell_type'].split('.')[-1], (row['PC1'], row['PC2']),
                       fontsize=7, alpha=0.7)

ax.set_xlabel(f'PC1 ({100*pca.explained_variance_ratio_[0]:.1f}%)')
ax.set_ylabel(f'PC2 ({100*pca.explained_variance_ratio_[1]:.1f}%)')
ax.set_title('Pre-TDP Map: Celltype × TDP Stage (PCA)')
ax.legend(title='TDP Stage')
ax.axhline(0, color='gray', linestyle='--', alpha=0.3)
ax.axvline(0, color='gray', linestyle='--', alpha=0.3)
plt.tight_layout()
plt.savefig(OUT_DIR / 'figures/preTDP_map_pca_celltype.png', dpi=150)
plt.close()
print(f"Saved: preTDP_map_pca_celltype.png")
print(f"  PC1: {100*pca.explained_variance_ratio_[0]:.1f}%, PC2: {100*pca.explained_variance_ratio_[1]:.1f}%")

# Save PCA loadings
loadings = pd.DataFrame(pca.components_.T, index=feature_cols,
                        columns=[f'PC{i+1}' for i in range(pca.n_components_)])
loadings.to_csv(OUT_DIR / 'tables/pca_loadings_celltype.csv')

# ============================================================================
# STEP 5: Donor-level tests for Trigger Gate (Mode A/B)
# ============================================================================
print("\n" + "=" * 70)
print("STEP 5: Donor-level statistical tests")
print("=" * 70)

# Aggregate by donor x celltype
donor_agg = df.groupby(['sample', 'condition', 'cell_type']).agg(
    trigger_gate_frac=('trigger_gate', 'mean'),
    DAMP_mean=('DAMP_det', 'mean'),
    NSA_mean=('NSA_det', 'mean'),
    DDR_mean=('DDR_det', 'mean'),
    n_cells=('cell_id', 'count'),
    n_genes_mean=('n_genes_detected', 'mean')
).reset_index()

donor_agg.to_csv(OUT_DIR / 'tables/donor_level_scores.csv', index=False)
print(f"Saved: donor_level_scores.csv")

# Statistical tests
print("\nTrigger Gate: Control vs TDPneg (donor-level)")
print("-" * 60)

tests = []
mode_a_types = counts_pivot[counts_pivot['analysis_mode'] == 'A'].index.tolist()
mode_b_types = counts_pivot[counts_pivot['analysis_mode'] == 'B'].index.tolist()

for ct in donor_agg['cell_type'].unique():
    ct_data = donor_agg[donor_agg['cell_type'] == ct]
    mode = 'A' if ct in mode_a_types else ('B' if ct in mode_b_types else 'C')

    ctrl = ct_data[ct_data['condition'] == 'Control']['trigger_gate_frac'].dropna()
    tdpneg = ct_data[ct_data['condition'] == 'TDPneg']['trigger_gate_frac'].dropna()
    tdpmed = ct_data[ct_data['condition'] == 'TDPmed']['trigger_gate_frac'].dropna()
    tdphigh = ct_data[ct_data['condition'] == 'TDPhigh']['trigger_gate_frac'].dropna()

    if len(ctrl) < 2 or len(tdpneg) < 2:
        continue

    # Mann-Whitney U test
    try:
        stat, pval = mannwhitneyu(ctrl, tdpneg, alternative='two-sided')
    except:
        pval = np.nan

    # Effect size (Cliff's delta)
    delta = cliffs_delta(tdpneg.values, ctrl.values)

    # Trend test
    all_vals, all_stages = [], []
    for stage_idx, (stage, vals) in enumerate([('Control', ctrl), ('TDPneg', tdpneg),
                                                ('TDPmed', tdpmed), ('TDPhigh', tdphigh)]):
        all_vals.extend(vals.tolist())
        all_stages.extend([stage_idx] * len(vals))

    if len(all_vals) > 5:
        rho, trend_p = spearmanr(all_stages, all_vals)
    else:
        rho, trend_p = np.nan, np.nan

    tests.append({
        'cell_type': ct,
        'analysis_mode': mode,
        'ctrl_mean': ctrl.mean(),
        'tdpneg_mean': tdpneg.mean(),
        'tdpmed_mean': tdpmed.mean() if len(tdpmed) > 0 else np.nan,
        'tdphigh_mean': tdphigh.mean() if len(tdphigh) > 0 else np.nan,
        'ctrl_vs_neg_pval': pval,
        'cliffs_delta': delta,
        'trend_rho': rho,
        'trend_pval': trend_p,
        'n_ctrl_donors': len(ctrl),
        'n_neg_donors': len(tdpneg)
    })

    sig = '*' if pval < 0.05 else ''
    if mode in ['A', 'B'] and (pval < 0.1 or abs(delta) > 0.3):
        pct = 100 * (tdpneg.mean() - ctrl.mean()) / ctrl.mean() if ctrl.mean() > 0 else np.nan
        print(f"  [{mode}] {ct:30s}: {pct:+6.1f}% (p={pval:.4f}{sig}, δ={delta:+.2f})")

tests_df = pd.DataFrame(tests)
tests_df = tests_df.sort_values('ctrl_vs_neg_pval')
tests_df.to_csv(OUT_DIR / 'tables/trend_tests_control_tdpneg_med_high.csv', index=False)

# ============================================================================
# STEP 6: Threshold sweep stability
# ============================================================================
print("\n" + "=" * 70)
print("STEP 6: Threshold sweep stability analysis")
print("=" * 70)

nsa_qs = [0.70, 0.75, 0.80, 0.85]
ddr_qs = [0.40, 0.50, 0.60]
damp_qs = [0.40, 0.50, 0.60]

sweep_results = []
reference_ranking = None

for nsa_q in nsa_qs:
    for ddr_q in ddr_qs:
        for damp_q in damp_qs:
            df_sweep = define_trigger_gate(df.copy(), nsa_q=nsa_q, ddr_q=ddr_q, damp_q=damp_q)

            # Get ranking by Control->TDPneg increase
            sweep_frac = df_sweep.groupby(['cell_type', 'condition'])['trigger_gate'].mean().reset_index()
            sweep_pivot = sweep_frac.pivot(index='cell_type', columns='condition', values='trigger_gate')
            if 'Control' in sweep_pivot.columns and 'TDPneg' in sweep_pivot.columns:
                sweep_pivot['diff'] = sweep_pivot['TDPneg'] - sweep_pivot['Control']
                ranking = sweep_pivot.sort_values('diff', ascending=False).index.tolist()

                if reference_ranking is None:
                    reference_ranking = ranking

                # Calculate top-5 overlap with reference
                top5_overlap = len(set(ranking[:5]) & set(reference_ranking[:5])) / 5

                sweep_results.append({
                    'nsa_q': nsa_q,
                    'ddr_q': ddr_q,
                    'damp_q': damp_q,
                    'gate_pct': 100 * df_sweep['trigger_gate'].mean(),
                    'top5_stability': top5_overlap,
                    'top5_celltypes': ', '.join(ranking[:5])
                })

sweep_df = pd.DataFrame(sweep_results)
sweep_df.to_csv(OUT_DIR / 'tables/trigger_gate_threshold_sweep_stability.csv', index=False)
print(f"Saved: trigger_gate_threshold_sweep_stability.csv")
print(f"  Mean top-5 stability: {sweep_df['top5_stability'].mean():.2f}")

# ============================================================================
# STEP 7: DAMP vs QC bias check
# ============================================================================
print("\n" + "=" * 70)
print("STEP 7: DAMP vs QC bias check")
print("=" * 70)

# Correlation at donor x celltype level
qc_check = donor_agg.groupby('cell_type').apply(
    lambda x: pd.Series({
        'damp_ngenes_corr': x['DAMP_mean'].corr(x['n_genes_mean']),
        'n_obs': len(x)
    })
).reset_index()
qc_check.to_csv(OUT_DIR / 'tables/qc_bias_check.csv', index=False)
print(f"Saved: qc_bias_check.csv")

# Plot DAMP vs n_genes by stage
fig, axes = plt.subplots(1, 4, figsize=(16, 4))
for i, stage in enumerate(TDP_STAGES):
    ax = axes[i]
    stage_data = donor_agg[donor_agg['condition'] == stage]
    ax.scatter(stage_data['n_genes_mean'], stage_data['DAMP_mean'], alpha=0.5, s=30)
    r = stage_data['n_genes_mean'].corr(stage_data['DAMP_mean'])
    ax.set_xlabel('Mean n_genes detected')
    ax.set_ylabel('DAMP (mean)')
    ax.set_title(f'{stage} (r={r:.2f})')
plt.tight_layout()
plt.savefig(OUT_DIR / 'figures/damp_vs_qc_scatter.png', dpi=150)
plt.close()
print(f"Saved: damp_vs_qc_scatter.png")

# ============================================================================
# STEP 8: Additional visualizations
# ============================================================================
print("\n" + "=" * 70)
print("STEP 8: Additional visualizations")
print("=" * 70)

# Trigger gate heatmap
fig, ax = plt.subplots(figsize=(8, 10))
gate_heatmap = gate_pivot[[f'{s}_frac' for s in TDP_STAGES]].copy()
gate_heatmap.columns = TDP_STAGES
gate_heatmap = gate_heatmap.sort_values('TDPneg', ascending=False)
sns.heatmap(gate_heatmap, cmap='YlOrRd', ax=ax, annot=True, fmt='.2f', cbar_kws={'label': 'Trigger Gate Fraction'})
ax.set_title('Trigger Gate Fraction by Cell Type × TDP Stage')
plt.tight_layout()
plt.savefig(OUT_DIR / 'figures/trigger_gate_heatmap.png', dpi=150)
plt.close()
print(f"Saved: trigger_gate_heatmap.png")

# Trigger gate stage bars (top cell types)
fig, ax = plt.subplots(figsize=(12, 6))
top_cts = gate_pivot.head(10).index.tolist()
plot_data = gate_frac[gate_frac['cell_type'].isin(top_cts)]
plot_pivot = plot_data.pivot(index='cell_type', columns='condition', values='gate_fraction')
if TDP_STAGES[0] in plot_pivot.columns:
    plot_pivot = plot_pivot[TDP_STAGES]
plot_pivot.plot(kind='bar', ax=ax, width=0.8, color=[colors[s] for s in TDP_STAGES])
ax.set_ylabel('Trigger Gate Fraction')
ax.set_title('Trigger Gate by TDP Stage (Top Cell Types)')
ax.legend(title='TDP Stage')
plt.xticks(rotation=45, ha='right')
plt.tight_layout()
plt.savefig(OUT_DIR / 'figures/trigger_gate_stage_bars.png', dpi=150)
plt.close()
print(f"Saved: trigger_gate_stage_bars.png")

# Donor PCA
donor_features = df.groupby(['sample', 'condition'])[feature_cols].mean().reset_index()
X_donor = donor_features[feature_cols].values
X_donor_scaled = scaler.fit_transform(X_donor)
pca_donor = PCA(n_components=min(5, len(feature_cols)))
pca_donor_result = pca_donor.fit_transform(X_donor_scaled)
donor_features['PC1'] = pca_donor_result[:, 0]
donor_features['PC2'] = pca_donor_result[:, 1]

fig, ax = plt.subplots(figsize=(10, 8))
for stage in TDP_STAGES:
    stage_data = donor_features[donor_features['condition'] == stage]
    ax.scatter(stage_data['PC1'], stage_data['PC2'], c=colors[stage],
               label=f'{stage} (n={len(stage_data)})', alpha=0.7, s=100)
ax.set_xlabel(f'PC1 ({100*pca_donor.explained_variance_ratio_[0]:.1f}%)')
ax.set_ylabel(f'PC2 ({100*pca_donor.explained_variance_ratio_[1]:.1f}%)')
ax.set_title('Pre-TDP Map: Donor Signature (PCA)')
ax.legend(title='TDP Stage')
ax.axhline(0, color='gray', linestyle='--', alpha=0.3)
ax.axvline(0, color='gray', linestyle='--', alpha=0.3)
plt.tight_layout()
plt.savefig(OUT_DIR / 'figures/preTDP_map_pca_donor.png', dpi=150)
plt.close()
print(f"Saved: preTDP_map_pca_donor.png")

# ============================================================================
# STEP 9: Generate Report
# ============================================================================
print("\n" + "=" * 70)
print("STEP 9: Generating report")
print("=" * 70)

# Get top findings
mode_a_tests = tests_df[tests_df['analysis_mode'] == 'A'].head(5)
mode_b_tests = tests_df[tests_df['analysis_mode'] == 'B'].head(5)

# QC bias summary
mean_damp_corr = qc_check['damp_ngenes_corr'].mean()

report = f"""# GSE212630 Trigger Stage Analysis Report

## Key Premise
- TDP staging (Control/TDPneg/TDPmed/TDPhigh) = **OUTCOME axis**
- True trigger events may exist **before TDP-43 pathology**
- Goal: Identify cell types where "Trigger Gate" enriches in pre-pathology (TDPneg)

## Data Summary
- Total cells: {len(df):,}
- TDP stages: Control {(df['condition']=='Control').sum():,}, TDPneg {(df['condition']=='TDPneg').sum():,}, TDPmed {(df['condition']=='TDPmed').sum():,}, TDPhigh {(df['condition']=='TDPhigh').sum():,}
- Cell types: {df['cell_type'].nunique()}

## Analysis Mode Distribution
- Mode A (main, n≥500 cells, ≥5 donors): {(counts_pivot['analysis_mode']=='A').sum()} cell types
- Mode B (secondary, n≥100 cells): {(counts_pivot['analysis_mode']=='B').sum()} cell types
- Mode C (exploratory, n<100): {(counts_pivot['analysis_mode']=='C').sum()} cell types

## Trigger Gate Definition
- |NSA| ≥ q75 (nucleic acid sensing imbalance)
- DDR ≥ q50 (DNA damage response active)
- DAMP ≤ q50 (not yet in death/damage state)
- Total: {df['trigger_gate'].sum():,} cells ({100*df['trigger_gate'].mean():.1f}%)

## Key Findings

### Top Cell Types (Control → TDPneg, Mode A)
"""

for _, row in mode_a_tests.iterrows():
    sig = '*' if row['ctrl_vs_neg_pval'] < 0.05 else ''
    pct = 100 * (row['tdpneg_mean'] - row['ctrl_mean']) / row['ctrl_mean'] if row['ctrl_mean'] > 0 else np.nan
    report += f"\n- **{row['cell_type']}**: {pct:+.1f}% (p={row['ctrl_vs_neg_pval']:.4f}{sig}, δ={row['cliffs_delta']:+.2f})"

report += f"""

### Top Cell Types (Control → TDPneg, Mode B)
"""

for _, row in mode_b_tests.iterrows():
    sig = '*' if row['ctrl_vs_neg_pval'] < 0.05 else ''
    pct = 100 * (row['tdpneg_mean'] - row['ctrl_mean']) / row['ctrl_mean'] if row['ctrl_mean'] > 0 else np.nan
    report += f"\n- **{row['cell_type']}**: {pct:+.1f}% (p={row['ctrl_vs_neg_pval']:.4f}{sig}, δ={row['cliffs_delta']:+.2f})"

report += f"""

### Threshold Stability
- Mean top-5 stability across threshold sweep: {sweep_df['top5_stability'].mean():.2f}
- Most stable top cell types: {sweep_df.iloc[0]['top5_celltypes']}

### DAMP vs QC Bias Check
- Mean DAMP-nGenes correlation across cell types: r = {mean_damp_corr:.2f}
- Interpretation: {"Potential QC bias (positive corr suggests DAMP tracks RNA quality)" if mean_damp_corr > 0.3 else "Low QC bias (DAMP not strongly correlated with RNA quality)"}

## PCA Summary
- Celltype × Stage PCA: PC1 {100*pca.explained_variance_ratio_[0]:.1f}%, PC2 {100*pca.explained_variance_ratio_[1]:.1f}%
- Donor PCA: PC1 {100*pca_donor.explained_variance_ratio_[0]:.1f}%, PC2 {100*pca_donor.explained_variance_ratio_[1]:.1f}%

## Caveats
1. TDP staging is derived from pathological burden, not temporal sequence
2. "Trigger Gate" is a hypothesis-driven definition, not ground truth
3. Small cell counts in Mode B/C limit statistical power
4. Cross-validate with Phase13 (ALS/Control) recommended

## Output Files

### Tables
- `celltype_counts_by_stage.csv`
- `trigger_gate_fraction_by_celltype_stage.csv`
- `donor_level_scores.csv`
- `trend_tests_control_tdpneg_med_high.csv`
- `trigger_gate_threshold_sweep_stability.csv`
- `qc_bias_check.csv`
- `pca_loadings_celltype.csv`

### Figures
- `preTDP_map_pca_celltype.png`
- `preTDP_map_pca_donor.png`
- `trigger_gate_heatmap.png`
- `trigger_gate_stage_bars.png`
- `damp_vs_qc_scatter.png`
"""

with open(OUT_DIR / 'reports/trigger_stage_summary.md', 'w') as f:
    f.write(report)
print(f"Saved: trigger_stage_summary.md")

print("\n" + "=" * 70)
print("COMPLETE")
print("=" * 70)
print(f"Output: {OUT_DIR}")
