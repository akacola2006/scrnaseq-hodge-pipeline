#!/usr/bin/env python3
"""
Phase PTv2 - Step 3: Compare 3 PTs and Test Robustness
=======================================================

Purpose:
- Compare PT_imes, PT_dpt, PT_stress across cell types and subtypes
- Test if key findings are robust across different PT definitions:
  1. VAT1L timing differs by subtype
  2. Microglia inverse correlation
  3. Cell type ordering (Oligo, Upper, VAT1L, Micro)

Input:
- cell_level_features_ALL_with_PT_dpt_stress.csv
- patient_clusters.csv (from Phase 5′)

Output:
- celltype_PT_summary.csv
- VAT1L_PT_by_subtype_comparison.csv
- PT_robustness_analysis.csv
- PT_comparison_plots.png
"""

import pandas as pd
import numpy as np
from scipy import stats
from scipy.stats import kendalltau
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

print("="*80)
print("Phase PTv2 - Step 3: Comparing 3 PTs and Testing Robustness")
print("="*80)

# Paths
base_dir = Path('/home/akaco/als/motor_cortex_analysis')
ids_dir = base_dir / 'ids_causal_analysis'
input_file = ids_dir / 'results' / 'PTv2_robustness' / 'cell_level_features_ALL_with_PT_dpt_stress.csv'
clusters_file = ids_dir / 'results' / 'patient_stratified' / 'patient_clusters.csv'
output_dir = ids_dir / 'results' / 'PTv2_robustness'

# ============================================================================
# Load Data
# ============================================================================
print("\n[Loading data...]")

df = pd.read_csv(input_file)
print(f"  Total cells: {len(df):,}")

# Load patient clusters
clusters = pd.read_csv(clusters_file)
print(f"  Patient clusters: {len(clusters)}")

# Define subtype mapping
subtype_map = {
    0: 'Upper-layer',
    1: 'Pure Oligo',
    2: 'Oligo-Inflammation',
    3: 'Upper-layer',
    4: 'Pure Oligo',
    5: 'Oligo-Inflammation'
}
clusters['subtype'] = clusters['cluster_id'].map(subtype_map)

# Merge with cell data
df = df.merge(clusters[['patient_id', 'cluster_id', 'subtype']], on='patient_id', how='left')

# Filter ALS only
df_als = df[df['condition'] == 'ALS'].copy()
print(f"  ALS cells: {len(df_als):,}")
print(f"  Subtypes: {df_als['subtype'].value_counts().to_dict()}")

# ============================================================================
# Step 3-1: Cell-Type Level PT Summary
# ============================================================================
print("\n[Step 3-1] Calculating cell-type-level PT summaries...")

PT_versions = ['PT_imes', 'PT_dpt', 'PT_stress']

celltype_summary_list = []

for celltype in df_als['cell_type'].unique():
    celltype_cells = df_als[df_als['cell_type'] == celltype]

    if len(celltype_cells) < 5:  # Skip rare cell types
        continue

    summary = {
        'cell_type': celltype,
        'n_cells': len(celltype_cells)
    }

    for pt_version in PT_versions:
        summary[f'{pt_version}_median'] = celltype_cells[pt_version].median()
        summary[f'{pt_version}_mean'] = celltype_cells[pt_version].mean()
        summary[f'{pt_version}_std'] = celltype_cells[pt_version].std()

    celltype_summary_list.append(summary)

celltype_summary = pd.DataFrame(celltype_summary_list)
celltype_summary = celltype_summary.sort_values('PT_imes_median')

output_file = output_dir / 'celltype_PT_summary.csv'
celltype_summary.to_csv(output_file, index=False)
print(f"  Saved: {output_file}")
print(f"  Cell types analyzed: {len(celltype_summary)}")

# ============================================================================
# Step 3-2: VAT1L PT by Subtype (for all 3 PTs)
# ============================================================================
print("\n[Step 3-2] Analyzing VAT1L PT by subtype...")

vat1l_cells = df_als[df_als['cell_type'].str.contains('VAT1L', na=False)]
print(f"  VAT1L cells: {len(vat1l_cells)}")

vat1l_summary_list = []

for pt_version in PT_versions:
    print(f"\n  Analyzing {pt_version}...")

    for subtype in ['Pure Oligo', 'Oligo-Inflammation', 'Upper-layer']:
        subtype_cells = vat1l_cells[vat1l_cells['subtype'] == subtype]

        if len(subtype_cells) == 0:
            continue

        summary = {
            'PT_version': pt_version,
            'subtype': subtype,
            'n_cells': len(subtype_cells),
            'n_patients': subtype_cells['patient_id'].nunique(),
            'PT_median': subtype_cells[pt_version].median(),
            'PT_mean': subtype_cells[pt_version].mean(),
            'PT_std': subtype_cells[pt_version].std(),
            'PT_min': subtype_cells[pt_version].min(),
            'PT_max': subtype_cells[pt_version].max()
        }

        vat1l_summary_list.append(summary)
        print(f"    {subtype}: PT={summary['PT_mean']:.3f} ± {summary['PT_std']:.3f} (n={summary['n_cells']})")

vat1l_summary = pd.DataFrame(vat1l_summary_list)

output_file = output_dir / 'VAT1L_PT_by_subtype_comparison.csv'
vat1l_summary.to_csv(output_file, index=False)
print(f"\n  Saved: {output_file}")

# Statistical tests for each PT version
print("\n  Statistical tests (VAT1L PT across subtypes):")

for pt_version in PT_versions:
    print(f"\n    {pt_version}:")

    oligo_pure = vat1l_cells[vat1l_cells['subtype'] == 'Pure Oligo'][pt_version]
    oligo_inflam = vat1l_cells[vat1l_cells['subtype'] == 'Oligo-Inflammation'][pt_version]
    upper = vat1l_cells[vat1l_cells['subtype'] == 'Upper-layer'][pt_version]

    if len(oligo_pure) > 0 and len(oligo_inflam) > 0:
        t, p = stats.ttest_ind(oligo_pure, oligo_inflam)
        print(f"      Pure Oligo vs Oligo-Inflammation: p={p:.4f}")

    if len(oligo_pure) > 0 and len(upper) > 0:
        t, p = stats.ttest_ind(oligo_pure, upper)
        print(f"      Pure Oligo vs Upper-layer: p={p:.4f}")

    if len(oligo_inflam) > 0 and len(upper) > 0:
        t, p = stats.ttest_ind(oligo_inflam, upper)
        print(f"      Oligo-Inflammation vs Upper-layer: p={p:.4f}")

# ============================================================================
# Step 3-3: Microglia PT by Subtype (test inverse correlation)
# ============================================================================
print("\n[Step 3-3] Analyzing Microglia PT by subtype...")

micro_cells = df_als[df_als['cell_type'].str.contains('Glia.Micro', na=False)]
print(f"  Microglia cells: {len(micro_cells)}")

micro_summary_list = []

for pt_version in PT_versions:
    for subtype in ['Pure Oligo', 'Oligo-Inflammation', 'Upper-layer']:
        subtype_cells = micro_cells[micro_cells['subtype'] == subtype]

        if len(subtype_cells) == 0:
            continue

        summary = {
            'PT_version': pt_version,
            'subtype': subtype,
            'n_cells': len(subtype_cells),
            'PT_median': subtype_cells[pt_version].median(),
            'PT_mean': subtype_cells[pt_version].mean(),
            'PT_std': subtype_cells[pt_version].std()
        }

        micro_summary_list.append(summary)

micro_summary = pd.DataFrame(micro_summary_list)

output_file = output_dir / 'Microglia_PT_by_subtype_comparison.csv'
micro_summary.to_csv(output_file, index=False)
print(f"  Saved: {output_file}")

# ============================================================================
# Step 3-4: Kendall's Tau for Cell Type Ordering
# ============================================================================
print("\n[Step 3-4] Calculating Kendall's tau for cell type ordering...")

# Key cell types for ordering comparison
key_celltypes = [
    'Ex.L5.VAT1L.EYA4',
    'Ex.L5.VAT1L.THSD4',
    'Glia.Oligo',
    'Glia.Astro.GFAP-neg',
    'Ex.L2.L3.CUX2.RASGRF2',
    'Ex.L3.L5.CUX2.RORB',
    'Glia.Micro',
    'In.PV.PVALB.CEMIP',
    'Ex.L4.L6.RORB.LRRK1',
    'Ex.L5.L6.THEMIS.TMEM233'
]

# Filter to available cell types
key_celltypes = [ct for ct in key_celltypes if ct in celltype_summary['cell_type'].values]
print(f"  Key cell types: {len(key_celltypes)}")

# Get rankings for each PT version
rankings = {}

for pt_version in PT_versions:
    # Sort cell types by median PT
    sorted_celltypes = celltype_summary.sort_values(f'{pt_version}_median')

    # Get rank of key cell types
    ranking = []
    for ct in key_celltypes:
        if ct in sorted_celltypes['cell_type'].values:
            rank = sorted_celltypes[sorted_celltypes['cell_type'] == ct].index[0]
            ranking.append(rank)
        else:
            ranking.append(np.nan)

    rankings[pt_version] = ranking

# Calculate pairwise Kendall's tau
kendall_results = []

for i, pt1 in enumerate(PT_versions):
    for j, pt2 in enumerate(PT_versions):
        if i < j:  # Only upper triangle
            # Remove NaN pairs
            mask = ~(np.isnan(rankings[pt1]) | np.isnan(rankings[pt2]))
            rank1 = np.array(rankings[pt1])[mask]
            rank2 = np.array(rankings[pt2])[mask]

            tau, p = kendalltau(rank1, rank2)
            kendall_results.append({
                'PT1': pt1,
                'PT2': pt2,
                'kendall_tau': tau,
                'p_value': p,
                'n_celltypes': len(rank1)
            })

            print(f"  {pt1} vs {pt2}: τ={tau:.3f}, p={p:.4f}")

kendall_df = pd.DataFrame(kendall_results)

output_file = output_dir / 'PT_kendall_tau.csv'
kendall_df.to_csv(output_file, index=False)
print(f"\n  Saved: {output_file}")

# ============================================================================
# Step 3-5: Comprehensive Robustness Report
# ============================================================================
print("\n[Step 3-5] Creating robustness summary...")

robustness_summary = {
    'VAT1L_subtype_ordering': {},
    'Microglia_inverse_pattern': {},
    'Cell_type_ordering_correlation': {}
}

# VAT1L ordering
for pt_version in PT_versions:
    pt_data = vat1l_summary[vat1l_summary['PT_version'] == pt_version]

    oligo_pt = pt_data[pt_data['subtype'] == 'Pure Oligo']['PT_mean'].values[0] if len(pt_data[pt_data['subtype'] == 'Pure Oligo']) > 0 else np.nan
    inflam_pt = pt_data[pt_data['subtype'] == 'Oligo-Inflammation']['PT_mean'].values[0] if len(pt_data[pt_data['subtype'] == 'Oligo-Inflammation']) > 0 else np.nan
    upper_pt = pt_data[pt_data['subtype'] == 'Upper-layer']['PT_mean'].values[0] if len(pt_data[pt_data['subtype'] == 'Upper-layer']) > 0 else np.nan

    ordering = sorted([
        ('Pure Oligo', oligo_pt),
        ('Oligo-Inflammation', inflam_pt),
        ('Upper-layer', upper_pt)
    ], key=lambda x: x[1] if not np.isnan(x[1]) else 999)

    robustness_summary['VAT1L_subtype_ordering'][pt_version] = [x[0] for x in ordering if not np.isnan(x[1])]

print("\n  VAT1L ordering across PTs:")
for pt_version, ordering in robustness_summary['VAT1L_subtype_ordering'].items():
    print(f"    {pt_version}: {' < '.join(ordering)}")

# Microglia inverse pattern
for pt_version in PT_versions:
    pt_data = micro_summary[micro_summary['PT_version'] == pt_version]

    oligo_pt = pt_data[pt_data['subtype'] == 'Pure Oligo']['PT_mean'].values[0] if len(pt_data[pt_data['subtype'] == 'Pure Oligo']) > 0 else np.nan
    inflam_pt = pt_data[pt_data['subtype'] == 'Oligo-Inflammation']['PT_mean'].values[0] if len(pt_data[pt_data['subtype'] == 'Oligo-Inflammation']) > 0 else np.nan

    inverse_pattern = inflam_pt > oligo_pt if not (np.isnan(inflam_pt) or np.isnan(oligo_pt)) else None

    robustness_summary['Microglia_inverse_pattern'][pt_version] = inverse_pattern

print("\n  Microglia inverse pattern (Oligo-Inflammation > Pure Oligo):")
for pt_version, pattern in robustness_summary['Microglia_inverse_pattern'].items():
    print(f"    {pt_version}: {pattern}")

# Cell type ordering correlation
robustness_summary['Cell_type_ordering_correlation'] = kendall_df.to_dict('records')

# Save robustness summary
import json
output_file = output_dir / 'PT_robustness_summary.json'
with open(output_file, 'w') as f:
    json.dump(robustness_summary, f, indent=2)
print(f"\n  Saved: {output_file}")

# ============================================================================
# Visualization
# ============================================================================
print("\n[Creating comparison plots...]")

fig = plt.figure(figsize=(20, 16))
gs = fig.add_gridspec(4, 3, hspace=0.3, wspace=0.3)

# Plot 1: PT correlations (3x3 scatter matrix)
for i, pt1 in enumerate(PT_versions):
    for j, pt2 in enumerate(PT_versions):
        ax = fig.add_subplot(gs[i, j])

        if i == j:
            # Diagonal: distribution
            ax.hist(df_als[pt1], bins=50, alpha=0.7, color='blue', edgecolor='black')
            ax.set_ylabel('Frequency' if j == 0 else '')
            ax.set_title(pt1, fontsize=10, fontweight='bold')
        else:
            # Off-diagonal: scatter
            ax.scatter(df_als[pt2], df_als[pt1], alpha=0.05, s=0.5, c='black')
            corr = np.corrcoef(df_als[pt2].dropna(), df_als[pt1].dropna())[0, 1]
            ax.text(0.05, 0.95, f'r={corr:.3f}', transform=ax.transAxes,
                   fontsize=10, verticalalignment='top',
                   bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

        if i == len(PT_versions) - 1:
            ax.set_xlabel(pt2, fontsize=10)
        if j == 0:
            ax.set_ylabel(pt1, fontsize=10)

# Plot 4: VAT1L PT by subtype (all 3 PTs)
ax = fig.add_subplot(gs[3, :])

x_pos = []
labels = []
colors_list = []
positions = {'Pure Oligo': 0, 'Oligo-Inflammation': 1, 'Upper-layer': 2}
offset = {'PT_imes': -0.2, 'PT_dpt': 0, 'PT_stress': 0.2}
colors = {'PT_imes': 'blue', 'PT_dpt': 'green', 'PT_stress': 'red'}

for pt_version in PT_versions:
    pt_data = vat1l_summary[vat1l_summary['PT_version'] == pt_version]

    for _, row in pt_data.iterrows():
        pos = positions[row['subtype']] + offset[pt_version]
        x_pos.append(pos)
        labels.append(f"{row['subtype']}\n{pt_version}")
        colors_list.append(colors[pt_version])

        ax.bar(pos, row['PT_mean'], width=0.18, color=colors[pt_version],
              alpha=0.7, label=pt_version if pos == 0 + offset[pt_version] else "")
        ax.errorbar(pos, row['PT_mean'], yerr=row['PT_std'], fmt='none',
                   color='black', alpha=0.5, capsize=3)

ax.set_xticks([0, 1, 2])
ax.set_xticklabels(['Pure Oligo', 'Oligo-Inflammation', 'Upper-layer'])
ax.set_ylabel('VAT1L PT (mean ± std)', fontsize=12)
ax.set_title('VAT1L PT by Subtype Across 3 PT Versions', fontsize=14, fontweight='bold')
ax.legend()
ax.grid(True, alpha=0.3, axis='y')

plt.suptitle('Phase PTv2: 3-PT Comparison & Robustness Analysis', fontsize=16, fontweight='bold')

plot_file = output_dir / 'PT_comparison_plots.png'
plt.savefig(plot_file, dpi=200, bbox_inches='tight')
plt.close()
print(f"  Saved: {plot_file}")

print("\n" + "="*80)
print("Step 3 completed!")
print(f"Robustness analysis complete for 3 PT versions")
print(f"Key findings tested across all PTs")
print("="*80)
