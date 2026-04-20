#!/usr/bin/env python3
"""
Phase TDP NVU Landscape Analysis
=================================

Purpose: Map PT_TDP (TDP-43 pathology pseudotime) onto Phase 13 NVU framework
to reveal TDP-43 abnormality distribution across NVU nodes and disease progression.

Uses same methodology as Phase 13:
- NVU nodes: Vascular, Glia, Upper, VAT1L
- PT_dpt bins: ~30 bins along disease pseudotime
- Analysis: Mean PT_TDP per node × PT_dpt bin
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from scipy import stats

# =====================================================
# Configuration
# =====================================================

BASE_DIR = Path("/home/akaco/als/motor_cortex_analysis")
IDS_DIR = BASE_DIR / "ids_causal_analysis"
RESULTS_DIR = IDS_DIR / "results" / "phaseTDP_NVU"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

# NVU Node definitions - LAYER-SPECIFIC
NVU_NODES = {
    'Vascular': [
        'Vasc.Endo.Capillary', 'Vasc.Endo.Venous', 'Vasc.Endo.Arterial',
        'Vasc.Mural.Pericyte', 'Vasc.Mural.SMC', 'Vasc.Fibro.CLMP.PDGFRA'
    ],
    'Glia': [
        'Glia.Oligo', 'Glia.Astro.GFAP-neg', 'Glia.Astro.GFAP-pos',
        'Glia.OPC', 'Glia.Micro'
    ],
    'L2/L3': [
        'Ex.L2.L3.CUX2.RASGRF2'
    ],
    'L3-L6': [
        'Ex.L3.L5.CUX2.RORB', 'Ex.L4.L5.RORB.FOXO1',
        'Ex.L4.L6.RORB.LRRK1'
    ],
    'VAT1L': [
        'Ex.L5.VAT1L.EYA4', 'Ex.L5.VAT1L.THSD4'
    ]
}

print("="*60)
print("Phase TDP NVU Landscape Analysis")
print("="*60)

# =====================================================
# Load Data
# =====================================================

print("\nLoading PT_TDP results...")
tdp_file = IDS_DIR / "results" / "phaseTDP_proper" / "cell_level_features_with_PT_TDP_proper.csv"
df = pd.read_csv(tdp_file)

print(f"  Loaded {len(df):,} cells")
print(f"  Groups: {df['group'].value_counts().to_dict()}")

# Filter to ALS only with valid PT_TDP and PT_dpt
df_als = df[df['group'] == 'ALS'].copy()
df_als = df_als[df_als['PT_TDP'].notna() & df_als['PT_dpt'].notna()]

print(f"  ALS cells with PT_TDP & PT_dpt: {len(df_als):,}")

# Assign NVU nodes
df_als['NVU_node'] = 'Other'

for node_name, cell_types in NVU_NODES.items():
    mask = df_als['cell_type'].isin(cell_types)
    df_als.loc[mask, 'NVU_node'] = node_name

# Count cells per node
node_counts = df_als['NVU_node'].value_counts()
print("\n  Cells per NVU node:")
for node, count in node_counts.items():
    print(f"    {node}: {count:,}")

# Remove "Other" cells
df_als = df_als[df_als['NVU_node'] != 'Other']

print(f"\n  Total NVU cells for analysis: {len(df_als):,}")

# =====================================================
# PT_dpt Binning (Phase 13 style)
# =====================================================

print("\nBinning PT_dpt...")

n_bins = 30
pt_min = df_als['PT_dpt'].min()
pt_max = df_als['PT_dpt'].max()

print(f"  PT_dpt range: [{pt_min:.3f}, {pt_max:.3f}]")
print(f"  Number of bins: {n_bins}")

# Create bins
pt_bins = np.linspace(pt_min, pt_max, n_bins + 1)
bin_centers = (pt_bins[:-1] + pt_bins[1:]) / 2

df_als['PT_bin'] = pd.cut(df_als['PT_dpt'], bins=pt_bins, labels=False, include_lowest=True)

# =====================================================
# Compute PT_TDP profiles per NVU node
# =====================================================

print("\nComputing PT_TDP(PT_dpt) profiles per NVU node...")

profiles = []

for node in ['Vascular', 'Glia', 'L2/L3', 'L3-L6', 'VAT1L']:
    df_node = df_als[df_als['NVU_node'] == node]

    if len(df_node) == 0:
        continue

    for bin_idx in range(n_bins):
        df_bin = df_node[df_node['PT_bin'] == bin_idx]

        if len(df_bin) < 5:  # Minimum 5 cells per bin
            continue

        mean_pt_tdp = df_bin['PT_TDP'].mean()
        std_pt_tdp = df_bin['PT_TDP'].std()
        n_cells = len(df_bin)
        pt_center = bin_centers[bin_idx]

        profiles.append({
            'NVU_node': node,
            'PT_bin': bin_idx,
            'PT_dpt_center': pt_center,
            'mean_PT_TDP': mean_pt_tdp,
            'std_PT_TDP': std_pt_tdp,
            'n_cells': n_cells
        })

df_profiles = pd.DataFrame(profiles)

# Save profiles
profile_file = RESULTS_DIR / "PT_TDP_NVU_profiles.csv"
df_profiles.to_csv(profile_file, index=False)
print(f"  Saved: {profile_file.name}")

# =====================================================
# Visualization 1: PT_TDP trajectories per NVU node
# =====================================================

print("\nGenerating PT_TDP trajectory plot...")

fig, ax = plt.subplots(figsize=(12, 8))

colors = {'Vascular': '#e74c3c', 'Glia': '#f39c12', 'L2/L3': '#3498db', 'L3-L6': '#9b59b6', 'VAT1L': '#2ecc71'}
markers = {'Vascular': 'o', 'Glia': 's', 'L2/L3': '^', 'L3-L6': 'v', 'VAT1L': 'D'}

for node in ['Vascular', 'Glia', 'L2/L3', 'L3-L6', 'VAT1L']:
    df_node = df_profiles[df_profiles['NVU_node'] == node]

    if len(df_node) == 0:
        continue

    ax.plot(df_node['PT_dpt_center'], df_node['mean_PT_TDP'],
            marker=markers[node], linewidth=2.5, markersize=8,
            label=f'{node} (n={node_counts[node]:,})',
            color=colors[node], alpha=0.8)

    # Error bars
    ax.fill_between(df_node['PT_dpt_center'],
                    df_node['mean_PT_TDP'] - df_node['std_PT_TDP'],
                    df_node['mean_PT_TDP'] + df_node['std_PT_TDP'],
                    alpha=0.2, color=colors[node])

ax.set_xlabel('PT_dpt (Disease Pseudotime)', fontsize=14, fontweight='bold')
ax.set_ylabel('PT_TDP (TDP-43 Pathology Level)', fontsize=14, fontweight='bold')
ax.set_title('TDP-43 Pathology Distribution Across NVU Nodes\n(Phase 13 Framework)',
             fontsize=16, fontweight='bold')
ax.legend(fontsize=11, loc='best', framealpha=0.9)
ax.grid(alpha=0.3, linestyle='--')

plt.tight_layout()
plt.savefig(RESULTS_DIR / 'Fig_PT_TDP_NVU_trajectories.png', dpi=300, bbox_inches='tight')
plt.close()

print(f"  Saved: Fig_PT_TDP_NVU_trajectories.png")

# =====================================================
# Visualization 2: Heatmap (NVU nodes × PT_dpt bins)
# =====================================================

print("\nGenerating PT_TDP heatmap...")

# Pivot for heatmap
pivot_data = df_profiles.pivot(index='NVU_node', columns='PT_bin', values='mean_PT_TDP')

# Reorder rows
row_order = ['Vascular', 'Glia', 'L2/L3', 'L3-L6', 'VAT1L']
pivot_data = pivot_data.reindex(row_order)

fig, ax = plt.subplots(figsize=(16, 6))

sns.heatmap(pivot_data, cmap='RdYlBu_r', center=0.15,
            cbar_kws={'label': 'Mean PT_TDP (TDP-43 Pathology)'},
            linewidths=0.5, linecolor='white',
            ax=ax, vmin=0, vmax=0.3)

ax.set_xlabel('PT_dpt Bin (Disease Progression →)', fontsize=13, fontweight='bold')
ax.set_ylabel('NVU Node', fontsize=13, fontweight='bold')
ax.set_title('TDP-43 Pathology Landscape: NVU Nodes × Disease Pseudotime\n(Phase 13 Framework)',
             fontsize=15, fontweight='bold')

plt.tight_layout()
plt.savefig(RESULTS_DIR / 'Fig_PT_TDP_NVU_heatmap.png', dpi=300, bbox_inches='tight')
plt.close()

print(f"  Saved: Fig_PT_TDP_NVU_heatmap.png")

# =====================================================
# Statistical Summary
# =====================================================

print("\nComputing statistical summary...")

summary_stats = []

for node in ['Vascular', 'Glia', 'L2/L3', 'L3-L6', 'VAT1L']:
    df_node = df_als[df_als['NVU_node'] == node]

    if len(df_node) == 0:
        continue

    mean_pt_tdp = df_node['PT_TDP'].mean()
    median_pt_tdp = df_node['PT_TDP'].median()
    std_pt_tdp = df_node['PT_TDP'].std()
    min_pt_tdp = df_node['PT_TDP'].min()
    max_pt_tdp = df_node['PT_TDP'].max()

    # Correlation with PT_dpt
    r_corr, p_corr = stats.pearsonr(df_node['PT_dpt'], df_node['PT_TDP'])

    # Early vs late PT_dpt
    early_mask = df_node['PT_dpt'] < 0.3
    late_mask = df_node['PT_dpt'] > 0.5

    early_pt_tdp = df_node[early_mask]['PT_TDP'].mean() if early_mask.sum() > 0 else np.nan
    late_pt_tdp = df_node[late_mask]['PT_TDP'].mean() if late_mask.sum() > 0 else np.nan

    summary_stats.append({
        'NVU_node': node,
        'n_cells': len(df_node),
        'mean_PT_TDP': mean_pt_tdp,
        'median_PT_TDP': median_pt_tdp,
        'std_PT_TDP': std_pt_tdp,
        'min_PT_TDP': min_pt_tdp,
        'max_PT_TDP': max_pt_tdp,
        'corr_PT_dpt_PT_TDP': r_corr,
        'corr_pval': p_corr,
        'early_PT_TDP': early_pt_tdp,
        'late_PT_TDP': late_pt_tdp
    })

df_summary = pd.DataFrame(summary_stats)

# Save summary
summary_file = RESULTS_DIR / "PT_TDP_NVU_summary_stats.csv"
df_summary.to_csv(summary_file, index=False)

print("\n  PT_TDP Statistics per NVU Node:")
print(df_summary.to_string(index=False))
print(f"\n  Saved: {summary_file.name}")

# =====================================================
# Interpretation
# =====================================================

print("\n" + "="*60)
print("Key Findings:")
print("="*60)

# Rank by mean PT_TDP
df_summary_sorted = df_summary.sort_values('mean_PT_TDP', ascending=False)

print("\n1. TDP-43 Pathology Severity (by mean PT_TDP):")
for idx, row in df_summary_sorted.iterrows():
    print(f"   {row['NVU_node']}: {row['mean_PT_TDP']:.3f} ± {row['std_PT_TDP']:.3f}")

print("\n2. TDP-43 Pathology vs Disease Progression Correlation:")
for idx, row in df_summary_sorted.iterrows():
    sig = "***" if row['corr_pval'] < 0.001 else "**" if row['corr_pval'] < 0.01 else "*" if row['corr_pval'] < 0.05 else "ns"
    print(f"   {row['NVU_node']}: r = {row['corr_PT_dpt_PT_TDP']:.3f} ({sig})")

print("\n3. Early vs Late Disease Stage:")
for idx, row in df_summary_sorted.iterrows():
    if not np.isnan(row['early_PT_TDP']) and not np.isnan(row['late_PT_TDP']):
        change = row['late_PT_TDP'] - row['early_PT_TDP']
        direction = "↑" if change > 0 else "↓"
        print(f"   {row['NVU_node']}: {row['early_PT_TDP']:.3f} → {row['late_PT_TDP']:.3f} ({direction}{abs(change):.3f})")

print("\n" + "="*60)
print("Phase TDP NVU Landscape Complete!")
print(f"Results saved to: {RESULTS_DIR}")
print("="*60)
