#!/usr/bin/env python3
"""
Phase 9c+9d: Vascular Time-Lag and Correlation Analysis

Purpose:
  - Analyze vascular heterogeneity using PT_dpt × modules
  - Calculate time-lag between Vasc, Upper, Glia, VAT1L (Phase 8 style)
  - Calculate module correlations between cell types (Phase 7 style)
  - Identify if vascular modules onset before Upper/Glia modules

Author: Claude Code
Date: 2025-11-24
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from scipy.stats import pearsonr, spearmanr
import warnings
warnings.filterwarnings('ignore')

# Configuration
DATA_FILE = 'results/phase9_vascular/cell_level_features_ALL_with_PTdpt.csv'
OUTPUT_DIR = 'results/phase9_vascular'
Path(OUTPUT_DIR).mkdir(parents=True, exist_ok=True)

# Cell type definitions
VASCULAR_TYPES = ['Vasc.Endo.Arterial', 'Vasc.Endo.Capillary', 'Vasc.Endo.Venous',
                  'Vasc.Mural.Pericyte', 'Vasc.Mural.SMC', 'Vasc.Fibro.CLMP.PDGFRA']
UPPER_TYPES = ['Ex.L2.L3.CUX2.RASGRF2']
GLIA_TYPES = ['Glia.Astro.GFAP-neg', 'Glia.Astro.GFAP-pos', 'Glia.Micro', 'Glia.Oligo', 'Glia.OPC']
VAT1L_TYPES = ['Ex.L5.VAT1L.EYA4', 'Ex.L5.VAT1L.THSD4']

def assign_group(cell_type):
    """Assign cell type to group"""
    if cell_type in VAT1L_TYPES:
        return 'VAT1L'
    elif cell_type in UPPER_TYPES:
        return 'Upper'
    elif cell_type in GLIA_TYPES:
        return 'Glia'
    elif cell_type in VASCULAR_TYPES:
        return 'Vascular'
    else:
        return 'Other'

def calculate_onset_for_group(df, module_names, num_bins=25, threshold_multiplier=0.5):
    """
    Calculate module onset PT for each group (Phase 8 style)
    """
    # Determine PT_dpt range
    pt_min = df['PT_dpt'].min()
    pt_99th = df['PT_dpt'].quantile(0.99)
    pt_max = pt_99th * 1.05

    # Create bins
    bin_edges = np.linspace(pt_min, pt_max, num_bins + 1)
    df_copy = df.copy()
    df_copy['pt_bin'] = pd.cut(df_copy['PT_dpt'], bins=bin_edges, labels=False, include_lowest=True)

    # Calculate bin centers
    def get_bin_center(bin_idx):
        if pd.isna(bin_idx):
            return np.nan
        bin_idx = int(bin_idx)
        if bin_idx < 0 or bin_idx >= len(bin_edges) - 1:
            return np.nan
        return (bin_edges[bin_idx] + bin_edges[bin_idx + 1]) / 2

    df_copy['pt_bin_center'] = df_copy['pt_bin'].apply(get_bin_center)

    # Calculate module profiles by group and bin
    profiles_list = []
    for group in ['Upper', 'Glia', 'VAT1L', 'Vascular']:
        group_df = df_copy[df_copy['group'] == group]
        if len(group_df) == 0:
            continue

        for bin_idx in sorted(group_df['pt_bin'].dropna().unique()):
            bin_df = group_df[group_df['pt_bin'] == bin_idx]
            if len(bin_df) == 0:
                continue

            row = {
                'group': group,
                'pt_bin': bin_idx,
                'pt_center': bin_df['pt_bin_center'].iloc[0],
                'n_cells': len(bin_df)
            }

            for mod in module_names:
                row[mod] = bin_df[f'module_{mod}'].mean()

            profiles_list.append(row)

    profiles_df = pd.DataFrame(profiles_list)

    # Calculate onset PT for each group × module
    onset_results = []

    for group in ['Upper', 'Glia', 'VAT1L', 'Vascular']:
        group_profile = profiles_df[profiles_df['group'] == group].copy()
        if len(group_profile) == 0:
            continue

        group_profile = group_profile.sort_values('pt_center')

        for mod in module_names:
            values = group_profile[mod].values
            pt_centers = group_profile['pt_center'].values

            # Calculate threshold
            mean_val = np.nanmean(values)
            std_val = np.nanstd(values)
            onset_threshold = mean_val + threshold_multiplier * std_val

            # Find onset (first bin above threshold)
            onset_idx = None
            for i, val in enumerate(values):
                if not np.isnan(val) and val > onset_threshold:
                    onset_idx = i
                    break

            onset_pt = pt_centers[onset_idx] if onset_idx is not None else np.nan

            onset_results.append({
                'group': group,
                'module': mod,
                'onset_pt': onset_pt,
                'onset_threshold': onset_threshold,
                'onset_bin_idx': onset_idx
            })

    return pd.DataFrame(onset_results), profiles_df

print("="*80)
print("Phase 9c+9d: Vascular Time-Lag and Correlation Analysis")
print("="*80)

# ============================================================================
# Step 1: Load data
# ============================================================================
print("\n[1] Loading data...")
df = pd.read_csv(DATA_FILE)
print(f"Total cells: {len(df):,}")

# Assign groups
df['group'] = df['cell_type'].apply(assign_group)

# Filter to key groups
df = df[df['group'].isin(['Upper', 'Glia', 'VAT1L', 'Vascular'])]
print(f"\nFiltered to {len(df):,} cells in 4 groups:")
for group in ['Vascular', 'Upper', 'Glia', 'VAT1L']:
    print(f"  {group}: {len(df[df['group'] == group]):,} cells")

# Get module names
module_cols = [c for c in df.columns if c.startswith('module_')]
module_names = [c.replace('module_', '') for c in module_cols]
print(f"\nFound {len(module_names)} modules")

# ============================================================================
# Step 2: Calculate onset for each group (Phase 8 style)
# ============================================================================
print("\n[2] Calculating module onset for each group...")
onset_df, profiles_df = calculate_onset_for_group(df, module_names, num_bins=25, threshold_multiplier=0.5)

# Save onset results
onset_df.to_csv(f'{OUTPUT_DIR}/module_onset_by_group.csv', index=False)
print(f"Saved: {OUTPUT_DIR}/module_onset_by_group.csv")

# ============================================================================
# Step 3: Calculate time-lags (Phase 8 style)
# ============================================================================
print("\n[3] Calculating time-lags...")

# Pivot to wide format
onset_wide = onset_df.pivot(index='module', columns='group', values='onset_pt')

# Calculate time-lags
time_lags = onset_wide.copy()
time_lags['Delta_Vasc_minus_Upper'] = onset_wide['Vascular'] - onset_wide['Upper']
time_lags['Delta_Vasc_minus_Glia'] = onset_wide['Vascular'] - onset_wide['Glia']
time_lags['Delta_Vasc_minus_VAT1L'] = onset_wide['Vascular'] - onset_wide['VAT1L']
time_lags['Delta_Upper_minus_Glia'] = onset_wide['Upper'] - onset_wide['Glia']
time_lags['Delta_Upper_minus_VAT1L'] = onset_wide['Upper'] - onset_wide['VAT1L']

time_lags.to_csv(f'{OUTPUT_DIR}/module_time_lags.csv', index=True)
print(f"Saved: {OUTPUT_DIR}/module_time_lags.csv")

# Print summary
print("\nTime-lag summary:")
for comparison in ['Delta_Vasc_minus_Upper', 'Delta_Vasc_minus_Glia', 'Delta_Upper_minus_Glia']:
    delta = time_lags[comparison].dropna()
    if len(delta) > 0:
        mean_delta = delta.mean()
        median_delta = delta.median()
        n_first = (delta < 0).sum()
        n_second = (delta > 0).sum()
        print(f"\n  {comparison}:")
        print(f"    Mean ΔPT: {mean_delta:.6f}")
        print(f"    Median ΔPT: {median_delta:.6f}")
        print(f"    First-onset (Δ<0): {n_first} modules")
        print(f"    Second-onset (Δ>0): {n_second} modules")

# ============================================================================
# Step 4: Module correlations by patient (Phase 7 style)
# ============================================================================
print("\n[4] Calculating module correlations by patient...")

# Check if patient_id exists
if 'patient_id' not in df.columns:
    print("  Warning: patient_id not found, skipping patient-level correlations")
    patient_module_corr = None
else:
    # Calculate mean module scores by patient × cell_type
    patient_means = df.groupby(['patient_id', 'group'])[module_cols].mean().reset_index()

    # Pivot to get separate dataframes for each group
    vasc_means = patient_means[patient_means['group'] == 'Vascular'].set_index('patient_id')[module_cols]
    upper_means = patient_means[patient_means['group'] == 'Upper'].set_index('patient_id')[module_cols]
    glia_means = patient_means[patient_means['group'] == 'Glia'].set_index('patient_id')[module_cols]
    vat1l_means = patient_means[patient_means['group'] == 'VAT1L'].set_index('patient_id')[module_cols]

    # Calculate correlations
    correlations = []

    for vasc_mod in module_names:
        vasc_col = f'module_{vasc_mod}'

        for target_mod in module_names:
            target_col = f'module_{target_mod}'

            # Vasc vs Upper
            common_pts = vasc_means.index.intersection(upper_means.index)
            if len(common_pts) >= 5:
                r, p = pearsonr(vasc_means.loc[common_pts, vasc_col],
                               upper_means.loc[common_pts, target_col])
                correlations.append({
                    'vasc_module': vasc_mod,
                    'target_group': 'Upper',
                    'target_module': target_mod,
                    'r': r,
                    'p': p,
                    'n': len(common_pts)
                })

            # Vasc vs Glia
            common_pts = vasc_means.index.intersection(glia_means.index)
            if len(common_pts) >= 5:
                r, p = pearsonr(vasc_means.loc[common_pts, vasc_col],
                               glia_means.loc[common_pts, target_col])
                correlations.append({
                    'vasc_module': vasc_mod,
                    'target_group': 'Glia',
                    'target_module': target_mod,
                    'r': r,
                    'p': p,
                    'n': len(common_pts)
                })

    patient_module_corr = pd.DataFrame(correlations)
    patient_module_corr.to_csv(f'{OUTPUT_DIR}/vasc_module_correlations.csv', index=False)
    print(f"Saved: {OUTPUT_DIR}/vasc_module_correlations.csv")

    # Print top correlations
    print("\nTop positive correlations (Vasc → Upper):")
    top_upper = patient_module_corr[patient_module_corr['target_group'] == 'Upper'].nlargest(5, 'r')
    for _, row in top_upper.iterrows():
        print(f"  Vasc.{row['vasc_module']:20s} → Upper.{row['target_module']:20s}  r={row['r']:.3f}, p={row['p']:.3f}")

    print("\nTop positive correlations (Vasc → Glia):")
    top_glia = patient_module_corr[patient_module_corr['target_group'] == 'Glia'].nlargest(5, 'r')
    for _, row in top_glia.iterrows():
        print(f"  Vasc.{row['vasc_module']:20s} → Glia.{row['target_module']:20s}  r={row['r']:.3f}, p={row['p']:.3f}")

# ============================================================================
# Step 5: Vascular subtype analysis
# ============================================================================
print("\n[5] Analyzing vascular subtype heterogeneity...")

vasc_df = df[df['group'] == 'Vascular'].copy()

# Calculate mean modules by vascular subtype
vasc_subtype_means = vasc_df.groupby('cell_type')[module_cols + ['PT_dpt', 'stress_total']].mean()
vasc_subtype_means.to_csv(f'{OUTPUT_DIR}/vascular_subtype_module_means.csv', index=True)
print(f"Saved: {OUTPUT_DIR}/vascular_subtype_module_means.csv")

print("\nVascular subtype PT_dpt and stress:")
for ct in VASCULAR_TYPES:
    ct_df = vasc_df[vasc_df['cell_type'] == ct]
    if len(ct_df) > 0:
        print(f"  {ct:30s}  PT={ct_df['PT_dpt'].mean():.3f}  stress={ct_df['stress_total'].mean():.3f}  n={len(ct_df):4d}")

# ============================================================================
# Step 6: Visualizations
# ============================================================================
print("\n[6] Creating visualizations...")

# Figure 1: Time-lag comparison (Phase 8 style)
fig, axes = plt.subplots(1, 3, figsize=(18, 8))

comparisons = [
    ('Delta_Vasc_minus_Upper', 'Vascular vs Upper'),
    ('Delta_Vasc_minus_Glia', 'Vascular vs Glia'),
    ('Delta_Upper_minus_Glia', 'Upper vs Glia')
]

for ax_idx, (col, title) in enumerate(comparisons):
    ax = axes[ax_idx]

    plot_data = time_lags[[col]].dropna().sort_values(col)

    y_pos = np.arange(len(plot_data))
    colors = ['red' if x < 0 else 'blue' for x in plot_data[col]]

    ax.barh(y_pos, plot_data[col], color=colors, alpha=0.7)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(plot_data.index, fontsize=8)
    ax.set_xlabel('ΔPT (Onset Time Difference)', fontsize=10)
    ax.set_title(title, fontsize=12, fontweight='bold')
    ax.axvline(x=0, color='black', linestyle='--', alpha=0.5)
    ax.grid(axis='x', alpha=0.3)

    # Add legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='red', alpha=0.7, label='First-onset (left cell type earlier)'),
        Patch(facecolor='blue', alpha=0.7, label='Second-onset (right cell type earlier)')
    ]
    ax.legend(handles=legend_elements, loc='best', fontsize=8)

plt.tight_layout()
plt.savefig(f'{OUTPUT_DIR}/Fig_module_time_lags.png', dpi=300, bbox_inches='tight')
print(f"Saved: {OUTPUT_DIR}/Fig_module_time_lags.png")
plt.close()

# Figure 2: Correlation heatmaps
if patient_module_corr is not None:
    fig, axes = plt.subplots(1, 2, figsize=(16, 10))

    for ax_idx, target_group in enumerate(['Upper', 'Glia']):
        ax = axes[ax_idx]

        # Pivot to heatmap format
        corr_subset = patient_module_corr[patient_module_corr['target_group'] == target_group]
        corr_matrix = corr_subset.pivot(index='vasc_module', columns='target_module', values='r')

        # Plot heatmap
        sns.heatmap(corr_matrix, cmap='RdBu_r', center=0, vmin=-0.8, vmax=0.8,
                   annot=False, cbar_kws={'label': 'Pearson r'}, ax=ax)
        ax.set_title(f'Vascular → {target_group} Module Correlations', fontsize=12, fontweight='bold')
        ax.set_xlabel(f'{target_group} Modules', fontsize=10)
        ax.set_ylabel('Vascular Modules', fontsize=10)

    plt.tight_layout()
    plt.savefig(f'{OUTPUT_DIR}/Fig_module_correlations_heatmap.png', dpi=300, bbox_inches='tight')
    print(f"Saved: {OUTPUT_DIR}/Fig_module_correlations_heatmap.png")
    plt.close()

# Figure 3: Vascular subtype module profiles
fig, axes = plt.subplots(2, 3, figsize=(18, 12))
axes = axes.flatten()

key_modules = ['Complement', 'Angiogenesis', 'Mitochondria', 'Inflammation',
               'Oxidative_Stress', 'Calcium_Signaling']

for idx, mod in enumerate(key_modules):
    ax = axes[idx]

    subtype_values = []
    subtype_labels = []
    for ct in VASCULAR_TYPES:
        ct_df = vasc_df[vasc_df['cell_type'] == ct]
        if len(ct_df) > 0:
            subtype_values.append(ct_df[f'module_{mod}'].values)
            subtype_labels.append(ct.replace('Vasc.', ''))

    ax.violinplot(subtype_values, positions=range(len(subtype_labels)), showmeans=True)
    ax.set_xticks(range(len(subtype_labels)))
    ax.set_xticklabels(subtype_labels, rotation=45, ha='right', fontsize=8)
    ax.set_ylabel(f'{mod} Score', fontsize=10)
    ax.set_title(f'{mod} by Vascular Subtype', fontsize=10, fontweight='bold')
    ax.grid(axis='y', alpha=0.3)

plt.tight_layout()
plt.savefig(f'{OUTPUT_DIR}/Fig_vascular_subtype_modules.png', dpi=300, bbox_inches='tight')
print(f"Saved: {OUTPUT_DIR}/Fig_vascular_subtype_modules.png")
plt.close()

print("\n" + "="*80)
print("Phase 9c+9d Complete!")
print("="*80)
print(f"\nKey outputs:")
print(f"  - {OUTPUT_DIR}/module_onset_by_group.csv")
print(f"  - {OUTPUT_DIR}/module_time_lags.csv")
print(f"  - {OUTPUT_DIR}/vasc_module_correlations.csv")
print(f"  - {OUTPUT_DIR}/vascular_subtype_module_means.csv")
print(f"  - {OUTPUT_DIR}/Fig_module_time_lags.png")
print(f"  - {OUTPUT_DIR}/Fig_module_correlations_heatmap.png")
print(f"  - {OUTPUT_DIR}/Fig_vascular_subtype_modules.png")
print(f"\nNext: Phase 9e - Integrated report")
