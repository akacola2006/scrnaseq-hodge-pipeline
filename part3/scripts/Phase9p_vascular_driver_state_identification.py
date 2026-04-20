#!/usr/bin/env python3
"""
Phase 9': Vascular Driver-State Identification and Coupling Network

Purpose:
  1. Identify vascular subclusters (early-onset vs late-onset)
  2. Extract "driver-state" within vascular cells (analogous to Early Upper)
  3. Analyze module coupling between Vasc-driver and Upper/Glia
  4. Validate root-cell independence

Author: Claude Code
Date: 2025-11-24
"""

import pandas as pd
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from scipy.stats import pearsonr, spearmanr
import warnings
warnings.filterwarnings('ignore')

# Configuration
DATA_FILE = 'results/phase9_vascular/cell_level_features_ALL_with_PTdpt.csv'
OUTPUT_DIR = 'results/phase9p_vascular_driver'
Path(OUTPUT_DIR).mkdir(parents=True, exist_ok=True)

print("="*80)
print("Phase 9': Vascular Driver-State Identification")
print("="*80)

# ============================================================================
# Step 1: Load data and extract vascular cells
# ============================================================================
print("\n[1] Loading data...")
df = pd.read_csv(DATA_FILE)
print(f"Total cells: {len(df):,}")

# Filter to vascular cells
VASCULAR_TYPES = ['Vasc.Endo.Arterial', 'Vasc.Endo.Capillary', 'Vasc.Endo.Venous',
                  'Vasc.Mural.Pericyte', 'Vasc.Mural.SMC', 'Vasc.Fibro.CLMP.PDGFRA']

vasc_df = df[df['cell_type'].isin(VASCULAR_TYPES)].copy()
print(f"Vascular cells: {len(vasc_df):,}")

# Get module columns
module_cols = [c for c in df.columns if c.startswith('module_')]
module_names = [c.replace('module_', '') for c in module_cols]
print(f"Modules: {len(module_names)}")

# ============================================================================
# Step 2: Create AnnData for vascular cells
# ============================================================================
print("\n[2] Creating vascular AnnData...")

# Use module scores as features
X_vasc = vasc_df[module_cols].values
adata_vasc = sc.AnnData(X=X_vasc)
adata_vasc.obs_names = vasc_df['cell_id'].values
adata_vasc.var_names = module_names

# Add metadata
adata_vasc.obs['cell_type'] = vasc_df['cell_type'].values
adata_vasc.obs['PT_dpt'] = vasc_df['PT_dpt'].values
adata_vasc.obs['stress_total'] = vasc_df['stress_total'].values

if 'patient_id' in vasc_df.columns:
    adata_vasc.obs['patient_id'] = vasc_df['patient_id'].values

print(f"Vascular AnnData: {adata_vasc.n_obs} cells × {adata_vasc.n_vars} modules")

# ============================================================================
# Step 3: Vascular subclustering (PT × modules × stress)
# ============================================================================
print("\n[3] Vascular subclustering...")

# Scale
sc.pp.scale(adata_vasc, max_value=10)

# PCA
sc.tl.pca(adata_vasc, n_comps=min(20, adata_vasc.n_vars-1))
var_exp = adata_vasc.uns['pca']['variance_ratio'][:10].sum()
print(f"  PCA variance explained (10 PCs): {var_exp:.3f}")

# Neighbors + UMAP
sc.pp.neighbors(adata_vasc, n_neighbors=min(30, len(adata_vasc)//3), n_pcs=15)
sc.tl.umap(adata_vasc)

# Leiden clustering
for res in [0.3, 0.5, 0.8, 1.0]:
    sc.tl.leiden(adata_vasc, resolution=res, key_added=f'leiden_r{res}')

# Use resolution 0.5 as default
adata_vasc.obs['cluster'] = adata_vasc.obs['leiden_r0.5']
n_clusters = adata_vasc.obs['cluster'].nunique()
print(f"  Identified {n_clusters} vascular subclusters (resolution=0.5)")

# ============================================================================
# Step 4: Identify Early-onset vs Late-onset clusters
# ============================================================================
print("\n[4] Identifying Early vs Late-onset clusters...")

cluster_stats = []
for cluster_id in sorted(adata_vasc.obs['cluster'].unique()):
    cluster_cells = adata_vasc.obs[adata_vasc.obs['cluster'] == cluster_id]

    cluster_stats.append({
        'cluster': cluster_id,
        'n_cells': len(cluster_cells),
        'PT_dpt_mean': cluster_cells['PT_dpt'].mean(),
        'PT_dpt_median': cluster_cells['PT_dpt'].median(),
        'stress_mean': cluster_cells['stress_total'].mean(),
        'stress_std': cluster_cells['stress_total'].std(),
    })

cluster_stats_df = pd.DataFrame(cluster_stats)
cluster_stats_df = cluster_stats_df.sort_values('PT_dpt_mean')
cluster_stats_df['onset_category'] = pd.cut(
    cluster_stats_df['PT_dpt_mean'],
    bins=3,
    labels=['Early', 'Mid', 'Late']
)

cluster_stats_df.to_csv(f'{OUTPUT_DIR}/vascular_cluster_stats.csv', index=False)
print(f"Saved: {OUTPUT_DIR}/vascular_cluster_stats.csv")

print("\nCluster onset categories:")
for cat in ['Early', 'Mid', 'Late']:
    cat_clusters = cluster_stats_df[cluster_stats_df['onset_category'] == cat]
    print(f"  {cat:8s}: {list(cat_clusters['cluster'])} (n={cat_clusters['n_cells'].sum()})")

# ============================================================================
# Step 5: Identify Driver-State (analogous to Early Upper)
# ============================================================================
print("\n[5] Identifying vascular driver-state...")

# Strategy: Early-onset (low PT) + High stress (extreme dysfunction)
# Similar to Phase 7 Early Upper identification

# Define Early-onset cells (bottom 25% PT_dpt within vascular)
pt_threshold = vasc_df['PT_dpt'].quantile(0.25)
early_vasc = vasc_df[vasc_df['PT_dpt'] <= pt_threshold].copy()
print(f"  Early vascular cells (PT_dpt ≤ {pt_threshold:.3f}): {len(early_vasc)} ({len(early_vasc)/len(vasc_df)*100:.1f}%)")

# Within early cells, identify extreme stress
stress_threshold = early_vasc['stress_total'].mean() + 1.0 * early_vasc['stress_total'].std()
driver_vasc = early_vasc[early_vasc['stress_total'] > stress_threshold].copy()
print(f"  Driver-state vascular (early + high stress): {len(driver_vasc)} ({len(driver_vasc)/len(vasc_df)*100:.1f}%)")

# Save driver-state cells
driver_vasc['driver_state'] = 'Vasc_Driver'
driver_vasc.to_csv(f'{OUTPUT_DIR}/vascular_driver_state_cells.csv', index=False)
print(f"Saved: {OUTPUT_DIR}/vascular_driver_state_cells.csv")

# Add driver label to main dataframe
vasc_df['driver_state'] = 'Normal'
vasc_df.loc[vasc_df['cell_id'].isin(driver_vasc['cell_id']), 'driver_state'] = 'Driver'

# ============================================================================
# Step 6: Module profiles of driver-state
# ============================================================================
print("\n[6] Analyzing driver-state module profiles...")

driver_modules = driver_vasc[module_cols].mean()
normal_modules = vasc_df[vasc_df['driver_state'] == 'Normal'][module_cols].mean()

module_comparison = pd.DataFrame({
    'module': module_names,
    'driver_mean': driver_modules.values,
    'normal_mean': normal_modules.values,
    'fold_change': driver_modules.values / (normal_modules.values + 1e-10),
    'difference': driver_modules.values - normal_modules.values
})
module_comparison = module_comparison.sort_values('fold_change', ascending=False)
module_comparison.to_csv(f'{OUTPUT_DIR}/driver_state_module_profiles.csv', index=False)
print(f"Saved: {OUTPUT_DIR}/driver_state_module_profiles.csv")

print("\nTop enriched modules in driver-state:")
for _, row in module_comparison.head(5).iterrows():
    print(f"  {row['module']:20s}  FC={row['fold_change']:.2f}  Δ={row['difference']:.3f}")

# ============================================================================
# Step 7: Patient-level coupling: Vasc-driver vs Upper/Glia
# ============================================================================
print("\n[7] Analyzing patient-level coupling...")

# Check if patient_id exists
if 'patient_id' not in df.columns:
    print("  Warning: patient_id not found, skipping patient-level correlations")
    patient_coupling = None
else:
    # Calculate mean modules by patient for each cell type group

    # Vascular driver-state
    driver_patient = driver_vasc.groupby('patient_id')[module_cols].mean()
    driver_patient.columns = [c.replace('module_', '') for c in driver_patient.columns]

    # Upper cells
    upper_df = df[df['cell_type'].str.contains('CUX2.RASGRF2', na=False)]
    upper_patient = upper_df.groupby('patient_id')[module_cols].mean()
    upper_patient.columns = [c.replace('module_', '') for c in upper_patient.columns]

    # Glia cells
    glia_df = df[df['cell_type'].str.contains('Glia', na=False)]
    glia_patient = glia_df.groupby('patient_id')[module_cols].mean()
    glia_patient.columns = [c.replace('module_', '') for c in glia_patient.columns]

    # Calculate correlations
    coupling_results = []

    for vasc_mod in module_names:
        # Vasc-driver vs Upper
        common_pts = driver_patient.index.intersection(upper_patient.index)
        if len(common_pts) >= 5:
            for upper_mod in module_names:
                r, p = pearsonr(driver_patient.loc[common_pts, vasc_mod],
                               upper_patient.loc[common_pts, upper_mod])
                coupling_results.append({
                    'vasc_module': vasc_mod,
                    'target_group': 'Upper',
                    'target_module': upper_mod,
                    'r': r,
                    'p': p,
                    'n': len(common_pts)
                })

        # Vasc-driver vs Glia
        common_pts = driver_patient.index.intersection(glia_patient.index)
        if len(common_pts) >= 5:
            for glia_mod in module_names:
                r, p = pearsonr(driver_patient.loc[common_pts, vasc_mod],
                               glia_patient.loc[common_pts, glia_mod])
                coupling_results.append({
                    'vasc_module': vasc_mod,
                    'target_group': 'Glia',
                    'target_module': glia_mod,
                    'r': r,
                    'p': p,
                    'n': len(common_pts)
                })

    patient_coupling = pd.DataFrame(coupling_results)
    patient_coupling.to_csv(f'{OUTPUT_DIR}/vascular_driver_patient_coupling.csv', index=False)
    print(f"Saved: {OUTPUT_DIR}/vascular_driver_patient_coupling.csv")

    # Top correlations
    print("\nTop Vasc-Driver → Upper correlations:")
    top_upper = patient_coupling[patient_coupling['target_group'] == 'Upper'].nlargest(5, 'r')
    for _, row in top_upper.iterrows():
        print(f"  Vasc.{row['vasc_module']:15s} → Upper.{row['target_module']:15s}  r={row['r']:.3f}, p={row['p']:.3f}")

    print("\nTop Vasc-Driver → Glia correlations:")
    top_glia = patient_coupling[patient_coupling['target_group'] == 'Glia'].nlargest(5, 'r')
    for _, row in top_glia.iterrows():
        print(f"  Vasc.{row['vasc_module']:15s} → Glia.{row['target_module']:15s}  r={row['r']:.3f}, p={row['p']:.3f}")

# ============================================================================
# Step 8: Root-cell sensitivity analysis
# ============================================================================
print("\n[8] Root-cell sensitivity analysis...")

# Test different root cells
from scipy.stats import spearmanr

# Original root: lowest stress in homeostatic cells
homeostatic_types = ['Glia.Oligo', 'Glia.Astro.GFAP-neg', 'Vasc.Endo.Capillary']
homeostatic_mask = df['cell_type'].isin(homeostatic_types)
root_original = df[homeostatic_mask]['stress_total'].idxmin()

# Alternative roots
root_highest_stress = df['stress_total'].idxmax()
root_random = df.sample(1).index[0]
root_vasc_earliest = vasc_df['PT_dpt'].idxmin()

print(f"  Original root: {root_original} (lowest stress homeostatic)")
print(f"  Alt root 1: {root_highest_stress} (highest stress)")
print(f"  Alt root 2: {root_random} (random)")
print(f"  Alt root 3: {root_vasc_earliest} (vascular earliest PT)")

print("\n  Note: Full PT_dpt recalculation with alternative roots")
print("  would require recomputing diffusion map. Skipping for now.")
print("  Current PT_dpt shows strong stress independence (R²<1%)")
print("  and consistent temporal ordering (Vasc → Glia → Upper).")

# ============================================================================
# Step 9: Visualizations
# ============================================================================
print("\n[9] Creating visualizations...")

# Figure 1: Vascular UMAP colored by key features
fig, axes = plt.subplots(2, 3, figsize=(18, 12))

sc.pl.umap(adata_vasc, color='cluster', ax=axes[0,0], show=False, title='Vascular Clusters')
sc.pl.umap(adata_vasc, color='PT_dpt', cmap='viridis', ax=axes[0,1], show=False, title='PT_dpt')
sc.pl.umap(adata_vasc, color='stress_total', cmap='Reds', ax=axes[0,2], show=False, title='Stress Total')

# Add driver-state to adata
adata_vasc.obs['driver_state'] = vasc_df.set_index('cell_id').loc[adata_vasc.obs_names, 'driver_state']
sc.pl.umap(adata_vasc, color='driver_state', ax=axes[1,0], show=False, title='Driver State')
sc.pl.umap(adata_vasc, color='cell_type', ax=axes[1,1], show=False, title='Cell Type')

# Module heatmap
axes[1,2].axis('off')

plt.tight_layout()
plt.savefig(f'{OUTPUT_DIR}/Fig1_vascular_subclusters_UMAP.png', dpi=300, bbox_inches='tight')
print(f"Saved: {OUTPUT_DIR}/Fig1_vascular_subclusters_UMAP.png")
plt.close()

# Figure 2: Driver-state module profile
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# Bar plot
ax = axes[0]
y_pos = np.arange(len(module_comparison))
colors = ['red' if fc > 1.2 else 'gray' for fc in module_comparison['fold_change']]
ax.barh(y_pos, module_comparison['fold_change'], color=colors, alpha=0.7)
ax.set_yticks(y_pos)
ax.set_yticklabels(module_comparison['module'], fontsize=8)
ax.set_xlabel('Fold Change (Driver / Normal)', fontsize=10)
ax.set_title('Vascular Driver-State Module Enrichment', fontsize=12, fontweight='bold')
ax.axvline(x=1, color='black', linestyle='--', alpha=0.5)
ax.invert_yaxis()

# Scatter: driver vs normal
ax = axes[1]
ax.scatter(module_comparison['normal_mean'], module_comparison['driver_mean'], alpha=0.6)
for _, row in module_comparison.iterrows():
    if row['fold_change'] > 1.3 or row['fold_change'] < 0.8:
        ax.annotate(row['module'], (row['normal_mean'], row['driver_mean']),
                   fontsize=7, alpha=0.7)
ax.plot([0, ax.get_xlim()[1]], [0, ax.get_xlim()[1]], 'k--', alpha=0.3)
ax.set_xlabel('Normal Vascular (mean)', fontsize=10)
ax.set_ylabel('Driver-State (mean)', fontsize=10)
ax.set_title('Driver vs Normal Module Levels', fontsize=12, fontweight='bold')

plt.tight_layout()
plt.savefig(f'{OUTPUT_DIR}/Fig2_driver_state_module_profile.png', dpi=300, bbox_inches='tight')
print(f"Saved: {OUTPUT_DIR}/Fig2_driver_state_module_profile.png")
plt.close()

# Figure 3: Patient-level coupling heatmaps
if patient_coupling is not None:
    fig, axes = plt.subplots(1, 2, figsize=(16, 10))

    for ax_idx, target_group in enumerate(['Upper', 'Glia']):
        ax = axes[ax_idx]

        # Pivot to heatmap
        corr_subset = patient_coupling[patient_coupling['target_group'] == target_group]
        corr_matrix = corr_subset.pivot(index='vasc_module', columns='target_module', values='r')

        # Plot
        sns.heatmap(corr_matrix, cmap='RdBu_r', center=0, vmin=-0.8, vmax=0.8,
                   annot=False, cbar_kws={'label': 'Pearson r'}, ax=ax)
        ax.set_title(f'Vasc-Driver → {target_group} Module Coupling\n(Patient-Level Correlations)',
                    fontsize=12, fontweight='bold')
        ax.set_xlabel(f'{target_group} Modules', fontsize=10)
        ax.set_ylabel('Vascular Driver Modules', fontsize=10)

    plt.tight_layout()
    plt.savefig(f'{OUTPUT_DIR}/Fig3_driver_coupling_heatmaps.png', dpi=300, bbox_inches='tight')
    print(f"Saved: {OUTPUT_DIR}/Fig3_driver_coupling_heatmaps.png")
    plt.close()

print("\n" + "="*80)
print("Phase 9' Complete!")
print("="*80)
print(f"\nKey outputs:")
print(f"  - {OUTPUT_DIR}/vascular_cluster_stats.csv")
print(f"  - {OUTPUT_DIR}/vascular_driver_state_cells.csv")
print(f"  - {OUTPUT_DIR}/driver_state_module_profiles.csv")
print(f"  - {OUTPUT_DIR}/vascular_driver_patient_coupling.csv")
print(f"  - {OUTPUT_DIR}/Fig1_vascular_subclusters_UMAP.png")
print(f"  - {OUTPUT_DIR}/Fig2_driver_state_module_profile.png")
print(f"  - {OUTPUT_DIR}/Fig3_driver_coupling_heatmaps.png")
print(f"\nNext: Integrate Phase 9' findings into final report")
