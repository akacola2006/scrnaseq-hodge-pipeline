#!/usr/bin/env python3
"""
Phase 9a: Recalculate PT_dpt for ALL Cell Types (including Vascular)

Purpose:
  - Load cell_level_features_ALL.csv with all cell types
  - Calculate PT_dpt using diffusion pseudotime (stress-independent)
  - Save results to cell_level_features_ALL_with_PTdpt.csv

Rationale:
  - PT_IDS/PT_imes are stress-dependent (R²=62.5% for VAT1L)
  - Need fair comparison of vascular vs Upper/Glia temporal ordering
  - Use module scores (23D) as feature space for diffusion map

Author: Claude Code
Date: 2025-11-24
"""

import pandas as pd
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# Configuration
DATA_FILE = 'results/cell_state_causality/cell_level_features_ALL.csv'
OUTPUT_DIR = 'results/phase9_vascular'
OUTPUT_FILE = f'{OUTPUT_DIR}/cell_level_features_ALL_with_PTdpt.csv'
Path(OUTPUT_DIR).mkdir(parents=True, exist_ok=True)

print("="*80)
print("Phase 9a: Recalculate PT_dpt for ALL Cell Types")
print("="*80)

# ============================================================================
# Step 1: Load data
# ============================================================================
print("\n[1] Loading data...")
df = pd.read_csv(DATA_FILE)
print(f"Total cells: {len(df):,}")
print(f"Columns: {len(df.columns)}")

# Get module columns
module_cols = [c for c in df.columns if c.startswith('module_')]
print(f"Found {len(module_cols)} modules")

# Check cell types
print(f"\nCell types: {df['cell_type'].nunique()}")
print("\nTop cell types:")
print(df['cell_type'].value_counts().head(10))

# Identify vascular cells
vascular_cts = df['cell_type'][df['cell_type'].str.contains('Vasc', na=False)].unique()
print(f"\nVascular cell types ({len(vascular_cts)}):")
for ct in sorted(vascular_cts):
    n = len(df[df['cell_type'] == ct])
    print(f"  {ct}: {n:,} cells")

# ============================================================================
# Step 2: Create AnnData object from module scores
# ============================================================================
print("\n[2] Creating AnnData object...")

# Extract module scores as X matrix
X = df[module_cols].values
print(f"Feature matrix shape: {X.shape}")

# Create AnnData
adata = sc.AnnData(X=X)
adata.obs_names = df['cell_id'].astype(str).values
adata.var_names = [c.replace('module_', '') for c in module_cols]

# Add metadata
adata.obs['cell_type'] = df['cell_type'].values
adata.obs['condition'] = df['condition'].values
adata.obs['PT_IDS'] = df['PT_IDS'].values
adata.obs['stress_total'] = df['stress_total'].values

# Add patient_id if exists
if 'patient_id' in df.columns:
    adata.obs['patient_id'] = df['patient_id'].values

print(f"AnnData: {adata.n_obs} cells × {adata.n_vars} modules")

# ============================================================================
# Step 3: Preprocessing for diffusion map
# ============================================================================
print("\n[3] Preprocessing...")

# Normalize module scores (already 0-1 range, but ensure proper scaling)
sc.pp.scale(adata, max_value=10)

# PCA
print("  Running PCA...")
sc.tl.pca(adata, svd_solver='full', n_comps=20)  # Use 20 PCs from 23 modules

# Calculate explained variance
var_explained = adata.uns['pca']['variance_ratio'][:10]
print(f"  PCA variance explained (first 10 PCs): {var_explained.sum():.3f}")

# ============================================================================
# Step 4: Compute neighbor graph
# ============================================================================
print("\n[4] Computing neighbor graph...")
sc.pp.neighbors(adata, n_neighbors=30, n_pcs=20, method='umap')

# UMAP for visualization
print("  Computing UMAP...")
sc.tl.umap(adata)

# ============================================================================
# Step 5: Diffusion map and PT_dpt
# ============================================================================
print("\n[5] Computing diffusion map...")
sc.tl.diffmap(adata, n_comps=15)

# Determine root cell for PT_dpt
# Strategy: Use cell with lowest stress_total in homeostatic cell types
# (Glia, Endo, Fibro - representing normal homeostatic state)
print("\n[6] Determining root cell for PT_dpt...")

homeostatic_types = ['Glia.Oligo', 'Glia.Astro.GFAP-neg', 'Glia.OPC',
                     'Vasc.Endo.Capillary', 'Vasc.Fibro.CLMP.PDGFRA']
homeostatic_mask = adata.obs['cell_type'].isin(homeostatic_types)

if homeostatic_mask.sum() > 0:
    homeostatic_cells = adata.obs[homeostatic_mask]
    root_idx = homeostatic_cells['stress_total'].idxmin()
    root_cell_idx = adata.obs_names.get_loc(root_idx)
    root_cell_type = adata.obs.loc[root_idx, 'cell_type']
    root_stress = adata.obs.loc[root_idx, 'stress_total']
    print(f"  Root cell: {root_idx}")
    print(f"  Cell type: {root_cell_type}")
    print(f"  Stress: {root_stress:.4f}")
else:
    # Fallback: lowest stress overall
    root_idx = adata.obs['stress_total'].idxmin()
    root_cell_idx = adata.obs_names.get_loc(root_idx)
    root_cell_type = adata.obs.loc[root_idx, 'cell_type']
    root_stress = adata.obs.loc[root_idx, 'stress_total']
    print(f"  Root cell (fallback): {root_idx}")
    print(f"  Cell type: {root_cell_type}")
    print(f"  Stress: {root_stress:.4f}")

# Compute DPT
print("\n[7] Computing diffusion pseudotime (PT_dpt)...")
adata.uns['iroot'] = root_cell_idx
sc.tl.dpt(adata, n_dcs=10)

# Extract PT_dpt
pt_dpt = adata.obs['dpt_pseudotime'].values
print(f"  PT_dpt range: {pt_dpt.min():.4f} - {pt_dpt.max():.4f}")
print(f"  PT_dpt mean: {pt_dpt.mean():.4f}")
print(f"  PT_dpt median: {np.median(pt_dpt):.4f}")

# ============================================================================
# Step 8: Check PT_dpt vs stress correlation
# ============================================================================
print("\n[8] Validating PT_dpt independence from stress...")

from scipy.stats import pearsonr, spearmanr

# Overall correlation
r_pearson, p_pearson = pearsonr(pt_dpt, adata.obs['stress_total'])
r_spearman, p_spearman = spearmanr(pt_dpt, adata.obs['stress_total'])
r2 = r_pearson ** 2

print(f"  Overall PT_dpt vs stress_total:")
print(f"    Pearson r: {r_pearson:.4f} (R² = {r2:.4f}, p = {p_pearson:.2e})")
print(f"    Spearman ρ: {r_spearman:.4f} (p = {p_spearman:.2e})")

# Per cell type
print(f"\n  PT_dpt vs stress_total by cell type:")
key_celltypes = ['Ex.L2.L3.CUX2.RASGRF2', 'Glia.Oligo', 'Glia.Astro.GFAP-neg',
                 'Ex.L5.VAT1L.EYA4', 'Ex.L5.VAT1L.THSD4',
                 'Vasc.Endo.Capillary', 'Vasc.Mural.Pericyte']

for ct in key_celltypes:
    ct_mask = adata.obs['cell_type'] == ct
    if ct_mask.sum() > 10:
        ct_pt = pt_dpt[ct_mask]
        ct_stress = adata.obs.loc[ct_mask, 'stress_total']
        r, p = pearsonr(ct_pt, ct_stress)
        r2_ct = r ** 2
        print(f"    {ct:30s} r={r:6.3f} (R²={r2_ct:.3f}, n={ct_mask.sum():4d})")

# ============================================================================
# Step 9: Save results
# ============================================================================
print("\n[9] Saving results...")

# Add PT_dpt to original dataframe
df['PT_dpt'] = pt_dpt

# Save full dataset
df.to_csv(OUTPUT_FILE, index=False)
print(f"Saved: {OUTPUT_FILE}")
print(f"  {len(df):,} cells × {len(df.columns)} columns")

# Save summary statistics
summary = []
for ct in sorted(df['cell_type'].unique()):
    ct_df = df[df['cell_type'] == ct]
    summary.append({
        'cell_type': ct,
        'n_cells': len(ct_df),
        'PT_dpt_mean': ct_df['PT_dpt'].mean(),
        'PT_dpt_std': ct_df['PT_dpt'].std(),
        'PT_dpt_median': ct_df['PT_dpt'].median(),
        'stress_total_mean': ct_df['stress_total'].mean(),
        'PT_stress_corr': pearsonr(ct_df['PT_dpt'], ct_df['stress_total'])[0] if len(ct_df) > 10 else np.nan
    })

summary_df = pd.DataFrame(summary)
summary_df.to_csv(f'{OUTPUT_DIR}/PT_dpt_summary_by_celltype.csv', index=False)
print(f"Saved: {OUTPUT_DIR}/PT_dpt_summary_by_celltype.csv")

# ============================================================================
# Step 10: Visualizations
# ============================================================================
print("\n[10] Creating visualizations...")

# Figure 1: UMAP colored by PT_dpt
fig, axes = plt.subplots(1, 3, figsize=(18, 5))

sc.pl.umap(adata, color='dpt_pseudotime', cmap='viridis',
           title='PT_dpt (Diffusion Pseudotime)', ax=axes[0], show=False)

sc.pl.umap(adata, color='stress_total', cmap='Reds',
           title='Stress Total', ax=axes[1], show=False)

# Cell type (simplified labels)
adata.obs['cell_group'] = adata.obs['cell_type'].apply(lambda x:
    'Upper' if 'CUX2.RASGRF2' in x else
    'VAT1L' if 'VAT1L' in x else
    'Glia' if 'Glia' in x else
    'Vascular' if 'Vasc' in x else
    'Other'
)
sc.pl.umap(adata, color='cell_group', title='Cell Groups', ax=axes[2], show=False)

plt.tight_layout()
plt.savefig(f'{OUTPUT_DIR}/Fig1_PT_dpt_UMAP.png', dpi=300, bbox_inches='tight')
print(f"Saved: {OUTPUT_DIR}/Fig1_PT_dpt_UMAP.png")
plt.close()

# Figure 2: PT_dpt distribution by cell group
fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# Violin plot
ax = axes[0, 0]
cell_groups = ['Upper', 'Glia', 'VAT1L', 'Vascular']
plot_data = []
for cg in cell_groups:
    cg_pt = df[df['cell_type'].str.contains(
        'CUX2.RASGRF2' if cg == 'Upper' else
        'VAT1L' if cg == 'VAT1L' else
        'Glia' if cg == 'Glia' else
        'Vasc', na=False
    )]['PT_dpt']
    plot_data.append(cg_pt)

ax.violinplot(plot_data, positions=range(len(cell_groups)), showmeans=True)
ax.set_xticks(range(len(cell_groups)))
ax.set_xticklabels(cell_groups, rotation=45)
ax.set_ylabel('PT_dpt')
ax.set_title('PT_dpt Distribution by Cell Group')
ax.grid(axis='y', alpha=0.3)

# Vascular subtypes
ax = axes[0, 1]
vasc_types = sorted(df[df['cell_type'].str.contains('Vasc', na=False)]['cell_type'].unique())
vasc_data = [df[df['cell_type'] == vt]['PT_dpt'].values for vt in vasc_types]
ax.violinplot(vasc_data, positions=range(len(vasc_types)), showmeans=True)
ax.set_xticks(range(len(vasc_types)))
ax.set_xticklabels([vt.replace('Vasc.', '') for vt in vasc_types], rotation=45, ha='right')
ax.set_ylabel('PT_dpt')
ax.set_title('PT_dpt Distribution - Vascular Subtypes')
ax.grid(axis='y', alpha=0.3)

# PT_dpt vs stress scatter (key celltypes)
ax = axes[1, 0]
key_ct_labels = {
    'Ex.L2.L3.CUX2.RASGRF2': 'Upper',
    'Glia.Oligo': 'Oligo',
    'Vasc.Endo.Capillary': 'Endo.Cap',
    'Vasc.Mural.Pericyte': 'Pericyte'
}
for ct, label in key_ct_labels.items():
    ct_df = df[df['cell_type'] == ct]
    if len(ct_df) > 0:
        ax.scatter(ct_df['PT_dpt'], ct_df['stress_total'],
                   alpha=0.3, s=5, label=label)

ax.set_xlabel('PT_dpt')
ax.set_ylabel('Stress Total')
ax.set_title('PT_dpt vs Stress (Key Cell Types)')
ax.legend()
ax.grid(alpha=0.3)

# Correlation summary
ax = axes[1, 1]
corr_summary = summary_df[['cell_type', 'PT_stress_corr']].dropna()
corr_summary = corr_summary.sort_values('PT_stress_corr')
y_pos = np.arange(min(15, len(corr_summary)))
if len(corr_summary) > 15:
    # Show top 7 and bottom 7
    top7 = corr_summary.tail(7)
    bottom7 = corr_summary.head(7)
    plot_df = pd.concat([bottom7, top7])
else:
    plot_df = corr_summary

ax.barh(range(len(plot_df)), plot_df['PT_stress_corr'], alpha=0.7)
ax.set_yticks(range(len(plot_df)))
ax.set_yticklabels([ct.replace('Ex.', '').replace('Glia.', '').replace('Vasc.', '')
                     for ct in plot_df['cell_type']], fontsize=8)
ax.set_xlabel('Pearson r (PT_dpt vs stress)')
ax.set_title('PT_dpt-Stress Correlation by Cell Type')
ax.axvline(x=0, color='red', linestyle='--', alpha=0.5)
ax.grid(axis='x', alpha=0.3)

plt.tight_layout()
plt.savefig(f'{OUTPUT_DIR}/Fig2_PT_dpt_distributions.png', dpi=300, bbox_inches='tight')
print(f"Saved: {OUTPUT_DIR}/Fig2_PT_dpt_distributions.png")
plt.close()

print("\n" + "="*80)
print("Phase 9a Complete!")
print("="*80)
print(f"\nKey outputs:")
print(f"  - {OUTPUT_FILE}")
print(f"  - {OUTPUT_DIR}/PT_dpt_summary_by_celltype.csv")
print(f"  - {OUTPUT_DIR}/Fig1_PT_dpt_UMAP.png")
print(f"  - {OUTPUT_DIR}/Fig2_PT_dpt_distributions.png")
print(f"\nNext: Phase 9b - Define vascular-specific modules")
