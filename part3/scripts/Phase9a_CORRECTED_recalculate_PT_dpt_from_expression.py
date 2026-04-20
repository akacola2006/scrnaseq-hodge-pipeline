#!/usr/bin/env python3
"""
Phase 9a CORRECTED: Recalculate PT_dpt from Gene Expression Data

Purpose:
  - Load original gene expression data for ALL cell types
  - Calculate PT_dpt using high-variable genes (HVG)
  - Produce accurate, stress-independent pseudotime

Author: Claude Code
Date: 2025-11-24
"""

import pandas as pd
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import glob
import gzip
import warnings
warnings.filterwarnings('ignore')

# Configuration
EXPRESSION_DIR = '/home/akaco/als/motor_cortex_analysis'
OUTPUT_DIR = 'results/phase9_vascular_corrected'
Path(OUTPUT_DIR).mkdir(parents=True, exist_ok=True)

print("="*80)
print("Phase 9a CORRECTED: PT_dpt from Gene Expression Data")
print("="*80)

# ============================================================================
# Step 1: Find all expression files
# ============================================================================
print("\n[1] Finding expression files...")

# Find both .csv and .csv.gz files
csv_files = glob.glob(f'{EXPRESSION_DIR}/motor_cortex_*_expression.csv')
gz_files = glob.glob(f'{EXPRESSION_DIR}/motor_cortex_*_expression.csv.gz')

all_files = csv_files + gz_files
print(f"Found {len(all_files)} expression files")
print(f"  .csv files: {len(csv_files)}")
print(f"  .csv.gz files: {len(gz_files)}")

# ============================================================================
# Step 2: Load expression data incrementally
# ============================================================================
print("\n[2] Loading expression data...")

adatas = []
total_cells = 0

for idx, filepath in enumerate(sorted(all_files)):
    cell_type = Path(filepath).stem.replace('motor_cortex_', '').replace('_expression', '')

    print(f"\n  [{idx+1}/{len(all_files)}] Loading {cell_type}...")

    try:
        # Load expression matrix
        if filepath.endswith('.gz'):
            df = pd.read_csv(filepath, compression='gzip', index_col=0)
        else:
            df = pd.read_csv(filepath, index_col=0)

        print(f"    Shape: {df.shape[0]:,} cells × {df.shape[1]:,} genes")

        # Create AnnData
        adata = sc.AnnData(X=df.values, obs=pd.DataFrame(index=df.index))
        adata.var_names = df.columns
        adata.obs['cell_type'] = cell_type

        adatas.append(adata)
        total_cells += adata.n_obs

        print(f"    Loaded: {adata.n_obs:,} cells")

    except Exception as e:
        print(f"    ERROR loading {cell_type}: {e}")
        continue

print(f"\nTotal cells loaded: {total_cells:,}")
print(f"Successfully loaded {len(adatas)} cell types")

# ============================================================================
# Step 3: Concatenate all cell types
# ============================================================================
print("\n[3] Concatenating all cell types...")

adata_combined = sc.concat(adatas, join='outer', fill_value=0)
print(f"Combined AnnData: {adata_combined.n_obs:,} cells × {adata_combined.n_vars:,} genes")

# Save cell type mapping
adata_combined.obs.to_csv(f'{OUTPUT_DIR}/cell_metadata.csv')
print(f"Saved: {OUTPUT_DIR}/cell_metadata.csv")

# ============================================================================
# Step 4: Load existing module scores and metadata
# ============================================================================
print("\n[4] Loading existing module scores and metadata...")

# Load the module scores and metadata from Phase 9a
existing_data = pd.read_csv('results/phase9_vascular/cell_level_features_ALL_with_PTdpt.csv')
print(f"Loaded existing data: {len(existing_data):,} cells")

# Match cells by cell_id
cell_ids = adata_combined.obs_names
matched_data = existing_data.set_index('cell_id').reindex(cell_ids)

# Add metadata to adata
adata_combined.obs['condition'] = matched_data['condition'].values
adata_combined.obs['stress_total'] = matched_data['stress_total'].values
adata_combined.obs['PT_IDS'] = matched_data['PT_IDS'].values

if 'patient_id' in matched_data.columns:
    adata_combined.obs['patient_id'] = matched_data['patient_id'].values

# Add module scores
module_cols = [c for c in matched_data.columns if c.startswith('module_')]
for col in module_cols:
    adata_combined.obs[col] = matched_data[col].values

print(f"Added {len(module_cols)} module scores and metadata")

# ============================================================================
# Step 5: Preprocessing
# ============================================================================
print("\n[5] Preprocessing...")

# Basic filtering
print("  Filtering genes...")
sc.pp.filter_genes(adata_combined, min_cells=10)
print(f"  After filtering: {adata_combined.n_vars:,} genes")

# Normalize
print("  Normalizing...")
sc.pp.normalize_total(adata_combined, target_sum=1e4)
sc.pp.log1p(adata_combined)

# Identify highly variable genes
print("  Identifying highly variable genes...")
sc.pp.highly_variable_genes(adata_combined, n_top_genes=3000, flavor='seurat_v3')
n_hvg = adata_combined.var['highly_variable'].sum()
print(f"  Selected {n_hvg} highly variable genes")

# ============================================================================
# Step 6: PCA on HVG
# ============================================================================
print("\n[6] Running PCA...")
sc.pp.scale(adata_combined, max_value=10)
sc.tl.pca(adata_combined, n_comps=50, use_highly_variable=True, svd_solver='arpack')

var_explained = adata_combined.uns['pca']['variance_ratio'][:10].sum()
print(f"  PCA variance explained (first 10 PCs): {var_explained:.3f}")

# ============================================================================
# Step 7: Compute neighbor graph
# ============================================================================
print("\n[7] Computing neighbor graph...")
sc.pp.neighbors(adata_combined, n_neighbors=30, n_pcs=50, method='umap')

# UMAP for visualization
print("  Computing UMAP...")
sc.tl.umap(adata_combined)

# ============================================================================
# Step 8: Diffusion map and PT_dpt
# ============================================================================
print("\n[8] Computing diffusion map...")
sc.tl.diffmap(adata_combined, n_comps=15)

# Determine root cell
print("\n[9] Determining root cell for PT_dpt...")

homeostatic_types = ['Glia.Oligo', 'Glia.Astro.GFAP-neg', 'Glia.OPC',
                     'Vasc.Endo.Capillary', 'Vasc.Fibro.CLMP.PDGFRA']
homeostatic_mask = adata_combined.obs['cell_type'].isin(homeostatic_types)

if homeostatic_mask.sum() > 0:
    homeostatic_cells = adata_combined.obs[homeostatic_mask]
    root_idx = homeostatic_cells['stress_total'].idxmin()
    root_cell_idx = adata_combined.obs_names.get_loc(root_idx)
    root_cell_type = adata_combined.obs.loc[root_idx, 'cell_type']
    root_stress = adata_combined.obs.loc[root_idx, 'stress_total']
    print(f"  Root cell: {root_idx}")
    print(f"  Cell type: {root_cell_type}")
    print(f"  Stress: {root_stress:.4f}")
else:
    root_idx = adata_combined.obs['stress_total'].idxmin()
    root_cell_idx = adata_combined.obs_names.get_loc(root_idx)
    root_cell_type = adata_combined.obs.loc[root_idx, 'cell_type']
    root_stress = adata_combined.obs.loc[root_idx, 'stress_total']
    print(f"  Root cell (fallback): {root_idx}")

# Compute DPT
print("\n[10] Computing diffusion pseudotime (PT_dpt)...")
adata_combined.uns['iroot'] = root_cell_idx
sc.tl.dpt(adata_combined, n_dcs=10)

pt_dpt = adata_combined.obs['dpt_pseudotime'].values
print(f"  PT_dpt range: {pt_dpt.min():.4f} - {pt_dpt.max():.4f}")
print(f"  PT_dpt mean: {pt_dpt.mean():.4f}")
print(f"  PT_dpt median: {np.median(pt_dpt):.4f}")

# ============================================================================
# Step 11: Validate PT_dpt vs stress independence
# ============================================================================
print("\n[11] Validating PT_dpt independence from stress...")

from scipy.stats import pearsonr, spearmanr

r_pearson, p_pearson = pearsonr(pt_dpt, adata_combined.obs['stress_total'])
r_spearman, p_spearman = spearmanr(pt_dpt, adata_combined.obs['stress_total'])
r2 = r_pearson ** 2

print(f"  Overall PT_dpt vs stress_total:")
print(f"    Pearson r: {r_pearson:.4f} (R² = {r2:.4f}, p = {p_pearson:.2e})")
print(f"    Spearman ρ: {r_spearman:.4f} (p = {p_spearman:.2e})")

# Per cell type
print(f"\n  PT_dpt vs stress_total by key cell types:")
key_celltypes = ['Ex.L2.L3.CUX2.RASGRF2', 'Glia.Oligo', 'Glia.Astro.GFAP-neg',
                 'Ex.L5.VAT1L.EYA4', 'Ex.L5.VAT1L.THSD4',
                 'Vasc.Endo.Arterial', 'Vasc.Endo.Capillary', 'Vasc.Mural.Pericyte']

for ct in key_celltypes:
    ct_mask = adata_combined.obs['cell_type'] == ct
    if ct_mask.sum() > 10:
        ct_pt = pt_dpt[ct_mask]
        ct_stress = adata_combined.obs.loc[ct_mask, 'stress_total']
        r, p = pearsonr(ct_pt, ct_stress)
        r2_ct = r ** 2
        print(f"    {ct:35s} r={r:6.3f} (R²={r2_ct:.3f}, n={ct_mask.sum():5d})")

# ============================================================================
# Step 12: Save results
# ============================================================================
print("\n[12] Saving results...")

# Create output dataframe
output_df = adata_combined.obs.copy()
output_df['PT_dpt_corrected'] = pt_dpt
output_df['cell_id'] = output_df.index

# Reorder columns
cols = ['cell_id', 'cell_type', 'condition', 'PT_dpt_corrected', 'PT_IDS', 'stress_total']
module_cols = [c for c in output_df.columns if c.startswith('module_')]
if 'patient_id' in output_df.columns:
    cols.insert(3, 'patient_id')
cols.extend(module_cols)

output_df = output_df[cols]
output_df.to_csv(f'{OUTPUT_DIR}/cell_level_features_ALL_with_PTdpt_CORRECTED.csv', index=False)
print(f"Saved: {OUTPUT_DIR}/cell_level_features_ALL_with_PTdpt_CORRECTED.csv")
print(f"  {len(output_df):,} cells × {len(output_df.columns)} columns")

# Save summary statistics
summary = []
for ct in sorted(adata_combined.obs['cell_type'].unique()):
    ct_df = output_df[output_df['cell_type'] == ct]
    if len(ct_df) > 0:
        summary.append({
            'cell_type': ct,
            'n_cells': len(ct_df),
            'PT_dpt_mean': ct_df['PT_dpt_corrected'].mean(),
            'PT_dpt_std': ct_df['PT_dpt_corrected'].std(),
            'PT_dpt_median': ct_df['PT_dpt_corrected'].median(),
            'stress_total_mean': ct_df['stress_total'].mean(),
            'PT_stress_corr': pearsonr(ct_df['PT_dpt_corrected'], ct_df['stress_total'])[0] if len(ct_df) > 10 else np.nan
        })

summary_df = pd.DataFrame(summary)
summary_df.to_csv(f'{OUTPUT_DIR}/PT_dpt_summary_by_celltype_CORRECTED.csv', index=False)
print(f"Saved: {OUTPUT_DIR}/PT_dpt_summary_by_celltype_CORRECTED.csv")

# Save AnnData object (compressed)
print("\n[13] Saving AnnData object...")
adata_combined.write_h5ad(f'{OUTPUT_DIR}/adata_all_cells_with_PTdpt.h5ad', compression='gzip')
print(f"Saved: {OUTPUT_DIR}/adata_all_cells_with_PTdpt.h5ad")

# ============================================================================
# Step 14: Visualizations
# ============================================================================
print("\n[14] Creating visualizations...")

# Figure 1: UMAP
fig, axes = plt.subplots(1, 3, figsize=(18, 5))

sc.pl.umap(adata_combined, color='dpt_pseudotime', cmap='viridis',
           title='PT_dpt (Corrected from Expression)', ax=axes[0], show=False)

sc.pl.umap(adata_combined, color='stress_total', cmap='Reds',
           title='Stress Total', ax=axes[1], show=False)

# Cell groups
adata_combined.obs['cell_group'] = adata_combined.obs['cell_type'].apply(lambda x:
    'Upper' if 'CUX2.RASGRF2' in x else
    'VAT1L' if 'VAT1L' in x else
    'Glia' if 'Glia' in x else
    'Vascular' if 'Vasc' in x else
    'Other'
)
sc.pl.umap(adata_combined, color='cell_group', title='Cell Groups', ax=axes[2], show=False)

plt.tight_layout()
plt.savefig(f'{OUTPUT_DIR}/Fig1_PT_dpt_CORRECTED_UMAP.png', dpi=300, bbox_inches='tight')
print(f"Saved: {OUTPUT_DIR}/Fig1_PT_dpt_CORRECTED_UMAP.png")
plt.close()

# Figure 2: PT_dpt distributions
fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# Violin plot
ax = axes[0, 0]
cell_groups = ['Upper', 'Glia', 'VAT1L', 'Vascular']
plot_data = []
for cg in cell_groups:
    if cg == 'Upper':
        mask = output_df['cell_type'].str.contains('CUX2.RASGRF2', na=False)
    elif cg == 'VAT1L':
        mask = output_df['cell_type'].str.contains('VAT1L', na=False)
    elif cg == 'Glia':
        mask = output_df['cell_type'].str.contains('Glia', na=False)
    elif cg == 'Vascular':
        mask = output_df['cell_type'].str.contains('Vasc', na=False)

    cg_pt = output_df[mask]['PT_dpt_corrected']
    plot_data.append(cg_pt)

ax.violinplot(plot_data, positions=range(len(cell_groups)), showmeans=True)
ax.set_xticks(range(len(cell_groups)))
ax.set_xticklabels(cell_groups, rotation=45)
ax.set_ylabel('PT_dpt (Corrected)')
ax.set_title('PT_dpt Distribution by Cell Group')
ax.grid(axis='y', alpha=0.3)

# Vascular subtypes
ax = axes[0, 1]
vasc_types = sorted(output_df[output_df['cell_type'].str.contains('Vasc', na=False)]['cell_type'].unique())
vasc_data = [output_df[output_df['cell_type'] == vt]['PT_dpt_corrected'].values for vt in vasc_types]
if len(vasc_data) > 0:
    ax.violinplot(vasc_data, positions=range(len(vasc_types)), showmeans=True)
    ax.set_xticks(range(len(vasc_types)))
    ax.set_xticklabels([vt.replace('Vasc.', '') for vt in vasc_types], rotation=45, ha='right')
    ax.set_ylabel('PT_dpt (Corrected)')
    ax.set_title('PT_dpt - Vascular Subtypes')
    ax.grid(axis='y', alpha=0.3)

# PT vs stress scatter
ax = axes[1, 0]
for ct, label in [('Ex.L2.L3.CUX2.RASGRF2', 'Upper'),
                  ('Glia.Oligo', 'Oligo'),
                  ('Vasc.Endo.Capillary', 'Endo.Cap'),
                  ('Vasc.Mural.Pericyte', 'Pericyte')]:
    ct_df = output_df[output_df['cell_type'] == ct]
    if len(ct_df) > 0:
        ax.scatter(ct_df['PT_dpt_corrected'], ct_df['stress_total'],
                   alpha=0.3, s=5, label=label)

ax.set_xlabel('PT_dpt (Corrected)')
ax.set_ylabel('Stress Total')
ax.set_title('PT_dpt vs Stress (Key Cell Types)')
ax.legend()
ax.grid(alpha=0.3)

# Correlation summary
ax = axes[1, 1]
corr_summary = summary_df[['cell_type', 'PT_stress_corr']].dropna()
corr_summary = corr_summary.sort_values('PT_stress_corr')
if len(corr_summary) > 15:
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
plt.savefig(f'{OUTPUT_DIR}/Fig2_PT_dpt_CORRECTED_distributions.png', dpi=300, bbox_inches='tight')
print(f"Saved: {OUTPUT_DIR}/Fig2_PT_dpt_CORRECTED_distributions.png")
plt.close()

print("\n" + "="*80)
print("Phase 9a CORRECTED Complete!")
print("="*80)
print(f"\nKey outputs:")
print(f"  - {OUTPUT_DIR}/cell_level_features_ALL_with_PTdpt_CORRECTED.csv")
print(f"  - {OUTPUT_DIR}/PT_dpt_summary_by_celltype_CORRECTED.csv")
print(f"  - {OUTPUT_DIR}/adata_all_cells_with_PTdpt.h5ad")
print(f"  - {OUTPUT_DIR}/Fig1_PT_dpt_CORRECTED_UMAP.png")
print(f"  - {OUTPUT_DIR}/Fig2_PT_dpt_CORRECTED_distributions.png")
print(f"\nNext: Re-run Phase 9c+9d with corrected PT_dpt")
