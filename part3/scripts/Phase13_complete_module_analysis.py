#!/usr/bin/env python3
"""
Phase 13 Complete Module Analysis
==================================
Comprehensive analysis of all 23 Phase 13 functional modules across NVU nodes

Modules analyzed:
1. Metabolic Axis: Mitochondria, ER_Stress, Protein_Homeostasis, Metabolism
2. Hyperexcitability Axis: Calcium_Signaling, Ion_Transport
3. Inflammatory Axis: Oxidative_Stress, Inflammation, Complement
4. Structural/ECM: ECM, Cytoskeleton, Myelination
5. Vascular: Angiogenesis
6. Cellular Processes: Apoptosis, Autophagy, Cell_Cycle, DNA_Repair, Epigenetic
7. Signaling: Growth_Factors, Synaptic
8. Transcriptional: Transcription, RNA_Processing, lncRNA

This script:
- Computes module φ (energy) for all 23 modules across 5 NVU nodes
- Identifies which axis each module belongs to
- Analyzes temporal dynamics (correlation with PT_dpt)
- Compares Control vs ALS for each module
- Integrates with PT_STMN2 analysis
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

# Configuration
OUTPUT_DIR = Path("results/phase13_complete_modules")
OUTPUT_DIR.mkdir(exist_ok=True, parents=True)

# NVU nodes (5-node framework from Phase 9)
# Updated with actual cell type names (using dots, not underscores)
NVU_NODES = {
    'Vascular': ['Vasc.Endo', 'Vasc.Fibro', 'Vasc.Mural'],
    'Glia': ['Glia.Astro', 'Glia.Oligo', 'Glia.OPC', 'Glia.Micro'],
    'L2/L3': ['Ex.L2.L3', 'Ex.L3.L5.CUX2'],
    'L3-L6': ['Ex.L4', 'Ex.L5.L6', 'Ex.L6', 'Ex.L5.PCP4'],
    'VAT1L': ['Ex.L5.VAT1L']
}

# Module classification by axis (from Phase 13)
MODULE_AXES = {
    'Metabolic': ['Mitochondria', 'ER_Stress', 'Protein_Homeostasis', 'Metabolism'],
    'Hyperexcitability': ['Calcium_Signaling', 'Ion_Transport'],
    'Inflammatory': ['Oxidative_Stress', 'Inflammation', 'Complement'],
    'Structural_ECM': ['ECM', 'Cytoskeleton', 'Myelination'],
    'Vascular': ['Angiogenesis'],
    'Cellular_Processes': ['Apoptosis', 'Autophagy', 'Cell_Cycle', 'DNA_Repair', 'Epigenetic'],
    'Signaling': ['Growth_Factors', 'Synaptic'],
    'Transcriptional': ['Transcription', 'RNA_Processing', 'lncRNA']
}

# All 23 modules
ALL_MODULES = [
    'Angiogenesis', 'Apoptosis', 'Autophagy', 'Calcium_Signaling', 'Cell_Cycle',
    'Complement', 'Cytoskeleton', 'DNA_Repair', 'ECM', 'ER_Stress', 'Epigenetic',
    'Growth_Factors', 'Inflammation', 'Ion_Transport', 'Metabolism', 'Mitochondria',
    'Myelination', 'Oxidative_Stress', 'Protein_Homeostasis', 'RNA_Processing',
    'Synaptic', 'Transcription', 'lncRNA'
]

print("=" * 80)
print("PHASE 13: COMPLETE MODULE ANALYSIS")
print("=" * 80)
print(f"\nAnalyzing {len(ALL_MODULES)} functional modules across 5 NVU nodes")
print(f"Output directory: {OUTPUT_DIR}")

# Load Phase 9 data with all modules and PT_dpt
print("\n" + "=" * 80)
print("STEP 1: Loading cell-level features with all modules")
print("=" * 80)

data_path = "results/phase9_vascular/cell_level_features_ALL_with_PTdpt.csv"
print(f"Loading: {data_path}")

df = pd.read_csv(data_path)
print(f"Loaded {len(df):,} cells")

# Check available modules
available_modules = [col.replace('module_', '') for col in df.columns if col.startswith('module_')]
print(f"\nAvailable modules: {len(available_modules)}")

missing_modules = set(ALL_MODULES) - set(available_modules)
if missing_modules:
    print(f"WARNING: Missing modules: {missing_modules}")

# Map cell types to NVU nodes
def map_to_nvu(cell_type):
    for nvu_node, cell_types in NVU_NODES.items():
        if any(ct in cell_type for ct in cell_types):
            return nvu_node
    return 'Other'

df['NVU_node'] = df['cell_type'].apply(map_to_nvu)

print("\nNVU node distribution:")
print(df['NVU_node'].value_counts())

# Condition distribution
print("\nCondition distribution:")
print(df['condition'].value_counts())

# Compute module φ (energy) for each module across NVU nodes
print("\n" + "=" * 80)
print("STEP 2: Computing module φ (energy) across NVU nodes")
print("=" * 80)

nvu_module_summary = []

for module in ALL_MODULES:
    module_col = f'module_{module}'

    if module_col not in df.columns:
        print(f"WARNING: {module_col} not found, skipping...")
        continue

    # Compute mean φ for each NVU node and condition
    for nvu_node in ['Vascular', 'Glia', 'L2/L3', 'L3-L6', 'VAT1L']:
        subset = df[df['NVU_node'] == nvu_node]

        if len(subset) == 0:
            continue

        for condition in ['Control', 'ALS']:
            cond_subset = subset[subset['condition'] == condition]

            if len(cond_subset) == 0:
                continue

            phi_mean = cond_subset[module_col].mean()
            phi_std = cond_subset[module_col].std()
            n_cells = len(cond_subset)

            nvu_module_summary.append({
                'Module': module,
                'NVU_node': nvu_node,
                'Condition': condition,
                'phi_mean': phi_mean,
                'phi_std': phi_std,
                'n_cells': n_cells
            })

# Convert to DataFrame
df_summary = pd.DataFrame(nvu_module_summary)
print(f"\nComputed φ for {len(df_summary)} module-NVU-condition combinations")

# Save summary
df_summary.to_csv(OUTPUT_DIR / "module_phi_NVU_summary.csv", index=False)
print(f"Saved: {OUTPUT_DIR / 'module_phi_NVU_summary.csv'}")

# Compute Control vs ALS differences
print("\n" + "=" * 80)
print("STEP 3: Computing Control vs ALS differences for each module")
print("=" * 80)

module_diffs = []

for module in ALL_MODULES:
    module_col = f'module_{module}'

    if module_col not in df.columns:
        continue

    for nvu_node in ['Vascular', 'Glia', 'L2/L3', 'L3-L6', 'VAT1L']:
        subset = df[df['NVU_node'] == nvu_node]

        if len(subset) == 0:
            continue

        control = subset[subset['condition'] == 'Control'][module_col].dropna()
        als = subset[subset['condition'] == 'ALS'][module_col].dropna()

        if len(control) > 10 and len(als) > 10:
            # T-test
            t_stat, p_val = stats.ttest_ind(control, als)

            # Effect size (Cohen's d)
            cohens_d = (als.mean() - control.mean()) / np.sqrt((als.var() + control.var()) / 2)

            module_diffs.append({
                'Module': module,
                'NVU_node': nvu_node,
                'Control_mean': control.mean(),
                'ALS_mean': als.mean(),
                'Difference': als.mean() - control.mean(),
                'Fold_change': als.mean() / control.mean() if control.mean() > 0 else np.nan,
                't_stat': t_stat,
                'p_value': p_val,
                'cohens_d': cohens_d,
                'n_control': len(control),
                'n_als': len(als)
            })

df_diffs = pd.DataFrame(module_diffs)
print(f"\nComputed differences for {len(df_diffs)} module-NVU combinations")

# FDR correction
from statsmodels.stats.multitest import multipletests
df_diffs['p_adj_fdr'] = multipletests(df_diffs['p_value'], method='fdr_bh')[1]

# Save differences
df_diffs.to_csv(OUTPUT_DIR / "module_Control_vs_ALS_differences.csv", index=False)
print(f"Saved: {OUTPUT_DIR / 'module_Control_vs_ALS_differences.csv'}")

# Identify significantly different modules
sig_modules = df_diffs[df_diffs['p_adj_fdr'] < 0.05].sort_values('p_adj_fdr')
print(f"\nSignificantly different modules (FDR < 0.05): {len(sig_modules)}")
print("\nTop 10 most significant:")
print(sig_modules[['Module', 'NVU_node', 'Difference', 'cohens_d', 'p_adj_fdr']].head(10))

# Temporal analysis: Correlation with PT_dpt
print("\n" + "=" * 80)
print("STEP 4: Temporal analysis - correlation with PT_dpt")
print("=" * 80)

module_pt_corr = []

for module in ALL_MODULES:
    module_col = f'module_{module}'

    if module_col not in df.columns:
        continue

    for nvu_node in ['Vascular', 'Glia', 'L2/L3', 'L3-L6', 'VAT1L']:
        subset = df[df['NVU_node'] == nvu_node].dropna(subset=[module_col, 'PT_dpt'])

        if len(subset) < 20:
            continue

        r, p = stats.pearsonr(subset[module_col], subset['PT_dpt'])

        module_pt_corr.append({
            'Module': module,
            'NVU_node': nvu_node,
            'r_with_PT_dpt': r,
            'p_value': p,
            'n_cells': len(subset)
        })

df_pt_corr = pd.DataFrame(module_pt_corr)

# FDR correction
df_pt_corr['p_adj_fdr'] = multipletests(df_pt_corr['p_value'], method='fdr_bh')[1]

# Save correlations
df_pt_corr.to_csv(OUTPUT_DIR / "module_PT_dpt_correlations.csv", index=False)
print(f"Saved: {OUTPUT_DIR / 'module_PT_dpt_correlations.csv'}")

# Identify modules that progress with disease
sig_corr = df_pt_corr[df_pt_corr['p_adj_fdr'] < 0.05].sort_values('r_with_PT_dpt', key=abs, ascending=False)
print(f"\nModules significantly correlated with PT_dpt (FDR < 0.05): {len(sig_corr)}")
print("\nTop 10 strongest correlations:")
print(sig_corr[['Module', 'NVU_node', 'r_with_PT_dpt', 'p_adj_fdr']].head(10))

# Visualizations
print("\n" + "=" * 80)
print("STEP 5: Creating visualizations")
print("=" * 80)

# 5.1: Heatmap of module φ across NVU nodes (ALS only)
print("\n[1/4] Creating heatmap of module φ across NVU nodes...")

pivot_data = df_summary[df_summary['Condition'] == 'ALS'].pivot(
    index='Module', columns='NVU_node', values='phi_mean'
)

# Reorder columns
col_order = ['Vascular', 'Glia', 'L2/L3', 'L3-L6', 'VAT1L']
pivot_data = pivot_data[[col for col in col_order if col in pivot_data.columns]]

fig, ax = plt.subplots(figsize=(12, 16))
sns.heatmap(pivot_data, annot=True, fmt='.2f', cmap='YlOrRd', ax=ax,
            cbar_kws={'label': 'Module φ (Energy)'}, linewidths=0.5)
ax.set_title('Module φ (Energy) Across NVU Nodes (ALS Condition)', fontsize=16, fontweight='bold', pad=20)
ax.set_xlabel('NVU Node', fontsize=14)
ax.set_ylabel('Functional Module', fontsize=14)
plt.tight_layout()
plt.savefig(OUTPUT_DIR / "Fig1_module_phi_heatmap_ALS.png", dpi=300, bbox_inches='tight')
print(f"Saved: {OUTPUT_DIR / 'Fig1_module_phi_heatmap_ALS.png'}")
plt.close()

# 5.2: Heatmap of ALS vs Control fold change
print("\n[2/4] Creating heatmap of ALS vs Control fold change...")

pivot_fc = df_diffs.pivot(index='Module', columns='NVU_node', values='Fold_change')
pivot_fc = pivot_fc[[col for col in col_order if col in pivot_fc.columns]]

fig, ax = plt.subplots(figsize=(12, 16))
sns.heatmap(pivot_fc, annot=True, fmt='.2f', cmap='RdBu_r', center=1.0,
            vmin=0.5, vmax=1.5, ax=ax, cbar_kws={'label': 'ALS/Control Fold Change'},
            linewidths=0.5)
ax.set_title('Module Dysregulation: ALS vs Control Fold Change Across NVU', fontsize=16, fontweight='bold', pad=20)
ax.set_xlabel('NVU Node', fontsize=14)
ax.set_ylabel('Functional Module', fontsize=14)
plt.tight_layout()
plt.savefig(OUTPUT_DIR / "Fig2_module_fold_change_heatmap.png", dpi=300, bbox_inches='tight')
print(f"Saved: {OUTPUT_DIR / 'Fig2_module_fold_change_heatmap.png'}")
plt.close()

# 5.3: Heatmap of correlation with PT_dpt
print("\n[3/4] Creating heatmap of PT_dpt correlations...")

pivot_corr = df_pt_corr.pivot(index='Module', columns='NVU_node', values='r_with_PT_dpt')
pivot_corr = pivot_corr[[col for col in col_order if col in pivot_corr.columns]]

fig, ax = plt.subplots(figsize=(12, 16))
sns.heatmap(pivot_corr, annot=True, fmt='.3f', cmap='RdBu_r', center=0,
            vmin=-0.6, vmax=0.6, ax=ax, cbar_kws={'label': 'Pearson r with PT_dpt'},
            linewidths=0.5)
ax.set_title('Module Temporal Dynamics: Correlation with Disease Progression (PT_dpt)', fontsize=16, fontweight='bold', pad=20)
ax.set_xlabel('NVU Node', fontsize=14)
ax.set_ylabel('Functional Module', fontsize=14)
plt.tight_layout()
plt.savefig(OUTPUT_DIR / "Fig3_module_PT_dpt_correlation_heatmap.png", dpi=300, bbox_inches='tight')
print(f"Saved: {OUTPUT_DIR / 'Fig3_module_PT_dpt_correlation_heatmap.png'}")
plt.close()

# 5.4: Bar plot showing modules by axis
print("\n[4/4] Creating axis classification bar plot...")

# Assign axis to each module
module_axis_map = {}
for axis, modules in MODULE_AXES.items():
    for module in modules:
        module_axis_map[module] = axis

df_diffs['Axis'] = df_diffs['Module'].map(module_axis_map)

# Compute mean fold change by axis and NVU node
axis_summary = df_diffs.groupby(['Axis', 'NVU_node'])['Fold_change'].mean().reset_index()

fig, ax = plt.subplots(figsize=(14, 8))
pivot_axis = axis_summary.pivot(index='Axis', columns='NVU_node', values='Fold_change')
pivot_axis = pivot_axis[[col for col in col_order if col in pivot_axis.columns]]
pivot_axis.plot(kind='bar', ax=ax, width=0.8)
ax.set_title('Module Dysregulation by Functional Axis Across NVU', fontsize=16, fontweight='bold')
ax.set_ylabel('ALS/Control Fold Change', fontsize=14)
ax.set_xlabel('Functional Axis', fontsize=14)
ax.axhline(y=1.0, color='black', linestyle='--', linewidth=1, alpha=0.5)
ax.legend(title='NVU Node', bbox_to_anchor=(1.05, 1), loc='upper left')
ax.tick_params(axis='x', rotation=45)
plt.tight_layout()
plt.savefig(OUTPUT_DIR / "Fig4_axis_dysregulation_barplot.png", dpi=300, bbox_inches='tight')
print(f"Saved: {OUTPUT_DIR / 'Fig4_axis_dysregulation_barplot.png'}")
plt.close()

print("\n" + "=" * 80)
print("PHASE 13 COMPLETE MODULE ANALYSIS FINISHED!")
print("=" * 80)
print(f"\nOutput files saved to: {OUTPUT_DIR}/")
print("\nKey findings:")
print(f"  - Analyzed {len(ALL_MODULES)} functional modules")
print(f"  - {len(sig_modules)} modules significantly dysregulated (FDR < 0.05)")
print(f"  - {len(sig_corr)} modules significantly correlated with disease progression")
print(f"  - Results organized by {len(MODULE_AXES)} functional axes")
