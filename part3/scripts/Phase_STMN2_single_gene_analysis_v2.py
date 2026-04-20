#!/usr/bin/env python3
"""
Phase STMN2 Single Gene Analysis (Version 2 - Fixed)
====================================================
Analyze STMN2 SINGLE GENE expression (not pathway) to create PT_STMN2_single

This script:
1. Extracts STMN2 (ENSG00000104435) expression from all cell-type expression files
2. Computes φ_STMN2_single = [(STMN2_expr - μ_Control) / σ_Control]²
3. Creates PT_STMN2_single for all cells
4. Compares PT_STMN2_single vs PT_STMN2_pathway vs PT_TDP
5. Integrates with Phase 13 axes framework

Key Difference:
- STMN2_pathway (9 genes): STMN1-4, MAPT, MAP2, MAP1B, DPYSL2-3
- STMN2_single (1 gene): STMN2 only (TDP-43's direct target)
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from scipy import stats
import glob
import warnings
warnings.filterwarnings('ignore')

# Configuration
STMN2_ENSEMBL_ID = "ENSG00000104435"
OUTPUT_DIR = Path("results/phaseTDP_STMN2_single")
OUTPUT_DIR.mkdir(exist_ok=True, parents=True)

# Find all expression files
EXPRESSION_DIR = Path("/home/akaco/als/motor_cortex_analysis")
expr_files = glob.glob(str(EXPRESSION_DIR / "*expression*.csv"))

print("=" * 80)
print("PHASE STMN2 SINGLE GENE ANALYSIS (v2)")
print("=" * 80)
print(f"\nSTMN2 Ensembl ID: {STMN2_ENSEMBL_ID}")
print(f"Output directory: {OUTPUT_DIR}")
print(f"\nFound {len(expr_files)} expression files")

# Step 1: Extract STMN2 expression from all cell types
print("\n" + "=" * 80)
print("STEP 1: Extracting STMN2 expression from all expression files")
print("=" * 80)

stmn2_data_list = []

for expr_file in sorted(expr_files):
    # Extract cell type from filename
    # Format: motor_cortex_<CellType>_expression.csv
    filename = Path(expr_file).name
    cell_type = filename.replace('motor_cortex_', '').replace('_expression.csv', '')

    print(f"\n[{cell_type}]")
    print(f"  File: {filename}")

    try:
        # Load expression file (genes as rows, cells as columns)
        expr_df = pd.read_csv(expr_file, index_col=0)
        print(f"  Shape: {expr_df.shape} (genes x cells)")

        # Check if STMN2 is present
        if STMN2_ENSEMBL_ID in expr_df.index:
            stmn2_expr = expr_df.loc[STMN2_ENSEMBL_ID, :]
            print(f"  ✓ STMN2 found! Expression range: [{stmn2_expr.min():.3f}, {stmn2_expr.max():.3f}]")
            print(f"    Mean: {stmn2_expr.mean():.3f}, Median: {stmn2_expr.median():.3f}")

            # Create dataframe with cell IDs and STMN2 expression
            cell_data = pd.DataFrame({
                'cell_id': stmn2_expr.index,
                'cell_type': cell_type,
                'STMN2_expr': stmn2_expr.values
            })
            stmn2_data_list.append(cell_data)
            print(f"    Added {len(cell_data)} cells")
        else:
            print(f"  ✗ STMN2 NOT FOUND in this cell type")

    except Exception as e:
        print(f"  ERROR: {e}")

# Combine all STMN2 data
if len(stmn2_data_list) == 0:
    raise ValueError("STMN2 not found in any cell type!")

stmn2_all = pd.concat(stmn2_data_list, ignore_index=True)
print(f"\n{'=' * 80}")
print(f"STMN2 EXTRACTION SUMMARY")
print(f"{'=' * 80}")
print(f"Total cells with STMN2 data: {len(stmn2_all):,}")
print(f"Cell types with STMN2: {stmn2_all['cell_type'].nunique()}")
print(f"\nCell type distribution:")
print(stmn2_all['cell_type'].value_counts())
print(f"\nSTMN2 expression statistics:")
print(stmn2_all['STMN2_expr'].describe())

# Step 2: Load metadata to get Control/ALS labels
print("\n" + "=" * 80)
print("STEP 2: Loading cell metadata (Control/ALS labels)")
print("=" * 80)

# Load metadata from Phase 9 data
metadata_path = "results/phase9_vascular/cell_level_features_ALL_with_PTdpt.csv"
print(f"Loading: {metadata_path}")

df_meta = pd.read_csv(metadata_path)
print(f"Loaded {len(df_meta):,} cells from metadata")
print(f"Columns in metadata: {list(df_meta.columns[:10])}...")

# Check condition distribution
print(f"\nCondition distribution in metadata:")
print(df_meta['condition'].value_counts())

# Step 3: Merge STMN2 expression with metadata
print("\n" + "=" * 80)
print("STEP 3: Merging STMN2 expression with metadata")
print("=" * 80)

# Merge on cell_id
print(f"\nMerging {len(stmn2_all)} STMN2 cells with {len(df_meta)} metadata cells...")
merged = pd.merge(
    stmn2_all,
    df_meta[['cell_id', 'cell_type', 'condition']],
    on='cell_id',
    how='left',
    suffixes=('_expr', '_meta')
)

print(f"Merged data: {len(merged):,} cells")
print(f"Cells with condition info: {merged['condition'].notna().sum():,}")
print(f"Cells without condition info: {merged['condition'].isna().sum():,}")

if merged['condition'].isna().sum() > 0:
    print("\nWARNING: Some cells don't have condition labels")
    print("Sample cells without condition:")
    print(merged[merged['condition'].isna()]['cell_id'].head())

# Use cell_type from expression file (more specific)
merged['cell_type_final'] = merged['cell_type_expr']

# Step 4: Compute φ_STMN2_single (Control-normalized energy)
print("\n" + "=" * 80)
print("STEP 4: Computing φ_STMN2_single (Control-normalized energy)")
print("=" * 80)

# Remove cells without condition labels
merged_clean = merged[merged['condition'].notna()].copy()
print(f"\nWorking with {len(merged_clean):,} cells with condition labels")

# Compute Control mean and std for each cell type
control_stats = merged_clean[merged_clean['condition'] == 'Control'].groupby('cell_type_final')['STMN2_expr'].agg(['mean', 'std']).reset_index()
control_stats.columns = ['cell_type_final', 'control_mean', 'control_std']

print("\nControl STMN2 statistics by cell type:")
for idx, row in control_stats.iterrows():
    print(f"  {row['cell_type_final']}: μ={row['control_mean']:.3f}, σ={row['control_std']:.3f}")

# Merge control stats
merged_clean = pd.merge(merged_clean, control_stats, on='cell_type_final', how='left')

# Compute φ = [(expr - μ_Control) / σ_Control]²
merged_clean['STMN2_zscore'] = (merged_clean['STMN2_expr'] - merged_clean['control_mean']) / merged_clean['control_std']
merged_clean['phi_STMN2_single'] = merged_clean['STMN2_zscore'] ** 2

print(f"\nφ_STMN2_single statistics:")
print(merged_clean['phi_STMN2_single'].describe())

# Step 5: Create PT_STMN2_single
print("\n" + "=" * 80)
print("STEP 5: Creating PT_STMN2_single pseudotime")
print("=" * 80)

# PT_STMN2_single = φ_STMN2_single (single gene, no averaging needed)
merged_clean['PT_STMN2_single'] = merged_clean['phi_STMN2_single']

print(f"\nPT_STMN2_single statistics:")
print(merged_clean['PT_STMN2_single'].describe())

# Compare Control vs ALS
print("\n" + "-" * 80)
print("PT_STMN2_single by condition:")
print("-" * 80)
for condition in ['Control', 'ALS']:
    subset = merged_clean[merged_clean['condition'] == condition]
    if len(subset) > 0:
        print(f"{condition:8s}: {subset['PT_STMN2_single'].mean():8.4f} ± {subset['PT_STMN2_single'].std():8.4f} (n={len(subset):,})")

# Statistical test
control_pt = merged_clean[merged_clean['condition'] == 'Control']['PT_STMN2_single'].dropna()
als_pt = merged_clean[merged_clean['condition'] == 'ALS']['PT_STMN2_single'].dropna()

if len(control_pt) > 0 and len(als_pt) > 0:
    t_stat, p_val = stats.ttest_ind(control_pt, als_pt)
    print(f"\nT-test: t={t_stat:.4f}, p={p_val:.2e}")
    print(f"Effect size (Cohen's d): {(als_pt.mean() - control_pt.mean()) / np.sqrt((als_pt.var() + control_pt.var()) / 2):.4f}")

# Save PT_STMN2_single data
output_file = OUTPUT_DIR / "PT_STMN2_single_all_cells.csv"
merged_clean.to_csv(output_file, index=False)
print(f"\nSaved: {output_file}")

# Step 6: Compare with PT_STMN2_pathway and PT_TDP
print("\n" + "=" * 80)
print("STEP 6: Comparing PT_STMN2_single vs PT_STMN2_pathway vs PT_TDP")
print("=" * 80)

# Load PT_STMN2_pathway data
pt_pathway_path = "results/phaseTDP_STMN2/cell_level_features_ALL_with_PT_STMN2.csv"
print(f"\nLoading PT_STMN2_pathway: {pt_pathway_path}")

if Path(pt_pathway_path).exists():
    df_pathway = pd.read_csv(pt_pathway_path)
    print(f"Loaded {len(df_pathway):,} cells with PT_STMN2_pathway")

    # Merge on cell_id
    comparison = pd.merge(
        merged_clean[['cell_id', 'cell_type_final', 'condition', 'STMN2_expr', 'PT_STMN2_single']],
        df_pathway[['cell_id', 'PT_STMN2', 'PT_TDP', 'PT_dpt']],
        on='cell_id',
        how='inner'
    )

    comparison.rename(columns={'PT_STMN2': 'PT_STMN2_pathway'}, inplace=True)

    print(f"\nComparison data: {len(comparison):,} cells")

    # Compute correlations
    print("\n" + "=" * 80)
    print("CORRELATION ANALYSIS")
    print("=" * 80)

    corr_cols = ['PT_STMN2_single', 'PT_STMN2_pathway', 'PT_TDP', 'PT_dpt']
    corr_data = comparison[corr_cols].dropna()

    corr_matrix = corr_data.corr()
    print("\nCorrelation matrix:")
    print(corr_matrix)

    # Key comparisons
    print("\n" + "-" * 80)
    print("KEY CORRELATIONS:")
    print("-" * 80)
    print(f"PT_STMN2_single vs PT_STMN2_pathway: r = {corr_matrix.loc['PT_STMN2_single', 'PT_STMN2_pathway']:.4f}")
    print(f"PT_STMN2_single vs PT_TDP:           r = {corr_matrix.loc['PT_STMN2_single', 'PT_TDP']:.4f}")
    print(f"PT_STMN2_single vs PT_dpt:           r = {corr_matrix.loc['PT_STMN2_single', 'PT_dpt']:.4f}")
    print(f"PT_STMN2_pathway vs PT_TDP:          r = {corr_matrix.loc['PT_STMN2_pathway', 'PT_TDP']:.4f}")

    # Save correlation matrix
    corr_matrix.to_csv(OUTPUT_DIR / "PT_comparison_correlations.csv")

    # Visualize correlations
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    axes = axes.flatten()

    comparisons_list = [
        ('PT_STMN2_single', 'PT_STMN2_pathway', 'Single Gene vs Pathway (9 genes)'),
        ('PT_STMN2_single', 'PT_TDP', 'Single Gene vs TDP-43 Pathology'),
        ('PT_STMN2_single', 'PT_dpt', 'Single Gene vs Disease Progression'),
        ('PT_STMN2_pathway', 'PT_TDP', 'Pathway (9 genes) vs TDP-43'),
        ('PT_STMN2_pathway', 'PT_dpt', 'Pathway (9 genes) vs Disease'),
        ('PT_TDP', 'PT_dpt', 'TDP-43 Pathology vs Disease')
    ]

    for idx, (x_col, y_col, title) in enumerate(comparisons_list):
        ax = axes[idx]

        # Scatter plot with condition coloring
        for condition, color in [('Control', 'blue'), ('ALS', 'red')]:
            subset = comparison[comparison['condition'] == condition]
            if len(subset) > 0:
                ax.scatter(subset[x_col], subset[y_col], alpha=0.3, s=10, c=color, label=condition)

        # Compute correlation
        r, p = stats.pearsonr(comparison[x_col].dropna(), comparison[y_col].dropna())

        ax.set_xlabel(x_col.replace('_', ' '), fontsize=11)
        ax.set_ylabel(y_col.replace('_', ' '), fontsize=11)
        ax.set_title(f'{title}\nr = {r:.3f}, p = {p:.2e}', fontsize=12, fontweight='bold')
        ax.legend()
        ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / "Fig1_PT_comparisons_scatter.png", dpi=300, bbox_inches='tight')
    print(f"\nSaved: {OUTPUT_DIR / 'Fig1_PT_comparisons_scatter.png'}")
    plt.close()

    # Heatmap of correlations
    fig, ax = plt.subplots(figsize=(10, 8))
    sns.heatmap(corr_matrix, annot=True, fmt='.3f', cmap='RdBu_r', center=0,
                vmin=-1, vmax=1, square=True, ax=ax, cbar_kws={'label': 'Pearson r'},
                linewidths=1, linecolor='gray')
    ax.set_title('Correlation Matrix: PT Metrics Comparison\nSTMN2 Single Gene vs Pathway vs TDP-43 vs Disease',
                 fontsize=14, fontweight='bold', pad=20)
    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / "Fig2_PT_correlation_heatmap.png", dpi=300, bbox_inches='tight')
    print(f"Saved: {OUTPUT_DIR / 'Fig2_PT_correlation_heatmap.png'}")
    plt.close()

    # Save comparison data
    comparison.to_csv(OUTPUT_DIR / "PT_comparison_all_cells.csv", index=False)
    print(f"Saved: {OUTPUT_DIR / 'PT_comparison_all_cells.csv'}")

else:
    print(f"WARNING: {pt_pathway_path} not found. Skipping comparison.")

# Step 7: Visualize STMN2 single gene expression
print("\n" + "=" * 80)
print("STEP 7: Visualizing STMN2 expression and PT_STMN2_single")
print("=" * 80)

fig, axes = plt.subplots(2, 2, figsize=(16, 12))

# 7.1: STMN2 expression by condition
ax = axes[0, 0]
sns.violinplot(data=merged_clean, x='condition', y='STMN2_expr', ax=ax, palette=['blue', 'red'])
sns.swarmplot(data=merged_clean.sample(min(1000, len(merged_clean))),
              x='condition', y='STMN2_expr', ax=ax, color='black', alpha=0.3, size=2)
ax.set_title('STMN2 Gene Expression by Condition', fontsize=14, fontweight='bold')
ax.set_ylabel('STMN2 Expression (log-normalized)', fontsize=12)
ax.set_xlabel('Condition', fontsize=12)

# 7.2: PT_STMN2_single by condition
ax = axes[0, 1]
sns.violinplot(data=merged_clean, x='condition', y='PT_STMN2_single', ax=ax, palette=['blue', 'red'])
ax.set_title('PT_STMN2_single by Condition', fontsize=14, fontweight='bold')
ax.set_ylabel('PT_STMN2_single (φ energy)', fontsize=12)
ax.set_xlabel('Condition', fontsize=12)
ax.set_ylim(0, merged_clean['PT_STMN2_single'].quantile(0.99))

# 7.3: STMN2 expression by cell type
ax = axes[1, 0]
cell_type_order = merged_clean.groupby('cell_type_final')['STMN2_expr'].median().sort_values(ascending=False).index
sns.boxplot(data=merged_clean, x='cell_type_final', y='STMN2_expr', ax=ax, order=cell_type_order)
ax.set_title('STMN2 Expression by Cell Type', fontsize=14, fontweight='bold')
ax.set_ylabel('STMN2 Expression', fontsize=12)
ax.set_xlabel('Cell Type', fontsize=12)
ax.tick_params(axis='x', rotation=45, labelsize=10)

# 7.4: PT_STMN2_single by cell type
ax = axes[1, 1]
pt_order = merged_clean.groupby('cell_type_final')['PT_STMN2_single'].median().sort_values(ascending=False).index
sns.boxplot(data=merged_clean, x='cell_type_final', y='PT_STMN2_single', ax=ax, order=pt_order)
ax.set_title('PT_STMN2_single by Cell Type', fontsize=14, fontweight='bold')
ax.set_ylabel('PT_STMN2_single', fontsize=12)
ax.set_xlabel('Cell Type', fontsize=12)
ax.tick_params(axis='x', rotation=45, labelsize=10)

plt.tight_layout()
plt.savefig(OUTPUT_DIR / "Fig3_STMN2_single_overview.png", dpi=300, bbox_inches='tight')
print(f"Saved: {OUTPUT_DIR / 'Fig3_STMN2_single_overview.png'}")
plt.close()

print("\n" + "=" * 80)
print("PHASE STMN2 SINGLE GENE ANALYSIS COMPLETE!")
print("=" * 80)
print(f"\nOutput files saved to: {OUTPUT_DIR}/")
print("\nKEY FINDINGS:")
print(f"  - Total cells analyzed: {len(merged_clean):,}")
print(f"  - Cell types with STMN2: {merged_clean['cell_type_final'].nunique()}")
print(f"  - PT_STMN2_single range: [{merged_clean['PT_STMN2_single'].min():.4f}, {merged_clean['PT_STMN2_single'].max():.4f}]")
if len(control_pt) > 0 and len(als_pt) > 0:
    print(f"  - Control: {control_pt.mean():.4f} ± {control_pt.std():.4f}")
    print(f"  - ALS:     {als_pt.mean():.4f} ± {als_pt.std():.4f}")
    print(f"  - Difference: {als_pt.mean() - control_pt.mean():.4f} (p={p_val:.2e})")
