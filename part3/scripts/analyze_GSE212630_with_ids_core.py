#!/usr/bin/env python3
"""
GSE_212630 Dataset Analysis with IDS Core
==========================================

This script analyzes the GSE_212630 raw expression data using the IDS Core library.

Dataset characteristics:
- 21 cell types (Glia, Excitatory neurons, Inhibitory neurons, Vascular)
- 36,601 genes
- 4 conditions: Control, TDPneg, TDPmed, TDPhigh (TDP-43 pathology stages)
- Total cells: ~60,000+

Analysis pipeline:
1. Load expression and metadata for all cell types
2. Compute ALS stress components
3. Apply IDS Core φ-flow analysis
4. Generate summary reports and visualizations

Author: Claude Code + User
Date: 2025-12-01
"""

import pandas as pd
import numpy as np
import sys
from pathlib import Path
import warnings
import gc
warnings.filterwarnings('ignore')

# Add script directory to path
SCRIPT_DIR = Path(__file__).parent
sys.path.insert(0, str(SCRIPT_DIR))

# Import project modules
import ids_core
sys.path.insert(0, str(SCRIPT_DIR.parent / 'config'))
from als_stress_config import ALS_STRESS_COMPONENTS, compute_als_stress

print("=" * 80)
print("GSE_212630 Dataset Analysis with IDS Core")
print("=" * 80)

# ============================================================================
# CONFIGURATION
# ============================================================================

DATA_DIR = Path(__file__).parent.parent / 'GSE_212630_raw_expression_transposed'
OUTPUT_DIR = Path(__file__).parent.parent / 'results' / 'GSE212630_ids_analysis'
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# TDP-43 staging to pseudotime mapping
TDP_STAGE_TO_PT = {
    'Control': 0.0,
    'TDPneg': 0.33,
    'TDPmed': 0.67,
    'TDPhigh': 1.0
}

# Cell type groupings for NVU analysis
NVU_MAPPING = {
    'Vascular': ['Vasc.Capillary', 'Vasc.Endo', 'Vasc.Fibro', 'Vasc.Pericyte', 'Vasc.SMC', 'Vasc.Unknown'],
    'Glia': ['Glia.Oligo', 'Glia.Astro.GFAP.neg', 'Glia.Astro.GFAP.pos', 'Glia.OPC', 'Glia.Micro'],
    'Excitatory': ['Ex.L2/3', 'Ex.L4', 'Ex.L5', 'Ex.L6', 'Ex.Unknown'],
    'Inhibitory': ['In.SST', 'In.VIP', 'In.DISC1', 'In.PVALB', 'In.Unknown']
}

# Minimum cells threshold
MIN_CELLS = 50

# Maximum cells per cell type (to avoid memory issues)
MAX_CELLS_PER_TYPE = 1000

# Random seed for reproducibility
np.random.seed(42)

# ============================================================================
# STEP 1: LOAD DATA
# ============================================================================

print("\n" + "=" * 80)
print("STEP 1: Loading data")
print("=" * 80)

# Read extraction summary
summary_df = pd.read_csv(DATA_DIR / 'extraction_summary.csv')
print(f"\nAvailable cell types: {len(summary_df)}")
print(summary_df[['cell_type', 'n_cells', 'n_Control', 'n_TDPneg', 'n_TDPmed', 'n_TDPhigh']].to_string(index=False))

# Filter cell types with sufficient cells, sort by size (smallest first to test)
# Exclude very large cell types to avoid memory issues
valid_df = summary_df[(summary_df['n_cells'] >= MIN_CELLS) & (summary_df['n_cells'] <= 10000)].sort_values('n_cells')
valid_cell_types = valid_df['cell_type'].tolist()
print(f"\nCell types with {MIN_CELLS}-10000 cells: {len(valid_cell_types)} (processing smallest first)")

# Load all data
all_metadata = []
all_stress_scores = []

print("\nLoading expression data for each cell type...")

for ct in valid_cell_types:
    # Find files
    ct_safe = ct.replace('/', '_').replace('.', '_')
    expr_file = None
    meta_file = None

    for f in DATA_DIR.glob('*_expression_transposed.csv.gz'):
        # Match cell type
        fname = f.stem.replace('_expression_transposed.csv', '')
        if ct_safe.lower() in fname.lower() or fname.replace('_', '.').lower() == ct.lower().replace('/', '_'):
            expr_file = f
            meta_file = DATA_DIR / f.name.replace('_expression_transposed.csv.gz', '_metadata.csv')
            break

    if expr_file is None:
        # Try direct match
        for f in DATA_DIR.glob('*_expression_transposed.csv.gz'):
            # Build expected filename
            ct_fname = ct.replace('.', '_').replace('/', '_')
            if ct_fname in f.name:
                expr_file = f
                meta_file = DATA_DIR / f.name.replace('_expression_transposed.csv.gz', '_metadata.csv')
                break

    if expr_file is None or not meta_file.exists():
        print(f"  [SKIP] {ct}: files not found")
        continue

    print(f"\n  Loading {ct}...")

    try:
        # Load metadata
        meta_df = pd.read_csv(meta_file)
        meta_df['cell_type'] = ct

        # Map condition to pseudotime
        meta_df['PT_TDP'] = meta_df['condition'].map(TDP_STAGE_TO_PT)

        # Map to NVU node
        nvu_node = 'Other'
        for node, cell_types in NVU_MAPPING.items():
            if ct in cell_types:
                nvu_node = node
                break
        meta_df['nvu_node'] = nvu_node

        # Sample if too many cells
        if len(meta_df) > MAX_CELLS_PER_TYPE:
            print(f"    Sampling {MAX_CELLS_PER_TYPE} cells from {len(meta_df)}...")
            # Stratified sampling by condition
            sampled_idx = []
            for cond in meta_df['condition'].unique():
                cond_idx = meta_df[meta_df['condition'] == cond].index.tolist()
                n_sample = min(len(cond_idx), MAX_CELLS_PER_TYPE // meta_df['condition'].nunique())
                sampled_idx.extend(np.random.choice(cond_idx, size=n_sample, replace=False).tolist())
            meta_df = meta_df.loc[sampled_idx].reset_index(drop=True)

        # Build column names: {cell_id}_{sample}_{condition}
        meta_df['expr_col'] = meta_df['cell_id'] + '_' + meta_df['sample'].astype(str) + '_' + meta_df['condition']
        cells_to_load = meta_df['expr_col'].tolist()

        # Load expression (all columns, then filter)
        print(f"    Loading expression file...")
        expr_df = pd.read_csv(expr_file, index_col=0)
        print(f"    Raw expression shape: {expr_df.shape}")

        # Filter to needed columns
        available_cols = [c for c in cells_to_load if c in expr_df.columns]
        if len(available_cols) == 0:
            print(f"    [SKIP] No matching columns found")
            continue
        if len(available_cols) < len(cells_to_load) * 0.5:
            print(f"    [WARN] Only {len(available_cols)}/{len(cells_to_load)} columns found")

        expr_df = expr_df[available_cols]
        # Update metadata to match
        meta_df = meta_df[meta_df['expr_col'].isin(available_cols)].reset_index(drop=True)

        print(f"    Filtered expression shape: {expr_df.shape} (genes × cells)")
        print(f"    Filtered metadata shape: {meta_df.shape}")

        # Transpose to cells × genes, align with metadata
        expr_T = expr_df.T
        expr_T = expr_T.loc[meta_df['expr_col'].values]

        # Verify alignment
        if len(meta_df) != expr_T.shape[0]:
            print(f"    [SKIP] Shape mismatch after alignment: meta={len(meta_df)}, expr={expr_T.shape[0]}")
            continue

        # Compute ALS stress scores
        gene_names = expr_T.columns.tolist()
        X = expr_T.values

        print(f"    Computing stress scores...")
        S_total, S_components = compute_als_stress(X, gene_names, return_components=True)

        # Add to metadata
        meta_df['stress_total'] = S_total
        for comp_name, comp_scores in S_components.items():
            meta_df[f'stress_{comp_name}'] = comp_scores

        all_metadata.append(meta_df)

        print(f"    Done: {len(meta_df)} cells, mean stress = {S_total.mean():.4f}")

        # Memory cleanup
        del expr_df, expr_T, X
        gc.collect()

    except Exception as e:
        print(f"    [ERROR] {e}")
        gc.collect()
        continue

# Combine all metadata
print("\n" + "-" * 40)
print("Combining all cell types...")

if len(all_metadata) == 0:
    print("ERROR: No data loaded!")
    sys.exit(1)

combined_df = pd.concat(all_metadata, ignore_index=True)
print(f"\nTotal cells loaded: {len(combined_df):,}")
print(f"Cell types: {combined_df['cell_type'].nunique()}")
print(f"Conditions: {combined_df['condition'].unique().tolist()}")

# Save combined metadata
combined_df.to_csv(OUTPUT_DIR / 'combined_metadata_with_stress.csv', index=False)
print(f"\nSaved: combined_metadata_with_stress.csv")

# ============================================================================
# STEP 2: APPLY IDS CORE φ-FLOW ANALYSIS
# ============================================================================

print("\n" + "=" * 80)
print("STEP 2: Computing φ-flow analysis with IDS Core")
print("=" * 80)

# Use PT_TDP as pseudotime (0 = Control, 1 = TDPhigh)
pt_col = 'PT_TDP'
control_mask = combined_df['condition'] == 'Control'

print(f"\nControl cells: {control_mask.sum():,}")
print(f"Disease cells: {(~control_mask).sum():,}")

# Compute φ for stress_total
print("\nComputing φ for stress_total...")
phi_stress, z_stress = ids_core.compute_phi_from_scores(
    scores=pd.Series(combined_df['stress_total'].values, index=combined_df.index),
    control_mask=pd.Series(control_mask.values, index=combined_df.index),
    return_z=True
)
combined_df['phi_stress_total'] = phi_stress.values
combined_df['z_stress_total'] = z_stress.values

# Compute φ for each stress component
stress_components = [col for col in combined_df.columns if col.startswith('stress_') and col != 'stress_total']
for comp in stress_components:
    phi_comp = ids_core.compute_phi_from_scores(
        scores=pd.Series(combined_df[comp].values, index=combined_df.index),
        control_mask=pd.Series(control_mask.values, index=combined_df.index)
    )
    combined_df[f'phi_{comp}'] = phi_comp.values

print(f"Computed φ for {len(stress_components) + 1} stress components")

# ============================================================================
# STEP 3: AGGREGATE BY NVU NODE AND PT BIN
# ============================================================================

print("\n" + "=" * 80)
print("STEP 3: Aggregating φ by NVU node and PT bin")
print("=" * 80)

# Bin by PT_TDP
n_bins = 4  # Match TDP stages
pt_series = pd.Series(combined_df['PT_TDP'].values, index=combined_df.index)
bin_index, bin_edges, bin_centers = ids_core.bin_by_pt(pt_series, n_bins=n_bins)
combined_df['pt_bin'] = bin_index.values

print(f"PT bins: {n_bins}")
print(f"Bin edges: {bin_edges}")
print(f"Bin centers: {bin_centers}")

# Summarize φ by NVU node and PT bin
phi_stress_series = pd.Series(combined_df['phi_stress_total'].values, index=combined_df.index)
group_series = pd.Series(combined_df['nvu_node'].values, index=combined_df.index)

profiles = ids_core.summarize_phi_by_group_and_bin(
    phi=phi_stress_series,
    group=group_series,
    bin_index=bin_index,
    group_order=['Vascular', 'Glia', 'Excitatory', 'Inhibitory']
)

print("\nφ profiles by NVU node and PT bin:")
print(profiles.to_string(index=False))

# Detect onset and peak
flow_summary = ids_core.detect_onset_and_peak(
    phi_mean_df=profiles[['group', 'bin', 'phi_mean']].copy(),
    bin_centers=bin_centers,
    threshold=0.5
)

print("\n📊 Flow summary (onset and peak detection):")
print(flow_summary.to_string(index=False))

# Save results
profiles.to_csv(OUTPUT_DIR / 'phi_profiles_by_nvu_and_pt.csv', index=False)
flow_summary.to_csv(OUTPUT_DIR / 'flow_summary_by_nvu.csv', index=False)

# ============================================================================
# STEP 4: DETAILED ANALYSIS BY CELL TYPE
# ============================================================================

print("\n" + "=" * 80)
print("STEP 4: Detailed analysis by cell type")
print("=" * 80)

# Aggregate by cell_type
celltype_summary = combined_df.groupby(['cell_type', 'condition']).agg({
    'stress_total': ['mean', 'std', 'count'],
    'phi_stress_total': ['mean', 'std'],
    'PT_TDP': 'first'
}).round(4)

celltype_summary.columns = ['_'.join(col).strip() for col in celltype_summary.columns]
celltype_summary = celltype_summary.reset_index()

print("\nCell type × Condition summary:")
print(celltype_summary.head(20).to_string(index=False))

celltype_summary.to_csv(OUTPUT_DIR / 'celltype_condition_summary.csv', index=False)

# ============================================================================
# STEP 5: VISUALIZATION
# ============================================================================

print("\n" + "=" * 80)
print("STEP 5: Generating visualizations")
print("=" * 80)

try:
    import matplotlib.pyplot as plt
    import seaborn as sns

    # Set style
    plt.style.use('default')
    sns.set_palette('husl')

    # Figure 1: φ trajectories by NVU node
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Panel A: φ trajectory
    ax = axes[0]
    for node in ['Vascular', 'Glia', 'Excitatory', 'Inhibitory']:
        node_data = profiles[profiles['group'] == node]
        ax.plot(node_data['bin'], node_data['phi_mean'], 'o-', label=node, linewidth=2, markersize=8)

    ax.set_xlabel('TDP-43 Stage (0=Control, 3=TDPhigh)', fontsize=12)
    ax.set_ylabel('φ (Stress Energy)', fontsize=12)
    ax.set_title('φ Trajectory by NVU Node', fontsize=14)
    ax.legend(loc='upper left')
    ax.set_xticks([0, 1, 2, 3])
    ax.set_xticklabels(['Control', 'TDPneg', 'TDPmed', 'TDPhigh'])
    ax.grid(True, alpha=0.3)

    # Panel B: Heatmap
    ax = axes[1]
    pivot_df = profiles.pivot(index='group', columns='bin', values='phi_mean')
    pivot_df = pivot_df.reindex(['Vascular', 'Glia', 'Excitatory', 'Inhibitory'])

    sns.heatmap(pivot_df, annot=True, fmt='.2f', cmap='YlOrRd', ax=ax,
                cbar_kws={'label': 'φ (Stress Energy)'})
    ax.set_xlabel('TDP-43 Stage', fontsize=12)
    ax.set_ylabel('NVU Node', fontsize=12)
    ax.set_title('φ Heatmap by NVU Node × TDP Stage', fontsize=14)
    ax.set_xticklabels(['Control', 'TDPneg', 'TDPmed', 'TDPhigh'])

    plt.tight_layout()
    fig.savefig(OUTPUT_DIR / 'Fig1_phi_trajectories_by_nvu.png', dpi=150, bbox_inches='tight')
    print(f"  Saved: Fig1_phi_trajectories_by_nvu.png")
    plt.close()

    # Figure 2: Stress components comparison
    fig, ax = plt.subplots(figsize=(12, 6))

    # Compute mean φ by component and condition
    phi_cols = [col for col in combined_df.columns if col.startswith('phi_stress_') and col != 'phi_stress_total']

    comp_summary = []
    for col in phi_cols:
        comp_name = col.replace('phi_stress_', '')
        for cond in ['Control', 'TDPneg', 'TDPmed', 'TDPhigh']:
            mean_phi = combined_df[combined_df['condition'] == cond][col].mean()
            comp_summary.append({
                'Component': comp_name,
                'Condition': cond,
                'phi_mean': mean_phi
            })

    comp_df = pd.DataFrame(comp_summary)
    comp_pivot = comp_df.pivot(index='Component', columns='Condition', values='phi_mean')
    comp_pivot = comp_pivot[['Control', 'TDPneg', 'TDPmed', 'TDPhigh']]

    comp_pivot.plot(kind='bar', ax=ax, width=0.8)
    ax.set_xlabel('Stress Component', fontsize=12)
    ax.set_ylabel('Mean φ', fontsize=12)
    ax.set_title('φ by Stress Component × TDP-43 Stage', fontsize=14)
    ax.legend(title='Condition', loc='upper right')
    ax.tick_params(axis='x', rotation=45)

    plt.tight_layout()
    fig.savefig(OUTPUT_DIR / 'Fig2_phi_by_stress_component.png', dpi=150, bbox_inches='tight')
    print(f"  Saved: Fig2_phi_by_stress_component.png")
    plt.close()

    # Figure 3: Cell type detailed heatmap
    fig, ax = plt.subplots(figsize=(14, 10))

    # Compute mean φ by cell type and condition
    ct_phi = combined_df.groupby(['cell_type', 'condition'])['phi_stress_total'].mean().reset_index()
    ct_pivot = ct_phi.pivot(index='cell_type', columns='condition', values='phi_stress_total')
    ct_pivot = ct_pivot[['Control', 'TDPneg', 'TDPmed', 'TDPhigh']]

    sns.heatmap(ct_pivot, annot=True, fmt='.2f', cmap='YlOrRd', ax=ax,
                cbar_kws={'label': 'φ (Stress Energy)'})
    ax.set_xlabel('TDP-43 Stage', fontsize=12)
    ax.set_ylabel('Cell Type', fontsize=12)
    ax.set_title('φ Heatmap by Cell Type × TDP-43 Stage', fontsize=14)

    plt.tight_layout()
    fig.savefig(OUTPUT_DIR / 'Fig3_phi_by_celltype.png', dpi=150, bbox_inches='tight')
    print(f"  Saved: Fig3_phi_by_celltype.png")
    plt.close()

except Exception as e:
    print(f"  [WARN] Visualization error: {e}")

# ============================================================================
# STEP 6: SUMMARY REPORT
# ============================================================================

print("\n" + "=" * 80)
print("STEP 6: Generating summary report")
print("=" * 80)

# Key findings
print("\n📊 KEY FINDINGS:")
print("-" * 40)

# 1. Cell count by condition
cond_counts = combined_df['condition'].value_counts()
print("\n1. Cell counts by condition:")
for cond in ['Control', 'TDPneg', 'TDPmed', 'TDPhigh']:
    if cond in cond_counts.index:
        print(f"   {cond}: {cond_counts[cond]:,} cells")

# 2. Mean φ by NVU node and condition
print("\n2. Mean φ by NVU node (Control vs TDPhigh):")
for node in ['Vascular', 'Glia', 'Excitatory', 'Inhibitory']:
    ctrl_phi = combined_df[(combined_df['nvu_node'] == node) & (combined_df['condition'] == 'Control')]['phi_stress_total'].mean()
    high_phi = combined_df[(combined_df['nvu_node'] == node) & (combined_df['condition'] == 'TDPhigh')]['phi_stress_total'].mean()
    fold_change = high_phi / (ctrl_phi + 1e-8)
    print(f"   {node:12s}: Control={ctrl_phi:.3f}, TDPhigh={high_phi:.3f}, FC={fold_change:.2f}x")

# 3. Temporal ordering
print("\n3. Temporal ordering (onset PT):")
for _, row in flow_summary.iterrows():
    onset = row['onset_PT'] if pd.notna(row['onset_PT']) else 'N/A'
    peak = row['peak_PT'] if pd.notna(row['peak_PT']) else 'N/A'
    print(f"   {row['group']:12s}: onset_PT={onset}, peak_PT={peak}, peak_φ={row['peak_phi']:.3f}")

# 4. Most affected cell types
print("\n4. Top 5 most affected cell types (highest φ at TDPhigh):")
ct_tdphigh = combined_df[combined_df['condition'] == 'TDPhigh'].groupby('cell_type')['phi_stress_total'].mean().sort_values(ascending=False)
for ct, phi in ct_tdphigh.head(5).items():
    print(f"   {ct}: φ = {phi:.3f}")

# Save summary to file
summary_text = f"""
# GSE_212630 IDS Core Analysis Summary
Generated: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}

## Dataset Overview
- Total cells analyzed: {len(combined_df):,}
- Cell types: {combined_df['cell_type'].nunique()}
- Conditions: Control, TDPneg, TDPmed, TDPhigh

## Cell Counts by Condition
{cond_counts.to_string()}

## φ Summary by NVU Node
{profiles.to_string(index=False)}

## Flow Summary (Onset and Peak)
{flow_summary.to_string(index=False)}

## Top Affected Cell Types (TDPhigh)
{ct_tdphigh.head(10).to_string()}

## Files Generated
- combined_metadata_with_stress.csv
- phi_profiles_by_nvu_and_pt.csv
- flow_summary_by_nvu.csv
- celltype_condition_summary.csv
- Fig1_phi_trajectories_by_nvu.png
- Fig2_phi_by_stress_component.png
- Fig3_phi_by_celltype.png
"""

with open(OUTPUT_DIR / 'ANALYSIS_SUMMARY.md', 'w') as f:
    f.write(summary_text)

print(f"\n  Saved: ANALYSIS_SUMMARY.md")

print("\n" + "=" * 80)
print("✅ ANALYSIS COMPLETE")
print(f"Output directory: {OUTPUT_DIR}")
print("=" * 80)
