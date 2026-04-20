#!/usr/bin/env python3
"""
Phase 5a: Prepare Cell-Level Features for Cell-State Causality
===============================================================

Prepares integrated cell-level features across multiple cell types
to enable cell-state causal hierarchy analysis.

This script leverages existing PT-integrated data from the main pipeline:
- PT_continuous values (already computed via diffusion map)
- Module scores (I_m for all 23 modules)
- φ-metrics (delta_phi_log, SR_log for stress indicators)

Steps:
1. Select key cell types
2. For each cell type:
   - Load PT-integrated module data
   - Pivot to cell-level format (one row per cell)
   - Calculate composite stress from φ-metrics
3. Integrate into single cell_level_features_ALL.csv

Output:
- results/cell_state_causality/cell_level_features_ALL.csv
- results/cell_state_causality/cell_type_summary.csv
"""

import pandas as pd
import numpy as np
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# Paths
base_dir = Path('/home/akaco/als/motor_cortex_analysis')
ids_dir = base_dir / 'ids_causal_analysis'
results_dir = ids_dir / 'results'
cell_state_dir = results_dir / 'cell_state_causality'
cell_state_dir.mkdir(parents=True, exist_ok=True)

# Pre-computed PT data directory
pt_data_dir = base_dir / 'complete_ids_results_withPTcontinuous'

print("="*80)
print("Phase 5a: Prepare Cell-Level Features")
print("="*80)

# ============================================================================
# 1. Select key cell types
# ============================================================================
print("\n[1] Selecting key cell types...")

# Get ALL available cell types from PT data directory
import glob
pt_files = glob.glob(str(pt_data_dir / 'cell_module_*_synFull_withPTcontinuous.csv'))
key_celltypes = []
for f in pt_files:
    # Extract cell type name from filename
    basename = Path(f).name
    celltype_underscore = basename.replace('cell_module_', '').replace('_synFull_withPTcontinuous.csv', '')
    # Convert back to dot notation
    celltype = celltype_underscore.replace('_', '.')
    # Handle special cases
    celltype = celltype.replace('GFAP.neg', 'GFAP-neg').replace('GFAP.pos', 'GFAP-pos')
    key_celltypes.append(celltype)

key_celltypes = sorted(key_celltypes)

print(f"  Selected {len(key_celltypes)} cell types:")
for ct in key_celltypes:
    print(f"    - {ct}")

# ============================================================================
# 2. Process each cell type
# ============================================================================
print("\n[2] Processing cell types...")

all_cells_data = []
cell_type_summary = []

for celltype in key_celltypes:
    print(f"\n  Processing: {celltype}")

    # Convert to underscore format for file naming
    celltype_file = celltype.replace('.', '_').replace('-', '_')

    # Load PT-integrated data
    pt_file = pt_data_dir / f'cell_module_{celltype_file}_synFull_withPTcontinuous.csv'

    if not pt_file.exists():
        print(f"    WARNING: File not found: {pt_file.name}, skipping")
        continue

    # Load data (format: each row is cell × module combination)
    df = pd.read_csv(pt_file)
    print(f"    Loaded PT data: {df.shape[0]} rows (cell × module combinations)")

    # Get unique cells
    unique_cells = df['cell_id'].unique()
    print(f"    Unique cells: {len(unique_cells)}")

    # Pivot to cell-level format (one row per cell)
    # Each module becomes a feature column
    print(f"    Pivoting to cell-level format...")

    # Extract module scores (I_m)
    module_pivot = df.pivot_table(
        index='cell_id',
        columns='module',
        values='I_m',
        aggfunc='first'
    )
    module_pivot.columns = [f'module_{col}' for col in module_pivot.columns]

    # Extract PT and condition (same for all modules of a cell)
    cell_info = df.groupby('cell_id').first()[['Disease', 'PT_continuous']].copy()
    cell_info.rename(columns={'Disease': 'condition', 'PT_continuous': 'PT_IDS'}, inplace=True)

    # Calculate composite stress from φ-metrics
    # Using delta_phi_log as stress indicator (deviation from golden ratio)
    if 'delta_phi_log' in df.columns:
        stress_pivot = df.pivot_table(
            index='cell_id',
            columns='module',
            values='delta_phi_log',
            aggfunc='first'
        )
        # Composite stress = mean absolute deviation across stress-related modules
        stress_modules = ['ER_Stress', 'Oxidative_Stress', 'Mitochondria', 'Protein_Homeostasis']
        available_stress_modules = [m for m in stress_modules if m in stress_pivot.columns]

        if available_stress_modules:
            stress_total = stress_pivot[available_stress_modules].abs().mean(axis=1)
        else:
            # Fallback: use all modules
            stress_total = stress_pivot.abs().mean(axis=1)

        cell_info['stress_total'] = stress_total
    else:
        print(f"    WARNING: delta_phi_log not found, setting stress to 0")
        cell_info['stress_total'] = 0

    # Combine all features
    cell_features = pd.concat([cell_info, module_pivot], axis=1)
    cell_features['cell_type'] = celltype
    cell_features.reset_index(inplace=True)

    # Append to master list
    all_cells_data.append(cell_features)

    # Summary statistics
    n_cells = len(cell_features)
    n_als = (cell_features['condition'] == 'ALS').sum()
    n_control = (cell_features['condition'].isin(['Control', 'CTRL'])).sum()
    mean_stress = cell_features['stress_total'].mean()
    mean_pt = cell_features['PT_IDS'].mean()

    cell_type_summary.append({
        'cell_type': celltype,
        'n_cells': n_cells,
        'n_ALS': n_als,
        'n_Control': n_control,
        'mean_stress': mean_stress,
        'mean_PT': mean_pt
    })

    print(f"    Completed: {n_cells} cells")
    print(f"      ALS: {n_als}, Control: {n_control}")
    print(f"      Mean PT: {mean_pt:.3f}, Mean stress: {mean_stress:.3f}")

# ============================================================================
# 3. Create integrated dataframe
# ============================================================================
print("\n[3] Creating integrated cell-level features...")

if len(all_cells_data) == 0:
    print("  ERROR: No data loaded! Exiting.")
    exit(1)

all_cells_df = pd.concat(all_cells_data, ignore_index=True)
print(f"  Total cells: {len(all_cells_df)}")
print(f"  Cell types: {all_cells_df['cell_type'].nunique()}")
print(f"  Features: {len(all_cells_df.columns)}")
print(f"  Columns: {list(all_cells_df.columns[:10])}... (showing first 10)")

# Reorder columns for clarity
col_order = ['cell_id', 'cell_type', 'condition', 'PT_IDS', 'stress_total']
module_cols = [c for c in all_cells_df.columns if c.startswith('module_')]
col_order.extend(module_cols)
all_cells_df = all_cells_df[col_order]

# Save
output_file = cell_state_dir / 'cell_level_features_ALL.csv'
all_cells_df.to_csv(output_file, index=False)
print(f"\n  Saved: {output_file}")

# Save summary
summary_df = pd.DataFrame(cell_type_summary)
summary_file = cell_state_dir / 'cell_type_summary.csv'
summary_df.to_csv(summary_file, index=False)
print(f"  Saved: {summary_file}")

print("\n  Cell type summary:")
print(summary_df.to_string(index=False))

print("\n" + "="*80)
print("Phase 5a completed!")
print("="*80)
