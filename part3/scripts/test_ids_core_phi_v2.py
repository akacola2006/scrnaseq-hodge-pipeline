#!/usr/bin/env python3
"""
Test: ids_core.compute_phi_from_scores() - v2
==============================================
Test with cell_type-specific Control statistics (matching original implementation).
"""

import pandas as pd
import numpy as np
from pathlib import Path
import sys

# Import ids_core
sys.path.insert(0, str(Path(__file__).parent))
import ids_core

print("=" * 80)
print("TEST: ids_core.compute_phi_from_scores() - cell_type-specific")
print("=" * 80)

# Load existing STMN2_single data
data_path = Path("results/phaseTDP_STMN2_single/PT_STMN2_single_all_cells.csv")

if not data_path.exists():
    print(f"ERROR: {data_path} not found.")
    sys.exit(1)

print(f"\nLoading: {data_path}")
df = pd.read_csv(data_path)
print(f"Loaded {len(df):,} cells")

# Compute φ using ids_core, cell_type by cell_type
print("\nComputing φ using ids_core (cell_type-specific Control stats)...")

phi_ids_core_list = []

for cell_type in df['cell_type_final'].unique():
    # Subset for this cell type
    subset = df[df['cell_type_final'] == cell_type].copy()

    # Extract scores and control mask
    scores = subset.set_index('cell_id')['STMN2_expr']
    control_mask = subset.set_index('cell_id')['condition'] == 'Control'

    # Compute φ using ids_core
    phi = ids_core.compute_phi_from_scores(scores, control_mask, eps=1e-8)

    # Convert back to DataFrame
    phi_df = phi.reset_index()
    phi_df.columns = ['cell_id', 'phi_ids_core']
    phi_df['cell_type_final'] = cell_type

    phi_ids_core_list.append(phi_df)

# Combine
df_phi_core = pd.concat(phi_ids_core_list, ignore_index=True)
print(f"Computed φ for {len(df_phi_core):,} cells across {len(df['cell_type_final'].unique())} cell types")

# Merge with original data
df_merged = df.merge(df_phi_core, on=['cell_id', 'cell_type_final'], how='inner')

# Compare
phi_original = df_merged['phi_STMN2_single']
phi_ids_core = df_merged['phi_ids_core']

diff = phi_ids_core - phi_original

print("\n" + "=" * 80)
print("COMPARISON")
print("=" * 80)

print(f"\nMatched cells: {len(df_merged):,}")

print(f"\n  Original φ stats:")
print(f"    Mean: {phi_original.mean():.6f}")
print(f"    Std:  {phi_original.std():.6f}")
print(f"    Min:  {phi_original.min():.6f}")
print(f"    Max:  {phi_original.max():.6f}")

print(f"\n  ids_core φ stats:")
print(f"    Mean: {phi_ids_core.mean():.6f}")
print(f"    Std:  {phi_ids_core.std():.6f}")
print(f"    Min:  {phi_ids_core.min():.6f}")
print(f"    Max:  {phi_ids_core.max():.6f}")

print(f"\n  Difference (ids_core - original):")
print(f"    Mean: {diff.mean():.2e}")
print(f"    Std:  {diff.std():.2e}")
print(f"    Max abs: {diff.abs().max():.2e}")

# Correlation
corr = phi_ids_core.corr(phi_original)
print(f"\n  Pearson correlation: r = {corr:.10f}")

# Success criterion
mean_diff = abs(diff.mean())
max_diff = diff.abs().max()

print("\n" + "=" * 80)
print("VERDICT")
print("=" * 80)

if mean_diff < 1e-6 and max_diff < 1e-4 and corr > 0.9999999:
    print("\n✅ TEST PASSED")
    print(f"   ids_core.compute_phi_from_scores() produces identical results")
    print(f"   when applied cell_type-specifically.")
    print(f"   Mean diff: {mean_diff:.2e} < 1e-6")
    print(f"   Max diff:  {max_diff:.2e} < 1e-4 (floating-point precision)")
    print(f"   Correlation: {corr:.10f} > 0.9999999")
    sys.exit(0)
else:
    print("\n❌ TEST FAILED")
    print(f"   Results differ from original implementation.")
    print(f"   Mean diff: {mean_diff:.2e}")
    print(f"   Max diff:  {max_diff:.2e}")
    print(f"   Correlation: {corr:.10f}")

    # Show top differences
    print("\n  Top 10 differences:")
    top_diff_idx = diff.abs().nlargest(10).index
    for i, idx in enumerate(top_diff_idx, 1):
        cell_id = df_merged.loc[idx, 'cell_id']
        cell_type = df_merged.loc[idx, 'cell_type_final']
        print(f"    {i}. cell {cell_id} ({cell_type}): orig={phi_original.iloc[idx]:.6f}, core={phi_ids_core.iloc[idx]:.6f}, diff={diff.iloc[idx]:.6e}")

    sys.exit(1)
