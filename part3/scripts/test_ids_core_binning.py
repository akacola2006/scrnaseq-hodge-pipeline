#!/usr/bin/env python3
"""
Test: ids_core PT binning and φ aggregation
=============================================
Verify that bin_by_pt() and summarize_phi_by_group_and_bin() produce
identical results to Phase13c_NVU_energy_complete.py
"""

import pandas as pd
import numpy as np
from pathlib import Path
import sys

# Import ids_core
sys.path.insert(0, str(Path(__file__).parent))
import ids_core

print("=" * 80)
print("TEST: ids_core PT binning and φ aggregation")
print("=" * 80)

# Configuration (matching Phase13c)
INPUT_FILE = 'results/phase9_vascular/cell_level_features_ALL_with_PTdpt.csv'
REFERENCE_FILE = 'results/phase13_nvu_energy/nvu_phi_profiles_by_node_and_bin.csv'
TEST_MODULE = 'Mitochondria'  # Use this module for testing
NUM_BINS = 30

# NVU node definitions (matching Phase13c)
NVU_NODES = {
    'Vascular': [
        'Vasc.Endo.Capillary',
        'Vasc.Endo.Venous',
        'Vasc.Endo.Arterial',
        'Vasc.Mural.Pericyte',
        'Vasc.Mural.SMC',
        'Vasc.Fibro.CLMP.PDGFRA'
    ],
    'Glia': [
        'Glia.Oligo',
        'Glia.Astro.GFAP-neg',
        'Glia.Micro'
    ],
    'Upper': [
        'Ex.L2.L3.CUX2.RASGRF2'
    ],
    'VAT1L': [
        'Ex.L5.VAT1L.EYA4',
        'Ex.L5.VAT1L.THSD4'
    ]
}

# ========================================================================
# LOAD DATA
# ========================================================================

print(f"\nLoading: {INPUT_FILE}")
if not Path(INPUT_FILE).exists():
    print(f"ERROR: {INPUT_FILE} not found")
    sys.exit(1)

df = pd.read_csv(INPUT_FILE)
print(f"Loaded {len(df):,} cells")

# Create NVU node labels
def assign_nvu_node(cell_type):
    for node, types in NVU_NODES.items():
        if cell_type in types:
            return node
    return None

df['nvu_node'] = df['cell_type'].apply(assign_nvu_node)
df_nvu = df[df['nvu_node'].notna()].copy()

# Check module column exists
col_name = f'module_{TEST_MODULE}'
if col_name not in df_nvu.columns:
    print(f"ERROR: {col_name} not found in data")
    sys.exit(1)

# Split Control and ALS
df_ctrl = df_nvu[df_nvu['condition'] == 'Control'].copy()
df_als = df_nvu[df_nvu['condition'] == 'ALS'].copy()

print(f"\nNVU cells: {len(df_nvu):,}")
print(f"Control cells: {len(df_ctrl):,}")
print(f"ALS cells: {len(df_als):,}")

# ========================================================================
# COMPUTE φ USING ORIGINAL METHOD
# ========================================================================

print(f"\nComputing φ_{TEST_MODULE} using original method...")

# Control statistics (global)
ctrl_vals = df_ctrl[col_name].dropna()
mu = ctrl_vals.mean()
sigma = ctrl_vals.std()  # default ddof=1

Z = (df_nvu[col_name] - mu) / sigma
df_nvu[f'phi_{TEST_MODULE}'] = Z ** 2

# Filter to ALS
df_als = df_nvu[df_nvu['condition'] == 'ALS'].copy()

print(f"  μ_Control = {mu:.4f}")
print(f"  σ_Control = {sigma:.4f}")
print(f"  φ range: [{df_als[f'phi_{TEST_MODULE}'].min():.4f}, {df_als[f'phi_{TEST_MODULE}'].max():.4f}]")

# ========================================================================
# BIN AND AGGREGATE USING ids_core
# ========================================================================

print(f"\nBinning and aggregating using ids_core...")

# Extract PT and φ values (ALS only)
pt = df_als.set_index('cell_id')['PT_dpt']
phi = df_als.set_index('cell_id')[f'phi_{TEST_MODULE}']
group = df_als.set_index('cell_id')['nvu_node']

# Use ids_core to bin PT
bin_index, bin_edges, bin_centers = ids_core.bin_by_pt(pt, n_bins=NUM_BINS)

print(f"  Created {NUM_BINS} bins")
print(f"  PT range: [{pt.min():.4f}, {pt.max():.4f}]")
print(f"  Bin edges: [{bin_edges[0]:.4f}, ..., {bin_edges[-1]:.4f}]")

# Use ids_core to aggregate
summary = ids_core.summarize_phi_by_group_and_bin(
    phi=phi,
    group=group,
    bin_index=bin_index,
    group_order=['Vascular', 'Glia', 'Upper', 'VAT1L']
)

print(f"  Aggregated to {len(summary)} (group, bin) combinations")

# Add bin centers for comparison
summary['pt_bin_center'] = summary['bin'].apply(lambda b: bin_centers[int(b)])

# Rename for consistency with original
summary = summary.rename(columns={'group': 'node', 'bin': 'pt_bin'})

# ========================================================================
# LOAD REFERENCE DATA
# ========================================================================

print(f"\nLoading reference: {REFERENCE_FILE}")
if not Path(REFERENCE_FILE).exists():
    print(f"ERROR: {REFERENCE_FILE} not found")
    sys.exit(1)

ref = pd.read_csv(REFERENCE_FILE)
ref_module = ref[ref['module'] == TEST_MODULE].copy()

print(f"Loaded reference with {len(ref_module)} (node, bin) combinations for {TEST_MODULE}")

# ========================================================================
# COMPARE
# ========================================================================

print("\n" + "=" * 80)
print("COMPARISON")
print("=" * 80)

# Merge on (node, pt_bin)
merged = summary.merge(
    ref_module[['node', 'pt_bin', 'pt_bin_center', 'phi_mean', 'n_cells']],
    on=['node', 'pt_bin'],
    how='inner',
    suffixes=('_ids_core', '_original')
)

print(f"\nMatched {len(merged)} (node, bin) combinations")

# Compare phi_mean
phi_ids_core = merged['phi_mean_ids_core']
phi_original = merged['phi_mean_original']

diff = phi_ids_core - phi_original

print(f"\n  Original φ_mean stats:")
print(f"    Mean: {phi_original.mean():.6f}")
print(f"    Std:  {phi_original.std():.6f}")
print(f"    Min:  {phi_original.min():.6f}")
print(f"    Max:  {phi_original.max():.6f}")

print(f"\n  ids_core φ_mean stats:")
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

# Compare PT bin centers
pt_center_ids_core = merged['pt_bin_center_ids_core']
pt_center_original = merged['pt_bin_center_original']
pt_center_diff = pt_center_ids_core - pt_center_original

print(f"\n  PT bin center difference:")
print(f"    Mean: {pt_center_diff.mean():.2e}")
print(f"    Max abs: {pt_center_diff.abs().max():.2e}")

# Compare n_cells
n_cells_ids_core = merged['n_cells_ids_core']
n_cells_original = merged['n_cells_original']
n_cells_diff = n_cells_ids_core - n_cells_original

print(f"\n  Cell count difference:")
print(f"    Mean: {n_cells_diff.mean():.2e}")
print(f"    Max abs: {n_cells_diff.abs().max():.2e}")
print(f"    Non-zero diffs: {(n_cells_diff != 0).sum()}")

# ========================================================================
# VERDICT
# ========================================================================

mean_diff = abs(diff.mean())
max_diff = diff.abs().max()

print("\n" + "=" * 80)
print("VERDICT")
print("=" * 80)

if mean_diff < 1e-6 and max_diff < 1e-4 and corr > 0.9999999 and n_cells_diff.abs().max() == 0:
    print("\n✅ TEST PASSED")
    print(f"   ids_core binning and aggregation produce identical results.")
    print(f"   Mean φ diff: {mean_diff:.2e} < 1e-6")
    print(f"   Max φ diff:  {max_diff:.2e} < 1e-4")
    print(f"   Correlation: {corr:.10f} > 0.9999999")
    print(f"   Cell counts match exactly")
    sys.exit(0)
else:
    print("\n❌ TEST FAILED")
    print(f"   Results differ from original implementation.")
    print(f"   Mean φ diff: {mean_diff:.2e}")
    print(f"   Max φ diff:  {max_diff:.2e}")
    print(f"   Correlation: {corr:.10f}")
    print(f"   Max cell count diff: {n_cells_diff.abs().max()}")

    # Show top differences
    print("\n  Top 10 φ_mean differences:")
    top_diff_idx = diff.abs().nlargest(10).index
    for i, idx in enumerate(top_diff_idx, 1):
        node = merged.loc[idx, 'node']
        pt_bin = merged.loc[idx, 'pt_bin']
        orig = phi_original.iloc[idx]
        core = phi_ids_core.iloc[idx]
        d = diff.iloc[idx]
        print(f"    {i}. {node} bin {pt_bin}: orig={orig:.6f}, core={core:.6f}, diff={d:.6e}")

    sys.exit(1)
