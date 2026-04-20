#!/usr/bin/env python3
"""
Phase 13c: Complete NVU Energy Landscape - REFACTORED with ids_core
====================================================================

**THIS IS THE REFACTORED VERSION using ids_core.py**

Extends Phase 13/13b φ-flow analysis to remaining modules:
  - Angiogenesis, Apoptosis, Autophagy, Cell_Cycle, DNA_Repair
  - ECM, Growth_Factors, Ion_Transport, Metabolism, Myelination

Key changes from original Phase13c_NVU_energy_complete.py:
  - φ computation → ids_core.compute_phi_from_scores()
  - PT binning → ids_core.bin_by_pt()
  - Aggregation → ids_core.summarize_phi_by_group_and_bin()
  - Onset/peak detection → ids_core.detect_onset_and_peak()

Author: Claude Code + User
Date: 2025-11-30 (IDS Core v1 - Refactored)
"""

import numpy as np
import pandas as pd
from pathlib import Path
import sys

# Import ids_core
sys.path.insert(0, str(Path(__file__).parent))
import ids_core

# ========================================================================
# CONFIGURATION
# ========================================================================

INPUT_FILE = 'results/phase9_vascular/cell_level_features_ALL_with_PTdpt.csv'
OUTPUT_DIR = 'results/phase13c_nvu_energy_refactored'

# Final 10 modules to analyze (completing the 23-module landscape)
FINAL_MODULES = [
    'Angiogenesis',
    'Apoptosis',
    'Autophagy',
    'Cell_Cycle',
    'DNA_Repair',
    'ECM',
    'Growth_Factors',
    'Ion_Transport',
    'Metabolism',
    'Myelination'
]

# NVU node definitions (MUST match Phase 13/13b for consistency!)
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
        'Ex.L2.L3.CUX2.RASGRF2'  # L2/L3 ONLY (as specified by user)
    ],
    'VAT1L': [
        'Ex.L5.VAT1L.EYA4',
        'Ex.L5.VAT1L.THSD4'
    ]
}

NUM_BINS = 30
PHI_THRESHOLD = 0.5

# Create output directory
Path(OUTPUT_DIR).mkdir(parents=True, exist_ok=True)

print("=" * 80)
print("Phase 13c: NVU Energy Landscape - REFACTORED with ids_core")
print("=" * 80)

# ========================================================================
# LOAD DATA
# ========================================================================

print("\nLoading cell-level features with PT_dpt...")
df = pd.read_csv(INPUT_FILE)
print(f"Loaded {len(df):,} cells, {len(df.columns)} features")

# Create NVU node labels using ids_core.map_cell_types_to_groups()
# (Note: the original used a custom function, but we can use ids_core for consistency)
# For now, keep the original assign_nvu_node since it's simple
def assign_nvu_node(cell_type):
    for node, types in NVU_NODES.items():
        if cell_type in types:
            return node
    return None

df['nvu_node'] = df['cell_type'].apply(assign_nvu_node)
df_nvu = df[df['nvu_node'].notna()].copy()
print(f"NVU cells: {len(df_nvu):,}")

# Check which modules actually exist in the data
available_modules = []
for module in FINAL_MODULES:
    col_name = f'module_{module}'
    if col_name in df_nvu.columns:
        available_modules.append(module)
    else:
        print(f"WARNING: {col_name} not found in data")

print(f"\nAnalyzing {len(available_modules)} final modules: {available_modules}")

# Split Control and ALS
df_ctrl = df_nvu[df_nvu['condition'] == 'Control'].copy()
df_als_full = df_nvu[df_nvu['condition'] == 'ALS'].copy()

print(f"Control cells: {len(df_ctrl):,}")
print(f"ALS cells: {len(df_als_full):,}")

# ========================================================================
# COMPUTE φ_m = Z² FROM CONTROL USING ids_core
# ========================================================================

print("\n" + "=" * 80)
print("COMPUTING φ USING ids_core.compute_phi_from_scores()")
print("=" * 80)

for module in available_modules:
    col_name = f'module_{module}'

    print(f"\nModule: {module}")

    # Extract scores and Control mask
    scores = df_nvu.set_index('cell_id')[col_name]
    control_mask = df_nvu.set_index('cell_id')['condition'] == 'Control'

    # Check for zero variance
    ctrl_vals = df_ctrl[col_name].dropna()
    if ctrl_vals.std() == 0:
        print(f"  WARNING: σ=0 in Control, skipping")
        continue

    # Compute φ using ids_core
    phi = ids_core.compute_phi_from_scores(scores, control_mask, eps=1e-8)

    # Add to DataFrame (align by index)
    # phi has cell_id as index, so we need to set df_nvu's index temporarily
    df_nvu = df_nvu.set_index('cell_id')
    df_nvu[f'phi_{module}'] = phi
    df_nvu = df_nvu.reset_index()

    print(f"  μ_Control = {ctrl_vals.mean():.4f}")
    print(f"  σ_Control = {ctrl_vals.std():.4f}")
    print(f"  φ range: [{phi.min():.4f}, {phi.max():.4f}]")

print(f"\n✓ Computed φ for {len(available_modules)} modules using ids_core")

# Filter to ALS cells (after calculating φ columns)
df_als = df_nvu[df_nvu['condition'] == 'ALS'].copy()

# ========================================================================
# BIN BY PT_dpt USING ids_core
# ========================================================================

print("\n" + "=" * 80)
print(f"BINNING BY PT_dpt USING ids_core.bin_by_pt() ({NUM_BINS} bins)")
print("=" * 80)

# Extract PT values (ALS only)
pt = df_als.set_index('cell_id')['PT_dpt']

# Use ids_core to bin PT
bin_index, bin_edges, bin_centers = ids_core.bin_by_pt(pt, n_bins=NUM_BINS)

print(f"\n✓ Created {NUM_BINS} bins")
print(f"  PT range: [{pt.min():.4f}, {pt.max():.4f}]")
print(f"  Bin edges: [{bin_edges[0]:.4f}, ..., {bin_edges[-1]:.4f}]")

# Add bin assignments back to DataFrame
# Important: bin_index has cell_id as index, so we need to align properly
df_als = df_als.set_index('cell_id')
df_als['PT_bin'] = bin_index
df_als = df_als.reset_index()

# Debug: check if φ columns exist
phi_cols_check = [f'phi_{m}' for m in available_modules]
print(f"\nDebug: Checking φ columns in df_als:")
for col in phi_cols_check:
    if col in df_als.columns:
        non_nan = df_als[col].notna().sum()
        print(f"  {col}: exists, {non_nan}/{len(df_als)} non-NaN values")
    else:
        print(f"  {col}: MISSING!")

# ========================================================================
# AGGREGATE φ BY (NVU_NODE, PT_BIN) USING ids_core
# ========================================================================

print("\n" + "=" * 80)
print("AGGREGATING φ BY (NVU_NODE, PT_BIN) USING ids_core")
print("=" * 80)

# Aggregate each module separately and save in long format
all_profiles = []

for module in available_modules:
    phi_col = f'phi_{module}'

    # Create a temporary DataFrame with required columns
    temp_df = df_als[['cell_id', 'nvu_node', 'PT_bin', phi_col]].copy()
    temp_df = temp_df.dropna(subset=[phi_col])

    print(f"  Module {module}: {len(temp_df)} cells with non-NaN φ")

    # Extract φ, group (nvu_node), and bin_index for this module
    phi = temp_df.set_index('cell_id')[phi_col]
    group = temp_df.set_index('cell_id')['nvu_node']
    bin_idx = temp_df.set_index('cell_id')['PT_bin']

    # Use ids_core to aggregate
    summary = ids_core.summarize_phi_by_group_and_bin(
        phi=phi,
        group=group,
        bin_index=bin_idx,
        group_order=['Vascular', 'Glia', 'Upper', 'VAT1L']
    )

    print(f"  → Aggregated to {len(summary)} (node, bin) combinations")

    # Add module name and PT bin center
    summary['module'] = module
    summary['pt_bin_center'] = summary['bin'].apply(lambda b: bin_centers[int(b)])

    # Rename columns to match original format
    summary = summary.rename(columns={
        'group': 'node',
        'bin': 'pt_bin'
    })

    # Reorder columns
    summary = summary[['node', 'module', 'pt_bin', 'pt_bin_center', 'phi_mean', 'phi_std', 'n_cells']]

    all_profiles.append(summary)

# Combine all profiles
if len(all_profiles) > 0:
    df_long = pd.concat(all_profiles, ignore_index=True)
    print(f"\n✓ Created {len(df_long)} aggregated (node, module, bin) profiles")
else:
    print(f"\nERROR: No profiles created!")
    sys.exit(1)

# Save profiles
output_file = f'{OUTPUT_DIR}/nvu_phi_profiles_final_by_node_and_bin.csv'
df_long.to_csv(output_file, index=False)
print(f"\n✓ Saved: {output_file}")

# ========================================================================
# DETERMINE FLOW ORDER USING ids_core.detect_onset_and_peak()
# ========================================================================

print("\n" + "=" * 80)
print("DETERMINING FLOW ORDER USING ids_core.detect_onset_and_peak()")
print("=" * 80)

flow_records = []

for module in available_modules:
    print(f"\nModule: {module}")

    # Get phi_mean profile for this module
    df_module = df_long[df_long['module'] == module].copy()

    # Prepare data for detect_onset_and_peak
    phi_mean_df = df_module[['node', 'pt_bin', 'phi_mean']].copy()
    phi_mean_df = phi_mean_df.rename(columns={'node': 'group', 'pt_bin': 'bin'})

    # Use ids_core to detect onset and peak
    flow_summary = ids_core.detect_onset_and_peak(
        phi_mean_df=phi_mean_df,
        bin_centers=bin_centers,
        threshold=PHI_THRESHOLD,
        min_bins=1
    )

    # Add module name
    flow_summary['module'] = module

    # Rename columns to match original format
    flow_summary = flow_summary.rename(columns={
        'group': 'node',
        'onset_PT': 'onset_pt',
        'peak_PT': 'peak_pt'
    })

    # Reorder columns
    flow_summary = flow_summary[['module', 'node', 'onset_pt', 'peak_pt', 'peak_phi']]

    flow_records.append(flow_summary)

    print(f"  ✓ Detected onset/peak for {len(flow_summary)} nodes")

# Combine all flow summaries
df_flow = pd.concat(flow_records, ignore_index=True)

# Compute ranks within each module (matching original implementation)
df_flow['onset_rank'] = df_flow.groupby('module')['onset_pt'].rank()
df_flow['peak_rank'] = df_flow.groupby('module')['peak_pt'].rank()

print(f"\n✓ Detected onset/peak for {len(df_flow)} (module, node) combinations")

# Save flow summary
flow_summary_file = f'{OUTPUT_DIR}/module_flow_summary_final.csv'
df_flow.to_csv(flow_summary_file, index=False)
print(f"✓ Saved: {flow_summary_file}")

# ========================================================================
# CREATE FLOW ORDER SUMMARY
# ========================================================================

print("\n" + "=" * 80)
print("CREATING FLOW ORDER SUMMARY")
print("=" * 80)

flow_order_records = []

for module in available_modules:
    df_mod = df_flow[df_flow['module'] == module].copy()

    # Sort by onset_rank
    df_mod_sorted = df_mod.sort_values('onset_rank')
    flow_onset = ' → '.join(df_mod_sorted['node'].values)

    # Sort by peak_rank
    df_mod_sorted = df_mod.sort_values('peak_rank')
    flow_peak = ' → '.join(df_mod_sorted['node'].values)

    flow_order_records.append({
        'module': module,
        'flow_order_onset': flow_onset,
        'flow_order_peak': flow_peak
    })

df_flow_order = pd.DataFrame(flow_order_records)
flow_order_file = f'{OUTPUT_DIR}/module_flow_order_final.csv'
df_flow_order.to_csv(flow_order_file, index=False)
print(f"✓ Saved: {flow_order_file}")

print(f"\n{df_flow_order.to_string(index=False)}")

# ========================================================================
# SUMMARY
# ========================================================================

print("\n" + "=" * 80)
print("REFACTORED PHASE 13c COMPLETE")
print("=" * 80)

print(f"\n✅ All core φ-flow operations performed using ids_core:")
print(f"   - compute_phi_from_scores(): {len(available_modules)} modules")
print(f"   - bin_by_pt(): {NUM_BINS} bins")
print(f"   - summarize_phi_by_group_and_bin(): {len(df_long)} profiles")
print(f"   - detect_onset_and_peak(): {len(df_flow)} flow records")

print(f"\n📁 Output files:")
print(f"   - {output_file}")
print(f"   - {flow_summary_file}")
print(f"   - {flow_order_file}")

print(f"\n🔬 Next step: Compare with original Phase13c outputs to verify identical results")
