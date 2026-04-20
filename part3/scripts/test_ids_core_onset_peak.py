#!/usr/bin/env python3
"""
Test: ids_core onset/peak detection
====================================
Verify that detect_onset_and_peak() produces identical results
to Phase13c_NVU_energy_complete.py
"""

import pandas as pd
import numpy as np
from pathlib import Path
import sys

# Import ids_core
sys.path.insert(0, str(Path(__file__).parent))
import ids_core

print("=" * 80)
print("TEST: ids_core onset/peak detection")
print("=" * 80)

# Configuration
PROFILES_FILE = 'results/phase13_nvu_energy/nvu_phi_profiles_by_node_and_bin.csv'
REFERENCE_FILE = 'results/phase13_nvu_energy/module_flow_summary.csv'
TEST_MODULE = 'Mitochondria'
PHI_THRESHOLD = 0.5

# ========================================================================
# LOAD DATA
# ========================================================================

print(f"\nLoading aggregated profiles: {PROFILES_FILE}")
if not Path(PROFILES_FILE).exists():
    print(f"ERROR: {PROFILES_FILE} not found")
    sys.exit(1)

profiles = pd.read_csv(PROFILES_FILE)
print(f"Loaded {len(profiles)} (node, module, bin) combinations")

# Filter to test module
df_module = profiles[profiles['module'] == TEST_MODULE].copy()
print(f"Filtering to {TEST_MODULE}: {len(df_module)} rows")

# ========================================================================
# EXTRACT BIN CENTERS
# ========================================================================

# Extract bin centers from the aggregated data
# (In real usage, this would come from bin_by_pt() output)
n_bins = df_module['pt_bin'].nunique()
bin_centers = df_module.groupby('pt_bin')['pt_bin_center'].first().sort_index().values

print(f"\nExtracted {n_bins} bin centers")
print(f"  PT range: [{bin_centers[0]:.4f}, {bin_centers[-1]:.4f}]")

# ========================================================================
# DETECT ONSET AND PEAK USING ids_core
# ========================================================================

print(f"\nDetecting onset and peak using ids_core...")

# Prepare data in the format expected by detect_onset_and_peak()
# The function expects: phi_mean_df with columns ['group', 'bin', 'phi_mean']
phi_mean_df = df_module[['node', 'pt_bin', 'phi_mean']].copy()
phi_mean_df = phi_mean_df.rename(columns={'node': 'group', 'pt_bin': 'bin'})

# Call ids_core function
flow_summary = ids_core.detect_onset_and_peak(
    phi_mean_df=phi_mean_df,
    bin_centers=bin_centers,
    threshold=PHI_THRESHOLD,
    min_bins=1
)

print(f"Detected onset/peak for {len(flow_summary)} groups:")
print(flow_summary[['group', 'onset_PT', 'peak_PT', 'peak_phi']])

# ========================================================================
# LOAD REFERENCE DATA
# ========================================================================

print(f"\nLoading reference: {REFERENCE_FILE}")
if not Path(REFERENCE_FILE).exists():
    print(f"ERROR: {REFERENCE_FILE} not found")
    sys.exit(1)

ref = pd.read_csv(REFERENCE_FILE)
ref_module = ref[ref['module'] == TEST_MODULE].copy()

print(f"Loaded reference with {len(ref_module)} groups for {TEST_MODULE}")

# ========================================================================
# COMPARE
# ========================================================================

print("\n" + "=" * 80)
print("COMPARISON")
print("=" * 80)

# Merge on group (node)
merged = flow_summary.merge(
    ref_module[['node', 'onset_pt', 'peak_pt', 'peak_phi']],
    left_on='group',
    right_on='node',
    how='inner',
    suffixes=('_ids_core', '_original')
)

print(f"\nMatched {len(merged)} groups")

# Compare onset_PT
onset_ids_core = merged['onset_PT']
onset_original = merged['onset_pt']

# Handle NaN for onset (some groups may not exceed threshold)
onset_diff = (onset_ids_core - onset_original).abs()
onset_diff_valid = onset_diff.dropna()

print(f"\n  Onset PT comparison:")
print(f"    ids_core NaN count: {onset_ids_core.isna().sum()}")
print(f"    original NaN count: {onset_original.isna().sum()}")
if len(onset_diff_valid) > 0:
    print(f"    Mean abs diff (non-NaN): {onset_diff_valid.mean():.2e}")
    print(f"    Max abs diff (non-NaN):  {onset_diff_valid.max():.2e}")

# Compare peak_PT
peak_pt_ids_core = merged['peak_PT']
peak_pt_original = merged['peak_pt']

peak_pt_diff = (peak_pt_ids_core - peak_pt_original).abs()

print(f"\n  Peak PT comparison:")
print(f"    Mean abs diff: {peak_pt_diff.mean():.2e}")
print(f"    Max abs diff:  {peak_pt_diff.max():.2e}")

# Compare peak_phi
peak_phi_ids_core = merged['peak_phi_ids_core']
peak_phi_original = merged['peak_phi_original']

peak_phi_diff = (peak_phi_ids_core - peak_phi_original).abs()

print(f"\n  Peak φ comparison:")
print(f"    Mean abs diff: {peak_phi_diff.mean():.2e}")
print(f"    Max abs diff:  {peak_phi_diff.max():.2e}")

# Correlation
corr_peak_pt = peak_pt_ids_core.corr(peak_pt_original)
corr_peak_phi = peak_phi_ids_core.corr(peak_phi_original)

print(f"\n  Correlations:")
print(f"    Peak PT:  r = {corr_peak_pt:.10f}")
print(f"    Peak φ:   r = {corr_peak_phi:.10f}")

# ========================================================================
# DETAILED COMPARISON
# ========================================================================

print("\n  Detailed comparison:")
print(merged[['group', 'onset_PT', 'onset_pt', 'peak_PT', 'peak_pt', 'peak_phi_ids_core', 'peak_phi_original']].to_string(index=False))

# ========================================================================
# VERDICT
# ========================================================================

print("\n" + "=" * 80)
print("VERDICT")
print("=" * 80)

# Check if onset NaN patterns match
onset_nan_match = (onset_ids_core.isna() == onset_original.isna()).all()

# Check numerical agreement
onset_ok = len(onset_diff_valid) == 0 or (onset_diff_valid.max() < 1e-10)
peak_pt_ok = peak_pt_diff.max() < 1e-10
peak_phi_ok = peak_phi_diff.max() < 1e-10

if onset_nan_match and onset_ok and peak_pt_ok and peak_phi_ok:
    print("\n✅ TEST PASSED")
    print(f"   ids_core.detect_onset_and_peak() produces identical results.")
    print(f"   Onset NaN patterns match: {onset_nan_match}")
    if len(onset_diff_valid) > 0:
        print(f"   Max onset PT diff:  {onset_diff_valid.max():.2e} < 1e-10")
    print(f"   Max peak PT diff:   {peak_pt_diff.max():.2e} < 1e-10")
    print(f"   Max peak φ diff:    {peak_phi_diff.max():.2e} < 1e-10")
    sys.exit(0)
else:
    print("\n❌ TEST FAILED")
    print(f"   Results differ from original implementation.")
    print(f"   Onset NaN patterns match: {onset_nan_match}")
    if len(onset_diff_valid) > 0:
        print(f"   Max onset PT diff:  {onset_diff_valid.max():.2e}")
    print(f"   Max peak PT diff:   {peak_pt_diff.max():.2e}")
    print(f"   Max peak φ diff:    {peak_phi_diff.max():.2e}")
    sys.exit(1)
