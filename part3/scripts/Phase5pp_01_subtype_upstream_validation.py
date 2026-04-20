#!/usr/bin/env python3
"""
Phase 5″ (Double-Prime): Subtype-Specific Upstream/Downstream Validation

Purpose:
- Use FIXED subtype labels (from Phase 5' IMES-based clustering)
- Re-evaluate upstream/downstream relationships WITHIN each subtype
- Use PT_dpt (stress-independent) instead of PT_imes
- Validate if subtype-specific causal patterns are robust

Key Innovation:
- Subtype definition: IMES-based (fixed)
- Ordering validation: PT_dpt-based (independent)
- This separates "subtype discovery" from "ordering validation"
"""

import pandas as pd
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

print("="*80)
print("Phase 5″: Subtype-Specific Upstream/Downstream Validation")
print("="*80)

# Paths
base_dir = Path('/home/akaco/als/motor_cortex_analysis')
ids_dir = base_dir / 'ids_causal_analysis'
output_dir = ids_dir / 'results' / 'phase5pp_subtype_validation'
output_dir.mkdir(parents=True, exist_ok=True)

# ============================================================================
# Step 1: Load and Merge Data
# ============================================================================
print("\n[Step 1] Loading data and merging subtype labels...")

# Load cell-level data (use quick test data for now)
cell_data_file = ids_dir / 'results' / 'PTv2_robustness' / 'PTv2_quick_key_celltypes_with_PTs.csv'

if not cell_data_file.exists():
    # Fallback to full dataset
    cell_data_file = ids_dir / 'results' / 'cell_state_causality' / 'cell_level_features_ALL.csv'
    print(f"  Using full dataset: {cell_data_file.name}")
else:
    print(f"  Using quick test data: {cell_data_file.name}")

df = pd.read_csv(cell_data_file)
df = df.rename(columns={'PT_IDS': 'PT_imes'})

# Extract patient_id if not present
if 'patient_id' not in df.columns:
    df['patient_id'] = df['cell_id'].str.split('-').str[1]

print(f"  Loaded {len(df):,} cells")
print(f"  Columns: {df.columns.tolist()}")

# Load patient clusters (subtype labels)
clusters = pd.read_csv(ids_dir / 'results' / 'patient_stratified' / 'patient_clusters.csv')
subtype_map = {
    0: 'Upper-layer',
    1: 'Pure Oligo',
    2: 'Oligo-Inflammation',
    3: 'Upper-layer',
    4: 'Pure Oligo',
    5: 'Oligo-Inflammation'
}
clusters['subtype'] = clusters['cluster_id'].map(subtype_map)

print(f"\n  Patient subtypes:")
print(clusters[['patient_id', 'cluster_id', 'subtype']].to_string(index=False))

# Merge (check if subtype already exists)
if 'subtype' in df.columns:
    print(f"\n  Subtype already in data, skipping merge")
    print(f"  After merge: {len(df):,} cells")
    print(f"  Cells with subtype: {df['subtype'].notna().sum():,}")
else:
    df = df.merge(clusters[['patient_id', 'subtype']], on='patient_id', how='left')
    print(f"\n  After merge: {len(df):,} cells")
    print(f"  Cells with subtype: {df['subtype'].notna().sum():,}")

# Filter to ALS only
df_als = df[df['condition'] == 'ALS'].copy()
print(f"  ALS cells: {len(df_als):,}")

# Check PT columns
pt_cols = [col for col in df_als.columns if col.startswith('PT_')]
print(f"\n  Available PT columns: {pt_cols}")

# Ensure we have PT_dpt
if 'PT_dpt' not in df_als.columns:
    print("\n  ERROR: PT_dpt not found! Need to run PTv2 pipeline first.")
    exit(1)

# ============================================================================
# Step 2: Define Key Cell Types
# ============================================================================
print("\n[Step 2] Defining key cell types...")

key_cell_types = [
    'Ex.L5.VAT1L.EYA4',
    'Ex.L5.VAT1L.THSD4',
    'Glia.Oligo',
    'Glia.Astro.GFAP-neg',
    'Glia.Micro',
    'Ex.L2.L3.CUX2.RASGRF2'
]

# Check availability
available_cts = [ct for ct in key_cell_types if ct in df_als['cell_type'].values]
print(f"  Requested: {len(key_cell_types)}")
print(f"  Available: {len(available_cts)}")
for ct in available_cts:
    n = len(df_als[df_als['cell_type'] == ct])
    print(f"    {ct}: {n:,} cells")

# Filter to key cell types
df_key = df_als[df_als['cell_type'].isin(available_cts)].copy()
print(f"\n  Key cell types subset: {len(df_key):,} cells")

# Save per-subtype data
for subtype in df_key['subtype'].unique():
    if pd.notna(subtype):
        subtype_df = df_key[df_key['subtype'] == subtype]
        safe_name = subtype.replace(' ', '_').replace('-', '_')
        output_file = output_dir / f'subtype_{safe_name}_cells.csv'
        subtype_df.to_csv(output_file, index=False)
        print(f"  Saved: {output_file.name} ({len(subtype_df):,} cells)")

# ============================================================================
# Step 3: Calculate PT Medians by Subtype and Cell Type
# ============================================================================
print("\n[Step 3] Calculating PT medians by subtype and cell type...")

PT_versions = ['PT_imes']
if 'PT_dpt' in df_key.columns:
    PT_versions.append('PT_dpt')
if 'PT_stress' in df_key.columns:
    PT_versions.append('PT_stress')

print(f"  PT versions: {PT_versions}")

# Calculate for each subtype
all_results = []

for subtype in ['Pure Oligo', 'Oligo-Inflammation', 'Upper-layer']:
    subtype_cells = df_key[df_key['subtype'] == subtype]

    if len(subtype_cells) == 0:
        continue

    print(f"\n  {subtype} ({len(subtype_cells):,} cells):")

    for ct in available_cts:
        ct_cells = subtype_cells[subtype_cells['cell_type'] == ct]

        if len(ct_cells) >= 5:  # Minimum threshold
            result = {
                'subtype': subtype,
                'cell_type': ct,
                'n_cells': len(ct_cells)
            }

            for pt in PT_versions:
                result[f'median_{pt}'] = ct_cells[pt].median()
                result[f'mean_{pt}'] = ct_cells[pt].mean()
                result[f'std_{pt}'] = ct_cells[pt].std()

            all_results.append(result)

            # Print summary
            pt_str = ', '.join([f"{pt}={result[f'median_{pt}']:.3f}"
                               for pt in PT_versions])
            print(f"    {ct}: {pt_str} (n={len(ct_cells)})")

# Save comprehensive summary
summary_df = pd.DataFrame(all_results)
summary_df.to_csv(output_dir / 'subtype_celltype_PT_summary_all.csv', index=False)
print(f"\n  Saved: subtype_celltype_PT_summary_all.csv")

# ============================================================================
# Step 4: Pairwise Ordering Tests
# ============================================================================
print("\n[Step 4] Pairwise ordering tests within subtypes...")

# Define key pairs to test
pairs = [
    ('Ex.L2.L3.CUX2.RASGRF2', 'Ex.L5.VAT1L.EYA4', 'Upper vs VAT1L'),
    ('Ex.L2.L3.CUX2.RASGRF2', 'Ex.L5.VAT1L.THSD4', 'Upper vs VAT1L'),
    ('Glia.Oligo', 'Ex.L5.VAT1L.EYA4', 'Oligo vs VAT1L'),
    ('Glia.Oligo', 'Ex.L5.VAT1L.THSD4', 'Oligo vs VAT1L'),
    ('Glia.Micro', 'Ex.L5.VAT1L.EYA4', 'Micro vs VAT1L'),
    ('Glia.Micro', 'Ex.L5.VAT1L.THSD4', 'Micro vs VAT1L'),
]

pairwise_results = []

for subtype in ['Pure Oligo', 'Oligo-Inflammation', 'Upper-layer']:
    subtype_cells = df_key[df_key['subtype'] == subtype]

    if len(subtype_cells) == 0:
        continue

    print(f"\n  {subtype}:")

    for ct1, ct2, pair_name in pairs:
        # Check if both cell types exist
        cells1 = subtype_cells[subtype_cells['cell_type'] == ct1]
        cells2 = subtype_cells[subtype_cells['cell_type'] == ct2]

        if len(cells1) < 5 or len(cells2) < 5:
            continue

        # Test for each PT version
        for pt in PT_versions:
            # t-test
            t, p = stats.ttest_ind(cells1[pt], cells2[pt])

            # Effect size (Cohen's d)
            pooled_std = np.sqrt((cells1[pt].std()**2 + cells2[pt].std()**2) / 2)
            cohen_d = (cells1[pt].mean() - cells2[pt].mean()) / pooled_std

            # Median difference
            median_diff = cells1[pt].median() - cells2[pt].median()

            # Ordering
            if median_diff < 0:
                ordering = f"{ct1.split('.')[-1]} < {ct2.split('.')[-1]}"
            else:
                ordering = f"{ct1.split('.')[-1]} > {ct2.split('.')[-1]}"

            pairwise_results.append({
                'subtype': subtype,
                'pair_name': pair_name,
                'cell_type_1': ct1,
                'cell_type_2': ct2,
                'PT_version': pt,
                'median_diff': median_diff,
                'mean_diff': cells1[pt].mean() - cells2[pt].mean(),
                'cohen_d': cohen_d,
                't_statistic': t,
                'p_value': p,
                'ordering': ordering,
                'n1': len(cells1),
                'n2': len(cells2)
            })

            sig = '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else 'ns'
            print(f"    {pair_name:20s} ({pt:12s}): {ordering:30s} p={p:.4f} {sig}")

pairwise_df = pd.DataFrame(pairwise_results)
pairwise_df.to_csv(output_dir / 'subtype_pairwise_order_tests.csv', index=False)
print(f"\n  Saved: subtype_pairwise_order_tests.csv")

# ============================================================================
# Step 5: Cross-Subtype Comparison
# ============================================================================
print("\n[Step 5] Cross-subtype comparison...")

print("\n  VAT1L position across subtypes (using PT_dpt):")

if 'PT_dpt' in PT_versions:
    for ct in ['Ex.L5.VAT1L.EYA4', 'Ex.L5.VAT1L.THSD4']:
        if ct not in available_cts:
            continue

        print(f"\n  {ct}:")

        for subtype in ['Pure Oligo', 'Oligo-Inflammation', 'Upper-layer']:
            subtype_vat1l = df_key[
                (df_key['subtype'] == subtype) &
                (df_key['cell_type'] == ct)
            ]

            if len(subtype_vat1l) >= 5:
                median_pt = subtype_vat1l['PT_dpt'].median()
                mean_pt = subtype_vat1l['PT_dpt'].mean()
                print(f"    {subtype:20s}: PT_dpt = {median_pt:.4f} (mean={mean_pt:.4f}, n={len(subtype_vat1l)})")

# ============================================================================
# Summary
# ============================================================================
print("\n" + "="*80)
print("Phase 5″ Step 1-3 Completed!")
print("="*80)
print(f"\nGenerated files:")
print(f"  - subtype_celltype_PT_summary_all.csv")
print(f"  - subtype_pairwise_order_tests.csv")
print(f"  - subtype_*_cells.csv (per-subtype data)")
print(f"\nNext steps:")
print(f"  - Visualization")
print(f"  - State-level validation (optional)")
print(f"  - Comprehensive report")
