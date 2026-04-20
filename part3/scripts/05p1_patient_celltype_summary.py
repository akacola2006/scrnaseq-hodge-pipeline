#!/usr/bin/env python3
"""
Phase 5′-1: Patient × Cell Type Summary
==========================================

Creates patient-level summaries to identify disease progression patterns
and potential clinical subtypes.

Output:
- patient_celltype_summary.csv: Patient × Cell type aggregated features
- patient_progression_matrix.csv: Patient × Key cell types PT matrix
"""

import pandas as pd
import numpy as np
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# Paths
base_dir = Path('/home/akaco/als/motor_cortex_analysis')
ids_dir = base_dir / 'ids_causal_analysis'
results_dir = ids_dir / 'results' / 'patient_stratified'
cell_state_dir = ids_dir / 'results' / 'cell_state_causality'

print("="*80)
print("Phase 5′-1: Patient × Cell Type Summary")
print("="*80)

# ============================================================================
# 1. Load Cell-Level Features
# ============================================================================
print("\n[1] Loading cell-level features...")

cell_features_file = cell_state_dir / 'cell_level_features_ALL.csv'
df = pd.read_csv(cell_features_file)

print(f"  Loaded: {len(df):,} cells")
print(f"  Columns: {len(df.columns)}")

# ============================================================================
# 2. Extract Patient ID from Cell ID
# ============================================================================
print("\n[2] Extracting patient IDs...")

# Extract patient_id from cell_id (format: BARCODE-PATIENTID)
df['patient_id'] = df['cell_id'].str.split('-').str[1]

n_patients = df['patient_id'].nunique()
n_cell_types = df['cell_type'].nunique()

print(f"  Patients: {n_patients}")
print(f"  Cell types: {n_cell_types}")

# Count cells per patient
patient_counts = df.groupby('patient_id').size().sort_values(ascending=False)
print(f"\n  Top 5 patients by cell count:")
for patient, count in patient_counts.head(5).items():
    condition = df[df['patient_id'] == patient]['condition'].iloc[0]
    print(f"    {patient}: {count:,} cells ({condition})")

# ============================================================================
# 3. Identify Module Columns
# ============================================================================
print("\n[3] Identifying module columns...")

module_cols = [col for col in df.columns if col.startswith('module_')]
print(f"  Found {len(module_cols)} modules")

# Key modules for analysis
key_modules = [
    'module_Angiogenesis',
    'module_Synaptic',
    'module_Mitochondria',
    'module_Apoptosis',
    'module_ER_Stress',
    'module_Oxidative_Stress',
    'module_Protein_Homeostasis',
    'module_Inflammation'
]
key_modules = [m for m in key_modules if m in module_cols]
print(f"  Key modules for summary: {len(key_modules)}")

# ============================================================================
# 4. Create Patient × Cell Type Summary
# ============================================================================
print("\n[4] Creating patient × cell type summary...")

# Group by patient and cell type
grouped = df.groupby(['patient_id', 'cell_type'])

# Calculate aggregated metrics
summary_list = []

for (patient, celltype), group in grouped:
    n_cells = len(group)

    summary = {
        'patient_id': patient,
        'cell_type': celltype,
        'condition': group['condition'].iloc[0],
        'n_cells': n_cells,
        'PT_IDS_median': group['PT_IDS'].median(),
        'PT_IDS_mean': group['PT_IDS'].mean(),
        'PT_IDS_std': group['PT_IDS'].std(),
        'stress_total_mean': group['stress_total'].mean(),
        'stress_total_median': group['stress_total'].median(),
        'stress_total_std': group['stress_total'].std()
    }

    # Add key module means
    for module in key_modules:
        module_name = module.replace('module_', '')
        summary[f'{module_name}_mean'] = group[module].mean()

    summary_list.append(summary)

# Create summary DataFrame
summary_df = pd.DataFrame(summary_list)

print(f"  Created summary: {len(summary_df)} combinations")
print(f"    Patients: {summary_df['patient_id'].nunique()}")
print(f"    Cell types: {summary_df['cell_type'].nunique()}")

# Save patient × cell type summary
output_file = results_dir / 'patient_celltype_summary.csv'
summary_df.to_csv(output_file, index=False)
print(f"\n  Saved: {output_file}")

# ============================================================================
# 5. Create Patient Progression Matrix
# ============================================================================
print("\n[5] Creating patient progression matrix...")

# Select key cell types for progression matrix
key_cell_types = [
    'Ex.L5.VAT1L.EYA4',
    'Ex.L5.VAT1L.THSD4',
    'Glia.Oligo',
    'Glia.Astro.GFAP-neg',
    'Glia.Astro.GFAP-pos',
    'Glia.Micro',
    'Ex.L2.L3.CUX2.RASGRF2',
    'Ex.L3.L5.CUX2.RORB',
    'Ex.L4.L6.RORB.LRRK1',
    'Ex.L5.L6.THEMIS.TMEM233',
    'In.SOM.SST.GALNT14',
    'In.PV.PVALB.CEMIP'
]

# Filter summary for key cell types
key_summary = summary_df[summary_df['cell_type'].isin(key_cell_types)].copy()

# Pivot to create patient × cell type matrix
progression_matrix = key_summary.pivot_table(
    index='patient_id',
    columns='cell_type',
    values='PT_IDS_median',
    aggfunc='first'
)

# Add overall patient metrics
patient_overall = df.groupby('patient_id').agg({
    'PT_IDS': 'median',
    'stress_total': 'mean',
    'condition': 'first',
    'cell_id': 'count'
}).rename(columns={
    'PT_IDS': 'PT_IDS_global_median',
    'stress_total': 'stress_global_mean',
    'cell_id': 'total_cells'
})

# Merge with progression matrix
progression_matrix = progression_matrix.join(patient_overall)

# Reorder columns: metadata first, then cell types
metadata_cols = ['condition', 'total_cells', 'PT_IDS_global_median', 'stress_global_mean']
celltype_cols = [col for col in progression_matrix.columns if col not in metadata_cols]
progression_matrix = progression_matrix[metadata_cols + celltype_cols]

print(f"  Progression matrix: {progression_matrix.shape}")
print(f"    Patients: {len(progression_matrix)}")
print(f"    Cell types: {len(celltype_cols)}")

# Save progression matrix
output_file = results_dir / 'patient_progression_matrix.csv'
progression_matrix.to_csv(output_file)
print(f"\n  Saved: {output_file}")

# ============================================================================
# 6. Quick Summary Statistics
# ============================================================================
print("\n[6] Summary statistics:")

# Overall patient stats
als_patients = summary_df[summary_df['condition'] == 'ALS']['patient_id'].unique()
ctrl_patients = summary_df[summary_df['condition'] == 'Control']['patient_id'].unique()

print(f"\n  ALS patients: {len(als_patients)}")
print(f"  Control patients: {len(ctrl_patients)}")

# PT distribution by condition
print(f"\n  Global PT_IDS distribution:")
for condition in ['ALS', 'Control']:
    pt_vals = df[df['condition'] == condition]['PT_IDS']
    print(f"    {condition}: median={pt_vals.median():.3f}, "
          f"mean={pt_vals.mean():.3f}, "
          f"std={pt_vals.std():.3f}")

# Cell type coverage across patients
print(f"\n  Cell type coverage:")
celltype_coverage = summary_df.groupby('cell_type')['patient_id'].nunique().sort_values(ascending=False)
print(f"    Top 5 cell types by patient coverage:")
for celltype, n_patients in celltype_coverage.head(5).items():
    print(f"      {celltype}: {n_patients}/{n_patients} patients")

print("\n" + "="*80)
print("Phase 5′-1 completed!")
print(f"Output files in: {results_dir}")
print("="*80)
