#!/usr/bin/env python3
"""
GSE212630: Validate PT_dpt Hypothesis
=====================================

Hypothesis: PT_dpt shows low correlation with TDP-43 stages because:
1. TDP-43 pathology is a downstream convergence point
2. PT_dpt captures upstream functional changes
3. Upstream cells (Vascular/Glia) show different correlation patterns than neurons

Validation approach:
1. Compute PT_dpt vs TDP-43 stage correlation BY CELL TYPE
2. Compute module score vs TDP-43 stage correlation BY MODULE
3. Compare with existing project onset order

Author: Claude Code + User
Date: 2025-12-09
"""

import pandas as pd
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

print("=" * 80)
print("GSE212630: Validate PT_dpt Hypothesis")
print("=" * 80)

# ============================================================================
# CONFIGURATION
# ============================================================================

DATA_DIR = Path('/home/akaco/als/motor_cortex_analysis/ids_causal_analysis')
GSE_RESULTS = DATA_DIR / 'results' / 'GSE212630_ids_analysis'
OUTPUT_DIR = GSE_RESULTS / 'hypothesis_validation'
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# TDP-43 stage mapping
TDP_STAGE_TO_NUMERIC = {'Control': 0, 'TDPneg': 1, 'TDPmed': 2, 'TDPhigh': 3}

# NVU grouping
NVU_MAPPING = {
    'Excitatory': ['Ex.RORB', 'Ex.L2.3.RORB', 'Ex.THEMIS.UNC5C', 'Ex.OPRK1',
                   'Ex.L5.6.THEMIS', 'Ex.L3.5.RORB', 'Ex.L5.6.RORB'],
    'Inhibitory': ['Inh.VIP', 'Inh.PVALB', 'Inh.SST', 'Inh.LAMP5', 'Inh.PAX6'],
    'Glia': ['Astro', 'Oligo', 'OPC', 'Micro'],
    'Vascular': ['Endo', 'VLMC', 'Peri', 'Immune']
}

def get_nvu_group(cell_type):
    for group, types in NVU_MAPPING.items():
        for t in types:
            if t in cell_type:
                return group
    return 'Other'

# ============================================================================
# STEP 1: LOAD DATA
# ============================================================================

print("\n[Step 1] Loading data...")

# Load combined metadata with modules and PT
data_file = GSE_RESULTS / 'combined_metadata_with_modules_and_PT.csv'
if not data_file.exists():
    data_file = GSE_RESULTS / 'combined_metadata_with_all_PTs.csv'

df = pd.read_csv(data_file)
print(f"  Total cells: {len(df):,}")
print(f"  Columns: {list(df.columns)}")

# Add TDP-43 numeric stage
df['TDP_stage_numeric'] = df['condition'].map(TDP_STAGE_TO_NUMERIC)

# Add NVU group
df['NVU_group'] = df['cell_type'].apply(get_nvu_group)

print(f"\n  NVU groups: {df['NVU_group'].value_counts().to_dict()}")
print(f"  TDP stages: {df['condition'].value_counts().to_dict()}")

# ============================================================================
# STEP 2: PT_dpt vs TDP-43 CORRELATION BY CELL TYPE / NVU GROUP
# ============================================================================

print("\n[Step 2] Computing PT_dpt vs TDP-43 correlation by cell type...")

# Check which PT columns exist
pt_cols = [c for c in df.columns if c.startswith('PT_')]
print(f"  Available PT columns: {pt_cols}")

# Use PT_dpt if available, otherwise PT_dpt_module
pt_col = 'PT_dpt_module' if 'PT_dpt_module' in df.columns else 'PT_dpt'
if pt_col not in df.columns:
    print(f"  [ERROR] No PT_dpt column found!")
    exit(1)

print(f"  Using: {pt_col}")

# Correlation by NVU group
print("\n  === PT_dpt vs TDP-43 Stage Correlation by NVU Group ===")
results_nvu = []

for group in ['Vascular', 'Glia', 'Excitatory', 'Inhibitory']:
    subset = df[df['NVU_group'] == group]
    if len(subset) < 10:
        continue

    # Spearman correlation
    rho, pval = stats.spearmanr(subset[pt_col], subset['TDP_stage_numeric'])

    results_nvu.append({
        'NVU_group': group,
        'n_cells': len(subset),
        'rho': rho,
        'p_value': pval,
        'significant': pval < 0.05
    })

    print(f"    {group}: ρ = {rho:.4f}, p = {pval:.2e}, n = {len(subset):,}")

results_nvu_df = pd.DataFrame(results_nvu)

# Correlation by individual cell type
print("\n  === PT_dpt vs TDP-43 Stage Correlation by Cell Type ===")
results_ct = []

for ct in df['cell_type'].unique():
    subset = df[df['cell_type'] == ct]
    if len(subset) < 30:
        continue

    rho, pval = stats.spearmanr(subset[pt_col], subset['TDP_stage_numeric'])
    nvu = get_nvu_group(ct)

    results_ct.append({
        'cell_type': ct,
        'NVU_group': nvu,
        'n_cells': len(subset),
        'rho': rho,
        'p_value': pval
    })

results_ct_df = pd.DataFrame(results_ct).sort_values('rho', ascending=False)
print(results_ct_df.to_string())

# ============================================================================
# STEP 3: MODULE SCORE vs TDP-43 CORRELATION
# ============================================================================

print("\n[Step 3] Computing module score vs TDP-43 correlation...")

# Find module columns
module_cols = [c for c in df.columns if c.startswith('module_')]
print(f"  Found {len(module_cols)} module columns")

if len(module_cols) == 0:
    print("  [WARN] No module columns found, skipping module analysis")
    results_module_df = pd.DataFrame()
else:
    results_module = []

    for mod_col in module_cols:
        module_name = mod_col.replace('module_', '')

        # Overall correlation
        rho, pval = stats.spearmanr(df[mod_col].dropna(),
                                     df.loc[df[mod_col].notna(), 'TDP_stage_numeric'])

        results_module.append({
            'module': module_name,
            'rho': rho,
            'p_value': pval,
            'significant': pval < 0.05
        })

    results_module_df = pd.DataFrame(results_module).sort_values('rho', ascending=False)

    print("\n  === Module vs TDP-43 Stage Correlation (Top 10 positive) ===")
    print(results_module_df.head(10).to_string())

    print("\n  === Module vs TDP-43 Stage Correlation (Top 10 negative) ===")
    print(results_module_df.tail(10).to_string())

# ============================================================================
# STEP 4: COMPARE WITH EXISTING PROJECT ONSET ORDER
# ============================================================================

print("\n[Step 4] Comparing with existing project onset order...")

# Load existing project onset data
existing_onset_file = DATA_DIR / 'results' / 'phase13_nvu_energy' / 'module_flow_summary_final.csv'
if existing_onset_file.exists():
    existing_df = pd.read_csv(existing_onset_file)

    # Get mean onset PT by module
    existing_onset = existing_df.groupby('module')['onset_pt'].mean().sort_values()

    print("\n  === Existing Project: Module Onset Order (earliest first) ===")
    print(existing_onset.head(10))

    # Compare: early onset modules should have LOW correlation with TDP-43 stages
    if len(results_module_df) > 0:
        # Merge
        comparison = results_module_df.merge(
            existing_onset.reset_index().rename(columns={'onset_pt': 'existing_onset_pt'}),
            left_on='module',
            right_on='module',
            how='inner'
        )

        if len(comparison) > 0:
            # Correlation between onset order and TDP-43 correlation
            rho_onset_tdp, pval = stats.spearmanr(comparison['existing_onset_pt'],
                                                   comparison['rho'])

            print(f"\n  === Key Finding ===")
            print(f"  Correlation between existing onset PT and TDP-43 correlation: ρ = {rho_onset_tdp:.4f}, p = {pval:.4f}")
            print(f"  Interpretation: {'Early onset modules show LOWER TDP-43 correlation (supports hypothesis)' if rho_onset_tdp > 0 else 'No clear pattern'}")

            comparison.to_csv(OUTPUT_DIR / 'module_onset_vs_TDP_correlation.csv', index=False)
else:
    print("  [WARN] Existing onset file not found")

# ============================================================================
# STEP 5: GENERATE FIGURES
# ============================================================================

print("\n[Step 5] Generating figures...")

# Figure 1: PT_dpt vs TDP-43 correlation by NVU group
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# Bar plot of correlations
ax = axes[0]
colors = {'Vascular': 'red', 'Glia': 'orange', 'Excitatory': 'green', 'Inhibitory': 'blue'}
if len(results_nvu_df) > 0:
    bars = ax.bar(results_nvu_df['NVU_group'], results_nvu_df['rho'],
                  color=[colors.get(g, 'gray') for g in results_nvu_df['NVU_group']])
    ax.axhline(y=0, color='black', linestyle='-', linewidth=0.5)
    ax.set_xlabel('NVU Group', fontsize=12)
    ax.set_ylabel('Spearman ρ (PT_dpt vs TDP-43 Stage)', fontsize=12)
    ax.set_title('PT_dpt vs TDP-43 Stage Correlation by NVU Group', fontsize=12)
    ax.grid(True, alpha=0.3, axis='y')

    # Add significance markers
    for i, (_, row) in enumerate(results_nvu_df.iterrows()):
        if row['p_value'] < 0.001:
            ax.text(i, row['rho'] + 0.01, '***', ha='center', fontsize=10)
        elif row['p_value'] < 0.01:
            ax.text(i, row['rho'] + 0.01, '**', ha='center', fontsize=10)
        elif row['p_value'] < 0.05:
            ax.text(i, row['rho'] + 0.01, '*', ha='center', fontsize=10)

# Scatter plot by cell type
ax = axes[1]
if len(results_ct_df) > 0:
    for nvu, color in colors.items():
        subset = results_ct_df[results_ct_df['NVU_group'] == nvu]
        ax.scatter(subset['n_cells'], subset['rho'], c=color, label=nvu, alpha=0.7, s=50)
    ax.axhline(y=0, color='black', linestyle='-', linewidth=0.5)
    ax.set_xlabel('Number of Cells', fontsize=12)
    ax.set_ylabel('Spearman ρ (PT_dpt vs TDP-43 Stage)', fontsize=12)
    ax.set_title('PT_dpt vs TDP-43 Correlation by Cell Type', fontsize=12)
    ax.legend()
    ax.grid(True, alpha=0.3)

plt.tight_layout()
fig.savefig(OUTPUT_DIR / 'Fig_PT_dpt_TDP_correlation_by_NVU.png', dpi=150, bbox_inches='tight')
plt.close()
print(f"  Saved: Fig_PT_dpt_TDP_correlation_by_NVU.png")

# Figure 2: Module correlation with TDP-43 stage
if len(results_module_df) > 0:
    fig, ax = plt.subplots(figsize=(12, 8))

    # Sort by correlation
    sorted_df = results_module_df.sort_values('rho')

    # Color by significance
    colors = ['red' if p < 0.05 and r < 0 else 'green' if p < 0.05 and r > 0 else 'gray'
              for r, p in zip(sorted_df['rho'], sorted_df['p_value'])]

    ax.barh(sorted_df['module'], sorted_df['rho'], color=colors)
    ax.axvline(x=0, color='black', linestyle='-', linewidth=0.5)
    ax.set_xlabel('Spearman ρ (Module Score vs TDP-43 Stage)', fontsize=12)
    ax.set_title('Module Score vs TDP-43 Stage Correlation\n(Green=positive sig., Red=negative sig., Gray=n.s.)', fontsize=12)
    ax.grid(True, alpha=0.3, axis='x')

    plt.tight_layout()
    fig.savefig(OUTPUT_DIR / 'Fig_module_TDP_correlation.png', dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: Fig_module_TDP_correlation.png")

# Save results
results_nvu_df.to_csv(OUTPUT_DIR / 'PT_dpt_TDP_correlation_by_NVU.csv', index=False)
results_ct_df.to_csv(OUTPUT_DIR / 'PT_dpt_TDP_correlation_by_celltype.csv', index=False)
if len(results_module_df) > 0:
    results_module_df.to_csv(OUTPUT_DIR / 'module_TDP_correlation.csv', index=False)

# ============================================================================
# SUMMARY
# ============================================================================

print("\n" + "=" * 80)
print("HYPOTHESIS VALIDATION SUMMARY")
print("=" * 80)

print("""
HYPOTHESIS: PT_dpt shows low overall correlation with TDP-43 stages because
TDP-43 pathology is a downstream convergence point, while PT_dpt captures
upstream functional changes.

PREDICTIONS:
1. Neurons (downstream) should show HIGHER PT_dpt vs TDP-43 correlation
2. Vascular/Glia (upstream) should show LOWER or NEGATIVE correlation
3. Upstream modules (early onset) should show lower TDP-43 correlation

RESULTS:
""")

print("1. PT_dpt vs TDP-43 Correlation by NVU Group:")
for _, row in results_nvu_df.iterrows():
    print(f"   {row['NVU_group']}: ρ = {row['rho']:.4f} {'(sig)' if row['significant'] else '(n.s.)'}")

if len(results_module_df) > 0:
    top_pos = results_module_df.head(3)
    top_neg = results_module_df.tail(3)
    print("\n2. Modules with HIGHEST TDP-43 correlation (downstream?):")
    for _, row in top_pos.iterrows():
        print(f"   {row['module']}: ρ = {row['rho']:.4f}")
    print("\n3. Modules with LOWEST TDP-43 correlation (upstream?):")
    for _, row in top_neg.iterrows():
        print(f"   {row['module']}: ρ = {row['rho']:.4f}")

print(f"\nOutput: {OUTPUT_DIR}")
print("=" * 80)
