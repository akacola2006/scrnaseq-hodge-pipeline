#!/usr/bin/env python3
"""
GSE212630: Upstream Module Analysis
====================================

Check if upstream modules (Ion_Transport, DNA_Repair, Calcium_Signaling)
show early changes specifically in Vascular cells.

Key questions:
1. Which NVU group shows the earliest change in upstream modules?
2. Is there a cell-type specific pattern for upstream vs downstream modules?
3. Can we trace the causal flow?

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
print("GSE212630: Upstream Module Analysis")
print("=" * 80)

# ============================================================================
# CONFIGURATION
# ============================================================================

DATA_DIR = Path('/home/akaco/als/motor_cortex_analysis/ids_causal_analysis')
GSE_RESULTS = DATA_DIR / 'results' / 'GSE212630_ids_analysis'
OUTPUT_DIR = GSE_RESULTS / 'upstream_analysis'
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# TDP-43 stage mapping
TDP_STAGE_TO_NUMERIC = {'Control': 0, 'TDPneg': 1, 'TDPmed': 2, 'TDPhigh': 3}

# NVU grouping - more inclusive
NVU_MAPPING = {
    'Excitatory': ['Ex.L2/3', 'Ex.L4', 'Ex.L5', 'Ex.L6', 'Ex.RORB', 'Ex.THEMIS',
                   'Ex.OPRK1', 'Ex.Unknown'],
    'Inhibitory': ['In.VIP', 'In.PVALB', 'In.SST', 'In.LAMP5', 'In.PAX6',
                   'In.DISC1', 'In.Unknown'],
    'Glia': ['Glia.Astro', 'Glia.Oligo', 'Glia.OPC', 'Glia.Micro', 'Astro', 'Oligo', 'OPC', 'Micro'],
    'Vascular': ['Vasc.Endo', 'Vasc.Pericyte', 'Vasc.Fibro', 'Vasc.Capillary',
                 'Endo', 'VLMC', 'Peri', 'Immune']
}

def get_nvu_group(cell_type):
    for group, types in NVU_MAPPING.items():
        for t in types:
            if t in cell_type:
                return group
    return 'Other'

# Upstream vs Downstream modules based on previous analysis
UPSTREAM_MODULES = ['Ion_Transport', 'DNA_Repair', 'Calcium_Signaling', 'RNA_Processing']
DOWNSTREAM_MODULES = ['Mitochondria', 'Protein_Homeostasis', 'ER_Stress', 'Cytoskeleton']

# ============================================================================
# STEP 1: LOAD DATA
# ============================================================================

print("\n[Step 1] Loading data...")

data_file = GSE_RESULTS / 'combined_metadata_with_modules_and_PT.csv'
df = pd.read_csv(data_file)

# Add TDP-43 numeric stage
df['TDP_stage_numeric'] = df['condition'].map(TDP_STAGE_TO_NUMERIC)

# Add NVU group
df['NVU_group'] = df['cell_type'].apply(get_nvu_group)

print(f"  Total cells: {len(df):,}")
print(f"  NVU groups: {df['NVU_group'].value_counts().to_dict()}")

# Get module columns
module_cols = [c for c in df.columns if c.startswith('module_')]
print(f"  Module columns: {len(module_cols)}")

# ============================================================================
# STEP 2: MODULE CHANGE BY TDP STAGE AND NVU GROUP
# ============================================================================

print("\n[Step 2] Analyzing module changes by TDP stage and NVU group...")

# For each module, compute mean score at each TDP stage for each NVU group
results = []

for mod_col in module_cols:
    module_name = mod_col.replace('module_', '')

    for nvu in ['Vascular', 'Glia', 'Excitatory', 'Inhibitory']:
        subset = df[df['NVU_group'] == nvu]
        if len(subset) < 30:
            continue

        # Mean at each stage
        stage_means = subset.groupby('condition')[mod_col].mean()

        # Control vs TDPhigh difference (normalized)
        if 'Control' in stage_means.index and 'TDPhigh' in stage_means.index:
            control_mean = stage_means['Control']
            tdphigh_mean = stage_means['TDPhigh']

            # Effect size (Cohen's d approximation)
            control_std = subset[subset['condition'] == 'Control'][mod_col].std()
            if control_std > 0:
                effect_size = (tdphigh_mean - control_mean) / control_std
            else:
                effect_size = np.nan

            # When does the change first appear? (TDPneg vs Control)
            if 'TDPneg' in stage_means.index:
                tdpneg_mean = stage_means['TDPneg']
                early_change = (tdpneg_mean - control_mean) / control_std if control_std > 0 else np.nan
            else:
                early_change = np.nan

            # Correlation with TDP stage
            rho, pval = stats.spearmanr(subset[mod_col], subset['TDP_stage_numeric'])

            results.append({
                'module': module_name,
                'NVU_group': nvu,
                'n_cells': len(subset),
                'control_mean': control_mean,
                'tdphigh_mean': tdphigh_mean,
                'effect_size': effect_size,
                'early_change': early_change,
                'rho_TDP': rho,
                'p_value': pval
            })

results_df = pd.DataFrame(results)

# ============================================================================
# STEP 3: IDENTIFY WHICH NVU GROUP SHOWS EARLIEST CHANGE
# ============================================================================

print("\n[Step 3] Which NVU group shows earliest change in upstream modules?")

print("\n  === UPSTREAM MODULES (negative TDP correlation = potential cause) ===")
for mod in UPSTREAM_MODULES:
    mod_data = results_df[results_df['module'] == mod].sort_values('early_change')
    print(f"\n  {mod}:")
    if len(mod_data) > 0:
        for _, row in mod_data.iterrows():
            direction = "↑" if row['early_change'] > 0 else "↓"
            print(f"    {row['NVU_group']}: early_change={row['early_change']:.3f}{direction}, rho={row['rho_TDP']:.3f}")

print("\n  === DOWNSTREAM MODULES (positive TDP correlation = consequence) ===")
for mod in DOWNSTREAM_MODULES:
    mod_data = results_df[results_df['module'] == mod].sort_values('early_change', ascending=False)
    print(f"\n  {mod}:")
    if len(mod_data) > 0:
        for _, row in mod_data.iterrows():
            direction = "↑" if row['early_change'] > 0 else "↓"
            print(f"    {row['NVU_group']}: early_change={row['early_change']:.3f}{direction}, rho={row['rho_TDP']:.3f}")

# ============================================================================
# STEP 4: ONSET TIMING ANALYSIS
# ============================================================================

print("\n[Step 4] Computing onset timing (when does module exceed threshold?)...")

# Define onset as: module score > Control mean + 1 SD
onset_results = []

for mod_col in module_cols:
    module_name = mod_col.replace('module_', '')

    # Global Control threshold
    control_data = df[df['condition'] == 'Control'][mod_col]
    threshold = control_data.mean() + control_data.std()

    for nvu in ['Vascular', 'Glia', 'Excitatory', 'Inhibitory']:
        subset = df[df['NVU_group'] == nvu]
        if len(subset) < 30:
            continue

        # Fraction exceeding threshold at each stage
        for stage, stage_num in TDP_STAGE_TO_NUMERIC.items():
            stage_data = subset[subset['condition'] == stage][mod_col]
            if len(stage_data) > 0:
                frac_above = (stage_data > threshold).mean()
                onset_results.append({
                    'module': module_name,
                    'NVU_group': nvu,
                    'stage': stage,
                    'stage_num': stage_num,
                    'frac_above_threshold': frac_above
                })

onset_df = pd.DataFrame(onset_results)

# Find earliest onset (first stage where frac > 10% above Control baseline)
print("\n  === Earliest Onset by Module and NVU Group ===")

earliest_onset = []
for mod in UPSTREAM_MODULES + DOWNSTREAM_MODULES:
    mod_onset = onset_df[onset_df['module'] == mod]

    for nvu in ['Vascular', 'Glia', 'Excitatory', 'Inhibitory']:
        nvu_onset = mod_onset[mod_onset['NVU_group'] == nvu].sort_values('stage_num')

        if len(nvu_onset) == 0:
            continue

        # Control baseline
        control_frac = nvu_onset[nvu_onset['stage'] == 'Control']['frac_above_threshold'].values
        if len(control_frac) == 0:
            continue
        control_frac = control_frac[0]

        # Find first stage with >10% increase
        onset_stage = None
        for _, row in nvu_onset.iterrows():
            if row['frac_above_threshold'] > control_frac + 0.1:
                onset_stage = row['stage']
                break

        earliest_onset.append({
            'module': mod,
            'NVU_group': nvu,
            'onset_stage': onset_stage if onset_stage else 'None',
            'control_frac': control_frac
        })

earliest_df = pd.DataFrame(earliest_onset)

# Pivot for easier viewing
if len(earliest_df) > 0:
    pivot = earliest_df.pivot(index='module', columns='NVU_group', values='onset_stage')
    print(pivot)

# ============================================================================
# STEP 5: GENERATE FIGURES
# ============================================================================

print("\n[Step 5] Generating figures...")

# Figure 1: Upstream module changes by NVU group across TDP stages
fig, axes = plt.subplots(2, 2, figsize=(14, 12))

for i, mod in enumerate(UPSTREAM_MODULES):
    ax = axes[i // 2, i % 2]
    mod_col = f'module_{mod}'

    if mod_col not in df.columns:
        continue

    # Plot mean and error by NVU group and TDP stage
    for nvu, color in [('Vascular', 'red'), ('Glia', 'orange'),
                       ('Excitatory', 'green'), ('Inhibitory', 'blue')]:
        subset = df[df['NVU_group'] == nvu]
        if len(subset) < 30:
            continue

        means = subset.groupby('TDP_stage_numeric')[mod_col].mean()
        sems = subset.groupby('TDP_stage_numeric')[mod_col].sem()

        ax.errorbar(means.index, means.values, yerr=sems.values,
                   label=nvu, color=color, marker='o', capsize=3)

    ax.set_xlabel('TDP-43 Stage (0=Control, 3=TDPhigh)', fontsize=11)
    ax.set_ylabel(f'{mod} Score', fontsize=11)
    ax.set_title(f'{mod} (UPSTREAM - negative TDP corr)', fontsize=12)
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_xticks([0, 1, 2, 3])
    ax.set_xticklabels(['Control', 'TDPneg', 'TDPmed', 'TDPhigh'])

plt.tight_layout()
fig.savefig(OUTPUT_DIR / 'Fig_upstream_modules_by_NVU.png', dpi=150, bbox_inches='tight')
plt.close()
print(f"  Saved: Fig_upstream_modules_by_NVU.png")

# Figure 2: Downstream module changes by NVU group
fig, axes = plt.subplots(2, 2, figsize=(14, 12))

for i, mod in enumerate(DOWNSTREAM_MODULES):
    ax = axes[i // 2, i % 2]
    mod_col = f'module_{mod}'

    if mod_col not in df.columns:
        continue

    for nvu, color in [('Vascular', 'red'), ('Glia', 'orange'),
                       ('Excitatory', 'green'), ('Inhibitory', 'blue')]:
        subset = df[df['NVU_group'] == nvu]
        if len(subset) < 30:
            continue

        means = subset.groupby('TDP_stage_numeric')[mod_col].mean()
        sems = subset.groupby('TDP_stage_numeric')[mod_col].sem()

        ax.errorbar(means.index, means.values, yerr=sems.values,
                   label=nvu, color=color, marker='o', capsize=3)

    ax.set_xlabel('TDP-43 Stage (0=Control, 3=TDPhigh)', fontsize=11)
    ax.set_ylabel(f'{mod} Score', fontsize=11)
    ax.set_title(f'{mod} (DOWNSTREAM - positive TDP corr)', fontsize=12)
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_xticks([0, 1, 2, 3])
    ax.set_xticklabels(['Control', 'TDPneg', 'TDPmed', 'TDPhigh'])

plt.tight_layout()
fig.savefig(OUTPUT_DIR / 'Fig_downstream_modules_by_NVU.png', dpi=150, bbox_inches='tight')
plt.close()
print(f"  Saved: Fig_downstream_modules_by_NVU.png")

# Figure 3: Heatmap of early changes
fig, ax = plt.subplots(figsize=(10, 12))

# Pivot early_change data
pivot_data = results_df.pivot(index='module', columns='NVU_group', values='early_change')

# Reorder modules: upstream first, then downstream
module_order = UPSTREAM_MODULES + DOWNSTREAM_MODULES + \
               [m for m in pivot_data.index if m not in UPSTREAM_MODULES + DOWNSTREAM_MODULES]
pivot_data = pivot_data.reindex([m for m in module_order if m in pivot_data.index])

sns.heatmap(pivot_data, cmap='RdBu_r', center=0, annot=True, fmt='.2f',
            ax=ax, cbar_kws={'label': 'Early Change (TDPneg vs Control)'})
ax.set_title('Module Early Change by NVU Group\n(Negative = decrease early, Positive = increase early)', fontsize=12)
ax.set_ylabel('Module', fontsize=11)
ax.set_xlabel('NVU Group', fontsize=11)

# Add separator line between upstream and downstream
n_upstream = len([m for m in UPSTREAM_MODULES if m in pivot_data.index])
ax.axhline(y=n_upstream, color='black', linewidth=2)

plt.tight_layout()
fig.savefig(OUTPUT_DIR / 'Fig_early_change_heatmap.png', dpi=150, bbox_inches='tight')
plt.close()
print(f"  Saved: Fig_early_change_heatmap.png")

# Figure 4: Causal flow diagram (simplified)
fig, ax = plt.subplots(figsize=(14, 8))

# Get key data for visualization
key_modules = UPSTREAM_MODULES + DOWNSTREAM_MODULES
nvu_groups = ['Vascular', 'Glia', 'Excitatory', 'Inhibitory']

# Create a summary: which NVU shows earliest change for each module
summary_data = []
for mod in key_modules:
    mod_results = results_df[results_df['module'] == mod]

    # Find NVU with most negative early_change (earliest decrease) for upstream
    # or most positive early_change (earliest increase) for downstream
    if mod in UPSTREAM_MODULES:
        earliest = mod_results.loc[mod_results['early_change'].idxmin()] if len(mod_results) > 0 else None
    else:
        earliest = mod_results.loc[mod_results['early_change'].idxmax()] if len(mod_results) > 0 else None

    if earliest is not None:
        summary_data.append({
            'module': mod,
            'earliest_NVU': earliest['NVU_group'],
            'early_change': earliest['early_change'],
            'type': 'upstream' if mod in UPSTREAM_MODULES else 'downstream'
        })

summary_df = pd.DataFrame(summary_data)

# Plot as a flow diagram
y_positions = {'upstream': 1, 'downstream': 0}
nvu_colors = {'Vascular': 'red', 'Glia': 'orange', 'Excitatory': 'green', 'Inhibitory': 'blue'}

for i, row in summary_df.iterrows():
    x = i
    y = y_positions[row['type']]
    color = nvu_colors.get(row['earliest_NVU'], 'gray')

    ax.scatter(x, y, s=500, c=color, alpha=0.7, edgecolors='black', linewidths=2)
    ax.annotate(row['module'], (x, y), textcoords="offset points",
                xytext=(0, 15), ha='center', fontsize=9, rotation=45)
    ax.annotate(f"{row['earliest_NVU']}\n({row['early_change']:.2f})", (x, y),
                textcoords="offset points", xytext=(0, -25), ha='center', fontsize=8)

# Add labels
ax.text(-0.5, 1, 'UPSTREAM\n(TDP-43 negative)', fontsize=12, fontweight='bold', va='center')
ax.text(-0.5, 0, 'DOWNSTREAM\n(TDP-43 positive)', fontsize=12, fontweight='bold', va='center')

# Add legend
for nvu, color in nvu_colors.items():
    ax.scatter([], [], c=color, s=100, label=nvu)
ax.legend(loc='upper right', title='Earliest Change In:')

ax.set_xlim(-1, len(summary_df))
ax.set_ylim(-0.5, 1.5)
ax.axis('off')
ax.set_title('Which NVU Group Shows Earliest Change in Each Module?', fontsize=14)

plt.tight_layout()
fig.savefig(OUTPUT_DIR / 'Fig_causal_flow_summary.png', dpi=150, bbox_inches='tight')
plt.close()
print(f"  Saved: Fig_causal_flow_summary.png")

# Save results
results_df.to_csv(OUTPUT_DIR / 'module_changes_by_NVU.csv', index=False)
summary_df.to_csv(OUTPUT_DIR / 'earliest_change_summary.csv', index=False)

# ============================================================================
# SUMMARY
# ============================================================================

print("\n" + "=" * 80)
print("UPSTREAM MODULE ANALYSIS SUMMARY")
print("=" * 80)

print("""
QUESTION: Do upstream modules (Ion_Transport, DNA_Repair, Calcium_Signaling)
show early changes specifically in Vascular cells?

RESULTS:
""")

print("=== Which NVU shows EARLIEST change in UPSTREAM modules? ===")
for mod in UPSTREAM_MODULES:
    mod_results = results_df[results_df['module'] == mod].sort_values('early_change')
    if len(mod_results) > 0:
        earliest = mod_results.iloc[0]
        print(f"  {mod}: {earliest['NVU_group']} (early_change = {earliest['early_change']:.3f})")

print("\n=== Which NVU shows EARLIEST change in DOWNSTREAM modules? ===")
for mod in DOWNSTREAM_MODULES:
    mod_results = results_df[results_df['module'] == mod].sort_values('early_change', ascending=False)
    if len(mod_results) > 0:
        earliest = mod_results.iloc[0]
        print(f"  {mod}: {earliest['NVU_group']} (early_change = {earliest['early_change']:.3f})")

print(f"\nOutput: {OUTPUT_DIR}")
print("=" * 80)
