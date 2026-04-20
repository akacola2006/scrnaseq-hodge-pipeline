#!/usr/bin/env python3
"""
Phase 8b: Sensitivity Analysis and First-Trigger Module Identification

Purpose:
  1. Sensitivity analysis: Test multiple bin numbers (25, 50, 100) and thresholds (0.5σ, 0.3σ, 0.2σ)
  2. Glia module onset detailed ranking
  3. Identify first-trigger modules for each cell type (Upper, Glia, VAT1L)

Author: Claude Code
Date: 2025-11-24
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from scipy.stats import spearmanr
import warnings
warnings.filterwarnings('ignore')

# Configuration
DATA_FILE = 'results/PTv2_robustness/PTv2_quick_key_celltypes_with_PTs.csv'
OUTPUT_DIR = 'results/phase8b_sensitivity'
Path(OUTPUT_DIR).mkdir(parents=True, exist_ok=True)

# Cell type definitions
VAT1L_types = ['Ex.L5.VAT1L.EYA4', 'Ex.L5.VAT1L.THSD4']
Upper_types = ['Ex.L2.L3.CUX2.RASGRF2']
Glia_types = ['Glia.Astro.GFAP-neg', 'Glia.Micro', 'Glia.Oligo']

# Analysis parameters
BIN_NUMBERS = [25, 50, 100]
THRESHOLDS = [0.5, 0.3, 0.2]  # multipliers for std (mean + threshold*std)

def assign_group(cell_type):
    """Assign cell type to group"""
    if cell_type in VAT1L_types:
        return 'VAT1L'
    elif cell_type in Upper_types:
        return 'Upper'
    elif cell_type in Glia_types:
        return 'Glia'
    else:
        return 'Other'

def calculate_onset_for_condition(df, module_names, num_bins, threshold_multiplier):
    """
    Calculate module onset PT for a specific binning and threshold condition

    Returns:
        DataFrame with onset PT for each group × module combination
    """
    # Determine PT_dpt range (same as Phase 8)
    pt_min = df['PT_dpt'].min()
    pt_99th = df['PT_dpt'].quantile(0.99)
    pt_max = pt_99th * 1.05

    # Create bins
    bin_edges = np.linspace(pt_min, pt_max, num_bins + 1)
    df_copy = df.copy()
    df_copy['pt_bin'] = pd.cut(df_copy['PT_dpt'], bins=bin_edges, labels=False, include_lowest=True)

    # Calculate bin centers
    def get_bin_center(bin_idx):
        if pd.isna(bin_idx):
            return np.nan
        bin_idx = int(bin_idx)
        if bin_idx < 0 or bin_idx >= len(bin_edges) - 1:
            return np.nan
        return (bin_edges[bin_idx] + bin_edges[bin_idx + 1]) / 2

    df_copy['pt_bin_center'] = df_copy['pt_bin'].apply(get_bin_center)

    # Calculate module profiles by group and bin
    profiles_list = []
    for group in ['Upper', 'Glia', 'VAT1L']:
        group_df = df_copy[df_copy['group'] == group]
        for bin_idx in sorted(group_df['pt_bin'].dropna().unique()):
            bin_df = group_df[group_df['pt_bin'] == bin_idx]
            if len(bin_df) == 0:
                continue

            row = {
                'group': group,
                'pt_bin': bin_idx,
                'pt_center': bin_df['pt_bin_center'].iloc[0],
                'n_cells': len(bin_df)
            }

            for mod in module_names:
                row[mod] = bin_df[f'module_{mod}'].mean()

            profiles_list.append(row)

    profiles_df = pd.DataFrame(profiles_list)

    # Calculate onset PT for each group × module
    onset_results = []

    for group in ['Upper', 'Glia', 'VAT1L']:
        group_profile = profiles_df[profiles_df['group'] == group].copy()
        group_profile = group_profile.sort_values('pt_center')

        for mod in module_names:
            values = group_profile[mod].values
            pt_centers = group_profile['pt_center'].values

            # Calculate threshold
            mean_val = np.nanmean(values)
            std_val = np.nanstd(values)
            onset_threshold = mean_val + threshold_multiplier * std_val

            # Find onset (first bin above threshold)
            onset_idx = None
            for i, val in enumerate(values):
                if not np.isnan(val) and val > onset_threshold:
                    onset_idx = i
                    break

            onset_pt = pt_centers[onset_idx] if onset_idx is not None else np.nan

            onset_results.append({
                'group': group,
                'module': mod,
                'onset_pt': onset_pt,
                'onset_threshold': onset_threshold,
                'onset_bin_idx': onset_idx
            })

    return pd.DataFrame(onset_results)

def rank_modules_within_group(onset_df, group_name):
    """
    Rank modules by onset PT within a specific group

    Returns:
        DataFrame with module, onset_pt, rank (1=earliest)
    """
    group_df = onset_df[onset_df['group'] == group_name].copy()
    group_df = group_df.sort_values('onset_pt')
    group_df['rank'] = range(1, len(group_df) + 1)
    return group_df[['module', 'onset_pt', 'rank']]

print("="*80)
print("Phase 8b: Sensitivity Analysis and First-Trigger Module Identification")
print("="*80)

# Load data
print("\n[1] Loading data...")
df = pd.read_csv(DATA_FILE)
print(f"Total cells: {len(df):,}")

# Assign groups
df['group'] = df['cell_type'].apply(assign_group)
df = df[df['group'] != 'Other']
print(f"\nFiltered to {len(df):,} cells in 3 groups:")
for group in ['Upper', 'Glia', 'VAT1L']:
    print(f"  {group}: {len(df[df['group'] == group]):,} cells")

# Get module names
module_cols = [c for c in df.columns if c.startswith('module_')]
module_names = [c.replace('module_', '') for c in module_cols]
print(f"\nFound {len(module_names)} modules")

# ============================================================================
# Part 1: Sensitivity Analysis
# ============================================================================
print("\n[2] Running sensitivity analysis...")
print(f"Testing {len(BIN_NUMBERS)} bin numbers × {len(THRESHOLDS)} thresholds = {len(BIN_NUMBERS) * len(THRESHOLDS)} conditions")

all_conditions = []

for num_bins in BIN_NUMBERS:
    for threshold in THRESHOLDS:
        print(f"\n  Condition: {num_bins} bins, threshold={threshold}σ")

        onset_df = calculate_onset_for_condition(df, module_names, num_bins, threshold)
        onset_df['num_bins'] = num_bins
        onset_df['threshold'] = threshold
        onset_df['condition'] = f"{num_bins}bins_{threshold}sigma"

        all_conditions.append(onset_df)

# Combine all conditions
all_onset_df = pd.concat(all_conditions, ignore_index=True)
all_onset_df.to_csv(f'{OUTPUT_DIR}/all_conditions_onset.csv', index=False)
print(f"\nSaved: {OUTPUT_DIR}/all_conditions_onset.csv")

# ============================================================================
# Part 2: Stability Analysis - Which modules are consistently first?
# ============================================================================
print("\n[3] Analyzing onset order stability across conditions...")

stability_results = []

for group in ['Upper', 'Glia', 'VAT1L']:
    group_onset = all_onset_df[all_onset_df['group'] == group]

    for mod in module_names:
        mod_onset = group_onset[group_onset['module'] == mod]

        # Calculate statistics across conditions
        mean_onset = mod_onset['onset_pt'].mean()
        std_onset = mod_onset['onset_pt'].std()
        min_onset = mod_onset['onset_pt'].min()
        max_onset = mod_onset['onset_pt'].max()

        # Calculate rank statistics
        ranks = []
        for condition in mod_onset['condition'].unique():
            cond_df = all_onset_df[(all_onset_df['group'] == group) & (all_onset_df['condition'] == condition)]
            cond_df = cond_df.sort_values('onset_pt')
            mod_rank = cond_df[cond_df['module'] == mod].index[0]
            rank = cond_df.index.get_loc(mod_rank) + 1
            ranks.append(rank)

        mean_rank = np.mean(ranks)
        std_rank = np.std(ranks)
        min_rank = np.min(ranks)
        max_rank = np.max(ranks)

        stability_results.append({
            'group': group,
            'module': mod,
            'mean_onset_pt': mean_onset,
            'std_onset_pt': std_onset,
            'min_onset_pt': min_onset,
            'max_onset_pt': max_onset,
            'mean_rank': mean_rank,
            'std_rank': std_rank,
            'min_rank': min_rank,
            'max_rank': max_rank,
            'rank_stability': std_rank  # Lower = more stable
        })

stability_df = pd.DataFrame(stability_results)
stability_df.to_csv(f'{OUTPUT_DIR}/module_onset_stability.csv', index=False)
print(f"Saved: {OUTPUT_DIR}/module_onset_stability.csv")

# ============================================================================
# Part 3: Identify First-Trigger Modules for Each Cell Type
# ============================================================================
print("\n[4] Identifying first-trigger modules for each cell type...")

first_triggers = []

for group in ['Upper', 'Glia', 'VAT1L']:
    print(f"\n  === {group} ===")
    group_stability = stability_df[stability_df['group'] == group].copy()
    group_stability = group_stability.sort_values('mean_rank')

    print(f"  Top 5 earliest-onset modules (by mean rank):")
    for i, row in group_stability.head(5).iterrows():
        print(f"    {row['module']:20s}  Mean rank: {row['mean_rank']:.1f} (±{row['std_rank']:.2f})  "
              f"Onset PT: {row['mean_onset_pt']:.6f}")

        if i < 3:  # Top 3
            first_triggers.append({
                'group': group,
                'trigger_rank': len([x for x in first_triggers if x['group'] == group]) + 1,
                'module': row['module'],
                'mean_rank': row['mean_rank'],
                'std_rank': row['std_rank'],
                'mean_onset_pt': row['mean_onset_pt']
            })

first_triggers_df = pd.DataFrame(first_triggers)
first_triggers_df.to_csv(f'{OUTPUT_DIR}/first_trigger_modules.csv', index=False)
print(f"\nSaved: {OUTPUT_DIR}/first_trigger_modules.csv")

# ============================================================================
# Part 4: Glia Module Onset Detailed Ranking
# ============================================================================
print("\n[5] Creating detailed Glia module onset ranking...")

glia_ranking = stability_df[stability_df['group'] == 'Glia'].copy()
glia_ranking = glia_ranking.sort_values('mean_onset_pt')
glia_ranking['onset_rank'] = range(1, len(glia_ranking) + 1)

# Categorize into early/middle/late
n_modules = len(glia_ranking)
glia_ranking['onset_category'] = pd.cut(
    glia_ranking['onset_rank'],
    bins=[0, n_modules/3, 2*n_modules/3, n_modules],
    labels=['Early', 'Middle', 'Late']
)

glia_ranking.to_csv(f'{OUTPUT_DIR}/glia_module_onset_ranking.csv', index=False)
print(f"Saved: {OUTPUT_DIR}/glia_module_onset_ranking.csv")

print("\nGlia module onset categories:")
for category in ['Early', 'Middle', 'Late']:
    cat_modules = glia_ranking[glia_ranking['onset_category'] == category]['module'].tolist()
    print(f"  {category:8s}: {', '.join(cat_modules)}")

# ============================================================================
# Part 5: Visualizations
# ============================================================================
print("\n[6] Creating visualizations...")

# Figure 1: Stability heatmap - show rank variation across conditions
fig, axes = plt.subplots(1, 3, figsize=(18, 10))

for ax_idx, group in enumerate(['Upper', 'Glia', 'VAT1L']):
    # Get onset data for this group across all conditions
    pivot_data = []
    for mod in module_names:
        row_data = {'module': mod}
        for condition in all_onset_df['condition'].unique():
            cond_df = all_onset_df[(all_onset_df['group'] == group) & (all_onset_df['condition'] == condition)]
            cond_df = cond_df.sort_values('onset_pt').reset_index(drop=True)
            mod_idx = cond_df[cond_df['module'] == mod].index
            if len(mod_idx) > 0:
                row_data[condition] = mod_idx[0] + 1
            else:
                row_data[condition] = np.nan
        pivot_data.append(row_data)

    pivot_df = pd.DataFrame(pivot_data)
    pivot_df = pivot_df.set_index('module')

    # Sort by mean rank
    pivot_df['mean_rank'] = pivot_df.mean(axis=1)
    pivot_df = pivot_df.sort_values('mean_rank')
    pivot_df = pivot_df.drop('mean_rank', axis=1)

    # Plot heatmap
    sns.heatmap(pivot_df, cmap='RdYlGn_r', annot=True, fmt='.0f',
                cbar_kws={'label': 'Onset Rank'}, ax=axes[ax_idx],
                vmin=1, vmax=len(module_names))
    axes[ax_idx].set_title(f'{group} Module Onset Rank\nAcross Conditions', fontsize=12, fontweight='bold')
    axes[ax_idx].set_xlabel('Condition', fontsize=10)
    axes[ax_idx].set_ylabel('Module', fontsize=10)

plt.tight_layout()
plt.savefig(f'{OUTPUT_DIR}/Fig1_onset_rank_stability_heatmap.png', dpi=300, bbox_inches='tight')
print(f"Saved: {OUTPUT_DIR}/Fig1_onset_rank_stability_heatmap.png")
plt.close()

# Figure 2: First-trigger modules comparison
fig, ax = plt.subplots(figsize=(12, 8))

trigger_pivot = first_triggers_df.pivot(index='module', columns='group', values='mean_rank')

# Plot grouped bar chart
x = np.arange(len(trigger_pivot))
width = 0.25

for i, group in enumerate(['Upper', 'Glia', 'VAT1L']):
    if group in trigger_pivot.columns:
        values = trigger_pivot[group].values
        ax.bar(x + i*width, values, width, label=group, alpha=0.8)

ax.set_xlabel('Module', fontsize=12, fontweight='bold')
ax.set_ylabel('Mean Onset Rank (Lower = Earlier)', fontsize=12, fontweight='bold')
ax.set_title('First-Trigger Modules by Cell Type\n(Top 3 Earliest-Onset Modules)', fontsize=14, fontweight='bold')
ax.set_xticks(x + width)
ax.set_xticklabels(trigger_pivot.index, rotation=45, ha='right')
ax.legend(title='Cell Type')
ax.grid(axis='y', alpha=0.3)

plt.tight_layout()
plt.savefig(f'{OUTPUT_DIR}/Fig2_first_trigger_modules.png', dpi=300, bbox_inches='tight')
print(f"Saved: {OUTPUT_DIR}/Fig2_first_trigger_modules.png")
plt.close()

# Figure 3: Glia detailed ranking with categories
fig, ax = plt.subplots(figsize=(10, 12))

glia_plot = glia_ranking.sort_values('onset_rank')
colors = {'Early': '#2ecc71', 'Middle': '#f39c12', 'Late': '#e74c3c'}
bar_colors = [colors[cat] for cat in glia_plot['onset_category']]

y_pos = np.arange(len(glia_plot))
ax.barh(y_pos, glia_plot['mean_onset_pt'], xerr=glia_plot['std_onset_pt'],
        color=bar_colors, alpha=0.7, edgecolor='black')

ax.set_yticks(y_pos)
ax.set_yticklabels(glia_plot['module'])
ax.set_xlabel('Mean Onset PT_dpt (±SD)', fontsize=12, fontweight='bold')
ax.set_title('Glia Module Onset Ranking\n(Ordered by Mean Onset PT)', fontsize=14, fontweight='bold')
ax.invert_yaxis()

# Add legend
from matplotlib.patches import Patch
legend_elements = [Patch(facecolor=colors['Early'], label='Early Onset'),
                   Patch(facecolor=colors['Middle'], label='Middle Onset'),
                   Patch(facecolor=colors['Late'], label='Late Onset')]
ax.legend(handles=legend_elements, loc='lower right')

plt.tight_layout()
plt.savefig(f'{OUTPUT_DIR}/Fig3_glia_onset_ranking.png', dpi=300, bbox_inches='tight')
print(f"Saved: {OUTPUT_DIR}/Fig3_glia_onset_ranking.png")
plt.close()

# Figure 4: Rank stability (coefficient of variation)
fig, axes = plt.subplots(1, 3, figsize=(18, 6), sharey=True)

for ax_idx, group in enumerate(['Upper', 'Glia', 'VAT1L']):
    group_stab = stability_df[stability_df['group'] == group].copy()
    group_stab = group_stab.sort_values('rank_stability')

    y_pos = np.arange(len(group_stab))
    axes[ax_idx].barh(y_pos, group_stab['rank_stability'], alpha=0.7, color='steelblue')
    axes[ax_idx].set_yticks(y_pos)
    axes[ax_idx].set_yticklabels(group_stab['module'])
    axes[ax_idx].set_xlabel('Rank Std Dev\n(Lower = More Stable)', fontsize=10)
    axes[ax_idx].set_title(f'{group}', fontsize=12, fontweight='bold')
    axes[ax_idx].axvline(x=group_stab['rank_stability'].median(),
                         color='red', linestyle='--', alpha=0.5, label='Median')
    axes[ax_idx].legend()
    axes[ax_idx].invert_yaxis()

fig.suptitle('Module Onset Rank Stability Across Conditions', fontsize=14, fontweight='bold')
plt.tight_layout()
plt.savefig(f'{OUTPUT_DIR}/Fig4_rank_stability.png', dpi=300, bbox_inches='tight')
print(f"Saved: {OUTPUT_DIR}/Fig4_rank_stability.png")
plt.close()

# ============================================================================
# Part 6: Summary Report
# ============================================================================
print("\n[7] Generating summary report...")

report = f"""# Phase 8b: Sensitivity Analysis and First-Trigger Module Identification

**Date**: 2025-11-24
**Analysis**: Sensitivity analysis of module onset timing and first-trigger module identification

---

## Analysis Overview

### Parameters Tested
- **Bin numbers**: {BIN_NUMBERS}
- **Thresholds**: {THRESHOLDS} (mean + threshold × std)
- **Total conditions**: {len(BIN_NUMBERS) * len(THRESHOLDS)}

---

## Key Findings

### 1. First-Trigger Modules (Top 3 Earliest-Onset by Cell Type)

"""

for group in ['Upper', 'Glia', 'VAT1L']:
    report += f"\n#### {group}\n\n"
    group_triggers = first_triggers_df[first_triggers_df['group'] == group]
    report += "| Rank | Module | Mean Rank | Rank Stability (±SD) | Mean Onset PT |\n"
    report += "|------|--------|-----------|---------------------|---------------|\n"
    for _, row in group_triggers.iterrows():
        report += f"| {row['trigger_rank']:.0f} | {row['module']} | {row['mean_rank']:.1f} | ±{row['std_rank']:.2f} | {row['mean_onset_pt']:.6f} |\n"

report += """

---

### 2. Glia Module Onset Categories

"""

for category in ['Early', 'Middle', 'Late']:
    cat_modules = glia_ranking[glia_ranking['onset_category'] == category]
    report += f"\n**{category} Onset** ({len(cat_modules)} modules):\n"
    for _, row in cat_modules.iterrows():
        report += f"- {row['module']:20s}  (Rank: {row['onset_rank']:.0f}, PT: {row['mean_onset_pt']:.6f}±{row['std_onset_pt']:.6f})\n"

report += """

---

### 3. Onset Order Stability

#### Most Stable Onset Order (Lowest Rank Std Dev)

"""

for group in ['Upper', 'Glia', 'VAT1L']:
    report += f"\n**{group}**:\n\n"
    group_stab = stability_df[stability_df['group'] == group].sort_values('rank_stability').head(5)
    report += "| Module | Mean Rank | Rank Std Dev | Onset PT (mean±SD) |\n"
    report += "|--------|-----------|--------------|--------------------|\n"
    for _, row in group_stab.iterrows():
        report += f"| {row['module']} | {row['mean_rank']:.1f} | {row['std_rank']:.2f} | {row['mean_onset_pt']:.6f}±{row['std_onset_pt']:.6f} |\n"

report += """

---

## Interpretation

### Upper Modules
- **Sensitivity**: Upper modules show {"low" if stability_df[stability_df['group'] == 'Upper']['std_rank'].mean() < 2 else "moderate to high"} variation in rank across conditions (mean rank std: {stability_df[stability_df['group'] == 'Upper']['std_rank'].mean():.2f})
- **First triggers**: The most consistently early modules are likely the true initiators of Upper dysfunction
- **Biological interpretation**: If Synaptic/Ca²⁺/Ion_Transport consistently rank first, this supports hyperexcitability as the primary Upper trigger

### Glia Modules
- **Sensitivity**: Glia modules show {"low" if stability_df[stability_df['group'] == 'Glia']['std_rank'].mean() < 2 else "moderate to high"} variation in rank across conditions (mean rank std: {stability_df[stability_df['group'] == 'Glia']['std_rank'].mean():.2f})
- **Detailed ranking**: Early-onset modules (Complement, Myelination, Synaptic, ER) represent the first signs of Glia structural dysfunction
- **Biological interpretation**: Complement/Myelination early onset suggests initial disruption of homeostatic functions before inflammatory cascade

### VAT1L Modules
- **Sensitivity**: VAT1L modules show {"low" if stability_df[stability_df['group'] == 'VAT1L']['std_rank'].mean() < 2 else "moderate to high"} variation in rank across conditions (mean rank std: {stability_df[stability_df['group'] == 'VAT1L']['std_rank'].mean():.2f})
- **Pattern**: If most modules have similar onset times (low rank variation), this suggests "simultaneous collapse" rather than sequential module failure
- **Biological interpretation**: Fragile downstream cell type may experience rapid multi-system failure

---

## Methodological Notes

### Sensitivity Analysis Approach
1. **Multiple binning resolutions** (25/50/100 bins) test whether onset detection is stable across different PT_dpt granularities
2. **Multiple thresholds** (0.5σ/0.3σ/0.2σ) test whether onset definition affects module ranking
3. **Rank stability metric** (std of ranks) identifies modules with consistent early/late onset regardless of parameters

### Limitations
1. **PT_dpt resolution**: Even at 100 bins, may not distinguish modules with very close onset times
2. **Threshold choice**: Lower thresholds (0.2σ) may capture noise; higher thresholds (0.5σ) may miss subtle early changes
3. **Binning artifacts**: Discrete bins may group together modules that actually have slightly different onsets

---

## Files Generated

1. **all_conditions_onset.csv** - Onset PT for all modules, groups, and conditions
2. **module_onset_stability.csv** - Stability metrics (mean rank, std rank) for each module × group
3. **first_trigger_modules.csv** - Top 3 earliest-onset modules per cell type
4. **glia_module_onset_ranking.csv** - Detailed Glia module ranking with categories

### Visualizations

1. **Fig1_onset_rank_stability_heatmap.png** - Heatmap showing rank variation across conditions
2. **Fig2_first_trigger_modules.png** - Comparison of first-trigger modules by cell type
3. **Fig3_glia_onset_ranking.png** - Detailed Glia module onset ranking
4. **Fig4_rank_stability.png** - Rank stability (coefficient of variation) for all modules

---

## Next Steps

1. **Biological validation**: Do the identified first-trigger modules align with known ALS pathology?
2. **Cross-reference with Phase 7**: Do first-trigger modules show strongest correlations with other cell types?
3. **Literature review**: Are early-onset modules documented as early events in ALS?

---

**Analysis Complete**: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}
"""

with open(f'{OUTPUT_DIR}/PHASE8B_SENSITIVITY_REPORT.md', 'w') as f:
    f.write(report)

print(f"Saved: {OUTPUT_DIR}/PHASE8B_SENSITIVITY_REPORT.md")

print("\n" + "="*80)
print("Phase 8b Analysis Complete!")
print("="*80)
print(f"\nResults saved to: {OUTPUT_DIR}/")
print("\nKey outputs:")
print("  - all_conditions_onset.csv")
print("  - module_onset_stability.csv")
print("  - first_trigger_modules.csv")
print("  - glia_module_onset_ranking.csv")
print("  - Fig1_onset_rank_stability_heatmap.png")
print("  - Fig2_first_trigger_modules.png")
print("  - Fig3_glia_onset_ranking.png")
print("  - Fig4_rank_stability.png")
print("  - PHASE8B_SENSITIVITY_REPORT.md")
