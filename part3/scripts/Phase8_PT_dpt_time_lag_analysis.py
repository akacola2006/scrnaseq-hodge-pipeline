#!/usr/bin/env python3
"""
Phase 8: PT_dpt-based Time-Lag Analysis
========================================

Goal: Determine which cell type (Upper, Glia, VAT1L) shows module changes FIRST
      along PT_dpt pseudotime.

Key Question: Do Upper cell modules precede Glia/VAT1L changes,
              or does Glia lead, or VAT1L?

Method:
  1. Bin cells by PT_dpt
  2. Calculate module averages for each group × bin
  3. Identify "onset" PT_dpt (when module first rises above threshold)
  4. Calculate time-lags: ΔPT(Upper vs Glia), ΔPT(Upper vs VAT1L)
  5. Quantify how many modules show Upper-first vs Glia-first patterns

Output:
  - Module profiles by group and bin (CSV)
  - Onset/peak PT_dpt for each group × module (CSV)
  - Time-lag statistics (ΔPT) (CSV)
  - Visualization plots (PNG)
  - Summary report (MD)
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

# Create output directory
import os
output_dir = 'results/phase8_time_lag'
os.makedirs(output_dir, exist_ok=True)

print("="*80)
print("Phase 8: PT_dpt Time-Lag Analysis")
print("="*80)

# ==============================================================================
# STEP 1: Load Data and Define Groups
# ==============================================================================

print("\n[1] Loading data...")
df = pd.read_csv('results/PTv2_robustness/PTv2_quick_key_celltypes_with_PTs.csv')

print(f"Total cells: {len(df):,}")
print(f"Unique cell types: {df['cell_type'].nunique()}")

# Check cell types
print("\nCell types in data:")
for ct in sorted(df['cell_type'].unique()):
    n = (df['cell_type'] == ct).sum()
    print(f"  {ct}: {n:,} cells")

# Define cell type groups
VAT1L_types = ['Ex.L5.VAT1L.EYA4', 'Ex.L5.VAT1L.THSD4']
Upper_types = ['Ex.L2.L3.CUX2.RASGRF2']  # Only this one in PTv2 data
Glia_types = ['Glia.Astro.GFAP-neg', 'Glia.Micro', 'Glia.Oligo']

# Verify all types exist
for group_name, types in [('VAT1L', VAT1L_types), ('Upper', Upper_types), ('Glia', Glia_types)]:
    existing = [t for t in types if t in df['cell_type'].unique()]
    missing = [t for t in types if t not in df['cell_type'].unique()]
    print(f"\n{group_name} group:")
    print(f"  Existing types: {existing}")
    if missing:
        print(f"  Missing types: {missing}")

# Add group column
def assign_group(cell_type):
    if cell_type in VAT1L_types:
        return 'VAT1L'
    elif cell_type in Upper_types:
        return 'Upper'
    elif cell_type in Glia_types:
        return 'Glia'
    else:
        return 'Other'

df['group'] = df['cell_type'].apply(assign_group)

# Filter to only groups of interest
df = df[df['group'].isin(['Upper', 'Glia', 'VAT1L'])].copy()

print(f"\nFiltered to {len(df):,} cells in 3 groups:")
for group in ['Upper', 'Glia', 'VAT1L']:
    n = (df['group'] == group).sum()
    print(f"  {group}: {n:,} cells")

# ==============================================================================
# STEP 2: PT_dpt Binning and Module Averaging
# ==============================================================================

print("\n[2] Binning by PT_dpt...")

# Check PT_dpt range
pt_min = df['PT_dpt'].min()
pt_max_full = df['PT_dpt'].max()
pt_99th = df['PT_dpt'].quantile(0.99)

print(f"PT_dpt full range: {pt_min:.6f} - {pt_max_full:.6f}")
print(f"PT_dpt 99th percentile: {pt_99th:.6f}")

# Use 99th percentile as max to avoid extreme outliers
# Most cells (99%) are in 0-0.024 range
pt_max = pt_99th * 1.05  # Slightly above 99th percentile

print(f"Using range for binning: {pt_min:.6f} - {pt_max:.6f}")
print(f"Cells in this range: {(df['PT_dpt'] <= pt_max).sum():,} ({100*(df['PT_dpt'] <= pt_max).sum()/len(df):.1f}%)")

# Create bins (25 bins for finer resolution)
num_bins = 25
bin_edges = np.linspace(pt_min, pt_max, num_bins + 1)
df['pt_bin'] = pd.cut(df['PT_dpt'], bins=bin_edges, labels=False, include_lowest=True)

# Calculate bin centers
def get_bin_center(bin_idx):
    if pd.isna(bin_idx):
        return np.nan
    bin_idx = int(bin_idx)
    if bin_idx < 0 or bin_idx >= len(bin_edges) - 1:
        return np.nan
    return (bin_edges[bin_idx] + bin_edges[bin_idx + 1]) / 2

df['pt_bin_center'] = df['pt_bin'].apply(get_bin_center)

print(f"Created {num_bins} bins")

# Get module columns
module_cols = [col for col in df.columns if col.startswith('module_')]
print(f"\nFound {len(module_cols)} modules:")
for col in module_cols[:5]:
    print(f"  {col}")
print(f"  ... and {len(module_cols)-5} more")

# Rename modules (remove 'module_' prefix for easier handling)
module_names = []
for col in module_cols:
    mod_name = col.replace('module_', '')
    df[mod_name] = df[col]
    module_names.append(mod_name)

print(f"\nProcessing {len(module_names)} modules:")
print(module_names)

# Calculate group × bin × module averages
print("\n[2.1] Calculating module profiles by group and bin...")

profiles = []
for group in ['Upper', 'Glia', 'VAT1L']:
    group_cells = df[df['group'] == group]

    for bin_idx in range(num_bins):
        bin_cells = group_cells[group_cells['pt_bin'] == bin_idx]

        if len(bin_cells) == 0:
            continue

        pt_center = bin_edges[bin_idx] + (bin_edges[bin_idx+1] - bin_edges[bin_idx]) / 2

        profile = {
            'group': group,
            'bin_index': bin_idx,
            'pt_center': pt_center,
            'n_cells': len(bin_cells)
        }

        # Module means
        for mod in module_names:
            if mod in bin_cells.columns:
                profile[mod] = bin_cells[mod].mean()

        profiles.append(profile)

profiles_df = pd.DataFrame(profiles)
print(f"Created {len(profiles_df)} group × bin combinations")

# Save
profiles_df.to_csv(f'{output_dir}/module_profiles_by_group_and_bin.csv', index=False)
print(f"Saved: {output_dir}/module_profiles_by_group_and_bin.csv")

# ==============================================================================
# STEP 3: Identify Onset and Peak PT_dpt for Each Group × Module
# ==============================================================================

print("\n[3] Identifying onset and peak PT_dpt...")

# For each group × module, find:
# - Onset: First bin where module > (mean + 0.5*std) across all bins
# - Peak: Bin with maximum module value

onset_peak_results = []

for group in ['Upper', 'Glia', 'VAT1L']:
    group_profile = profiles_df[profiles_df['group'] == group].copy()

    # Sort by PT
    group_profile = group_profile.sort_values('pt_center')

    for mod in module_names:
        if mod not in group_profile.columns:
            continue

        # Get module values
        values = group_profile[mod].values
        pt_centers = group_profile['pt_center'].values

        if len(values) == 0 or np.all(np.isnan(values)):
            continue

        # Calculate mean and std
        mean_val = np.nanmean(values)
        std_val = np.nanstd(values)

        # Onset threshold: mean + 0.5*std
        onset_threshold = mean_val + 0.5 * std_val

        # Find onset (first bin above threshold)
        onset_idx = None
        for i, val in enumerate(values):
            if not np.isnan(val) and val > onset_threshold:
                onset_idx = i
                break

        onset_pt = pt_centers[onset_idx] if onset_idx is not None else np.nan

        # Find peak (max value)
        peak_idx = np.nanargmax(values)
        peak_pt = pt_centers[peak_idx]
        peak_val = values[peak_idx]

        onset_peak_results.append({
            'group': group,
            'module': mod,
            'onset_PT': onset_pt,
            'peak_PT': peak_pt,
            'peak_value': peak_val,
            'mean_value': mean_val,
            'std_value': std_val,
            'onset_threshold': onset_threshold
        })

onset_peak_df = pd.DataFrame(onset_peak_results)
print(f"Calculated onset/peak for {len(onset_peak_df)} group × module combinations")

# Save
onset_peak_df.to_csv(f'{output_dir}/module_onset_peak_PT_by_group.csv', index=False)
print(f"Saved: {output_dir}/module_onset_peak_PT_by_group.csv")

# ==============================================================================
# STEP 4: Calculate Time-Lags (ΔPT)
# ==============================================================================

print("\n[4] Calculating time-lags (ΔPT)...")

# Pivot to get onset/peak PT for each group
onset_pivot = onset_peak_df.pivot(index='module', columns='group', values='onset_PT')
peak_pivot = onset_peak_df.pivot(index='module', columns='group', values='peak_PT')

# Calculate ΔPT
time_lags = pd.DataFrame({
    'module': onset_pivot.index,
    'Upper_onset_PT': onset_pivot['Upper'],
    'Glia_onset_PT': onset_pivot['Glia'],
    'VAT1L_onset_PT': onset_pivot['VAT1L'],
    'Upper_peak_PT': peak_pivot['Upper'],
    'Glia_peak_PT': peak_pivot['Glia'],
    'VAT1L_peak_PT': peak_pivot['VAT1L']
})

# ΔPT (positive = Upper is later, negative = Upper is earlier)
time_lags['Delta_onset_Upper_minus_Glia'] = time_lags['Upper_onset_PT'] - time_lags['Glia_onset_PT']
time_lags['Delta_onset_Upper_minus_VAT1L'] = time_lags['Upper_onset_PT'] - time_lags['VAT1L_onset_PT']
time_lags['Delta_peak_Upper_minus_Glia'] = time_lags['Upper_peak_PT'] - time_lags['Glia_peak_PT']
time_lags['Delta_peak_Upper_minus_VAT1L'] = time_lags['Upper_peak_PT'] - time_lags['VAT1L_peak_PT']

# Save
time_lags.to_csv(f'{output_dir}/module_time_lags.csv', index=False)
print(f"Saved: {output_dir}/module_time_lags.csv")

# Print summary
print("\nTime-lag summary (onset):")
print(f"  Upper vs Glia:")
print(f"    Mean ΔPT: {time_lags['Delta_onset_Upper_minus_Glia'].mean():.6f}")
print(f"    Median ΔPT: {time_lags['Delta_onset_Upper_minus_Glia'].median():.6f}")
print(f"    Upper-first (Δ<0): {(time_lags['Delta_onset_Upper_minus_Glia'] < 0).sum()} modules")
print(f"    Glia-first (Δ>0): {(time_lags['Delta_onset_Upper_minus_Glia'] > 0).sum()} modules")

print(f"\n  Upper vs VAT1L:")
print(f"    Mean ΔPT: {time_lags['Delta_onset_Upper_minus_VAT1L'].mean():.6f}")
print(f"    Median ΔPT: {time_lags['Delta_onset_Upper_minus_VAT1L'].median():.6f}")
print(f"    Upper-first (Δ<0): {(time_lags['Delta_onset_Upper_minus_VAT1L'] < 0).sum()} modules")
print(f"    VAT1L-first (Δ>0): {(time_lags['Delta_onset_Upper_minus_VAT1L'] > 0).sum()} modules")

# ==============================================================================
# STEP 5: Visualizations
# ==============================================================================

print("\n[5] Creating visualizations...")

# Figure 1: Time-lag bar plot (onset)
fig, axes = plt.subplots(2, 1, figsize=(14, 10))

# Upper vs Glia
ax = axes[0]
data = time_lags.sort_values('Delta_onset_Upper_minus_Glia')
colors = ['red' if x < 0 else 'blue' for x in data['Delta_onset_Upper_minus_Glia']]
ax.barh(range(len(data)), data['Delta_onset_Upper_minus_Glia'], color=colors, alpha=0.7)
ax.set_yticks(range(len(data)))
ax.set_yticklabels(data['module'], fontsize=8)
ax.axvline(0, color='black', linestyle='--', linewidth=1)
ax.set_xlabel('ΔPT (Upper onset - Glia onset)', fontsize=10)
ax.set_title('Module Onset Time-Lags: Upper vs Glia\n(Negative = Upper first, Positive = Glia first)', fontsize=12, fontweight='bold')
ax.grid(axis='x', alpha=0.3)

# Upper vs VAT1L
ax = axes[1]
data = time_lags.sort_values('Delta_onset_Upper_minus_VAT1L')
colors = ['red' if x < 0 else 'green' for x in data['Delta_onset_Upper_minus_VAT1L']]
ax.barh(range(len(data)), data['Delta_onset_Upper_minus_VAT1L'], color=colors, alpha=0.7)
ax.set_yticks(range(len(data)))
ax.set_yticklabels(data['module'], fontsize=8)
ax.axvline(0, color='black', linestyle='--', linewidth=1)
ax.set_xlabel('ΔPT (Upper onset - VAT1L onset)', fontsize=10)
ax.set_title('Module Onset Time-Lags: Upper vs VAT1L\n(Negative = Upper first, Positive = VAT1L first)', fontsize=12, fontweight='bold')
ax.grid(axis='x', alpha=0.3)

plt.tight_layout()
plt.savefig(f'{output_dir}/Fig_module_onset_time_lags.png', dpi=300, bbox_inches='tight')
print(f"Saved: {output_dir}/Fig_module_onset_time_lags.png")
plt.close()

# Figure 2: Module trajectories (selected modules)
selected_modules = ['Synaptic', 'Calcium_Signaling', 'ER_Stress', 'Inflammation',
                    'Oxidative_Stress', 'Apoptosis']

fig, axes = plt.subplots(2, 3, figsize=(18, 10))
axes = axes.flatten()

for i, mod in enumerate(selected_modules):
    if mod not in profiles_df.columns:
        continue

    ax = axes[i]

    for group in ['Upper', 'Glia', 'VAT1L']:
        group_data = profiles_df[profiles_df['group'] == group].sort_values('pt_center')

        if mod in group_data.columns:
            color = 'red' if group == 'Upper' else 'blue' if group == 'Glia' else 'green'
            ax.plot(group_data['pt_center'], group_data[mod],
                   marker='o', label=group, color=color, linewidth=2, markersize=4, alpha=0.7)

    ax.set_xlabel('PT_dpt', fontsize=9)
    ax.set_ylabel(f'{mod} expression', fontsize=9)
    ax.set_title(mod, fontsize=11, fontweight='bold')
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

plt.tight_layout()
plt.savefig(f'{output_dir}/Fig_module_trajectories.png', dpi=300, bbox_inches='tight')
print(f"Saved: {output_dir}/Fig_module_trajectories.png")
plt.close()

# ==============================================================================
# STEP 6: Statistical Summary and Report
# ==============================================================================

print("\n[6] Generating report...")

# Count modules by pattern
onset_stats = {
    'Upper_first_vs_Glia': (time_lags['Delta_onset_Upper_minus_Glia'] < 0).sum(),
    'Glia_first_vs_Upper': (time_lags['Delta_onset_Upper_minus_Glia'] > 0).sum(),
    'Upper_first_vs_VAT1L': (time_lags['Delta_onset_Upper_minus_VAT1L'] < 0).sum(),
    'VAT1L_first_vs_Upper': (time_lags['Delta_onset_Upper_minus_VAT1L'] > 0).sum(),
    'Total_modules': len(time_lags)
}

# Key modules for analysis
key_modules = {
    'Hyperexcitability_related': ['Synaptic', 'Calcium_Signaling', 'Ion_Transport'],
    'Stress_related': ['ER_Stress', 'Oxidative_Stress', 'Protein_Homeostasis'],
    'Inflammatory': ['Inflammation', 'Complement'],
    'Death_related': ['Apoptosis', 'Autophagy']
}

# Generate report
report = f"""# Phase 8: PT_dpt Time-Lag Analysis Report

**Date**: {pd.Timestamp.now().strftime('%Y-%m-%d')}
**Analysis**: Module onset/peak timing along PT_dpt pseudotime
**Dataset**: {len(df):,} cells across 3 groups (Upper, Glia, VAT1L)

---

## Executive Summary

**Goal**: Determine which cell type shows module changes FIRST along PT_dpt.

**Key Question**: Do Upper modules precede Glia/VAT1L changes (supporting Upper-first hypothesis),
or does Glia lead (supporting Glia-first), or VAT1L?

**Method**:
1. Binned cells by PT_dpt ({num_bins} bins from {pt_min:.6f} to {pt_max:.6f})
2. Calculated module averages for each group × bin
3. Identified "onset" PT (first bin where module > mean + 0.5*std)
4. Calculated time-lags: ΔPT (Upper onset - Glia/VAT1L onset)

**Interpretation**:
- **Negative ΔPT**: Upper module onset precedes Glia/VAT1L (Upper-first)
- **Positive ΔPT**: Glia/VAT1L module onset precedes Upper (Glia/VAT1L-first)

---

## Overall Statistics

### Upper vs Glia Timing

**Onset time-lags**:
- Mean ΔPT: {time_lags['Delta_onset_Upper_minus_Glia'].mean():.6f}
- Median ΔPT: {time_lags['Delta_onset_Upper_minus_Glia'].median():.6f}

**Module counts**:
- **Upper-first** (Δ<0): {onset_stats['Upper_first_vs_Glia']} modules ({100*onset_stats['Upper_first_vs_Glia']/onset_stats['Total_modules']:.1f}%)
- **Glia-first** (Δ>0): {onset_stats['Glia_first_vs_Upper']} modules ({100*onset_stats['Glia_first_vs_Upper']/onset_stats['Total_modules']:.1f}%)

---

### Upper vs VAT1L Timing

**Onset time-lags**:
- Mean ΔPT: {time_lags['Delta_onset_Upper_minus_VAT1L'].mean():.6f}
- Median ΔPT: {time_lags['Delta_onset_Upper_minus_VAT1L'].median():.6f}

**Module counts**:
- **Upper-first** (Δ<0): {onset_stats['Upper_first_vs_VAT1L']} modules ({100*onset_stats['Upper_first_vs_VAT1L']/onset_stats['Total_modules']:.1f}%)
- **VAT1L-first** (Δ>0): {onset_stats['VAT1L_first_vs_Upper']} modules ({100*onset_stats['VAT1L_first_vs_Upper']/onset_stats['Total_modules']:.1f}%)

---

## Key Module Time-Lags

### Hyperexcitability-Related Modules

"""

# Add key module details
for category, mods in key_modules.items():
    report += f"### {category}\n\n"
    report += "| Module | Upper Onset | Glia Onset | VAT1L Onset | Δ(U-G) | Δ(U-V) | Pattern |\n"
    report += "|--------|-------------|------------|-------------|--------|--------|----------|\n"

    for mod in mods:
        if mod not in time_lags['module'].values:
            continue

        row = time_lags[time_lags['module'] == mod].iloc[0]

        delta_ug = row['Delta_onset_Upper_minus_Glia']
        delta_uv = row['Delta_onset_Upper_minus_VAT1L']

        pattern_ug = 'U→G' if delta_ug < 0 else 'G→U'
        pattern_uv = 'U→V' if delta_uv < 0 else 'V→U'
        pattern = f"{pattern_ug}, {pattern_uv}"

        report += f"| {mod} | {row['Upper_onset_PT']:.5f} | {row['Glia_onset_PT']:.5f} | {row['VAT1L_onset_PT']:.5f} | {delta_ug:+.5f} | {delta_uv:+.5f} | {pattern} |\n"

    report += "\n"

report += """---

## Interpretation

### Patterns Observed

1. **Upper-first modules**: Modules where Upper onset precedes Glia/VAT1L
   - These support "Upper functional driver" hypothesis
   - Suggest Upper changes temporally precede Glia/VAT1L changes

2. **Glia-first modules**: Modules where Glia onset precedes Upper/VAT1L
   - These support "Glia metabolic upstream" hypothesis
   - Suggest Glia dysfunction precedes neuronal changes

3. **Mixed patterns**: Some modules show Upper→Glia but Glia→VAT1L
   - Suggest multi-step cascade

### Biological Interpretation

**If Hyperexcitability (Synaptic, Ca²⁺, Ion) shows Upper-first**:
→ Supports Upper neurons as early functional drivers
→ Excitotoxic mechanisms may propagate to Glia/VAT1L

**If ER Stress/Inflammation shows Glia-first**:
→ Supports Glia metabolic dysfunction as initiating event
→ Glia stress propagates to neurons

**If Mixed patterns**:
→ Supports dual-upstream model
→ Both Glia (metabolic) and Upper (functional) contribute

---

## Methodological Considerations

### Strengths

1. **PT_dpt is stress-independent**: R²=1.9% with stress_total
   - Avoids circular logic (unlike PT_imes with R²=62.5%)
2. **Quantitative time-lag measurement**: ΔPT provides numerical comparison
3. **Module-specific analysis**: Tests pathway-specific timing, not just aggregate stress

### Limitations

1. **Pseudotime ≠ real time**:
   - PT_dpt represents trajectory, not absolute time
   - "Onset" timing is relative, not absolute
   - Cannot distinguish fast vs slow processes

2. **Cross-sectional data**:
   - Infers dynamics from snapshot
   - True temporal precedence requires longitudinal data

3. **Binning artifacts**:
   - Onset detection depends on bin size and threshold
   - Sensitivity analysis needed (vary threshold 0.3-1.0σ)

4. **Cell density variation**:
   - Some bins have few cells (low n)
   - May introduce noise in low-density regions

5. **Interpretation caveats**:
   - "Upper-first" onset ≠ "Upper causes Glia changes"
   - Requires integration with Phase 7 correlational evidence
   - Spatial co-localization still needed for causality

---

## Integration with Phase 7 Findings

### Convergent Evidence for Upper Functional Driver

**Phase 7 (Correlational)**:
- Upper ER → Glia ER: r=0.664, p=0.010
- Upper Synaptic → Glia Inflammation: r=0.534, p=0.049
- Upper heterogeneity → VAT1L stress: r=0.569, p=0.034

**Phase 8 (Temporal)**:
- If ER Stress shows Upper-first onset → temporal precedence
- If Synaptic shows Upper-first onset → temporal precedence
- Combines with Phase 7 correlations → stronger driver hypothesis support

### Complementarity

- **Phase 7**: Tests "if Upper high, then Glia/VAT1L high" (correlation)
- **Phase 8**: Tests "Upper changes before Glia/VAT1L" (temporal sequence)
- **Together**: Both correlation AND temporal precedence → stronger causal inference

---

## Conclusions

### Main Findings

1. **Onset timing statistics**:
   - Upper vs Glia: {onset_stats['Upper_first_vs_Glia']}/{onset_stats['Total_modules']} modules show Upper-first
   - Upper vs VAT1L: {onset_stats['Upper_first_vs_VAT1L']}/{onset_stats['Total_modules']} modules show Upper-first

2. **Key module patterns**: (See tables above)

3. **Interpretation**:
   - Pattern suggests [TO BE DETERMINED after reviewing module-specific results]
   - Integration with Phase 7 correlational evidence needed

### Next Steps

1. **Sensitivity analysis**:
   - Test different onset thresholds (0.3σ, 0.5σ, 1.0σ)
   - Test different bin numbers (10, 20, 30)
   - Assess robustness of findings

2. **Module grouping**:
   - Group modules by functional category
   - Test if categories show consistent patterns

3. **Integration with correlations**:
   - Modules with both Upper-first onset AND strong Phase 7 correlation
   - Highest confidence driver candidates

4. **Validation**:
   - Longitudinal data (if available)
   - Spatial transcriptomics (test local co-progression)
   - Perturbation experiments (causal validation)

---

**Generated**: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}
**Scripts**: `scripts/Phase8_PT_dpt_time_lag_analysis.py`
**Data**: `results/phase8_time_lag/`

**Key Files**:
- `module_profiles_by_group_and_bin.csv`: Group × bin × module averages
- `module_onset_peak_PT_by_group.csv`: Onset/peak PT for each module
- `module_time_lags.csv`: ΔPT statistics
- `Fig_module_onset_time_lags.png`: Time-lag visualization
- `Fig_module_trajectories.png`: Module expression trajectories
"""

# Save report
with open(f'{output_dir}/PHASE8_TIME_LAG_REPORT.md', 'w') as f:
    f.write(report)

print(f"Saved: {output_dir}/PHASE8_TIME_LAG_REPORT.md")

print("\n" + "="*80)
print("Phase 8 Analysis Complete!")
print("="*80)
print(f"\nResults saved to: {output_dir}/")
print("\nKey outputs:")
print(f"  - module_profiles_by_group_and_bin.csv")
print(f"  - module_onset_peak_PT_by_group.csv")
print(f"  - module_time_lags.csv")
print(f"  - Fig_module_onset_time_lags.png")
print(f"  - Fig_module_trajectories.png")
print(f"  - PHASE8_TIME_LAG_REPORT.md")
