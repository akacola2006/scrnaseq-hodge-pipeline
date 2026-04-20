#!/usr/bin/env python3
"""
GSE_212630: φ-Flow Analysis with PT_dpt (Diffusion Pseudotime)
===============================================================

This script performs φ-flow analysis using PT_dpt for comparison with PT_pca.
PT_dpt is based on stress component diffusion, while PT_pca uses key ALS genes.

Author: Claude Code + User
Date: 2025-12-07
"""

import pandas as pd
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import sys
import warnings
warnings.filterwarnings('ignore')

# Import IDS Core
SCRIPT_DIR = Path(__file__).parent
sys.path.insert(0, str(SCRIPT_DIR))
import ids_core

print("=" * 80)
print("GSE_212630: φ-Flow Analysis with PT_dpt (Diffusion Pseudotime)")
print("=" * 80)

# ============================================================================
# CONFIGURATION
# ============================================================================

INPUT_FILE = SCRIPT_DIR.parent / 'results' / 'GSE212630_ids_analysis' / 'combined_metadata_with_PT_dpt.csv'
OUTPUT_DIR = SCRIPT_DIR.parent / 'results' / 'GSE212630_ids_analysis'

# Use PT_dpt as pseudotime
PT_COL = 'PT_dpt'

# NVU node mapping
NVU_MAPPING = {
    'Vascular': ['Vasc.Capillary', 'Vasc.Endo', 'Vasc.Fibro', 'Vasc.Pericyte', 'Vasc.SMC', 'Vasc.Unknown'],
    'Glia': ['Glia.Oligo', 'Glia.Astro.GFAP.neg', 'Glia.Astro.GFAP.pos', 'Glia.OPC', 'Glia.Micro'],
    'Excitatory': ['Ex.L2/3', 'Ex.L4', 'Ex.L5', 'Ex.L6', 'Ex.Unknown'],
    'Inhibitory': ['In.SST', 'In.VIP', 'In.DISC1', 'In.PVALB', 'In.Unknown']
}

NVU_ORDER = ['Vascular', 'Glia', 'Excitatory', 'Inhibitory']

# Analysis parameters
N_BINS = 10
PHI_THRESHOLD = 0.5

# ============================================================================
# STEP 1: LOAD DATA
# ============================================================================

print("\n[Step 1] Loading data...")
df = pd.read_csv(INPUT_FILE)
print(f"  Total cells: {len(df):,}")

# Check PT_dpt
print(f"\n  PT_dpt range: [{df[PT_COL].min():.4f}, {df[PT_COL].max():.4f}]")
print(f"  PT_dpt mean: {df[PT_COL].mean():.4f}")

# Map to NVU nodes
def map_to_nvu(cell_type):
    for node, types in NVU_MAPPING.items():
        if cell_type in types:
            return node
    return 'Other'

df['nvu_node'] = df['cell_type'].apply(map_to_nvu)
print(f"\n  NVU node distribution:")
for node in NVU_ORDER:
    n = (df['nvu_node'] == node).sum()
    print(f"    {node}: {n:,}")

# ============================================================================
# STEP 2: COMPUTE φ FOR STRESS COMPONENTS
# ============================================================================

print("\n[Step 2] Computing φ for stress components...")

control_mask = pd.Series((df['condition'] == 'Control').values, index=df.index)

stress_cols = ['stress_total'] + [col for col in df.columns if col.startswith('stress_') and col != 'stress_total']
print(f"  Stress components: {stress_cols}")

for col in stress_cols:
    scores = pd.Series(df[col].values, index=df.index)
    phi = ids_core.compute_phi_from_scores(scores, control_mask)
    df[f'phi_{col}'] = phi.values

print(f"  Computed φ for {len(stress_cols)} stress components")

# ============================================================================
# STEP 3: BIN BY PT_dpt AND COMPUTE φ PROFILES
# ============================================================================

print("\n[Step 3] Binning by PT_dpt and computing φ profiles...")

pt_series = pd.Series(df[PT_COL].values, index=df.index)
bin_index, bin_edges, bin_centers = ids_core.bin_by_pt(pt_series, n_bins=N_BINS)
df['pt_bin'] = bin_index.values

print(f"  Number of bins: {N_BINS}")
print(f"  Bin centers: {bin_centers[:3]}... {bin_centers[-3:]}")

# Compute φ profiles for stress_total by NVU node
phi_series = pd.Series(df['phi_stress_total'].values, index=df.index)
group_series = pd.Series(df['nvu_node'].values, index=df.index)

profiles = ids_core.summarize_phi_by_group_and_bin(
    phi=phi_series,
    group=group_series,
    bin_index=bin_index,
    group_order=NVU_ORDER
)

print(f"  Profile records: {len(profiles)}")
profiles.to_csv(OUTPUT_DIR / 'phi_profiles_PT_dpt.csv', index=False)
print(f"  Saved: phi_profiles_PT_dpt.csv")

# ============================================================================
# STEP 4: DETECT ONSET AND PEAK
# ============================================================================

print("\n[Step 4] Detecting onset and peak...")

flow_summary = ids_core.detect_onset_and_peak(
    phi_mean_df=profiles[['group', 'bin', 'phi_mean']].copy(),
    bin_centers=bin_centers,
    threshold=PHI_THRESHOLD
)

print("\n  Flow Summary (PT_dpt):")
print(flow_summary.to_string(index=False))

flow_summary.to_csv(OUTPUT_DIR / 'flow_summary_PT_dpt.csv', index=False)
print(f"\n  Saved: flow_summary_PT_dpt.csv")

# ============================================================================
# STEP 5: CELL TYPE ANALYSIS
# ============================================================================

print("\n[Step 5] Cell type-level analysis...")

late_pt_mask = df['pt_bin'] >= N_BINS - 3
late_phi = df[late_pt_mask].groupby('cell_type')['phi_stress_total'].mean().sort_values(ascending=False)

print("\n  Top 5 cell types with highest φ at late PT_dpt:")
for ct, phi in late_phi.head(5).items():
    print(f"    {ct}: φ = {phi:.3f}")

# ============================================================================
# STEP 6: GENERATE FIGURES
# ============================================================================

print("\n[Step 6] Generating figures...")

fig, axes = plt.subplots(2, 2, figsize=(14, 12))

colors = {'Vascular': '#e74c3c', 'Glia': '#f39c12', 'Excitatory': '#27ae60', 'Inhibitory': '#3498db'}

# Panel A: φ trajectory
ax = axes[0, 0]
for node in NVU_ORDER:
    node_data = profiles[profiles['group'] == node]
    ax.plot(node_data['bin'], node_data['phi_mean'], 'o-',
            label=node, color=colors[node], linewidth=2, markersize=8)
    ax.fill_between(node_data['bin'],
                    node_data['phi_mean'] - node_data['phi_std']/2,
                    node_data['phi_mean'] + node_data['phi_std']/2,
                    alpha=0.2, color=colors[node])

ax.set_xlabel('PT_dpt Bin (Diffusion Pseudotime →)', fontsize=12)
ax.set_ylabel('φ (Stress Energy)', fontsize=12)
ax.set_title('φ Trajectory by NVU Node\n(PT_dpt: Stress-based Diffusion Pseudotime)', fontsize=12)
ax.legend(loc='upper left')
ax.grid(True, alpha=0.3)

# Panel B: Heatmap
ax = axes[0, 1]
pivot = profiles.pivot(index='group', columns='bin', values='phi_mean')
pivot = pivot.reindex(NVU_ORDER)
sns.heatmap(pivot, cmap='YlOrRd', ax=ax, annot=True, fmt='.2f',
            cbar_kws={'label': 'φ (Stress Energy)'})
ax.set_xlabel('PT_dpt Bin', fontsize=12)
ax.set_ylabel('NVU Node', fontsize=12)
ax.set_title('φ Heatmap by NVU Node × PT_dpt Bin', fontsize=12)

# Panel C: Onset timing
ax = axes[1, 0]
onset_data = flow_summary[flow_summary['onset_PT'].notna()].copy()
onset_data = onset_data.sort_values('onset_PT')
bars = ax.barh(onset_data['group'], onset_data['onset_PT'],
               color=[colors.get(g, 'gray') for g in onset_data['group']])
ax.set_xlabel('Onset PT_dpt', fontsize=12)
ax.set_ylabel('NVU Node', fontsize=12)
ax.set_title('Onset Timing by NVU Node (PT_dpt)', fontsize=12)
ax.grid(True, alpha=0.3, axis='x')
for i, (_, row) in enumerate(onset_data.iterrows()):
    ax.text(row['onset_PT'] + 0.01, i, f"{row['onset_PT']:.3f}", va='center', fontsize=10)

# Panel D: Peak φ
ax = axes[1, 1]
peak_data = flow_summary.copy()
peak_data = peak_data.sort_values('peak_phi', ascending=True)
bars = ax.barh(peak_data['group'], peak_data['peak_phi'],
               color=[colors.get(g, 'gray') for g in peak_data['group']])
ax.set_xlabel('Peak φ', fontsize=12)
ax.set_ylabel('NVU Node', fontsize=12)
ax.set_title('Peak φ by NVU Node (PT_dpt)', fontsize=12)
ax.grid(True, alpha=0.3, axis='x')
for i, (_, row) in enumerate(peak_data.iterrows()):
    ax.text(row['peak_phi'] + 0.05, i, f"{row['peak_phi']:.2f}", va='center', fontsize=10)

plt.tight_layout()
fig.savefig(OUTPUT_DIR / 'Fig_phi_flow_PT_dpt_main.png', dpi=150, bbox_inches='tight')
plt.close()
print(f"  Saved: Fig_phi_flow_PT_dpt_main.png")

# Cell type heatmap
fig, ax = plt.subplots(figsize=(14, 10))
ct_pivot = df.groupby(['cell_type', 'pt_bin'])['phi_stress_total'].mean().unstack()
sns.heatmap(ct_pivot, cmap='YlOrRd', ax=ax, annot=True, fmt='.2f',
            cbar_kws={'label': 'φ (Stress Energy)'})
ax.set_xlabel('PT_dpt Bin', fontsize=12)
ax.set_ylabel('Cell Type', fontsize=12)
ax.set_title('φ Heatmap by Cell Type × PT_dpt Bin', fontsize=14)
plt.tight_layout()
fig.savefig(OUTPUT_DIR / 'Fig_phi_celltype_heatmap_PT_dpt.png', dpi=150, bbox_inches='tight')
plt.close()
print(f"  Saved: Fig_phi_celltype_heatmap_PT_dpt.png")

# ============================================================================
# STEP 7: COMPARISON WITH PT_pca
# ============================================================================

print("\n[Step 7] Loading PT_pca results for comparison...")

try:
    pca_flow = pd.read_csv(OUTPUT_DIR / 'flow_summary_PT_pca_all_components.csv')
    pca_flow = pca_flow[pca_flow['stress_component'] == 'stress_total']

    print("\n  COMPARISON: PT_pca vs PT_dpt")
    print("  " + "-" * 60)
    print(f"  {'NVU Node':<12} | {'PT_pca Onset':<12} | {'PT_dpt Onset':<12} | {'PT_pca Peak':<10} | {'PT_dpt Peak':<10}")
    print("  " + "-" * 60)

    for node in NVU_ORDER:
        pca_row = pca_flow[pca_flow['group'] == node]
        dpt_row = flow_summary[flow_summary['group'] == node]

        pca_onset = pca_row['onset_PT'].values[0] if len(pca_row) > 0 else np.nan
        dpt_onset = dpt_row['onset_PT'].values[0] if len(dpt_row) > 0 else np.nan
        pca_peak = pca_row['peak_phi'].values[0] if len(pca_row) > 0 else np.nan
        dpt_peak = dpt_row['peak_phi'].values[0] if len(dpt_row) > 0 else np.nan

        print(f"  {node:<12} | {pca_onset:<12.3f} | {dpt_onset:<12.3f} | {pca_peak:<10.2f} | {dpt_peak:<10.2f}")

    # Create comparison figure
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Onset comparison
    ax = axes[0]
    x = np.arange(len(NVU_ORDER))
    width = 0.35

    pca_onsets = [pca_flow[pca_flow['group'] == node]['onset_PT'].values[0] if len(pca_flow[pca_flow['group'] == node]) > 0 else 0 for node in NVU_ORDER]
    dpt_onsets = [flow_summary[flow_summary['group'] == node]['onset_PT'].values[0] if len(flow_summary[flow_summary['group'] == node]) > 0 else 0 for node in NVU_ORDER]

    ax.bar(x - width/2, pca_onsets, width, label='PT_pca', color='#3498db')
    ax.bar(x + width/2, dpt_onsets, width, label='PT_dpt', color='#e74c3c')
    ax.set_ylabel('Onset PT')
    ax.set_title('Onset Timing Comparison')
    ax.set_xticks(x)
    ax.set_xticklabels(NVU_ORDER)
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Peak φ comparison
    ax = axes[1]
    pca_peaks = [pca_flow[pca_flow['group'] == node]['peak_phi'].values[0] if len(pca_flow[pca_flow['group'] == node]) > 0 else 0 for node in NVU_ORDER]
    dpt_peaks = [flow_summary[flow_summary['group'] == node]['peak_phi'].values[0] if len(flow_summary[flow_summary['group'] == node]) > 0 else 0 for node in NVU_ORDER]

    ax.bar(x - width/2, pca_peaks, width, label='PT_pca', color='#3498db')
    ax.bar(x + width/2, dpt_peaks, width, label='PT_dpt', color='#e74c3c')
    ax.set_ylabel('Peak φ')
    ax.set_title('Peak φ Comparison')
    ax.set_xticks(x)
    ax.set_xticklabels(NVU_ORDER)
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    fig.savefig(OUTPUT_DIR / 'Fig_PT_pca_vs_PT_dpt_comparison.png', dpi=150, bbox_inches='tight')
    plt.close()
    print(f"\n  Saved: Fig_PT_pca_vs_PT_dpt_comparison.png")

except Exception as e:
    print(f"  [WARN] Could not load PT_pca results: {e}")

# ============================================================================
# SUMMARY
# ============================================================================

print("\n" + "=" * 80)
print("φ-FLOW ANALYSIS (PT_dpt) COMPLETE")
print("=" * 80)

print(f"""
RESULTS (PT_dpt - Stress-based Diffusion Pseudotime):

1. ONSET TIMING:
""")
for _, row in flow_summary.iterrows():
    onset = f"{row['onset_PT']:.3f}" if pd.notna(row['onset_PT']) else 'N/A'
    print(f"   - {row['group']}: onset PT = {onset}")

print(f"""
2. PEAK φ VALUES:
""")
for _, row in flow_summary.iterrows():
    print(f"   - {row['group']}: peak φ = {row['peak_phi']:.3f} at PT = {row['peak_PT']:.3f}")

print(f"""
3. TOP AFFECTED CELL TYPES (Late PT_dpt):
""")
for ct, phi in late_phi.head(5).items():
    print(f"   - {ct}: φ = {phi:.3f}")

print(f"""
Output: {OUTPUT_DIR}
""")
print("=" * 80)
