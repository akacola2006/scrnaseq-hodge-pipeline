#!/usr/bin/env python3
"""
Phase 13: NVU Energy Landscape per Module (φ-flow)

For each functional module (Mitochondria, ER_Stress, Synaptic, Calcium_Signaling, etc.),
compute an "energy curve" φ_m(node, PT_dpt) for each NVU node (Vascular, Glia, Upper, VAT1L)
and visualize the disease-space landscape.

Goals:
  - See which modules are activated earliest in which node
  - Detect shared flow directions (e.g. Mito: Vasc→Glia→Upper→VAT1L)
  - Identify modules where NVU nodes behave differently

Author: Claude Code
Date: 2025-11-24
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# ============================================================
# Configuration
# ============================================================

# Input file
INPUT_FILE = 'results/phase9_vascular/cell_level_features_ALL_with_PTdpt.csv'

# Output directory
OUTPUT_DIR = 'results/phase13_nvu_energy'
Path(OUTPUT_DIR).mkdir(parents=True, exist_ok=True)

# NVU node definitions
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

# Key energy modules to visualize
ENERGY_MODULES = [
    'Mitochondria',
    'ER_Stress',
    'Synaptic',
    'Calcium_Signaling',
    'Complement',
    'Oxidative_Stress',
    'Inflammation',
    'Protein_Homeostasis'
]

# PT_dpt binning parameters
NUM_BINS = 30
PHI_THRESHOLD = 0.5  # Threshold for determining onset

# Colors for NVU nodes
NODE_COLORS = {
    'Vascular': '#87CEEB',  # Sky blue
    'Glia': '#98D8C8',      # Mint green
    'Upper': '#FFB347',     # Orange
    'VAT1L': '#DDA0DD'      # Plum
}

print("=" * 80)
print("Phase 13: NVU Energy Landscape per Module (φ-flow)")
print("=" * 80)
print()

# ============================================================
# Step 1: Load data and define NVU nodes
# ============================================================

print("[1] Loading data and defining NVU nodes...")

df = pd.read_csv(INPUT_FILE)
print(f"  Loaded {len(df)} cells")

# Add NVU node column
def assign_nvu_node(cell_type):
    for node, subtypes in NVU_NODES.items():
        if cell_type in subtypes:
            return node
    return None

df['nvu_node'] = df['cell_type'].apply(assign_nvu_node)

# Filter to cells with NVU node assignment
df_nvu = df[df['nvu_node'].notna()].copy()
print(f"  Cells in NVU nodes: {len(df_nvu)}")

# Count by node and condition
node_counts = df_nvu.groupby(['nvu_node', 'condition']).size().unstack(fill_value=0)
print("\n  NVU node counts:")
print(node_counts)
print()

# ============================================================
# Step 2: Calculate φ_m (energy) for each module
# ============================================================

print("[2] Calculating φ_m (energy) from Control deviation...")

# Check which modules exist in data
available_modules = []
for module in ENERGY_MODULES:
    col_name = f'module_{module}'
    if col_name in df_nvu.columns:
        available_modules.append(module)
    else:
        print(f"  Warning: {col_name} not found in data")

print(f"  Available modules: {len(available_modules)}")
print(f"    {', '.join(available_modules)}")
print()

# Calculate φ_m = Z² from Control for each module
module_stats = {}

# First get Control reference
df_ctrl = df_nvu[df_nvu['condition'] == 'Control'].copy()
print(f"  Control cells in NVU: {len(df_ctrl)}")

for module in available_modules:
    col_name = f'module_{module}'

    # Compute Control mean and std
    ctrl_vals = df_ctrl[col_name].dropna()
    mu = ctrl_vals.mean()
    sigma = ctrl_vals.std()

    if sigma == 0 or np.isnan(sigma):
        print(f"  Warning: {module} has zero/NaN std in Control, skipping")
        continue

    # Store stats
    module_stats[module] = {'mu': mu, 'sigma': sigma}

    # Calculate Z-score and φ = Z² for all cells
    Z = (df_nvu[col_name] - mu) / sigma
    df_nvu[f'phi_{module}'] = Z ** 2

    print(f"  {module:25s}: μ={mu:.3f}, σ={sigma:.3f}")

print()

# Now filter to ALS cells (after calculating φ columns)
df_als = df_nvu[df_nvu['condition'] == 'ALS'].copy()
print(f"  ALS cells in NVU: {len(df_als)}")
print()

# ============================================================
# Step 3: Aggregate φ_m by NVU node × PT_dpt bin
# ============================================================

print("[3] Aggregating φ_m by NVU node × PT_dpt bins...")

# Define PT_dpt bins (ALS only)
pt_min = df_als['PT_dpt'].min()
pt_max = df_als['PT_dpt'].max()
bins = np.linspace(pt_min, pt_max, NUM_BINS + 1)
bin_centers = (bins[:-1] + bins[1:]) / 2

print(f"  PT_dpt range (ALS): [{pt_min:.3f}, {pt_max:.3f}]")
print(f"  Number of bins: {NUM_BINS}")
print()

# Assign bin index to ALS cells
df_als['pt_bin'] = pd.cut(df_als['PT_dpt'], bins=bins, labels=False, include_lowest=True)

# Aggregate φ_m by node × bin × module
phi_profiles = []

for module in available_modules:
    if module not in module_stats:
        continue

    phi_col = f'phi_{module}'

    for node in ['Vascular', 'Glia', 'Upper', 'VAT1L']:
        df_node = df_als[df_als['nvu_node'] == node].copy()

        if len(df_node) == 0:
            continue

        # Aggregate by bin
        for bin_idx in range(NUM_BINS):
            df_bin = df_node[df_node['pt_bin'] == bin_idx]

            if len(df_bin) == 0:
                continue

            phi_mean = df_bin[phi_col].mean()
            phi_std = df_bin[phi_col].std()
            n_cells = len(df_bin)
            pt_center = bin_centers[bin_idx]

            phi_profiles.append({
                'node': node,
                'module': module,
                'pt_bin': bin_idx,
                'pt_bin_center': pt_center,
                'phi_mean': phi_mean,
                'phi_std': phi_std if not np.isnan(phi_std) else 0,
                'n_cells': n_cells
            })

df_phi = pd.DataFrame(phi_profiles)
print(f"  Generated {len(df_phi)} node × bin × module profiles")

# Save
output_file = f'{OUTPUT_DIR}/nvu_phi_profiles_by_node_and_bin.csv'
df_phi.to_csv(output_file, index=False)
print(f"  Saved: {output_file}")
print()

# ============================================================
# Step 4: Visualize φ_m(PT) curves for each module
# ============================================================

print("[4] Visualizing φ_m(PT) energy curves...")

for module in available_modules:
    if module not in module_stats:
        continue

    fig, ax = plt.subplots(figsize=(12, 6))

    # Plot each node
    for node in ['Vascular', 'Glia', 'Upper', 'VAT1L']:
        df_node = df_phi[(df_phi['module'] == module) & (df_phi['node'] == node)]

        if len(df_node) == 0:
            continue

        # Sort by PT
        df_node = df_node.sort_values('pt_bin_center')

        # Plot line
        ax.plot(df_node['pt_bin_center'], df_node['phi_mean'],
               label=node, color=NODE_COLORS[node], linewidth=2.5, marker='o', markersize=4)

        # Optional: shaded area for ±1 std
        ax.fill_between(df_node['pt_bin_center'],
                       df_node['phi_mean'] - df_node['phi_std'],
                       df_node['phi_mean'] + df_node['phi_std'],
                       color=NODE_COLORS[node], alpha=0.2)

    # Add threshold line
    ax.axhline(y=PHI_THRESHOLD, color='gray', linestyle='--', linewidth=1, alpha=0.5, label='Onset threshold')

    ax.set_xlabel('PT_dpt (Disease-space coordinate)', fontsize=12)
    ax.set_ylabel(f'φ_{{{module}}} (Energy: Z² from Control)', fontsize=12)
    ax.set_title(f'Phase 13: Energy Flow for {module} Module Across NVU Nodes',
                fontsize=14, fontweight='bold')
    ax.legend(loc='best', frameon=True, fancybox=True, shadow=True)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    output_file = f'{OUTPUT_DIR}/Fig_phi_curve_{module}.png'
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    plt.close()

    print(f"  Saved: Fig_phi_curve_{module}.png")

print()

# ============================================================
# Step 5: Compute flow order and onset/peak ranks
# ============================================================

print("[5] Computing flow order (onset and peak ranks)...")

flow_summary = []

for module in available_modules:
    if module not in module_stats:
        continue

    node_metrics = []

    for node in ['Vascular', 'Glia', 'Upper', 'VAT1L']:
        df_node = df_phi[(df_phi['module'] == module) & (df_phi['node'] == node)]

        if len(df_node) == 0:
            continue

        df_node = df_node.sort_values('pt_bin_center')

        # Find onset_PT (first bin where φ_mean > threshold)
        onset_bins = df_node[df_node['phi_mean'] > PHI_THRESHOLD]
        if len(onset_bins) > 0:
            onset_pt = onset_bins.iloc[0]['pt_bin_center']
        else:
            onset_pt = np.nan

        # Find peak_PT (bin with maximum φ_mean)
        peak_idx = df_node['phi_mean'].idxmax()
        peak_pt = df_node.loc[peak_idx, 'pt_bin_center']
        peak_phi = df_node.loc[peak_idx, 'phi_mean']

        node_metrics.append({
            'node': node,
            'onset_pt': onset_pt,
            'peak_pt': peak_pt,
            'peak_phi': peak_phi
        })

    # Rank nodes by onset_pt and peak_pt
    if len(node_metrics) > 0:
        df_metrics = pd.DataFrame(node_metrics)

        # Onset rank (lower PT = rank 1)
        df_metrics['onset_rank'] = df_metrics['onset_pt'].rank(method='min', na_option='bottom')

        # Peak rank
        df_metrics['peak_rank'] = df_metrics['peak_pt'].rank(method='min')

        # Add module
        df_metrics['module'] = module

        flow_summary.append(df_metrics)

df_flow = pd.concat(flow_summary, ignore_index=True)

# Reorder columns
df_flow = df_flow[['module', 'node', 'onset_pt', 'peak_pt', 'peak_phi', 'onset_rank', 'peak_rank']]

# Save
output_file = f'{OUTPUT_DIR}/module_flow_summary.csv'
df_flow.to_csv(output_file, index=False)
print(f"  Saved: {output_file}")

# Create flow order table (ordered list by onset_pt)
flow_order = []

for module in available_modules:
    if module not in module_stats:
        continue

    df_mod = df_flow[df_flow['module'] == module].copy()
    df_mod = df_mod.sort_values('onset_rank')

    # Create ordered list
    order_str = ' → '.join(df_mod['node'].values)

    flow_order.append({
        'module': module,
        'flow_order_by_onset': order_str,
        'earliest_node': df_mod.iloc[0]['node'] if len(df_mod) > 0 else None,
        'earliest_onset_PT': df_mod.iloc[0]['onset_pt'] if len(df_mod) > 0 else None
    })

df_order = pd.DataFrame(flow_order)

output_file = f'{OUTPUT_DIR}/module_flow_order.csv'
df_order.to_csv(output_file, index=False)
print(f"  Saved: {output_file}")
print()

# Print summary
print("  Flow order by module (onset-based):")
for _, row in df_order.iterrows():
    print(f"    {row['module']:25s}: {row['flow_order_by_onset']}")

print()

# ============================================================
# Step 6: Create heatmap of onset ranks
# ============================================================

print("[6] Creating heatmap visualization...")

# Pivot onset_rank
pivot_onset = df_flow.pivot(index='module', columns='node', values='onset_rank')
pivot_onset = pivot_onset[['Vascular', 'Glia', 'Upper', 'VAT1L']]  # Order columns

# Create heatmap
fig, ax = plt.subplots(figsize=(10, 8))

sns.heatmap(pivot_onset, annot=True, fmt='.0f', cmap='RdYlGn_r',
           cbar_kws={'label': 'Onset Rank (1=earliest)'},
           linewidths=1, linecolor='gray', ax=ax)

ax.set_title('Phase 13: Module-Specific Flow Order\n(Onset Rank by NVU Node)',
            fontsize=14, fontweight='bold')
ax.set_xlabel('NVU Node', fontsize=12)
ax.set_ylabel('Functional Module', fontsize=12)

plt.tight_layout()
output_file = f'{OUTPUT_DIR}/Fig_onset_rank_heatmap.png'
plt.savefig(output_file, dpi=150, bbox_inches='tight')
plt.close()

print(f"  Saved: {output_file}")
print()

# ============================================================
# Completion
# ============================================================

print("=" * 80)
print("Phase 13 Complete!")
print("=" * 80)
print()
print("Key outputs:")
print(f"  - {OUTPUT_DIR}/nvu_phi_profiles_by_node_and_bin.csv")
print(f"  - {OUTPUT_DIR}/module_flow_summary.csv")
print(f"  - {OUTPUT_DIR}/module_flow_order.csv")
print(f"  - {OUTPUT_DIR}/Fig_phi_curve_{{module}}.png (×{len(available_modules)})")
print(f"  - {OUTPUT_DIR}/Fig_onset_rank_heatmap.png")
print()
print("Next: Create Phase 13 comprehensive report")
