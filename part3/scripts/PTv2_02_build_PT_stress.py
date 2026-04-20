#!/usr/bin/env python3
"""
Phase PTv2 - Step 2: Stress-based Pseudotime (PT_stress) Construction
=====================================================================

Purpose:
- Build a NEW pseudotime (PT_stress) based on STRESS POTENTIAL GRADIENT
- Independent of both PT_imes and PT_dpt
- Use directed graph where edges flow from low→high stress

Input:
- cell_level_features_ALL_with_PT_dpt.csv

Output:
- cell_level_features_ALL_with_PT_dpt_stress.csv
- PT_stress_QC_plots.png
"""

import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.neighbors import NearestNeighbors
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import shortest_path
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

print("="*80)
print("Phase PTv2 - Step 2: Building PT_stress (Stress-based Pseudotime)")
print("="*80)

# Paths
base_dir = Path('/home/akaco/als/motor_cortex_analysis')
ids_dir = base_dir / 'ids_causal_analysis'
input_file = ids_dir / 'results' / 'PTv2_robustness' / 'cell_level_features_ALL_with_PT_dpt.csv'
output_dir = ids_dir / 'results' / 'PTv2_robustness'

# ============================================================================
# Load Data
# ============================================================================
print("\n[Loading data...]")

df = pd.read_csv(input_file)
print(f"  Total cells: {len(df):,}")
print(f"  Columns: {len(df.columns)}")

# ============================================================================
# Step 2-1: Define Composite Stress Potential
# ============================================================================
print("\n[Step 2-1] Defining composite stress potential...")

# Define stress-related modules with weights
stress_modules = {
    'module_Angiogenesis': 1.5,        # Vascular stress
    'module_Oxidative_Stress': 2.0,    # High weight - core stress
    'module_Mitochondria': 1.8,        # Metabolic stress
    'module_ER_Stress': 2.0,           # High weight - core stress
    'module_Inflammation': 1.5,        # Immune stress
    'module_Cytoskeleton': 1.0,        # Structural stress
    'module_Apoptosis': 1.2,           # Death pathway
    'module_Autophagy': 0.8,           # Protective (inverse?)
    'module_DNA_Repair': 1.0,          # Genomic stress
    'module_Protein_Homeostasis': 1.5  # Proteostasis stress
}

# Check which modules are available
available_stress_modules = {k: v for k, v in stress_modules.items() if k in df.columns}
print(f"  Using {len(available_stress_modules)} stress modules")

# Calculate composite stress score
S_stress_base = np.zeros(len(df))

for module, weight in available_stress_modules.items():
    S_stress_base += weight * df[module].values
    print(f"    {module}: weight={weight}")

# Normalize
S_stress_base = (S_stress_base - S_stress_base.min()) / (S_stress_base.max() - S_stress_base.min())

df['S_stress_base'] = S_stress_base

print(f"\n  S_stress_base statistics:")
print(f"    Mean: {S_stress_base.mean():.4f}")
print(f"    Std: {S_stress_base.std():.4f}")
print(f"    Range: [{S_stress_base.min():.4f}, {S_stress_base.max():.4f}]")

# Correlation with existing stress_total
corr_stress = np.corrcoef(S_stress_base, df['stress_total'])[0, 1]
print(f"    Correlation with stress_total: {corr_stress:.3f}")

# ============================================================================
# Step 2-2: Build Stress-based Directed Graph
# ============================================================================
print("\n[Step 2-2] Building stress-based directed graph...")

# 2-2a. Build kNN graph in feature space (same as PT_dpt)
print("\n  [2-2a] Building kNN graph...")

module_cols = [col for col in df.columns if col.startswith('module_')]
top_modules = [
    'module_Angiogenesis',
    'module_Synaptic',
    'module_Mitochondria',
    'module_Apoptosis',
    'module_ER_Stress',
    'module_Oxidative_Stress',
    'module_Inflammation',
    'module_Cytoskeleton',
    'module_Protein_Homeostasis',
    'module_RNA_Processing',
    'module_Transcription',
    'module_Metabolism'
]
top_modules = [m for m in top_modules if m in module_cols]

feature_cols = top_modules + ['stress_total']
X = df[feature_cols].values

# Standardize
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

# Build kNN
K = 30
nbrs = NearestNeighbors(n_neighbors=K+1, metric='euclidean', n_jobs=-1)
nbrs.fit(X_scaled)
distances, indices = nbrs.kneighbors(X_scaled)

print(f"      kNN graph: K={K}, {len(df):,} nodes")

# 2-2b. Create directed edges based on stress gradient
print("\n  [2-2b] Creating directed edges (low→high stress)...")

n_cells = len(df)
row_idx = []
col_idx = []
edge_weights = []

n_forward = 0   # Edges where stress increases (i→j: S_j > S_i)
n_reverse = 0   # Edges where stress decreases (i→j: S_j < S_i)
n_kept = 0      # Edges kept after filtering

for i in range(n_cells):
    S_i = S_stress_base[i]

    for j_idx, j in enumerate(indices[i, 1:]):  # Skip self
        S_j = S_stress_base[j]
        delta_S = S_j - S_i

        # Count direction
        if delta_S > 0:
            n_forward += 1
        else:
            n_reverse += 1

        # Only keep edges where stress increases (forward direction)
        # OR allow small reverse edges with reduced weight
        if delta_S > 0:
            # Forward edge: full weight
            weight = 1.0 / (distances[i, j_idx + 1] + 1e-6)
            row_idx.append(i)
            col_idx.append(j)
            edge_weights.append(weight)
            n_kept += 1
        elif delta_S > -0.05:  # Allow small reverse with very low weight
            # Small reverse: reduced weight
            weight = 0.1 / (distances[i, j_idx + 1] + 1e-6)
            row_idx.append(i)
            col_idx.append(j)
            edge_weights.append(weight)
            n_kept += 1

print(f"      Forward edges (ΔS > 0): {n_forward:,}")
print(f"      Reverse edges (ΔS < 0): {n_reverse:,}")
print(f"      Kept edges: {n_kept:,} ({n_kept/(n_forward+n_reverse)*100:.1f}%)")

# Create directed graph (sparse matrix)
W_directed = csr_matrix((edge_weights, (row_idx, col_idx)), shape=(n_cells, n_cells))
print(f"      Directed graph: {W_directed.shape}, {W_directed.nnz:,} non-zero entries")

# 2-2c. Calculate geodesic distance from root
print("\n  [2-2c] Calculating geodesic distances...")

# Root: lowest stress cell
root_idx = S_stress_base.argmin()
root_cell = df.iloc[root_idx]
print(f"      Root cell: {root_cell['cell_id']}")
print(f"      Root cell_type: {root_cell['cell_type']}")
print(f"      Root S_stress_base: {S_stress_base[root_idx]:.4f}")

# Convert weights to distances (distance = 1/weight)
# Create distance matrix for shortest path
distances_matrix = W_directed.copy()
distances_matrix.data = 1.0 / (distances_matrix.data + 1e-10)

# Compute shortest paths from root
print(f"      Computing shortest paths from root...")
dist_from_root = shortest_path(distances_matrix, directed=True, indices=root_idx, return_predecessors=False)

# Handle unreachable nodes (inf distances)
n_unreachable = np.isinf(dist_from_root).sum()
if n_unreachable > 0:
    print(f"      WARNING: {n_unreachable} unreachable cells from root")
    # Set unreachable to max finite distance + 1
    max_finite = dist_from_root[np.isfinite(dist_from_root)].max()
    dist_from_root[np.isinf(dist_from_root)] = max_finite + 1

# Normalize to [0, 1]
PT_stress = (dist_from_root - dist_from_root.min()) / (dist_from_root.max() - dist_from_root.min())

print(f"\n  PT_stress statistics:")
print(f"    Mean: {PT_stress.mean():.4f}")
print(f"    Std: {PT_stress.std():.4f}")
print(f"    Range: [{PT_stress.min():.4f}, {PT_stress.max():.4f}]")

# Add to dataframe
df['PT_stress'] = PT_stress

# ============================================================================
# Save Results
# ============================================================================
print("\n[Saving results...]")

output_file = output_dir / 'cell_level_features_ALL_with_PT_dpt_stress.csv'
df.to_csv(output_file, index=False)
print(f"  Saved: {output_file}")

# ============================================================================
# QC Plots
# ============================================================================
print("\n[Creating QC plots...]")

fig, axes = plt.subplots(2, 3, figsize=(18, 12))

# Plot 1: PT_stress vs PT_imes
ax = axes[0, 0]
scatter = ax.scatter(df['PT_imes'], df['PT_stress'], alpha=0.1, s=1, c=df['S_stress_base'], cmap='viridis')
ax.set_xlabel('PT_imes (IMES-derived)', fontsize=12)
ax.set_ylabel('PT_stress (Stress-based)', fontsize=12)
corr = np.corrcoef(df['PT_imes'], df['PT_stress'])[0, 1]
ax.set_title(f'PT_imes vs PT_stress\nCorrelation: {corr:.3f}', fontsize=12)
plt.colorbar(scatter, ax=ax, label='S_stress_base')

# Plot 2: PT_stress vs PT_dpt
ax = axes[0, 1]
scatter = ax.scatter(df['PT_dpt'], df['PT_stress'], alpha=0.1, s=1, c=df['stress_total'], cmap='plasma')
ax.set_xlabel('PT_dpt (Diffusion-based)', fontsize=12)
ax.set_ylabel('PT_stress (Stress-based)', fontsize=12)
corr = np.corrcoef(df['PT_dpt'], df['PT_stress'])[0, 1]
ax.set_title(f'PT_dpt vs PT_stress\nCorrelation: {corr:.3f}', fontsize=12)
plt.colorbar(scatter, ax=ax, label='stress_total')

# Plot 3: PT_stress distribution
ax = axes[0, 2]
ax.hist(df['PT_stress'], bins=50, alpha=0.7, color='green', edgecolor='black')
ax.set_xlabel('PT_stress', fontsize=12)
ax.set_ylabel('Frequency', fontsize=12)
ax.set_title(f'PT_stress Distribution\nMean: {PT_stress.mean():.3f}, Std: {PT_stress.std():.3f}', fontsize=12)

# Plot 4: PT_stress by key cell types
ax = axes[1, 0]
key_celltypes = [
    'Ex.L5.VAT1L.EYA4',
    'Glia.Oligo',
    'Glia.Astro.GFAP-neg',
    'Ex.L2.L3.CUX2.RASGRF2',
    'Glia.Micro'
]
key_celltypes = [ct for ct in key_celltypes if ct in df['cell_type'].values]

data = [df[df['cell_type'] == ct]['PT_stress'].values for ct in key_celltypes]
bp = ax.boxplot(data, labels=[ct.split('.')[-1] if '.' in ct else ct for ct in key_celltypes], patch_artist=True)
for patch in bp['boxes']:
    patch.set_facecolor('lightgreen')
ax.set_ylabel('PT_stress', fontsize=12)
ax.set_title('PT_stress by Key Cell Types', fontsize=12)
ax.tick_params(axis='x', rotation=45)

# Plot 5: S_stress_base vs stress_total
ax = axes[1, 1]
scatter = ax.scatter(df['stress_total'], df['S_stress_base'], alpha=0.1, s=1, c=df['PT_stress'], cmap='coolwarm')
ax.set_xlabel('stress_total (existing)', fontsize=12)
ax.set_ylabel('S_stress_base (composite)', fontsize=12)
corr = np.corrcoef(df['stress_total'], df['S_stress_base'])[0, 1]
ax.set_title(f'Stress Measures\nCorrelation: {corr:.3f}', fontsize=12)
plt.colorbar(scatter, ax=ax, label='PT_stress')

# Plot 6: PT_stress vs S_stress_base
ax = axes[1, 2]
scatter = ax.scatter(df['S_stress_base'], df['PT_stress'], alpha=0.1, s=1, c=df['PT_imes'], cmap='viridis')
ax.set_xlabel('S_stress_base', fontsize=12)
ax.set_ylabel('PT_stress', fontsize=12)
corr = np.corrcoef(df['S_stress_base'], df['PT_stress'])[0, 1]
ax.set_title(f'Stress Potential vs PT_stress\nCorrelation: {corr:.3f}', fontsize=12)
plt.colorbar(scatter, ax=ax, label='PT_imes')

plt.tight_layout()
plot_file = output_dir / 'PT_stress_QC_plots.png'
plt.savefig(plot_file, dpi=150, bbox_inches='tight')
plt.close()
print(f"  Saved: {plot_file}")

print("\n" + "="*80)
print("Step 2 completed!")
print(f"PT_stress successfully constructed for {len(df):,} cells")
print(f"Correlation PT_imes vs PT_stress: {np.corrcoef(df['PT_imes'], df['PT_stress'])[0,1]:.3f}")
print(f"Correlation PT_dpt vs PT_stress: {np.corrcoef(df['PT_dpt'], df['PT_stress'])[0,1]:.3f}")
print("="*80)
