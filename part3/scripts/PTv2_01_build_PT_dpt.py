#!/usr/bin/env python3
"""
Phase PTv2 - Step 0 & 1: Diffusion Pseudotime (PT_dpt) Construction
====================================================================

Purpose:
- Build a NEW pseudotime (PT_dpt) that is INDEPENDENT of the original IMES-derived PT
- Use only modules + stress features (NO PT_imes!)
- Use Diffusion Pseudotime (DPT) algorithm

Input:
- cell_level_features_ALL.csv (111,837 cells)

Output:
- cell_level_features_ALL_with_PT_dpt.csv
- PT_dpt_QC_plots.png
"""

import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.neighbors import NearestNeighbors
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import eigsh
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

print("="*80)
print("Phase PTv2 - Step 1: Building PT_dpt (Diffusion Pseudotime)")
print("="*80)

# Paths
base_dir = Path('/home/akaco/als/motor_cortex_analysis')
ids_dir = base_dir / 'ids_causal_analysis'
input_file = ids_dir / 'results' / 'cell_state_causality' / 'cell_level_features_ALL.csv'
output_dir = ids_dir / 'results' / 'PTv2_robustness'
output_dir.mkdir(exist_ok=True, parents=True)

# ============================================================================
# Step 0: Load Data and Prepare PT_imes Baseline
# ============================================================================
print("\n[Step 0] Loading data...")

df = pd.read_csv(input_file)
print(f"  Total cells: {len(df):,}")
print(f"  Columns: {len(df.columns)}")

# Rename PT_IDS to PT_imes for clarity
df = df.rename(columns={'PT_IDS': 'PT_imes'})
print(f"  Renamed PT_IDS → PT_imes (IMES-derived baseline)")

# Extract patient_id
df['patient_id'] = df['cell_id'].str.split('-').str[1]
print(f"  Unique patients: {df['patient_id'].nunique()}")

# ============================================================================
# Step 1: Build PT_dpt - Diffusion Pseudotime
# ============================================================================
print("\n[Step 1] Building PT_dpt...")

# 1a. Define feature space (NO PT_imes!)
print("\n  [1a] Defining feature space...")

module_cols = [col for col in df.columns if col.startswith('module_')]
print(f"      Available modules: {len(module_cols)}")

# Select top modules for DPT
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
    'module_Metabolism',
    'module_Autophagy',
    'module_DNA_Repair',
    'module_Myelination'
]

# Filter to available modules
top_modules = [m for m in top_modules if m in module_cols]
print(f"      Using {len(top_modules)} modules for PT_dpt")

# Add stress_total
feature_cols = top_modules + ['stress_total']
print(f"      Total features: {len(feature_cols)}")

# Check for missing values
missing = df[feature_cols].isnull().sum().sum()
if missing > 0:
    print(f"      WARNING: {missing} missing values found, filling with 0")
    df[feature_cols] = df[feature_cols].fillna(0)

# Extract features
X = df[feature_cols].values
print(f"      Feature matrix shape: {X.shape}")

# Standardize
print("\n  [1b] Standardizing features...")
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)
print(f"      Scaled features: mean={X_scaled.mean():.6f}, std={X_scaled.std():.6f}")

# 1c. Build kNN graph
print("\n  [1c] Building kNN graph...")
K = 30  # Number of neighbors
print(f"      K = {K} neighbors")

nbrs = NearestNeighbors(n_neighbors=K+1, metric='euclidean', n_jobs=-1)
nbrs.fit(X_scaled)
distances, indices = nbrs.kneighbors(X_scaled)

print(f"      kNN graph built: {len(df)} nodes, ~{len(df)*K} edges")

# 1d. Build affinity matrix (Gaussian kernel)
print("\n  [1d] Building affinity matrix...")

# Estimate sigma from median distance
sigma = np.median(distances[:, 1:])  # Exclude self (distance=0)
print(f"      Sigma (median distance): {sigma:.4f}")

# Build sparse affinity matrix
n_cells = len(df)
row_idx = []
col_idx = []
affinity_values = []

for i in range(n_cells):
    for j_idx, j in enumerate(indices[i, 1:]):  # Skip self
        dist = distances[i, j_idx + 1]
        affinity = np.exp(- (dist**2) / (2 * sigma**2))
        row_idx.append(i)
        col_idx.append(j)
        affinity_values.append(affinity)

# Create sparse matrix
W = csr_matrix((affinity_values, (row_idx, col_idx)), shape=(n_cells, n_cells))
print(f"      Affinity matrix: {W.shape}, {W.nnz:,} non-zero entries")

# Make symmetric
W = (W + W.T) / 2

# 1e. Compute diffusion operator
print("\n  [1e] Computing diffusion operator...")

# Degree matrix
D = np.array(W.sum(axis=1)).flatten()
D_inv_sqrt = 1.0 / np.sqrt(D + 1e-10)
D_inv_sqrt_sparse = csr_matrix((D_inv_sqrt, (range(n_cells), range(n_cells))))

# Normalized Laplacian: L = D^(-1/2) W D^(-1/2)
L = D_inv_sqrt_sparse @ W @ D_inv_sqrt_sparse

print(f"      Normalized Laplacian computed")

# 1f. Compute diffusion components via eigendecomposition
print("\n  [1f] Computing diffusion components...")

n_components = 10
print(f"      Computing top {n_components} eigenvectors...")

# Use eigsh for symmetric sparse matrix
eigenvalues, eigenvectors = eigsh(L, k=n_components, which='LM')

print(f"      Eigenvalues: {eigenvalues}")
print(f"      Eigenvectors shape: {eigenvectors.shape}")

# Sort by eigenvalue (descending)
idx = eigenvalues.argsort()[::-1]
eigenvalues = eigenvalues[idx]
eigenvectors = eigenvectors[:, idx]

# Store diffusion components
diffusion_components = eigenvectors
print(f"      Diffusion components: {diffusion_components.shape}")

# 1g. Compute Diffusion Pseudotime (DPT)
print("\n  [1g] Computing Diffusion Pseudotime (DPT)...")

# Select root cell: lowest stress_total
root_idx = df['stress_total'].argmin()
root_cell = df.iloc[root_idx]
print(f"      Root cell: {root_cell['cell_id']}")
print(f"      Root cell_type: {root_cell['cell_type']}")
print(f"      Root stress_total: {root_cell['stress_total']:.4f}")

# DPT: Diffusion distance from root
# Using first few diffusion components (exclude first constant eigenvector)
dpt_components = diffusion_components[:, 1:6]  # Use components 2-6

# Calculate Euclidean distance in diffusion space from root
root_vector = dpt_components[root_idx, :]
dpt_distances = np.linalg.norm(dpt_components - root_vector, axis=1)

# Normalize to [0, 1]
PT_dpt = (dpt_distances - dpt_distances.min()) / (dpt_distances.max() - dpt_distances.min())

print(f"      PT_dpt range: [{PT_dpt.min():.4f}, {PT_dpt.max():.4f}]")
print(f"      PT_dpt mean: {PT_dpt.mean():.4f}")
print(f"      PT_dpt std: {PT_dpt.std():.4f}")

# Add to dataframe
df['PT_dpt'] = PT_dpt

# ============================================================================
# Save Results
# ============================================================================
print("\n[Saving results...]")

output_file = output_dir / 'cell_level_features_ALL_with_PT_dpt.csv'
df.to_csv(output_file, index=False)
print(f"  Saved: {output_file}")

# ============================================================================
# QC Plots
# ============================================================================
print("\n[Creating QC plots...]")

fig, axes = plt.subplots(2, 3, figsize=(18, 12))

# Plot 1: PT_imes vs PT_dpt correlation
ax = axes[0, 0]
scatter = ax.scatter(df['PT_imes'], df['PT_dpt'], alpha=0.1, s=1, c=df['stress_total'], cmap='viridis')
ax.set_xlabel('PT_imes (IMES-derived)', fontsize=12)
ax.set_ylabel('PT_dpt (Diffusion-based)', fontsize=12)
ax.set_title(f'PT_imes vs PT_dpt\nCorrelation: {np.corrcoef(df["PT_imes"], df["PT_dpt"])[0,1]:.3f}', fontsize=12)
plt.colorbar(scatter, ax=ax, label='stress_total')

# Plot 2: PT_dpt distribution
ax = axes[0, 1]
ax.hist(df['PT_dpt'], bins=50, alpha=0.7, color='blue', edgecolor='black')
ax.set_xlabel('PT_dpt', fontsize=12)
ax.set_ylabel('Frequency', fontsize=12)
ax.set_title(f'PT_dpt Distribution\nMean: {PT_dpt.mean():.3f}, Std: {PT_dpt.std():.3f}', fontsize=12)

# Plot 3: PT_dpt by condition
ax = axes[0, 2]
als_dpt = df[df['condition'] == 'ALS']['PT_dpt']
ctrl_dpt = df[df['condition'] == 'Control']['PT_dpt']
ax.hist([als_dpt, ctrl_dpt], bins=50, alpha=0.7, label=['ALS', 'Control'], color=['red', 'blue'], edgecolor='black')
ax.set_xlabel('PT_dpt', fontsize=12)
ax.set_ylabel('Frequency', fontsize=12)
ax.set_title('PT_dpt by Condition', fontsize=12)
ax.legend()

# Plot 4: PT_dpt by key cell types
ax = axes[1, 0]
key_celltypes = [
    'Ex.L5.VAT1L.EYA4',
    'Glia.Oligo',
    'Glia.Astro.GFAP-neg',
    'Ex.L2.L3.CUX2.RASGRF2',
    'Glia.Micro'
]
key_celltypes = [ct for ct in key_celltypes if ct in df['cell_type'].values]

data = [df[df['cell_type'] == ct]['PT_dpt'].values for ct in key_celltypes]
bp = ax.boxplot(data, labels=[ct.split('.')[-1] if '.' in ct else ct for ct in key_celltypes], patch_artist=True)
for patch in bp['boxes']:
    patch.set_facecolor('lightblue')
ax.set_ylabel('PT_dpt', fontsize=12)
ax.set_title('PT_dpt by Key Cell Types', fontsize=12)
ax.tick_params(axis='x', rotation=45)

# Plot 5: Stress vs PT_dpt
ax = axes[1, 1]
scatter = ax.scatter(df['PT_dpt'], df['stress_total'], alpha=0.1, s=1, c=df['PT_imes'], cmap='plasma')
ax.set_xlabel('PT_dpt', fontsize=12)
ax.set_ylabel('stress_total', fontsize=12)
ax.set_title(f'PT_dpt vs Stress\nCorrelation: {np.corrcoef(df["PT_dpt"], df["stress_total"])[0,1]:.3f}', fontsize=12)
plt.colorbar(scatter, ax=ax, label='PT_imes')

# Plot 6: Eigenvalue spectrum
ax = axes[1, 2]
ax.plot(range(1, len(eigenvalues)+1), eigenvalues, 'o-', color='black')
ax.set_xlabel('Component', fontsize=12)
ax.set_ylabel('Eigenvalue', fontsize=12)
ax.set_title('Diffusion Eigenvalue Spectrum', fontsize=12)
ax.grid(True, alpha=0.3)

plt.tight_layout()
plot_file = output_dir / 'PT_dpt_QC_plots.png'
plt.savefig(plot_file, dpi=150, bbox_inches='tight')
plt.close()
print(f"  Saved: {plot_file}")

print("\n" + "="*80)
print("Step 1 completed!")
print(f"PT_dpt successfully constructed for {len(df):,} cells")
print(f"Correlation PT_imes vs PT_dpt: {np.corrcoef(df['PT_imes'], df['PT_dpt'])[0,1]:.3f}")
print("="*80)
