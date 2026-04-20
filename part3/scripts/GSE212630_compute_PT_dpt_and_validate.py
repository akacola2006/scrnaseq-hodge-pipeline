#!/usr/bin/env python3
"""
GSE_212630: Compute Diffusion Pseudotime (PT_dpt) and Validate vs TDP-43 Stages
================================================================================

This script:
1. Loads the combined metadata with stress scores
2. Computes PT_dpt using Diffusion Pseudotime algorithm
3. Validates whether PT_dpt correlates with TDP-43 disease stages
4. Generates validation figures

Author: Claude Code + User
Date: 2025-12-07
"""

import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.neighbors import NearestNeighbors
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import eigsh
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

print("=" * 80)
print("GSE_212630: Compute PT_dpt and Validate vs TDP-43 Stages")
print("=" * 80)

# ============================================================================
# CONFIGURATION
# ============================================================================

INPUT_DIR = Path(__file__).parent.parent / 'results' / 'GSE212630_ids_analysis'
OUTPUT_DIR = INPUT_DIR  # Same directory
INPUT_FILE = INPUT_DIR / 'combined_metadata_with_stress.csv'

# TDP-43 stage ordering (expected disease progression)
TDP_STAGE_ORDER = ['Control', 'TDPneg', 'TDPmed', 'TDPhigh']
TDP_STAGE_TO_NUMERIC = {'Control': 0, 'TDPneg': 1, 'TDPmed': 2, 'TDPhigh': 3}

# DPT parameters
K_NEIGHBORS = 30
N_COMPONENTS = 10

# ============================================================================
# STEP 1: LOAD DATA
# ============================================================================

print("\n[Step 1] Loading data...")
df = pd.read_csv(INPUT_FILE)
print(f"  Total cells: {len(df):,}")
print(f"  Columns: {list(df.columns)}")

# Check conditions
print(f"\n  Conditions: {df['condition'].unique().tolist()}")
print(f"  Cells per condition:")
for cond in TDP_STAGE_ORDER:
    n = (df['condition'] == cond).sum()
    print(f"    {cond}: {n:,}")

# ============================================================================
# STEP 2: PREPARE FEATURE SPACE
# ============================================================================

print("\n[Step 2] Preparing feature space...")

# Use stress components as features (similar to module scores in original project)
stress_cols = [col for col in df.columns if col.startswith('stress_') and col != 'stress_total']
print(f"  Stress components: {stress_cols}")

# Add stress_total
feature_cols = stress_cols + ['stress_total']
print(f"  Total features: {len(feature_cols)}")

# Check for missing values
missing = df[feature_cols].isnull().sum().sum()
if missing > 0:
    print(f"  WARNING: {missing} missing values found, filling with 0")
    df[feature_cols] = df[feature_cols].fillna(0)

# Extract and scale features
X = df[feature_cols].values
print(f"  Feature matrix shape: {X.shape}")

scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)
print(f"  Scaled features: mean={X_scaled.mean():.6f}, std={X_scaled.std():.6f}")

# ============================================================================
# STEP 3: BUILD KNN GRAPH AND AFFINITY MATRIX
# ============================================================================

print("\n[Step 3] Building kNN graph...")
print(f"  K = {K_NEIGHBORS} neighbors")

nbrs = NearestNeighbors(n_neighbors=K_NEIGHBORS + 1, metric='euclidean', n_jobs=-1)
nbrs.fit(X_scaled)
distances, indices = nbrs.kneighbors(X_scaled)

print(f"  kNN graph built: {len(df)} nodes")

# Build affinity matrix
print("\n  Building affinity matrix...")
sigma = np.median(distances[:, 1:])
print(f"  Sigma (median distance): {sigma:.4f}")

n_cells = len(df)
row_idx = []
col_idx = []
affinity_values = []

for i in range(n_cells):
    for j_idx, j in enumerate(indices[i, 1:]):
        dist = distances[i, j_idx + 1]
        affinity = np.exp(-(dist ** 2) / (2 * sigma ** 2))
        row_idx.append(i)
        col_idx.append(j)
        affinity_values.append(affinity)

W = csr_matrix((affinity_values, (row_idx, col_idx)), shape=(n_cells, n_cells))
W = (W + W.T) / 2  # Make symmetric
print(f"  Affinity matrix: {W.shape}, {W.nnz:,} non-zero entries")

# ============================================================================
# STEP 4: COMPUTE DIFFUSION OPERATOR AND EIGENDECOMPOSITION
# ============================================================================

print("\n[Step 4] Computing diffusion operator...")

# Degree matrix
D = np.array(W.sum(axis=1)).flatten()
D_inv_sqrt = 1.0 / np.sqrt(D + 1e-10)
D_inv_sqrt_sparse = csr_matrix((D_inv_sqrt, (range(n_cells), range(n_cells))))

# Normalized Laplacian
L = D_inv_sqrt_sparse @ W @ D_inv_sqrt_sparse
print(f"  Normalized Laplacian computed")

# Eigendecomposition
print(f"\n  Computing top {N_COMPONENTS} eigenvectors...")
eigenvalues, eigenvectors = eigsh(L, k=N_COMPONENTS, which='LM')

# Sort by eigenvalue (descending)
idx = eigenvalues.argsort()[::-1]
eigenvalues = eigenvalues[idx]
eigenvectors = eigenvectors[:, idx]

print(f"  Eigenvalues: {eigenvalues[:5]}...")
print(f"  Eigenvectors shape: {eigenvectors.shape}")

# ============================================================================
# STEP 5: COMPUTE DIFFUSION PSEUDOTIME (DPT)
# ============================================================================

print("\n[Step 5] Computing Diffusion Pseudotime (DPT)...")

# Select root cell: Control cell with lowest stress
control_mask = df['condition'] == 'Control'
control_indices = df[control_mask].index.tolist()

if len(control_indices) == 0:
    print("  ERROR: No Control cells found!")
    exit(1)

# Find Control cell with lowest stress_total
control_stress = df.loc[control_indices, 'stress_total']
root_idx = control_stress.idxmin()
root_cell = df.loc[root_idx]

print(f"  Root cell index: {root_idx}")
print(f"  Root cell_type: {root_cell['cell_type']}")
print(f"  Root condition: {root_cell['condition']}")
print(f"  Root stress_total: {root_cell['stress_total']:.4f}")

# DPT: Diffusion distance from root
# Use components 2-6 (skip first constant eigenvector)
dpt_components = eigenvectors[:, 1:6]

# Get root position in array (not DataFrame index)
root_array_idx = df.index.get_loc(root_idx)
root_vector = dpt_components[root_array_idx, :]

# Calculate Euclidean distance in diffusion space
dpt_distances = np.linalg.norm(dpt_components - root_vector, axis=1)

# Normalize to [0, 1]
PT_dpt = (dpt_distances - dpt_distances.min()) / (dpt_distances.max() - dpt_distances.min())

print(f"  PT_dpt range: [{PT_dpt.min():.4f}, {PT_dpt.max():.4f}]")
print(f"  PT_dpt mean: {PT_dpt.mean():.4f}")
print(f"  PT_dpt std: {PT_dpt.std():.4f}")

# Add to dataframe
df['PT_dpt'] = PT_dpt

# ============================================================================
# STEP 6: VALIDATE PT_dpt vs TDP-43 STAGES
# ============================================================================

print("\n[Step 6] Validating PT_dpt vs TDP-43 stages...")

# Convert TDP stages to numeric
df['TDP_stage_numeric'] = df['condition'].map(TDP_STAGE_TO_NUMERIC)

# Correlation test
corr, p_value = stats.spearmanr(df['TDP_stage_numeric'], df['PT_dpt'])
print(f"\n  Spearman correlation (TDP stage vs PT_dpt):")
print(f"    rho = {corr:.4f}")
print(f"    p-value = {p_value:.2e}")

# Mean PT_dpt by stage
print(f"\n  Mean PT_dpt by TDP-43 stage:")
stage_means = df.groupby('condition')['PT_dpt'].agg(['mean', 'std', 'count'])
for stage in TDP_STAGE_ORDER:
    if stage in stage_means.index:
        row = stage_means.loc[stage]
        print(f"    {stage:10s}: mean={row['mean']:.4f} ± {row['std']:.4f} (n={int(row['count'])})")

# ANOVA test
print(f"\n  Kruskal-Wallis test (non-parametric ANOVA):")
groups = [df[df['condition'] == stage]['PT_dpt'].values for stage in TDP_STAGE_ORDER]
groups = [g for g in groups if len(g) > 0]
h_stat, kw_p = stats.kruskal(*groups)
print(f"    H-statistic = {h_stat:.4f}")
print(f"    p-value = {kw_p:.2e}")

# Trend test (Jonckheere-Terpstra style using Spearman)
print(f"\n  Monotonic trend test:")
if corr > 0:
    print(f"    Direction: INCREASING (as expected for disease progression)")
else:
    print(f"    Direction: DECREASING (opposite to expectation)")

if p_value < 0.001:
    print(f"    Significance: HIGHLY SIGNIFICANT (p < 0.001)")
elif p_value < 0.05:
    print(f"    Significance: SIGNIFICANT (p < 0.05)")
else:
    print(f"    Significance: NOT SIGNIFICANT (p >= 0.05)")

# ============================================================================
# STEP 7: SAVE RESULTS
# ============================================================================

print("\n[Step 7] Saving results...")

output_file = OUTPUT_DIR / 'combined_metadata_with_PT_dpt.csv'
df.to_csv(output_file, index=False)
print(f"  Saved: {output_file}")

# ============================================================================
# STEP 8: GENERATE VALIDATION FIGURES
# ============================================================================

print("\n[Step 8] Generating validation figures...")

fig, axes = plt.subplots(2, 3, figsize=(18, 12))

# Plot 1: PT_dpt distribution by TDP stage (boxplot)
ax = axes[0, 0]
order = TDP_STAGE_ORDER
bp_data = [df[df['condition'] == stage]['PT_dpt'].values for stage in order]
bp = ax.boxplot(bp_data, labels=order, patch_artist=True)
colors = ['#2ecc71', '#f1c40f', '#e67e22', '#e74c3c']
for patch, color in zip(bp['boxes'], colors):
    patch.set_facecolor(color)
    patch.set_alpha(0.7)
ax.set_xlabel('TDP-43 Stage', fontsize=12)
ax.set_ylabel('PT_dpt', fontsize=12)
ax.set_title(f'PT_dpt by TDP-43 Stage\nSpearman ρ = {corr:.3f}, p = {p_value:.2e}', fontsize=12)
ax.grid(True, alpha=0.3)

# Plot 2: PT_dpt distribution (histogram by condition)
ax = axes[0, 1]
for stage, color in zip(order, colors):
    data = df[df['condition'] == stage]['PT_dpt']
    ax.hist(data, bins=30, alpha=0.5, label=stage, color=color, density=True)
ax.set_xlabel('PT_dpt', fontsize=12)
ax.set_ylabel('Density', fontsize=12)
ax.set_title('PT_dpt Distribution by Condition', fontsize=12)
ax.legend()

# Plot 3: Mean PT_dpt trend
ax = axes[0, 2]
means = [df[df['condition'] == stage]['PT_dpt'].mean() for stage in order]
stds = [df[df['condition'] == stage]['PT_dpt'].std() for stage in order]
x = range(len(order))
ax.errorbar(x, means, yerr=stds, fmt='o-', capsize=5, capthick=2, linewidth=2, markersize=10, color='black')
ax.set_xticks(x)
ax.set_xticklabels(order)
ax.set_xlabel('TDP-43 Stage', fontsize=12)
ax.set_ylabel('Mean PT_dpt ± SD', fontsize=12)
ax.set_title('PT_dpt Trend Across TDP-43 Stages', fontsize=12)
ax.grid(True, alpha=0.3)

# Plot 4: PT_dpt vs stress_total
ax = axes[1, 0]
scatter = ax.scatter(df['PT_dpt'], df['stress_total'], c=df['TDP_stage_numeric'],
                     cmap='RdYlGn_r', alpha=0.3, s=5)
ax.set_xlabel('PT_dpt', fontsize=12)
ax.set_ylabel('stress_total', fontsize=12)
stress_corr = np.corrcoef(df['PT_dpt'], df['stress_total'])[0, 1]
ax.set_title(f'PT_dpt vs Stress Total\nPearson r = {stress_corr:.3f}', fontsize=12)
cbar = plt.colorbar(scatter, ax=ax)
cbar.set_label('TDP Stage')
cbar.set_ticks([0, 1, 2, 3])
cbar.set_ticklabels(['Control', 'TDPneg', 'TDPmed', 'TDPhigh'])

# Plot 5: PT_dpt by cell type (top 10)
ax = axes[1, 1]
ct_means = df.groupby('cell_type')['PT_dpt'].mean().sort_values(ascending=False)
top_cts = ct_means.head(10).index.tolist()
data = [df[df['cell_type'] == ct]['PT_dpt'].values for ct in top_cts]
bp = ax.boxplot(data, labels=[ct.replace('.', '\n') for ct in top_cts], patch_artist=True)
for patch in bp['boxes']:
    patch.set_facecolor('lightblue')
ax.set_ylabel('PT_dpt', fontsize=12)
ax.set_title('PT_dpt by Cell Type (Top 10)', fontsize=12)
ax.tick_params(axis='x', rotation=45, labelsize=8)

# Plot 6: Eigenvalue spectrum
ax = axes[1, 2]
ax.plot(range(1, len(eigenvalues) + 1), eigenvalues, 'o-', color='black', markersize=8)
ax.set_xlabel('Component', fontsize=12)
ax.set_ylabel('Eigenvalue', fontsize=12)
ax.set_title('Diffusion Eigenvalue Spectrum', fontsize=12)
ax.grid(True, alpha=0.3)

plt.tight_layout()
fig_file = OUTPUT_DIR / 'Fig_PT_dpt_validation.png'
plt.savefig(fig_file, dpi=150, bbox_inches='tight')
plt.close()
print(f"  Saved: {fig_file}")

# ============================================================================
# SUMMARY
# ============================================================================

print("\n" + "=" * 80)
print("VALIDATION SUMMARY")
print("=" * 80)

print(f"""
PT_dpt was calculated using Diffusion Pseudotime algorithm:
  - Features used: {len(feature_cols)} stress components
  - Root cell: Control cell with lowest stress
  - kNN neighbors: {K_NEIGHBORS}
  - Diffusion components: {N_COMPONENTS}

VALIDATION RESULTS:
  - Spearman correlation (TDP stage vs PT_dpt): ρ = {corr:.4f}, p = {p_value:.2e}
  - Kruskal-Wallis test: H = {h_stat:.4f}, p = {kw_p:.2e}
  - PT_dpt-Stress correlation: r = {stress_corr:.4f}

INTERPRETATION:
""")

if corr > 0 and p_value < 0.05:
    print("  ✓ PT_dpt INCREASES with TDP-43 disease stage (as expected)")
    print("  ✓ The computed pseudotime is CONSISTENT with labeled disease stages")
    print("  ✓ This validates that the Diffusion Pseudotime captures disease progression")
elif corr < 0 and p_value < 0.05:
    print("  ✗ PT_dpt DECREASES with TDP-43 stage (opposite to expectation)")
    print("  → Consider reversing PT_dpt or checking root cell selection")
else:
    print("  ? No significant correlation between PT_dpt and TDP-43 stages")
    print("  → The pseudotime may not align with disease staging")

print(f"""
FILES GENERATED:
  - combined_metadata_with_PT_dpt.csv
  - Fig_PT_dpt_validation.png
""")

print("=" * 80)
