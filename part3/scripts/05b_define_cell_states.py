#!/usr/bin/env python3
"""
Phase 5b: Define Cell States via Clustering
============================================

Defines "cell states" as (cell_type, cluster) pairs by:
1. For each cell type, cluster cells based on module features
2. Create state labels: "{cell_type}_State{i}"
3. Calculate state-level aggregated features

Clustering approach:
- Use module scores + PT_IDS as features
- Apply k-means clustering (k=3-5 per cell type)
- Validate clusters by PT distribution

Output:
- results/cell_state_causality/cell_states_labeled.csv
- results/cell_state_causality/state_features.csv
- results/cell_state_causality/state_summary.csv
"""

import pandas as pd
import numpy as np
from pathlib import Path
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler
import warnings
warnings.filterwarnings('ignore')

# Paths
ids_dir = Path('/home/akaco/als/motor_cortex_analysis/ids_causal_analysis')
cell_state_dir = ids_dir / 'results' / 'cell_state_causality'

print("="*80)
print("Phase 5b: Define Cell States via Clustering")
print("="*80)

# ============================================================================
# 1. Load cell-level features
# ============================================================================
print("\n[1] Loading cell-level features...")

features_file = cell_state_dir / 'cell_level_features_ALL.csv'
df = pd.read_csv(features_file)

print(f"  Loaded: {len(df)} cells, {len(df.columns)} features")
print(f"  Cell types: {df['cell_type'].nunique()}")

# ============================================================================
# 2. Cluster each cell type
# ============================================================================
print("\n[2] Clustering each cell type...")

# Clustering parameters - dynamic based on cell count
def get_n_clusters(n_cells):
    """Determine number of clusters based on cell count"""
    if n_cells < 100:
        return 2  # Very small: 2 states
    elif n_cells < 500:
        return 3  # Small: 3 states
    elif n_cells < 2000:
        return 4  # Medium: 4 states
    elif n_cells < 5000:
        return 5  # Large: 5 states
    else:
        return 6  # Very large: 6 states

all_state_labels = []
state_summary = []

for celltype in df['cell_type'].unique():
    print(f"\n  Clustering: {celltype}")

    # Filter cells of this type
    ct_df = df[df['cell_type'] == celltype].copy()
    n_cells = len(ct_df)

    # Get number of clusters based on cell count
    n_clusters = get_n_clusters(n_cells)
    print(f"    Cells: {n_cells}, Clusters: {n_clusters}")

    # Extract features for clustering
    module_cols = [c for c in ct_df.columns if c.startswith('module_')]
    feature_cols = module_cols + ['PT_IDS', 'stress_total']

    X = ct_df[feature_cols].values

    # Standardize features
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    # K-means clustering
    kmeans = KMeans(n_clusters=n_clusters, random_state=42, n_init=20)
    cluster_labels = kmeans.fit_predict(X_scaled)

    # Assign state labels
    ct_df['cluster'] = cluster_labels
    ct_df['cell_state'] = ct_df['cluster'].apply(lambda x: f"{celltype}_State{x}")

    all_state_labels.append(ct_df)

    # Analyze clusters
    for i in range(n_clusters):
        cluster_cells = ct_df[ct_df['cluster'] == i]
        n_cluster = len(cluster_cells)
        n_als = (cluster_cells['condition'] == 'ALS').sum()
        n_ctrl = (cluster_cells['condition'].isin(['Control', 'CTRL'])).sum()
        mean_pt = cluster_cells['PT_IDS'].mean()
        mean_stress = cluster_cells['stress_total'].mean()

        state_summary.append({
            'cell_type': celltype,
            'cluster': i,
            'cell_state': f"{celltype}_State{i}",
            'n_cells': n_cluster,
            'n_ALS': n_als,
            'n_Control': n_ctrl,
            'pct_ALS': 100 * n_als / n_cluster if n_cluster > 0 else 0,
            'mean_PT': mean_pt,
            'mean_stress': mean_stress
        })

        print(f"      State{i}: {n_cluster:5d} cells, PT={mean_pt:.3f}, stress={mean_stress:.3f}, ALS={100*n_als/n_cluster:.1f}%")

# ============================================================================
# 3. Combine and save labeled cells
# ============================================================================
print("\n[3] Saving cell states...")

df_labeled = pd.concat(all_state_labels, ignore_index=True)

output_file = cell_state_dir / 'cell_states_labeled.csv'
df_labeled.to_csv(output_file, index=False)
print(f"  Saved: {output_file}")
print(f"  Total cells: {len(df_labeled)}")
print(f"  Total states: {df_labeled['cell_state'].nunique()}")

# ============================================================================
# 4. Calculate state-level features
# ============================================================================
print("\n[4] Calculating state-level aggregated features...")

# Aggregate by state
feature_cols = ['PT_IDS', 'stress_total'] + [c for c in df_labeled.columns if c.startswith('module_')]

state_features = df_labeled.groupby('cell_state')[feature_cols].mean()
state_features['n_cells'] = df_labeled.groupby('cell_state').size()
state_features['n_ALS'] = df_labeled.groupby('cell_state').apply(
    lambda x: (x['condition'] == 'ALS').sum()
)
state_features['n_Control'] = df_labeled.groupby('cell_state').apply(
    lambda x: (x['condition'].isin(['Control', 'CTRL'])).sum()
)
state_features['pct_ALS'] = 100 * state_features['n_ALS'] / state_features['n_cells']

# Add cell_type and cluster info
state_features['cell_type'] = state_features.index.str.rsplit('_State', n=1).str[0]
state_features['cluster'] = state_features.index.str.rsplit('_State', n=1).str[1].astype(int)

state_features.reset_index(inplace=True)

# Save state features
state_features_file = cell_state_dir / 'state_features.csv'
state_features.to_csv(state_features_file, index=False)
print(f"  Saved: {state_features_file}")
print(f"  States: {len(state_features)}")
print(f"  Features per state: {len(feature_cols)}")

# ============================================================================
# 5. Save state summary
# ============================================================================
summary_df = pd.DataFrame(state_summary)
summary_file = cell_state_dir / 'state_summary.csv'
summary_df.to_csv(summary_file, index=False)
print(f"  Saved: {summary_file}")

print("\n  State summary (first 10):")
print(summary_df[['cell_state', 'n_cells', 'mean_PT', 'mean_stress', 'pct_ALS']].head(10).to_string(index=False))

print("\n" + "="*80)
print("Phase 5b completed!")
print(f"Total states defined: {df_labeled['cell_state'].nunique()}")
print("="*80)
