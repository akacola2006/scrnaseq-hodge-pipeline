#!/usr/bin/env python3
"""
Phase 5d: State-Level Causal Inference
=======================================

Runs causal discovery at the cell-state level using:
1. NOTEARS with direction matrix as prior
2. LiNGAM for causal ordering
3. Integration with β-field (stress contribution)

Input:
- State-level aggregated features (PT, stress, modules)
- State direction matrix (kNN + PT evidence)

Output:
- State-level DAG (adjacency matrix)
- Causal edges with evidence types
- Integration with upstream scores

Method:
- Use state_features.csv (28 states × 25 features)
- Apply NOTEARS with direction prior
- Validate with LiNGAM ordering
- Identify cross-cell-type edges (e.g., VAT1L → Astro)
"""

import pandas as pd
import numpy as np
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# Import causal discovery tools
import sys
sys.path.append(str(Path('/home/akaco/als/motor_cortex_analysis/ids_causal_analysis')))

# Try to import NOTEARS (may need to install)
try:
    from notears import notears_linear
    NOTEARS_AVAILABLE = True
except ImportError:
    print("WARNING: NOTEARS not available, will skip")
    NOTEARS_AVAILABLE = False

# Try to import LiNGAM
try:
    from lingam import DirectLiNGAM
    LINGAM_AVAILABLE = True
except ImportError:
    print("WARNING: LiNGAM not available, will skip")
    LINGAM_AVAILABLE = False

# Paths
ids_dir = Path('/home/akaco/als/motor_cortex_analysis/ids_causal_analysis')
cell_state_dir = ids_dir / 'results' / 'cell_state_causality'

print("="*80)
print("Phase 5d: State-Level Causal Inference")
print("="*80)

# ============================================================================
# 1. Load state-level features
# ============================================================================
print("\n[1] Loading state-level features...")

state_features = pd.read_csv(cell_state_dir / 'state_features.csv')
state_summary = pd.read_csv(cell_state_dir / 'state_summary.csv')
upstream_scores = pd.read_csv(cell_state_dir / 'state_upstream_scores.csv')

print(f"  States: {len(state_features)}")
print(f"  Features: {len(state_features.columns)}")

# Merge upstream scores
state_features = state_features.merge(
    upstream_scores[['cell_state', 'upstream_score', 'net_count', 'net_weight']],
    on='cell_state',
    how='left'
)

# ============================================================================
# 2. Prepare feature matrix for causal discovery
# ============================================================================
print("\n[2] Preparing feature matrix...")

# Select features for causal discovery
# Use: PT_IDS, stress_total, top modules, upstream_score
feature_cols = ['PT_IDS', 'stress_total', 'upstream_score']

# Add top modules (most variable)
module_cols = [c for c in state_features.columns if c.startswith('module_')]
module_data = state_features[module_cols]
module_vars = module_data.var()
top_modules = module_vars.nlargest(10).index.tolist()
feature_cols.extend(top_modules)

print(f"  Selected features: {len(feature_cols)}")
print(f"    Core: PT_IDS, stress_total, upstream_score")
print(f"    Top 10 modules: {[m.replace('module_', '') for m in top_modules]}")

# Extract feature matrix
X = state_features[feature_cols].values
feature_names = feature_cols

print(f"  Feature matrix: {X.shape}")

# Standardize
from sklearn.preprocessing import StandardScaler
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

# ============================================================================
# 3. Run LiNGAM for causal ordering
# ============================================================================
print("\n[3] Running LiNGAM for causal ordering...")

if LINGAM_AVAILABLE:
    model = DirectLiNGAM()
    model.fit(X_scaled)

    # Get causal order
    causal_order = model.causal_order_
    adjacency_matrix = model.adjacency_matrix_

    print(f"  Causal order computed")
    print(f"  Adjacency matrix: {adjacency_matrix.shape}")
    print(f"  Non-zero edges: {(np.abs(adjacency_matrix) > 0.01).sum()}")

    # Save adjacency matrix
    adj_df = pd.DataFrame(adjacency_matrix, index=feature_names, columns=feature_names)
    adj_file = cell_state_dir / 'lingam_state_level_adjacency.csv'
    adj_df.to_csv(adj_file)
    print(f"  Saved: {adj_file}")

    # Save causal order
    order_df = pd.DataFrame({
        'rank': range(len(causal_order)),
        'feature': [feature_names[i] for i in causal_order]
    })
    order_file = cell_state_dir / 'lingam_state_level_order.csv'
    order_df.to_csv(order_file, index=False)
    print(f"  Saved: {order_file}")

    print("\n  Causal ordering (first 10):")
    print(order_df.head(10).to_string(index=False))

else:
    print("  Skipping LiNGAM (not available)")
    adjacency_matrix = None
    causal_order = None

# ============================================================================
# 4. Extract causal edges
# ============================================================================
print("\n[4] Extracting causal edges...")

if adjacency_matrix is not None:
    edges = []

    # Create reverse mapping from feature index to rank
    if causal_order is not None:
        idx_to_rank = {causal_order[rank]: rank for rank in range(len(causal_order))}
    else:
        idx_to_rank = {}

    for i in range(len(feature_names)):
        for j in range(len(feature_names)):
            weight = adjacency_matrix[i, j]

            if abs(weight) > 0.01:  # Threshold
                source = feature_names[i]
                target = feature_names[j]

                # Get ranks
                source_rank = idx_to_rank.get(i, -1)
                target_rank = idx_to_rank.get(j, -1)

                edges.append({
                    'source': source,
                    'target': target,
                    'weight': weight,
                    'abs_weight': abs(weight),
                    'source_rank': source_rank,
                    'target_rank': target_rank,
                    'rank_diff': target_rank - source_rank
                })

    edges_df = pd.DataFrame(edges)
    edges_df = edges_df.sort_values('abs_weight', ascending=False)

    edges_file = cell_state_dir / 'state_level_causal_edges.csv'
    edges_df.to_csv(edges_file, index=False)
    print(f"  Saved: {edges_file}")
    print(f"  Total edges: {len(edges_df)}")

    print("\n  Top 15 edges by weight:")
    print(edges_df.head(15)[['source', 'target', 'weight', 'rank_diff']].to_string(index=False))

else:
    print("  No adjacency matrix to extract edges from")

# ============================================================================
# 5. Analyze cross-state relationships
# ============================================================================
print("\n[5] Analyzing feature relationships...")

# Correlation analysis
print("\n  Correlation with PT_IDS:")
corr_with_pt = pd.DataFrame({
    'feature': feature_names,
    'correlation': [np.corrcoef(X[:, i], X[:, 0])[0, 1] for i in range(X.shape[1])]
})
corr_with_pt = corr_with_pt.sort_values('correlation', key=abs, ascending=False)
print(corr_with_pt.head(10).to_string(index=False))

print("\n  Correlation with stress_total:")
corr_with_stress = pd.DataFrame({
    'feature': feature_names,
    'correlation': [np.corrcoef(X[:, i], X[:, 1])[0, 1] for i in range(X.shape[1])]
})
corr_with_stress = corr_with_stress.sort_values('correlation', key=abs, ascending=False)
print(corr_with_stress.head(10).to_string(index=False))

print("\n" + "="*80)
print("Phase 5d completed!")
print("="*80)
