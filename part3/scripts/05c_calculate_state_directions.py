#!/usr/bin/env python3
"""
Phase 5c: Calculate State-Level Direction Scores
=================================================

Calculates directional scores between cell states using:
1. kNN graph transitions: Count edges from state A → state B
2. PT ordering: If PT(B) > PT(A), stronger evidence for A → B
3. Upstream score: Combined metric

This identifies which states are "upstream" (sources) vs "downstream" (sinks)
in the cell-state causal hierarchy.

Method:
- Build kNN graph (k=15) in module feature space
- For each edge (cell_i, cell_j):
  - Record (state_i, state_j) transition
  - Weight by PT difference: if PT_j > PT_i, evidence for i→j
- Aggregate to state-level directional matrix

Output:
- results/cell_state_causality/state_direction_matrix.csv
- results/cell_state_causality/state_upstream_scores.csv
"""

import pandas as pd
import numpy as np
from pathlib import Path
from sklearn.neighbors import NearestNeighbors
from collections import defaultdict
import warnings
warnings.filterwarnings('ignore')

# Paths
ids_dir = Path('/home/akaco/als/motor_cortex_analysis/ids_causal_analysis')
cell_state_dir = ids_dir / 'results' / 'cell_state_causality'

print("="*80)
print("Phase 5c: Calculate State-Level Direction Scores")
print("="*80)

# ============================================================================
# 1. Load cell-level data with state labels
# ============================================================================
print("\n[1] Loading cell states...")

df = pd.read_csv(cell_state_dir / 'cell_states_labeled.csv')

print(f"  Loaded: {len(df)} cells")
print(f"  States: {df['cell_state'].nunique()}")

# ============================================================================
# 2. Build kNN graph
# ============================================================================
print("\n[2] Building kNN graph...")

# Extract features for kNN
module_cols = [c for c in df.columns if c.startswith('module_')]
feature_cols = module_cols + ['PT_IDS', 'stress_total']

X = df[feature_cols].values
print(f"  Feature matrix: {X.shape}")

# Build kNN graph (k=15)
K = 15
print(f"  Building kNN graph with k={K}...")

nbrs = NearestNeighbors(n_neighbors=K+1, metric='euclidean')  # k+1 because includes self
nbrs.fit(X)
distances, indices = nbrs.kneighbors(X)

print(f"  kNN graph built: {len(df)} nodes, ~{len(df)*K} edges")

# ============================================================================
# 3. Calculate state-level transitions
# ============================================================================
print("\n[3] Calculating state-level transitions...")

# Direction matrix: direction_counts[state_i][state_j] = count of i→j transitions
direction_counts = defaultdict(lambda: defaultdict(int))
pt_weighted_scores = defaultdict(lambda: defaultdict(float))

# Process each cell's neighbors
for i in range(len(df)):
    cell_i = df.iloc[i]
    state_i = cell_i['cell_state']
    pt_i = cell_i['PT_IDS']

    # Iterate over neighbors (skip first neighbor which is self)
    for j_idx in indices[i, 1:]:
        cell_j = df.iloc[j_idx]
        state_j = cell_j['cell_state']
        pt_j = cell_j['PT_IDS']

        # Skip self-loops (same state)
        if state_i == state_j:
            continue

        # Direction evidence from PT
        if pt_j > pt_i:
            # j is downstream of i (i → j)
            direction_counts[state_i][state_j] += 1
            pt_diff = pt_j - pt_i
            pt_weighted_scores[state_i][state_j] += pt_diff
        elif pt_i > pt_j:
            # i is downstream of j (j → i)
            direction_counts[state_j][state_i] += 1
            pt_diff = pt_i - pt_j
            pt_weighted_scores[state_j][state_i] += pt_diff
        else:
            # Equal PT: count both directions equally (weak evidence)
            direction_counts[state_i][state_j] += 0.5
            direction_counts[state_j][state_i] += 0.5

# ============================================================================
# 4. Create direction matrix
# ============================================================================
print("\n[4] Creating direction matrix...")

all_states = sorted(df['cell_state'].unique())
n_states = len(all_states)

# Initialize matrix
direction_matrix = np.zeros((n_states, n_states))
pt_weight_matrix = np.zeros((n_states, n_states))

state_to_idx = {state: i for i, state in enumerate(all_states)}

for state_i, targets in direction_counts.items():
    i_idx = state_to_idx[state_i]
    for state_j, count in targets.items():
        j_idx = state_to_idx[state_j]
        direction_matrix[i_idx, j_idx] = count
        pt_weight_matrix[i_idx, j_idx] = pt_weighted_scores[state_i][state_j]

# Convert to DataFrame
direction_df = pd.DataFrame(direction_matrix, index=all_states, columns=all_states)
pt_weight_df = pd.DataFrame(pt_weight_matrix, index=all_states, columns=all_states)

# Save
direction_file = cell_state_dir / 'state_direction_matrix.csv'
direction_df.to_csv(direction_file)
print(f"  Saved: {direction_file}")

pt_weight_file = cell_state_dir / 'state_pt_weight_matrix.csv'
pt_weight_df.to_csv(pt_weight_file)
print(f"  Saved: {pt_weight_file}")

print(f"  Direction matrix: {direction_df.shape}")
print(f"  Non-zero edges: {(direction_matrix > 0).sum()}")

# ============================================================================
# 5. Calculate upstream scores
# ============================================================================
print("\n[5] Calculating upstream scores...")

upstream_scores = []

for state in all_states:
    idx = state_to_idx[state]

    # Outgoing: state → others
    out_count = direction_matrix[idx, :].sum()
    out_weight = pt_weight_matrix[idx, :].sum()

    # Incoming: others → state
    in_count = direction_matrix[:, idx].sum()
    in_weight = pt_weight_matrix[:, idx].sum()

    # Net scores (positive = more upstream, negative = more downstream)
    net_count = out_count - in_count
    net_weight = out_weight - in_weight

    # Upstream score: weighted combination
    upstream_score = 0.5 * net_count + 0.5 * net_weight

    upstream_scores.append({
        'cell_state': state,
        'out_count': out_count,
        'in_count': in_count,
        'net_count': net_count,
        'out_weight': out_weight,
        'in_weight': in_weight,
        'net_weight': net_weight,
        'upstream_score': upstream_score
    })

upstream_df = pd.DataFrame(upstream_scores)
upstream_df = upstream_df.sort_values('upstream_score', ascending=False)

# Save
upstream_file = cell_state_dir / 'state_upstream_scores.csv'
upstream_df.to_csv(upstream_file, index=False)
print(f"  Saved: {upstream_file}")

# ============================================================================
# 6. Display top upstream and downstream states
# ============================================================================
print("\n[6] Top upstream states (global sources):")
print(upstream_df.head(10)[['cell_state', 'upstream_score', 'net_count', 'net_weight']].to_string(index=False))

print("\n[6] Top downstream states (global sinks):")
print(upstream_df.tail(10)[['cell_state', 'upstream_score', 'net_count', 'net_weight']].to_string(index=False))

print("\n" + "="*80)
print("Phase 5c completed!")
print(f"Direction matrix: {n_states} × {n_states} states")
print(f"Non-zero transitions: {(direction_matrix > 0).sum()}")
print("="*80)
