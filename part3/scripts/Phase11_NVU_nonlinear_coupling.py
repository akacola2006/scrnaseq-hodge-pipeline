#!/usr/bin/env python3
"""
Phase 11: NVU 4-node Nonlinear Coupling Analysis

Alternative to Phase 10 DAG approach.
Uses nonlinear regression, partial correlations, and severity adjustment
to identify coupling patterns among:
  - Vasc_driver
  - Glia_ER
  - Upper_hyper
  - VAT1L_Ca

Author: Claude Code
Date: 2025-11-24
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import spearmanr, pearsonr
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import PolynomialFeatures, StandardScaler
from sklearn.decomposition import PCA
from sklearn.metrics import r2_score
import networkx as nx
from itertools import combinations
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# ============================================================
# Configuration
# ============================================================

# Input files
NVU_DATA_FILE = 'results/phase10_nvu_dag/nvu_4node_data.csv'
PATHOLOGY_FILE = 'results/phase9p_vascular_driver/patient_pathology_signatures.csv'

# Output directory
OUTPUT_DIR = 'results/phase11_nvu_coupling'
PLOT_DIR = f'{OUTPUT_DIR}/plots'
Path(OUTPUT_DIR).mkdir(parents=True, exist_ok=True)
Path(PLOT_DIR).mkdir(parents=True, exist_ok=True)

# Variables
NVU_VARS = ['Vasc_driver', 'Glia_ER', 'Upper_hyper', 'VAT1L_Ca']

print("=" * 80)
print("Phase 11: NVU 4-node Nonlinear Coupling Analysis")
print("=" * 80)
print()

# ============================================================
# Step 1: Load and prepare data
# ============================================================

print("[1] Loading and preparing data...")

# Load NVU 4-node data
df_nvu = pd.read_csv(NVU_DATA_FILE)

# Load additional pathology signatures for severity calculation
df_path = pd.read_csv(PATHOLOGY_FILE)

# Merge, keeping NVU vars from df_nvu and adding only extra pathology columns
# Avoid column name collisions by suffixing duplicates
df = df_nvu.merge(df_path, on='patient_id', how='left', suffixes=('', '_path'))

print(f"  Total patients: {len(df)}")

# Split by group
df_all = df[df[NVU_VARS].notna().all(axis=1)].copy()
df_als = df_all[df_all['is_ALS']].copy()
df_ctrl = df_all[~df_all['is_ALS']].copy()

print(f"  All patients with complete data: n={len(df_all)}")
print(f"  ALS patients: n={len(df_als)}")
print(f"  Control patients: n={len(df_ctrl)}")
print()

# Basic statistics
print("  Variable statistics (all patients):")
for var in NVU_VARS:
    vals = df_all[var].values
    print(f"    {var:15s}: mean={np.mean(vals):.3f}, std={np.std(vals):.3f}, "
          f"min={np.min(vals):.3f}, max={np.max(vals):.3f}")

# Zero-inflation in Vasc_driver
n_zeros = (df_all['Vasc_driver'] == 0).sum()
print(f"  Vasc_driver zero-inflation: {n_zeros}/{len(df_all)} ({n_zeros/len(df_all):.1%})")
print()

# ============================================================
# Step 2: Construct Severity covariate
# ============================================================

print("[2] Constructing disease Severity covariate...")

# Define severity columns (z-score and average)
severity_cols = ['Upper_ER', 'Upper_Mito', 'Glia_ER', 'Glia_Mito', 'VAT1L_Stress']

# Check which columns are available
available_severity = [col for col in severity_cols if col in df_all.columns]
print(f"  Available severity columns: {available_severity}")

if len(available_severity) >= 2:
    # Z-score normalization
    scaler = StandardScaler()
    severity_z = scaler.fit_transform(df_all[available_severity])

    # Severity = mean of z-scores
    df_all['Severity'] = severity_z.mean(axis=1)

    # Apply same transformation to subsets
    df_als['Severity'] = scaler.transform(df_als[available_severity]).mean(axis=1)
    df_ctrl['Severity'] = scaler.transform(df_ctrl[available_severity]).mean(axis=1)

    print(f"  Severity constructed from {len(available_severity)} variables")
    print(f"  Severity range (all): [{df_all['Severity'].min():.3f}, {df_all['Severity'].max():.3f}]")
else:
    print("  Warning: Insufficient severity columns, using Upper_ER only")
    df_all['Severity'] = StandardScaler().fit_transform(df_all[['Upper_ER']])
    df_als['Severity'] = StandardScaler().fit_transform(df_als[['Upper_ER']])
    df_ctrl['Severity'] = StandardScaler().fit_transform(df_ctrl[['Upper_ER']])

print()

# ============================================================
# Step 3: Pairwise nonlinear analysis
# ============================================================

print("[3] Running pairwise nonlinear analysis...")

def compute_residuals(X, Y, covariate):
    """Compute residuals after regressing out covariate."""
    lr_x = LinearRegression()
    lr_y = LinearRegression()

    lr_x.fit(covariate.reshape(-1, 1), X)
    lr_y.fit(covariate.reshape(-1, 1), Y)

    X_res = X - lr_x.predict(covariate.reshape(-1, 1))
    Y_res = Y - lr_y.predict(covariate.reshape(-1, 1))

    return X_res, Y_res

def analyze_pair(df_sub, X_name, Y_name, subset_name):
    """Analyze relationship between two variables."""
    X = df_sub[X_name].values
    Y = df_sub[Y_name].values
    n = len(X)

    if n < 5:
        return None

    # Spearman correlation
    rho_spear, p_spear = spearmanr(X, Y)

    # Linear regression
    lr = LinearRegression()
    lr.fit(X.reshape(-1, 1), Y)
    y_pred_lin = lr.predict(X.reshape(-1, 1))
    r2_linear = r2_score(Y, y_pred_lin)

    # Cubic polynomial regression
    poly = PolynomialFeatures(degree=3)
    X_poly = poly.fit_transform(X.reshape(-1, 1))
    lr_poly = LinearRegression()
    lr_poly.fit(X_poly, Y)
    y_pred_poly = lr_poly.predict(X_poly)
    r2_cubic = r2_score(Y, y_pred_poly)
    r2_improv = r2_cubic - r2_linear

    # Partial correlation (severity-adjusted)
    if 'Severity' in df_sub.columns:
        S = df_sub['Severity'].values
        X_res, Y_res = compute_residuals(X, Y, S)
        rho_partial, p_partial = spearmanr(X_res, Y_res)
    else:
        rho_partial, p_partial = np.nan, np.nan

    return {
        'subset': subset_name,
        'X': X_name,
        'Y': Y_name,
        'n': n,
        'spearman_r': rho_spear,
        'spearman_p': p_spear,
        'R2_linear': r2_linear,
        'R2_cubic': r2_cubic,
        'R2_improv': r2_improv,
        'partial_spear_r': rho_partial,
        'partial_spear_p': p_partial
    }

# Analyze all pairs for each subset
results = []

for subset_name, df_sub in [('all', df_all), ('ALS', df_als), ('Control', df_ctrl)]:
    if len(df_sub) < 5:
        continue

    print(f"  Analyzing {subset_name} subset (n={len(df_sub)})...")

    for X_name, Y_name in combinations(NVU_VARS, 2):
        # X -> Y
        res_xy = analyze_pair(df_sub, X_name, Y_name, subset_name)
        if res_xy:
            results.append(res_xy)

        # Y -> X
        res_yx = analyze_pair(df_sub, Y_name, X_name, subset_name)
        if res_yx:
            results.append(res_yx)

df_pairwise = pd.DataFrame(results)

# Save
output_file = f'{OUTPUT_DIR}/nvu_pairwise_nonlin.csv'
df_pairwise.to_csv(output_file, index=False)
print(f"  Saved: {output_file}")
print()

# Print top relationships
print("  Top relationships by partial Spearman (|rho| > 0.3):")
top_partial = df_pairwise[df_pairwise['partial_spear_r'].abs() > 0.3].sort_values(
    'partial_spear_r', key=abs, ascending=False
)
for _, row in top_partial.head(10).iterrows():
    sig = '***' if row['partial_spear_p'] < 0.001 else '**' if row['partial_spear_p'] < 0.01 else '*' if row['partial_spear_p'] < 0.05 else ''
    print(f"    [{row['subset']:7s}] {row['X']:15s} → {row['Y']:15s}: "
          f"ρ_partial={row['partial_spear_r']:+.3f} (p={row['partial_spear_p']:.3f}){sig}, "
          f"R²_improv={row['R2_improv']:.3f}")

print()

# ============================================================
# Step 4: Direction scoring
# ============================================================

print("[4] Computing direction scores...")

def compute_direction_score(df_pw, X_name, Y_name, subset_name):
    """Compute direction score for X -> Y."""
    row_xy = df_pw[(df_pw['subset'] == subset_name) &
                   (df_pw['X'] == X_name) &
                   (df_pw['Y'] == Y_name)]

    if len(row_xy) == 0:
        return 0.0

    row = row_xy.iloc[0]

    # Score = weighted combination of evidence
    # Higher score = more evidence for X -> Y direction
    w1 = 1.0  # Spearman weight
    w2 = 1.0  # Partial correlation weight
    w3 = 0.5  # R2 improvement weight

    score = (w1 * abs(row['spearman_r']) +
             w2 * abs(row['partial_spear_r']) +
             w3 * max(0, row['R2_improv']))

    return score

# Compute direction scores for all pairs
direction_results = []

for subset_name in ['all', 'ALS', 'Control']:
    df_sub_pw = df_pairwise[df_pairwise['subset'] == subset_name]

    if len(df_sub_pw) == 0:
        continue

    for X_name, Y_name in combinations(NVU_VARS, 2):
        score_xy = compute_direction_score(df_pairwise, X_name, Y_name, subset_name)
        score_yx = compute_direction_score(df_pairwise, Y_name, X_name, subset_name)

        # Determine direction preference
        if score_xy > score_yx * 1.2:  # 20% threshold
            direction = f'{X_name}→{Y_name}'
        elif score_yx > score_xy * 1.2:
            direction = f'{Y_name}→{X_name}'
        else:
            direction = 'undetermined'

        # Get partial correlation for preferred direction
        if direction.startswith(X_name):
            partial_r = df_sub_pw[(df_sub_pw['X']==X_name) & (df_sub_pw['Y']==Y_name)]['partial_spear_r'].values
            r2_improv = df_sub_pw[(df_sub_pw['X']==X_name) & (df_sub_pw['Y']==Y_name)]['R2_improv'].values
        elif direction.startswith(Y_name):
            partial_r = df_sub_pw[(df_sub_pw['X']==Y_name) & (df_sub_pw['Y']==X_name)]['partial_spear_r'].values
            r2_improv = df_sub_pw[(df_sub_pw['X']==Y_name) & (df_sub_pw['Y']==X_name)]['R2_improv'].values
        else:
            partial_r = [np.nan]
            r2_improv = [np.nan]

        direction_results.append({
            'subset': subset_name,
            'node1': X_name,
            'node2': Y_name,
            'score_12': score_xy,
            'score_21': score_yx,
            'direction_preference': direction,
            'partial_r': partial_r[0] if len(partial_r) > 0 else np.nan,
            'R2_improv': r2_improv[0] if len(r2_improv) > 0 else np.nan
        })

df_direction = pd.DataFrame(direction_results)

# Save
output_file = f'{OUTPUT_DIR}/nvu_direction_scores.csv'
df_direction.to_csv(output_file, index=False)
print(f"  Saved: {output_file}")
print()

# Print direction candidates
print("  Direction candidates (ALS subset):")
df_dir_als = df_direction[df_direction['subset'] == 'ALS']
for _, row in df_dir_als.iterrows():
    if row['direction_preference'] != 'undetermined':
        print(f"    {row['direction_preference']:30s} (score ratio: "
              f"{max(row['score_12'], row['score_21']):.2f}/{min(row['score_12'], row['score_21']):.2f}, "
              f"ρ_partial={row['partial_r']:+.3f})")

print()

# ============================================================
# Step 5: Visualizations
# ============================================================

print("[5] Creating visualizations...")

# Heatmaps
def plot_correlation_heatmap(df_pw, subset_name, metric='spearman_r', title_suffix=''):
    """Plot correlation heatmap."""
    # Create matrix
    matrix = np.zeros((len(NVU_VARS), len(NVU_VARS)))

    for i, var1 in enumerate(NVU_VARS):
        for j, var2 in enumerate(NVU_VARS):
            if i == j:
                matrix[i, j] = 1.0
            else:
                row = df_pw[(df_pw['subset'] == subset_name) &
                           (df_pw['X'] == var1) &
                           (df_pw['Y'] == var2)]
                if len(row) > 0:
                    matrix[i, j] = row[metric].values[0]

    # Plot
    fig, ax = plt.subplots(figsize=(8, 7))
    sns.heatmap(matrix, annot=True, fmt='.3f', cmap='RdBu_r', center=0,
                vmin=-1, vmax=1, xticklabels=NVU_VARS, yticklabels=NVU_VARS,
                cbar_kws={'label': metric.replace('_', ' ').title()}, ax=ax)
    ax.set_title(f'{metric.replace("_", " ").title()} - {subset_name} patients {title_suffix}',
                fontsize=13, fontweight='bold')

    plt.tight_layout()
    filename = f'{OUTPUT_DIR}/nvu_{metric}_heatmap_{subset_name}.png'
    plt.savefig(filename, dpi=150, bbox_inches='tight')
    plt.close()

    return filename

# Plot heatmaps for each subset
for subset in ['all', 'ALS', 'Control']:
    df_sub = df_pairwise[df_pairwise['subset'] == subset]
    if len(df_sub) > 0:
        plot_correlation_heatmap(df_sub, subset, 'spearman_r')
        plot_correlation_heatmap(df_sub, subset, 'partial_spear_r', '(Severity-adjusted)')

print(f"  Saved correlation heatmaps")

# Direction candidate graph
def plot_direction_graph(df_dir, subset_name):
    """Plot direction candidate graph."""
    fig, ax = plt.subplots(figsize=(10, 8))

    G = nx.DiGraph()

    # Add nodes
    for var in NVU_VARS:
        G.add_node(var)

    # Add edges based on direction preference
    df_sub = df_dir[df_dir['subset'] == subset_name]

    for _, row in df_sub.iterrows():
        if row['direction_preference'] != 'undetermined':
            parts = row['direction_preference'].split('→')
            if len(parts) == 2:
                src, dst = parts
                weight = max(row['score_12'], row['score_21'])
                G.add_edge(src, dst, weight=weight, partial_r=row['partial_r'])

    # Layout
    pos = nx.spring_layout(G, k=2, iterations=50, seed=42)

    # Draw nodes
    nx.draw_networkx_nodes(G, pos, node_size=3000, node_color='lightblue',
                          edgecolors='black', linewidths=2, ax=ax)

    # Draw edges
    for (src, dst) in G.edges():
        weight = G[src][dst]['weight']
        partial_r = G[src][dst]['partial_r']
        edge_color = 'red' if partial_r > 0 else 'blue'
        edge_width = min(weight * 2, 5)
        nx.draw_networkx_edges(G, pos, [(src, dst)],
                              width=edge_width, edge_color=edge_color,
                              arrowsize=20, arrowstyle='->', ax=ax)

    # Draw labels
    nx.draw_networkx_labels(G, pos, font_size=10, font_weight='bold', ax=ax)

    # Edge labels
    edge_labels = {(src, dst): f'{G[src][dst]["partial_r"]:+.2f}'
                   for (src, dst) in G.edges()}
    nx.draw_networkx_edge_labels(G, pos, edge_labels, font_size=8, ax=ax)

    ax.set_title(f'Direction Candidates ({subset_name} patients)',
                fontsize=14, fontweight='bold')
    ax.axis('off')

    plt.tight_layout()
    filename = f'{OUTPUT_DIR}/nvu_direction_candidate_graph_{subset_name}.png'
    plt.savefig(filename, dpi=150, bbox_inches='tight')
    plt.close()

    return filename

# Plot direction graphs
for subset in ['all', 'ALS']:
    df_sub = df_direction[df_direction['subset'] == subset]
    if len(df_sub) > 0:
        plot_direction_graph(df_direction, subset)

print(f"  Saved direction candidate graphs")
print()

# ============================================================
# Completion
# ============================================================

print("=" * 80)
print("Phase 11 Complete!")
print("=" * 80)
print()
print("Key outputs:")
print(f"  - {OUTPUT_DIR}/nvu_pairwise_nonlin.csv")
print(f"  - {OUTPUT_DIR}/nvu_direction_scores.csv")
print(f"  - {OUTPUT_DIR}/nvu_*_heatmap_*.png")
print(f"  - {OUTPUT_DIR}/nvu_direction_candidate_graph_*.png")
print()
print("Next: Create comprehensive Phase 11 report")
