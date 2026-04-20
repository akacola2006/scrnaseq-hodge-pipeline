#!/usr/bin/env python3
"""
Phase 3c: LiNGAM因果方向の精緻化

NOTEARSで得た構造スケルトンに対して、LiNGAMによる
非ガウス性ベースの因果方向推定を行い、エッジを検証・分類する。

理論:
    LiNGAM (Linear Non-Gaussian Acyclic Model):
    x = Bx + e
    where B is strictly lower triangular (after permutation)
          e_i ~ non-Gaussian

    Non-Gaussianity allows identifying causal direction
    without interventions.

    DirectLiNGAM:
    - Directly estimates causal ordering via ICA
    - No need for prior structure
    - Robust to moderate sample sizes

Step 3a: Full DirectLiNGAM (24 modules, unconstrained)
Step 3b: Skeleton-constrained comparison with NOTEARS
Step 3c: Focused triplet/cluster tests

Input:
- results/causal_graphs/notears_data_VAT1L.npy
- results/causal_graphs/notears_adj_matrix_VAT1L.csv
- results/causal_graphs/notears_modules_VAT1L.txt

Output:
- results/causal_graphs/lingam_adj_matrix_VAT1L.csv
- results/causal_graphs/lingam_causal_order_VAT1L.csv
- results/causal_graphs/core_dag_VAT1L.csv (NOTEARS ∩ LiNGAM)
- results/causal_graphs/edge_classification_VAT1L.csv
- results/figures/lingam_comparison_VAT1L.png

Usage:
    python 03c_lingam_refinement.py
"""

import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx
from pathlib import Path
import json
from scipy import stats

# =====================================================
# パス設定
# =====================================================

BASE_DIR = Path("/home/akaco/als/motor_cortex_analysis")
IDS_DIR = BASE_DIR / "ids_causal_analysis"
GRAPH_DIR = IDS_DIR / "results" / "causal_graphs"
FIG_DIR = IDS_DIR / "results" / "figures"
RESULTS_DIR = IDS_DIR / "results" / "ids_pseudotime"

# =====================================================
# LiNGAM実装
# =====================================================

try:
    from lingam import DirectLiNGAM
    LINGAM_AVAILABLE = True
except ImportError:
    print("Warning: lingam package not available. Will use ICA-based fallback.")
    LINGAM_AVAILABLE = False


def check_non_gaussianity(X, module_names):
    """
    各モジュールの非ガウス性を検定

    Parameters
    ----------
    X : np.ndarray
        データ行列 (n_samples, n_features)
    module_names : list
        モジュール名

    Returns
    -------
    non_gaussian_scores : pd.DataFrame
        非ガウス性スコア
    """
    print(f"\n{'='*60}")
    print(f"Checking Non-Gaussianity")
    print(f"{'='*60}")

    results = []

    for i, module in enumerate(module_names):
        x = X[:, i]

        # Shapiro-Wilk test
        stat_sw, p_sw = stats.shapiro(x)

        # Jarque-Bera test
        stat_jb, p_jb = stats.jarque_bera(x)

        # Kurtosis (excess kurtosis)
        kurt = stats.kurtosis(x)

        # Skewness
        skew = stats.skew(x)

        results.append({
            'module': module,
            'shapiro_stat': stat_sw,
            'shapiro_p': p_sw,
            'jarque_bera_stat': stat_jb,
            'jarque_bera_p': p_jb,
            'kurtosis': kurt,
            'skewness': skew,
            'is_non_gaussian': (p_sw < 0.05) or (p_jb < 0.05)
        })

    df = pd.DataFrame(results)

    n_non_gaussian = df['is_non_gaussian'].sum()
    print(f"\n  Non-Gaussian modules: {n_non_gaussian} / {len(module_names)}")
    print(f"  Mean |kurtosis|: {np.abs(df['kurtosis']).mean():.3f}")
    print(f"  Mean |skewness|: {np.abs(df['skewness']).mean():.3f}")

    print(f"\n  Top 5 non-Gaussian modules (by |kurtosis|):")
    top_ng = df.nlargest(5, 'kurtosis', keep='all')
    for _, row in top_ng.iterrows():
        print(f"    {row['module']:25s}: kurt={row['kurtosis']:+.3f}, skew={row['skewness']:+.3f}")

    return df


def run_direct_lingam(X, module_names):
    """
    DirectLiNGAMによる因果順序推定

    Parameters
    ----------
    X : np.ndarray
        データ行列 (n_samples, n_features)
    module_names : list
        モジュール名

    Returns
    -------
    adjacency : np.ndarray
        隣接行列 (n_features, n_features)
    causal_order : np.ndarray
        因果順序 (permutation)
    """
    print(f"\n{'='*60}")
    print(f"Running DirectLiNGAM")
    print(f"{'='*60}")
    print(f"  Data shape: {X.shape}")

    if LINGAM_AVAILABLE:
        # Use official lingam package
        model = DirectLiNGAM()
        model.fit(X)

        adjacency = model.adjacency_matrix_
        causal_order = model.causal_order_

    else:
        # Fallback: simple ICA-based heuristic
        from sklearn.decomposition import FastICA

        print("  Using ICA-based fallback (not true DirectLiNGAM)")

        # ICA
        ica = FastICA(n_components=X.shape[1], random_state=42, max_iter=1000)
        S = ica.fit_transform(X)
        A = ica.mixing_

        # Estimate adjacency via regression
        from sklearn.linear_model import LassoCV

        adjacency = np.zeros((X.shape[1], X.shape[1]))

        for i in range(X.shape[1]):
            # Regress X[:, i] on all other X[:, j]
            mask = np.ones(X.shape[1], dtype=bool)
            mask[i] = False

            if X[:, mask].shape[1] > 0:
                model = LassoCV(cv=5, random_state=42)
                model.fit(X[:, mask], X[:, i])

                adjacency[mask, i] = model.coef_

        # Causal order: heuristic based on total out-degree
        out_degrees = np.sum(np.abs(adjacency) > 0.01, axis=1)
        causal_order = np.argsort(-out_degrees)  # High out-degree = upstream

    print(f"\n  Causal order (upstream → downstream):")
    for rank, idx in enumerate(causal_order[:10]):
        print(f"    {rank+1:2d}. {module_names[idx]:25s}")

    print(f"\n  Adjacency matrix:")
    print(f"    Non-zero edges: {np.sum(np.abs(adjacency) > 0.01)}")
    print(f"    Max weight: {np.abs(adjacency).max():.4f}")

    return adjacency, causal_order


def compare_with_notears(W_notears, W_lingam, module_names, threshold=0.01):
    """
    NOTEARSとLiNGAMのエッジを比較・分類

    Parameters
    ----------
    W_notears : np.ndarray
        NOTEARS隣接行列
    W_lingam : np.ndarray
        LiNGAM隣接行列
    module_names : list
        モジュール名
    threshold : float
        エッジ判定閾値

    Returns
    -------
    edge_classification : pd.DataFrame
        エッジ分類結果
    """
    print(f"\n{'='*60}")
    print(f"Comparing NOTEARS vs LiNGAM")
    print(f"{'='*60}")

    edges = []

    n = len(module_names)

    for i in range(n):
        for j in range(n):
            if i == j:
                continue

            w_nt = W_notears[i, j]
            w_lg = W_lingam[i, j]
            w_lg_rev = W_lingam[j, i]

            # NOTEARS edge exists?
            nt_edge = abs(w_nt) > threshold

            # LiNGAM edge exists (same direction)?
            lg_edge_same = abs(w_lg) > threshold

            # LiNGAM edge exists (reverse direction)?
            lg_edge_rev = abs(w_lg_rev) > threshold

            if nt_edge or lg_edge_same or lg_edge_rev:
                # Classify
                if nt_edge and lg_edge_same:
                    category = 'core'  # Both agree
                elif nt_edge and lg_edge_rev:
                    category = 'ambiguous_flip'  # Direction disagreement
                elif nt_edge and not lg_edge_same and not lg_edge_rev:
                    category = 'notears_only'  # NOTEARS artifact?
                elif lg_edge_same and not nt_edge:
                    category = 'lingam_only'  # Weak but real?
                elif lg_edge_rev and not nt_edge:
                    category = 'lingam_reverse_only'
                else:
                    category = 'other'

                edges.append({
                    'source': module_names[i],
                    'target': module_names[j],
                    'notears_weight': w_nt,
                    'lingam_weight': w_lg,
                    'lingam_weight_reverse': w_lg_rev,
                    'category': category
                })

    df = pd.DataFrame(edges)

    # Statistics
    print(f"\n  Edge classification:")
    for cat in ['core', 'ambiguous_flip', 'notears_only', 'lingam_only', 'lingam_reverse_only']:
        n_edges = (df['category'] == cat).sum()
        print(f"    {cat:25s}: {n_edges:3d} edges")

    print(f"\n  Core edges (Top 10 by NOTEARS weight):")
    core = df[df['category'] == 'core'].nlargest(10, 'notears_weight', keep='all')
    for _, row in core.iterrows():
        print(f"    {row['source']:20s} → {row['target']:20s}: NT={row['notears_weight']:.3f}, LG={row['lingam_weight']:.3f}")

    print(f"\n  Ambiguous edges (direction flip):")
    ambig = df[df['category'] == 'ambiguous_flip']
    for _, row in ambig.iterrows():
        print(f"    {row['source']:20s} ⇄ {row['target']:20s}: NT={row['notears_weight']:.3f}, LG_rev={row['lingam_weight_reverse']:.3f}")

    return df


def create_core_dag(edge_classification, module_names):
    """
    CoreエッジのみからDAGを構築

    Parameters
    ----------
    edge_classification : pd.DataFrame
        エッジ分類
    module_names : list
        モジュール名

    Returns
    -------
    W_core : np.ndarray
        Core隣接行列
    """
    print(f"\n{'='*60}")
    print(f"Creating Core DAG")
    print(f"{'='*60}")

    n = len(module_names)
    W_core = np.zeros((n, n))

    module_to_idx = {m: i for i, m in enumerate(module_names)}

    core_edges = edge_classification[edge_classification['category'] == 'core']

    for _, row in core_edges.iterrows():
        i = module_to_idx[row['source']]
        j = module_to_idx[row['target']]
        W_core[i, j] = row['notears_weight']

    # Check DAG property
    G = nx.DiGraph(W_core)
    is_dag = nx.is_directed_acyclic_graph(G)

    n_edges = np.sum(W_core != 0)

    print(f"\n  Core DAG:")
    print(f"    Edges: {n_edges}")
    print(f"    Is DAG: {is_dag}")

    # Source/sink analysis
    in_degrees = dict(G.in_degree())
    out_degrees = dict(G.out_degree())

    sources = [i for i in range(n) if out_degrees[i] > 0 and in_degrees[i] == 0]
    sinks = [i for i in range(n) if in_degrees[i] > 0 and out_degrees[i] == 0]

    print(f"\n  Source nodes: {len(sources)}")
    for i in sources:
        print(f"    {module_names[i]:25s} (out={out_degrees[i]})")

    print(f"\n  Sink nodes: {len(sinks)}")
    for i in sinks[:5]:
        print(f"    {module_names[i]:25s} (in={in_degrees[i]})")

    return W_core


def run_focused_triplet_tests(X, module_names, triplets):
    """
    特定のモジュール組み合わせでLiNGAMを実行

    Parameters
    ----------
    X : np.ndarray
        データ行列
    module_names : list
        モジュール名
    triplets : list of list
        テストするモジュール組み合わせ

    Returns
    -------
    results : dict
        各tripletの結果
    """
    print(f"\n{'='*60}")
    print(f"Running Focused Triplet Tests")
    print(f"{'='*60}")

    module_to_idx = {m: i for i, m in enumerate(module_names)}

    results = {}

    for triplet in triplets:
        print(f"\n  Triplet: {triplet}")

        # Extract indices
        indices = [module_to_idx[m] for m in triplet]
        X_sub = X[:, indices]

        # Run DirectLiNGAM
        if LINGAM_AVAILABLE:
            model = DirectLiNGAM()
            model.fit(X_sub)
            adj_sub = model.adjacency_matrix_
            order_sub = model.causal_order_
        else:
            # Simple fallback
            from sklearn.linear_model import LassoCV
            adj_sub = np.zeros((len(indices), len(indices)))
            for i in range(len(indices)):
                mask = np.ones(len(indices), dtype=bool)
                mask[i] = False
                if X_sub[:, mask].shape[1] > 0:
                    model = LassoCV(cv=5, random_state=42)
                    model.fit(X_sub[:, mask], X_sub[:, i])
                    adj_sub[mask, i] = model.coef_
            order_sub = np.argsort(-np.sum(np.abs(adj_sub) > 0.01, axis=1))

        # Report
        print(f"    Causal order:")
        for rank, idx in enumerate(order_sub):
            print(f"      {rank+1}. {triplet[idx]}")

        print(f"    Edges:")
        for i in range(len(indices)):
            for j in range(len(indices)):
                if abs(adj_sub[i, j]) > 0.01:
                    print(f"      {triplet[i]:20s} → {triplet[j]:20s}: {adj_sub[i, j]:.3f}")

        results[tuple(triplet)] = {
            'adjacency': adj_sub,
            'causal_order': [triplet[idx] for idx in order_sub],
            'modules': triplet
        }

    return results


# =====================================================
# 可視化
# =====================================================

def visualize_lingam_comparison(W_notears, W_lingam, W_core,
                                causal_order, module_names, coeffs,
                                edge_classification, output_file):
    """
    NOTEARS vs LiNGAM の比較可視化
    """
    fig = plt.figure(figsize=(20, 16))
    gs = fig.add_gridspec(3, 3, hspace=0.3, wspace=0.3)

    # 1. NOTEARS adjacency (sorted by β)
    ax1 = fig.add_subplot(gs[0, 0])
    sorted_idx = coeffs.sort_values(ascending=False).index
    module_to_idx = {m: i for i, m in enumerate(module_names)}
    sorted_indices = [module_to_idx[m] for m in sorted_idx]

    W_nt_sorted = W_notears[sorted_indices, :][:, sorted_indices]

    sns.heatmap(W_nt_sorted, cmap='RdYlBu_r', center=0, ax=ax1,
                xticklabels=sorted_idx, yticklabels=sorted_idx,
                cbar_kws={'label': 'Weight'})
    ax1.set_title('NOTEARS Adjacency\n(sorted by β)', fontsize=11)
    ax1.set_xticklabels(ax1.get_xticklabels(), rotation=90, fontsize=6)
    ax1.set_yticklabels(ax1.get_yticklabels(), rotation=0, fontsize=6)

    # 2. LiNGAM adjacency (sorted by causal order)
    ax2 = fig.add_subplot(gs[0, 1])

    W_lg_sorted = W_lingam[causal_order, :][:, causal_order]
    ordered_names = [module_names[i] for i in causal_order]

    sns.heatmap(W_lg_sorted, cmap='RdYlBu_r', center=0, ax=ax2,
                xticklabels=ordered_names, yticklabels=ordered_names,
                cbar_kws={'label': 'Weight'})
    ax2.set_title('LiNGAM Adjacency\n(sorted by causal order)', fontsize=11)
    ax2.set_xticklabels(ax2.get_xticklabels(), rotation=90, fontsize=6)
    ax2.set_yticklabels(ax2.get_yticklabels(), rotation=0, fontsize=6)

    # 3. Core DAG (sorted by β)
    ax3 = fig.add_subplot(gs[0, 2])

    W_core_sorted = W_core[sorted_indices, :][:, sorted_indices]

    sns.heatmap(W_core_sorted, cmap='RdYlBu_r', center=0, ax=ax3,
                xticklabels=sorted_idx, yticklabels=sorted_idx,
                cbar_kws={'label': 'Weight'})
    ax3.set_title('Core DAG (NOTEARS ∩ LiNGAM)\n(sorted by β)', fontsize=11)
    ax3.set_xticklabels(ax3.get_xticklabels(), rotation=90, fontsize=6)
    ax3.set_yticklabels(ax3.get_yticklabels(), rotation=0, fontsize=6)

    # 4. Edge category distribution
    ax4 = fig.add_subplot(gs[1, 0])

    cat_counts = edge_classification['category'].value_counts()
    colors = {'core': 'green', 'ambiguous_flip': 'orange',
              'notears_only': 'red', 'lingam_only': 'blue',
              'lingam_reverse_only': 'purple'}

    bars = ax4.bar(range(len(cat_counts)), cat_counts.values,
                   color=[colors.get(c, 'gray') for c in cat_counts.index],
                   alpha=0.7)
    ax4.set_xticks(range(len(cat_counts)))
    ax4.set_xticklabels(cat_counts.index, rotation=45, ha='right', fontsize=8)
    ax4.set_ylabel('Number of edges', fontsize=10)
    ax4.set_title('Edge Classification', fontsize=11)
    ax4.grid(alpha=0.3, axis='y')

    for i, v in enumerate(cat_counts.values):
        ax4.text(i, v, str(v), ha='center', va='bottom', fontsize=9)

    # 5. Causal order comparison (β rank vs LiNGAM rank)
    ax5 = fig.add_subplot(gs[1, 1])

    beta_rank = {m: i for i, m in enumerate(coeffs.sort_values(ascending=False).index)}
    lingam_rank = {module_names[causal_order[i]]: i for i in range(len(causal_order))}

    common_modules = [m for m in module_names if m in beta_rank and m in lingam_rank]

    x = [beta_rank[m] for m in common_modules]
    y = [lingam_rank[m] for m in common_modules]

    ax5.scatter(x, y, alpha=0.6, s=50)

    # Add labels for key modules
    key_modules = ['Transcription', 'Angiogenesis', 'Mitochondria', 'ECM', 'ER_Stress']
    for m in key_modules:
        if m in common_modules:
            ax5.annotate(m, (beta_rank[m], lingam_rank[m]), fontsize=7, alpha=0.8)

    ax5.plot([0, len(common_modules)], [0, len(common_modules)], 'r--', alpha=0.5, label='Perfect agreement')
    ax5.set_xlabel('β rank (high → low)', fontsize=10)
    ax5.set_ylabel('LiNGAM causal order', fontsize=10)
    ax5.set_title('Causal Order Comparison', fontsize=11)
    ax5.legend(fontsize=8)
    ax5.grid(alpha=0.3)

    # 6. Core edges network
    ax6 = fig.add_subplot(gs[1, 2])

    G_core = nx.DiGraph(W_core)

    # Only show edges
    if G_core.number_of_edges() > 0:
        # Position by β
        pos = {}
        for node in G_core.nodes():
            beta = coeffs.iloc[node]
            pos[node] = (np.random.randn()*0.3, beta)

        pos = nx.spring_layout(G_core, pos=pos, k=2, iterations=50)

        node_colors = [coeffs.iloc[node] for node in G_core.nodes()]

        nx.draw_networkx_nodes(G_core, pos, ax=ax6, node_color=node_colors,
                              cmap='RdYlBu_r', node_size=300, alpha=0.8)
        nx.draw_networkx_labels(G_core, pos, ax=ax6,
                               labels={i: module_names[i] for i in G_core.nodes()},
                               font_size=6)
        nx.draw_networkx_edges(G_core, pos, ax=ax6,
                              edge_color='gray', arrows=True,
                              arrowsize=10, alpha=0.6,
                              connectionstyle='arc3,rad=0.1')

    ax6.set_title(f'Core DAG Network ({G_core.number_of_edges()} edges)', fontsize=11)
    ax6.axis('off')

    # 7. Agreement score per module
    ax7 = fig.add_subplot(gs[2, :])

    # Count core edges per module (as source)
    module_agreement = []
    for i, module in enumerate(module_names):
        n_notears = np.sum(W_notears[i, :] != 0)
        n_core = np.sum(W_core[i, :] != 0)

        if n_notears > 0:
            agreement = n_core / n_notears
        else:
            agreement = 1.0 if n_core == 0 else 0.0

        module_agreement.append({
            'module': module,
            'notears_edges': n_notears,
            'core_edges': n_core,
            'agreement': agreement,
            'beta': coeffs[module]
        })

    df_agree = pd.DataFrame(module_agreement)
    df_agree = df_agree.sort_values('beta', ascending=False)

    x = range(len(df_agree))
    width = 0.35

    ax7.bar([i - width/2 for i in x], df_agree['notears_edges'], width,
            label='NOTEARS edges', alpha=0.7, color='coral')
    ax7.bar([i + width/2 for i in x], df_agree['core_edges'], width,
            label='Core edges', alpha=0.7, color='green')

    ax7.set_xticks(x)
    ax7.set_xticklabels(df_agree['module'], rotation=90, fontsize=7)
    ax7.set_ylabel('Number of outgoing edges', fontsize=10)
    ax7.set_title('NOTEARS vs Core Edges per Module (sorted by β)', fontsize=11)
    ax7.legend(fontsize=9)
    ax7.grid(alpha=0.3, axis='y')

    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"\n✓ Saved: {output_file}")
    plt.close()


# =====================================================
# メイン処理
# =====================================================

def main():
    # Load data
    print(f"\n{'='*60}")
    print(f"Loading data")
    print(f"{'='*60}")

    data_file = GRAPH_DIR / "notears_data_VAT1L.npy"
    notears_file = GRAPH_DIR / "notears_adj_matrix_VAT1L.csv"
    modules_file = GRAPH_DIR / "notears_modules_VAT1L.txt"
    coeffs_file = RESULTS_DIR / "local_causal_coeffs_VAT1L.csv"

    X = np.load(data_file)
    W_notears = pd.read_csv(notears_file, index_col=0).values

    with open(modules_file, 'r') as f:
        module_names = [line.strip() for line in f]

    coeffs = pd.read_csv(coeffs_file, index_col=0)['coefficient']

    print(f"  Data: {X.shape}")
    print(f"  NOTEARS adjacency: {W_notears.shape}")
    print(f"  Modules: {len(module_names)}")

    # Step 3a: Check non-Gaussianity
    non_gaussian_df = check_non_gaussianity(X, module_names)
    ng_file = GRAPH_DIR / "non_gaussianity_VAT1L.csv"
    non_gaussian_df.to_csv(ng_file, index=False)
    print(f"\n✓ Saved: {ng_file}")

    # Step 3a: Full DirectLiNGAM
    W_lingam, causal_order = run_direct_lingam(X, module_names)

    # Save LiNGAM results
    lingam_adj_file = GRAPH_DIR / "lingam_adj_matrix_VAT1L.csv"
    W_lingam_df = pd.DataFrame(W_lingam, index=module_names, columns=module_names)
    W_lingam_df.to_csv(lingam_adj_file)
    print(f"\n✓ Saved: {lingam_adj_file}")

    causal_order_file = GRAPH_DIR / "lingam_causal_order_VAT1L.csv"
    order_df = pd.DataFrame({
        'rank': range(1, len(causal_order) + 1),
        'module': [module_names[i] for i in causal_order],
        'beta': [coeffs.iloc[i] for i in causal_order]
    })
    order_df.to_csv(causal_order_file, index=False)
    print(f"\n✓ Saved: {causal_order_file}")

    # Step 3b: Compare with NOTEARS
    edge_classification = compare_with_notears(W_notears, W_lingam, module_names)

    edge_class_file = GRAPH_DIR / "edge_classification_VAT1L.csv"
    edge_classification.to_csv(edge_class_file, index=False)
    print(f"\n✓ Saved: {edge_class_file}")

    # Create core DAG
    W_core = create_core_dag(edge_classification, module_names)

    core_dag_file = GRAPH_DIR / "core_dag_VAT1L.csv"
    W_core_df = pd.DataFrame(W_core, index=module_names, columns=module_names)
    W_core_df.to_csv(core_dag_file)
    print(f"\n✓ Saved: {core_dag_file}")

    # Step 3c: Focused triplet tests
    triplets = [
        ['Transcription', 'lncRNA', 'Epigenetic', 'RNA_Processing'],
        ['Transcription', 'Cytoskeleton', 'ALS_Genes'],
        ['Angiogenesis', 'ECM', 'Synaptic', 'Ion_Transport'],
        ['Mitochondria', 'ER_Stress', 'Oxidative_Stress'],
        ['Transcription', 'Protein_Homeostasis', 'DNA_Repair']
    ]

    triplet_results = run_focused_triplet_tests(X, module_names, triplets)

    # Save triplet results
    triplet_file = GRAPH_DIR / "triplet_tests_VAT1L.json"
    triplet_export = {}
    for key, val in triplet_results.items():
        triplet_export[str(key)] = {
            'causal_order': val['causal_order'],
            'modules': val['modules'],
            'edges': []
        }
        adj = val['adjacency']
        for i in range(len(val['modules'])):
            for j in range(len(val['modules'])):
                if abs(adj[i, j]) > 0.01:
                    triplet_export[str(key)]['edges'].append({
                        'source': val['modules'][i],
                        'target': val['modules'][j],
                        'weight': float(adj[i, j])
                    })

    with open(triplet_file, 'w') as f:
        json.dump(triplet_export, f, indent=2)
    print(f"\n✓ Saved: {triplet_file}")

    # Visualization
    fig_file = FIG_DIR / "lingam_comparison_VAT1L.png"
    visualize_lingam_comparison(W_notears, W_lingam, W_core,
                                causal_order, module_names, coeffs,
                                edge_classification, fig_file)

    # Final summary
    print(f"\n{'='*60}")
    print(f"LiNGAM Refinement Complete")
    print(f"{'='*60}")

    print(f"\nCausal order (Top 5):")
    for i in range(min(5, len(causal_order))):
        idx = causal_order[i]
        print(f"  {i+1}. {module_names[idx]:25s} (β={coeffs.iloc[idx]:+.4f})")

    print(f"\nEdge classification:")
    for cat in ['core', 'ambiguous_flip', 'notears_only', 'lingam_only']:
        n = (edge_classification['category'] == cat).sum()
        print(f"  {cat:25s}: {n:3d} edges")

    print(f"\nCore DAG:")
    print(f"  Edges: {np.sum(W_core != 0)}")
    G_core = nx.DiGraph(W_core)
    print(f"  Is DAG: {nx.is_directed_acyclic_graph(G_core)}")

    # Source analysis
    in_deg = dict(G_core.in_degree())
    out_deg = dict(G_core.out_degree())
    sources = [i for i in range(len(module_names)) if out_deg[i] > 0 and in_deg[i] == 0]

    print(f"\n  Source nodes:")
    for i in sources:
        print(f"    {module_names[i]:25s} (β={coeffs.iloc[i]:+.4f}, out={out_deg[i]})")

    print(f"\n✓ Ready for RNA velocity integration (Step 4)")


if __name__ == "__main__":
    main()
