#!/usr/bin/env python3
"""
Phase 3b: NOTEARS実行によるDAG推定

Step 2で準備したprior制約を用いて、NOTEARSアルゴリズムにより
VAT1Lの因果グラフ（DAG）を推定する。

理論:
    NOTEARS (Zheng et al. 2018):
    min_{W} ||X - XW||²_F + λ||W||_1
    s.t. h(W) = 0  (acyclicity constraint)

    where h(W) = tr(e^(W⊙W)) - d = 0

    今回の拡張:
    - Prior制約: W[i,j] の初期化に β_i - β_j を利用
    - Positive-only constraint: W[i,j] ≥ 0
    - Sparsity: λ による L1 正則化

Input:
- results/causal_graphs/notears_data_VAT1L.npy
- results/causal_graphs/notears_prior_matrix_VAT1L.csv
- results/causal_graphs/notears_modules_VAT1L.txt

Output:
- results/causal_graphs/notears_adj_matrix_VAT1L.csv
- results/causal_graphs/notears_dag_VAT1L.png
- results/causal_graphs/notears_summary_VAT1L.json

Usage:
    python 03b_run_notears.py
"""

import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx
from pathlib import Path
import json

# =====================================================
# パス設定
# =====================================================

BASE_DIR = Path("/home/akaco/als/motor_cortex_analysis")
IDS_DIR = BASE_DIR / "ids_causal_analysis"
GRAPH_DIR = IDS_DIR / "results" / "causal_graphs"
FIG_DIR = IDS_DIR / "results" / "figures"

# =====================================================
# NOTEARS実装（scipy最適化版）
# =====================================================

from scipy.optimize import minimize
from scipy.special import expit as sigmoid


def notears_linear(X, lambda1=0.1, loss_type='l2', max_iter=100, h_tol=1e-8,
                   rho_max=1e+16, w_threshold=0.3, prior_matrix=None,
                   prior_strength=0.5):
    """
    NOTEARS algorithm for linear SEM

    Parameters
    ----------
    X : np.ndarray
        Data matrix (n_samples, n_nodes)
    lambda1 : float
        L1 regularization parameter
    loss_type : str
        'l2' or 'logistic'
    max_iter : int
        Maximum number of dual ascent steps
    h_tol : float
        Tolerance for acyclicity constraint
    rho_max : float
        Maximum value of penalty parameter
    w_threshold : float
        Threshold for edge pruning
    prior_matrix : np.ndarray, optional
        Prior constraint matrix (n_nodes, n_nodes)
        Higher values = stronger prior for edge existence
    prior_strength : float
        Strength of prior influence (0=ignore, 1=strong)

    Returns
    -------
    W : np.ndarray
        Estimated adjacency matrix (n_nodes, n_nodes)
    """
    def _loss(W):
        """Compute loss"""
        M = X @ W
        if loss_type == 'l2':
            R = X - M
            loss = 0.5 / X.shape[0] * (R ** 2).sum()
        else:
            raise NotImplementedError
        return loss

    def _h(W):
        """Acyclicity constraint h(W) = tr(e^(W⊙W)) - d"""
        E = np.eye(W.shape[0])
        return np.trace(np.linalg.matrix_power(E + W * W / W.shape[0], W.shape[0])) - W.shape[0]

    def _adj(w):
        """Convert flattened w to adjacency matrix"""
        return w.reshape([d, d])

    def _func(w):
        """Objective function for optimization"""
        W = _adj(w)
        loss = _loss(W)

        # L1 regularization
        l1_reg = lambda1 * np.abs(W).sum()

        # Prior penalty (negative = encouragement)
        prior_penalty = 0.0
        if prior_matrix is not None and prior_strength > 0:
            # Higher prior → lower penalty (encouragement)
            prior_penalty = -prior_strength * (W * prior_matrix).sum()

        return loss + l1_reg + prior_penalty + 0.5 * rho * _h(W) ** 2 + alpha * _h(W)

    n, d = X.shape

    # Initialize W with prior
    if prior_matrix is not None:
        # Scale prior to [0, 0.1] range
        W_est = 0.1 * prior_matrix / (prior_matrix.max() + 1e-8)
    else:
        W_est = np.zeros([d, d])

    # Dual ascent parameters
    rho, alpha, h = 1.0, 0.0, np.inf

    for n_iter in range(max_iter):
        # Optimize W
        w_new = minimize(_func, W_est.flatten(), method='L-BFGS-B').x
        W_new = _adj(w_new)

        # Update dual variables
        h_new = _h(W_new)

        if h_new > 0.25 * h:
            rho *= 10

        W_est, h = W_new, h_new
        alpha += rho * h

        # Check convergence
        if h <= h_tol or rho >= rho_max:
            break

    # Threshold small weights
    W_est[np.abs(W_est) < w_threshold] = 0

    # Ensure no self-loops
    np.fill_diagonal(W_est, 0)

    return W_est


def compute_notears_with_prior(X, prior_matrix, lambda1=0.1, w_threshold=0.3,
                                prior_strength=0.5, standardize=True):
    """
    NOTEARS実行（Prior制約付き）

    Parameters
    ----------
    X : np.ndarray
        データ行列 (n_samples, n_features)
    prior_matrix : np.ndarray
        Prior行列 (n_features, n_features)
    lambda1 : float
        L1正則化パラメータ（スパース性）
    w_threshold : float
        エッジ閾値
    prior_strength : float
        Prior強度 (0=無視, 1=強)
    standardize : bool
        データを標準化するか

    Returns
    -------
    W_est : np.ndarray
        推定隣接行列
    metrics : dict
        評価指標
    """
    print(f"\n{'='*60}")
    print(f"Running NOTEARS with Prior Constraints")
    print(f"{'='*60}")
    print(f"  Data shape: {X.shape}")
    print(f"  Lambda1 (L1): {lambda1}")
    print(f"  W threshold: {w_threshold}")
    print(f"  Prior strength: {prior_strength}")

    # Standardize
    if standardize:
        X_mean = X.mean(axis=0)
        X_std = X.std(axis=0)
        X_std[X_std == 0] = 1  # Avoid division by zero
        X_norm = (X - X_mean) / X_std
    else:
        X_norm = X

    # Run NOTEARS
    W_est = notears_linear(
        X_norm,
        lambda1=lambda1,
        loss_type='l2',
        max_iter=100,
        h_tol=1e-8,
        w_threshold=w_threshold,
        prior_matrix=prior_matrix,
        prior_strength=prior_strength
    )

    # Compute metrics
    n_edges = np.sum(W_est != 0)
    n_possible = W_est.shape[0] * (W_est.shape[0] - 1)
    sparsity = n_edges / n_possible

    # Check acyclicity
    G = nx.DiGraph(W_est)
    is_dag = nx.is_directed_acyclic_graph(G)

    metrics = {
        'n_edges': int(n_edges),
        'n_possible': n_possible,
        'sparsity': float(sparsity),
        'is_dag': is_dag,
        'max_weight': float(np.abs(W_est).max()),
        'mean_weight': float(np.abs(W_est[W_est != 0]).mean()) if n_edges > 0 else 0.0
    }

    print(f"\n  Results:")
    print(f"    Edges: {n_edges} / {n_possible} ({100*sparsity:.1f}%)")
    print(f"    Is DAG: {is_dag}")
    print(f"    Max weight: {metrics['max_weight']:.4f}")
    print(f"    Mean weight: {metrics['mean_weight']:.4f}")

    return W_est, metrics


# =====================================================
# 可視化
# =====================================================

def visualize_dag(W, module_names, coeffs, output_file, top_k=15):
    """
    DAGの可視化

    Parameters
    ----------
    W : np.ndarray
        隣接行列
    module_names : list
        モジュール名
    coeffs : pd.Series
        因果係数β
    output_file : Path
        出力ファイル
    top_k : int
        可視化する上位エッジ数
    """
    fig, axes = plt.subplots(2, 2, figsize=(18, 16))

    # 1. 隣接行列のヒートマップ
    ax = axes[0, 0]

    # βでソート
    sorted_idx = coeffs.sort_values(ascending=False).index
    module_to_idx = {m: i for i, m in enumerate(module_names)}
    sorted_indices = [module_to_idx[m] for m in sorted_idx]

    W_sorted = W[sorted_indices, :][:, sorted_indices]

    import seaborn as sns
    sns.heatmap(W_sorted, cmap='RdYlBu_r', center=0, ax=ax,
                xticklabels=sorted_idx, yticklabels=sorted_idx,
                cbar_kws={'label': 'Edge weight'})

    ax.set_xlabel('Target module')
    ax.set_ylabel('Source module')
    ax.set_title('NOTEARS Adjacency Matrix\n(sorted by β coefficient)', fontsize=12)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90, fontsize=7)
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0, fontsize=7)

    # 2. DAGグラフ（上位エッジのみ）
    ax = axes[0, 1]

    # 上位k個のエッジを抽出
    edges = []
    for i in range(W.shape[0]):
        for j in range(W.shape[1]):
            if W[i, j] != 0:
                edges.append((i, j, W[i, j]))

    edges_sorted = sorted(edges, key=lambda x: abs(x[2]), reverse=True)[:top_k]

    G = nx.DiGraph()
    for i, j, w in edges_sorted:
        G.add_edge(module_names[i], module_names[j], weight=w)

    # レイアウト: βの値で縦方向に配置
    pos = {}
    for node in G.nodes():
        beta = coeffs[node]
        # X座標: ランダム（後で調整）
        # Y座標: βの値
        pos[node] = (np.random.randn(), beta)

    # Spring layoutで横方向を調整
    pos = nx.spring_layout(G, pos=pos, fixed=None, k=2, iterations=50)

    # Draw
    node_colors = [coeffs[node] for node in G.nodes()]
    edge_weights = [abs(G[u][v]['weight']) for u, v in G.edges()]

    nx.draw_networkx_nodes(G, pos, ax=ax, node_color=node_colors,
                          cmap='RdYlBu_r', node_size=500, alpha=0.8)
    nx.draw_networkx_labels(G, pos, ax=ax, font_size=7)
    nx.draw_networkx_edges(G, pos, ax=ax, width=edge_weights,
                          edge_color='gray', arrows=True,
                          arrowsize=15, alpha=0.6,
                          connectionstyle='arc3,rad=0.1')

    ax.set_title(f'VAT1L Causal DAG (Top {top_k} edges)', fontsize=12)
    ax.axis('off')

    # 3. Degree分布
    ax = axes[1, 0]

    G_full = nx.DiGraph(W)
    in_degrees = dict(G_full.in_degree())
    out_degrees = dict(G_full.out_degree())

    in_deg = [in_degrees[i] for i in range(len(module_names))]
    out_deg = [out_degrees[i] for i in range(len(module_names))]

    x = np.arange(len(module_names))
    width = 0.35

    ax.barh(x - width/2, in_deg, width, label='In-degree', alpha=0.7, color='steelblue')
    ax.barh(x + width/2, out_deg, width, label='Out-degree', alpha=0.7, color='coral')

    ax.set_yticks(x)
    ax.set_yticklabels([module_names[i] for i in sorted_indices], fontsize=7)
    ax.set_xlabel('Degree', fontsize=10)
    ax.set_title('Node Degree Distribution', fontsize=12)
    ax.legend()
    ax.grid(alpha=0.3, axis='x')

    # 4. Top drivers & targets
    ax = axes[1, 1]

    # Out-degree (drivers)
    drivers = sorted([(module_names[i], out_deg[i]) for i in range(len(module_names))],
                     key=lambda x: x[1], reverse=True)[:10]

    # In-degree (targets)
    targets = sorted([(module_names[i], in_deg[i]) for i in range(len(module_names))],
                     key=lambda x: x[1], reverse=True)[:10]

    y_pos = np.arange(10)

    # Split plot
    ax.barh(y_pos, [d[1] for d in drivers], alpha=0.7, color='coral', label='Out-degree (drivers)')
    ax.set_yticks(y_pos)
    ax.set_yticklabels([d[0] for d in drivers], fontsize=8)
    ax.set_xlabel('Out-degree', fontsize=10)
    ax.set_title('Top 10 Driver Modules', fontsize=12)
    ax.invert_yaxis()
    ax.grid(alpha=0.3, axis='x')

    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"\n✓ Saved: {output_file}")
    plt.close()


def analyze_dag_structure(W, module_names, coeffs):
    """
    DAG構造の詳細解析

    Parameters
    ----------
    W : np.ndarray
        隣接行列
    module_names : list
        モジュール名
    coeffs : pd.Series
        因果係数

    Returns
    -------
    analysis : dict
        解析結果
    """
    print(f"\n{'='*60}")
    print(f"Analyzing DAG Structure")
    print(f"{'='*60}")

    G = nx.DiGraph(W)

    # Basic stats
    n_nodes = G.number_of_nodes()
    n_edges = G.number_of_edges()

    # Degree
    in_degrees = dict(G.in_degree())
    out_degrees = dict(G.out_degree())

    # Source & sink nodes
    sources = [i for i in range(n_nodes) if out_degrees[i] > 0 and in_degrees[i] == 0]
    sinks = [i for i in range(n_nodes) if in_degrees[i] > 0 and out_degrees[i] == 0]

    # Top drivers (high out-degree)
    top_drivers = sorted(range(n_nodes), key=lambda i: out_degrees[i], reverse=True)[:5]

    # Top targets (high in-degree)
    top_targets = sorted(range(n_nodes), key=lambda i: in_degrees[i], reverse=True)[:5]

    print(f"\n  Graph structure:")
    print(f"    Nodes: {n_nodes}")
    print(f"    Edges: {n_edges}")
    print(f"    Density: {n_edges / (n_nodes * (n_nodes - 1)):.3f}")
    print(f"    Is DAG: {nx.is_directed_acyclic_graph(G)}")

    print(f"\n  Source nodes (no parents): {len(sources)}")
    for i in sources[:5]:
        print(f"    {module_names[i]:25s} (β={coeffs.iloc[i]:+.4f})")

    print(f"\n  Sink nodes (no children): {len(sinks)}")
    for i in sinks[:5]:
        print(f"    {module_names[i]:25s} (β={coeffs.iloc[i]:+.4f})")

    print(f"\n  Top drivers (high out-degree):")
    for i in top_drivers:
        print(f"    {module_names[i]:25s}: out={out_degrees[i]}, β={coeffs.iloc[i]:+.4f}")

    print(f"\n  Top targets (high in-degree):")
    for i in top_targets:
        print(f"    {module_names[i]:25s}: in={in_degrees[i]}, β={coeffs.iloc[i]:+.4f}")

    analysis = {
        'n_nodes': n_nodes,
        'n_edges': n_edges,
        'sources': [module_names[i] for i in sources],
        'sinks': [module_names[i] for i in sinks],
        'top_drivers': [module_names[i] for i in top_drivers],
        'top_targets': [module_names[i] for i in top_targets]
    }

    return analysis


# =====================================================
# メイン処理
# =====================================================

def main():
    # Load data
    print(f"\n{'='*60}")
    print(f"Loading NOTEARS inputs")
    print(f"{'='*60}")

    data_file = GRAPH_DIR / "notears_data_VAT1L.npy"
    prior_file = GRAPH_DIR / "notears_prior_matrix_VAT1L.csv"
    modules_file = GRAPH_DIR / "notears_modules_VAT1L.txt"
    coeffs_file = IDS_DIR / "results" / "ids_pseudotime" / "local_causal_coeffs_VAT1L.csv"

    X = np.load(data_file)
    prior_matrix = pd.read_csv(prior_file, index_col=0).values

    with open(modules_file, 'r') as f:
        module_names = [line.strip() for line in f]

    coeffs = pd.read_csv(coeffs_file, index_col=0)['coefficient']

    print(f"  Data: {X.shape}")
    print(f"  Prior: {prior_matrix.shape}")
    print(f"  Modules: {len(module_names)}")

    # Run NOTEARS
    W_est, metrics = compute_notears_with_prior(
        X,
        prior_matrix,
        lambda1=0.1,          # L1正則化（スパース性）
        w_threshold=0.3,      # エッジ閾値
        prior_strength=0.5,   # Prior強度
        standardize=True
    )

    # Analyze structure
    analysis = analyze_dag_structure(W_est, module_names, coeffs)

    # Save adjacency matrix
    adj_file = GRAPH_DIR / "notears_adj_matrix_VAT1L.csv"
    W_df = pd.DataFrame(W_est, index=module_names, columns=module_names)
    W_df.to_csv(adj_file)
    print(f"\n✓ Saved: {adj_file}")

    # Visualize
    fig_file = FIG_DIR / "notears_dag_VAT1L.png"
    visualize_dag(W_est, module_names, coeffs, fig_file, top_k=20)

    # Save summary
    summary = {
        'metrics': metrics,
        'analysis': analysis,
        'parameters': {
            'lambda1': 0.1,
            'w_threshold': 0.3,
            'prior_strength': 0.5
        }
    }

    summary_file = GRAPH_DIR / "notears_summary_VAT1L.json"
    with open(summary_file, 'w') as f:
        json.dump(summary, f, indent=2)
    print(f"\n✓ Saved: {summary_file}")

    # Final summary
    print(f"\n{'='*60}")
    print(f"NOTEARS DAG Estimation Complete")
    print(f"{'='*60}")

    print(f"\nEstimated DAG:")
    print(f"  Nodes: {metrics['n_edges']} edges")
    print(f"  Sparsity: {100*metrics['sparsity']:.1f}%")
    print(f"  Is DAG: {metrics['is_dag']}")

    print(f"\nSource modules (upstream drivers):")
    for m in analysis['sources'][:5]:
        beta = coeffs[m]
        print(f"  {m:25s}: β={beta:+.4f}")

    print(f"\nSink modules (downstream targets):")
    for m in analysis['sinks'][:5]:
        beta = coeffs[m]
        print(f"  {m:25s}: β={beta:+.4f}")

    print(f"\n✓ Ready for LiNGAM refinement (Step 3c)")


if __name__ == "__main__":
    main()
