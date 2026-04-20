#!/usr/bin/env python3
"""
Phase 3a: NOTEARS用Prior制約行列の準備

Step 1で得られた因果係数βを用いて、NOTEARSに渡す
有向制約（prior）を生成する。

理論:
    β_i > 0 → モジュールiはストレス駆動（上流候補）
    β_i < 0 → モジュールiはストレス保護（下流候補）

    Prior edge weight:
    W_prior[i,j] ∝ sign(β_i - β_j) if β_i > β_j
                  = 0 otherwise

Input:
- results/ids_pseudotime/local_causal_coeffs_VAT1L.csv
- results/ids_pseudotime/module_scores_Ex.L5.VAT1L_EYA4.csv

Output:
- results/causal_graphs/notears_prior_matrix_VAT1L.csv
- results/causal_graphs/prior_graph_structure_VAT1L.png

Usage:
    python 03a_prepare_notears_prior.py
"""

import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

# =====================================================
# パス設定
# =====================================================

BASE_DIR = Path("/home/akaco/als/motor_cortex_analysis")
IDS_DIR = BASE_DIR / "ids_causal_analysis"
RESULTS_DIR = IDS_DIR / "results" / "ids_pseudotime"
GRAPH_DIR = IDS_DIR / "results" / "causal_graphs"
FIG_DIR = IDS_DIR / "results" / "figures"

GRAPH_DIR.mkdir(parents=True, exist_ok=True)

# =====================================================
# Prior制約の生成
# =====================================================

def load_causal_coefficients():
    """
    Step 1で得られた因果係数を読み込む

    Returns
    -------
    coeffs : pd.Series
        因果係数 β
    """
    print(f"\n{'='*60}")
    print(f"Loading causal coefficients")
    print(f"{'='*60}")

    coeffs_file = RESULTS_DIR / "local_causal_coeffs_VAT1L.csv"

    if not coeffs_file.exists():
        raise FileNotFoundError(f"Causal coefficients not found: {coeffs_file}")

    df = pd.read_csv(coeffs_file, index_col=0)

    coeffs = df['coefficient']

    print(f"  Loaded {len(coeffs)} module coefficients")
    print(f"  Range: [{coeffs.min():.4f}, {coeffs.max():.4f}]")

    return coeffs


def create_prior_matrix(coeffs, method='gradient', strength=1.0):
    """
    因果係数から prior 制約行列を生成

    Parameters
    ----------
    coeffs : pd.Series
        因果係数 β
    method : str
        'gradient' - β差分ベース（推奨）
        'binary' - β符号のみ
        'magnitude' - β絶対値ベース
    strength : float
        Prior の強さ（0=なし, 1=標準）

    Returns
    -------
    prior_matrix : pd.DataFrame
        Prior 制約行列 (modules × modules)
        prior[i,j] > 0 → i→j の edge を優先
        prior[i,j] = 0 → edge なし or 不明
    """
    print(f"\n{'='*60}")
    print(f"Creating prior constraint matrix")
    print(f"{'='*60}")
    print(f"  Method: {method}")
    print(f"  Strength: {strength}")

    modules = coeffs.index.tolist()
    n = len(modules)

    prior = np.zeros((n, n))

    if method == 'gradient':
        # β_i > β_j なら i→j のedge候補
        # 差分が大きいほど強い prior
        for i in range(n):
            for j in range(n):
                if i == j:
                    continue

                beta_diff = coeffs.iloc[i] - coeffs.iloc[j]

                # 正の差分のみ prior として使用
                if beta_diff > 0:
                    prior[i, j] = strength * beta_diff
                # 負の差分は逆方向を禁止
                elif beta_diff < 0:
                    prior[i, j] = 0  # j→i が正しい方向

    elif method == 'binary':
        # β符号のみ使用（簡易版）
        beta_positive = coeffs > 0
        beta_negative = coeffs < 0

        for i in range(n):
            for j in range(n):
                if i == j:
                    continue

                # Positive → Negative の方向を優先
                if beta_positive.iloc[i] and beta_negative.iloc[j]:
                    prior[i, j] = strength
                # Negative → Positive は禁止
                elif beta_negative.iloc[i] and beta_positive.iloc[j]:
                    prior[i, j] = 0

    elif method == 'magnitude':
        # β絶対値ベース（大→小の方向）
        abs_coeffs = np.abs(coeffs)

        for i in range(n):
            for j in range(n):
                if i == j:
                    continue

                mag_diff = abs_coeffs.iloc[i] - abs_coeffs.iloc[j]

                if mag_diff > 0:
                    prior[i, j] = strength * mag_diff

    else:
        raise ValueError(f"Unknown method: {method}")

    # DataFrameに変換
    prior_df = pd.DataFrame(prior, index=modules, columns=modules)

    # 統計
    n_nonzero = np.sum(prior > 0)
    n_total = n * (n - 1)  # 対角線除く

    print(f"\n  Prior statistics:")
    print(f"    Non-zero edges: {n_nonzero} / {n_total} ({100*n_nonzero/n_total:.1f}%)")
    print(f"    Mean prior strength: {prior[prior > 0].mean():.4f}")
    print(f"    Max prior strength: {prior.max():.4f}")

    return prior_df


def identify_causal_hierarchy(coeffs, n_levels=5):
    """
    因果係数に基づいて階層構造を同定

    Parameters
    ----------
    coeffs : pd.Series
        因果係数
    n_levels : int
        階層レベル数

    Returns
    -------
    hierarchy : pd.DataFrame
        階層情報
    """
    print(f"\n{'='*60}")
    print(f"Identifying causal hierarchy")
    print(f"{'='*60}")

    # βでソート
    sorted_coeffs = coeffs.sort_values(ascending=False)

    # 階層に分割（quantileベース）
    hierarchy = pd.DataFrame({
        'module': sorted_coeffs.index,
        'coefficient': sorted_coeffs.values,
        'rank': range(1, len(sorted_coeffs) + 1)
    })

    # レベル割り当て
    hierarchy['level'] = pd.qcut(
        hierarchy['rank'],
        q=n_levels,
        labels=list(range(1, n_levels + 1))
    )

    print(f"\n  Hierarchy levels:")
    for level in range(1, n_levels + 1):
        modules_in_level = hierarchy[hierarchy['level'] == level]['module'].tolist()
        print(f"    Level {level}: {len(modules_in_level)} modules")
        print(f"      {', '.join(modules_in_level[:3])}, ...")

    return hierarchy


def visualize_prior_structure(prior_matrix, coeffs, hierarchy, output_file):
    """
    Prior構造の可視化

    Parameters
    ----------
    prior_matrix : pd.DataFrame
        Prior行列
    coeffs : pd.Series
        因果係数
    hierarchy : pd.DataFrame
        階層情報
    output_file : Path
        出力ファイル
    """
    fig, axes = plt.subplots(2, 2, figsize=(16, 14))

    # 1. Prior行列のヒートマップ
    ax = axes[0, 0]

    # 係数でソート
    sorted_modules = coeffs.sort_values(ascending=False).index
    prior_sorted = prior_matrix.loc[sorted_modules, sorted_modules]

    sns.heatmap(prior_sorted, cmap='YlOrRd', ax=ax,
                xticklabels=True, yticklabels=True,
                cbar_kws={'label': 'Prior strength'})

    ax.set_xlabel('Target module')
    ax.set_ylabel('Source module')
    ax.set_title('NOTEARS Prior Matrix\n(sorted by β coefficient)', fontsize=12)

    # ラベルを斜めに
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90, fontsize=6)
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0, fontsize=6)

    # 2. Prior の強さ分布
    ax = axes[0, 1]

    prior_values = prior_matrix.values.flatten()
    prior_nonzero = prior_values[prior_values > 0]

    ax.hist(prior_nonzero, bins=50, edgecolor='black', alpha=0.7)
    ax.set_xlabel('Prior strength', fontsize=10)
    ax.set_ylabel('Frequency', fontsize=10)
    ax.set_title(f'Distribution of Prior Strengths\n(n={len(prior_nonzero)} edges)', fontsize=12)
    ax.grid(alpha=0.3)

    # 3. 階層構造（係数ベース）
    ax = axes[1, 0]

    hierarchy_sorted = hierarchy.sort_values('coefficient', ascending=False)
    colors = plt.cm.RdYlBu_r(hierarchy_sorted['level'].astype(int) / hierarchy['level'].max())

    bars = ax.barh(range(len(hierarchy_sorted)), hierarchy_sorted['coefficient'],
                   color=colors, alpha=0.8)

    ax.set_yticks(range(len(hierarchy_sorted)))
    ax.set_yticklabels(hierarchy_sorted['module'], fontsize=7)
    ax.axvline(0, color='black', linestyle='--', linewidth=1)
    ax.set_xlabel('Causal coefficient β', fontsize=10)
    ax.set_title('Causal Hierarchy (by β coefficient)', fontsize=12)
    ax.invert_yaxis()
    ax.grid(alpha=0.3, axis='x')

    # 4. レベル別の統計
    ax = axes[1, 1]

    level_stats = hierarchy.groupby('level')['coefficient'].agg(['mean', 'std', 'count'])
    level_stats.index = [f'Level {i}' for i in level_stats.index]

    x = range(len(level_stats))
    ax.bar(x, level_stats['mean'], yerr=level_stats['std'],
           alpha=0.7, capsize=5, color='steelblue')

    ax.set_xticks(x)
    ax.set_xticklabels(level_stats.index, rotation=45)
    ax.axhline(0, color='black', linestyle='--', linewidth=1)
    ax.set_ylabel('Mean coefficient β', fontsize=10)
    ax.set_title('Hierarchy Level Statistics', fontsize=12)
    ax.grid(alpha=0.3, axis='y')

    # 数値表示
    for i, (mean_val, count) in enumerate(zip(level_stats['mean'], level_stats['count'])):
        ax.text(i, mean_val, f'n={count}', ha='center',
                va='bottom' if mean_val > 0 else 'top', fontsize=9)

    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"\n✓ Saved: {output_file}")
    plt.close()


def export_for_notears(prior_matrix, module_scores_file, output_dir):
    """
    NOTEARS実行用のファイルをエクスポート

    Parameters
    ----------
    prior_matrix : pd.DataFrame
        Prior行列
    module_scores_file : Path
        モジュールスコアファイル
    output_dir : Path
        出力ディレクトリ

    Returns
    -------
    export_files : dict
        エクスポートしたファイルのパス
    """
    print(f"\n{'='*60}")
    print(f"Exporting for NOTEARS")
    print(f"{'='*60}")

    export_files = {}

    # 1. Prior行列
    prior_file = output_dir / "notears_prior_matrix_VAT1L.csv"
    prior_matrix.to_csv(prior_file)
    export_files['prior'] = prior_file
    print(f"  Prior matrix: {prior_file}")

    # 2. モジュールスコア（numpy形式）
    scores_df = pd.read_csv(module_scores_file, index_col=0)

    # CRITICAL: モジュール順序を prior_matrix と一致させる
    # prior_matrix の index/columns が正しい順序
    modules = prior_matrix.index.tolist()

    # この順序でデータを抽出
    X = scores_df[modules].values

    # numpy保存
    data_file = output_dir / "notears_data_VAT1L.npy"
    np.save(data_file, X)
    export_files['data'] = data_file
    print(f"  Data matrix: {data_file} (shape: {X.shape})")

    # モジュール名保存
    modules_file = output_dir / "notears_modules_VAT1L.txt"
    with open(modules_file, 'w') as f:
        for m in modules:
            f.write(f"{m}\n")
    export_files['modules'] = modules_file
    print(f"  Module names: {modules_file}")

    # 3. メタデータ
    import json
    meta = {
        'n_samples': X.shape[0],
        'n_modules': X.shape[1],
        'modules': modules,
        'prior_edges': int(np.sum(prior_matrix.values > 0)),
        'cell_type': 'Ex.L5.VAT1L_EYA4'
    }

    meta_file = output_dir / "notears_metadata_VAT1L.json"
    with open(meta_file, 'w') as f:
        json.dump(meta, f, indent=2)
    export_files['metadata'] = meta_file
    print(f"  Metadata: {meta_file}")

    return export_files


# =====================================================
# メイン処理
# =====================================================

def main():
    # 因果係数読み込み
    coeffs = load_causal_coefficients()

    # Prior行列生成
    prior_matrix = create_prior_matrix(
        coeffs,
        method='gradient',  # β差分ベース
        strength=1.0
    )

    # 階層同定
    hierarchy = identify_causal_hierarchy(coeffs, n_levels=5)

    # 可視化
    fig_file = FIG_DIR / "prior_graph_structure_VAT1L.png"
    visualize_prior_structure(prior_matrix, coeffs, hierarchy, fig_file)

    # NOTEARS用にエクスポート
    module_scores_file = RESULTS_DIR / "module_scores_Ex.L5.VAT1L_EYA4.csv"
    export_files = export_for_notears(prior_matrix, module_scores_file, GRAPH_DIR)

    # サマリー
    print(f"\n{'='*60}")
    print(f"NOTEARS Prior Preparation Complete")
    print(f"{'='*60}")

    print(f"\nPrior matrix:")
    print(f"  Shape: {prior_matrix.shape}")
    print(f"  Non-zero edges: {np.sum(prior_matrix.values > 0)}")
    print(f"  Mean strength: {prior_matrix.values[prior_matrix.values > 0].mean():.4f}")

    print(f"\nCausal hierarchy:")
    print(f"  Levels: {hierarchy['level'].nunique()}")
    print(f"  Top level (strongest drivers):")
    top_level = hierarchy[hierarchy['level'] == 1]
    for _, row in top_level.iterrows():
        print(f"    {row['module']:25s}: β={row['coefficient']:+.4f}")

    print(f"\nExported files:")
    for key, path in export_files.items():
        print(f"  {key}: {path}")

    print(f"\n✓ Ready for NOTEARS execution (Step 3b)")


if __name__ == "__main__":
    main()
