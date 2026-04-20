#!/usr/bin/env python3
"""
Phase 1: 複合ストレスポテンシャルの計算

既存の発現データから、ALS特化型ストレス成分を計算し、
各細胞のストレススコアを生成する。

Input:
- motor_cortex_{cell_type}_expression.csv(.gz)
- motor_cortex_{cell_type}_metadata.csv(.gz)

Output:
- results/stress_components/stress_scores_{cell_type}.csv
- results/stress_components/stress_components_{cell_type}.csv

Usage:
    python 01_composite_stress.py --cell_type Ex.L5.VAT1L_EYA4
    python 01_composite_stress.py --all  # 全細胞型
"""

import os
import sys
import argparse
import numpy as np
import pandas as pd
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from config.als_stress_config import (
    compute_als_stress,
    ALS_STRESS_COMPONENTS,
    plot_stress_components
)

# =====================================================
# パス設定
# =====================================================

BASE_DIR = Path("/home/akaco/als/motor_cortex_analysis")
IDS_DIR = BASE_DIR / "ids_causal_analysis"
RESULTS_DIR = IDS_DIR / "results" / "stress_components"

RESULTS_DIR.mkdir(parents=True, exist_ok=True)

# =====================================================
# ユーティリティ関数
# =====================================================

def load_expression_data(cell_type):
    """
    発現データとメタデータを読み込む（ENSG→Symbol変換付き）

    Parameters
    ----------
    cell_type : str
        細胞型名（例: "Ex.L5.VAT1L_EYA4"）

    Returns
    -------
    expr : pd.DataFrame
        発現データ (genes × cells), gene names as Symbol
    meta : pd.DataFrame
        メタデータ
    """
    print(f"\n{'='*60}")
    print(f"Loading data for: {cell_type}")
    print(f"{'='*60}")

    # ファイルパス
    expr_file_csv = BASE_DIR / f"motor_cortex_{cell_type}_expression.csv"
    expr_file_gz = BASE_DIR / f"motor_cortex_{cell_type}_expression.csv.gz"
    meta_file_csv = BASE_DIR / f"motor_cortex_{cell_type}_metadata.csv"
    meta_file_gz = BASE_DIR / f"motor_cortex_{cell_type}_metadata.csv.gz"

    # 発現データ読み込み
    if expr_file_csv.exists():
        expr_file = expr_file_csv
    elif expr_file_gz.exists():
        expr_file = expr_file_gz
    else:
        raise FileNotFoundError(f"Expression file not found for {cell_type}")

    print(f"  Reading expression: {expr_file.name}")
    expr = pd.read_csv(expr_file, index_col=0, compression='infer')

    # 転置チェック（遺伝子×細胞 にする）
    if expr.shape[0] < expr.shape[1]:
        print(f"  Transposing: {expr.shape} → ", end='')
        expr = expr.T
        print(f"{expr.shape}")

    print(f"  Expression shape: {expr.shape} (genes × cells)")

    # ENSG→Symbol変換
    print(f"  Converting ENSG IDs to gene symbols...")
    gene_mapping_file = BASE_DIR / "gene_mapping.json"

    if gene_mapping_file.exists():
        import json
        with open(gene_mapping_file, 'r') as f:
            gene_mapping_raw = json.load(f)

        # ensg_to_symbol辞書を作成
        if 'ensg_to_symbol' in gene_mapping_raw:
            ensg_to_symbol = gene_mapping_raw['ensg_to_symbol']
        else:
            # 旧形式の場合
            ensg_to_symbol = {}
            for k, v in gene_mapping_raw.items():
                if isinstance(v, str):
                    ensg_to_symbol[k] = v
                elif isinstance(v, dict) and 'symbol' in v:
                    ensg_to_symbol[k] = v['symbol']

        # 遺伝子名を変換
        new_index = []
        for ensg_id in expr.index:
            if ensg_id in ensg_to_symbol:
                new_index.append(ensg_to_symbol[ensg_id])
            else:
                new_index.append(ensg_id)  # 変換できない場合はそのまま

        expr.index = new_index
        n_converted = sum([1 for idx in expr.index if not idx.startswith('ENSG')])
        print(f"  Converted {n_converted} / {len(expr.index)} genes to symbols")
    else:
        print(f"  Warning: gene_mapping.json not found, using ENSG IDs")

    print(f"  Final expression shape: {expr.shape} (genes × cells)")

    # メタデータ読み込み
    if meta_file_csv.exists():
        meta_file = meta_file_csv
    elif meta_file_gz.exists():
        meta_file = meta_file_gz
    else:
        raise FileNotFoundError(f"Metadata file not found for {cell_type}")

    print(f"  Reading metadata: {meta_file.name}")
    meta = pd.read_csv(meta_file, compression='infer')

    print(f"  Metadata shape: {meta.shape}")

    # 列名の正規化
    if 'cell_id' not in meta.columns:
        meta.rename(columns={meta.columns[0]: 'cell_id'}, inplace=True)

    if 'condition' not in meta.columns and 'Disease' in meta.columns:
        meta.rename(columns={'Disease': 'condition'}, inplace=True)

    # 共通細胞の抽出
    common_cells = list(set(expr.columns) & set(meta['cell_id']))
    print(f"  Common cells: {len(common_cells)}")

    if len(common_cells) == 0:
        raise ValueError("No common cells found between expression and metadata!")

    # データフィルタリング
    expr = expr[common_cells]
    meta = meta[meta['cell_id'].isin(common_cells)].set_index('cell_id')
    meta = meta.loc[common_cells]  # 順序を揃える

    # Condition分布
    if 'condition' in meta.columns:
        print(f"\n  Condition distribution:")
        for cond, count in meta['condition'].value_counts().items():
            print(f"    {cond}: {count}")

    return expr, meta


def compute_stress_for_celltype(cell_type, save_plots=True):
    """
    特定の細胞型のストレススコアを計算

    Parameters
    ----------
    cell_type : str
        細胞型名
    save_plots : bool
        可視化を保存するか

    Returns
    -------
    results_df : pd.DataFrame
        ストレススコア（cell_id × components）
    """
    # データ読み込み
    expr, meta = load_expression_data(cell_type)

    # 遺伝子×細胞 → 細胞×遺伝子に転置
    X = expr.T.values
    gene_names = expr.index.tolist()
    cell_ids = expr.columns.tolist()

    print(f"\nComputing ALS stress components...")
    print(f"  Input shape: {X.shape} (cells × genes)")

    # ストレス計算
    S_total, S_components = compute_als_stress(
        X, gene_names,
        config=ALS_STRESS_COMPONENTS,
        return_components=True
    )

    print(f"\n  Total stress range: [{S_total.min():.3f}, {S_total.max():.3f}]")
    print(f"  Total stress mean: {S_total.mean():.3f}")

    # 結果をDataFrame化
    results_df = pd.DataFrame({
        'cell_id': cell_ids,
        'stress_total': S_total
    })

    # 成分ごとのスコアを追加
    for comp_name, comp_values in S_components.items():
        results_df[f'stress_{comp_name}'] = comp_values

    # メタデータを結合
    results_df = results_df.merge(
        meta[['condition']].reset_index(),
        on='cell_id',
        how='left'
    )

    # 保存
    output_file = RESULTS_DIR / f"stress_scores_{cell_type}.csv"
    results_df.to_csv(output_file, index=False)
    print(f"\n✓ Saved: {output_file}")

    # 成分別の保存
    components_df = pd.DataFrame(S_components, index=cell_ids)
    components_df.index.name = 'cell_id'
    comp_output_file = RESULTS_DIR / f"stress_components_{cell_type}.csv"
    components_df.to_csv(comp_output_file)
    print(f"✓ Saved: {comp_output_file}")

    # 可視化
    if save_plots:
        try:
            import matplotlib
            matplotlib.use('Agg')  # Non-interactive backend
            import matplotlib.pyplot as plt
            import seaborn as sns

            # Figure 1: ストレス成分の分布
            fig, axes = plt.subplots(2, 3, figsize=(18, 10))
            axes = axes.flatten()

            for i, (comp_name, comp_values) in enumerate(S_components.items()):
                if i >= 6:
                    break

                ax = axes[i]

                # Condition別のヒストグラム
                for cond in meta['condition'].unique():
                    mask = meta['condition'] == cond
                    ax.hist(comp_values[mask], bins=50, alpha=0.5,
                           label=cond, density=True)

                ax.set_xlabel('Stress Score')
                ax.set_ylabel('Density')
                ax.set_title(f'{comp_name}')
                ax.legend()
                ax.grid(alpha=0.3)

            plt.tight_layout()
            fig_file = RESULTS_DIR / f"stress_distributions_{cell_type}.png"
            plt.savefig(fig_file, dpi=300, bbox_inches='tight')
            plt.close()
            print(f"✓ Saved: {fig_file}")

            # Figure 2: 総合ストレスの比較
            fig, ax = plt.subplots(figsize=(10, 6))

            for cond in meta['condition'].unique():
                mask = meta['condition'] == cond
                ax.hist(S_total[mask], bins=50, alpha=0.5,
                       label=cond, density=True)

            ax.set_xlabel('Total Stress Score')
            ax.set_ylabel('Density')
            ax.set_title(f'Total ALS Stress: {cell_type}')
            ax.legend()
            ax.grid(alpha=0.3)

            plt.tight_layout()
            fig_file2 = RESULTS_DIR / f"stress_total_{cell_type}.png"
            plt.savefig(fig_file2, dpi=300, bbox_inches='tight')
            plt.close()
            print(f"✓ Saved: {fig_file2}")

            # Figure 3: 成分間相関
            fig, ax = plt.subplots(figsize=(10, 8))

            corr_matrix = components_df.corr()
            sns.heatmap(corr_matrix, annot=True, fmt='.2f',
                       cmap='coolwarm', center=0, ax=ax,
                       cbar_kws={'label': 'Correlation'})

            ax.set_title(f'Stress Component Correlations: {cell_type}')
            plt.tight_layout()
            fig_file3 = RESULTS_DIR / f"stress_correlations_{cell_type}.png"
            plt.savefig(fig_file3, dpi=300, bbox_inches='tight')
            plt.close()
            print(f"✓ Saved: {fig_file3}")

        except Exception as e:
            print(f"Warning: Visualization failed: {e}")

    # 統計サマリー
    print(f"\n{'='*60}")
    print(f"Stress Statistics Summary")
    print(f"{'='*60}")

    if 'condition' in meta.columns:
        for cond in meta['condition'].unique():
            mask = meta['condition'] == cond
            print(f"\n{cond}:")
            print(f"  Total stress: {S_total[mask].mean():.3f} ± {S_total[mask].std():.3f}")

            for comp_name in S_components.keys():
                mean_val = S_components[comp_name][mask].mean()
                print(f"  {comp_name:20s}: {mean_val:.3f}")

    return results_df


# =====================================================
# メイン処理
# =====================================================

def main():
    parser = argparse.ArgumentParser(
        description='Compute composite stress potential for ALS data'
    )
    parser.add_argument(
        '--cell_type',
        type=str,
        help='Cell type to process (e.g., Ex.L5.VAT1L_EYA4)'
    )
    parser.add_argument(
        '--all',
        action='store_true',
        help='Process all cell types'
    )
    parser.add_argument(
        '--no_plots',
        action='store_true',
        help='Skip visualization'
    )

    args = parser.parse_args()

    # 細胞型リスト
    TARGET_CELLTYPES = [
        "Ex.L5.VAT1L_EYA4",
        "Ex.L5.VAT1L_THSD4",
        "Ex.L2_L3.CUX2_RASGRF2",
        "Glia.Astro.GFAP-pos",
        "Glia.Micro",
        "In.PV.PVALB_PTHLH"
    ]

    if args.all:
        cell_types = TARGET_CELLTYPES
    elif args.cell_type:
        cell_types = [args.cell_type]
    else:
        print("Error: Specify --cell_type or --all")
        return

    # 処理
    print(f"\nProcessing {len(cell_types)} cell type(s)...")

    results = {}
    for ct in cell_types:
        try:
            results[ct] = compute_stress_for_celltype(
                ct,
                save_plots=not args.no_plots
            )
        except Exception as e:
            print(f"\n✗ Failed for {ct}: {e}")
            import traceback
            traceback.print_exc()
            continue

    # 全体サマリー
    print(f"\n{'='*60}")
    print(f"Summary: Processed {len(results)} / {len(cell_types)} cell types")
    print(f"{'='*60}")

    for ct in results.keys():
        n_cells = len(results[ct])
        print(f"  {ct}: {n_cells} cells")

    print(f"\nResults saved to: {RESULTS_DIR}")


if __name__ == "__main__":
    main()
