#!/usr/bin/env python3
"""
Phase 2a: 24モジュールスコアの計算

既存の24モジュール定義を使用して、各細胞のモジュールスコアを計算する。
これはストレスポテンシャルと統合され、局所因果ベクトル場推定の基盤となる。

Input:
- motor_cortex_{cell_type}_expression.csv(.gz)
- 24_functional_modules_fixed.json

Output:
- results/ids_pseudotime/module_scores_{cell_type}.csv

Usage:
    python 02a_compute_module_scores.py --cell_type Ex.L5.VAT1L_EYA4
"""

import os
import sys
import json
import argparse
import numpy as np
import pandas as pd
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

# =====================================================
# パス設定
# =====================================================

BASE_DIR = Path("/home/akaco/als/motor_cortex_analysis")
IDS_DIR = BASE_DIR / "ids_causal_analysis"
RESULTS_DIR = IDS_DIR / "results" / "ids_pseudotime"

RESULTS_DIR.mkdir(parents=True, exist_ok=True)

# =====================================================
# モジュールスコア計算
# =====================================================

def load_module_definitions():
    """
    既存の24モジュール定義を読み込む

    Returns
    -------
    modules : dict
        {module_name: [gene_symbols]}
    """
    module_file = BASE_DIR / "24_functional_modules_fixed.json"

    print(f"Loading module definitions from: {module_file.name}")

    with open(module_file, 'r') as f:
        data = json.load(f)

    modules = data.get('modules', data)

    print(f"  Loaded {len(modules)} modules")
    for name, genes in list(modules.items())[:5]:
        print(f"    {name}: {len(genes)} genes")

    return modules


def load_expression_with_conversion(cell_type):
    """
    発現データを読み込み、ENSG→Symbolに変換

    Parameters
    ----------
    cell_type : str
        細胞型名

    Returns
    -------
    expr : pd.DataFrame
        発現データ (genes × cells), gene names as Symbol
    meta : pd.DataFrame
        メタデータ
    """
    print(f"\n{'='*60}")
    print(f"Loading expression data for: {cell_type}")
    print(f"{'='*60}")

    # 発現ファイル
    expr_file_csv = BASE_DIR / f"motor_cortex_{cell_type}_expression.csv"
    expr_file_gz = BASE_DIR / f"motor_cortex_{cell_type}_expression.csv.gz"

    if expr_file_csv.exists():
        expr_file = expr_file_csv
    elif expr_file_gz.exists():
        expr_file = expr_file_gz
    else:
        raise FileNotFoundError(f"Expression file not found for {cell_type}")

    print(f"  Reading: {expr_file.name}")
    expr = pd.read_csv(expr_file, index_col=0, compression='infer')

    # 転置チェック
    if expr.shape[0] < expr.shape[1]:
        expr = expr.T

    print(f"  Expression shape: {expr.shape} (genes × cells)")

    # ENSG→Symbol変換
    print(f"  Converting ENSG IDs to gene symbols...")
    gene_mapping_file = BASE_DIR / "gene_mapping.json"

    if gene_mapping_file.exists():
        with open(gene_mapping_file, 'r') as f:
            gene_mapping_raw = json.load(f)

        if 'ensg_to_symbol' in gene_mapping_raw:
            ensg_to_symbol = gene_mapping_raw['ensg_to_symbol']
        else:
            ensg_to_symbol = {}
            for k, v in gene_mapping_raw.items():
                if isinstance(v, str):
                    ensg_to_symbol[k] = v
                elif isinstance(v, dict) and 'symbol' in v:
                    ensg_to_symbol[k] = v['symbol']

        # 変換
        new_index = []
        for ensg_id in expr.index:
            if ensg_id in ensg_to_symbol:
                new_index.append(ensg_to_symbol[ensg_id])
            else:
                new_index.append(ensg_id)

        expr.index = new_index
        n_converted = sum([1 for idx in expr.index if not idx.startswith('ENSG')])
        print(f"  Converted {n_converted} / {len(expr.index)} genes to symbols")

    # メタデータ
    meta_file_csv = BASE_DIR / f"motor_cortex_{cell_type}_metadata.csv"
    meta_file_gz = BASE_DIR / f"motor_cortex_{cell_type}_metadata.csv.gz"

    if meta_file_csv.exists():
        meta_file = meta_file_csv
    elif meta_file_gz.exists():
        meta_file = meta_file_gz
    else:
        raise FileNotFoundError(f"Metadata file not found for {cell_type}")

    print(f"  Reading metadata: {meta_file.name}")
    meta = pd.read_csv(meta_file, compression='infer')

    # 列名正規化
    if 'cell_id' not in meta.columns:
        meta.rename(columns={meta.columns[0]: 'cell_id'}, inplace=True)
    if 'condition' not in meta.columns and 'Disease' in meta.columns:
        meta.rename(columns={'Disease': 'condition'}, inplace=True)

    # 共通細胞
    common_cells = list(set(expr.columns) & set(meta['cell_id']))
    print(f"  Common cells: {len(common_cells)}")

    expr = expr[common_cells]
    meta = meta[meta['cell_id'].isin(common_cells)].set_index('cell_id')
    meta = meta.loc[common_cells]

    return expr, meta


def compute_module_scores(expr, modules, method='mean', center=True):
    """
    各細胞のモジュールスコアを計算

    Parameters
    ----------
    expr : pd.DataFrame
        発現データ (genes × cells)
    modules : dict
        モジュール定義 {module_name: [genes]}
    method : str
        'mean', 'median', 'pca'
    center : bool
        Control平均で中心化するか

    Returns
    -------
    module_scores : pd.DataFrame
        モジュールスコア (cells × modules)
    module_info : dict
        各モジュールの情報（使用遺伝子数など）
    """
    print(f"\nComputing module scores (method={method}, center={center})...")

    cell_ids = expr.columns.tolist()
    n_cells = len(cell_ids)

    # 遺伝子名を大文字化（マッチング用）
    expr_genes_upper = {g.upper(): g for g in expr.index}

    module_scores = pd.DataFrame(index=cell_ids)
    module_info = {}

    for module_name, module_genes in modules.items():
        # モジュール遺伝子をマッチング
        matched_genes = []
        for g in module_genes:
            g_upper = g.upper()
            if g_upper in expr_genes_upper:
                matched_genes.append(expr_genes_upper[g_upper])

        if len(matched_genes) == 0:
            print(f"  Warning: No genes found for module '{module_name}'")
            module_scores[module_name] = 0.0
            module_info[module_name] = {
                'n_genes_total': len(module_genes),
                'n_genes_found': 0,
                'genes_found': []
            }
            continue

        # モジュール発現を取得
        module_expr = expr.loc[matched_genes, :].T  # cells × genes

        # スコア計算
        if method == 'mean':
            scores = module_expr.mean(axis=1)
        elif method == 'median':
            scores = module_expr.median(axis=1)
        elif method == 'pca':
            from sklearn.decomposition import PCA
            pca = PCA(n_components=1)
            scores = pca.fit_transform(module_expr).flatten()
        else:
            raise ValueError(f"Unknown method: {method}")

        module_scores[module_name] = scores

        module_info[module_name] = {
            'n_genes_total': len(module_genes),
            'n_genes_found': len(matched_genes),
            'genes_found': matched_genes[:10]  # 最初の10個のみ保存
        }

        print(f"  {module_name:25s}: {len(matched_genes):4d}/{len(module_genes):4d} genes, "
              f"score range=[{scores.min():.2f}, {scores.max():.2f}]")

    # 中心化（オプション）
    if center:
        print(f"\n  Centering module scores...")
        module_scores = module_scores - module_scores.mean()

    return module_scores, module_info


def integrate_with_stress(module_scores, stress_file):
    """
    モジュールスコアとストレススコアを統合

    Parameters
    ----------
    module_scores : pd.DataFrame
        モジュールスコア (cells × modules)
    stress_file : Path
        ストレススコアファイル

    Returns
    -------
    integrated : pd.DataFrame
        統合データ (cells × [modules + stress_components])
    """
    print(f"\nIntegrating with stress scores...")

    if not stress_file.exists():
        print(f"  Warning: Stress file not found: {stress_file}")
        return module_scores

    stress_df = pd.read_csv(stress_file)
    stress_df = stress_df.set_index('cell_id')

    # ストレス成分カラムのみ抽出
    stress_cols = [col for col in stress_df.columns if col.startswith('stress_')]
    stress_data = stress_df[stress_cols]

    print(f"  Stress components: {len(stress_cols)}")

    # 共通細胞で統合
    common_cells = list(set(module_scores.index) & set(stress_data.index))
    print(f"  Common cells: {len(common_cells)}")

    integrated = pd.concat([
        module_scores.loc[common_cells],
        stress_data.loc[common_cells]
    ], axis=1)

    # Conditionも追加
    if 'condition' in stress_df.columns:
        integrated['condition'] = stress_df.loc[common_cells, 'condition']

    return integrated


# =====================================================
# メイン処理
# =====================================================

def main():
    parser = argparse.ArgumentParser(
        description='Compute module scores for causal inference'
    )
    parser.add_argument(
        '--cell_type',
        type=str,
        required=True,
        help='Cell type to process (e.g., Ex.L5.VAT1L_EYA4)'
    )
    parser.add_argument(
        '--method',
        type=str,
        default='mean',
        choices=['mean', 'median', 'pca'],
        help='Method to compute module scores'
    )
    parser.add_argument(
        '--no_center',
        action='store_true',
        help='Do not center module scores'
    )

    args = parser.parse_args()

    # モジュール定義読み込み
    modules = load_module_definitions()

    # 発現データ読み込み
    expr, meta = load_expression_with_conversion(args.cell_type)

    # モジュールスコア計算
    module_scores, module_info = compute_module_scores(
        expr, modules,
        method=args.method,
        center=not args.no_center
    )

    # ストレススコアと統合
    stress_file = IDS_DIR / "results" / "stress_components" / f"stress_scores_{args.cell_type}.csv"
    integrated = integrate_with_stress(module_scores, stress_file)

    # 保存
    output_file = RESULTS_DIR / f"module_scores_{args.cell_type}.csv"
    integrated.to_csv(output_file)
    print(f"\n✓ Saved: {output_file}")

    # モジュール情報も保存
    info_file = RESULTS_DIR / f"module_info_{args.cell_type}.json"
    with open(info_file, 'w') as f:
        json.dump(module_info, f, indent=2)
    print(f"✓ Saved: {info_file}")

    # サマリー
    print(f"\n{'='*60}")
    print(f"Module Scores Summary")
    print(f"{'='*60}")
    print(f"Cells: {len(integrated)}")
    print(f"Modules: {len(module_scores.columns)}")
    print(f"Stress components: {len([c for c in integrated.columns if c.startswith('stress_')])}")
    print(f"Total features: {len(integrated.columns)}")

    # 統計
    print(f"\nModule score statistics:")
    print(integrated[module_scores.columns].describe().loc[['mean', 'std', 'min', 'max']].T.head(10))


if __name__ == "__main__":
    main()
