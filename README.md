# scRNAseq Hodge Decomposition Pipeline

Single-cell RNA-seq解析パイプライン。離散Hodge分解を用いて、疾患や条件に関連する
**上流セルタイプ**と**上流遺伝子**を同定します。

## 概要

このパイプラインは、scRNAseqデータから以下を自動的に実行します：

1. **データ前処理** — QC、正規化 (log1p CPM)、回帰残差 (donor + nUMI + pct_mito)
2. **次元圧縮** — GPU加速ランダムPCA
3. **共分散構造推定** — Ledoit-Wolf shrinkage + SPD多様体上のLog-Euclidean距離
4. **疑似時間構築** — 拡散マップによるPT-B
5. **セルタイプレベル解析 (Lane A)** — 完全グラフ上の離散Hodge分解 → 上流セルタイプ同定
6. **ブートストラップ検証** — ドナーレベルクラスターブートストラップ（B=100）
7. **遺伝子レベル解析 (Lane B)** — 上流セルタイプ内の遺伝子phi スコア

## クイックスタート

### 1. 環境セットアップ
```bash
pip install -r requirements.txt
```

### 2. データ配置
```
data/
  h5ad/           <- h5adファイルを配置（サンプルごとに1ファイル）
  metadata/
    sample_info.csv  <- サンプルメタデータ（donor_id, file, condition）
```

### 3. 設定編集
`project_config.yaml` を編集して、セルタイプ名、条件ラベル等を設定

### 4. 実行
```bash
python run_pipeline.py
```

## ファイル構成

```
scrnaseq_hodge_pipeline/
├── CLAUDE.md              # Claude Code用ドキュメント
├── README.md              # このファイル
├── project_config.yaml    # パイプライン設定（ユーザー編集）
├── requirements.txt       # Python依存パッケージ
├── run_pipeline.py        # メイン実行スクリプト
├── data/                  # データ配置ディレクトリ
│   ├── h5ad/             # h5adファイル
│   ├── metadata/         # メタデータ
│   └── README.md         # データフォーマット説明
├── results/              # 出力ディレクトリ（自動生成）
└── scripts/              # パイプラインモジュール
    ├── config.py         # 設定読み込み
    ├── data_loader.py    # データ読み込み・フィルタリング
    ├── residuals.py      # GPU加速回帰残差
    ├── pca_engine.py     # GPU加速PCA
    ├── spd.py            # SPD共分散・Log-Euclidean演算
    ├── pseudotime.py     # 疑似時間構築
    ├── hodge.py          # 離散Hodge分解
    ├── lane_a.py         # セルタイプレベル上流解析
    ├── lane_b.py         # 上流PC同定・遺伝子抽出
    ├── gene_hodge.py     # 遺伝子レベルHodge分解
    ├── bootstrap.py      # ブートストラップ信頼区間
    ├── seed_utils.py     # 再現性のためのシード管理
    └── log_utils.py      # 失敗ログ
```

## 必要なデータ形式

### h5adファイル（サンプルごとに1ファイル）
- `.X` — 生カウント行列 (細胞 x 遺伝子)
- `.obs` — 細胞メタデータ（`CellType`列を含む）
- `.var_names` — 遺伝子識別子

### sample_info.csv
| donor_id | file         | condition |
|----------|--------------|-----------|
| D001     | sample1.h5ad | Disease   |
| D002     | sample2.h5ad | Control   |

詳細は `data/README.md` を参照してください。

## 数学的背景

### 離散Hodge分解
完全グラフ K_N 上のエッジフローを3つの直交成分に分解：
```
f = gradient + curl + harmonic
```
- **gradient** (B0 @ phi): ポテンシャル場（上流/下流の順序）
- **curl** (B1^T @ c): 循環的フロー
- **harmonic**: 残差

セルタイプの**phi（ポテンシャル）スコア**が高いほど、疾患進行の「上流」に位置します。

### Log-Euclidean距離
SPD (対称正定値) 共分散行列間の距離を行列対数空間のFrobeniusノルムで測定：
```
d_LE(A, B) = || log(A) - log(B) ||_F
```

## GPU加速

- 残差計算: PyTorch GPU-batched OLS
- PCA: `torch.pca_lowrank` (ランダムSVD)
- SPD batch log/exp: PyTorch batch eigendecomposition

GPUが利用できない場合は自動的にCPUにフォールバックします。

## Claude Codeでの使い方

このフォルダをClaude Codeで開くと、CLAUDE.mdに従って解析を実行できます。

1. `data/h5ad/` にh5adファイルを配置
2. `data/metadata/sample_info.csv` を作成
3. Claude Codeに「解析を実行して」と指示
4. Claude Codeが `project_config.yaml` を設定し、`run_pipeline.py` を実行
5. 結果の解釈をClaude Codeに依頼

## ライセンス

研究利用を目的としています。
