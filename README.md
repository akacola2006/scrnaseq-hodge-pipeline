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
7. **遺伝子レベル解析 (Lane B)** — 上流セルタイプ内の遺伝子phiスコア
8. **Dual-Mode Hodge** — K_N完全グラフ vs sparse k-NNグラフの比較解析
9. **方向性分解** — Δ⁺（相関増加）/ Δ⁻（相関減少）の分離Hodge
10. **Multi-transition統合** — 全窓遷移にわたるphi推定の統合（4重み付けスキーム）
11. **ランダム行列ベースライン** — 帰無分布からのGF有意性検定
12. **2軸モデル (TRS × MSS)** — 上流遺伝子の翻訳資源スコア × 構造形態スコア

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
├── CLAUDE.md                # Claude Code用ドキュメント
├── README.md                # このファイル
├── project_config.yaml      # パイプライン設定（ユーザー編集）
├── pyproject.toml           # Pythonパッケージ定義
├── requirements.txt         # Python依存パッケージ
├── run_pipeline.py          # メイン実行スクリプト
├── data/                    # データ配置ディレクトリ
│   ├── h5ad/               # h5adファイル
│   ├── metadata/           # メタデータ + 同梱リファレンス
│   │   ├── gene_annotation.csv       # 遺伝子アノテーション（60,403遺伝子）
│   │   └── functional_modules.json   # 24機能モジュール（エンリッチメント用）
│   └── README.md            # データフォーマット説明
├── results/                 # 出力ディレクトリ（自動生成）
└── scripts/                 # パイプラインモジュール
    ├── config.py            # 設定読み込み（project_config.yaml）
    ├── data_loader.py       # データ読み込み・フィルタリング
    ├── residuals.py         # GPU加速回帰残差
    ├── pca_engine.py        # GPU加速PCA（ランダムSVD）
    ├── spd.py               # SPD共分散・Log-Euclidean演算
    ├── pseudotime.py        # 疑似時間構築（拡散マップ）
    ├── hodge.py             # 離散Hodge分解（セルタイプレベル）
    ├── lane_a.py            # Lane A: セルタイプ上流解析
    ├── lane_b.py            # Lane B: 上流PC同定・遺伝子抽出
    ├── gene_hodge.py        # 遺伝子レベルHodge分解（K_N解析的解法）
    ├── bootstrap.py         # ブートストラップ信頼区間
    ├── sparse_hodge.py      # Dual-Mode: K_N vs sparse k-NN Hodge
    ├── directional.py       # 方向性分解（Δ⁺/Δ⁻ split）
    ├── multi_transition.py  # Multi-transition phi統合
    ├── random_baseline.py   # ランダム行列ベースライン（GF帰無分布）
    ├── two_axis.py          # 2軸モデル（TRS × MSS）
    ├── enrichment.py        # 機能モジュールエンリッチメント（Fisher検定）
    ├── seed_utils.py        # 再現性のためのシード管理
    └── log_utils.py         # 失敗ログ
```

## パイプラインステップ

```
python run_pipeline.py --list-steps

 1. validate        — 設定・データの検証
 2. load_data       — h5adファイル読み込み
 3. residuals       — 回帰残差計算（GPU加速）
 4. pca             — ランダムPCA（GPU加速）
 5. spd             — SPD共分散推定（Ledoit-Wolf）
 6. pseudotime      — 疑似時間構築（拡散マップ）
 7. lane_a          — セルタイプレベルHodge（上流セルタイプ同定）
 8. bootstrap       — ドナーブートストラップ（B=100）
 9. lane_b          — 遺伝子セット抽出（上流PC）
10. gene_hodge      — 遺伝子レベルHodge分解
11. dual_mode       — K_N vs sparse k-NN 比較
12. enrichment      — 機能モジュールエンリッチメント
```

追加の解析モジュール（`run_pipeline.py` から独立して使用可能）：

```python
from scripts.directional import run_directional_decomposition   # Δ⁺/Δ⁻ 分離
from scripts.multi_transition import run_multi_transition        # 全遷移統合
from scripts.random_baseline import run_random_baseline          # GF帰無検定
from scripts.two_axis import run_two_axis_model                  # TRS × MSS
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

### Dual-Mode: K_N vs Sparse k-NN
K_N（完全グラフ）ではGFが構造的ベースライン（~1/φ）に縮退する問題がある。
sparse k-NNグラフではこの縮退が破れ、局所的なcurl構造が可視化される。
パイプラインは両モードを自動比較し、phiの一貫性（Spearman ρ）を検証する。

### ランダム行列ベースライン
同サイズのランダム対称行列でGF帰無分布を生成し、観測GFが
構造的アーティファクトでないことを統計的に検証する。

### 方向性分解 (Δ⁺ / Δ⁻)
相関変化行列Δを：
- **Δ⁺**（相関が増加した要素）→ 協調的活性化の駆動源
- **Δ⁻**（相関が減少した要素）→ 脱共役の駆動源

に分離し、それぞれ独立にHodge分解。上流遺伝子が「何を駆動しているか」を解明する。

### 2軸モデル (TRS × MSS)
- **TRS** (Translation-Resource Score): High-phi遺伝子間のlog相関の平均 → 翻訳資源の協調度
- **MSS** (Morphogenesis-Structure Score): Medium-phi遺伝子のPC1 → 構造プログラム活性

窓ごとの（TRS, MSS）軌跡がTRS先行パターンを示す場合、supply-demand構造を示唆する。

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
