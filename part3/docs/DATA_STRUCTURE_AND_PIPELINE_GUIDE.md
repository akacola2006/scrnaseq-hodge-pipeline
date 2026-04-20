# ALS Motor Cortex 解析: データ構造とパイプライン完全ガイド

**作成日**: 2025-11-22
**目的**: プロジェクトの生データから最終解析結果までの全体像を把握

---

## 📋 目次

1. [データ構造概要](#1-データ構造概要)
2. [生データファイル](#2-生データファイル)
3. [解析パイプライン](#3-解析パイプライン)
4. [中間データ](#4-中間データ)
5. [最終結果データ](#5-最終結果データ)
6. [データフロー図](#6-データフロー図)
7. [ファイル命名規則](#7-ファイル命名規則)

---

## 1. データ構造概要

### プロジェクト規模
- **細胞型**: 43種類
- **遺伝子**: 36,503 (ENSG ID)
- **総細胞数**: >100,000
- **機能モジュール**: 23種類
- **サンプル**: ALS患者 vs Control

### データ階層
```
生データ (Raw)
  ↓
前処理済みメトリクス (Processed)
  ↓
モジュール別メトリクス (Module-level)
  ↓
Pseudotime付与 (With PT)
  ↓
最終解析結果 (Final Results)
```

---

## 2. 生データファイル

### 2.1 発現データ

**ファイル形式**:
```
motor_cortex_{細胞型}_expression.csv(.gz)
```

**例**:
```
motor_cortex_Ex.L2_L3.CUX2_RASGRF2_expression.csv
motor_cortex_Ex.L5.VAT1L_EYA4_expression.csv
motor_cortex_Glia.Astro.GFAP-pos_expression.csv.gz
```

**データ構造**:
```
形式: 遺伝子(行) × 細胞ID(列)
行数: 36,504 (ヘッダー + 36,503遺伝子)
列数: 数千〜1万 (細胞型により異なる)
```

**CSV形式**:
```csv
,AAAGGATTCACCTGTC-110MCX,AAAGGTACATACCAGT-110MCX,...
ENSG00000243485,0,0,...
ENSG00000237613,0,0,...
ENSG00000186092,0,0,...
```

**列の説明**:
- 第1列: 空(index用)
- 第2列以降: 細胞バーコードID
  - 形式: `{16文字バーコード}-{サンプルID}`
  - 例: `AAAGGATTCACCTGTC-110MCX`

**行の説明**:
- 第1行: ヘッダー(細胞ID)
- 第2行以降: ENSG ID(Ensembl遺伝子ID)
  - 形式: `ENSG00000######`
  - 例: `ENSG00000104435` (STMN2遺伝子)

**データ値**:
- 生カウント(raw counts)
- 整数値(0以上)

**ファイルサイズ例**:
```
Ex.L2_L3.CUX2_RASGRF2:  831MB (11,871細胞)
Ex.L3_L5.CUX2_RORB:     563MB (8,000細胞)
Ex.L5.VAT1L_EYA4:        23MB (1,300細胞)
In.PV.PVALB_PTHLH:      6.8MB (圧縮)
Vasc.T_Cell:            小規模
```

### 2.2 メタデータ

**ファイル形式**:
```
motor_cortex_{細胞型}_metadata.csv(.gz)
```

**データ構造**:
```csv
cell_id,condition
AAAGGATTCACCTGTC-110MCX,ALS
AAAGGTACATACCAGT-110MCX,ALS
AAAGTCCCAAGTCATC-110MCX,ALS
...
AATCGACCATCGAACT-311MCX,Control
TAGGTTGTCCAGTGTA-311MCX,Control
...
```

**列の説明**:
- `cell_id`: 細胞バーコードID(発現データの列名と一致)
- `condition`: 疾患状態
  - `ALS`: ALS患者由来
  - `Control`: 健常者由来

**サンプルID規則**:
- `###MCX`: サンプル番号
  - 100番台: ALS患者
  - 300番台: Control

**細胞数の例** (Ex.L2_L3.CUX2_RASGRF2):
```
ALS:     8,166 細胞 (68.8%)
Control: 3,705 細胞 (31.2%)
Total:  11,871 細胞
```

### 2.3 全細胞型リスト

#### 興奮性ニューロン (Ex) - 15種類
```
Ex.L2_L3.CUX2_RASGRF2      # Layer 2/3 上層
Ex.L3_L5.CUX2_RORB
Ex.L3_L5.SCN4B_NEFH
Ex.L4_L5.RORB_FOXO1
Ex.L4_L5.RORB_POU3F2
Ex.L4_L6.RORB_ADGRL4
Ex.L4_L6.RORB_LRRK1
Ex.L5.PCP4_NXPH2
Ex.L5.VAT1L_EYA4           # Layer 5 運動ニューロン (重要)
Ex.L5.VAT1L_THSD4          # Layer 5 運動ニューロン (重要)
Ex.L5_L6.THEMIS_NR4A2
Ex.L5_L6.THEMIS_TMEM233
Ex.L6.TLE4_CCBE1
Ex.L6.TLE4_MEGF11
Ex.L6.TLE4_SEMA3D
```

#### 抑制性ニューロン (In) - 16種類
```
# 5HT3aR陽性 (6種類)
In.5HT3aR.CDH4_CCK
In.5HT3aR.CDH4_SCGN
In.5HT3aR.DISC1_CCK
In.5HT3aR.DISC1_RELN
In.5HT3aR.VIP_CLSTN2
In.5HT3aR.VIP_HTR2C
In.5HT3aR.VIP_LAMA3

# PV陽性 (3種類)
In.PV.PVALB_CEMIP
In.PV.PVALB_MYBPC1
In.PV.PVALB_PTHLH

# SOM陽性 (4種類)
In.SOM.SST_ADAMTS19
In.SOM.SST_BRINP3
In.SOM.SST_GALNT14
In.SOM.SST_NPY

# Rosehip細胞 (2種類)
In.Rosehip.LAMP5_CA3
In.Rosehip.LAMP5_PMEPA1
```

#### グリア (Glia) - 4種類
```
Glia.Astro.GFAP-pos        # GFAP陽性アストロサイト
Glia.Astro.GFAP-neg        # GFAP陰性アストロサイト
Glia.Oligo                 # オリゴデンドロサイト
Glia.Micro                 # ミクログリア
Glia.OPC                   # オリゴ前駆細胞
```

#### 血管系 (Vasc) - 7種類
```
Vasc.Endo.Arterial         # 動脈内皮
Vasc.Endo.Venous           # 静脈内皮
Vasc.Endo.Capillary        # 毛細血管内皮
Vasc.Mural.SMC             # 平滑筋細胞
Vasc.Mural.Pericyte        # ペリサイト
Vasc.Fibro.CLMP_PDGFRA     # 線維芽細胞
Vasc.T_Cell                # T細胞
```

### 2.4 補助データファイル

#### モジュール定義
**ファイル**: `24_functional_modules_fixed.json`

```json
{
  "modules": {
    "Synaptic": ["ENSG00000...", "ENSG00000...", ...],
    "Mitochondria": [...],
    "Inflammation": [...],
    ...
  }
}
```

**23モジュールリスト**:
```
1.  Angiogenesis         - 血管新生
2.  Apoptosis            - アポトーシス
3.  Autophagy            - オートファジー
4.  Calcium_Signaling    - カルシウムシグナル
5.  Cell_Cycle           - 細胞周期
6.  Complement           - 補体系
7.  Cytoskeleton         - 細胞骨格
8.  DNA_Repair           - DNA修復
9.  ECM                  - 細胞外マトリクス
10. ER_Stress            - 小胞体ストレス
11. Epigenetic           - エピジェネティクス
12. Growth_Factors       - 成長因子
13. Inflammation         - 炎症
14. Ion_Transport        - イオン輸送
15. Metabolism           - 代謝
16. Mitochondria         - ミトコンドリア
17. Myelination          - ミエリン化
18. Oxidative_Stress     - 酸化ストレス
19. Protein_Homeostasis  - タンパク質恒常性
20. RNA_Processing       - RNA処理
21. Synaptic             - シナプス (重要)
22. Transcription        - 転写
23. lncRNA               - 長鎖非コードRNA
```

#### 遺伝子マッピング
**ファイル**: `gene_mapping.json`

```json
{
  "ensg_to_symbol": {
    "ENSG00000104435": "STMN2",
    "ENSG00000089280": "FUS",
    "ENSG00000142168": "SOD1",
    ...
  }
}
```

---

## 3. 解析パイプライン

### 3.1 ステップ1: モジュールメトリクス計算

**スクリプト**: `complete_batch_processor.py`

**入力**:
- `motor_cortex_{細胞型}_expression.csv(.gz)`
- `motor_cortex_{細胞型}_metadata.csv(.gz)`
- `24_functional_modules_fixed.json`
- `gene_mapping.json`

**実行例**:
```bash
python3 complete_batch_processor.py 8 synFull
# 8 = Ex.L5.VAT1L_EYA4 (細胞型インデックス: 1-43)
# synFull = バリアント名
```

**処理内容**:

1. **データ読み込み**
   - 発現データ読み込み(CSV/GZ自動判定)
   - メタデータ読み込み
   - 共通細胞の抽出

2. **Fast/Slow遺伝子分類**
   ```python
   # Control平均で中心化
   ctrl_mean = expr[ctrl_cells].mean()
   centered_expr = expr - ctrl_mean

   # 分散計算
   gene_vars = centered_expr.var()

   # パーセンタイル分類
   fast_threshold = np.percentile(gene_vars, 70)  # 上位30%
   slow_threshold = np.percentile(gene_vars, 30)  # 下位30%

   fast_genes = genes[gene_vars >= fast_threshold]
   slow_genes = genes[gene_vars <= slow_threshold]
   ```

3. **エネルギー計算**
   ```python
   # 各細胞、各モジュールについて
   EF = mean((centered_expr[fast_genes])²)  # Fast energy
   ES = mean((centered_expr[slow_genes])²)  # Slow energy
   ```

4. **メトリクス計算**
   ```python
   # Module Intensity
   I_m = √(EF + ES)

   # Scale Ratio
   SR = √(EF / ES)

   # φ偏差 (生値)
   φ = 1.618  # 黄金比
   delta_phi_m = |SR - φ|
   ```

**出力**:
```
complete_ids_results/cell_module_{細胞型}_synFull.csv
```

**出力形式**:
```csv
cell_id,cell_type,variant,module,Disease,delta_phi_m,I,SR,EF,ES,n_fast,n_slow,n_genes,source
AAAA...-110MCX,Ex.L5.VAT1L_EYA4,synFull,Synaptic,ALS,48.52,10.23,50.14,102.1,2.04,15,15,50,metrics_...csv
...
```

**列の説明**:
- `cell_id`: 細胞バーコードID
- `cell_type`: 細胞型名
- `variant`: 解析バリアント(synFull/synIEG)
- `module`: モジュール名
- `Disease`: ALS/Control/CTRL
- `delta_phi_m`: Δφ(生値) = |SR - φ|
- `I`: Module Intensity = √(EF + ES)
- `SR`: Scale Ratio = √(EF/ES)
- `EF`: Fast遺伝子のエネルギー
- `ES`: Slow遺伝子のエネルギー
- `n_fast`: Fast遺伝子数
- `n_slow`: Slow遺伝子数
- `n_genes`: モジュール内総遺伝子数

### 3.2 ステップ2: Pseudotime (PT) 計算

**スクリプト**: `make_continuous_pt_T_only.py`

**入力**:
```
complete_ids_results/cell_module_*.csv (全966ファイル)
```

**実行例**:
```bash
python3 make_continuous_pt_T_only.py \
  --tcell complete_ids_results/cell_module_Vasc_T_Cell_synFull.csv \
  --out pt_results/tcell_pt_continuous.csv \
  --k 15 \
  --bins 10
```

**処理内容**:

1. **特徴行列構築**
   ```python
   # 46次元特徴空間
   for each cell:
       features = [
           I_m[module_1], ..., I_m[module_23],      # 23次元
           Δφ_signed[module_1], ..., Δφ_signed[module_23]  # 23次元
       ]

   # Δφ_signed計算
   Δφ_signed = log(SR) - log(φ)
              = 0.5 * log(EF/ES) - log(1.618)
   ```

2. **kNNグラフ構築**
   ```python
   # k近傍探索 (k=15)
   nbrs = NearestNeighbors(n_neighbors=15, metric='euclidean')
   distances, indices = nbrs.kneighbors(features)

   # Gaussian kernel重み付き隣接行列
   weight[i,j] = exp(-distance²/(2*σ²))
   # 対称化
   adjacency[i,j] = adjacency[j,i] = weight
   ```

3. **拡散マップ計算**
   ```python
   # 度数行列
   D = diag(sum(adjacency, axis=1))

   # 正規化ラプラシアン
   D_inv_sqrt = diag(1/√(D))
   L_norm = D_inv_sqrt @ adjacency @ D_inv_sqrt

   # 固有値分解
   eigenvals, eigenvecs = eigsh(L_norm, k=10, which='LA')

   # 第2固有ベクトル = PT_raw
   pt_raw = eigenvecs[:, 1]
   ```

4. **堅牢な正規化**
   ```python
   # Controlの下位5%群のメジアン → PT = 0
   ctrl_anchor = median(pt_raw[ctrl_cells][pt_raw <= p05])

   # ALSの上位95%群のメジアン → PT = 1
   als_anchor = median(pt_raw[als_cells][pt_raw >= p95])

   # 正規化
   PT_continuous = (pt_raw - ctrl_anchor) / (als_anchor - ctrl_anchor)
   PT_continuous = clip(PT_continuous, 0, 1)
   ```

**出力**:
```csv
cell_id,condition,PT_continuous,PT_raw
AAAA...-110MCX,ALS,0.756,2.341
BBBB...-311MCX,Control,0.234,-1.123
...
```

### 3.3 ステップ3: PTの統合

**処理**: 各細胞型の全モジュールファイルにPT_continuousを結合

**出力ディレクトリ**:
```
complete_ids_results_withPTcontinuous/
```

**出力ファイル** (966ファイル):
```
cell_module_{細胞型}_synFull_withPTcontinuous.csv
```

**出力形式**:
```csv
cell_id,cell_type,variant,module,Disease,delta_phi_m,I_m,SR_m,SR_log,EF_m,ES_m,delta_phi_log,PT_continuous
AAAA...,Ex.L5.VAT1L_EYA4,synFull,Synaptic,ALS,48.52,10.23,50.14,3.915,102.1,2.04,3.434,0.756
...
```

**新規列**:
- `I_m`: Module Intensity (列名統一)
- `SR_m`: Scale Ratio (生値)
- `SR_log`: log(SR) = 0.5*log(EF/ES)
- `delta_phi_log`: Δφ_signed = SR_log - log(φ)
- `PT_continuous`: Pseudotime [0,1]

### 3.4 ステップ4: Trajectory集計

**スクリプト**: `create_trajectory_summary_table.py`

**入力**:
```
complete_ids_results_withPTcontinuous/cell_module_*_synFull_withPTcontinuous.csv
(全966ファイル)
```

**実行**:
```bash
python3 create_trajectory_summary_table.py
```

**処理内容**:

1. **PT区間分割**
   ```python
   PT_BINS = [0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0]
   PT_LABELS = ['0.50-0.55', '0.55-0.60', ..., '0.95-1.00']

   # 10区間に分割
   pt_bin = pd.cut(PT_continuous, bins=PT_BINS, labels=PT_LABELS)
   ```

2. **各区間の統計量計算**
   ```python
   for each (cell_type, module, pt_bin):
       module_intensity_mean = mean(I_m)
       module_intensity_sem = sem(I_m)  # 標準誤差
       n_cells = count(cells)
   ```

3. **遺伝子発現の統合** (STMN2, FUS)
   ```python
   # 発現データを読み込み
   expr = read_csv(motor_cortex_{cell_type}_expression.csv)

   # STMN2, FUSの発現量を抽出
   stmn2_expr = expr.loc['ENSG00000104435', cells]
   fus_expr = expr.loc['ENSG00000089280', cells]

   # PT区間ごとに平均
   for pt_bin:
       STMN2_mean = mean(stmn2_expr[cells_in_bin])
       FUS_mean = mean(fus_expr[cells_in_bin])
   ```

**出力**:
```
module_trajectory_summary_10bins.csv
```

**出力形式**:
```csv
cell_type,module,pt_bin,module_intensity_mean,module_intensity_sem,n_cells,STMN2_mean,FUS_mean
Ex.L5.VAT1L_EYA4,Synaptic,0.50-0.55,8.234,0.123,245,12.34,8.76
Ex.L5.VAT1L_EYA4,Synaptic,0.55-0.60,8.567,0.145,223,11.89,8.54
...
Ex.L2_L3.CUX2_RASGRF2,Synaptic,0.50-0.55,3.456,0.089,512,15.67,9.23
...
```

**データ量**:
```
42細胞型 × 23モジュール × 10 PT bins = 9,660行
```

---

## 4. 中間データ

### 4.1 モジュールメトリクス

**ディレクトリ**: `complete_ids_results/`

**ファイル数**: 42ファイル (細胞型ごと)

**ファイル名**:
```
cell_module_{細胞型}_synFull.csv
```

**内容**:
- 各細胞×23モジュールのメトリクス
- EF, ES, I, SR, Δφ
- 疾患ラベル

**用途**:
- PT計算の入力
- 初期品質確認

### 4.2 epsilon値ログ

**ファイル**: `epsilon_log.csv`

**内容**:
```csv
cell_type,module,variant,eps_value
Vasc.T_Cell,Synaptic,synFull,1e-8
Vasc.T_Cell,Mitochondria,synFull,1e-8
...
```

**用途**:
- ゼロ除算防止のための微小値
- SR計算時に使用: SR = √((EF + ε)/(ES + ε))

---

## 5. 最終結果データ

### 5.1 PT付きモジュールメトリクス

**ディレクトリ**: `complete_ids_results_withPTcontinuous/`

**ファイル数**: 966ファイル

**命名規則**:
```
cell_module_{細胞型}_synFull_withPTcontinuous.csv
```

**データ規模** (例: Ex.L5.VAT1L_EYA4):
```
細胞数: 1,234
モジュール数: 23
行数: 1,234 × 23 = 28,382行
```

**主要列**:
```
cell_id          - 細胞バーコードID
cell_type        - 細胞型名
module           - モジュール名
Disease          - ALS/Control
I_m              - Module Intensity
SR_m             - Scale Ratio (生値)
SR_log           - log(SR)
EF_m             - Fast energy
ES_m             - Slow energy
delta_phi_m      - Δφ (生値)
delta_phi_log    - Δφ_signed (対数)
PT_continuous    - Pseudotime [0,1]
```

**用途**:
- 個別細胞レベルの詳細解析
- Trajectory可視化
- 統計解析

### 5.2 Trajectory集計データ

**ファイル**: `module_trajectory_summary_10bins.csv`

**データ規模**:
```
行数: 9,660
  = 42細胞型 × 23モジュール × 10 PT bins
列数: 8
```

**列の詳細**:
```csv
cell_type                 - 細胞型名
module                    - モジュール名
pt_bin                    - PT区間 (0.50-0.55等)
module_intensity_mean     - I_m平均
module_intensity_sem      - I_m標準誤差
n_cells                   - 細胞数
STMN2_mean               - STMN2平均発現
FUS_mean                 - FUS平均発現
```

**用途**:
- Module Trajectory可視化
- 疾患進行パターン解析
- VAT1L vs L2/3比較
- 論文Figure作成

**重要なパターン例**:
```python
# VAT1L (減少パターン)
pt_bin      I_m_mean
0.50-0.55   8.23
0.55-0.60   8.56
0.60-0.65   9.12
0.65-0.70   9.82  # ピーク
0.70-0.75   9.45
0.75-0.80   8.67
0.80-0.85   7.89
0.85-0.90   6.91
0.90-0.95   5.92
0.95-1.00   5.01  # 低下

# L2/3 (増加パターン)
pt_bin      I_m_mean
0.50-0.55   2.94
0.55-0.60   3.12
...
0.90-0.95   5.01
0.95-1.00   5.42  # 増加
```

### 5.3 ペアワイズ因果性ネットワーク

**ディレクトリ**: `pairwise_batch/`

**主要ファイル**: `pairwise_summary.csv`

**データ規模**:
```
行数: 903ペア
  = 42細胞型 × (42-1) / 2 ≈ 903
```

**列の構成**:
```csv
pair                - ペア名 (A_vs_B)
majority_decision   - 因果方向 (A→B, B→A, Uncertain)
support_A_to_B      - A→Bのサポートスコア
support_B_to_A      - B→Aのサポートスコア
O1_result           - O1テスト結果
O2_result           - O2テスト結果
O3_result           - O3テスト結果
O1_value            - O1統計量
O2_value            - O2統計量
O3_value            - O3統計量
```

**3つの因果性テスト**:
```
O1: Lag asymmetry      - 時間遅れ非対称性
O2: Conditional asymmetry - 条件付き独立性
O3: Energy flow        - エネルギー流方向

判定: Majority vote (3つ中2つ以上の合意)
```

**用途**:
- 細胞型間因果ネットワーク構築
- 5層階層構造の抽出
- 疾患カスケードの解明

---

## 6. データフロー図

```
[生データレイヤー]
  motor_cortex_*_expression.csv.gz (43ファイル)
  motor_cortex_*_metadata.csv.gz (43ファイル)
  24_functional_modules_fixed.json
  gene_mapping.json
        ↓
        ↓ complete_batch_processor.py
        ↓   - ENSG→Symbol変換
        ↓   - Control平均中心化
        ↓   - Fast/Slow分類 (70/30パーセンタイル)
        ↓   - EF, ES, I, SR, Δφ計算
        ↓
[中間レイヤー1]
  complete_ids_results/
    cell_module_{細胞型}_synFull.csv (42ファイル)
        ↓
        ↓ make_continuous_pt_T_only.py
        ↓   - 46次元特徴空間構築
        ↓   - kNNグラフ (k=15)
        ↓   - 拡散マップ
        ↓   - PT正規化
        ↓
[中間レイヤー2]
  complete_ids_results_withPTcontinuous/
    cell_module_{細胞型}_synFull_withPTcontinuous.csv (966ファイル)
        ↓
        ↓ create_trajectory_summary_table.py
        ↓   - PT区間分割 (10 bins)
        ↓   - I_m平均・SEM計算
        ↓   - STMN2/FUS発現統合
        ↓
[最終レイヤー]
  module_trajectory_summary_10bins.csv
    - 9,660行 (42 × 23 × 10)
    - 論文Figure用データ
```

---

## 7. ファイル命名規則

### 7.1 細胞型名の表記

**生データファイル**: ドット区切り
```
motor_cortex_Ex.L5.VAT1L_EYA4_expression.csv
motor_cortex_Glia.Astro.GFAP-pos_metadata.csv
```

**処理済みファイル**: アンダースコア区切り
```
cell_module_Ex_L5_VAT1L_EYA4_synFull.csv
cell_module_Glia_Astro_GFAP_pos_synFull_withPTcontinuous.csv
```

**変換ルール**:
```python
normalized_name = cell_type.replace('.', '_').replace('-', '_')
```

### 7.2 ファイルサフィックス

**発現データ**:
```
_expression.csv      - 非圧縮CSV
_expression.csv.gz   - gzip圧縮
```

**メタデータ**:
```
_metadata.csv
_metadata.csv.gz
```

**処理済みデータ**:
```
_synFull.csv                      - バリアント名
_synFull_withPTcontinuous.csv     - PT付き
_synFull_withPTlog.csv            - 旧形式(PT対数値)
```

### 7.3 ディレクトリ構造

```
/home/akaco/als/motor_cortex_analysis/
│
├── motor_cortex_*_expression.csv(.gz)    # 生発現データ (43ファイル)
├── motor_cortex_*_metadata.csv(.gz)      # メタデータ (43ファイル)
├── 24_functional_modules_fixed.json      # モジュール定義
├── gene_mapping.json                     # 遺伝子マッピング
│
├── complete_ids_results/                 # 中間結果
│   └── cell_module_*_synFull.csv         # 42ファイル
│
├── complete_ids_results_withPTcontinuous/ # PT付き結果
│   └── cell_module_*_synFull_withPTcontinuous.csv  # 966ファイル
│
├── pairwise_batch/                       # ペアワイズ解析
│   ├── pairwise_summary.csv              # 集計結果
│   └── individual_results/               # 個別ファイル
│
├── module_trajectory_summary_10bins.csv  # Trajectory集計
│
├── complete_batch_processor.py           # スクリプト1
├── make_continuous_pt_T_only.py          # スクリプト2
├── create_trajectory_summary_table.py    # スクリプト3
│
└── [各種ドキュメント].md
```

---

## 8. データアクセス例

### 8.1 Python

#### 発現データ読み込み
```python
import pandas as pd

# CSV読み込み
expr = pd.read_csv('motor_cortex_Ex.L5.VAT1L_EYA4_expression.csv', index_col=0)

# GZ読み込み
expr = pd.read_csv('motor_cortex_In.PV.PVALB_PTHLH_expression.csv.gz',
                   index_col=0, compression='gzip')

# 形式: 遺伝子(行) × 細胞(列)
print(expr.shape)  # (36503, 1234)

# 特定遺伝子の発現
stmn2_expr = expr.loc['ENSG00000104435', :]

# 特定細胞の全遺伝子発現
cell_expr = expr.loc[:, 'AAAGGATTCACCTGTC-110MCX']
```

#### メタデータ読み込み
```python
meta = pd.read_csv('motor_cortex_Ex.L5.VAT1L_EYA4_metadata.csv')

# ALS細胞のみ
als_cells = meta[meta['condition'] == 'ALS']['cell_id'].tolist()

# Control細胞のみ
ctrl_cells = meta[meta['condition'] == 'Control']['cell_id'].tolist()

# 発現データフィルタリング
als_expr = expr.loc[:, als_cells]
```

#### Trajectoryデータ読み込み
```python
traj = pd.read_csv('module_trajectory_summary_10bins.csv')

# VAT1Lのみ
vat1l = traj[traj['cell_type'] == 'Ex_L5_VAT1L_EYA4']

# Synapticモジュールのみ
syn = vat1l[vat1l['module'] == 'Synaptic']

# PT区間ごとの平均
print(syn[['pt_bin', 'module_intensity_mean']])
```

#### PT付きデータ読み込み
```python
# 個別細胞データ
data = pd.read_csv('complete_ids_results_withPTcontinuous/'
                   'cell_module_Ex_L5_VAT1L_EYA4_synFull_withPTcontinuous.csv')

# ALS細胞のみ
als_data = data[data['Disease'] == 'ALS']

# 特定モジュール
syn_data = als_data[als_data['module'] == 'Synaptic']

# PT区間でフィルタ
late_stage = syn_data[syn_data['PT_continuous'] >= 0.8]
```

### 8.2 R

```r
# 発現データ
expr <- read.csv('motor_cortex_Ex.L5.VAT1L_EYA4_expression.csv', row.names=1)

# メタデータ
meta <- read.csv('motor_cortex_Ex.L5.VAT1L_EYA4_metadata.csv')

# Trajectoryデータ
traj <- read.csv('module_trajectory_summary_10bins.csv')

# VAT1L Synaptic
vat1l_syn <- traj[traj$cell_type == 'Ex_L5_VAT1L_EYA4' &
                   traj$module == 'Synaptic', ]

# プロット
plot(vat1l_syn$pt_bin, vat1l_syn$module_intensity_mean, type='b')
```

### 8.3 Bash

```bash
# 細胞数カウント
wc -l motor_cortex_Ex.L5.VAT1L_EYA4_metadata.csv

# ALS細胞のみ抽出
awk -F',' '$2=="ALS" {print $1}' motor_cortex_Ex.L5.VAT1L_EYA4_metadata.csv

# 遺伝子数カウント
head -1 motor_cortex_Ex.L5.VAT1L_EYA4_expression.csv | tr ',' '\n' | wc -l

# STMN2発現抽出
grep "ENSG00000104435" motor_cortex_Ex.L5.VAT1L_EYA4_expression.csv
```

---

## 9. 重要な技術的詳細

### 9.1 φ-アトラクター理論パラメータ

```python
# 定数
PHI = (1 + np.sqrt(5)) / 2  # 1.618... (黄金比)
FAST_PERCENTILE = 70         # Fast遺伝子の閾値
SLOW_PERCENTILE = 30         # Slow遺伝子の閾値
MIN_GENES = 20               # モジュール最小遺伝子数

# 計算式
EF = mean((expr[fast_genes] - ctrl_mean)²)
ES = mean((expr[slow_genes] - ctrl_mean)²)

I_m = √(EF + ES)
SR = √(EF / ES)
Δφ_signed = log(SR) - log(φ)
          = 0.5 * log(EF/ES) - log(1.618)
```

### 9.2 Pseudotime正規化

```python
# アンカーポイント
ctrl_anchor = median(PT_raw[ctrl_cells][PT_raw <= p05])
als_anchor = median(PT_raw[als_cells][PT_raw >= p95])

# 正規化
PT = (PT_raw - ctrl_anchor) / (als_anchor - ctrl_anchor)
PT = clip(PT, 0, 1)

# PT範囲
解析対象: PT ∈ [0.5, 1.0]  # 疾患進行期
除外: PT < 0.5              # 健常/初期(細胞数少)
```

### 9.3 統計量

```python
# 標準誤差 (SEM)
sem = std(I_m) / √(n_cells)

# Pearson相関
from scipy.stats import pearsonr
r, p = pearsonr(PT_bins, I_m_means)

# 例: VAT1L Synaptic
r = -0.919  # 強い負の相関
p < 0.001   # 高度有意
```

---

## 10. よくある質問

### Q1: なぜ43細胞型なのに42と書かれている？

**A**: プロジェクトドキュメントでは42細胞型と記載されていますが、実際のファイルは43個存在します。いずれかの細胞型が除外されたか、カウント方法の違いによるものと思われます。

### Q2: synFullとsynIEGの違いは？

**A**:
- `synFull`: Synapticモジュールの全遺伝子を使用
- `synIEG`: Immediate Early Genes(IEG)のみを使用

現在のプロジェクトでは主にsynFullが使用されています。

### Q3: PT < 0.5のデータはどこ?

**A**: PT範囲全体は[0,1]ですが、解析では疾患進行期の[0.5,1.0]のみを使用しています。PT < 0.5は健常または初期状態で細胞数が少ないためです。

### Q4: ファイルサイズがなぜ異なる?

**A**: 細胞数の違いによります:
- Ex.L2_L3.CUX2_RASGRF2: 11,871細胞 → 831MB
- Vasc.T_Cell: 数百細胞 → 数MB(圧縮)

### Q5: Zone.Identifierファイルとは?

**A**: Windows環境でダウンロードされたファイルに付与される代替データストリームです。無視して構いません。

---

## 11. トラブルシューティング

### メモリ不足エラー

```python
# チャンク読み込み
for chunk in pd.read_csv('large_file.csv', chunksize=10000):
    process(chunk)

# メモリ解放
del large_df
import gc
gc.collect()
```

### ファイルが見つからない

```bash
# ファイル検索
find . -name "*VAT1L*" -type f

# ディレクトリ確認
ls -lh complete_ids_results_withPTcontinuous/ | head
```

### 圧縮ファイル読み込みエラー

```python
# 自動判定
df = pd.read_csv('file.csv.gz', compression='infer', index_col=0)

# 明示的指定
df = pd.read_csv('file.csv.gz', compression='gzip', index_col=0)
```

---

## 12. 更新履歴

| 日付 | バージョン | 変更内容 |
|------|-----------|---------|
| 2025-11-22 | 1.0 | 初版作成 |

---

**このドキュメントは完全版です。次回セッション時もこのファイルを参照してください。**
