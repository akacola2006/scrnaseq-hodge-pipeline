# Motor Cortex 全細胞タイプリスト

## データ概要
- 総ファイル数: 174ファイル (CSV 21個 + GZ圧縮 66個 + Zone.Identifier 87個)
- 総データサイズ: 約3.9GB
- 細胞タイプカテゴリー: 7種類

## 細胞タイプ一覧

### 1. 興奮性神経細胞 (Excitatory Neurons) - Ex
- Ex.L2_L3.CUX2_RASGRF2
- Ex.L3_L5.CUX2_RORB ✓
- Ex.L3_L5.SCN4B_NEFH
- Ex.L4_L5.RORB_FOXO1
- Ex.L4_L5.RORB_POU3F2
- Ex.L4_L6.RORB_ADGRL4
- Ex.L4_L6.RORB_LRRK1
- Ex.L5.PCP4_NXPH2
- Ex.L5.VAT1L_EYA4
- Ex.L5.VAT1L_THSD4
- Ex.L5_L6.THEMIS_NR4A2
- Ex.L5_L6.THEMIS_TMEM233
- Ex.L6.TLE4_CCBE1
- Ex.L6.TLE4_MEGF11
- Ex.L6.TLE4_SEMA3D

### 2. 抑制性神経細胞 (Inhibitory Neurons) - In
#### 2a. 5HT3aR陽性
- In.5HT3aR.CDH4_CCK
- In.5HT3aR.CDH4_SCGN
- In.5HT3aR.DISC1_CCK
- In.5HT3aR.DISC1_RELN
- In.5HT3aR.VIP_CLSTN2
- In.5HT3aR.VIP_HTR2C
- In.5HT3aR.VIP_LAMA3

#### 2b. PV陽性
- In.PV.PVALB_CEMIP
- In.PV.PVALB_MYBPC1
- In.PV.PVALB_PTHLH

#### 2c. Rosehip
- In.Rosehip.LAMP5_CA3
- In.Rosehip.LAMP5_PMEPA1

#### 2d. SOM陽性
- In.SOM.SST_ADAMTS19
- In.SOM.SST_BRINP3
- In.SOM.SST_GALNT14
- In.SOM.SST_NPY

### 3. グリア細胞 (Glia)
- Glia.Astro.GFAP-neg ✓
- Glia.Astro.GFAP-pos
- Glia.Micro (ミクログリア)
- Glia.OPC (オリゴデンドロサイト前駆細胞)
- Glia.Oligo (オリゴデンドロサイト)

### 4. 血管系細胞 (Vascular)
#### 4a. 内皮細胞
- Vasc.Endo.Arterial (動脈内皮)
- Vasc.Endo.Capillary (毛細血管内皮) ✓
- Vasc.Endo.Venous (静脈内皮)

#### 4b. 壁細胞
- Vasc.Mural.Pericyte (ペリサイト)
- Vasc.Mural.SMC (平滑筋細胞)

#### 4c. その他
- Vasc.Fibro.CLMP_PDGFRA (線維芽細胞)
- Vasc.T_Cell (T細胞)

### 5. サマリーファイル
- motor_cortex_cell_type_summary.csv

## ファイル形式
- 非圧縮CSV: 大きな主要細胞タイプ
  - Ex.L2_L3, Ex.L3_L5, Astro.GFAP-neg, Oligo, Vasc.Endo.Capillary等

- GZ圧縮: 小さな細胞タイプ
  - 各種抑制性神経細胞サブタイプ
  - 希少な血管系細胞等

## 注記
✓ マークは現在の解析で主に使用している3細胞タイプを示す

## データ構造
各細胞タイプに対して2ファイル:
1. `*_expression.csv(.gz)`: 遺伝子×細胞の発現マトリックス
2. `*_metadata.csv(.gz)`: 細胞のメタデータ（condition等）