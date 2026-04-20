# 細胞タイプ別モジュール時系列ガイド
## Early → Mid → Late の変遷を読み解く

---

## 📊 生成されたファイル（52個）

### **タイムライン要約CSV**（全43タイプ）
```
cell_module_timeline/{cell_type}_timeline_summary.csv
```

**カラム構成**:
- `module`: モジュール名
- `activation_time`: 発火時刻（PT中央値、0-1）
- `convergence_score`: 収束スコア（正=収束、負=発散）
- `slope_abs_delta_phi`: |Δφ|の傾き（負=収束）
- `late_within_phi_fraction`: late期のφ近接率
- `rho_I_PT_{early/mid/late}`: 各期のI-PT相関
- `q_I_PT_{early/mid/late}`: I-PT相関のBY-FDR
- `rho_dPhi_PT_{early/mid/late}`: 各期の|Δφ|-PT相関
- `rho_I_dPhi_{early/mid/late}`: 各期のI-|Δφ|相関

### **可視化（Focus 4タイプ）**

#### 1. **パネル図** (`*_panel.png`)
- **Ex_L5_VAT1L_EYA4** (944KB) - VAT1L脆弱性の時系列
- **In_PV_PVALB_CEMIP** (825KB) - PVALB抑制性ニューロン
- **Glia_Astro_GFAP_neg** (807KB) - アストロサイト
- **Glia_OPC** (881KB) - OPC（オリゴデンドロサイト前駆細胞）

**パネル図の見方**:
- 各サブプロット = 1モジュールの時系列
- 左軸（青）: |Δφ_signed|（φからの逸脱）
- 右軸（赤）: Eφ median（φエネルギー）
- 縦の破線: Phase境界（0.3, 0.7）
- 枠色: 青=収束（convergence_score>0）、赤=発散（<0）

#### 2. **要約ヒートマップ** (`*_summary_heatmap.png`)
- モジュール × 主要指標のヒートマップ
- Z-scoreで正規化
- 赤=高、青=低

#### 3. **Focus 4比較図** (`focus4_overview.png`, 161KB)
- VAT1L, PVALB, Astro, OPC の Top 10モジュール比較
- Convergence score vs -Slope(|Δφ|)

---

## 🔬 読み方のルール

### **1. Activation（発火）時刻**
```
activation_time が小さい → 早期に逸脱開始
activation_time が大きい → 後期に逸脱開始
```

**例（VAT1L_EYA4）**:
- **Apoptosis**: 0.0028（最早発火）
- **Autophagy**: 0.0028
- **Calcium_Signaling**: 0.2775（最遅発火）

### **2. Convergence（収束 vs 発散）**
```
convergence_score > 0 → 収束（弱いアトラクタに沿う）
convergence_score < 0 → 発散（φから離れ続ける）

slope_abs_delta_phi < 0 → early→lateで|Δφ|が減少
slope_abs_delta_phi > 0 → early→lateで|Δφ|が増加
```

**例（VAT1L_EYA4）**:
- **収束**: Complement (1.0), Epigenetic (0.5)
- **発散**: Apoptosis (-0.77), Calcium_Signaling (-0.40)

### **3. Phase別相関パターン**

#### **early期**（PT 0-0.3）
- `rho_I_PT > 0` → 振幅が増大開始
- `rho_dPhi_PT > 0` → φからの逸脱が増大
- `rho_I_dPhi > 0` → 振幅と逸脱が連動（危険）

#### **mid期**（PT 0.3-0.7）
- `rho_dPhi_PT < 0` → φへの収束開始
- `rho_I_PT > 0` → 振幅は継続増大

#### **late期**（PT 0.7-1.0）
- `rho_dPhi_PT < 0` → 収束継続
- `rho_I_dPhi > 0` → 振幅と逸脱の残存連動

### **4. Late期のφ近接率**
```
late_within_phi_fraction → late期で|Δφ|<0.18の細胞割合
高い値（>0.5）→ late期に多くの細胞がφに戻る（安全）
低い値（<0.1）→ late期でも逸脱が残る（危険）
```

---

## 📖 実例：VAT1L_EYA4の読み解き

### **タイムライン要約から**
```csv
module,activation_time,convergence_score,slope_abs_delta_phi,late_within_phi_fraction
Synaptic,0.003,0.2,-0.01,0.0
Oxidative_Stress,0.003,-0.3,0.05,0.0
Calcium_Signaling,0.278,-0.40,0.07,0.0
Apoptosis,0.003,-0.77,0.08,0.0
```

**解釈**:
1. **Synaptic**: 早期発火（0.003）、収束傾向（0.2）→ 早期介入で修正可能
2. **Oxidative_Stress**: 早期発火、発散傾向（-0.3）→ 早期から悪化
3. **Calcium_Signaling**: 最遅発火（0.278）、発散（-0.40）→ 終末期破綻
4. **Apoptosis**: 早期発火、強発散（-0.77）→ 一貫して悪化

### **Phase相関から**
```csv
module,rho_I_PT_early,rho_dPhi_PT_mid,rho_I_dPhi_late
Synaptic,-0.25,0.09,0.92
Calcium_Signaling,-0.37,0.40,0.56
```

**解釈**:
- **Synaptic**: late期でI-|Δφ|が強相関（0.92）→ 振幅と逸脱が同期（危険増幅）
- **Calcium_Signaling**: mid期から|Δφ|-PTが正相関（0.40）→ 収束失敗

---

## 🎯 タイプ別の特徴（Focus 4）

### **1. VAT1L (Ex_L5_VAT1L_EYA4)**
**特徴**:
- Synaptic, Oxidative_Stress, Calcium_Signaling が**早期発火 + late発散**
- late期でもlate_within_phi_fraction = 0（完全非収束）
- **最も脆弱**な細胞タイプ

**危険モジュール**:
- Apoptosis (convergence -0.77)
- Calcium_Signaling (convergence -0.40)
- Oxidative_Stress (convergence -0.30)

### **2. PVALB (In_PV_PVALB_CEMIP)**
**特徴**:
- Synaptic, Ion_Transport のmid期でI-PT相関が正
- **ネットワーク拡散の媒介役**
- 比較的高いconvergence（Complement, Epigenetic）

**役割**:
- 初期破綻のバッファー → 中期での拡散媒介

### **3. Astro (Glia_Astro_GFAP_neg)**
**特徴**:
- Complement, Protein_Homeostasis が**高収束**（1.0+）
- Ion_Transport, Calcium_Signaling も比較的安定
- **ホメオスタシス維持の要**

**保護的モジュール**:
- Complement (convergence 1.17)
- Protein_Homeostasis (convergence 0.52)

### **4. OPC (Glia_OPC)**
**特徴**:
- Myelination, Cytoskeleton のlate期で微回復
- Inflammation, Apoptosis は発散傾向
- **再生 vs 炎症の拮抗**

**両刃性**:
- Myelination: late回復の可能性
- Inflammation: 持続悪化

---

## 💡 介入ポイントの読み取り

### **Early介入（PT 0-0.3）**
**標的**: `activation_time < 0.01` かつ `convergence_score < -0.3`
- RNA_Processing, Mitochondria, **Synaptic**, **Oxidative_Stress**
- **理由**: 早期発火 + 強発散 = 上流破綻の起点

### **Mid介入（PT 0.3-0.7）**
**標的**: `rho_I_PT_mid > 0.3` かつ `q < 0.10`
- Ion_Transport, Calcium_Signaling
- **理由**: mid期での振幅増大 = 拡散期

### **Late救済（PT 0.7-1.0）**
**標的**: `late_within_phi_fraction > 0.3`（収束可能性あり）
- Complement, Epigenetic, Protein_Homeostasis
- **理由**: late期での部分収束 = 緩衝ノード

---

## 📁 ファイル構成

```
cell_module_timeline/
├── Ex_L5_VAT1L_EYA4_timeline_summary.csv     ★
├── Ex_L5_VAT1L_EYA4_panel.png                ★
├── Ex_L5_VAT1L_EYA4_summary_heatmap.png      ★
│
├── In_PV_PVALB_CEMIP_timeline_summary.csv    ★
├── In_PV_PVALB_CEMIP_panel.png               ★
├── In_PV_PVALB_CEMIP_summary_heatmap.png     ★
│
├── Glia_Astro_GFAP_neg_timeline_summary.csv  ★
├── Glia_Astro_GFAP_neg_panel.png             ★
├── Glia_Astro_GFAP_neg_summary_heatmap.png   ★
│
├── Glia_OPC_timeline_summary.csv             ★
├── Glia_OPC_panel.png                        ★
├── Glia_OPC_summary_heatmap.png              ★
│
├── focus4_overview.png                        ★
│
└── {他39タイプ}_timeline_summary.csv
```

---

## 🔧 次のステップ（オプション）

### **1. 追加比較図**
- VAT1L vs PCP4 vs TLE4（Exタイプ比較）
- 5HT3aR vs SOM vs PV（Inhibitory比較）
- Astro vs OPC vs Oligo vs Micro（Glia比較）

### **2. 介入候補表の最終化**
Timeline要約CSVから、以下の条件で抽出：
```python
# Tier 1（最優先）
activation_time < 0.01 AND convergence_score < -0.3 AND cell_type == 'VAT1L'

# Tier 2（早期介入）
activation_time < 0.05 AND slope_abs_delta_phi > 0.01

# Tier 3（Late救済）
late_within_phi_fraction > 0.3 AND convergence_score > 0.5
```

### **3. 時期×効果×ドラッガブルの統合表**
```
module | cell_type | phase | effect_size | druggability_score | intervention
```

---

**生成日時**: 2025-10-06
**ツール**: visualize_cell_timelines.py
**入力**: module_trajectories/ + phase_correlations/ + module_convergence_scores.csv
**出力**: 52ファイル（43 CSV + 9 PNG）

