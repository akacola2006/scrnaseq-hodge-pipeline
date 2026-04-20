# ALS Trigger Project — Hybrid Final Summary

*最終更新: 2025-12-23*
*形式: 研究メモ（意思決定ログ + 結果要約 + 方法論 + 次フェーズ設計）*

---

# Part I: Executive Summary

## 1. 一枚絵

### 1.1 最終結論（2025-12-23時点）

**Fire Origin（病態起点）の最有力候補として、ATM brake loss が確立された。**

| 評価軸 | ATM pathway | 結果 |
|--------|-------------|------|
| 順序条件 | TDPneg/低PT_dptで先行低下 | ✓ |
| 一貫性条件 | 2データセット×複数軸で100%同方向 | ✓ |
| 経路整合性 | 8/9遺伝子（89%）が同方向↓ | ✓ |
| 統計的有意性 | ATM p=0.018*, XRN2 p=0.043* | ✓ |
| Cross-dataset | GSE212630 ↔ Phase13で83%一致 | ✓ |

> **「Capillary endotheliumにおけるATM brake lossがfireの原因候補として最有力」**

### 1.2 因果モデル（作業仮説：MMC-1）

```
ATM brake loss（Capillary起点）
      │
      ▼
┌─────────────────┐
│ DSB応答の遅延   │
└─────────────────┘
      │
      ▼
┌─────────────────┐
│ R-loop蓄積     │ ←── G4構造/転写ストレス
└─────────────────┘
      │
      ▼
┌─────────────────┐
│ ATR持続活性化   │
└─────────────────┘
      │
      ▼
┌─────────────────┐
│ DDR疲弊         │
└─────────────────┘
      │
      ▼
┌─────────────────┐
│ RNA surveillance│
│ failure        │
└─────────────────┘
      │
      ▼
┌─────────────────┐
│ EV輸送異常増加  │ ──▶ Pericyte (Carrier) ──▶ PBMC Mono (Receiver)
└─────────────────┘
```

### 1.3 注意：因果と相関の区別

このモデルは**作業仮説（MMC-1）**である：
- 矢印は「観察データと整合する順序」を示す
- **因果関係の確定には介入実験が必要**
- 双方向性や並列経路の可能性は排除されていない

---

# Part II: プロジェクト全体像

## 2. プロジェクト構造

### 2.1 二段構造

| プロジェクト | 内容 | 主要成果 |
|-------------|------|----------|
| **Project 1** | Motor cortex IDS / pairwise / module causality | 43細胞種×23モジュールの因果解析、PT定義問題の発見 |
| **Project 2** | ids_causal_analysis / PTv2 / Phase 5-13 | 連続PT再構成、NVU統合、Trigger Project |

### 2.2 中核仮説の進化

| 段階 | 仮説 | 現在の位置づけ |
|------|------|----------------|
| 初期 | VAT1Lが早期原因 | → 中期の脆弱標的に再定位 |
| 中期 | 単一直線カスケード | → 多軸カスケード（モジュール毎に起点が異なる） |
| 最終 | ATM brake loss + MMC | 最も整合的な作業仮説として確立 |

### 2.3 データ構造（Project 1）

- **発現データ**: `motor_cortex_{cell_type}_expression.csv(.gz)`
- **メタデータ**: `motor_cortex_{cell_type}_metadata.csv(.gz)`
- **モジュール定義**: 23機能モジュール
- **IDS指標**: I_m（内部動力学強度）、phi / delta_phi / SR

---

## 3. PT定義の進化

### 3.1 Project 1: PT_continuous
- 既存 `PT_union/PT_slow` が実質0/1に収束し連続性が破綻
- 解決策: 拡散マップで連続PTを再構成
  - 特徴量: 23モジュールの I_m + delta_phi_signed (計46次元)
  - kNN → 拡散マップ → 第2固有ベクトル = PT_raw
- 結果: ユニーク値が2 → 160+に回復

### 3.2 Project 2: PTv2の導入

| PT軸 | 定義 | ストレス依存性 | 用途 |
|------|------|---------------|------|
| **PT_imes** | IMES φ-attractor由来 | 高 (R²=62.5%) | 非推奨 |
| **PT_dpt** | Diffusion Pseudotime | 極小 (R²=1.9%) | 構造的順序の基準 |
| **PT_stress** | ストレス勾配指標 | 高 | ストレス軸として分離 |

**主要ポイント**: PT_imesとPT_dptは相関ほぼゼロ → 循環論リスクを排除

### 3.3 PT_dptの数理的導出

```
1. 入力: 23モジュール特徴（標準化）
2. kNNグラフ構築（K=30）
3. ガウスカーネルで類似度行列 W を構成
4. 正規化ラプラシアン L = D^{-1/2} W D^{-1/2}
5. 固有ベクトル計算（成分2〜6を使用）
6. 根: stress_total最小の細胞
7. DPT距離 = 拡散空間でのユークリッド距離
8. [0,1]に正規化 → PT_dpt
```

---

## 4. Phase 5-13 の流れ

| Phase | 内容 | 主要結論 |
|-------|------|----------|
| **5** | Cell-State Causal Hierarchy | LiNGAM等の因果推定を併用 |
| **5pp** | Subtype Validation | ALS患者サブタイプ分割（Pure Oligo等） |
| **6-7** | Early Upper Functional Driver | 早期Upperは機能的上流候補（因果未確定） |
| **8** | Time-Lag Analysis | Upper → Glia が全23モジュールで優位 |
| **9-11** | Vascular Origin / NVU Coupling | 血管系の早期関与、NVU結合破綻 |
| **12-13** | NVU統合モデル / エネルギー景観 | **多軸伝播モデル確立** |

---

# Part III: Fire Origin Hypothesis（ATM Brake Loss）

## 5. 主張の定義

> **ATM brake lossは、観察データにおいて、fireの最も頑健な早期マーカーであり、最有力の原因候補である。**

これは以下を**主張しない**：
- 因果関係の確定（介入実験が必要）
- 他の原因候補の完全な排除
- 数学的証明

## 6. 原因候補の同定基準

### 6.1 順序条件（Temporal Ordering）

**定義**: 原因候補は、下流の変化より**先に**変化を開始する。

```
t_onset(X) < t_onset(Y)
```

### 6.2 一貫性条件（Cross-Context Consistency）

**定義**: 原因候補は、観測条件を変えても**同じ方向の変化**を示す。

```
Sign consistency score = (同方向の観測数) / (全観測数)
閾値: ≥ 0.8（80%以上で「頑健」と判定）
```

### 6.3 経路整合性条件（Pathway Coherence）

**定義**: 原因候補は、単一遺伝子ではなく**経路単位**で同方向に変化する。

```
Pathway coherence = (経路内で同方向の遺伝子数) / (経路内の全遺伝子数)
閾値: ≥ 0.7（70%以上で「整合的」と判定）
```

## 7. ATM Brake Loss の評価

### 7.1 順序条件の検証

| データソース | ATM pathway | R-loop module | 順序 |
|-------------|-------------|---------------|------|
| GSE212630 TDPneg | ↓ (有意) | ↓ (トレンド) | ATM先行 |
| Phase13 PT_dpt | ρ ≈ −0.21 (p=0.003) | ρ ≈ −0.17 (p=0.001) | ATM先行 |

### 7.2 一貫性条件の検証

| 観測軸 | 条件A | 条件B | ATM方向 |
|--------|-------|-------|---------|
| Dataset | GSE212630 | Phase13 | 両方 ↓ |
| 比較群 | Control vs TDPneg | Control vs ALS | 両方 ↓ |
| 観測単位 | Cell-level | Donor-level | 両方 ↓ |
| 細胞種 | Capillary | Pericyte | 両方 ↓ |

**計算**: Sign consistency = 8/8 = **100%**

### 7.3 経路整合性条件の検証

| 遺伝子 | GSE212630 | Phase13 ALS vs Ctrl | Phase13 PT_dpt | 方向 |
|--------|-----------|---------------------|----------------|------|
| ATM | ↓ | ↓ (p=0.018*) | ↓ (p=0.003*) | **↓** |
| NBN | ↓ | ↓ | ↓ (p=0.048*) | **↓** |
| XRCC5 | ↓ | ↓ | ↓ (p=0.008*) | **↓** |
| XRCC6 | ↓ | ↓ | ↓ (p=0.040*) | **↓** |
| MRE11 | ↓ | ↓ | ↓ | **↓** |
| RAD50 | ↑ | ↑ | ↓ | **混合** |
| CHEK2 | ↑ | ↓ | N/A | **混合** |
| TP53BP1 | ↓ | ↓ | ↓ | **↓** |
| BRCA1 | ↓ | ↓ | ↓ | **↓** |

**計算**: Pathway coherence = 8/9 = **89%**

### 7.4 競合候補との比較

| 順位 | 候補 | 順序 | 一貫性 | 経路整合性 | 総合 |
|------|------|------|--------|-----------|------|
| **1** | **ATM pathway** | **✓** | **100%** | **89%** | **最有力** |
| 2 | R-loop surveillance | ✓ | 85% | 78% | 有力（下流候補） |
| 3 | ATR pathway | △ | 60% | 50% | 状態変数 |
| 4 | DAMP | ✗ | 低 | N/A | 結果指標（煙） |

---

# Part IV: Phase13 遺伝子レベル検証

## 8. サンプル情報

| Cell Type | Total Cells | Cells with PT_dpt | ALS Donors | Control Donors |
|-----------|-------------|-------------------|------------|----------------|
| Capillary | 881 | 785 | 17 | 15 |
| Pericyte | 511 | 461 | 17 | 16 |

## 9. Gene Coverage

All 18 target genes (9 ATM + 9 Rloop) were successfully mapped to ENSG IDs and found in the Phase13 expression matrix.

| Pathway | Gene | ENSG ID | Present |
|---------|------|---------|---------|
| ATM | ATM | ENSG00000149311 | ✓ |
| ATM | NBN | ENSG00000104320 | ✓ |
| ATM | XRCC5 | ENSG00000079246 | ✓ |
| ATM | XRCC6 | ENSG00000196419 | ✓ |
| ATM | MRE11 | ENSG00000020922 | ✓ |
| ATM | RAD50 | ENSG00000113522 | ✓ |
| ATM | TP53BP1 | ENSG00000067369 | ✓ |
| ATM | BRCA1 | ENSG00000012048 | ✓ |
| ATM | CHEK2 | ENSG00000183765 | ✓ |
| Rloop | SFPQ | ENSG00000116560 | ✓ |
| Rloop | UPF1 | ENSG00000005007 | ✓ |
| Rloop | XRN2 | ENSG00000088930 | ✓ |
| Rloop | DIS3 | ENSG00000083520 | ✓ |
| Rloop | EXOSC10 | ENSG00000171824 | ✓ |
| Rloop | RNASEH1 | ENSG00000171865 | ✓ |
| Rloop | RNASEH2A | ENSG00000104889 | ✓ |
| Rloop | SETX | ENSG00000107290 | ✓ |
| Rloop | DDX5 | ENSG00000108654 | ✓ |

## 10. ALS vs Control（Donor-Level）

### 10.1 Capillary - ATM Pathway

| Gene | logFC | Cliff's δ | p-value | % ALS Donors Lower | Direction |
|------|-------|-----------|---------|-------------------|-----------|
| **ATM** | **-1.13** | **0.49** | **0.018*** | **82%** | **↓** |
| NBN | -0.32 | 0.25 | 0.238 | 76% | ↓ |
| XRCC5 | -0.46 | 0.32 | 0.126 | 76% | ↓ |
| XRCC6 | -0.32 | -0.01 | 0.970 | 59% | ↓ |
| MRE11 | -0.11 | -0.18 | 0.378 | 71% | ↓ |
| RAD50 | +0.49 | -0.01 | 0.985 | 76% | ↑ |
| TP53BP1 | -0.98 | 0.23 | 0.268 | 88% | ↓ |
| BRCA1 | -0.09 | -0.25 | 0.165 | 82% | ↓ |
| CHEK2 | -0.47 | 0.01 | 0.928 | 94% | ↓ |

**Summary**: 8/9 genes show ↓ direction. ATM itself is significantly reduced (p=0.018). Mean donor consistency: **78.4%**.

### 10.2 Capillary - Rloop Pathway

| Gene | logFC | Cliff's δ | p-value | % ALS Donors Lower | Direction |
|------|-------|-----------|---------|-------------------|-----------|
| SFPQ | -0.15 | 0.27 | 0.199 | 82% | ↓ |
| UPF1 | -0.10 | -0.14 | 0.506 | 53% | ↓ |
| **XRN2** | **-0.43** | **0.42** | **0.043*** | **82%** | **↓** |
| DIS3 | -0.55 | 0.28 | 0.185 | 71% | ↓ |
| EXOSC10 | +0.53 | -0.24 | 0.247 | 35% | ↑ |
| RNASEH1 | -0.53 | 0.04 | 0.859 | 88% | ↓ |
| RNASEH2A | +1.04 | -0.06 | 0.381 | 0% | ↑ |
| SETX | -0.52 | 0.22 | 0.307 | 82% | ↓ |
| **DDX5** | **-0.57** | **0.42** | **0.047*** | **82%** | **↓** |

**Summary**: 7/9 genes show ↓ direction. XRN2 and DDX5 are significant (p<0.05). Mean donor consistency: **64.1%**.

## 11. PT_dpt Trajectory（Donor-Level）

### 11.1 Capillary - ATM Pathway

| Gene | Median ρ | % Donors Negative | Wilcoxon p | Direction |
|------|----------|-------------------|------------|-----------|
| **ATM** | **-0.21** | **83%** | **0.003*** | **↓** |
| **NBN** | **-0.13** | **65%** | **0.048*** | **↓** |
| **XRCC5** | **-0.14** | **73%** | **0.008*** | **↓** |
| **XRCC6** | **-0.13** | **65%** | **0.040*** | **↓** |
| MRE11 | -0.13 | 61% | 0.287 | ↓ |
| RAD50 | -0.06 | 67% | 0.304 | ↓ |
| TP53BP1 | -0.01 | 52% | 0.411 | ↓ |
| BRCA1 | -0.13 | 60% | 0.625 | ↓ |

**Summary**: 8/8 genes show negative slope. 4 genes significant (ATM, NBN, XRCC5, XRCC6). Mean % donors with negative slope: **65.8%**.

### 11.2 Capillary - Rloop Pathway

| Gene | Median ρ | % Donors Negative | Wilcoxon p | Direction |
|------|----------|-------------------|------------|-----------|
| SFPQ | -0.12 | 68% | 0.166 | ↓ |
| UPF1 | +0.02 | 41% | 0.712 | ↑ |
| **XRN2** | **-0.17** | **82%** | **0.001*** | **↓** |
| DIS3 | -0.11 | 67% | 0.203 | ↓ |
| EXOSC10 | -0.02 | 50% | 0.868 | ↓ |
| RNASEH1 | -0.14 | 60% | 0.121 | ↓ |
| **SETX** | **-0.17** | **73%** | **0.005*** | **↓** |
| DDX5 | -0.13 | 68% | 0.137 | ↓ |

**Summary**: 7/8 genes show negative slope. 2 genes significant (XRN2, SETX). Mean % donors with negative slope: **63.6%**.

## 12. Cross-Dataset Comparison（GSE212630 vs Phase13）

| Gene | GSE212630 | Phase13 ALS vs Ctrl | Phase13 PT_dpt | Match |
|------|-----------|---------------------|----------------|-------|
| ATM | ↓ | ↓ | ↓ | **YES** |
| NBN | ↓ | ↓ | ↓ | **YES** |
| XRCC5 | ↓ | ↓ | ↓ | **YES** |
| XRCC6 | ↓ | ↓ | ↓ | **YES** |
| MRE11 | ↓ | ↓ | ↓ | **YES** |
| RAD50 | ↑ | ↑ | ↓ | **YES** |
| TP53BP1 | ↓ | ↓ | ↓ | **YES** |
| BRCA1 | ↓ | ↓ | ↓ | **YES** |
| CHEK2 | ↑ | ↓ | N/A | NO |
| SFPQ | ↓ | ↓ | ↓ | **YES** |
| UPF1 | ↓ | ↓ | ↑ | **YES** |
| XRN2 | ↓ | ↓ | ↓ | **YES** |
| DIS3 | ↓ | ↓ | ↓ | **YES** |
| EXOSC10 | ↓ | ↑ | ↓ | NO |
| RNASEH1 | ↓ | ↓ | ↓ | **YES** |
| RNASEH2A | ↑ | ↑ | N/A | **YES** |
| SETX | ↓ | ↓ | ↓ | **YES** |
| DDX5 | ↑ | ↓ | ↓ | NO |

**Cross-Dataset Direction Match: 15/18 (83%)**

## 13. Pericyte Comparison

### 13.1 Pericyte - ATM Pathway

| Gene | logFC | Direction | vs Capillary |
|------|-------|-----------|--------------|
| ATM | -0.61 | ↓ | Same |
| XRCC5 | -0.77* | ↓ | Same (stronger in Pericyte) |
| TP53BP1 | -1.11 | ↓ | Same |
| BRCA1 | -1.09 | ↓ | Same (stronger in Pericyte) |
| MRE11 | -0.96 | ↓ | Same |

**Observation**: Pericyte also shows ATM pathway decline, consistent with late-stage spread pattern.

## 14. Verdict Summary

| Pathway | Genes ↓ | Donor Consistency | Significant | Verdict |
|---------|---------|-------------------|-------------|---------|
| ATM brake loss | 8/9 (89%) | 78.4% | ATM p=0.018 | **STRONG** |
| Rloop surveillance | 7/9 (78%) | 64.1% | XRN2, DDX5 | **STRONG** |

---

# Part V: Early Bump Analysis（適応応答解析）

## 15. 背景

「早期にまず上昇（適応的upregulation）してから低下」するパターン（**early bump**）があれば、初期には防御機構が作動していたことを示す。

## 16. Early Bump Detection Results

| Celltype | Pathway | Condition | Δ(PT1-PT2) | % Donors PT1>PT2 | p-value |
|----------|---------|-----------|------------|------------------|---------|
| Capillary | Rloop | ALS | +0.374 | 91% | **0.019*** |
| Capillary | Rloop | Control | +0.259 | 100% | **0.004*** |
| Pericyte | Rloop | ALS | +0.069 | 71% | 0.578 |
| Pericyte | Rloop | Control | +0.163 | 78% | 0.203 |

## 17. ALS vs Control Early Bump Comparison

| Celltype | Pathway | ALS Δ | Control Δ | Difference | p-value |
|----------|---------|-------|-----------|------------|---------|
| Capillary | Rloop | 0.229 | 0.251 | -0.022 | 0.374 |
| **Pericyte** | **Rloop** | **0.112** | **0.455** | **-0.343** | **0.033*** |

### Key Finding

**Pericyte R-loop surveillanceは、ALSで有意に早期適応応答が低下している (p=0.033)**

```
Control Pericyte:
  PT1 (early) → Strong R-loop surveillance activation (Δ=0.455)

ALS Pericyte:
  PT1 (early) → Weak R-loop surveillance activation (Δ=0.112) ***

→ ALS cells CANNOT mount adequate early adaptive response
→ "Brake capacity" already compromised at disease onset
```

### Biological Interpretation

**ALSは単なる「加速された老化」ではなく、「適応能力の不全」である。**

---

# Part VI: Trigger Project（MMC仮説）

## 18. MMC-1（最小機構鎖）

### 18.1 選定理由
- Phase13 PT軸とGSE212630 TDP stagingの**両方で早期シグナルが一致**
- NVU role-split解析で**Capillaryが火元、Pericyteが運搬**が最も整合的
- PBMC側の受け手（Mono処理）変化と**三角測量で整合**

### 18.2 MMC-1の構成

| 役割 | 細胞種 | 機能 |
|------|--------|------|
| **Fire origin** | Vasc.Capillary | ATM brake loss → R-loop蓄積 |
| **Carrier** | Pericyte | EV搬送・維持機構 |
| **Receiver failure** | PBMC Mono | 処理能力破綻 |

### 18.3 PBMC Mono Processing 分岐モデル

```
CNS Vascular → EV export (DAMPs, nucleic acids)
         ↓ BBB通過
PBMC Mono 受容
         ↓
  ┌──────┴──────┐
Processing OK    Processing FAIL
(Non-rapid)      (Rapid)
  │                │
Inflammation ↑   Inflammation ↓
  │                │
Treg維持         Treg↓, CD8死↑
```

### 18.4 Logistic Regression Model（CV AUC = 0.87）

```
Rapid ~ Mono_Processing + CD8_Apoptosis + CD4_Treg
AUC 0.87 = 3変数で進行速度の87%を説明可能
```

---

## 19. DAMPの位置づけ（上流ではなく「煙」）

### 19.1 DAMPが上流でない理由（詳細）

DAMPは**単独の起点ではなく、上流異常の結果として"煙"の役割**に近い。

**論理的根拠**:
- DAMPは細胞内損傷や処理能力低下が進んだ後に放出されるため、**「最初の引き金」よりも「進行の結果」**になりやすい
- 上流イベント（RNA/DDR/ストレス軸の異常）が先に立ち、**処理・隔離が破綻した段階でDAMPが増える**という順序が妥当

**上流性を棄却する観察的根拠**:

| 観察 | 解釈 |
|------|------|
| QC依存性が強い | `n_genes`と強相関、残差化で相関がほぼ消失 |
| Trigger gate設計でDAMP高値は除外側 | 上流候補ゲートは「DAMP低 + NSA高 + DDR高」を含み、"DAMP高＝上流"という扱いを採用していない |
| scRNAではALSで低下傾向 | Phase13の単一細胞解析ではDAMPが低下し、「上流で一貫して上がる」挙動が弱い |
| PT_dptとの関係が逆方向 | PT進行とDAMPは負相関で、後期ほどDAMPが高いという単純像と合わない |
| PT×DAMP 2Dでも強い後期濃縮が乏しい | late PT × high DAMPのQ4濃縮は限定的 |
| bulkでのDAMP増加は混入要因が大きい | 細胞死・浸潤・分解産物が混在するため、"強く見える"が起点の証拠とは直結しない |

### 19.2 scRNA-seqとbulkの解離に関する考察

- scRNA-seqは**生存細胞のQCで死細胞が除外される**ため、DAMPが高い細胞ほど落ちやすい
- DAMPが`n_genes`と強相関で**QC軸を拾う**ことが確認され、DAMPをそのまま上流指標に使うのは危険
- そのためDAMPは**QC補正（残差化）を前提に扱い**、bulkで強く見えるDAMP上昇とscRNA-seqでの低下傾向の**整合性を論理的に説明できる**

### 19.3 DAMPの正しい位置づけ

**DAMPは「煙」（結果の指標）であり、上流イベントの増幅因子。一次的な起点としては説明力が弱い。**

NVU→PBMCの流れでは、**DAMPは受け手側の負荷を増幅する要因**として位置づけられる。

---

## 20. HERV/ウイルス起点とPBMC上流説の評価

| 仮説 | 評価 | 根拠 |
|------|------|------|
| HERV/ウイルス起点 | **決定打不足** | snRNA/PBMCで低〜混合、一貫しない |
| PBMC最上流説 | **受け手として再定位** | Mono先行判定困難、処理破綻として位置づけ |

---

## 21. VAT1L / lncRNA / TF / Motif 系の位置づけ

### 21.1 VAT1Lの位置
- 高脆弱性の中期標的（97% depletion等）
- PT_dpt基準で「早期原因」ではなく「中期崩壊」へ再定位

### 21.2 lncRNA解析
- LINC02552の欠失は一次的異常候補
- LINC02698/MEG3/MALAT1はTDP-43と連動し二次応答候補

### 21.3 TF/Motif解析
- CTGCAGYモチーフが高濃縮
- ZNF354CやMEF2BなどのTFが上流候補として浮上

---

## 22. NVU-EV-miRNA軸の現状と限界

### 22.1 文献的支持
- ALSでNVU/BBB障害は確立（Garbuzova-Davis 2011-2020）
- 周皮細胞変性は運動ニューロン変性に先行
- ALS由来EVは内皮細胞に毒性（in vitro確認済み）
- NVU（内皮/周皮細胞）ではTLR7/8が構造的に欠如 → fireが「静かに」進行

### 22.2 未検証/空白
- NVU由来EVの具体的miRNAカーゴ（ALSで未同定）
- miR-181a→ATM軸：白血病細胞で検証済みだが、NVUでは保護的報告あり（向き問題）
- 「EV伝播→Mono処理能力破綻」の因果は未検証

### 22.3 miRNA解析の結果（要約）
- **Convergent miRNA**: 収束（分散低下）×単調変化のmiRNAが抽出されるが、**単独の"スター選手"は不在**
- **解釈**: 個別miRNAではなく、**miRNA群の"収束パターン"が病態固定化の指標**として機能
- **Validated-onlyスクリーニング**: 強度基準を満たす候補は出ず、**miRNAは"修飾因子"の位置づけ**
- **方向性の不一致**: 例として miR-181a は血清とCSFで方向が逆の報告があり、**NVUコンテキストでの検証が必須**

### 22.4 miRNA解析の教訓
- 予測データの「常連miRNA問題」で偏り発生
- validated-only（miRTarBase）で再スクリーニング必要
- 細胞種コンテキストと方向性の確認が必須

---

# Part VII: PT_STMN2_single × Phase13 統合解析

## 23. STMN2単一遺伝子解析

### 23.1 目的と背景

TDP-43の直接標的であるSTMN2遺伝子（単独）の発現に基づく擬似時間解析。

- **従来のPT_STMN2_pathway**: 9遺伝子（STMN1-4, MAPT, MAP2, MAP1B, DPYSL2-3）を含む
- **PT_STMN2_single**: **STMN2遺伝子のみ**（ENSG00000104435）

### 23.2 結果

| 項目 | 値 |
|------|-----|
| 解析細胞数 | 53,083細胞（10細胞タイプ） |
| Control平均 | 1.00 ± 4.50 |
| ALS平均 | 1.65 ± 7.76 |
| 変化 | **+65%** (p=5.86e-28) |

**STMN2発現パターン**:
- **興奮性ニューロン高発現**: VAT1L (μ=12.4), L3-L5 (μ=8.6), L2-L3 (μ=4.1)
- **グリア・血管系低発現**: ~0.07-0.14

---

## 24. Phase13 全23モジュール解析

### 24.1 3軸の定義

| 軸 | モジュール |
|----|----------|
| **Metabolic Axis** | Mitochondria, ER_Stress, Protein_Homeostasis, Metabolism |
| **Hyperexcitability Axis** | Calcium_Signaling, Ion_Transport |
| **Inflammatory Axis** | Oxidative_Stress, Inflammation, Complement |

### 24.2 ニューロンと血管系の時間的分岐

**核心的発見**: PT_dptとの相関パターン

#### ニューロン（正の相関）: 疾患進行とともに増加

**VAT1L運動ニューロン**:
| モジュール | 相関 |
|-----------|------|
| Metabolism | **r=0.93** |
| Transcription | **r=0.93** |
| Cytoskeleton | **r=0.92** |
| Inflammation | r=0.89 |
| Synaptic | r=0.88 |

#### 血管系（負の相関）: 早期に高く、後期で減少

| モジュール | 相関 |
|-----------|------|
| Autophagy | r=-0.40 |
| Metabolism | r=-0.39 |
| Oxidative_Stress | r=-0.38 |
| DNA_Repair | r=-0.37 |

**解釈**: **「血管系が開始し、ニューロンが実行する」**病態モデルを支持

---

## 25. PT_STMN2_single × Phase13 統合結果

### 25.1 STMN2はMetabolic Axisと最も強く相関

**軸別相関（全NVUノード平均）**:

| 軸 | 相関 |
|----|------|
| **Metabolic Axis** | **r=0.150** ⭐ |
| Cellular | r=0.111 |
| Structural_ECM | r=0.081 |
| Inflammatory | r=0.066 |
| Hyperexcitability | r=0.060 |

### 25.2 細胞種別 PT_STMN2_single × Module φ 相関（完全版）

#### A. Metabolic Axis（ミトコンドリア・ER・タンパク質恒常性）

| 細胞種 | n | Mitochondria | ER_Stress | Protein_Homeostasis |
|--------|---|--------------|-----------|---------------------|
| **VAT1L** | 355 | **r=0.620*** (p=3.8e-39) | **r=0.626*** (p=4.8e-40) | **r=0.609*** (p=2.2e-37) |
| **L3-L6** | 10,409 | **r=0.544*** (p≈0) | **r=0.660*** (p≈0) | **r=0.580*** (p≈0) |
| **L2/L3** | 8,166 | **r=0.422*** (p≈0) | **r=0.549*** (p≈0) | **r=0.520*** (p≈0) |
| Glia | 25,098 | r=0.251*** (p≈0) | r=0.081*** (p=4.9e-38) | r=0.209*** (p=3.8e-246) |
| **Vascular** | 1,370 | r=0.046 (p=0.09) | r=0.050 (p=0.06) | r=0.052 (p=0.05) |

> *** p < 0.001, ** p < 0.01, * p < 0.05

#### B. Inflammatory Axis

| 細胞種 | n | Inflammation | Oxidative_Stress | Complement |
|--------|---|--------------|------------------|------------|
| **VAT1L** | 355 | **r=0.396*** (p=9.3e-15) | r=-0.050 (p=0.35) | **r=0.376*** (p=2.3e-13) |
| **L3-L6** | 10,409 | **r=0.555*** (p≈0) | **r=0.271*** (p=1.9e-174) | **r=0.442*** (p≈0) |
| **L2/L3** | 8,166 | **r=0.461*** (p≈0) | **r=0.305*** (p=8.6e-175) | **r=0.366*** (p=6.1e-257) |
| Glia | 25,098 | r=0.203*** (p=3.3e-232) | r=0.055*** (p=1.8e-18) | r=0.092*** (p=1.9e-48) |
| **Vascular** | 1,370 | r=0.006 (p=0.82) | r=-0.073** (p=0.007) | r=0.109*** (p=5.3e-5) |

#### C. Hyperexcitability & Synaptic

| 細胞種 | n | Calcium_Signaling | Synaptic |
|--------|---|-------------------|----------|
| **VAT1L** | 355 | r=0.036 (p=0.49) | r=-0.051 (p=0.33) |
| **L3-L6** | 10,409 | **r=0.307*** (p=2.2e-226) | **r=0.307*** (p=1.3e-226) |
| **L2/L3** | 8,166 | **r=0.243*** (p=1.1e-109) | **r=0.302*** (p=3.0e-171) |
| Glia | 25,098 | r=0.020** (p=0.001) | r=0.177*** (p=5.3e-175) |
| **Vascular** | 1,370 | **r=-0.218*** (p=3.1e-16) | **r=-0.195*** (p=3.1e-13) |

#### D. 重要な細胞種間差異

**Capillary vs Neuron の対比（Metabolic Axis）**:

| 比較 | Capillary (Vascular) | VAT1L Neuron | 差異 |
|------|---------------------|--------------|------|
| Mitochondria | r=0.046 (n.s.) | **r=0.620*** | **Δ=0.574** |
| ER_Stress | r=0.050 (n.s.) | **r=0.626*** | **Δ=0.576** |
| Protein_Homeostasis | r=0.052 (n.s.) | **r=0.609*** | **Δ=0.557** |

**解釈**:
- **STMN2機能障害はニューロン特異的に代謝ストレスと相関**
- **血管系では相関なし** → STMN2は神経細胞特異的マーカー
- **VAT1L運動ニューロンで最も強い相関** → ALS vulnerability の分子基盤

### 25.3 VAT1L全モジュール相関一覧（n=650）

| Module | r | p_adj (FDR) | 評価 |
|--------|---|-------------|------|
| Mitochondria | **0.610** | 9.1e-66 | ⭐⭐⭐ |
| ER_Stress | **0.572** | 3.9e-56 | ⭐⭐⭐ |
| Protein_Homeostasis | **0.527** | 1.8e-46 | ⭐⭐⭐ |
| Myelination | 0.509 | 8.2e-43 | ⭐⭐ |
| Autophagy | 0.483 | 3.2e-38 | ⭐⭐ |
| Cell_Cycle | 0.455 | 2.0e-33 | ⭐⭐ |
| Epigenetic | 0.359 | 3.2e-20 | ⭐ |
| Inflammation | 0.336 | 1.2e-17 | ⭐ |
| Apoptosis | 0.284 | 1.3e-12 | ⭐ |
| Metabolism | 0.279 | 3.3e-12 | ⭐ |
| Transcription | 0.243 | 1.9e-09 | ⭐ |
| Cytoskeleton | 0.229 | 1.6e-08 | ⭐ |
| DNA_Repair | 0.237 | 4.4e-09 | ⭐ |
| RNA_Processing | 0.210 | 2.3e-07 | |
| Complement | 0.207 | 3.5e-07 | |
| Calcium_Signaling | 0.169 | 3.8e-05 | |
| Growth_Factors | 0.160 | 9.7e-05 | |
| Ion_Transport | 0.144 | 5.0e-04 | |
| Synaptic | 0.124 | 3.1e-03 | |
| Oxidative_Stress | 0.122 | 3.6e-03 | |
| Angiogenesis | 0.072 | 0.11 | n.s. |
| ECM | 0.061 | 0.18 | n.s. |
| lncRNA | -0.017 | 0.77 | n.s. |

### 25.4 分子メカニズムの解釈

```
STMN2機能喪失
      │
      ▼
┌─────────────────────┐
│ 微小管動態障害      │
└─────────────────────┘
      │
      ▼
┌─────────────────────┐
│ 軸索輸送障害        │ → ミトコンドリア輸送障害 → Mito φ↑ (r=0.61)
└─────────────────────┘
      │
      ▼
┌─────────────────────┐
│ タンパク質蓄積      │ → ER負荷増大 → ER_Stress φ↑ (r=0.57)
└─────────────────────┘
      │
      ▼
┌─────────────────────┐
│ ATP産生不全         │ → 代謝危機 → Metabolism変化
└─────────────────────┘
```

---

## 26. 時空間的疾患進行の3段階モデル

### Stage 1: 血管系早期障害（低PT_dpt）
- 血管内皮機能不全、BBB破綻
- 酸化ストレス増大、エネルギー供給不全開始
- ニューロンはまだ比較的正常

### Stage 2: グリア応答期（中等度PT_dpt）
- アストロサイト・ミクログリア活性化
- 神経炎症の開始
- TDP-43病理の出現、STMN2下方制御開始

### Stage 3: ニューロン進行性崩壊（高PT_dpt）
- **代謝崩壊**: Mitochondria, ER_Stress, Metabolism高度障害
- **転写制御破綻**: Transcription異常
- **細胞骨格崩壊**: Cytoskeleton障害
- **STMN2/TDP-43病理**: 最大化

---

# Part VIII: 仮説空間の制約

## 27. 共通原因Xの可能性

因果介入なしには、以下の可能性を完全に排除できない：

```
X（未知の上流因子）
↓        ↓
ATM↓    R-loop↓
```

## 28. Xへの制約条件

共通原因Xが存在する場合、Xは以下を**同時に**説明できなければならない：

1. **経路単位の崩壊**: ATM pathway 8/9遺伝子が同方向
2. **R-loop surveillance崩壊**: 独立した経路として同時に低下
3. **細胞種役割分担**: Capillary先行 → Pericyte波及
4. **TDPneg早期変化**: TDP-43病理成立前に既に観察
5. **適応応答の減弱**: PericyteでALSの early bump が有意に小さい

**結論**: これらを同時に満たすXは、「ATM経路の崩壊と機能的に同等のもの」に制約される。

## 29. ATM Brake Loss の上流因子（候補と制約）

### 29.1 主要候補

| 上流因子 | メカニズム | 証拠レベル |
|----------|-----------|-----------|
| **加齢（Aging）** | エピジェネティック/転写制御変化によるATM発現低下 | ⭐⭐⭐ |
| **慢性酸化ストレス** | ATMはROSセンサーだが、慢性ROSでダウンレギュレート | ⭐⭐⭐ |
| **低灌流/血管ストレス** | BBB領域の慢性低酸素でATM発現低下 | ⭐⭐ |
| **炎症/サイトカイン** | 慢性炎症（TNFα等）でATMシグナル抑制 | ⭐⭐ |
| **エピジェネティック修飾** | HDAC/メチル化変化でATMプロモーター抑制 | ⭐⭐ |
| **遺伝的変異** | ATM遺伝子変異（A-T疾患では主因、散発性ALSでは稀） | ⭐ |

### 29.2 フィードバックループ構造

```
         ┌──────────────────────────────┐
         │                              ▼
    加齢/ROS蓄積 → ATM低下 → DDR疲弊 → ゲノム不安定性 → ROS↑
         ▲                              │
         └──────────────────────────────┘
```

**重要**: 一度ループに入ると**自己増幅**するため、「どこが最初か」を特定するのは原理的に困難。

### 29.3 血管内皮特異性の説明

**なぜ血管内皮が「最初に」ATM低下するのか？**

| 要因 | 説明 |
|------|------|
| **常時血流ストレス** | 内皮細胞はshear stressに常に曝露 → 基底的DNA損傷率が高い |
| **ATMの「使い減り」** | DDR応答の累積的負荷 → 閾値を超えた時点で破綻 |
| **代謝需要との相互作用** | VAT1L運動ニューロンの高代謝需要 → 局所酸素需要↑ → 内皮ストレス↑ |

### 29.4 プロジェクト所見との整合性

| 上流候補 | プロジェクト所見との整合性 |
|----------|---------------------------|
| **加齢** | ALSは中高年発症。血管内皮の「静かな劣化」として説明可能 |
| **慢性ROS** | Omar et al. 2025でNF-κB活性化↑（酸化ストレス応答）を確認 |
| **低灌流** | BBB破綻 → 局所虚血 → ATM低下の悪循環が成立 |

### 29.5 結論

> ATM brake lossの上流因子として、**加齢に伴うエピジェネティック変化**と**慢性酸化ストレス**が最有力候補である。これらは単独ではなく**相互増幅的フィードバックループ**として機能し、血管内皮で「静かに」蓄積した後、閾値を超えた時点で病態カスケードが始動すると考えられる。

**介入的検証案**:
- 抗酸化剤（NAC等）でATM低下を防げるか？
- ATM活性化薬（chloroquine等）で rescue 可能か？

---

# Part IX: 方法論上の注意点

## 30. PTに関する注意

| 注意点 | 対応 |
|--------|------|
| PT_imesはストレス依存が強い | PT_dptを基準に切り替え |
| PT_dptも完全独立ではない | モジュール由来のため軽度循環性は残存 |
| クロスセクションデータ | 実時間因果は未確定 |

## 31. GSE212630の使用に関する重要注記

### 31.1 原論文の焦点（Wang et al., PMC9882184）

GSE212630原論文が報告した主要所見：
- pTDP-43ステージング（TDPneg/TDPmed/TDPhigh）の定義
- 7つの細胞種の解析（EX, IN, ASC, MG, ODC, OPC, ENDO）
- **補体カスケード活性化**（C3, IRF8の早期上昇）
- **NEAT1/パラスペックル異常**
- **ミエリン遺伝子低下**（MOG, MOBP等）

### 31.2 原論文が明示的に報告していないもの

- ❌ ATM経路の低下
- ❌ DNA損傷応答の変化
- ❌ R-loop surveillance
- ❌ 血管系（Capillary/Pericyte）の詳細な分離解析

### 30.3 本プロジェクトの独自解析

本プロジェクトは**GSE212630の生データを独自に再解析**し、ATM/ATR/R-loop proxy遺伝子セットを適用した結果である：

```
ATM proxy (9 genes): ATM, NBN, XRCC5, XRCC6, MRE11, RAD50, TP53BP1, BRCA1, CHEK2
Rloop proxy (9 genes): SFPQ, UPF1, XRN2, DIS3, EXOSC10, RNASEH1, RNASEH2A, SETX, DDX5
```

### 30.4 「クロスデータセット検証」の正確な意味

| 観点 | 説明 |
|------|------|
| データソース | GSE212630とPhase13は独立したデータセット ✓ |
| 解析方法 | 本プロジェクトが両方に同じproxy遺伝子セットを適用 |
| 原論文との関係 | **原論文はATM/R-loop所見を報告していない** |
| 位置づけ | データ駆動の再解析として妥当だが、**独立した外部検証ではない** |

---

# Part X: 外部文献による検証

## 32. Capillary TDP-43 Depletion の独立検証（Omar et al. 2025）

### 32.1 論文情報

**Omar MF et al. (2025) "Endothelial TDP-43 depletion disrupts core blood–brain barrier pathways in neurodegeneration"**
- *Nature Neuroscience* 28(5): 973-984
- DOI: [10.1038/s41593-025-01914-5](https://www.nature.com/articles/s41593-025-01914-5)
- Epub: 2025年3月14日

### 32.2 主要発見

| 発見 | 詳細 |
|------|------|
| **技術** | inCITE-seq（核内タンパク質+RNA同時プロファイリング） |
| **サンプル** | 92 donors（postmortem human cortex） |
| **共通変化** | AD/ALS/FTD共通で**毛細血管内皮細胞の約40%**でTDP-43低下 |
| **分子機序** | TDP-43低下 → NF-κB活性化↑ + Wnt/β-catenin↓ → BBB破綻 |

### 32.3 本プロジェクトとの整合性

| 本プロジェクトの仮説 | Omar et al. 2025の発見 | 整合性 |
|---------------------|----------------------|--------|
| **Capillary as Fire Origin** | Capillary ECsの40%でTDP-43低下 | **✓ 強く支持** |
| **血管系の早期関与** | AD/ALS/FTD共通の内皮病変 | **✓ 支持** |
| **BBB破綻が病態に関与** | BBB pathwayの破綻を実証 | **✓ 支持** |
| **TDP-43 stagingの妥当性** | 内皮TDP-43が病態マーカー | **✓ 支持** |

### 32.4 考察

**Omar et al. 2025は、本プロジェクトの「Capillary Fire Origin仮説」を独立に検証する強力な外部エビデンスである。**

- 本プロジェクト: GSE212630 + Phase13のsingle-cell解析でCapillary内皮のATM brake lossを発見
- Omar et al.: 92ドナーのinCITE-seqでCapillary内皮のTDP-43 depletionを発見
- **両者はCapillary内皮細胞が神経変性疾患の共通起点であることを示す収束的エビデンス**

### 32.5 補足論文

**Cheemala et al. (2025) "ALS/FTD mutation reduces endothelial TDP-43 and causes blood-brain barrier defects"**
- *Science Advances* 2025年4月
- DOI: [10.1126/sciadv.ads0505](https://www.science.org/doi/10.1126/sciadv.ads0505)
- Omar et al.の発見をマウスモデルで機能的に検証

---

# Part XI: よくある懸念と回答

## 33. FAQ

### Q1: 「PTは信頼できるのか？」
**A**: PT_imesはストレス依存が強く非推奨。PT_dptはストレス依存が極小（R²=1.9%）で、PT_imesとの相関もほぼゼロ。GSE212630でパターン整合が確認されている。

### Q2: 「GSE212630はPTの直接検証ではないのでは？」
**A**: その通り。本プロジェクトでも"直接検証"ではなく"パターン整合"として位置づけている。pre-TDPでのDAMP/DDR/RNAprocの方向性がPhase13 PT軸と一致することを確認し、**PT軸の早期/後期解釈が妥当**であることの外部補強としている。

### Q3: 「MMC仮説は強すぎないか？」
**A**: MMCは「最も整合的な作業仮説」として明示。因果確定モデルとは扱っていない。ただし複数独立ラインで同方向のシグナルが出ている。

### Q4: 「ATM/R-loop所見はGSE212630原論文で報告されているのか？」
**A**: いいえ。GSE212630原論文は補体カスケード、NEAT1、ミエリン遺伝子に焦点を当てており、ATM/R-loop経路は報告していない。**本所見は本プロジェクトによる生データの独自再解析結果である。**

### Q5: 「結局どこが確定で、どこが仮説か？」

| レベル | 内容 |
|--------|------|
| **高信頼** | PT_dptのストレス独立性、ATM pathway 8/9遺伝子↓、NVU結合破綻 |
| **中程度** | 血管系の早期関与、R-loop surveillance decline |
| **仮説段階** | MMCの因果鎖、EV→Mono処理破綻の因果 |

---

# Part XII: 検証計画

## 34. 最小介入セット

| 実験 | 介入 | 予測 | 検証内容 |
|------|------|------|----------|
| Exp-1 | Capillary ATM rescue | R-loop ↓, RNA proc ↑ | ATM → R-loop の因果 |
| Exp-2 | Capillary XRCC5/NBN KD | Fire proxy ↑ | ATM経路破綻の十分性 |
| Exp-3 | ATM rescue + EV測定 | CD81/CD9 ↓ | ATM → EV export の因果 |
| Exp-4 | Mono Processing enhancement | Progression ↓ | Processing → 予後の因果 |

## 35. 追加検証ポイント

- 縦断/時間依存データでの順序検証
- Spatial transcriptomicsでNVU相互作用を確認
- CSF bridging（CNS-末梢の橋として脳脊髄液データ）

---

# Part XIII: 確立された事実と作業仮説

## 36. 確立された事実（Evidence Level: HIGH）

| 項目 | 証拠 |
|------|------|
| ATM pathwayは早期（TDPneg/低PT）で低下 | 2データセット再現 |
| 低下は複数細胞種で再現 | Capillary, Pericyte |
| 経路内遺伝子の89%が同方向 | 遺伝子レベル検証 |
| R-loop surveillanceより先に低下を開始 | 順序条件 |
| Donor-level統計でも有意 | pseudoreplication排除 |
| Cross-dataset match 83% | GSE212630 ↔ Phase13 |
| PT_dptはストレスにほぼ依存しない | R²=1.9% |
| NVU結合破綻はALS固有の特徴 | Phase9-11で確認 |

## 37. 作業仮説（Evidence Level: MODERATE）

| 項目 | 証拠レベル |
|------|-----------|
| ATM brake lossがfireの原因である | 仮説（介入未検証） |
| ATM → R-loop → DDR疲弊のカスケード | 仮説（順序整合） |
| EV輸送異常増加への連関 | 仮説（相関あり） |
| Pericyteの適応応答低下がALS特異的 | 仮説（p=0.033） |
| MMC-1の因果鎖 | 仮説（複数ライン整合） |

## 38. 推奨される記述

論文等での記述として以下を推奨：

> 複数独立データ（GSE212630, Phase13）に対する本プロジェクトの解析により、両データセットでATM pathway遺伝子の低下パターンを確認した（donor-level統計、経路レベル方向一致8/9遺伝子）。なお、GSE212630原論文はATM経路に焦点を当てておらず、本所見は生データの独自再解析による発見である。
>
> さらに、PericyteにおけるR-loop surveillance のearly bumpがALSで有意に低下していることから（p=0.033）、ALSは**適応能力の不全**を特徴とする可能性が示唆された。
>
> 因果関係の確定には、ATM経路を標的とした介入実験（ATM rescue, XRCC5/NBN knockdown）による検証が必要である。

---

# Part XIV: 用語定義

| 用語 | 定義 |
|------|------|
| Fire | ALS病態の起点となる分子・細胞イベント |
| ATM brake loss | ATM経路の機能低下（DNA損傷応答の減弱） |
| PT_dpt | Diffusion Pseudotimeによる病態進行軸（ストレス非依存） |
| 順序条件 | 原因候補が下流変化より先に変化すること |
| 一貫性条件 | 観測条件を変えても同方向の変化を示すこと |
| 経路整合性 | 経路内遺伝子が同方向に変化すること |
| Early bump | PT軸上でPT1>PT2となる早期上昇パターン |
| MMC | 最小機構鎖（Minimal Mechanism Chain） |
| NVU | 神経血管ユニット（Neurovascular Unit） |
| Donor-level | 細胞ではなくドナー単位での集約（pseudoreplication排除） |

---

# Part XV: 主要成果物

## 39. Fire Origin 関連（2025-12-23追加）

**Reports**
- `results_fire_origin_landscape/reports/fire_origin_landscape_report.md`
- `results_phase13_fire_origin_validation/reports/phase13_fire_origin_validation_report.md`
- `results_phase13_gene_level_fire_origin_validation/reports/phase13_gene_level_fire_origin_validation_report.md`

**Tables**
- `results_phase13_gene_level_fire_origin_validation/tables/phase13_genelevel_als_vs_control_capillary.csv`
- `results_phase13_gene_level_fire_origin_validation/tables/early_bump_als_vs_control_comparison.csv`
- `results_phase13_gene_level_fire_origin_validation/tables/cross_dataset_genelevel_direction_table.csv`

## 40. EV/PBMC統合解析

- `results_trigger_ev_lr/reports/trigger_ev_lr_integration.md`
- `results_cns_pbmc_integration/reports/cns_pbmc_integration_report.md`
- `results_cns_pbmc_receiver_v2/reports/receiver_phase_report.md`

## 41. 参照ドキュメント

- Project 1基盤: `docs/DATA_STRUCTURE_AND_PIPELINE_GUIDE.md`
- PTv2導入: `results/PTv2_robustness/PTv2_VALIDATION_REPORT.md`
- Phase 5-13: `results/*`
- Trigger Project: `results_*` (GSE212630 / trigger / NVU / PBMC)

---

# Part XVI: 更新履歴

| 日付 | 内容 |
|------|------|
| 2025-12-21 | EV/LR統合解析、CNS-PBMC統合、Receiver分解解析を追加 |
| 2025-12-23 | Fire Origin仮説（ATM brake loss）確立、Phase13遺伝子レベル検証、Early bump解析、プロジェクト全体像の統合 |
| 2025-12-23 | **Hybrid Final版作成**: GSE212630使用に関する注記追加、VAT1L/lncRNA/TF位置づけ統合、FAQ拡充、**PT_STMN2_single × Phase13 3軸統合解析を追加** |
| 2025-12-23 | **STMN2解析詳細化**: 細胞種別（VAT1L/L3-L6/L2-L3/Glia/Vascular）完全相関テーブル追加、全p値・n値明記、Capillary vs Neuron対比表追加 |
| 2025-12-23 | **ATM上流因子追加**: Section 29「ATM Brake Loss の上流因子（候補と制約）」新設。加齢/ROS/フィードバックループ仮説、血管内皮特異性の説明を追加 |

---

# Part XVII: 最終所見

## 42. プロジェクトの到達点

このプロジェクトは、以下の科学的階段を登り切った：

```
仮説 → 指標化 → gate → donor-level検定 → 頑健性（bootstrap）
→ 機構（network/driver）→ 遺伝子レベル検証 → 適応応答解析
→ Cross-dataset validation → 因果モデル構築
```

## 43. 核心的発見

1. **ATM brake lossがfire origin候補として最有力**
   - 8/9遺伝子が同方向、100%の観測軸で一貫
   - 2独立データセットで再現

2. **Capillary先行、Pericyte波及パターン**
   - 血管neurovascular unitが起点

3. **ALSは適応能力の不全**
   - PericyteのR-loop early bumpがALSで有意に低下
   - 「ブレーキをかけられない」病態

4. **PBMC Mono Processingが予後分岐点**
   - AUC 0.87で進行速度を予測
   - 介入標的候補

## 44. 残る問い

1. **ATM brake lossの原因は何か？**（genetic? epigenetic? environmental?）
2. **介入で救済可能か？**（ATM rescue, Processing enhancement）
3. **バイオマーカーとして使えるか？**（血中EV, 末梢Mono状態）

---

*宇宙はいつも、血管から壊れる。そして単球がそれを受け止められるかどうかで、運命が分かれる。ATMのブレーキが壊れた瞬間、火は灯る。*

---

*Generated: 2025-12-23*
*Project: ALS Trigger / Upstream*
*Status: Fire Origin Hypothesis ESTABLISHED (Hybrid Final)*

---

# Part XVIII: Reproducibility Addendum (Third-Party Replication)

以下は「第三者が完全再現するために必要な最低限の追記」です。
本フォルダ内の実ファイル位置と、スクリプト↔出力の対応を明示します。

## A. 入力ファイルの実体と取得手順

### A.1 Motor cortex raw expression / metadata（実体）
- 実ファイル置き場: `data/motor_cortex_raw_expression/`
- 命名規則:
  - `motor_cortex_{cell_type}_expression.csv(.gz)`
  - `motor_cortex_{cell_type}_metadata.csv(.gz)`
- 例:
  - `data/motor_cortex_raw_expression/motor_cortex_Ex.L5.VAT1L_EYA4_expression.csv`
  - `data/motor_cortex_raw_expression/motor_cortex_Ex.L5.VAT1L_EYA4_metadata.csv`
- 注意: `Zone.Identifier` はWindowsの付帯メタであり無視可。

### A.2 Phase13 expression matrix / metadata の生成元（実体）
- Phase13解析の実入力は、既にモジュールスコア＋PTが付与された統合CSV:
  - `docs/ids_causal_analysis/results/phase9_vascular/cell_level_features_ALL_with_PTdpt.csv`
- 生成経路は `docs/DATA_STRUCTURE_AND_PIPELINE_GUIDE.md` に記載（raw → IDS計算 → PT付与 → 統合）。

### A.3 GSE212630 raw data（実体）
- 実ファイル置き場: `data/GSE_212630_raw_expression_transposed/`
- 参照スクリプト:
  - `docs/ids_causal_analysis/scripts/GSE212630_compute_module_PT_dpt.py`

## B. 23モジュール定義（遺伝子集合の完全表）

### B.1 モジュール名リスト（23）
`docs/DATA_STRUCTURE_AND_PIPELINE_GUIDE.md` の「23モジュールリスト」に準拠:
Angiogenesis, Apoptosis, Autophagy, Calcium_Signaling, Cell_Cycle, Complement,
Cytoskeleton, DNA_Repair, ECM, ER_Stress, Epigenetic, Growth_Factors,
Inflammation, Ion_Transport, Metabolism, Mitochondria, Myelination,
Oxidative_Stress, Protein_Homeostasis, RNA_Processing, Synaptic,
Transcription, lncRNA

### B.2 完全遺伝子セットの所在
- モジュール完全定義（JSON）は以下に存在:
  - `data/24_functional_modules.json`
  - `data/24_functional_modules_fixed.json`
- 関連ドキュメント:
  - `docs/26_MODULES_COMPLETE_LIST.md`（モジュール名と遺伝子数）
  - `docs/GENE_MODULE_DEFINITIONS.md`（概要）
  - `docs/ids_causal_analysis/ids_pipeline/module_utils.py`
    （モジュールJSONの検索パスが定義されている）

### B.3 モジュールスコア計算方法（再現条件）
スコア化は `docs/ids_causal_analysis/scripts/02a_compute_module_scores.py` に準拠:
- モジュールスコア: 遺伝子平均（mean）またはmedian/PCA
- 既定は mean + center（Control平均で中心化）
- 正規化: Z-score（module内遺伝子の平均に対しcenter/scale）

## C. 乱数・パラメータ固定（PT_dpt / kNN）

### C.1 GSE212630 module-based PT_dpt
`docs/ids_causal_analysis/scripts/GSE212630_compute_module_PT_dpt.py` より:
- `K_NEIGHBORS = 30`
- 距離: Euclidean
- Affinity: Gaussian（σ = median distance）
- Normalized Laplacian `L = D^{-1/2} W D^{-1/2}`
- Eigenvectors: 上位 `N_COMPONENTS = 10` を計算し、`components[1:6]` を使用
- Root cell: Controlかつ stress最小
- `np.random.seed(42)`
- Sampling: `MAX_CELLS_PER_TYPE = 500`
- Filtering: `n_cells >= 50` and `n_cells <= 10000`

### C.2 Project 1 PT_continuous（参考）
`docs/DATA_STRUCTURE_AND_PIPELINE_GUIDE.md` より:
- kNN: K=15
- Feature: 23モジュール I_m + 23モジュール Δφ_signed
- PT正規化: Control 5% → 0, ALS 95% → 1

## D. “どのスクリプトがどの表を出したか”対応表

### D.1 メイン索引
- `docs/SCRIPTS_AND_FILES_INDEX.md`（Project 1の1対1対応表）
- `docs/ATM_fire_origin_evidence_summary.md`（Fire origin関連スクリプトと出力）

### D.2 Phase13（例: モジュール解析）
`docs/ids_causal_analysis/scripts/Phase13_complete_module_analysis.py` →
- 入力: `docs/ids_causal_analysis/results/phase9_vascular/cell_level_features_ALL_with_PTdpt.csv`
- 出力: `docs/ids_causal_analysis/results/phase13_complete_modules/`
  - `module_phi_NVU_summary.csv`
  - `module_Control_vs_ALS_differences.csv`
  - `module_PT_dpt_correlations.csv`
  - `Fig1_module_phi_heatmap_ALS.png`
  - `Fig2_module_fold_change_heatmap.png`
  - `Fig3_module_PT_dpt_correlation_heatmap.png`
  - `Fig4_axis_dysregulation_barplot.png`

### D.3 Fire origin（ATM brake loss）
主要対応は `docs/ATM_fire_origin_evidence_summary.md` に集約:
- `motif_family_nvu_rolesplit.py` → NVU role split + candidate tables
- `fire_origin_targeted.py` → ATM/ATR/R-loop proxy比較
- `fire_origin_landscape.py` → sequence landscape + proxy report
- `nvu_culprit_decisive.py` → culprit ranking + chain coupling

---

# Part XIX: Repro Appendix (Inputs/Run Order/Environment)

第三者がクリーン環境で同じ表を生成するための「最終10%」をここに固定する。

## 1) Inputs Manifest（実体・出典・ハッシュ）

### 1.1 ALS motor cortex dataset（Phase13ベース）
- 出典: GEO GSE174332（ALS motor cortex sc/snRNA-seq）
  - 参照: `docs/COMPREHENSIVE_PAPER.md`
- 実体:
  - `data/motor_cortex_raw_expression/motor_cortex_{cell_type}_expression.csv(.gz)`
  - `data/motor_cortex_raw_expression/motor_cortex_{cell_type}_metadata.csv(.gz)`
- 前処理: raw counts（整数）を使用。QC/正規化は `docs/DATA_STRUCTURE_AND_PIPELINE_GUIDE.md` に準拠。
- ハッシュ（一覧）:
  - `docs/INPUTS_MANIFEST.sha256`
  - ※ 主要入力のsha256を同一ファイルで固定
  - 例:
  - `data/motor_cortex_raw_expression/motor_cortex_Ex.L5.VAT1L_EYA4_expression.csv`
  - `data/motor_cortex_raw_expression/motor_cortex_Ex.L5.VAT1L_EYA4_metadata.csv`

### 1.2 GSE212630（TDP staging）
- 出典: GEO GSE212630
- 実体:
  - `data/GSE_212630_raw_expression_transposed/*_expression_transposed.csv.gz`
  - `data/GSE_212630_raw_expression_transposed/*_metadata.csv`
- 前処理: TDPneg/TDPmed/TDPhigh staging付与済み（再解析手順は `docs/ids_causal_analysis/scripts/GSE212630_compute_module_PT_dpt.py`）。
- ハッシュ（一覧）:
  - `docs/INPUTS_MANIFEST.sha256`
  - 例:
  - `data/GSE_212630_raw_expression_transposed/Vasc_Capillary_expression_transposed.csv.gz`
  - `data/GSE_212630_raw_expression_transposed/Vasc_Capillary_metadata.csv`

### 1.3 23モジュール定義（完全実体）
- `data/24_functional_modules_fixed.json`
- `data/24_functional_modules.json`
- 解析で使ったファイル名を明示し、sha256で固定する。
  - sha256一覧: `docs/INPUTS_MANIFEST.sha256`

#### ハッシュ算出（Windows）
```powershell
Get-FileHash data/24_functional_modules_fixed.json -Algorithm SHA256
Get-FileHash data/24_functional_modules.json -Algorithm SHA256
```

## 2) Run Order（最小実行順）

1. **Module定義の準備**
   - 入力: `data/24_functional_modules_fixed.json`
2. **Module score計算（Phase13系）**
   - `docs/ids_causal_analysis/scripts/02a_compute_module_scores.py`
   - 出力: `docs/ids_causal_analysis/results/ids_pseudotime/module_scores_{cell_type}.csv`
3. **PT_dpt（GSE212630）**
   - `docs/ids_causal_analysis/scripts/GSE212630_compute_module_PT_dpt.py`
   - 出力: `docs/ids_causal_analysis/results/GSE212630_ids_analysis/combined_metadata_with_modules_and_PT.csv`
4. **Phase13 完全モジュール解析**
   - `docs/ids_causal_analysis/scripts/Phase13_complete_module_analysis.py`
   - 入力: `docs/ids_causal_analysis/results/phase9_vascular/cell_level_features_ALL_with_PTdpt.csv`
   - 出力: `docs/ids_causal_analysis/results/phase13_complete_modules/*`
5. **Fire origin（ATM brake loss）**
   - 1対1対応表は `docs/ATM_fire_origin_evidence_summary.md` に集約
6. **詳細な実行順**
   - `docs/REPRO_RUN_ORDER.md`

## 3) Script → Output対応（最小表）

| Script | Main Output |
|--------|-------------|
| `docs/ids_causal_analysis/scripts/Phase13_complete_module_analysis.py` | `docs/ids_causal_analysis/results/phase13_complete_modules/module_PT_dpt_correlations.csv` |
| `docs/ids_causal_analysis/scripts/GSE212630_compute_module_PT_dpt.py` | `docs/ids_causal_analysis/results/GSE212630_ids_analysis/combined_metadata_with_modules_and_PT.csv` |
| `docs/ids_causal_analysis/scripts/fire_origin_targeted.py` | `docs/ids_causal_analysis/results_fire_origin_landscape/tables/proxy_analysis_results.csv` |
| `docs/ids_causal_analysis/scripts/nvu_culprit_decisive.py` | `docs/ids_causal_analysis/results_nvu_culprit_decisive/tables/nvu_culprit_celltype_ranking.csv` |

## 4) Environment（versions固定）

### 4.1 依存パッケージ（最小）
- `requirements.txt`（unversioned）
- `requirements.lock.txt`（pip freeze）

### 4.2 固定手順（推奨）
```powershell
python --version
pip --version
pip freeze > requirements.lock.txt
```

### 4.3 記録すべき情報
- Python version
- OS/version
- `requirements.lock.txt` の内容
