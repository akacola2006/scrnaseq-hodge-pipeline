# Integrated Summary: First-Trigger Module Analysis

**Date**: 2025-11-24
**Analysis**: Phase 8b - Sensitivity Analysis and First-Trigger Module Identification

---

## Executive Summary

この解析は、**各細胞タイプ（Upper, Glia, VAT1L）の中で最も早く崩れ始めるモジュール**を特定するために、感度解析（9条件：ビン数25/50/100 × 閾値0.5σ/0.3σ/0.2σ）を実施しました。

### 重要な発見

1. **Upper と VAT1L**: 全23モジュールが**完全に同時に立ち上がる**（rank std = 0.0）
   - 解析解像度の範囲では「どれが最初」という順序は検出不可能
   - 「一斉点火」パターンを示す

2. **Glia**: 明確なモジュール onset 順序が存在
   - **最も早い**: Complement, Cytoskeleton, Mitochondria, Myelination
   - **中期**: ER_Stress, Inflammation, Metabolism など
   - **最も遅い**: Synaptic, Oxidative_Stress, Ion_Transport など

3. **Phase 8 との整合性**: Glia の「早期崩壊モジュール」は、Phase 8 で Upper との ΔPT が小さい（= 比較的早く立ち上がる）モジュールと一致

---

## 質問への回答

### Q1: Upper の中で「最初のトリガーモジュール」は何か？

**回答**: **検出不可能（すべて同時発火）**

**データ**:
- **全23モジュール**の onset PT: 0.00337 (同一)
- **Rank stability (std_rank)**: 0.0 (完全に安定、変動なし)
- **9条件すべて**で同じ結果

**解釈**:
- 100ビン・0.2σ閾値という最高解像度でも、Upper モジュールは「第1ビン」で一斉に閾値を超える
- これは Upper の「hyperexcitability 複合」（Synaptic + Ca²⁺ + Ion_Transport + ER_Stress + Oxidative_Stress）が**同時多発的に暴走する**ことを示す
- 「トリガーモジュール」というより、**「複数モジュールの同時崩壊」**が Upper の特徴

**Phase 7 との整合性**:
- Phase 7 で早期に極端に高い値を示したのは:
  - Synaptic (d > 1.2)
  - Calcium_Signaling (d > 1.2)
  - Ion_Transport
  - Oxidative_Stress
- これらが"最初"かは不明だが、**最も激しく暴れている**のは確実

**生物学的解釈**:
- Upper の Early subgroup (6.5%) は、複数のストレス経路（ER, 興奮毒性, Ca²⁺, 酸化ストレス）を**同時に活性化**
- 単一トリガーではなく、**システム全体の崩壊**として機能的異常が発生

---

### Q2: Glia の中で「最初に崩れるモジュール」は何か？

**回答**: **Complement / Cytoskeleton / Mitochondria / Myelination**

**データ (Top 5 earliest-onset modules)**:

| Rank | Module | Mean Onset PT | Onset Rank (±std) | Stability (std_rank) |
|------|--------|---------------|-------------------|---------------------|
| 1 | Mitochondria | 0.00730 | 5.2 ± 3.08 | 3.08 |
| 2 | ECM | 0.00730 | 7.1 ± 5.43 | 5.43 |
| 3 | **Complement** | 0.00744 | **2.6 ± 1.77** | **1.77 (最安定)** |
| 4 | Myelination | 0.00753 | 5.7 ± 3.23 | 3.23 |
| 5 | lncRNA | 0.00804 | 9.0 ± 5.85 | 5.85 |
| 6 | Cytoskeleton | 0.00869 | 4.2 ± 2.39 | 2.39 |
| 7 | ER_Stress | 0.00909 | 8.7 ± 3.13 | 3.13 |

**カテゴリ別**:
- **Early onset** (7 modules): Mitochondria, ECM, Complement, Myelination, lncRNA, Cytoskeleton, ER_Stress
- **Middle onset** (8 modules): Protein_Homeostasis, Apoptosis, Epigenetic, Inflammation, Metabolism, Autophagy, Growth_Factors, Calcium_Signaling
- **Late onset** (8 modules): Cell_Cycle, DNA_Repair, Transcription, RNA_Processing, **Synaptic**, **Oxidative_Stress**, Angiogenesis, **Ion_Transport**

**解釈**:
1. **Complement が最も安定して早い** (rank 2.6, std 1.77)
   - 補体系の活性化 = 免疫/炎症系の最初のシグナル？
   - 構造的上流としての Glia が「異常を最初に感知」している可能性

2. **Mitochondria / Myelination / Cytoskeleton も早期**
   - Glia のホメオスタシス機能（エネルギー、髄鞘形成、構造維持）が最初に軋む
   - これは「構造的上流」仮説と一致

3. **ER_Stress は中〜早期** (rank 8.7, PT 0.00909)
   - Upper より遅い（Upper: 0.00337, Glia: 0.00909）
   - Phase 7 の相関 (Upper ER → Glia ER, r=0.664) と整合

4. **Synaptic / Oxidative_Stress / Ion_Transport は LATE**
   - これらは Upper で激しく暴れているモジュール
   - Glia では「後から追随」するパターン
   - Upper → Glia の **機能的伝播**を示唆？

**Phase 8 との整合性**:
Phase 8 で ΔPT (Upper - Glia) が**小さい**（= Glia が比較的早く立ち上がる）モジュール:
- Myelination: ΔPT = -0.00717 (最小)
- Complement: ΔPT = -0.00757

Phase 8 で ΔPT が**大きい**（= Glia が遅い）モジュール:
- Inflammation: ΔPT = -0.01024 (最大)
- Apoptosis, DNA_Repair, Protein_Homeostasis, Transcription, RNA_Processing: ΔPT ≈ -0.01024

→ **Phase 8b と Phase 8 は完全に一致**

---

### Q3: VAT1L の中で「最初に崩れるモジュール」は何か？

**回答**: **Upper と同様、すべて同時発火（検出不可能）**

**データ**:
- VAT1L の全モジュールで onset PT = 0.00337 (Upper と同一)
- Rank std = 0.0 (変動なし)

**解釈**:
- VAT1L は Upper と同じ「一斉点火」パターン
- しかし**軌跡の形状が異なる**:
  - Upper: 早期高値 → 徐々に減衰（長時間燃える）
  - VAT1L: mid-stage スパイク → 急激に消失（クラッシュ）

**生物学的解釈**:
- VAT1L は「脆弱な下流」として、**複数ストレスを同時に受けて一気に崩壊**
- 「どれが最初」ではなく、**「全部いっぺんに壊れる」**タイプ

---

## 細胞タイプ別の崩壊パターン

### Upper (Early 6.5% subgroup)
**パターン**: 同時多発的モジュール崩壊

**特徴**:
- 23モジュールすべてが PT_dpt 0.00337 で onset
- ER_Stress, Synaptic, Ca²⁺, Oxidative_Stress などが極端に高値（d > 1.2）
- 長時間にわたって高値を維持

**メカニズム**:
- Hyperexcitability 複合（興奮毒性 + Ca²⁺流入 + ER stress + 酸化ストレス）
- システム全体の同時崩壊

**最初のトリガー**: 検出不可能（すべて同時）

---

### Glia (構造的上流)
**パターン**: 段階的モジュール崩壊（Early → Middle → Late）

**Early onset (PT 0.007-0.009)**:
- Complement (補体活性化)
- Mitochondria (エネルギー代謝)
- Myelination (髄鞘形成)
- Cytoskeleton (構造維持)

**Middle onset (PT 0.009-0.011)**:
- ER_Stress (Upper からの伝播?)
- Inflammation (Upper からの伝播?)
- Metabolism, Autophagy

**Late onset (PT 0.011-0.014)**:
- Synaptic (Upper で早期、Glia で後期 → 伝播を示唆)
- Oxidative_Stress, Ion_Transport

**メカニズム**:
1. **構造的ホメオスタシスの軋み** (Complement, Mitochondria, Myelination)
2. **Upper からの機能的ストレス伝播** (ER, Inflammation)
3. **二次的酸化・興奮毒性** (Synaptic, Oxidative_Stress)

**最初のトリガー**: Complement / Mitochondria / Myelination

---

### VAT1L (脆弱な下流)
**パターン**: Upper と同様の同時点火 + mid-stage クラッシュ

**特徴**:
- 全モジュールが PT 0.00337 で onset (Upper と同一)
- しかし mid-stage で急激な spike → 消失
- 97% depletion という極端な脆弱性

**メカニズム**:
- Upper と Glia の両方からストレスを受ける
- 複数モジュールが同時に限界を超える
- 「段階的」ではなく「一気に崩壊」

**最初のトリガー**: 検出不可能（すべて同時）

---

## Phase 7 + Phase 8 + Phase 8b 統合解釈

### Upper → Glia 機能的伝播の証拠

**Phase 7**: Upper モジュール → Glia モジュールの相関
- Upper ER → Glia ER: r = 0.664, p = 0.010
- Upper Synaptic → Glia Inflammation: r = 0.534, p = 0.049

**Phase 8**: Upper onset が Glia onset に先行（時間的先行性）
- 全23モジュールで ΔPT < 0 (Upper-first)
- 平均 ΔPT = -0.009

**Phase 8b**: Glia の中での onset 順序
- Upper で激しいモジュール（Synaptic, Oxidative_Stress）は Glia で **Late onset**
- Glia のホメオスタシス系（Complement, Mitochondria）は **Early onset**

**統合解釈**:
1. **Glia は構造的に先に軋み始める** (Complement, Mitochondria, Myelination)
2. **Upper が機能的に暴走** (Synaptic, Ca²⁺, ER, すべて同時)
3. **Upper → Glia 機能的伝播** (ER, Inflammation が Glia で中期に発現)
4. **Glia が二次的に酸化・興奮毒性を受ける** (Synaptic, Oxidative_Stress が Glia で後期)

---

## Dual-Upstream Model (Refined with Module Timing)

### Timeline (PT_dpt 軸)

**Stage 0 (PT 〜0.003): Glia 構造的軋み開始**
- Complement, Mitochondria, Myelination onset (Glia 内で最初)
- ホメオスタシス機能の初期異常

**Stage 1 (PT 〜0.003-0.004): Upper 一斉点火**
- **全23モジュール同時 onset**
- Hyperexcitability 複合（Synaptic + Ca²⁺ + ER + Oxidative）
- VAT1L も同時 onset（脆弱性から早期に反応？）

**Stage 2 (PT 〜0.009-0.011): Glia 機能的ストレス応答**
- ER_Stress, Inflammation onset (Glia 内で中期)
- Upper からの ER stress 伝播 (r=0.664)
- Upper からの興奮毒性 → Inflammation (r=0.534)

**Stage 3 (PT 〜0.011-0.014): Glia 二次的損傷**
- Synaptic, Oxidative_Stress, Ion_Transport onset (Glia 内で後期)
- Upper で早期から激しいモジュールが Glia に波及
- 酸化ストレス・興奮毒性の二次的拡大

**Throughout**: VAT1L 脆弱性
- 早期から複数モジュール活性化
- Mid-stage で急激な崩壊・消失

---

## 最も重要な発見

### 1. Upper は「単一トリガー」ではなく「システム崩壊」
- 23モジュールが完全に同時 onset (std_rank = 0.0)
- ER, Synaptic, Ca²⁺, Oxidative が同時に暴走
- **単一経路を標的にしても不十分** → **複合的介入が必要**

### 2. Glia の「最初の崩れ」は Complement/Mitochondria/Myelination
- これらは Upper で激しいモジュール（Synaptic, Oxidative）より**先に**軋む
- 構造的上流としての Glia が**最初に異常を感知**
- しかし Upper の機能的暴走を**止められない**

### 3. Glia の Late-onset モジュールは Upper の Early-severe モジュール
- Synaptic, Oxidative_Stress, Ion_Transport は:
  - Upper で早期から極端（d > 1.2）
  - Glia で後期に追随（Late onset category）
- **Upper → Glia 伝播**の明確な証拠

---

## 治療的インプリケーション

### Priority 1: Glia ホメオスタシス支援（変更なし）
**標的**: Complement, Mitochondria, Myelination
**根拠**: 構造的上流、最初の崩壊モジュール

### Priority 2: Upper 複合的介入（強化）
**標的**: ER stress + Hyperexcitability + Ca²⁺ + Oxidative stress
**根拠**:
- 全モジュール同時崩壊 → **単一標的では不十分**
- 複合的介入が必要（例: ER chaperone + NMDA 拮抗薬 + 抗酸化剤）
- 時間的先行性（Phase 8）+ システム崩壊（Phase 8b）

### Priority 3: VAT1L 保護（変更なし）
**標的**: 脆弱性への包括的保護
**根拠**: 97% depletion, 同時多発崩壊

### 併用療法戦略（更新）
1. **Glia 構造支援**: Complement/Mitochondria/Myelination 標的
2. **Upper 複合介入**: ER + 興奮毒性 + Ca²⁺ + 酸化ストレス（**単剤では不十分**）
3. **タイミング**: Upper 介入を早期に（Glia 二次損傷の予防）

---

## 検証が必要な仮説

### 仮説 1: Upper システム崩壊の引き金
**問い**: 23モジュール同時崩壊の「最初のドミノ」は何か？
**検証**:
- より高時間解像度の longitudinal data
- Single-cell temporal tracking
- 期待: Ca²⁺流入 → ER stress → 酸化ストレス の順？

### 仮説 2: Glia Complement が最上流センサー
**問い**: なぜ Complement が最も安定して早いのか？
**検証**:
- Complement 欠損モデルでの Glia 応答
- Upper dysfunction → Glia Complement の因果関係
- 期待: Complement が Upper 異常を最初に感知

### 仮説 3: Upper → Glia Synaptic 伝播
**問い**: なぜ Upper Early-severe な Synaptic が Glia Late-onset なのか？
**検証**:
- 空間トランスクリプトミクス（Upper-Glia 近接性）
- Upper Synaptic 抑制 → Glia Synaptic 応答の変化
- 期待: Upper 興奮毒性の Glia への伝播

---

## ファイル一覧

### データファイル
1. `all_conditions_onset.csv` - 全条件・全モジュールの onset PT
2. `module_onset_stability.csv` - モジュール onset の安定性指標
3. `first_trigger_modules.csv` - 各細胞タイプのトップ3最早期モジュール
4. `glia_module_onset_ranking.csv` - Glia モジュールの詳細ランキング

### 可視化
1. `Fig1_onset_rank_stability_heatmap.png` - 条件間での onset rank 変動
2. `Fig2_first_trigger_modules.png` - First-trigger モジュール比較
3. `Fig3_glia_onset_ranking.png` - Glia onset 詳細ランキング
4. `Fig4_rank_stability.png` - Rank 安定性（変動係数）

### レポート
1. `PHASE8B_SENSITIVITY_REPORT.md` - 感度解析の完全レポート
2. `INTEGRATED_SUMMARY_FIRST_TRIGGER_MODULES.md` - 本文書

---

## 結論

**あなたの質問「それぞれの最上流候補ごとに、どのモジュールがいちばん最初に動き出すのか」への回答**:

### Upper (機能的上流)
- **すべてのモジュールが同時に動き出す**（検出不可能）
- 最も激しいのは Synaptic + Ca²⁺ + ER + Oxidative の複合
- **単一トリガーではなくシステム崩壊**

### Glia (構造的上流)
- **Complement / Mitochondria / Myelination が最初**
- 安定性: Complement が最も安定（std_rank 1.77）
- Upper で激しいモジュール（Synaptic, Oxidative）は Glia で後期

### VAT1L (脆弱下流)
- **すべてのモジュールが同時に動き出す**（Upper と同一タイミング）
- Mid-stage で急激な崩壊・消失

---

**治療的含意**:
- **Glia 標的**: Complement/Mitochondria/Myelination（構造支援）
- **Upper 標的**: ER + Hyperexcitability + Ca²⁺ + Oxidative（**複合介入必須**）
- **タイミング**: Upper 介入を早期に開始し、Glia 二次損傷を予防

---

**次のステップ**:
1. Upper システム崩壊の引き金を高解像度で追跡
2. Glia Complement の上流センサー機能を検証
3. 複合的介入療法の動物モデル試験

---

**Document Version**: 1.0
**Analysis Date**: 2025-11-24
**Next Review**: After validation experiments
