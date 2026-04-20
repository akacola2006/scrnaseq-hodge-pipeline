# 26機能モジュール完全リスト

## 概要
ALS研究で使用される26の機能モジュールとその遺伝子数です。
このモジュール定義は`26_functional_modules_split.json`に含まれています。

## 26モジュールと遺伝子数

| # | モジュール名 | 遺伝子数 | 主な機能 |
|---|-------------|----------|---------|
| 1 | **Synaptic** | 833 | シナプス伝達、神経伝達物質放出 |
| 2 | **Oxidative_Stress** | 1,649 | 酸化ストレス応答、ROS制御 |
| 3 | **Protein_Homeostasis** | 1,388 | タンパク質折り畳み、品質管理 |
| 4 | **RNA_Processing** | 1,834 | RNAスプライシング、RNA結合 |
| 5 | **Apoptosis** | 1,387 | プログラム細胞死、アポトーシス |
| 6 | **Autophagy** | 971 | オートファジー、細胞内分解 |
| 7 | **Inflammation** | 1,721 | 炎症応答、免疫反応 |
| 8 | **Calcium_Signaling** | 965 | カルシウムシグナリング |
| 9 | **Ion_Transport** | 1,033 | イオンチャネル、イオン輸送 |
| 10 | **Mitochondria** | 1,035 | ミトコンドリア機能、エネルギー代謝 |
| 11 | **ER_Stress** | 1,345 | 小胞体ストレス、UPR |
| 12 | **ALS_Genes** | 10 | ALS関連遺伝子（SOD1, FUS, TDP-43等） |
| 13 | **DNA_Repair** | 2,443 | DNA損傷修復、ゲノム安定性 |
| 14 | **Cell_Cycle** | 1,358 | 細胞周期制御 |
| 15 | **Cytoskeleton** | 2,638 | 細胞骨格、軸索輸送 |
| 16 | **Metabolism** | 3,797 | 代謝経路全般 |
| 17 | **ECM** | 1,051 | 細胞外マトリックス |
| 18 | **Angiogenesis** | 697 | 血管新生 |
| 19 | **Growth_Factors** | 742 | 成長因子シグナリング |
| 20 | **Transcription** | 4,218 | 転写制御 |
| 21 | **Epigenetic** | 842 | エピジェネティック制御 |
| 22 | **Complement** | 47 | 補体系 |
| 23 | **Myelination** | 102 | ミエリン形成 |
| 24 | **lncRNA** | 14,061 | 長鎖非コードRNA |
| **25** | **(未定義)** | - | - |
| **26** | **(未定義)** | - | - |

## 注記

- **合計24モジュール**: 現在のファイルには24個の機能モジュールが定義されています
- **lncRNA**: 最大のカテゴリで14,061遺伝子を含む
- **ALS_Genes**: 最小のカテゴリで10個のコア ALS遺伝子のみ
- **Complement**: 次に小さく47遺伝子

## ALS研究における重要モジュール

### 特に重要な5モジュール
1. **RNA_Processing** - TDP-43, FUS関連
2. **Protein_Homeostasis** - 異常タンパク質蓄積
3. **Mitochondria** - エネルギー代謝異常
4. **Cytoskeleton** - 軸索変性
5. **ALS_Genes** - 直接的な原因遺伝子

### データソース
- GO database: Enrichr GO 2023
- lncRNA source: ENSG_allgene_allcharacter.csv
- Creation date: 2025-07-15

## ファイル
- `26_functional_modules_split.json` - 完全な遺伝子リスト
- `24_functional_modules.json` - 24モジュール版
- `24_functional_modules_fixed.json` - 修正版