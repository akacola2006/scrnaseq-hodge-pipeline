# 遺伝子モジュール定義

## 概要
このプロジェクトで使用される遺伝子機能モジュールの定義です。
主にGO（Gene Ontology）とパスウェイベースで分類されています。

## モジュール定義ファイル

### 1. メインモジュール定義
- `src/ids_module_stress_analyzer.py`: GO/パスウェイベースの20モジュール
- `src/ids_go_module_stress_analyzer.py`: 拡張版モジュール定義
- `als_causal_inference_analysis/calculate_proper_cell_module_cie.py`: 標準24モジュール

## 定義されているモジュール（20カテゴリ）

### RNA代謝・スプライシング
1. **RNA_processing** (GO:0006396)
   - 遺伝子: TARDBP, FUS, HNRNPA1, HNRNPA2B1, TIA1, EWSR1, TAF15

2. **RNA_splicing** (GO:0008380)
   - 遺伝子: TARDBP, FUS, HNRNPA1, HNRNPA2B1, SFPQ, MATR3

### タンパク質品質管理
3. **Protein_folding** (GO:0006457)
   - 遺伝子: HSPA1A, HSPA8, HSPB1, DNAJB1, BAG3, HSPH1

4. **Ubiquitin_proteasome** (GO:0030433)
   - 遺伝子: UBQLN2, VCP, SQSTM1, OPTN, TBK1, UBE2A

5. **Autophagy** (GO:0006914)
   - 遺伝子: SQSTM1, OPTN, TBK1, VCP, C9orf72, UBQLN2

### ミトコンドリア機能
6. **Oxidative_phosphorylation** (GO:0006119)
   - 遺伝子: COX1-3, ND1-6

7. **Mitochondrial_transport** (GO:0006839)
   - 遺伝子: CHCHD10, CHCHD2, PINK1, PRKN, MFN2

### 軸索輸送
8. **Anterograde_axonal_transport** (GO:0008089)
   - 遺伝子: KIF5A, KIF5B, KIF5C, DCTN1, DYNC1H1

9. **Retrograde_axonal_transport** (GO:0008090)
   - 遺伝子: DCTN1, DYNC1H1, DYNLL1, DYNLL2

10. **Neuron_projection** (GO:0031175)
    - 遺伝子: NEFL, NEFM, NEFH, PRPH, INA

### シナプス機能
11. **Synaptic_transmission** (GO:0007268)
    - 遺伝子: SYT1, SYT2, SYN1, SYN2, SNAP25, STX1A

12. **Synaptic_signaling** (GO:0099536)
    - 遺伝子: GRIN1, GRIN2A, GRIA1, GRIA2, GABRA1

### 炎症・免疫応答
13. **Inflammatory_response** (GO:0006954)
    - 遺伝子: TNF, IL1B, IL6, NFKB1, RELA, TLR4

14. **Immune_system_process** (GO:0002376)
    - 遺伝子: TREM2, TYROBP, CD33, MS4A4A, MS4A6A

### 細胞死
15. **Cell_death** (GO:0008219)
    - 遺伝子: CASP3, CASP9, BCL2, BAX, BID, APAF1

16. **Apoptotic_signaling** (GO:0097190)
    - 遺伝子: TP53, PUMA, NOXA, BIM, BAD

### その他重要機能
17. **DNA_damage_response** (GO:0006974)
    - 遺伝子: ATM, ATR, CHEK1, CHEK2, BRCA1, BRCA2

18. **Cytoskeleton_organization** (GO:0007010)
    - 遺伝子: ACTB, TUBB, PFN1, CFL1, TUBA1A

19. **Transcription_regulation** (GO:0006355)
    - 遺伝子: HDAC1-3, HDAC6, SIRT1, EP300

20. **Ion_Homeostasis** (カスタム)
    - イオンチャネルとカルシウム制御関連遺伝子

## 標準24モジュール（calculate_proper_cell_module_cie.pyより）

上記20モジュールに加えて：

21. **Stress_Response** / **ROS_Stress**
22. **Protein_Aggregation** / **Protein_Homeostasis**
23. **Axon_Degeneration**
24. **Calcium_Signaling**

## 使用方法

```python
# モジュール定義の読み込み例
from src.ids_module_stress_analyzer import ModuleStressAnalyzer

analyzer = ModuleStressAnalyzer()
analyzer.load_go_annotations()

# モジュール遺伝子の取得
rna_genes = analyzer.module_genes['GO:0006396_RNA_processing']
```

## 注記
- これらのモジュールはALS研究に特化して選定されています
- GO IDは標準的なGene Ontologyデータベースに準拠
- 一部カスタムモジュールも含まれています（Ion_Homeostasis等）