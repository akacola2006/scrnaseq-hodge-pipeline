# scRNAseq Hodge Decomposition Pipeline

Companion code repository for:

> **Hodge decomposition of gene co-expression dynamics infers structural
> upstream organisation in sporadic ALS: from measurement instrument to
> integrated disease cascade**
> Kaneko S., Urushitani M. (Shiga University of Medical Science)

離散 Hodge 分解を用いて、**scRNAseq データから疾患の「上流細胞型」と
「上流遺伝子」を訓練なしで同定**する測定器 (Intrinsic Direction System; IDS)
の公式実装。論文で報告した全ての主要解析を再現可能な形で収録しています。

## Repository structure

```
scrnaseq_hodge_pipeline/
├── scripts/                        Part I: Core Hodge pipeline
│   ├── data_loader.py              データ読み込み (h5ad)
│   ├── residuals.py                GPU 加速回帰残差 (~ donor + nUMI + pct_mito)
│   ├── pca_engine.py               GPU ランダム PCA
│   ├── spd.py                      SPD 共分散 (Ledoit-Wolf + log-Euclidean)
│   ├── pseudotime.py               拡散マップ擬似時間
│   ├── hodge.py                    Hodge 分解 (セルタイプレベル)
│   ├── lane_a.py                   Lane A: 上流セルタイプ同定
│   ├── lane_b.py                   Lane B: 上流 PC 遺伝子抽出
│   ├── gene_hodge.py               遺伝子レベル K_N Hodge
│   ├── sparse_hodge.py             Dual-mode: K_N vs k-NN
│   ├── directional.py              Δ⁺/Δ⁻ 分離 Hodge
│   ├── multi_transition.py         全遷移統合 (stable-High)
│   ├── random_baseline.py          GF 帰無分布 (316 simulation)
│   ├── two_axis.py                 TRS × MSS 2軸モデル
│   ├── bootstrap.py                ドナーブートストラップ (B=100)
│   ├── enrichment.py               機能モジュール Fisher 検定
│   ├── whitening.py                ZCA whitening (摂動解析用)
│   ├── three_phi_residual.py       ★ 3φ residual framework (Section 4.1-4.2, 6.15)
│   ├── allgene_insertion.py        ★ All-gene 3φ insertion (Section 6.15, Appendix AS)
│   ├── nemf_screen.py              ★ NEMF genome-wide screen (Section 8.9)
│   ├── de_sals_only.py             ★ SALS-only DESeq2 + GSEA (Section 6.16)
│   ├── config.py / seed_utils.py / log_utils.py / file_checker.py
│
├── benchmark/                      External validation benchmarks
│   ├── README.md                   全ベンチマーク索引
│   ├── norman_perturb_seq/         Section 2.4 (IDS Top-1 = 66.4 %)
│   ├── shen_crispr/                Section 2.5 (mean 21.8 %ile, p = 0.002)
│   ├── timevault/                  Section 2.3 (TAS AUC = 0.9475)
│   ├── glioma/                     Section 2.6 (Null GF = 0.4248)
│   └── mendelian_randomization/    Section 4.5 / 8.7 (OR = 0.988, p = 0.159)
│
├── part3/                          Module-level causal analysis (Part III)
│   ├── README.md                   Section 9-11 完全実装（104 ファイル）
│   ├── ALS_TRIGGER_METHODOLOGY.md  legacy narrative 方法論
│   ├── ids_pipeline/               ★ IDSPipeline API + CLI (Python)
│   │   ├── pipeline.py             IDSPipeline class (orchestrator)
│   │   ├── run_pipeline.py         CLI entry point
│   │   ├── config.py, gene_utils.py, module_utils.py
│   │   └── templates/default_config.yaml
│   ├── scripts/                    ★ 52 analysis scripts
│   │   ├── ids_core.py             Core library (phi, PT binning, axis discovery)
│   │   ├── Phase8_*.py             Time-lag analysis
│   │   ├── Phase9*.py              Vascular (capillary + pericyte) analysis
│   │   ├── Phase10/11/12/13_*.py   NVU DAG, energy landscape
│   │   ├── 03a/b/c_*.py            NOTEARS + LiNGAM causal inference
│   │   ├── GSE212630_*.py (9)      Cross-dataset validation (Wang 2023)
│   │   ├── PTv2_*.py               PT_dpt robustness
│   │   ├── fire_origin_*.py        Fire origin identification
│   │   └── test_ids_core_*.py      Unit tests
│   ├── reports/ (24)               Phase reports + module definitions JSON
│   └── docs/ (18)                  Project-level narrative documentation
│
├── data/
│   ├── README.md
│   ├── h5ad/                       h5ad file 配置先
│   └── metadata/
│       ├── gene_annotation.csv              60,403 genes
│       ├── functional_modules.json          23 functional modules (Part III)
│       ├── als_causative_genes.json         ALS_Genes (Part II 参照, 10 genes)
│       └── functional_modules_full24.json   24 modules (backup: Part III + ALS_Genes)
│
├── README.md                       (this file)
├── CLAUDE.md                       Claude Code 用ドキュメント
├── project_config.yaml             パイプライン設定
├── requirements.txt                Python 依存パッケージ
├── pyproject.toml                  パッケージ定義
├── run_pipeline.py                 メイン実行スクリプト (Part I)
├── app.py                          Streamlit GUI
└── launcher.py                     起動ランチャー
```

## Paper Sections ↔ Scripts 対応表

| 論文 Section | 解析内容 | スクリプト |
|---|---|---|
| **Part I — The Instrument** | | |
| 2.1, 6.1-6.5 | Core Hodge pipeline (Stage 1-5) | `scripts/data_loader.py`, `residuals.py`, `pca_engine.py`, `spd.py`, `pseudotime.py`, `hodge.py`, `lane_a.py`, `lane_b.py`, `gene_hodge.py` |
| 2.2, 6.6, Supp Note 4 | Null GF baseline | `scripts/random_baseline.py` |
| 2.3, 6.11 | TimeVault TAS validation | `benchmark/timevault/run_tas.py` |
| 2.4, 6.8 | Norman Perturb-seq + whitening | `benchmark/norman_perturb_seq/run_ids_norman.py`, `scripts/whitening.py` |
| 2.5 | Shen CRISPR directional validation | `benchmark/shen_crispr/run_shen_ids.py` |
| 2.6 | Glioma cross-disease | `benchmark/glioma/*.py` |
| 3.1 | Oligo upstream + LODO/LOCO robustness | `scripts/lane_a.py`, `bootstrap.py`, `gene_hodge.py` |
| 3.3, 6.5 | Gradient / curl dichotomy | `scripts/sparse_hodge.py`, `directional.py` |
| 3.4, 6.10 | Two-axis TRS × MSS | `scripts/two_axis.py` |
| **Part I — Disease findings** | | |
| 4.1, 6.15 | Manifold preservation R² | `scripts/three_phi_residual.py` ★ |
| 4.2, 6.15 | Rewiring / collapse residual z | `scripts/three_phi_residual.py` ★ |
| 4.2, 6.16 | SALS-only DE + GSEA rank test | `scripts/de_sals_only.py` ★ |
| 4.2, 6.17 | Metascape rewiring enrichment | (external: [metascape.org](https://metascape.org); inputs from `three_phi_residual.py`) |
| 4.5, 6.12, 8.7 | Mendelian Randomization | `benchmark/mendelian_randomization/*.R` |
| 6.18 | LODO / LOCO robustness | `scripts/bootstrap.py` |
| **Part II** | | |
| 7.1-7.7 | Detection-rate gene-level φ landscape | `scripts/gene_hodge.py` + detection-rate variant (see `project_config.yaml`) |
| 8.9, Appendix AV | NEMF genome-wide screen | `scripts/nemf_screen.py` ★ |
| 8.9.5, Figure 6B | Cross-CT co-collapse / compensatory axes | `scripts/nemf_screen.py` ★ |
| **Part III** | | |
| 9.1 | 23 modules + PT_dpt + Phase 8 time-lag | `part3/scripts/Phase8_*.py`, `Phase9a_*.py`, `02a_compute_module_scores.py`, `PTv2_01_build_PT_dpt.py` |
| 9.1 | NOTEARS + LiNGAM causal inference | `part3/scripts/03a_prepare_notears_prior.py`, `03b_run_notears.py`, `03c_lingam_refinement.py` |
| 10.2 | ATM pathway 3-criteria evaluation | `part3/scripts/Phase9p_vascular_driver_state_identification.py`, `Phase9cd_vascular_timelag_correlation_analysis.py` |
| 10 | NVU DAG + MMC-1 cascade | `part3/scripts/Phase10_NVU_causal_DAG.py`, `Phase13_*.py`, `final_closure_mmc.py` |
| 10 | GSE212630 cross-dataset 83% concordance | `part3/scripts/GSE212630_*.py` (9 files), `analyze_GSE212630_with_ids_core.py` |
| ALL | IDSPipeline Python API | `part3/ids_pipeline/pipeline.py` + `run_pipeline.py` CLI |

★ = 本リリースで新規追加

## Quick start

### 1. 環境セットアップ

```bash
pip install -r requirements.txt

# Optional: GPU acceleration (for all-gene insertion)
pip install torch --index-url https://download.pytorch.org/whl/cu121

# Optional: DESeq2 (for scripts/de_sals_only.py)
pip install pydeseq2==0.5.0

# Optional: MR analysis (R)
# install.packages(c("TwoSampleMR", "coloc", "ieugwasr"))
```

### 2. データ配置

```
data/
  h5ad/                      # h5ad files (one per donor/sample)
    donor1.h5ad
    donor2.h5ad
    ...
  metadata/
    sample_info.csv          # donor_id, file, condition, sex
    gene_annotation.csv      # (pre-included)
    functional_modules.json  # (pre-included, 23 modules)
```

### 3. 設定編集

`project_config.yaml` を編集：
- `cell_types.include`: 解析する細胞型リスト
- `cell_types.pseudotime_set`: 擬似時間構築用（4–6 細胞型推奨）
- `conditions.column_name`, `control_label`, `disease_labels`
- `hardware.use_gpu`

### 4. 実行

```bash
# Part I core pipeline
python run_pipeline.py

# 3φ residual framework (Table 1 R², Rewiring/Collapse)
python -m scripts.three_phi_residual

# All-gene 3φ insertion (NEMF screen input)
python -m scripts.allgene_insertion --cell-type Oligo --n-base 3500

# NEMF screen (Section 8.9)
python -m scripts.nemf_screen \
    --zscore-matrix results/allgene_3phi/verification/zscore_matrix_wide.csv \
    --target-gene NEMF

# SALS-only DE + GSEA rank test
python -m scripts.de_sals_only

# External benchmarks (see benchmark/ and part3/ READMEs)
python benchmark/timevault/run_tas.py ...
python benchmark/norman_perturb_seq/run_ids_norman.py ...
```

## 再現可能性についての注意

本リポジトリは論文で報告した **全ての central findings** を再現可能な形で
実装しています：

- ✅ **Oligo φ = 0.900** (Section 3.1) — `scripts/lane_a.py`
- ✅ **LOCO Jaccard = 0.954, ρ = 0.981** (Section 3.1) — `scripts/bootstrap.py`
- ✅ **Null GF = 0.425** (Section 2.2) — `scripts/random_baseline.py`
- ✅ **TimeVault TAS AUC = 0.9475** (Section 2.3) — `benchmark/timevault/run_tas.py`
- ✅ **Norman Top-1 = 66.4 %** (Section 2.4) — `benchmark/norman_perturb_seq/run_ids_norman.py`
- ✅ **Shen CRISPR p = 0.002** (Section 2.5) — `benchmark/shen_crispr/run_shen_ids.py`
- ✅ **Glioma Null GF = 0.4248** (Section 2.6) — `benchmark/glioma/run_glioma_pseudotime_insertion.py`
- ✅ **NEMF z = −6.95 in L3_L5, Bonferroni p = 3.3×10⁻⁸** (Section 8.9) — `scripts/nemf_screen.py`
- ✅ **MR OR = 0.988, p = 0.159** (Section 4.5) — `benchmark/mendelian_randomization/02_run_mr.R`
- ✅ **ATM pathway p = 0.018, 8/9 concordant** (Section 10) — `part3/scripts/Phase9p_vascular_driver_state_identification.py`
- ✅ **GSE212630 83% cross-dataset concordance** (Section 10) — `part3/scripts/analyze_GSE212630_with_ids_core.py`

各数値は対応するスクリプトの実行により再現されます。必要に応じて
`benchmark/<name>/README.md` および `part3/README.md` の
「Reproducing the paper's numbers」を参照。

## Mathematical background

### 離散 Hodge 分解 (完全グラフ K_N)

エッジフローを 3 つの直交成分に分解：

```
f = grad(φ) + curl(ψ) + harmonic
```

- **gradient** ∇φ : 不可逆カスケード (上流→下流の順序)
- **curl** δψ : 循環フィードバック
- **harmonic** : K_N では自明に 0

Hodge potential φ が高い遺伝子ほど「構造的上流」。論文は下記の 3 つの φ を
用いる：

- `φ_static` : 健常ドナーのみから構築 → 健常時ネットワーク位相
- `φ_disease` : 疾患ドナーのみから構築 → 疾患時ネットワーク構造
- `φ_condition` : 両群の直接差分 → case-control 差

疾患特異的残差は `φ_disease` を `φ_static` の 3 次多項式で回帰した残差 z-score。

### Null GF baseline (Section 2.2, Supp Note 4)

```
GF_0 = 2 / [3 (1 + CV²)]
```

where CV = edge-weight coefficient of variation. Random-matrix null gives
CV ≈ 0.76 → GF_0 ≈ 0.425. sALS observed CV ≈ 0.79 → predicted GF ≈ 0.412;
observed GF = 0.398 indicates biological structure shifts flow energy toward
curl (|δGF| ≈ 0.014).

### Edge-weight flow (Methods 6.4, Supp Note 3)

```
f(i, j) = |Δ_ij| × sign(d_i − d_j),   d_i = ||Δ_{i,·}||_2
```

公理系 A1–A5 (antisymmetry, magnitude fidelity, node coherence, permutation
invariance, energy fidelity) の下で、この構造が唯一の energy-preserving
antisymmetric flow maximally aligned with the discrete gradient of d_i
(Supp Note 3, Theorem 3.3)。

## Data availability

- **sALS snRNA-seq**: NYGC ALS Consortium (Pineda et al. 2024, *Cell*)
  Synapse syn51105515, BioProject PRJNA1073234
- **Norman Perturb-seq**: GSE133344 / scPerturb
- **Shen spatial CRISPR**: GSE274058
- **TCGA + GTEx glioma**: UCSC Xena TOIL
- **TimeVault PC9**: [github.com/thechenlab/TimeVault](https://github.com/thechenlab/TimeVault)
- **Bryois eQTL**: Zenodo 7276971
- **van Rheenen ALS GWAS**: GCST90027164
- **GSE212630 (Part III cross-dataset)**: Wang, H.L.V. et al. 2023, *bioRxiv*
  2023.01.12.523820

## Citing

If you use this code, please cite:

```
Kaneko, S. & Urushitani, M. (2026). Hodge decomposition of gene co-expression
dynamics infers structural upstream organisation in sporadic ALS. [journal TBD].
```

And the companion theoretical manuscript:

```
Kaneko, A. (in preparation). IDS theory: stress geometry, swap-aggregation
dynamics, and principle-theoretic foundations of covariance manifold analysis.
```

## License

Research use.
