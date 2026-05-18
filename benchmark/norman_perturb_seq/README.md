# Norman Perturb-seq Benchmark

Reproduces paper Section 2.4: IDS perturbation target identification on
Norman et al. (2019) CRISPRa Perturb-seq in K562 cells.

## Dataset

- **Source**: [scPerturb repository](https://scperturb.org/) / GEO GSE133344
- **Cells**: 111,122 K562 cells (test set: 34,624)
- **Perturbations**: 235 CRISPR perturbation classes
  (232 groups evaluated after excluding 3 classes with n_cells < 3)
- **Candidate targets**: 102 genes

## Paper claims (Section 2.4 Table)

| Method                    | Top-1     | Median Rank | Training | Perturbation Labels |
|---------------------------|-----------|-------------|----------|---------------------|
| Random                    | 0.98 %    | 51.5        | —        | —                   |
| GRNBoost2 (best of 3)     | 7.3 %     | —           | None     | Not used            |
| d_corr (whitened)         | 57.3 %    | 1.0         | None     | Not used            |
| **IDS (whitened + Hodge)**| **66.4 %**| **1.0**     | **None** | **Not used**        |
| DE ranking                | 81.0 %    | 1.0         | None     | Required            |

- IDS vs d_corr: McNemar p = 1.1 × 10⁻⁵ (+9.1 pp improvement)
- IDS vs DE: McNemar p = 2.1 × 10⁻⁶
- Union IDS + DE Top-1: 85.3 %

## Scripts

### `run_ids_norman.py`
Main IDS benchmark on Norman et al. 2019 K562 Perturb-seq, paper test
split (`norman_test.h5ad`, 34,624 cells, 235 categories). Configuration:
- **R2b mode**: global whitening (W = Σ_ctrl^{-1/2}) + bootstrap (B = 100)
- Per-group analytical K_N Hodge phi (edge-weight flow on the candidate
  gene complete graph)
- Dynamic candidate construction from the perturbation labels (tokens of
  `cat.split("_")` that match `adata.var_names`), matching the
  NormanCellNaviDataset loader.
- Reproduces the **66.4 % Top-1** reference value (Section 2.4).

### `run_ids_replogle.py`
Independent external replication of the same R2b pipeline on Replogle et
al. 2022 K562 essential Perturb-seq
(`K562_essential_raw_singlecell_01.h5ad`, ~310k cells). Reported as a
robustness check that the IDS / whitening / Hodge pipeline is not
dataset-specific. Same algorithm as `run_ids_norman.py`; only the data
loader and (large, ~2k) candidate gene set differ.

### `run_grnboost2_norman.py`
GRNBoost2 baseline using `sklearn.ensemble.GradientBoostingRegressor`
with the same hyperparameters as the SCENIC pipeline
(n_estimators=100, max_features='sqrt', learning_rate=0.05, subsample=0.9, max_depth=3).
Three aggregation strategies are tested; best strategy ("max") yields Top-1 = 7.3 %.

### Requires `scripts/whitening.py`
ZCA whitening implementation from `scripts/whitening.py`:
    W = V D^{-1/2} V^T  where V, D are eigendecomposition of Ledoit-Wolf Σ_ctrl
    (eigenvalues clamped at 10⁻¹⁰ for numerical stability).

## Reproducing the paper's numbers

```bash
# 1. Download Norman 2019 data (requires scPerturb account or GEO download)
# See https://scperturb.org/ ; extract the paper test split as
# data/external/norman/processed/norman_test.h5ad

# 2. Point the script at the data (default location works if the file is
#    at the path above; otherwise override with ALS_NORMAN_DATA_DIR or
#    ALS_NORMAN_TEST_H5AD).
export ALS_NORMAN_TEST_H5AD=/path/to/norman_test.h5ad
export ALS_NORMAN_OUTPUT_DIR=results/norman_benchmark

# 3. Run IDS benchmark on Norman
python benchmark/norman_perturb_seq/run_ids_norman.py

# Expected output:
#   $ALS_NORMAN_OUTPUT_DIR/ids_norman.json           # summary including top1_rate ~ 0.66
#   $ALS_NORMAN_OUTPUT_DIR/ids_norman_per_group.csv  # 232 rows, one per perturbation

# 4. Run GRNBoost2 baseline (optional comparison)
python benchmark/norman_perturb_seq/run_grnboost2_norman.py

# 5. Optional: Replogle external replication (independent dataset)
export ALS_REPLOGLE_DATA_DIR=/path/to/replogle
python benchmark/norman_perturb_seq/run_ids_replogle.py
```

## Key results (from `champion_tournament_final_summary.md`)

- **IDS R2b Top-1: 66.4 %** (Section 2.4 canonical value)
  Equivalent alternative run (232-group per-target evaluation) yields 65.9 %;
  both numbers refer to the same IDS R2b configuration with minor
  bootstrap/control-subsample differences.
- Failure modes (~34 % of groups):
  - Diffuse effectors (MAPK1)
  - Indirect effectors (ARID1A)
  - Non-transcriptional effectors (ATL1)
  - Epistatic masking (AHR in AHR_KLF1)

## Robustness checks (Section 2.4)

- Null test (random 102 genes): GF = 0.000 (pipeline has no bias)
- Dilution test (102 → 300 candidates): median rank = 1.0 retained
- Open-world whitening (2,102 genes): Top-1 = 68.5 %
- Strict pseudotime evaluation: 74.8 % at rank ≤ 5 (15.3-fold enrichment)
- Pseudotime concordance: 95.1 %

## Source

Migrated from:
- `cell navi/benchmark/run_centrality_baselines_norman_v3.py`  (BV.61 v3,
  Norman test split + dynamic candidate construction)
- `cell navi/benchmark/run_ids_replogle_r2b.py`  (R2b whitening + Hodge
  phi pipeline)
  → combined into `run_ids_norman.py`
- `cell navi/benchmark/run_ids_replogle_r2b.py` (Replogle data loader)
  → `run_ids_replogle.py`
- `cell navi/benchmark/run_grnboost2_norman.py` → `run_grnboost2_norman.py`
- `sals_analysis_frozen_20260211/scripts/whitening.py` → `scripts/whitening.py`

The previous repository state placed the R2b Replogle script under
`run_ids_norman.py`; that script has been renamed to `run_ids_replogle.py`
and a true Norman evaluation script has been added under the original
filename.
