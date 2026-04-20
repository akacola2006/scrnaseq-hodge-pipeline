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
Main IDS benchmark. Configuration:
- **R2b mode**: global whitening (W = Σ_ctrl^{-1/2}) + bootstrap (B = 100)
- Per-group Hodge decomposition on K₁₀₂ graph (102 candidate genes)
- Reproduces the **66.4 % Top-1** reference value.

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
# See https://scperturb.org/

# 2. Run IDS benchmark
python benchmark/norman_perturb_seq/run_ids_norman.py \
    --data-dir /path/to/norman_2019 \
    --output-dir results/norman_benchmark \
    --mode R2b  # whitening + bootstrap

# Expected output:
#   results/norman_benchmark/top1_accuracy.txt  → 66.4 %
#   results/norman_benchmark/per_group_results.csv  (232 rows)
#   results/norman_benchmark/mcnemar_vs_d_corr.json  → p = 1.1×10⁻⁵

# 3. Run GRNBoost2 baseline (optional comparison)
python benchmark/norman_perturb_seq/run_grnboost2_norman.py \
    --data-dir /path/to/norman_2019 \
    --output-dir results/norman_benchmark
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
- `cell navi/benchmark/run_ids_replogle_r2b.py` → `run_ids_norman.py`
- `cell navi/benchmark/run_grnboost2_norman.py` → `run_grnboost2_norman.py`
- `sals_analysis_frozen_20260211/scripts/whitening.py` → `scripts/whitening.py`
