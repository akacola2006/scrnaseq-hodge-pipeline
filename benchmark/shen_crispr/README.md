# Shen Spatial CRISPR Perturb-seq Benchmark

Reproduces paper Section 2.5: Gene Hodge applied within each KO condition of
Shen et al. (2026) spatial Perturb-seq in mouse hippocampal neurons.

## Dataset

- **Source**: GEO GSE274058 (Shen et al. 2026, *Nat. Commun.* 17, 3018)
- **Cells**: 41,749 mouse hippocampal cells
- **Perturbations**: 18 conditions (17 single-gene KOs + 1 mSafe negative control)
- **Evaluated**: 16 targets (Ndufaf2 excluded for n_KO = 12; mSafe excluded as no target)

## Paper claims (Section 2.5)

- **Mean percentile**: 21.8 th (lower = more upstream)
- **Wilcoxon signed-rank p = 0.002**
- **ALS targets**: 12.8 th percentile (most upstream)
- **PD targets**: 18.2 th
- **AD targets**: 61.1 th (driven by Clu 78.5 th, Trem2 91.1 th — microglial effectors)
- **Sparse k=30**: p = 0.020
- **Re-execution with different seed**: mean 35.1 th, p = 0.031

## Design

For each KO condition's cells (n_KO ≥ 15), Gene Hodge is run independently on
the 316-gene edge-weight flow network. The KO target gene's φ percentile is
measured within its own perturbation network. This tests whether KO targets
occupy structurally foundational positions — i.e. whether the genes that
were removed anchor the post-perturbation co-expression structure.

## Scripts

### `perturbation_loader.py`
Loads GSE274058 data, assigns cells to perturbation barcodes, and handles
mouse → human ortholog mapping via `ortholog_map.py`.

### `run_shen_ids.py`
Track 2: per-perturbation Gene Hodge with Ledoit-Wolf pooled covariance,
edge-weight flow, bootstrap B=50, permutation B=500.

### `run_shen_ids_sparse.py`
Track 2 sparse-graph variant (k-NN with k ∈ {17, 30, 50}) for sensitivity analysis.
Produces the k=30 p=0.020 number cited in the paper.

### `track2_gene_hodge.py`
Core gene-level Hodge decomposition for perturbation data (different from
observational `scripts/gene_hodge.py` in that it uses pooled LW covariance
rather than per-donor).

## Reproducing the paper's numbers

```bash
# 1. Download GSE274058
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE274nnn/GSE274058/suppl/...

# 2. Preprocess (cell × gene matrix, barcode assignment)
python benchmark/shen_crispr/perturbation_loader.py \
    --input /path/to/GSE274058 \
    --output results/shen_preprocessed

# 3. Run per-perturbation IDS (main result)
python benchmark/shen_crispr/run_shen_ids.py \
    --data-dir results/shen_preprocessed \
    --output-dir results/shen_benchmark \
    --n-hvg 300 \
    --n-hvg-pool 5000 \
    --flow-mode edge_weight \
    --n-bootstrap 50 \
    --n-perm 500 \
    --min-ko-cells 15 \
    --ctrl-ratio 8

# Expected output:
#   results/shen_benchmark/track2_results.json
#   results/shen_benchmark/track2_per_perturbation/*.json (×17)
#   Summary: mean_percentile = 78.2 %  (→ 21.8 th in paper convention)
#   Wilcoxon p = 0.00222

# 4. Sparse variant for k=30 comparison
python benchmark/shen_crispr/run_shen_ids_sparse.py \
    --data-dir results/shen_preprocessed \
    --k-values 17 30 50 \
    --output-dir results/shen_benchmark
```

## Key result interpretation

**Percentile convention mismatch**: The pipeline outputs high = upstream;
the paper uses low = upstream (consistent throughout). A 78.2 % pipeline
output corresponds to 100 − 78.2 = **21.8 th percentile** in paper notation.

## Per-target results (from P4_TRACK2_RESULTS_REPORT.md)

| Gene     | Category     | n_KO | Paper %ile | Pipeline %ile | Tier   |
|----------|--------------|------|------------|---------------|--------|
| Tbk1     | ALS          | 27   | 2.5        | 97.5          | High   |
| C9orf72  | ALS          | 33   | 6.3        | 93.7          | Medium |
| Cfap410  | ALS          | 29   | 16.1       | 83.9          | Medium |
| Dpp6     | ALS          | 39   | 26.3       | 73.7          | Medium |
| Flcn     | Control      | 45   | 4.1        | 95.9          | Medium |
| Olig2    | Control      | 25   | 5.7        | 94.3          | Medium |
| Gfap     | Control      | 51   | 7.0        | 93.0          | Medium |
| Lrrk2    | PD           | 22   | 8.9        | 91.1          | Medium |
| Sh3gl2   | PD           | 22   | 10.4       | 89.6          | Medium |
| Srf      | Control      | 18   | 10.4       | 89.6          | Medium |
| Stk39    | PD           | 61   | 35.4       | 64.6          | Medium |
| Rraga    | Control      | 56   | 18.7       | 81.3          | Medium |
| Fasn     | Control      | 45   | 13.6       | 86.4          | Medium |
| Rbfox1   | AD           | 45   | 13.6       | 86.4          | Medium |
| Clu      | AD           | 75   | 78.5       | 21.5          | Low    |
| Trem2    | AD           | 38   | 91.1       | 8.9           | Low    |
| Ndufaf2  | AD           | 12   | excluded   | —             | —      |
| mSafe    | Ctrl-neg     | 20   | excluded   | —             | —      |

**Category means (paper convention, low = upstream)**:
- ALS: 12.8 %    |  PD: 18.2 %    |  Control: 9.9 %    |  AD: 61.1 %

## Source

Migrated from `sals_analysis_frozen_20260211/p4_perturbation_validation/scripts/`
