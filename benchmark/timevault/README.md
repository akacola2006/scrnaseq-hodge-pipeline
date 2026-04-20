# TimeVault TAS Benchmark

Reproduces paper Section 2.3: Temporal Asymmetry Score (TAS) on TimeVault
PC9 persister cell data (Chao et al. 2026, *Science*).

## Dataset

- **Source**: [TimeVault GitHub](https://github.com/thechenlab/TimeVault) Figure 3
- **Cells**: PC9 lung adenocarcinoma (EGFR exon 19 deletion)
- **Treatment**: Osimertinib (EGFR inhibitor) 100 nM, 4 days
- **Samples**: n = 12 bulk RNA-seq (4 conditions × 3 replicates)
  - NDc, ND (Non-drug Recorded, Non-drug Present)
  - Dc, D (Drug-treated Recorded, Drug-treated Present)
- **Genes**: 17,543 after filtering

## Ground truth classification (DEG-based)

| Category    | Count | Definition                            |
|-------------|-------|---------------------------------------|
| UPSTREAM    | 5,458 | Significant in Recorded only (Dc vs NDc) |
| DOWNSTREAM  |   547 | Significant in Present only (D vs ND) |
| BOTH        |   315 | Significant in both                   |
| NEITHER     |11,223 | Neither (background)                  |

## Paper claims (Section 2.3)

- **TAS ROC-AUC = 0.9475** (UPSTREAM vs DOWNSTREAM discrimination)
- **Cohen's d = 2.06** (very large effect)
- **F1-score = 0.93**
- **PCA temporal axis recovery AUC = 1.0** (Recorded vs Present) on n = 12
- **Treatment-effect p = 1.0** (drug vs control remains interleaved along the
  temporal axis, confirming temporal specificity)

## Temporal Asymmetry Score (TAS) definition

TAS is a per-gene univariate score (Section 6.11):

```
TAS_i = (|dmean_recorded_i| + |dvar_recorded_i|)
      − (|dmean_present_i| + |dvar_present_i|)

where:
  dmean_recorded = mean(Dc) − mean(NDc)   [drug effect in past RNA]
  dvar_recorded  = var(Dc)  − var(NDc)
  dmean_present  = mean(D)  − mean(ND)    [drug effect in present RNA]
  dvar_present   = var(D)   − var(ND)
```

**Interpretation**:
- TAS > 0 → Recorded (past) variance > Present variance → UPSTREAM candidate
- TAS < 0 → Present variance > Recorded variance → DOWNSTREAM candidate

TAS is NOT part of the IDS-Hodge pipeline. It is a premise-level readout
demonstrating that **variance asymmetry encodes temporal precedence** at the
per-gene univariate level. The IDS-Hodge pipeline uses this same operating
principle at the multivariate gene-gene network level.

## Scripts

### `run_tas.py`
Main TAS analysis: data loading, DEG analysis (t-test + BH-FDR), gene
classification, TAS computation, ROC-AUC / Cohen's d evaluation, and
4-panel validation figure.

## Reproducing the paper's numbers

```bash
# 1. Download TimeVault bulk RNA-seq data
curl -L -o metadata_PC9.csv \
    "https://raw.githubusercontent.com/thechenlab/TimeVault/main/Figure3/metadata_PC9.csv"
curl -L -o all_gex_matrix_genename.txt \
    "https://raw.githubusercontent.com/thechenlab/TimeVault/main/Figure3/all_gex_matrix_genename.txt"

# 2. Run TAS validation
python benchmark/timevault/run_tas.py \
    --counts all_gex_matrix_genename.txt \
    --metadata metadata_PC9.csv \
    --output-dir results/timevault_benchmark

# Expected output:
#   results/timevault_benchmark/brake_loss_validation_results.csv  (17,543 genes)
#   results/timevault_benchmark/gene_classification_timevault.csv  (4 categories)
#   results/timevault_benchmark/brake_loss_validation_figure.png   (4 panels)
#
#   Printed: ROC-AUC: 0.9475, Cohen's d: 2.06, F1: 0.93
```

## Interpretation notes (Section 2.3)

- **TAS validates the premise, not the IDS instrument itself.** Section 2.3
  explicitly clarifies: "Hodge decomposition of gene-gene covariance is not
  involved in its computation."
- **Biological context**: PC9 cancer cells + osimertinib is fundamentally
  different from sALS post-mortem brain. The result is a proof of principle
  that SPD-based pseudotime can recover true temporal ordering, not a direct
  validation of the instrument in the ALS context.
- **Statistical range limit**: AUC = 1.0 on n = 12 has limited power to
  distinguish "perfect" from "near-perfect" recovery.
- **Future experiment** (Section 2.3 end): direct comparison of IDS φ rankings
  with TimeVault temporal labels on matched longitudinal single-cell data.

## Source

Migrated from: `TimeVault/brake_loss_validation_pipeline.py`
