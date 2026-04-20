# Glioma (TCGA + GTEx) Cross-Disease Benchmark

Reproduces paper Section 2.6: Same IDS pipeline applied to TCGA glioma + GTEx
bulk RNA-seq, confirming null GF invariance across diseases and data modalities.

## Dataset

- **Source**: TCGA + GTEx via UCSC Xena TOIL (Vivian et al. 2017, Goldman et al. 2020)
- **Preprocessing**: log₂(TPM + 1)
- **Samples**: **723** after QC
  - W0: 105 samples (GTEx Cortex - Normal)
  - W1: 221 samples (TCGA LGG IDH-mut G2)
  - W2: 187 samples (TCGA LGG IDH-mut G3)
  - W3:  92 samples (TCGA LGG IDH-wt)
  - W4: 118 samples (TCGA GBM)
- **Genes**: 25,431 total, base set of 2,000 top-variance genes
- **Windows**: 5 clinical grades (Normal → LGG G2 → LGG G3 → LGG IDH-wt → GBM)

## Paper claims (Section 2.6)

- **Null GF invariance**: GF_null = **0.4248 ± 0.0004** (n = 2,000)
  matching sALS GF_null = 0.4250 (n = 922)
- **Data-driven pseudotime** vs clinical grade: Spearman ρ = **−0.7554** (reversed direction), AUC = **0.998**
- **φ agreement** between pseudotime-based and clinical-grade-based windows: ρ = **0.891**
- **Transition-specific upstream programmes** match established glioma biology:
  - Normal → LGG: ECM remodelling
  - G2 → G3: proliferation (cell cycle)
  - G3 → IDH-wt: synaptic acquisition (cf. Venkatesh et al. 2019)
  - IDH-wt → GBM: immune remodelling
- **Top upstream programmes**:
  - Cell cycle (G2 → G3)
  - Neuronal system (G3 → IDH-wt)
  - Vesicle transport (IDH-wt → GBM)
- **Myelin**: most upstream module overall (mean φ = +0.30)
- **Cell-cycle**: feedback-dominated internal structure
- **Caveat**: Major ribosomal genes (RPL3, RPL5, RPS6, …) absent from the
  2,000-gene base set due to bulk RNA-seq variance filtering; so the absence
  of translation from glioma GSEA reflects gene-set composition, not
  anti-artefact.

## Scripts

### `step8_pseudotime.py`
Build glioma pseudotime from 30-component PCA + diffusion map.
Root = highest-density donor. Output: `pt.csv` with per-sample PT value.

### `step9_gene_insertion.py`
All-gene insertion φ ranking with pseudotime-based windows.
Complementary to `step8_pseudotime.py`.

### `run_glioma_pseudotime_insertion.py`
Full GPU-accelerated all-gene PT insertion (integrates the above steps).
Output: `glioma_pt_allgene_insertion.csv` with 25,431 gene-level φ values.

### `run_idh_stratified_insertion.py`
IDH-mut (W2) vs IDH-wt (W3) condition-φ analysis.
Isolates the driver-mutation effect on co-expression geometry.

## Reproducing the paper's numbers

```bash
# 1. Download TCGA + GTEx expression data from UCSC Xena Toil
# URL: https://xenabrowser.net/datapages/
#   - TCGA LGG: TCGA.LGG.sampleMap/HiSeqV2
#   - TCGA GBM: TCGA.GBM.sampleMap/HiSeqV2
#   - GTEx Cortex: gtex_RSEM_Hugo_norm_count.txt.gz (Cortex samples)

# 2. Run pseudotime construction (requires combined expression matrix)
python benchmark/glioma/step8_pseudotime.py \
    --input combined_glioma_expression.csv \
    --output-dir results/glioma_pseudotime

# Expected:
#   results/glioma_pseudotime/pt.csv  (723 samples, PT values)
#   Spearman(PT, clinical_grade) = -0.7554

# 3. All-gene insertion
python benchmark/glioma/run_glioma_pseudotime_insertion.py \
    --expression combined_glioma_expression.csv \
    --pt results/glioma_pseudotime/pt.csv \
    --base-size 2000 \
    --output-dir results/glioma_insertion

# Expected:
#   results/glioma_insertion/glioma_pt_allgene_insertion.csv  (25,431 genes)
#   Null GF = 0.4248 ± 0.0004

# 4. IDH-stratified comparison
python benchmark/glioma/run_idh_stratified_insertion.py \
    --expression combined_glioma_expression.csv \
    --output-dir results/glioma_insertion
```

## Sample count reconciliation

The paper states **723 samples** (matching the actual run log). An earlier
draft of the validation report cites 721 due to a stricter (2-sample)
intermediate QC filter; the 723 figure corresponds to the published analysis.

## Source

Migrated from:
- `sals_analysis_frozen_20260211/scripts/glioma_step8_pseudotime.py`
- `sals_analysis_frozen_20260211/scripts/glioma_step9_gene_insertion.py`
- `sals_analysis_frozen_20260211/glioma_pseudotime_insertion/run_pseudotime_allgene_insertion_gpu.py`
- `sals_analysis_frozen_20260211/glioma_idh_stratified/run_idh_stratified_insertion_gpu.py`
