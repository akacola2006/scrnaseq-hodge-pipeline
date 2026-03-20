# scRNAseq Hodge Decomposition Pipeline

## Overview
This is a generalized single-cell RNA-seq (scRNAseq) analysis pipeline that uses
**discrete Hodge decomposition** to identify upstream cell types and upstream genes
driving disease or condition-associated changes in gene co-expression structure.

The pipeline was originally developed for ALS research and has been generalized to
work with any scRNAseq dataset comparing conditions (e.g., disease vs. control).

## How to Run This Pipeline

### Step 1: Prepare Data
Place your data in the `data/` folder:
- `data/h5ad/` — One h5ad file per sample/donor (see Data Format below)
- `data/metadata/sample_info.csv` — Sample metadata (see below)
- `data/metadata/gene_annotation.csv` — (Optional) Gene annotation with chromosome info

### Step 2: Edit Configuration
Edit `project_config.yaml` to set:
- Cell type column name and list of cell types
- Condition labels (disease vs. control)
- Hardware settings (GPU/CPU)
- Analysis parameters

### Step 3: Run the Pipeline
```bash
python run_pipeline.py
```
Or run individual steps:
```bash
python run_pipeline.py --step residuals
python run_pipeline.py --step pca
python run_pipeline.py --step pseudotime
python run_pipeline.py --step lane_a
python run_pipeline.py --step bootstrap
python run_pipeline.py --step lane_b
python run_pipeline.py --step gene_hodge
python run_pipeline.py --step enrichment
python run_pipeline.py --step dual_mode
```

Additional analysis modules (import and use directly in Python):
```python
from scripts.directional import run_directional_decomposition   # Δ⁺/Δ⁻ split
from scripts.multi_transition import run_multi_transition        # All-transition integration
from scripts.random_baseline import run_random_baseline          # GF null test
from scripts.two_axis import run_two_axis_model                  # TRS × MSS trajectory
```

## Data Format Requirements

### h5ad Files (one per sample/donor)
Each h5ad file must contain:
- `.X` — Raw count matrix (cells x genes), sparse or dense
- `.obs` — Cell-level metadata with at minimum:
  - A cell type column (name configurable in `project_config.yaml`, default: `CellType`)
  - (Optional) `Condition` column per cell
- `.var` or `.var_names` — Gene identifiers (ENSEMBL IDs or gene symbols)

**Naming convention**: Files should be named consistently. The mapping from
filename to donor/sample ID is defined in `data/metadata/sample_info.csv`.

### sample_info.csv (Required)
CSV with columns:
| Column     | Description                           | Example          |
|------------|---------------------------------------|------------------|
| donor_id   | Unique sample/donor identifier        | "D001"           |
| file       | h5ad filename (basename only)         | "sample1.h5ad"   |
| condition  | Condition label (disease/control)     | "ALS" or "Control"|
| (optional) | Any additional metadata columns       |                  |

### gene_annotation.csv (Optional)
CSV with columns:
| Column     | Description                | Example            |
|------------|----------------------------|--------------------|
| gene_id    | Gene identifier            | "ENSG00000198888"  |
| chromosome | Chromosome location        | "chrM"             |
| gene_name  | Gene symbol (optional)     | "MT-ND1"           |

If not provided, mitochondrial gene detection defaults to gene names starting with "MT-".

### Pre-included Data Files

The following reference files are already included in `data/metadata/`:

- **`gene_annotation.csv`** — Comprehensive gene annotation (~60,000 genes)
  - Columns: gene_id (ENSEMBL), gene_symbol, gene_type, chromosome,
    transcript_length, neuron_or_glia, neuron (expression), glia (expression),
    gc_ratio, g4_score, tdp43_prob, ftt_prob, mirna_score, mirna_hit_no
  - Used for: mitochondrial gene detection, ENSEMBL-to-symbol mapping,
    neuron/glia classification, gene-level annotation of results

- **`functional_modules.json`** — 24 functional gene modules for enrichment analysis
  - Modules: Synaptic, Oxidative_Stress, Protein_Homeostasis, RNA_Processing,
    Apoptosis, Autophagy, Inflammation, Calcium_Signaling, Ion_Transport,
    Mitochondria, ER_Stress, ALS_Genes, DNA_Repair, Cell_Cycle, Cytoskeleton,
    Metabolism, ECM, Angiogenesis, Growth_Factors, Transcription, Epigenetic,
    Complement, Myelination, lncRNA
  - Used for: Fisher's exact test enrichment of upstream gene sets

## Pipeline Architecture

```
h5ad files (per donor)
  |
  v
[1. Data Loading & QC] — Filter cells by type, compute nUMI, pct_mito
  |
  v
[2. Normalization] — log1p(CPM)
  |
  v
[3. Residual Computation] — Regress out: donor + nUMI + pct_mito (GPU-accelerated)
  |
  v
[4. PCA] — Randomized SVD, k components (GPU-accelerated)
  |
  v
[5. SPD Covariance] — Per-donor Ledoit-Wolf covariance on PC scores
  |
  v
[6. Pseudotime (PT-B)] — Diffusion map on SPD covariance structure
  |
  v
[7. Window Assignment] — Partition donors into K time windows
  |
  v
=== LANE A: Cell-type Level ===
[8. Window Representatives] — Log-Euclidean mean SPD per window
[9. Distance Decomposition] — d_cov, d_corr, d_var between consecutive windows
[10. Pairwise Flow] — Build edge flow on cell-type complete graph
[11. Hodge Decomposition] — Decompose flow into gradient + curl + harmonic
[12. Cell-type phi] — Upstream score (gradient potential per cell type)
  |
  v
=== BOOTSTRAP (100 iterations) ===
[13. Donor Resampling] — Cluster bootstrap at donor level
[14. Stability Assessment] — Top-3 reproducibility, phi 95% CI
  |
  v
=== LANE B: Gene Level (within upstream cell type) ===
[15. Upstream PC Identification] — Top PCs by d_corr contribution
[16. Gene Set Extraction] — Top 200 genes per upstream PC
[17. Gene-level Hodge] — Correlation-change flow on gene complete graph
[18. Gene phi Scores] — Upstream gene ranking
[19. Permutation Testing] — Statistical significance
  |
  v
[20. Rank Shift Analysis] — Disease vs. control gene ranking comparison
  |
  v
=== DUAL-MODE ===
[21. K_N vs Sparse k-NN] — Compare complete graph vs k-NN sparsified Hodge
[22. Curl Hub Analysis] — Triangle curl scores on sparse graph
  |
  v
=== ENRICHMENT ===
[23. Functional Module Enrichment] — Fisher's exact test against 24 modules
[24. Gene Annotation] — ENSEMBL/symbol mapping, neuron/glia, TDP-43 prob
```

### Additional Analysis Modules (standalone)

These can be run independently after the main pipeline completes:

- **Directional Decomposition** (`scripts/directional.py`)
  Split Δ into Δ⁺ (correlation gain) and Δ⁻ (correlation loss), run Hodge on
  each separately. Reveals whether upstream genes drive coordinated activation
  or decoupling.

- **Multi-Transition Integration** (`scripts/multi_transition.py`)
  Compute phi across ALL window transitions (not just the primary), then
  integrate with 4 weighting schemes (uniform, GF-weighted, rank-median, best).
  Concordance matrix shows stability across the pseudotime axis.

- **Random Matrix Baseline** (`scripts/random_baseline.py`)
  Generate null GF distribution from random symmetric matrices (same N).
  Tests whether observed GF is genuine signal or structural artifact of K_N.
  Critical for interpretation: if observed GF ≤ null p95, the gradient signal
  is not significant.

- **2-Axis Model** (`scripts/two_axis.py`)
  TRS (Translation-Resource Score) × MSS (Morphogenesis-Structure Score) per
  window. TRS = coordination among High-phi genes; MSS = PC1 of Medium-phi genes.
  If TRS peaks before MSS → supply-demand cascade structure.

## Key Parameters (project_config.yaml)

| Parameter               | Default | Description                                    |
|-------------------------|---------|------------------------------------------------|
| pca_k                   | 100     | Number of PCA components                       |
| n_windows               | 8       | Number of pseudotime windows                   |
| min_cells_per_donor     | 100     | Minimum cells per donor per cell type           |
| gene_detection_min_frac | 0.01    | Minimum fraction of cells detecting a gene     |
| bootstrap_n             | 100     | Number of bootstrap iterations                 |
| bootstrap_top_k         | 3       | Top-K cell types for reproducibility check     |
| bootstrap_stable_thresh | 0.80    | Stability threshold (>=80%)                    |
| gene_hodge_perm_n       | 1000    | Gene Hodge permutation count                   |
| control_label           | varies  | Label for control/reference condition           |

## Dependencies
```
anndata>=0.10
scanpy>=1.9
numpy>=1.24
scipy>=1.10
pandas>=2.0
scikit-learn>=1.3
statsmodels>=0.14
torch>=2.0          # Optional: GPU acceleration
pyyaml>=6.0
```

## Output Structure
```
results/
  run_YYYYMMDD_HHMM/
    logs/                    — Pipeline logs and failure journal
    pt/                      — Pseudotime assignments
    laneA/                   — Cell-type level Hodge results
      hodge_quality/         — Gradient fraction, permutation tests
    laneB/                   — Gene-level results
      gene_hodge/            — Gene phi scores per cell type
      rank_shift/            — Condition comparison results
      enrichment/            — Functional module enrichment results
      dual_mode/             — K_N vs sparse k-NN comparison
      directional/           — Δ⁺/Δ⁻ directional decomposition
      multi_transition/      — All-transition phi integration
      random_baseline/       — Null GF distribution
      two_axis/              — TRS × MSS trajectory
    bootstrap/               — Bootstrap confidence estimates
    figures/                 — Visualization outputs
    qc/                      — Quality control reports
    summary.json             — Run summary with all key results
```

## For Claude Code Users
When the user asks you to analyze scRNAseq data with this pipeline:
1. Check that `data/h5ad/` contains h5ad files
2. Check that `data/metadata/sample_info.csv` exists and is correctly formatted
3. Read `project_config.yaml` and verify/update settings for the user's data
4. Run `python run_pipeline.py` to execute the full pipeline
5. Interpret results from `results/` and summarize findings

When troubleshooting:
- Check `results/*/logs/failures.jsonl` for error details
- Verify h5ad files can be loaded: `python -c "import anndata; print(anndata.read_h5ad('data/h5ad/FILE.h5ad'))"`
- Verify GPU availability: `python -c "import torch; print(torch.cuda.is_available())"`
