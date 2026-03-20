# Data Directory

Place your scRNAseq data here for analysis.

## Required Structure

```
data/
  h5ad/                      <- Place h5ad files here (one per sample/donor)
    sample1.h5ad
    sample2.h5ad
    ...
  metadata/
    sample_info.csv           <- Required: sample metadata
    gene_annotation.csv       <- Optional: gene annotations
```

## h5ad File Format

Each h5ad file should contain data for one sample (donor/patient):

- **`.X`** — Raw count matrix (cells x genes)
  - Can be sparse (CSR/CSC) or dense
  - Integer counts preferred (not normalized)
- **`.obs`** — Cell-level metadata DataFrame
  - Must include a cell type column (default name: `CellType`)
  - Cell type labels should be consistent across all samples
- **`.var`** or **`.var_names`** — Gene identifiers
  - ENSEMBL IDs (e.g., `ENSG00000198888`) or gene symbols (e.g., `MT-ND1`)
  - Must be consistent across all samples

### Example h5ad structure:
```python
import anndata as ad

adata = ad.read_h5ad("sample1.h5ad")
print(adata)
# AnnData object with n_obs x n_vars = 5000 x 30000
#     obs: 'CellType', 'nCount_RNA', ...
#     var: 'gene_name', ...

print(adata.obs["CellType"].value_counts())
# Oligo     1200
# Astro      800
# Micro      600
# ...
```

## sample_info.csv Format

| donor_id | file           | condition |
|----------|----------------|-----------|
| D001     | sample1.h5ad   | Disease   |
| D002     | sample2.h5ad   | Disease   |
| D003     | sample3.h5ad   | Control   |
| ...      | ...            | ...       |

- **donor_id**: Unique identifier for each sample (string)
- **file**: Filename of the h5ad file (basename only, not full path)
- **condition**: Condition label (must match `conditions.control_label` and
  `conditions.disease_labels` in `project_config.yaml`)
- Additional columns are allowed and will be preserved

## gene_annotation.csv Format (Optional)

| gene_id           | chromosome | gene_name |
|-------------------|------------|-----------|
| ENSG00000198888   | chrM       | MT-ND1    |
| ENSG00000198763   | chrM       | MT-ND2    |
| ENSG00000141510   | chr17      | TP53      |

If not provided, mitochondrial genes are identified by names starting with "MT-".

## Recommended Dataset Properties

- **Minimum samples**: 10+ donors (20+ recommended for robust pseudotime)
- **Minimum cells per donor per cell type**: 100 (configurable)
- **Cell types**: 3+ cell types for meaningful Hodge decomposition
- **Conditions**: At least 2 conditions (disease + control)
- **Genes**: Standard whole-transcriptome (~20,000-30,000 genes)

## Common Data Sources

This pipeline works with scRNAseq data from:
- 10x Genomics Chromium (Cell Ranger output -> convert to h5ad)
- Drop-seq
- Smart-seq2
- snRNA-seq (single-nucleus)

### Converting Cell Ranger output to h5ad:
```python
import scanpy as sc

# From filtered_feature_bc_matrix
adata = sc.read_10x_h5("filtered_feature_bc_matrix.h5")
# or
adata = sc.read_10x_mtx("filtered_feature_bc_matrix/")

# Add cell type annotations (from your clustering/annotation)
adata.obs["CellType"] = your_celltype_labels

# Save
adata.write_h5ad("sample1.h5ad")
```
