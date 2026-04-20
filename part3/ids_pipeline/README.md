# IDS Pipeline

**Unified φ-Flow Energy Analysis Pipeline for ALS Motor Cortex Data**

## Quick Start

```python
from ids_pipeline import IDSPipeline

# Initialize pipeline (auto-detects config files)
pipeline = IDSPipeline()

# Run full analysis
results = pipeline.run(
    expression_file='expression.csv',
    metadata_file='metadata.csv',
    output_dir='results/'
)
```

## Command Line Usage

```bash
# Basic usage
python run_pipeline.py expression.csv metadata.csv -o results/

# With custom config
python run_pipeline.py expression.csv metadata.csv -o results/ -c config.yaml

# With column specifications
python run_pipeline.py expression.csv metadata.csv -o results/ \
    --condition-column condition \
    --celltype-column cell_type \
    --pt-column PT_dpt

# Show configuration
python run_pipeline.py --show-config

# Generate config template
python run_pipeline.py --generate-config my_config.yaml
```

## Features

### Automatic Configuration
- Auto-detects gene mapping files (`gene_mapping_enhanced.json`)
- Auto-detects module files (`24_functional_modules_fixed.json`)
- Auto-maps cell types to NVU groups (Vascular, Glia, Excitatory, Inhibitory, VAT1L)

### Gene Mapping
```python
from ids_pipeline import GeneMapper

mapper = GeneMapper()

# Convert single gene
mapper.symbol_to_ensg('STMN2')  # -> 'ENSG00000104435'
mapper.ensg_to_symbol('ENSG00000104435')  # -> 'STMN2'

# Convert lists
mapper.convert_genes(['STMN2', 'TARDBP'], to='ensg')
```

### Module Scoring
```python
from ids_pipeline import ModuleScorer

scorer = ModuleScorer()

# Score expression matrix
scores = scorer.score_expression(expression_df)

# Get module genes
genes = scorer.get_module_genes('Mitochondria', as_ensg=True)

# List available modules
modules = scorer.list_modules()
```

### φ-Flow Analysis
```python
# Step-by-step analysis
pipeline = IDSPipeline()
pipeline.load_data('expression.csv', 'metadata.csv')
pipeline.compute_modules()
pipeline.compute_phi()
pipeline.analyze_flow()
pipeline.save_results('results/')
```

## Configuration

### YAML Configuration File

```yaml
# config.yaml
n_bins: 25
phi_threshold: 0.5
min_cells: 10

nvu_mapping:
  Vascular:
    - Vasc.Endo
    - Vasc.Fibro
  Glia:
    - Glia.Astro
    - Glia.Oligo
  # ...

stress_components:
  TDP43_targets:
    weight: 0.30
    genes:
      - STMN2
      - UNC13A
      # ...
```

### Python Configuration

```python
from ids_pipeline import IDSConfig

# Load from YAML
config = IDSConfig.from_yaml('config.yaml')

# Or create programmatically
config = IDSConfig(n_bins=30, phi_threshold=0.3)

# Validate
issues = config.validate()

# Print summary
config.print_summary()
```

## File Structure

```
ids_pipeline/
├── __init__.py          # Package init
├── config.py            # Configuration management
├── gene_utils.py        # Gene ID conversion
├── module_utils.py      # Module scoring
├── pipeline.py          # Main pipeline class
├── run_pipeline.py      # CLI runner
├── README.md            # This file
└── templates/
    └── default_config.yaml  # Config template
```

## Input Data Format

### Expression Matrix
- CSV or gzipped CSV
- Genes as rows (index), Cells as columns
- Gene IDs can be ENSG or Symbols (auto-detected)

### Metadata
- CSV file
- Required columns:
  - `cell_id`: Cell identifier (matching expression columns)
  - `condition`: Control/ALS status
  - `cell_type`: Cell type annotation
- Optional columns:
  - `PT_dpt` or similar: Pre-computed pseudotime

## Output Files

- `module_scores.csv`: Module expression scores per cell
- `phi_values.csv`: φ energy values per cell per module
- `flow_profiles.csv`: φ values aggregated by NVU group and PT bin
- `flow_summary.csv`: Onset and peak timing for each module
- `config.yaml`: Configuration used for analysis

## 23 Functional Modules

1. **Metabolic Axis**: Mitochondria, ER_Stress, Protein_Homeostasis, Metabolism
2. **Hyperexcitability Axis**: Calcium_Signaling, Ion_Transport, Synaptic
3. **Inflammatory Axis**: Oxidative_Stress, Inflammation, Complement
4. **Structural**: ECM, Cytoskeleton, Myelination, Angiogenesis
5. **Cellular Processes**: Apoptosis, Autophagy, Cell_Cycle, DNA_Repair
6. **Regulatory**: Epigenetic, Growth_Factors, Transcription, RNA_Processing, lncRNA

## Dependencies

- Python 3.8+
- pandas
- numpy
- PyYAML
- scipy (optional, for some analyses)
- scikit-learn (optional, for axis discovery)

## Related Files

- `scripts/ids_core.py`: Core IDS functions (φ computation, binning, etc.)
- `config/24_functional_modules_fixed.json`: Module gene definitions
- `gene_mapping_enhanced.json`: ENSG↔Symbol mapping

## Author

Claude Code + User
Date: 2025-12-11
