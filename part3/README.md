# Part III: Module-Level Causal Analysis — ATM Brake Loss / NVU Trigger Hypothesis

Reproduces paper **Sections 9–11** (Part III): an independent analytical
pipeline applied to the same NYGC ALS Consortium snRNA-seq dataset used in
Parts I–II, extending to 43 cell types and 23 functional modules with
diffusion pseudotime (PT_dpt), NOTEARS + LiNGAM causal inference, and
cross-dataset concordance against GSE212630 (Wang et al. 2023).

## Directory structure

```
part3/
├── README.md                               ← this file
├── ALS_TRIGGER_METHODOLOGY.md              ← narrative methodology (older summary)
│
├── ids_pipeline/                           ← Main Python API
│   ├── __init__.py
│   ├── config.py                           Auto-detect config file
│   ├── gene_utils.py                       Gene symbol ↔ ENSG mapping
│   ├── module_utils.py                     Module scoring (23 / 24 modules)
│   ├── pipeline.py                         `IDSPipeline` class (orchestrator)
│   ├── run_pipeline.py                     CLI: `python run_pipeline.py expr.csv meta.csv`
│   ├── templates/default_config.yaml       Config template
│   └── README.md                           Pipeline usage
│
├── scripts/                                ← Analysis scripts (52 files)
│   ├── ids_core.py                         ★ Core library (phi energy, PT binning,
│   │                                         onset/peak detection, axis discovery)
│   ├── 01_composite_stress.py              Composite stress score
│   ├── 02a_compute_module_scores.py        Module score computation
│   ├── 03a_prepare_notears_prior.py        NOTEARS causal inference (preparation)
│   ├── 03b_run_notears.py                  NOTEARS execution
│   ├── 03c_lingam_refinement.py            ★ LiNGAM refinement (paper Section 9.1)
│   ├── 05a–05d_*.py                        Cell state / state-level causal inference
│   ├── 05p1–05p4_*.py                      Patient stratification
│   ├── Phase5pp_01_subtype_upstream_validation.py
│   ├── Phase6_ext_driver_validation.py
│   ├── Phase7_functional_driver_advanced_validation.py
│   ├── Phase8_PT_dpt_time_lag_analysis.py  ★ Phase 8 time-lag (paper)
│   ├── Phase8b_sensitivity_and_trigger_module_analysis.py
│   ├── Phase9a_*.py                        PT_dpt recalculation
│   ├── Phase9cd_vascular_timelag_correlation_analysis.py ★ Vascular time-lag
│   ├── Phase9p_vascular_driver_state_identification.py ★ Vascular driver state
│   ├── Phase9pp_patient_vascular_driver_coupling.py
│   ├── Phase10_NVU_causal_DAG.py           ★ NVU causal DAG (Section 10)
│   ├── Phase11_NVU_nonlinear_coupling.py
│   ├── Phase12_create_NVU_integrated_figure.py
│   ├── Phase13_NVU_energy_landscape.py     ★ Phase 13 NVU energy landscape
│   ├── Phase13_complete_module_analysis.py
│   ├── Phase13c_NVU_energy_REFACTORED.py
│   ├── Phase_STMN2_*.py                    STMN2 single-gene analyses
│   ├── Phase_TDP_*.py                      TDP-43 analyses
│   ├── PTv2_01_build_PT_dpt.py             PT_dpt construction (robustness)
│   ├── PTv2_02_build_PT_stress.py          PT_stress (alternative axis)
│   ├── PTv2_03_compare_robustness.py
│   ├── GSE212630_*.py (9 files)            ★ Cross-dataset validation (Wang 2023)
│   ├── gse212630_triggerstage_v1/v2/v3.py  ★ Trigger-stage analysis
│   ├── analyze_GSE212630_with_ids_core.py
│   ├── fire_origin_landscape.py            ★ Fire origin identification
│   ├── fire_origin_targeted.py
│   ├── final_closure_mmc.py                MMC cascade closure
│   ├── trigger_network_discovery.py
│   ├── trigger_gene_network.py
│   └── test_ids_core_*.py                  Unit tests (3 files)
│
├── reports/                                ← Phase reports (24 files)
│   ├── 24_functional_modules.json          Module definitions (canonical, 24)
│   ├── 24_functional_modules_fixed.json    Corrected version
│   ├── PHASE9_VASCULAR_ORIGIN_REPORT.md    ★ Core capillary-endothelium origin report
│   ├── PHASE8_TIME_LAG_REPORT.md
│   ├── PHASE8_EXECUTIVE_SUMMARY.md
│   ├── PHASE8B_SENSITIVITY_REPORT.md
│   ├── PHASE11_NVU_NONLINEAR_COUPLING_REPORT_CORRECTED.md
│   ├── PHASE12_NVU_INTEGRATED_MODEL.md
│   ├── PHASE13B_EXTENDED_LANDSCAPE.md
│   ├── PHASE13C_COMPLETE_LANDSCAPE.md
│   ├── PHASE13_NVU_ENERGY_LANDSCAPE.md
│   ├── PHASE6_FUNCTIONAL_DRIVER_FINAL_REPORT.md
│   ├── PHASE7_ADVANCED_VALIDATION_REPORT.md
│   ├── PHASE9PP_PATIENT_COUPLING_REPORT.md
│   ├── SUBTYPE_UPSTREAM_DOWNSTREAM_VALIDATION.md
│   ├── GSE212630_ANALYSIS_SUMMARY.md       ★ Cross-dataset summary
│   ├── GSE212630_PHI_FLOW_ANALYSIS_REPORT.md
│   ├── PATIENT_STRATIFIED_INTEGRATED_REPORT.md
│   ├── PTv2_VALIDATION_REPORT.md           PT_dpt v2 robustness validation
│   ├── PT_DPT_STRESS_TRAJECTORY_REPORT.md
│   ├── COMPLETE_ANALYSIS_SUMMARY.md
│   ├── FINAL_INTEGRATED_ANALYSIS_SUMMARY_REVISED.md
│   ├── REVISION_SUMMARY.md
│   └── INTEGRATED_SUMMARY_FIRST_TRIGGER_MODULES.md
│
└── docs/                                   ← Project-level narrative (18 files)
    ├── ATM_fire_origin_evidence_summary.md ★ Complete ATM evidence
    ├── CELL_STATE_CAUSAL_HIERARCHY.md
    ├── 26_MODULES_COMPLETE_LIST.md
    ├── ALS_Cortex_IDS_Network_Paper.md     ★ Full paper narrative draft
    ├── ALS_TRIGGER_PROJECT_HYBRID_FINAL.md
    ├── COMPREHENSIVE_PAPER.md
    ├── COMPREHENSIVE_PAPER_DETAILED.md
    ├── COMPREHENSIVE_PROJECT_REPORT.md
    ├── CELL_TIMELINE_GUIDE.md
    ├── ALL_CELLTYPES_LIST.md               All 43 cell types listed
    ├── CONTINUOUS_PT_ANALYSIS_SUMMARY.md
    ├── CROSS_VALIDATION_WITH_ORIGINAL_PROJECT.md
    ├── GENE_MODULE_DEFINITIONS.md          Module definitions narrative
    ├── DATA_STRUCTURE_AND_PIPELINE_GUIDE.md
    ├── METHODOLOGY_COMPARISON_CORRECTED.md
    ├── DATA_DRIVEN_INTERPRETATION_ALS_vulnerability.md
    ├── PHASE5_COMPLETE_SUMMARY.md
    └── LAYERED_NETWORK_REPORT.md
```

## Paper claims and corresponding scripts

| Paper claim | Primary script | Supporting report |
|---|---|---|
| 43 cell types × 23 modules × PT_dpt (Section 9.1) | `Phase9a_recalculate_PT_dpt_all_cells.py` | `PHASE9_VASCULAR_ORIGIN_REPORT.md` |
| Stress-independent PT_dpt (R² = 1.9%) | `PTv2_01_build_PT_dpt.py` + `PTv2_03_compare_robustness.py` | `PTv2_VALIDATION_REPORT.md`, `PT_DPT_STRESS_TRAJECTORY_REPORT.md` |
| LiNGAM causal inference (Section 9.1) | `03a/b/c_*.py` | `PHASE11_NVU_NONLINEAR_COUPLING_REPORT_CORRECTED.md` |
| Phase 8 time-lag analysis | `Phase8_PT_dpt_time_lag_analysis.py` | `PHASE8_TIME_LAG_REPORT.md` |
| ATM pathway 3-criteria evaluation (Section 10.2) | `Phase9p_vascular_driver_state_identification.py` | `PHASE9_VASCULAR_ORIGIN_REPORT.md`, `ATM_fire_origin_evidence_summary.md` |
| ATM p = 0.018, 89% pathway coherence | `Phase9cd_vascular_timelag_correlation_analysis.py` | same |
| Capillary 881 cells + Pericyte 511 cells (Section 10.6) | `Phase9a_recalculate_PT_dpt_all_cells.py` | same |
| GSE212630 cross-dataset 83% concordance (Section 10) | `GSE212630_compute_PT_dpt_and_validate.py`, `analyze_GSE212630_with_ids_core.py` | `GSE212630_ANALYSIS_SUMMARY.md`, `GSE212630_PHI_FLOW_ANALYSIS_REPORT.md` |
| NVU cascade + MMC-1 (Section 10.5, 12) | `Phase10_NVU_causal_DAG.py`, `Phase13_NVU_energy_landscape.py`, `final_closure_mmc.py` | `PHASE10_*.md`, `PHASE12_NVU_INTEGRATED_MODEL.md`, `PHASE13_*.md` |
| Fire origin identification | `fire_origin_landscape.py`, `fire_origin_targeted.py` | `ATM_fire_origin_evidence_summary.md` |

## Quick start

### Python API

```python
from part3.ids_pipeline import IDSPipeline

# Auto-detect config and modules file
pipeline = IDSPipeline()

# Run full analysis
results = pipeline.run(
    expression_file='expression.csv',
    metadata_file='metadata.csv',
    output_dir='results/'
)
```

### Command-line

```bash
cd scrnaseq_hodge_pipeline/part3
python ids_pipeline/run_pipeline.py \
    /path/to/expression.csv \
    /path/to/metadata.csv \
    -o results/part3_analysis/ \
    --condition-column condition \
    --celltype-column cell_type \
    --pt-column PT_dpt
```

## Reproducing specific paper claims

### Phase 9 vascular origin (Section 10)

```bash
# Capillary endothelium + Pericyte PT_dpt gene trends
python scripts/Phase9a_CORRECTED_recalculate_PT_dpt_from_expression.py \
    --input <motor_cortex_expression> \
    --output-dir results/part3_phase9/

# Vascular driver state identification (ATM pathway 3-criteria)
python scripts/Phase9p_vascular_driver_state_identification.py

# Time-lag correlation across modules
python scripts/Phase9cd_vascular_timelag_correlation_analysis.py
```

Expected output: ATM pathway ATM p = 0.018, 8/9 gene concordance (89%),
satisfies all 3 criteria (temporal ordering, cross-context ≥ 80%, pathway
coherence ≥ 70%). See `reports/PHASE9_VASCULAR_ORIGIN_REPORT.md` for
detailed expected values.

### Phase 10 NVU causal DAG

```bash
python scripts/03a_prepare_notears_prior.py
python scripts/03b_run_notears.py
python scripts/03c_lingam_refinement.py
python scripts/Phase10_NVU_causal_DAG.py
```

### GSE212630 cross-dataset validation (Section 10, 83% concordance)

```bash
# PT_dpt recomputation on Wang 2023 dlPFC BA9 data
python scripts/GSE212630_compute_PT_dpt_and_validate.py \
    --input <GSE212630_raw_expression> \
    --output-dir results/part3_gse212630/

# Phi-flow analysis with PT_dpt
python scripts/GSE212630_phi_flow_with_PT_dpt.py

# Upstream module analysis for cross-dataset concordance
python scripts/GSE212630_upstream_module_analysis.py
```

Expected output: 83% direction concordance across 18 genes (ATM + R-loop
pathways). See `reports/GSE212630_ANALYSIS_SUMMARY.md` for exact values.

### Phase 13 NVU energy landscape

```bash
python scripts/Phase13_NVU_energy_landscape.py
python scripts/Phase13_complete_module_analysis.py
python scripts/Phase13c_NVU_energy_REFACTORED.py
```

## Data dependencies

- **NYGC ALS Consortium snRNA-seq** (same as Parts I–II)
  Primary motor cortex BA4; 43 cell-type annotations required
- **GSE212630** (Wang et al. 2023, *bioRxiv* 2023.01.12.523820)
  Dorsolateral prefrontal cortex BA9; 7 controls + 19 C9orf72 ALS/FTD donors
  Used for cross-dataset concordance analysis
- **`reports/24_functional_modules.json`** — canonical module definitions
  (Paper uses 23 of these; ALS_Genes excluded. See
  `../data/metadata/als_causative_genes.json` for that subset.)

## Python package dependencies

Beyond the core Part I pipeline requirements:

```
# LiNGAM causal inference
lingam>=1.8

# NOTEARS (if using the NOTEARS prior)
notears>=0.1  # or equivalent

# scanpy for diffusion pseudotime
scanpy>=1.10

# Standard scientific stack
scikit-learn>=1.3  # for axis discovery (PCA, NMF, KMeans)
networkx>=3.0      # for DAG manipulation
```

## Citations

- **Part III methodology**: Kaneko S., Urushitani M. (2026) — paper Sections 9–11
- **Cross-dataset data**: Wang, H.L.V. et al. (2023) *bioRxiv* 2023.01.12.523820
- **ids_core library**: extracted from Phase 13 + TDP-43/STMN2 analyses
  (2025-11-27 v1, 2025-11-30 v1.2 with Axis Discovery Engine)

## Important caveats (from paper Section 11)

1. **Limited confidence level**: Sample sizes for vascular cell types are
   hundreds of cells (881 Capillary, 511 Pericyte), one to two orders of
   magnitude smaller than major neuronal/glial populations.

2. **Cross-region + cross-disease comparison**: NYGC (sALS motor cortex BA4)
   ↔ GSE212630 (C9orf72 ALS/FTD dorsolateral prefrontal cortex BA9). The 83%
   concordance spans both region and disease boundaries. Generalisability
   to sALS motor cortex and other brain regions remains to be established.

3. **LiNGAM optimism under small sample**: LiNGAM estimates can be
   optimistic when sample sizes are small relative to variable
   dimensionality. The 23-module × hundreds-of-cells regime applies.

4. **mRNA ≠ protein**: ATM kinase activity is regulated by
   post-translational modifications invisible to transcriptomic analysis.

5. **Subtype composition confound**: If ALS alters vascular cell subtype
   composition, the observed ATM decline could reflect compositional
   change rather than per-cell downregulation.

## Source

Migrated from the author's private research directory:
`C:\Users\akaco\OneDrive\デスクトップ\projects\codex\docs\ids_causal_analysis\`

Total: **104 files** (52 scripts + 8 ids_pipeline + 24 reports + 18 docs +
1 methodology + this README).
