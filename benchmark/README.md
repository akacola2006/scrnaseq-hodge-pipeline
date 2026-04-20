# External Benchmarks

Reproduction scripts for the four external validation benchmarks reported
in the paper: Norman Perturb-seq (Section 2.4), Shen CRISPR (Section 2.5),
TimeVault (Section 2.3), Glioma (Section 2.6), and Mendelian Randomization
(Section 4.5 / 8.7).

Each subdirectory contains a `README.md` documenting the dataset, paper
claims, reproduction commands, and source references.

## Directory structure

```
benchmark/
├── README.md                       ← this file
├── norman_perturb_seq/             ← Section 2.4 (IDS Top-1 = 66.4%)
│   ├── README.md
│   ├── run_ids_norman.py           ← R2b (whitening + bootstrap)
│   └── run_grnboost2_norman.py     ← GRNBoost2 baseline
├── shen_crispr/                    ← Section 2.5 (mean 21.8 %ile, p = 0.002)
│   ├── README.md
│   ├── perturbation_loader.py
│   ├── ortholog_map.py
│   ├── run_shen_ids.py             ← Track 2: per-perturbation Gene Hodge
│   ├── run_shen_ids_sparse.py      ← k-NN variant (k = 30: p = 0.020)
│   └── track2_gene_hodge.py
├── timevault/                      ← Section 2.3 (TAS AUC = 0.9475)
│   ├── README.md
│   └── run_tas.py                  ← Brake-loss pipeline, full TAS evaluation
├── glioma/                         ← Section 2.6 (Null GF = 0.4248)
│   ├── README.md
│   ├── step8_pseudotime.py
│   ├── step9_gene_insertion.py
│   ├── run_glioma_pseudotime_insertion.py
│   └── run_idh_stratified_insertion.py
└── mendelian_randomization/        ← Section 4.5 / 8.7 (Track B OR = 0.988)
    ├── README.md
    ├── 00_inspect_data.R           ← Bryois exposure data QC
    ├── 00b_inspect_als.R           ← van Rheenen outcome QC
    ├── 02_run_mr.R                 ← Main MR
    ├── 03_run_mr_clumped.R         ← LD clumping
    ├── 03b_run_mr_api_clump.R
    ├── 04_relaxed_threshold_mr.R   ← Track B (p < 5e-6)
    ├── 05_smr_heidi.R
    ├── 06_final_figures.R
    ├── 07_rqc_pathway_mr.R         ← Track F (HBS1L)
    ├── 08_mam_pathway_mr.R         ← Track G (MAM pathway)
    ├── set_opengwas_token.R
    └── setup_local_clumping.R
```

## Summary table (paper benchmarks → scripts)

| Paper Section | Benchmark           | Dataset                | Key number         | Script                                           |
|---------------|---------------------|------------------------|--------------------|--------------------------------------------------|
| 2.3           | TimeVault TAS       | Chao 2026 (PC9)        | AUC = 0.9475       | `timevault/run_tas.py`                           |
| 2.4           | Norman Perturb-seq  | Norman 2019 (K562)     | Top-1 = 66.4 %     | `norman_perturb_seq/run_ids_norman.py`           |
| 2.5           | Shen CRISPR         | Shen 2026 (hippocampus)| p = 0.002          | `shen_crispr/run_shen_ids.py`                    |
| 2.6           | Glioma              | TCGA + GTEx            | Null GF = 0.4248   | `glioma/run_glioma_pseudotime_insertion.py`      |
| 4.5, 8.7      | MR                  | Bryois × van Rheenen   | OR = 0.988, p=0.159| `mendelian_randomization/02_run_mr.R`            |

## Note on internal consistency

These scripts are migrated from the original research directories
(`cell navi/`, `sals_analysis_frozen_20260211/`, `TimeVault/`). They rely on
the pre-paper analysis directory structure and will require minor path edits
for your environment. Each `README.md` provides a step-by-step reproduction
guide with expected numerical outputs matching the paper.

For the core Hodge pipeline applied to sALS, see `scripts/` (Part I) and
the 3φ residual framework scripts (`scripts/three_phi_residual.py`,
`scripts/allgene_insertion.py`, `scripts/nemf_screen.py`,
`scripts/de_sals_only.py`).
