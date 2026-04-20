# Mendelian Randomization (MR) Analysis

Reproduces paper Section 4.5 / Methods 6.12: Two-sample MR using
oligodendrocyte cis-eQTLs (Bryois et al. 2022) as instruments for
exposure, ALS GWAS (van Rheenen et al. 2021) as outcome, across seven
gene set tracks (Tracks A–G).

## Data sources

- **Exposure**: Bryois et al. 2022, *Nat Neurosci* 25, 1104–1112
  [Zenodo 7276971](https://doi.org/10.5281/zenodo.7276971)
  Brain cell-type cis-eQTLs, n = 196 individuals, 8 cell types (GRCh37/38)

- **Outcome**: van Rheenen et al. 2021, *Nat Genet* 53, 1779–1788
  GWAS Catalog GCST90027164, N_eff = 80,713 (~28,000 cases)

- **LD reference**: 1000 Genomes Phase 3 EUR, n = 503 (GRCh37)

## Paper claims (Section 4.5, 8.7, Methods 6.12)

| Track | Gene Set | N IVs | IVW OR (95% CI) | p |
|-------|----------|-------|-----------------|----|
| A | All Oligo eQTL | 576 | 1.000 (0.998–1.002) | 0.766 |
| **B** | **Stable-High (ribosomal)** | **15** | **0.988 (0.972–1.005)** | **0.159** |
| D | Darkgrey (Rho/Ras GTPase) | 74 | 1.001 (0.995–1.008) | 0.660 |
| E | RTK/Trophic pathway | 8 | 0.997 (0.974–1.020) | 0.795 |
| **F** | **RQC (HBS1L)** | **1** | **0.957 (0.913–1.004)** | **0.072** |

- All four MR methods (IVW, MR-Egger, Weighted Median, Weighted Mode) show
  protective direction (OR < 1) in Tracks B and F.
- Egger intercept p = 0.791 (no directional pleiotropy)
- Cochran's Q p = 0.493 (no heterogeneity)
- **Power**: min detectable OR at 80% power = **1.012** → observed 0.988 at
  the boundary of detectability (power-limited null, not effect-absent null)

- **Colocalization (coloc.abf)**: KANSL1 at chr17q21.31 PP.H4 = **0.794**
  (only gene with moderate evidence of shared causal variant)

## Scripts

### R scripts (TwoSampleMR 0.7.0 + coloc)

| Script | Purpose |
|--------|---------|
| `00_inspect_data.R` | QC and preview of Bryois exposure data |
| `00b_inspect_als.R` | QC of van Rheenen outcome data |
| `02_run_mr.R` | Main MR execution (strict p < 5e-8) |
| `03_run_mr_clumped.R` | MR after LD clumping |
| `03b_run_mr_api_clump.R` | LD clumping via IEU OpenGWAS API |
| `04_relaxed_threshold_mr.R` | Relaxed p < 5e-6 threshold analysis |
| `05_smr_heidi.R` | SMR + HEIDI test |
| `06_final_figures.R` | Publication figures |
| `07_rqc_pathway_mr.R` | **Track F: RQC pathway MR (HBS1L)** |
| `08_mam_pathway_mr.R` | **Track G: MAM pathway MR** |
| `set_opengwas_token.R` | JWT token configuration |
| `setup_local_clumping.R` | Local PLINK + 1000G EUR setup |

## Reproducing

```bash
# Requires R 4.4.3+, TwoSampleMR 0.7.0, coloc, genetics.binaRies
Rscript benchmark/mendelian_randomization/00_inspect_data.R
Rscript benchmark/mendelian_randomization/00b_inspect_als.R
Rscript benchmark/mendelian_randomization/02_run_mr.R
Rscript benchmark/mendelian_randomization/03_run_mr_clumped.R
Rscript benchmark/mendelian_randomization/04_relaxed_threshold_mr.R   # Main Track B analysis
Rscript benchmark/mendelian_randomization/07_rqc_pathway_mr.R          # Main Track F analysis
Rscript benchmark/mendelian_randomization/08_mam_pathway_mr.R          # Main Track G analysis
```

## Key interpretations (Section 4.5, 8.7)

- **Track B (ribosomal 135 genes)**: All 4 MR methods show β < 0 (protective);
  Egger intercept / Cochran's Q confirm the signal is not pleiotropy or
  heterogeneity. This is the **cleanest MR signal** across all tracks.

- **Track F (HBS1L)**: Only 1 MR-testable gene (HBS1L) among 14 RQC-related
  genes because core RQC factors (LTN1, NEMF, TCF25, ZNF598) are housekeeping
  genes with tightly constrained expression and no significant cis-eQTLs.
  The Track B (ribosomal) + Track F (NGD) protective direction spans two
  biological layers (translation structure + stalled-ribosome clearance)
  with **zero gene overlap** — a convergence that strengthens the causal
  interpretation.

- **Colocalization KANSL1 (PP.H4 = 0.794)**: Moderate evidence of shared
  causal variant at chr17q21.31 inversion (top ALS GWAS locus). KANSL1 is
  NOT in the stable-High set, suggesting a chromatin/epigenetic pathway
  distinct from ribosomal disruption identified by IDS.

## Source

Migrated from:
- `sals_analysis_frozen_20260211/scripts/mr_validation/` (all R scripts)
- `sals_analysis_frozen_20260211/scripts/08_mam_pathway_mr.R`
- Full documentation: `sals_analysis_frozen_20260211/MR_VALIDATION_REPORT.md`
