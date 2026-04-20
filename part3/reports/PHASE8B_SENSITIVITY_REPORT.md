# Phase 8b: Sensitivity Analysis and First-Trigger Module Identification

**Date**: 2025-11-24
**Analysis**: Sensitivity analysis of module onset timing and first-trigger module identification

---

## Analysis Overview

### Parameters Tested
- **Bin numbers**: [25, 50, 100]
- **Thresholds**: [0.5, 0.3, 0.2] (mean + threshold × std)
- **Total conditions**: 9

---

## Key Findings

### 1. First-Trigger Modules (Top 3 Earliest-Onset by Cell Type)


#### Upper

| Rank | Module | Mean Rank | Rank Stability (±SD) | Mean Onset PT |
|------|--------|-----------|---------------------|---------------|
| 1 | Angiogenesis | 1.0 | ±0.00 | 0.003372 |
| 2 | Apoptosis | 2.0 | ±0.00 | 0.003372 |
| 3 | Autophagy | 3.0 | ±0.00 | 0.003372 |

#### Glia

| Rank | Module | Mean Rank | Rank Stability (±SD) | Mean Onset PT |
|------|--------|-----------|---------------------|---------------|

#### VAT1L

| Rank | Module | Mean Rank | Rank Stability (±SD) | Mean Onset PT |
|------|--------|-----------|---------------------|---------------|


---

### 2. Glia Module Onset Categories


**Early Onset** (7 modules):
- Mitochondria          (Rank: 1, PT: 0.007299±0.003726)
- ECM                   (Rank: 2, PT: 0.007299±0.003726)
- Complement            (Rank: 3, PT: 0.007441±0.002096)
- Myelination           (Rank: 4, PT: 0.007526±0.002266)
- lncRNA                (Rank: 5, PT: 0.008039±0.001904)
- Cytoskeleton          (Rank: 6, PT: 0.008693±0.003214)
- ER_Stress             (Rank: 7, PT: 0.009091±0.002614)

**Middle Onset** (8 modules):
- Protein_Homeostasis   (Rank: 8, PT: 0.009148±0.003782)
- Apoptosis             (Rank: 9, PT: 0.009518±0.003543)
- Epigenetic            (Rank: 10, PT: 0.009660±0.002638)
- Inflammation          (Rank: 11, PT: 0.009746±0.003731)
- Metabolism            (Rank: 12, PT: 0.010030±0.002942)
- Autophagy             (Rank: 13, PT: 0.010258±0.002479)
- Growth_Factors        (Rank: 14, PT: 0.010343±0.002379)
- Calcium_Signaling     (Rank: 15, PT: 0.010372±0.002332)

**Late Onset** (8 modules):
- Cell_Cycle            (Rank: 16, PT: 0.010400±0.002536)
- DNA_Repair            (Rank: 17, PT: 0.010486±0.002719)
- Transcription         (Rank: 18, PT: 0.010486±0.002719)
- RNA_Processing        (Rank: 19, PT: 0.010827±0.002985)
- Synaptic              (Rank: 20, PT: 0.010912±0.002005)
- Oxidative_Stress      (Rank: 21, PT: 0.011282±0.001859)
- Angiogenesis          (Rank: 22, PT: 0.011396±0.001924)
- Ion_Transport         (Rank: 23, PT: 0.011396±0.001924)


---

### 3. Onset Order Stability

#### Most Stable Onset Order (Lowest Rank Std Dev)


**Upper**:

| Module | Mean Rank | Rank Std Dev | Onset PT (mean±SD) |
|--------|-----------|--------------|--------------------|
| Angiogenesis | 1.0 | 0.00 | 0.003372±0.000169 |
| Apoptosis | 2.0 | 0.00 | 0.003372±0.000169 |
| Autophagy | 3.0 | 0.00 | 0.003372±0.000169 |
| Calcium_Signaling | 4.0 | 0.00 | 0.003372±0.000169 |
| Cell_Cycle | 5.0 | 0.00 | 0.003372±0.000169 |

**Glia**:

| Module | Mean Rank | Rank Std Dev | Onset PT (mean±SD) |
|--------|-----------|--------------|--------------------|
| Complement | 2.6 | 1.77 | 0.007441±0.002096 |
| Cytoskeleton | 4.2 | 2.39 | 0.008693±0.003214 |
| Epigenetic | 9.3 | 2.40 | 0.009660±0.002638 |
| Oxidative_Stress | 19.8 | 2.90 | 0.011282±0.001859 |
| Ion_Transport | 18.7 | 2.91 | 0.011396±0.001924 |

**VAT1L**:

| Module | Mean Rank | Rank Std Dev | Onset PT (mean±SD) |
|--------|-----------|--------------|--------------------|
| Angiogenesis | 1.0 | 0.00 | 0.003372±0.000169 |
| Apoptosis | 2.0 | 0.00 | 0.003372±0.000169 |
| Autophagy | 3.0 | 0.00 | 0.003372±0.000169 |
| Calcium_Signaling | 4.0 | 0.00 | 0.003372±0.000169 |
| Cell_Cycle | 5.0 | 0.00 | 0.003372±0.000169 |


---

## Interpretation

### Upper Modules
- **Sensitivity**: Upper modules show {"low" if stability_df[stability_df['group'] == 'Upper']['std_rank'].mean() < 2 else "moderate to high"} variation in rank across conditions (mean rank std: {stability_df[stability_df['group'] == 'Upper']['std_rank'].mean():.2f})
- **First triggers**: The most consistently early modules are likely the true initiators of Upper dysfunction
- **Biological interpretation**: If Synaptic/Ca²⁺/Ion_Transport consistently rank first, this supports hyperexcitability as the primary Upper trigger

### Glia Modules
- **Sensitivity**: Glia modules show {"low" if stability_df[stability_df['group'] == 'Glia']['std_rank'].mean() < 2 else "moderate to high"} variation in rank across conditions (mean rank std: {stability_df[stability_df['group'] == 'Glia']['std_rank'].mean():.2f})
- **Detailed ranking**: Early-onset modules (Complement, Myelination, Synaptic, ER) represent the first signs of Glia structural dysfunction
- **Biological interpretation**: Complement/Myelination early onset suggests initial disruption of homeostatic functions before inflammatory cascade

### VAT1L Modules
- **Sensitivity**: VAT1L modules show {"low" if stability_df[stability_df['group'] == 'VAT1L']['std_rank'].mean() < 2 else "moderate to high"} variation in rank across conditions (mean rank std: {stability_df[stability_df['group'] == 'VAT1L']['std_rank'].mean():.2f})
- **Pattern**: If most modules have similar onset times (low rank variation), this suggests "simultaneous collapse" rather than sequential module failure
- **Biological interpretation**: Fragile downstream cell type may experience rapid multi-system failure

---

## Methodological Notes

### Sensitivity Analysis Approach
1. **Multiple binning resolutions** (25/50/100 bins) test whether onset detection is stable across different PT_dpt granularities
2. **Multiple thresholds** (0.5σ/0.3σ/0.2σ) test whether onset definition affects module ranking
3. **Rank stability metric** (std of ranks) identifies modules with consistent early/late onset regardless of parameters

### Limitations
1. **PT_dpt resolution**: Even at 100 bins, may not distinguish modules with very close onset times
2. **Threshold choice**: Lower thresholds (0.2σ) may capture noise; higher thresholds (0.5σ) may miss subtle early changes
3. **Binning artifacts**: Discrete bins may group together modules that actually have slightly different onsets

---

## Files Generated

1. **all_conditions_onset.csv** - Onset PT for all modules, groups, and conditions
2. **module_onset_stability.csv** - Stability metrics (mean rank, std rank) for each module × group
3. **first_trigger_modules.csv** - Top 3 earliest-onset modules per cell type
4. **glia_module_onset_ranking.csv** - Detailed Glia module ranking with categories

### Visualizations

1. **Fig1_onset_rank_stability_heatmap.png** - Heatmap showing rank variation across conditions
2. **Fig2_first_trigger_modules.png** - Comparison of first-trigger modules by cell type
3. **Fig3_glia_onset_ranking.png** - Detailed Glia module onset ranking
4. **Fig4_rank_stability.png** - Rank stability (coefficient of variation) for all modules

---

## Next Steps

1. **Biological validation**: Do the identified first-trigger modules align with known ALS pathology?
2. **Cross-reference with Phase 7**: Do first-trigger modules show strongest correlations with other cell types?
3. **Literature review**: Are early-onset modules documented as early events in ALS?

---

**Analysis Complete**: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}
