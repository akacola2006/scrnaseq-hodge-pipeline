# Phase 8: PT_dpt Time-Lag Analysis Report

**Date**: 2025-11-24
**Analysis**: Module onset/peak timing along PT_dpt pseudotime
**Dataset**: 27,170 cells across 3 groups (Upper, Glia, VAT1L)

---

## Executive Summary

**Goal**: Determine which cell type shows module changes FIRST along PT_dpt.

**Key Question**: Do Upper modules precede Glia/VAT1L changes (supporting Upper-first hypothesis),
or does Glia lead (supporting Glia-first), or VAT1L?

**Method**:
1. Binned cells by PT_dpt (25 bins from 0.000000 to 0.025609)
2. Calculated module averages for each group × bin
3. Identified "onset" PT (first bin where module > mean + 0.5*std)
4. Calculated time-lags: ΔPT (Upper onset - Glia/VAT1L onset)

**Interpretation**:
- **Negative ΔPT**: Upper module onset precedes Glia/VAT1L (Upper-first)
- **Positive ΔPT**: Glia/VAT1L module onset precedes Upper (Glia/VAT1L-first)

---

## Overall Statistics

### Upper vs Glia Timing

**Onset time-lags**:
- Mean ΔPT: -0.008863
- Median ΔPT: -0.009219

**Module counts**:
- **Upper-first** (Δ<0): 23 modules (100.0%)
- **Glia-first** (Δ>0): 0 modules (0.0%)

---

### Upper vs VAT1L Timing

**Onset time-lags**:
- Mean ΔPT: 0.000000
- Median ΔPT: 0.000000

**Module counts**:
- **Upper-first** (Δ<0): 0 modules (0.0%)
- **VAT1L-first** (Δ>0): 0 modules (0.0%)

---

## Key Module Time-Lags

### Hyperexcitability-Related Modules

### Hyperexcitability_related

| Module | Upper Onset | Glia Onset | VAT1L Onset | Δ(U-G) | Δ(U-V) | Pattern |
|--------|-------------|------------|-------------|--------|--------|----------|
| Synaptic | 0.00359 | 0.01178 | 0.00359 | -0.00820 | +0.00000 | U→G, V→U |
| Calcium_Signaling | 0.00359 | 0.01178 | 0.00359 | -0.00820 | +0.00000 | U→G, V→U |
| Ion_Transport | 0.00359 | 0.01280 | 0.00359 | -0.00922 | +0.00000 | U→G, V→U |

### Stress_related

| Module | Upper Onset | Glia Onset | VAT1L Onset | Δ(U-G) | Δ(U-V) | Pattern |
|--------|-------------|------------|-------------|--------|--------|----------|
| ER_Stress | 0.00359 | 0.01178 | 0.00359 | -0.00820 | +0.00000 | U→G, V→U |
| Oxidative_Stress | 0.00359 | 0.01280 | 0.00359 | -0.00922 | +0.00000 | U→G, V→U |
| Protein_Homeostasis | 0.00359 | 0.01383 | 0.00359 | -0.01024 | +0.00000 | U→G, V→U |

### Inflammatory

| Module | Upper Onset | Glia Onset | VAT1L Onset | Δ(U-G) | Δ(U-V) | Pattern |
|--------|-------------|------------|-------------|--------|--------|----------|
| Inflammation | 0.00359 | 0.01383 | 0.00359 | -0.01024 | +0.00000 | U→G, V→U |
| Complement | 0.00359 | 0.00973 | 0.00359 | -0.00615 | +0.00000 | U→G, V→U |

### Death_related

| Module | Upper Onset | Glia Onset | VAT1L Onset | Δ(U-G) | Δ(U-V) | Pattern |
|--------|-------------|------------|-------------|--------|--------|----------|
| Apoptosis | 0.00359 | 0.01383 | 0.00359 | -0.01024 | +0.00000 | U→G, V→U |
| Autophagy | 0.00359 | 0.01178 | 0.00359 | -0.00820 | +0.00000 | U→G, V→U |

---

## Interpretation

### Patterns Observed

1. **Upper-first modules**: Modules where Upper onset precedes Glia/VAT1L
   - These support "Upper functional driver" hypothesis
   - Suggest Upper changes temporally precede Glia/VAT1L changes

2. **Glia-first modules**: Modules where Glia onset precedes Upper/VAT1L
   - These support "Glia metabolic upstream" hypothesis
   - Suggest Glia dysfunction precedes neuronal changes

3. **Mixed patterns**: Some modules show Upper→Glia but Glia→VAT1L
   - Suggest multi-step cascade

### Biological Interpretation

**If Hyperexcitability (Synaptic, Ca²⁺, Ion) shows Upper-first**:
→ Supports Upper neurons as early functional drivers
→ Excitotoxic mechanisms may propagate to Glia/VAT1L

**If ER Stress/Inflammation shows Glia-first**:
→ Supports Glia metabolic dysfunction as initiating event
→ Glia stress propagates to neurons

**If Mixed patterns**:
→ Supports dual-upstream model
→ Both Glia (metabolic) and Upper (functional) contribute

---

## Methodological Considerations

### Strengths

1. **PT_dpt is stress-independent**: R²=1.9% with stress_total
   - Avoids circular logic (unlike PT_imes with R²=62.5%)
2. **Quantitative time-lag measurement**: ΔPT provides numerical comparison
3. **Module-specific analysis**: Tests pathway-specific timing, not just aggregate stress

### Limitations

1. **Pseudotime ≠ real time**:
   - PT_dpt represents trajectory, not absolute time
   - "Onset" timing is relative, not absolute
   - Cannot distinguish fast vs slow processes

2. **Cross-sectional data**:
   - Infers dynamics from snapshot
   - True temporal precedence requires longitudinal data

3. **Binning artifacts**:
   - Onset detection depends on bin size and threshold
   - Sensitivity analysis needed (vary threshold 0.3-1.0σ)

4. **Cell density variation**:
   - Some bins have few cells (low n)
   - May introduce noise in low-density regions

5. **Interpretation caveats**:
   - "Upper-first" onset ≠ "Upper causes Glia changes"
   - Requires integration with Phase 7 correlational evidence
   - Spatial co-localization still needed for causality

---

## Integration with Phase 7 Findings

### Convergent Evidence for Upper Functional Driver

**Phase 7 (Correlational)**:
- Upper ER → Glia ER: r=0.664, p=0.010
- Upper Synaptic → Glia Inflammation: r=0.534, p=0.049
- Upper heterogeneity → VAT1L stress: r=0.569, p=0.034

**Phase 8 (Temporal)**:
- If ER Stress shows Upper-first onset → temporal precedence
- If Synaptic shows Upper-first onset → temporal precedence
- Combines with Phase 7 correlations → stronger driver hypothesis support

### Complementarity

- **Phase 7**: Tests "if Upper high, then Glia/VAT1L high" (correlation)
- **Phase 8**: Tests "Upper changes before Glia/VAT1L" (temporal sequence)
- **Together**: Both correlation AND temporal precedence → stronger causal inference

---

## Conclusions

### Main Findings

1. **Onset timing statistics**:
   - Upper vs Glia: {onset_stats['Upper_first_vs_Glia']}/{onset_stats['Total_modules']} modules show Upper-first
   - Upper vs VAT1L: {onset_stats['Upper_first_vs_VAT1L']}/{onset_stats['Total_modules']} modules show Upper-first

2. **Key module patterns**: (See tables above)

3. **Interpretation**:
   - Pattern suggests [TO BE DETERMINED after reviewing module-specific results]
   - Integration with Phase 7 correlational evidence needed

### Next Steps

1. **Sensitivity analysis**:
   - Test different onset thresholds (0.3σ, 0.5σ, 1.0σ)
   - Test different bin numbers (10, 20, 30)
   - Assess robustness of findings

2. **Module grouping**:
   - Group modules by functional category
   - Test if categories show consistent patterns

3. **Integration with correlations**:
   - Modules with both Upper-first onset AND strong Phase 7 correlation
   - Highest confidence driver candidates

4. **Validation**:
   - Longitudinal data (if available)
   - Spatial transcriptomics (test local co-progression)
   - Perturbation experiments (causal validation)

---

**Generated**: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}
**Scripts**: `scripts/Phase8_PT_dpt_time_lag_analysis.py`
**Data**: `results/phase8_time_lag/`

**Key Files**:
- `module_profiles_by_group_and_bin.csv`: Group × bin × module averages
- `module_onset_peak_PT_by_group.csv`: Onset/peak PT for each module
- `module_time_lags.csv`: ΔPT statistics
- `Fig_module_onset_time_lags.png`: Time-lag visualization
- `Fig_module_trajectories.png`: Module expression trajectories
