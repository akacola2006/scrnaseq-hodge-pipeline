# Original vs IDS-PT Methodology: Corrected Analysis
**Created**: 2025-11-23
**Status**: Corrected interpretation based on actual implementation

---

## Executive Summary

After careful examination of the actual implementation code, both Original pairwise and IDS-PT state-level analyses measure **disease progression timing** (progressive causality), not regulatory control. However, they differ fundamentally in:
1. **What they compare**: Modules vs whole cells
2. **PT space**: Per-cell-type vs global
3. **Decision mechanism**: O1/O2/O3 majority vote vs pure PT-ordering
4. **Resolution**: Cell-type level vs cell-state level

The low correlation (r=-0.041) reflects these fundamental methodological differences, not a failure of either method.

---

## 1. Original Pairwise Analysis (O1-O3)

### Implementation Review

#### O1: Lag Asymmetry (pairwise_o1o3_analysis.py:51-140)
```python
# For each module m shared by cell types A and B:
for lag in range(-lagwin, lagwin+1):
    if lag > 0:  # A leads B
        corr_forward.append(spearmanr(I_A[lag:], I_B[:]))
    else:  # B leads A
        corr_backward.append(spearmanr(I_B[lag:], I_A[:]))

stat_obs = mean(corr_forward) - mean(corr_backward)
# If stat_obs > 0 → A→B
```

**What it measures**: Temporal precedence in module activity correlation
**Interpretation**: If module m in cell A correlates with future values of module m in cell B, then A→B
**Causality type**: Temporal correlation (NOT causal, can be confounded by common drivers)
**Signal**: Weak - only 0.08 modules/pair on average show significance

---

#### O2: Conditional Asymmetry (pairwise_o1o3_analysis.py:142-196)
```python
# For each module m:
# Calculate Δφ variance in ALS vs Control
drop_A = var(Δφ_A^Control) - var(Δφ_A^ALS)
drop_B = var(Δφ_B^Control) - var(Δφ_B^ALS)

# If drop_A > drop_B → A→B
```

**What it measures**: Differential loss of homeostasis
**Interpretation**: If cell type A loses Δφ variance more than B (Control→ALS), then A→B
**Causality type**: Differential vulnerability/response (NOT directional causality)
**Signal**: Strong - 20.78 modules/pair on average show significance
**Problem**: This measures "which fails first" not "which controls which"

---

#### O3: Energy Flow (pairwise_o1o3_analysis.py:198-252)
```python
# For each module m:
slope_A = linregress(PT_A, Δφ_A).slope
slope_B = linregress(PT_B, Δφ_B).slope

# If slope_A < slope_B → A→B
```

**What it measures**: Rate of stress accumulation along pseudotime
**Interpretation**:
- If slope_A < slope_B (A's stress increases slower than B's):
  - A reaches high stress earlier in PT (low PT, high Δφ)
  - B reaches high stress later in PT (high PT, high Δφ)
  - Therefore A→B in disease progression

**Causality type**: **Progressive causality** (disease trajectory timing)
**Signal**: Strong - 16.36 modules/pair on average
**This is valid trajectory timing**, but uses per-cell-type module PT, not global cell PT

---

#### Majority Vote Decision
```python
decisions = [O1_decision, O2_decision, O3_decision]
final_decision = most_common(decisions)
```

**Problem**: Combines three different concepts:
1. O1: Temporal correlation (weak signal, 0.08 modules/pair)
2. O2: Differential vulnerability (strong signal, 20.78 modules/pair)
3. O3: Progressive timing (strong signal, 16.36 modules/pair)

**Result bias**: O2 and O3 dominate (both ~90% sensitivity), so final decision primarily reflects "which cell type is more vulnerable" mixed with "which degenerates earlier".

---

## 2. IDS-PT State-Level Analysis

### Implementation Review (scripts/05c_calculate_state_directions.py)

```python
# Build kNN graph on feature space
features = [PT_IDS, stress_total, module_1, ..., module_10]
nbrs = NearestNeighbors(n_neighbors=15)
nbrs.fit(features)

# For each cell i and its neighbors j:
for j in neighbors(i):
    if PT_j > PT_i:
        # j is downstream of i
        direction_count[state_i][state_j] += 1
        pt_weight[state_i][state_j] += (PT_j - PT_i)

# Upstream score
upstream_score = (out_count - in_count) + (out_weight - in_weight)
```

**What it measures**: Disease progression ordering among similar cells
**Interpretation**:
- Cells similar in feature space (stress + modules) are neighbors
- Among neighbors, cells with lower PT are upstream, higher PT are downstream
- Aggregated to state level, then to cell-type level

**Causality type**: **Progressive causality** (disease trajectory timing)
**PT space**: **Global cell-level PT_IDS** (DPT across all cells)
**Resolution**: Cell-state level (170 states from 42 types)

---

## 3. Why Are They Measuring Different Things?

### Key Differences

| Aspect | Original O1-O3 | IDS-PT State-Level |
|--------|----------------|-------------------|
| **Unit of comparison** | Module activity between cell types | Whole cells in feature space |
| **PT space** | Per-cell-type module PT | Global cell-level PT_IDS |
| **What defines "neighbors"** | All pairs compared | kNN (k=15) in feature space |
| **Direction signal** | O1(correlation) + O2(vulnerability) + O3(timing) | Pure PT ordering |
| **Resolution** | Cell-type level (43 types) | Cell-state level (170 states → 42 types) |
| **Aggregation** | Majority vote across modules | Sum across states |

---

### Mechanistic Interpretation

#### Original O3 (Progressive Component Only)
For module m in cell types A and B:
1. Calculate slope of PT vs Δφ for each cell type independently
2. If slope_A < slope_B:
   - A reaches high stress earlier in its own PT trajectory
   - B reaches high stress later in its own PT trajectory
   - **Conclusion**: A→B (A degenerates before B in disease progression)

**PT space**: Each cell type has its own module-level PT
**Comparison**: Between-cell-type, within-module

---

#### IDS-PT State-Level
For cells i and j in a shared kNN neighborhood:
1. Both cells are similar in [PT, stress, modules]
2. If PT_j > PT_i:
   - Cell j appears later in global disease trajectory
   - Cell i appears earlier in global disease trajectory
   - **Conclusion**: i→j (i's state precedes j's state)

**PT space**: Global cell-level PT_IDS (DPT across all 111,837 cells)
**Comparison**: Between-cell, within-neighborhood

---

## 4. Why Correlation is Low (r=-0.041)

### Hypothesis 1: Different PT Spaces
- **Original**: Uses per-cell-type PT for each module independently
  - Ex_L5_VAT1L has its own PT for each module
  - Glia_Oligo has its own PT for each module
  - Comparison: "Does VAT1L's Angiogenesis module reach stress earlier in VAT1L's PT than Oligo's Angiogenesis module reaches stress in Oligo's PT?"

- **IDS-PT**: Uses global PT_IDS across all cells
  - All cells exist in the same PT space
  - Comparison: "Do VAT1L cells appear earlier or later than Oligo cells in the global disease trajectory?"

**These are fundamentally different questions.**

---

### Hypothesis 2: O2 Contamination
- Original's majority vote includes O2 (differential vulnerability)
- O2 has highest sensitivity (20.78 modules/pair) and dominates decisions
- O2 measures "which loses homeostasis more" not "which appears earlier in PT"
- This shifts Original's results away from pure progressive timing

**Evidence**:
- O2 alone: 20.78 significant modules/pair
- O3 alone: 16.36 significant modules/pair
- O1 alone: 0.08 significant modules/pair
- Final decision heavily weighted by O2

---

### Hypothesis 3: Aggregation Artifacts
- IDS-PT finds biphasic patterns (e.g., Glia_Oligo):
  - State4: +15,356 (strongest upstream)
  - State0: -16,346 (strongest downstream)
  - **Aggregated**: +3,242 (appears upstream)

- Original sees single directionality:
  - Glia_Oligo: -15 (downstream)

**The aggregation masks state-level heterogeneity.**

---

### Hypothesis 4: Module vs Cell-Level Features
- **Original**: Compares individual modules between cell types
  - "Is Angiogenesis in A earlier than Angiogenesis in B?"
  - 23 modules × 43 cell types = 989 module×cell combinations
  - Each module analyzed independently

- **IDS-PT**: Uses aggregated module features
  - 10 top-variance modules combined with PT and stress
  - kNN groups cells with similar overall profiles
  - Direction based on global PT, not per-module timing

**Different feature spaces lead to different neighborhoods and directions.**

---

## 5. Which Method is "Correct"?

### Both are valid for different questions:

#### Original O3 (isolating the progressive component):
**Question**: "For a given module, which cell type reaches stress earlier in its own disease trajectory?"
**Use case**: Identifying cell-type-specific module vulnerability patterns
**Limitation**:
- Uses per-cell-type PT (cannot compare timing across types directly)
- O2 contamination in majority vote
- No state-level resolution

---

#### IDS-PT State-Level:
**Question**: "In the global disease trajectory, which cell states appear earlier and which appear later?"
**Use case**: Mapping the overall progression across all cell types and states
**Limitation**:
- Aggregation to cell-type level masks biphasic patterns
- kNN threshold (k=15) affects neighborhood definition
- Does not capture per-module timing differences

---

## 6. Reconciling the Two Approaches

### Case Study: VAT1L Motor Neurons

#### Original Pairwise Result:
- Ex_L5_VAT1L_EYA4: net_directionality = -1 (slightly downstream)
- Ex_L5_VAT1L_THSD4: net_directionality = -5 (downstream)
- **Interpretation**: VAT1L is a victim, not initiator (O2+O3 consensus)

#### IDS-PT Result:
- Ex_L5_VAT1L_EYA4 aggregated: upstream_score = -1,793 (strong downstream)
- Ex_L5_VAT1L_THSD4 aggregated: upstream_score = -576 (downstream)
- **State-level breakdown**:
  - EYA4_State0: -8,442 (extreme sink)
  - EYA4_State2: -319 (moderate sink)
  - THSD4_State1: +1,021 (moderate source)
  - THSD4_State2: -1,597 (sink)

**Consensus**: Both agree VAT1L is downstream/victim
**Difference**: IDS-PT shows stronger downstream signal and reveals state heterogeneity

---

### Case Study: Oligodendrocytes (Biphasic)

#### Original Pairwise Result:
- Glia_Oligo: net_directionality = -15 (downstream)
- **Interpretation**: Victim of myelination stress (O2+O3)

#### IDS-PT Result:
- Glia_Oligo aggregated: upstream_score = +3,242 (upstream)
- **State-level breakdown**:
  - State4: +15,356 (strongest global source)
  - State0: -16,346 (strongest global sink)
  - States 1,2,3,5: Mixed

**Interpretation**:
- Oligo has **both** upstream (early, State4) and downstream (late, State0) roles
- Original captures the dominant downstream signal (averaged across states)
- IDS-PT reveals biphasic pattern: oligodendrocytes both initiate and suffer

**This is NOT a contradiction** - it's complementary information at different resolutions.

---

## 7. Corrected Conclusions

### What We Learned:

1. **Both methods measure progressive causality** (disease trajectory timing), not regulatory control
   - Original O3: Per-cell-type module PT slope
   - IDS-PT: Global cell-level PT ordering

2. **Original's low O1 sensitivity suggests weak temporal correlation** between cell types
   - Module activities don't show strong lag correlations
   - Temporal precedence signal is weak

3. **Original's high O2 sensitivity suggests strong differential vulnerability**
   - Cell types show different rates of homeostasis loss
   - This dominates the majority vote, shifting results toward vulnerability rather than timing

4. **IDS-PT reveals state-level heterogeneity** that aggregated cell-type analysis misses
   - Biphasic patterns (upstream early states, downstream late states)
   - Important for understanding disease dynamics within cell types

5. **Low correlation (r=-0.041) is explained by**:
   - Different PT spaces (per-type vs global)
   - O2 contamination in Original
   - Aggregation artifacts in IDS-PT
   - Module-level vs cell-level comparisons

---

## 8. Recommendations for Future Analysis

### To test these hypotheses:

1. **Re-run Original pairwise using O3 only** (exclude O1/O2)
   - This isolates progressive timing signal
   - Compare with IDS-PT: should correlate better

2. **Re-run IDS-PT without state aggregation**
   - Report state-level upstream scores separately
   - Don't aggregate to cell-type level

3. **Build global module-level PT**
   - Use DPT on concatenated module data across all cell types
   - Compare module timing in shared PT space
   - This would match IDS-PT's global approach

4. **Analyze biphasic patterns systematically**
   - Identify all cell types with both upstream and downstream states
   - Investigate transition mechanisms
   - This is unique insight from state-level resolution

---

## 9. Final Answer to the User's Question

**Question**: "ほんとうにオリジナルのpairwise解析は制御の方向性をみているんですか？"
**Translation**: "Does the original pairwise analysis really look at regulatory directionality?"

**Answer**: **No, it does not measure regulatory causality.**

- **O1** measures temporal correlation precedence (weak signal, likely confounded)
- **O2** measures differential vulnerability/homeostasis loss (strong signal, but not directional causality)
- **O3** measures disease progression timing (strong signal, **this IS progressive causality**)

The majority vote combines these three, with O2 (vulnerability) dominating due to highest sensitivity. The final "direction" primarily indicates **which cell type is more vulnerable** mixed with **which degenerates earlier in its own PT trajectory**.

**IDS-PT measures the same progressive causality as O3**, but using global cell-level PT instead of per-cell-type module PT. The low correlation reflects fundamental methodological differences, not measurement of different types of causality (my original interpretation was incorrect).

**Both are valid disease trajectory analyses, not regulatory control analyses.**

---

**End of Corrected Analysis**
