
# GSE_212630 φ-Flow Analysis with Validated PT_pca
Generated: 2025-12-07 13:18:41

## Analysis Overview
- Pseudotime: PT_pca (validated against TDP-43 stages, ρ = 0.153, p < 1e-32)
- Total cells analyzed: 6,044
- PT bins: 10
- φ threshold for onset: 0.5

## Flow Summary by NVU Node (stress_total)
     group  onset_bin  onset_PT  peak_bin  peak_PT  peak_phi
  Vascular          1      0.15         7     0.75  6.184494
      Glia          0      0.05         8     0.85  4.657305
Excitatory          0      0.05         9     0.95  3.006091
Inhibitory          0      0.05         4     0.45  3.830967

## Key Findings

### 1. Temporal Ordering (Onset PT)
- Vascular: 0.1500
- Glia: 0.0500
- Excitatory: 0.0500

### 2. Peak φ by NVU Node
     group  peak_phi  peak_PT
  Vascular  6.184494     0.75
      Glia  4.657305     0.85
Excitatory  3.006091     0.95
Inhibitory  3.830967     0.45

### 3. Top Affected Cell Types (Late PT)
cell_type
Glia.Astro.GFAP.neg    11.557707
Vasc.Capillary          6.184494
Ex.L4                   6.150463
Ex.L5                   4.842111
Ex.L2/3                 3.782956

## Interpretation
- **Vascular onset**: PT = 0.1500 (earliest)
- **Glia onset**: PT = 0.0500
- **Excitatory onset**: PT = 0.0500

This pattern suggests: Vascular → Glia → Excitatory temporal ordering of stress accumulation.

## Files Generated
- phi_profiles_PT_pca_detailed.csv
- flow_summary_PT_pca_all_components.csv
- timelag_analysis_PT_pca.csv
- phi_profiles_by_celltype_PT_pca.csv
- Fig_phi_flow_PT_pca_main.png
- Fig_phi_flow_by_component.png
- Fig_phi_celltype_heatmap_PT_pca.png
