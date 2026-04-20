
# GSE_212630 IDS Core Analysis Summary
Generated: 2025-12-07 12:43:36

## Dataset Overview
- Total cells analyzed: 9,766
- Cell types: 17
- Conditions: Control, TDPneg, TDPmed, TDPhigh

## Cell Counts by Condition
condition
TDPmed     3331
Control    2594
TDPneg     2507
TDPhigh    1334

## φ Summary by NVU Node
     group  bin  phi_mean  phi_std  n_cells
  Vascular    0  2.043567 2.475044       48
  Vascular    1  1.610104 2.266267      173
  Vascular    2  1.605161 2.211282      119
  Vascular    3  1.875748 2.296836      142
      Glia    0  0.974230 1.460082     1000
      Glia    1  1.194875 1.605576     1000
      Glia    2  0.941208 1.359135     1000
      Glia    3  1.002052 1.520844     1000
Excitatory    0  0.934419 1.261351      734
Excitatory    1  1.475666 2.007148      773
Excitatory    2  1.059133 1.404018     1284
Excitatory    3  1.797498 2.508792       70
Inhibitory    0  1.028097 1.446304      812
Inhibitory    1  1.088951 1.535263      561
Inhibitory    2  1.155861 1.616003      928
Inhibitory    3  0.834793 1.255534      122

## Flow Summary (Onset and Peak)
     group  onset_bin  onset_PT  peak_bin  peak_PT  peak_phi
  Vascular          0     0.125         0    0.125  2.043567
      Glia          0     0.125         1    0.375  1.194875
Excitatory          0     0.125         3    0.875  1.797498
Inhibitory          0     0.125         2    0.625  1.155861

## Top Affected Cell Types (TDPhigh)
cell_type
Ex.L5                  3.691016
Ex.Unknown             3.190997
Vasc.Endo              2.605419
Ex.L6                  2.219988
Vasc.Fibro             1.960559
Vasc.Pericyte          1.738652
Ex.L2/3                1.700312
Glia.Micro             1.308148
In.Unknown             0.965274
Glia.Astro.GFAP.pos    0.942935

## Files Generated
- combined_metadata_with_stress.csv
- phi_profiles_by_nvu_and_pt.csv
- flow_summary_by_nvu.csv
- celltype_condition_summary.csv
- Fig1_phi_trajectories_by_nvu.png
- Fig2_phi_by_stress_component.png
- Fig3_phi_by_celltype.png
