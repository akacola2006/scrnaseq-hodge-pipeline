#!/usr/bin/env python3
"""
Phase 6 Extended: Functional Driver Hypothesis Validation

Purpose:
- Test if Early Upper subgroup is a functional driver (causal) or just dying cells (victim)
- Patient-level correlation: Early Upper ratio vs Glia/VAT1L stress
- Cross-cell module coupling: Upper hyperexcitability vs VAT1L/Glia stress
- Mid-bin trajectory: Early → Mid → Late continuity
"""

import pandas as pd
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

print("="*80)
print("Phase 6 Extended: Functional Driver Hypothesis Validation")
print("="*80)

# Paths
base_dir = Path('/home/akaco/als/motor_cortex_analysis')
ids_dir = base_dir / 'ids_causal_analysis'
output_dir = ids_dir / 'results' / 'phase6_driver_validation'
output_dir.mkdir(parents=True, exist_ok=True)

# ============================================================================
# Step 1: Load Data
# ============================================================================
print("\n[Step 1] Loading data...")

data_file = ids_dir / 'results' / 'PTv2_robustness' / 'PTv2_quick_key_celltypes_with_PTs.csv'
df = pd.read_csv(data_file)

print(f"  Loaded: {data_file.name}")
print(f"  Total cells: {len(df):,}")

# Filter to ALS only
df_als = df[df['condition'] == 'ALS'].copy()
print(f"  ALS cells: {len(df_als):,}")

# Define PT_dpt groups
df_als['PT_dpt_group'] = pd.cut(
    df_als['PT_dpt'],
    bins=[0, 0.005, 0.010, 1.0],
    labels=['Early', 'Mid', 'Late'],
    include_lowest=True
)

# ============================================================================
# Step 2: Patient-Level Correlation Analysis
# ============================================================================
print("\n[Step 2] Patient-level correlation: Early Upper ratio vs Glia/VAT1L stress...")

# Define cell type groups
upper_types = ['Ex.L2.L3.CUX2.RASGRF2']
vat1l_types = ['Ex.L5.VAT1L.EYA4', 'Ex.L5.VAT1L.THSD4']
glia_types = ['Glia.Oligo', 'Glia.Astro.GFAP-neg', 'Glia.Micro']

# Calculate patient-level metrics
patient_metrics = []

for patient_id in df_als['patient_id'].unique():
    patient_cells = df_als[df_als['patient_id'] == patient_id]

    # Upper metrics
    upper_cells = patient_cells[patient_cells['cell_type'].isin(upper_types)]
    if len(upper_cells) >= 10:  # Minimum threshold
        upper_early = upper_cells[upper_cells['PT_dpt_group'] == 'Early']
        early_upper_ratio = len(upper_early) / len(upper_cells) * 100

        # Hyperexcitability score for Upper
        hyperexc_modules = ['module_Synaptic', 'module_Calcium_Signaling', 'module_Ion_Transport']
        available_hyperexc = [m for m in hyperexc_modules if m in upper_cells.columns]
        if len(available_hyperexc) > 0:
            upper_hyperexc = upper_cells[available_hyperexc].mean(axis=1).mean()
        else:
            upper_hyperexc = np.nan
    else:
        early_upper_ratio = np.nan
        upper_hyperexc = np.nan

    # VAT1L metrics
    vat1l_cells = patient_cells[patient_cells['cell_type'].isin(vat1l_types)]
    if len(vat1l_cells) >= 5:
        vat1l_stress = vat1l_cells['stress_total'].mean()

        # VAT1L ER stress
        if 'module_ER_Stress' in vat1l_cells.columns:
            vat1l_er_stress = vat1l_cells['module_ER_Stress'].mean()
        else:
            vat1l_er_stress = np.nan
    else:
        vat1l_stress = np.nan
        vat1l_er_stress = np.nan

    # Glia metrics
    glia_cells = patient_cells[patient_cells['cell_type'].isin(glia_types)]
    if len(glia_cells) >= 10:
        glia_stress = glia_cells['stress_total'].mean()

        # Glia Inflammation
        if 'module_Inflammation' in glia_cells.columns:
            glia_inflammation = glia_cells['module_Inflammation'].mean()
        else:
            glia_inflammation = np.nan
    else:
        glia_stress = np.nan
        glia_inflammation = np.nan

    patient_metrics.append({
        'patient_id': patient_id,
        'n_cells': len(patient_cells),
        'n_upper': len(upper_cells),
        'n_vat1l': len(vat1l_cells),
        'n_glia': len(glia_cells),
        'early_upper_ratio': early_upper_ratio,
        'upper_hyperexc': upper_hyperexc,
        'vat1l_stress': vat1l_stress,
        'vat1l_er_stress': vat1l_er_stress,
        'glia_stress': glia_stress,
        'glia_inflammation': glia_inflammation
    })

patient_df = pd.DataFrame(patient_metrics)
patient_df = patient_df.dropna(subset=['early_upper_ratio', 'vat1l_stress', 'glia_stress'])

print(f"\n  Patients with complete data: {len(patient_df)}")
print(f"  Early Upper ratio range: {patient_df['early_upper_ratio'].min():.1f}% - {patient_df['early_upper_ratio'].max():.1f}%")

# Save
patient_df.to_csv(output_dir / 'patient_level_metrics.csv', index=False)
print(f"  Saved: patient_level_metrics.csv")

# Correlations
print(f"\n  Patient-level correlations:")

# Early Upper ratio vs VAT1L stress
if len(patient_df) >= 3:
    r_vat, p_vat = stats.pearsonr(patient_df['early_upper_ratio'], patient_df['vat1l_stress'])
    print(f"    Early Upper ratio vs VAT1L stress:  r={r_vat:+.3f}, p={p_vat:.4f}")

    # Early Upper ratio vs Glia stress
    r_glia, p_glia = stats.pearsonr(patient_df['early_upper_ratio'], patient_df['glia_stress'])
    print(f"    Early Upper ratio vs Glia stress:   r={r_glia:+.3f}, p={p_glia:.4f}")

    # Upper hyperexcitability vs VAT1L stress
    patient_df_hyperexc = patient_df.dropna(subset=['upper_hyperexc', 'vat1l_stress'])
    if len(patient_df_hyperexc) >= 3:
        r_hyperexc_vat, p_hyperexc_vat = stats.pearsonr(
            patient_df_hyperexc['upper_hyperexc'],
            patient_df_hyperexc['vat1l_stress']
        )
        print(f"    Upper hyperexc vs VAT1L stress:     r={r_hyperexc_vat:+.3f}, p={p_hyperexc_vat:.4f}")

    # Upper hyperexcitability vs Glia stress
    patient_df_hyperexc2 = patient_df.dropna(subset=['upper_hyperexc', 'glia_stress'])
    if len(patient_df_hyperexc2) >= 3:
        r_hyperexc_glia, p_hyperexc_glia = stats.pearsonr(
            patient_df_hyperexc2['upper_hyperexc'],
            patient_df_hyperexc2['glia_stress']
        )
        print(f"    Upper hyperexc vs Glia stress:      r={r_hyperexc_glia:+.3f}, p={p_hyperexc_glia:.4f}")

# ============================================================================
# Step 3: Cross-Cell Module Coupling
# ============================================================================
print("\n[Step 3] Cross-cell module coupling analysis...")

# Calculate cell-level coupling (within each patient)
coupling_results = []

for patient_id in df_als['patient_id'].unique():
    patient_cells = df_als[df_als['patient_id'] == patient_id]

    upper_cells = patient_cells[patient_cells['cell_type'].isin(upper_types)]
    vat1l_cells = patient_cells[patient_cells['cell_type'].isin(vat1l_types)]
    glia_cells = patient_cells[patient_cells['cell_type'].isin(glia_types)]

    if len(upper_cells) >= 10 and len(vat1l_cells) >= 5 and len(glia_cells) >= 10:
        # Upper hyperexcitability distribution
        hyperexc_modules = ['module_Synaptic', 'module_Calcium_Signaling', 'module_Ion_Transport']
        available_hyperexc = [m for m in hyperexc_modules if m in upper_cells.columns]

        if len(available_hyperexc) > 0:
            upper_hyperexc_mean = upper_cells[available_hyperexc].mean(axis=1).mean()
            upper_hyperexc_std = upper_cells[available_hyperexc].mean(axis=1).std()

            # VAT1L stress
            vat1l_stress_mean = vat1l_cells['stress_total'].mean()

            # Glia stress
            glia_stress_mean = glia_cells['stress_total'].mean()

            coupling_results.append({
                'patient_id': patient_id,
                'upper_hyperexc_mean': upper_hyperexc_mean,
                'upper_hyperexc_std': upper_hyperexc_std,
                'vat1l_stress_mean': vat1l_stress_mean,
                'glia_stress_mean': glia_stress_mean
            })

coupling_df = pd.DataFrame(coupling_results)
coupling_df.to_csv(output_dir / 'cross_cell_coupling.csv', index=False)

if len(coupling_df) >= 3:
    r_coup, p_coup = stats.pearsonr(coupling_df['upper_hyperexc_mean'], coupling_df['vat1l_stress_mean'])
    print(f"\n  Cross-patient coupling (Upper hyperexc vs VAT1L stress):")
    print(f"    r = {r_coup:+.3f}, p = {p_coup:.4f}")

# ============================================================================
# Step 4: Mid-Bin Trajectory Analysis
# ============================================================================
print("\n[Step 4] Mid-bin trajectory analysis (Early → Mid → Late)...")

# Extract Upper cells by PT_dpt group
df_upper = df_als[df_als['cell_type'].isin(upper_types)].copy()

early_upper = df_upper[df_upper['PT_dpt_group'] == 'Early']
mid_upper = df_upper[df_upper['PT_dpt_group'] == 'Mid']
late_upper = df_upper[df_upper['PT_dpt_group'] == 'Late']

print(f"\n  Upper subgroups:")
print(f"    Early: {len(early_upper):,} cells ({len(early_upper)/len(df_upper)*100:.1f}%)")
print(f"    Mid:   {len(mid_upper):,} cells ({len(mid_upper)/len(df_upper)*100:.1f}%)")
print(f"    Late:  {len(late_upper):,} cells ({len(late_upper)/len(df_upper)*100:.1f}%)")

# Key signatures
signatures = {
    'Hyperexcitability': ['module_Synaptic', 'module_Calcium_Signaling', 'module_Ion_Transport'],
    'ER_Stress': ['module_ER_Stress', 'module_Protein_Homeostasis'],
    'Inflammation': ['module_Inflammation'],
    'Apoptosis': ['module_Apoptosis']
}

trajectory_results = []

for sig_name, modules in signatures.items():
    available_mods = [m for m in modules if m in df_upper.columns]

    if len(available_mods) > 0:
        early_score = early_upper[available_mods].mean(axis=1).mean()
        mid_score = mid_upper[available_mods].mean(axis=1).mean() if len(mid_upper) > 0 else np.nan
        late_score = late_upper[available_mods].mean(axis=1).mean()

        trajectory_results.append({
            'signature': sig_name,
            'early_mean': early_score,
            'mid_mean': mid_score,
            'late_mean': late_score,
            'early_to_late_diff': early_score - late_score,
            'early_to_late_pct': ((early_score - late_score) / late_score) * 100 if late_score > 0 else np.nan
        })

        print(f"\n  {sig_name}:")
        print(f"    Early: {early_score:.4f}")
        if not np.isnan(mid_score):
            print(f"    Mid:   {mid_score:.4f}")
        print(f"    Late:  {late_score:.4f}")
        print(f"    Trend: {'Early > Mid > Late' if not np.isnan(mid_score) and early_score > mid_score > late_score else 'Early > Late'}")

trajectory_df = pd.DataFrame(trajectory_results)
trajectory_df.to_csv(output_dir / 'trajectory_analysis.csv', index=False)
print(f"\n  Saved: trajectory_analysis.csv")

# ============================================================================
# Step 5: Visualization
# ============================================================================
print("\n[Step 5] Creating visualizations...")

# Figure 1: Patient-level correlations
fig, axes = plt.subplots(2, 2, figsize=(14, 12))

if len(patient_df) >= 3:
    # Early Upper ratio vs VAT1L stress
    ax = axes[0, 0]
    ax.scatter(patient_df['early_upper_ratio'], patient_df['vat1l_stress'], alpha=0.7, s=100)
    z = np.polyfit(patient_df['early_upper_ratio'], patient_df['vat1l_stress'], 1)
    p_fit = np.poly1d(z)
    x_fit = np.linspace(patient_df['early_upper_ratio'].min(), patient_df['early_upper_ratio'].max(), 100)
    ax.plot(x_fit, p_fit(x_fit), "r--", alpha=0.8)
    ax.set_xlabel('Early Upper Ratio (%)')
    ax.set_ylabel('VAT1L Stress (mean)')
    ax.set_title(f'Early Upper Ratio vs VAT1L Stress\nr={r_vat:.3f}, p={p_vat:.4f}')
    ax.grid(alpha=0.3)

    # Early Upper ratio vs Glia stress
    ax = axes[0, 1]
    ax.scatter(patient_df['early_upper_ratio'], patient_df['glia_stress'], alpha=0.7, s=100, color='purple')
    z = np.polyfit(patient_df['early_upper_ratio'], patient_df['glia_stress'], 1)
    p_fit = np.poly1d(z)
    x_fit = np.linspace(patient_df['early_upper_ratio'].min(), patient_df['early_upper_ratio'].max(), 100)
    ax.plot(x_fit, p_fit(x_fit), "r--", alpha=0.8)
    ax.set_xlabel('Early Upper Ratio (%)')
    ax.set_ylabel('Glia Stress (mean)')
    ax.set_title(f'Early Upper Ratio vs Glia Stress\nr={r_glia:.3f}, p={p_glia:.4f}')
    ax.grid(alpha=0.3)

    # Upper hyperexc vs VAT1L stress
    if len(patient_df_hyperexc) >= 3:
        ax = axes[1, 0]
        ax.scatter(patient_df_hyperexc['upper_hyperexc'], patient_df_hyperexc['vat1l_stress'],
                   alpha=0.7, s=100, color='green')
        z = np.polyfit(patient_df_hyperexc['upper_hyperexc'], patient_df_hyperexc['vat1l_stress'], 1)
        p_fit = np.poly1d(z)
        x_fit = np.linspace(patient_df_hyperexc['upper_hyperexc'].min(),
                            patient_df_hyperexc['upper_hyperexc'].max(), 100)
        ax.plot(x_fit, p_fit(x_fit), "r--", alpha=0.8)
        ax.set_xlabel('Upper Hyperexcitability Score')
        ax.set_ylabel('VAT1L Stress (mean)')
        ax.set_title(f'Upper Hyperexc vs VAT1L Stress\nr={r_hyperexc_vat:.3f}, p={p_hyperexc_vat:.4f}')
        ax.grid(alpha=0.3)

    # Upper hyperexc vs Glia stress
    if len(patient_df_hyperexc2) >= 3:
        ax = axes[1, 1]
        ax.scatter(patient_df_hyperexc2['upper_hyperexc'], patient_df_hyperexc2['glia_stress'],
                   alpha=0.7, s=100, color='orange')
        z = np.polyfit(patient_df_hyperexc2['upper_hyperexc'], patient_df_hyperexc2['glia_stress'], 1)
        p_fit = np.poly1d(z)
        x_fit = np.linspace(patient_df_hyperexc2['upper_hyperexc'].min(),
                            patient_df_hyperexc2['upper_hyperexc'].max(), 100)
        ax.plot(x_fit, p_fit(x_fit), "r--", alpha=0.8)
        ax.set_xlabel('Upper Hyperexcitability Score')
        ax.set_ylabel('Glia Stress (mean)')
        ax.set_title(f'Upper Hyperexc vs Glia Stress\nr={r_hyperexc_glia:.3f}, p={p_hyperexc_glia:.4f}')
        ax.grid(alpha=0.3)

plt.tight_layout()
plt.savefig(output_dir / 'Fig_patient_level_correlations.png', dpi=300, bbox_inches='tight')
print(f"  Saved: Fig_patient_level_correlations.png")

# Figure 2: Trajectory plot
fig, ax = plt.subplots(1, 1, figsize=(10, 6))

for idx, row in trajectory_df.iterrows():
    x = [0, 1, 2]  # Early, Mid, Late
    y = [row['early_mean'], row['mid_mean'], row['late_mean']]

    ax.plot(x, y, marker='o', markersize=8, label=row['signature'], linewidth=2)

ax.set_xticks([0, 1, 2])
ax.set_xticklabels(['Early', 'Mid', 'Late'])
ax.set_xlabel('PT_dpt Group')
ax.set_ylabel('Mean Signature Score')
ax.set_title('Functional Signature Trajectory: Early → Mid → Late Upper')
ax.legend()
ax.grid(alpha=0.3)

plt.tight_layout()
plt.savefig(output_dir / 'Fig_trajectory_signatures.png', dpi=300, bbox_inches='tight')
print(f"  Saved: Fig_trajectory_signatures.png")

# ============================================================================
# Summary
# ============================================================================
print("\n" + "="*80)
print("Phase 6 Extended: Driver Validation Completed!")
print("="*80)

print(f"\nKey Findings:")

if len(patient_df) >= 3:
    print(f"\n  1. Patient-level correlations:")
    print(f"     Early Upper ratio vs VAT1L stress:  r={r_vat:+.3f}, p={p_vat:.4f}")
    print(f"     Early Upper ratio vs Glia stress:   r={r_glia:+.3f}, p={p_glia:.4f}")

    if abs(r_vat) > 0.3 or abs(r_glia) > 0.3:
        print(f"\n     → Moderate-to-strong correlation detected")
        if (r_vat > 0.3 and p_vat < 0.05) or (r_glia > 0.3 and p_glia < 0.05):
            print(f"       Supports 'functional driver' hypothesis ✓")
    else:
        print(f"\n     → Weak correlation")
        print(f"       Does NOT strongly support 'functional driver' hypothesis")

print(f"\n  2. Trajectory analysis:")
for idx, row in trajectory_df.iterrows():
    trend = 'Decreasing' if row['early_mean'] > row['late_mean'] else 'Increasing'
    print(f"     {row['signature']:20s}: {trend} (Early: {row['early_mean']:.2f} → Late: {row['late_mean']:.2f})")

print(f"\nGenerated files:")
print(f"  - patient_level_metrics.csv")
print(f"  - cross_cell_coupling.csv")
print(f"  - trajectory_analysis.csv")
print(f"  - Fig_patient_level_correlations.png")
print(f"  - Fig_trajectory_signatures.png")
