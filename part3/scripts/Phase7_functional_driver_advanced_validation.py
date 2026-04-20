#!/usr/bin/env python3
"""
Phase 7: Advanced Functional Driver Validation
==============================================

Goal: Test Early Upper functional driver hypothesis with advanced methods:
1. Non-linear/threshold analysis
2. Module-specific coupling (not just stress_total)
3. Within-patient heterogeneity
4. Subpopulation refinement
5. Cross-cell-type module coupling

This addresses the critical limitation: stress_total alone cannot capture
excitatory network dysfunction (hyperexcitability, Ca²⁺, glutamate, etc.)
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from scipy.interpolate import UnivariateSpline
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import PolynomialFeatures
import warnings
warnings.filterwarnings('ignore')

# Create output directory
import os
output_dir = 'results/phase7_advanced_driver_validation'
os.makedirs(output_dir, exist_ok=True)

print("="*80)
print("Phase 7: Advanced Functional Driver Validation")
print("="*80)

# Load data
print("\n[1] Loading data...")
df = pd.read_csv('results/PTv2_robustness/PTv2_quick_key_celltypes_with_PTs.csv')

# Define cell types (based on actual data)
upper_cts = ['Ex.L2.L3.CUX2.RASGRF2']  # Only this Upper-layer type in PTv2 data
vat1l_cts = ['Ex.L5.VAT1L.EYA4', 'Ex.L5.VAT1L.THSD4']
glia_cts = ['Glia.Astro.GFAP-neg', 'Glia.Micro', 'Glia.Oligo']

print(f"Total cells: {len(df):,}")
print(f"Upper cells: {df['cell_type'].isin(upper_cts).sum():,}")
print(f"VAT1L cells: {df['cell_type'].isin(vat1l_cts).sum():,}")
print(f"Glia cells: {df['cell_type'].isin(glia_cts).sum():,}")

# Define modules (with module_ prefix)
module_names = [
    'Synaptic', 'Calcium_Signaling', 'Ion_Transport',
    'ER_Stress', 'Inflammation', 'Apoptosis', 'Oxidative_Stress',
    'Mitochondria', 'Metabolism', 'Autophagy', 'Protein_Homeostasis',
    'DNA_Repair', 'Cell_Cycle', 'Transcription', 'Epigenetic',
    'RNA_Processing', 'lncRNA', 'Growth_Factors', 'Cytoskeleton',
    'ECM', 'Angiogenesis', 'Myelination', 'Complement'
]

# Rename columns to remove module_ prefix for easier handling
for mod in module_names:
    if f'module_{mod}' in df.columns:
        df[mod] = df[f'module_{mod}']

# Hyperexcitability composite (if not present, create it)
if 'Hyperexcitability' not in df.columns:
    df['Hyperexcitability'] = (df['Synaptic'] + df['Calcium_Signaling'] + df['Ion_Transport']) / 3

print("\n[2] Binning Upper cells by PT_dpt...")
upper_cells = df[df['cell_type'].isin(upper_cts)].copy()

# PT_dpt bins
bins = [0, 0.005, 0.010, 1.0]
labels = ['Early', 'Mid', 'Late']
upper_cells['pt_bin'] = pd.cut(upper_cells['PT_dpt'], bins=bins, labels=labels)

print(f"\nUpper cell distribution:")
print(upper_cells['pt_bin'].value_counts().sort_index())

# Flag Early Upper cells
df['is_early_upper'] = (df['cell_type'].isin(upper_cts)) & (df['PT_dpt'] < 0.005)
print(f"\nEarly Upper cells: {df['is_early_upper'].sum():,} ({100*df['is_early_upper'].sum()/len(upper_cells):.1f}% of Upper)")

print("\n" + "="*80)
print("ANALYSIS 1: Non-linear and Threshold Effects")
print("="*80)

# Patient-level metrics
print("\n[1.1] Calculating patient-level metrics...")
patient_metrics = []

for patient_id in df['patient_id'].unique():
    patient_cells = df[df['patient_id'] == patient_id]

    # Cell counts
    n_cells = len(patient_cells)
    n_upper = patient_cells['cell_type'].isin(upper_cts).sum()
    n_vat1l = patient_cells['cell_type'].isin(vat1l_cts).sum()
    n_glia = patient_cells['cell_type'].isin(glia_cts).sum()
    n_early_upper = patient_cells['is_early_upper'].sum()

    if n_upper < 10:  # Skip patients with too few Upper cells
        continue

    # Early Upper ratio
    early_upper_ratio = 100 * n_early_upper / n_upper if n_upper > 0 else 0

    # Upper metrics
    upper_patient = patient_cells[patient_cells['cell_type'].isin(upper_cts)]
    upper_stress = upper_patient['stress_total'].mean()
    upper_hyperexc = upper_patient['Hyperexcitability'].mean()
    upper_ca = upper_patient['Calcium_Signaling'].mean()
    upper_synaptic = upper_patient['Synaptic'].mean()
    upper_er = upper_patient['ER_Stress'].mean()

    # VAT1L metrics
    if n_vat1l >= 5:
        vat1l_patient = patient_cells[patient_cells['cell_type'].isin(vat1l_cts)]
        vat1l_stress = vat1l_patient['stress_total'].mean()
        vat1l_er = vat1l_patient['ER_Stress'].mean()
        vat1l_ca = vat1l_patient['Calcium_Signaling'].mean()
        vat1l_apoptosis = vat1l_patient['Apoptosis'].mean()
    else:
        vat1l_stress = np.nan
        vat1l_er = np.nan
        vat1l_ca = np.nan
        vat1l_apoptosis = np.nan

    # Glia metrics
    if n_glia >= 10:
        glia_patient = patient_cells[patient_cells['cell_type'].isin(glia_cts)]
        glia_stress = glia_patient['stress_total'].mean()
        glia_inflammation = glia_patient['Inflammation'].mean()
        glia_er = glia_patient['ER_Stress'].mean()
    else:
        glia_stress = np.nan
        glia_inflammation = np.nan
        glia_er = np.nan

    patient_metrics.append({
        'patient_id': patient_id,
        'n_cells': n_cells,
        'n_upper': n_upper,
        'n_vat1l': n_vat1l,
        'n_glia': n_glia,
        'n_early_upper': n_early_upper,
        'early_upper_ratio': early_upper_ratio,
        'upper_stress': upper_stress,
        'upper_hyperexc': upper_hyperexc,
        'upper_ca': upper_ca,
        'upper_synaptic': upper_synaptic,
        'upper_er': upper_er,
        'vat1l_stress': vat1l_stress,
        'vat1l_er': vat1l_er,
        'vat1l_ca': vat1l_ca,
        'vat1l_apoptosis': vat1l_apoptosis,
        'glia_stress': glia_stress,
        'glia_inflammation': glia_inflammation,
        'glia_er': glia_er
    })

patient_df = pd.DataFrame(patient_metrics)
patient_df = patient_df.dropna(subset=['vat1l_stress', 'glia_stress'])  # Keep complete cases
print(f"Patients with complete data: {len(patient_df)}")

# Save patient metrics
patient_df.to_csv(f'{output_dir}/patient_level_advanced_metrics.csv', index=False)

print("\n[1.2] Testing threshold effects...")
# Test if Early Upper ratio >10% shows different VAT1L stress
patient_df['high_early_upper'] = patient_df['early_upper_ratio'] > 10

high_group = patient_df[patient_df['high_early_upper']]
low_group = patient_df[~patient_df['high_early_upper']]

print(f"\nHigh Early Upper (>10%): n={len(high_group)}")
print(f"  VAT1L stress: {high_group['vat1l_stress'].mean():.3f} ± {high_group['vat1l_stress'].std():.3f}")
print(f"  Glia stress: {high_group['glia_stress'].mean():.3f} ± {high_group['glia_stress'].std():.3f}")

print(f"\nLow Early Upper (≤10%): n={len(low_group)}")
print(f"  VAT1L stress: {low_group['vat1l_stress'].mean():.3f} ± {low_group['vat1l_stress'].std():.3f}")
print(f"  Glia stress: {low_group['glia_stress'].mean():.3f} ± {low_group['glia_stress'].std():.3f}")

if len(high_group) >= 3 and len(low_group) >= 3:
    t_vat1l, p_vat1l = stats.ttest_ind(high_group['vat1l_stress'], low_group['vat1l_stress'])
    t_glia, p_glia = stats.ttest_ind(high_group['glia_stress'], low_group['glia_stress'])
    print(f"\nThreshold test (>10% vs ≤10%):")
    print(f"  VAT1L stress: t={t_vat1l:.3f}, p={p_vat1l:.3f}")
    print(f"  Glia stress: t={t_glia:.3f}, p={p_glia:.3f}")

print("\n[1.3] Polynomial regression (quadratic, cubic)...")
# Test non-linear relationships
X = patient_df['early_upper_ratio'].values.reshape(-1, 1)

results_nonlinear = []

for target, target_name in [
    ('vat1l_stress', 'VAT1L stress'),
    ('vat1l_er', 'VAT1L ER stress'),
    ('glia_stress', 'Glia stress')
]:
    y = patient_df[target].values

    # Linear
    r_linear, p_linear = stats.pearsonr(patient_df['early_upper_ratio'], patient_df[target])

    # Quadratic
    poly2 = PolynomialFeatures(degree=2)
    X_poly2 = poly2.fit_transform(X)
    model2 = LinearRegression().fit(X_poly2, y)
    r2_poly2 = model2.score(X_poly2, y)

    # Cubic
    poly3 = PolynomialFeatures(degree=3)
    X_poly3 = poly3.fit_transform(X)
    model3 = LinearRegression().fit(X_poly3, y)
    r2_poly3 = model3.score(X_poly3, y)

    print(f"\n{target_name}:")
    print(f"  Linear: r={r_linear:.3f}, p={p_linear:.3f}, R²={r_linear**2:.3f}")
    print(f"  Quadratic: R²={r2_poly2:.3f} (Δ={r2_poly2 - r_linear**2:.3f})")
    print(f"  Cubic: R²={r2_poly3:.3f} (Δ={r2_poly3 - r_linear**2:.3f})")

    results_nonlinear.append({
        'target': target_name,
        'r_linear': r_linear,
        'p_linear': p_linear,
        'R2_linear': r_linear**2,
        'R2_quadratic': r2_poly2,
        'R2_cubic': r2_poly3,
        'improvement_quadratic': r2_poly2 - r_linear**2,
        'improvement_cubic': r2_poly3 - r_linear**2
    })

results_nonlinear_df = pd.DataFrame(results_nonlinear)
results_nonlinear_df.to_csv(f'{output_dir}/nonlinear_regression_results.csv', index=False)

print("\n" + "="*80)
print("ANALYSIS 2: Module-Specific Coupling (Beyond stress_total)")
print("="*80)

print("\n[2.1] Testing module-specific correlations...")
# Key hypothesis: Upper hyperexcitability → VAT1L Ca²⁺/ER stress
# (Not just stress_total, but specific functional pathways)

module_correlations = []

# Upper modules of interest
upper_modules = ['upper_hyperexc', 'upper_ca', 'upper_synaptic', 'upper_er']

# VAT1L/Glia targets
targets = [
    ('vat1l_stress', 'VAT1L stress_total'),
    ('vat1l_er', 'VAT1L ER stress'),
    ('vat1l_ca', 'VAT1L Ca²⁺'),
    ('vat1l_apoptosis', 'VAT1L Apoptosis'),
    ('glia_stress', 'Glia stress_total'),
    ('glia_er', 'Glia ER stress'),
    ('glia_inflammation', 'Glia Inflammation')
]

print("\nCorrelation matrix (Upper modules → VAT1L/Glia):")
print(f"{'Upper module':<20} {'Target':<25} {'r':>8} {'p':>8} {'n':>4}")
print("-" * 70)

for upper_mod in upper_modules:
    for target, target_name in targets:
        valid = patient_df[[upper_mod, target]].dropna()
        if len(valid) >= 5:
            r, p = stats.pearsonr(valid[upper_mod], valid[target])
            module_correlations.append({
                'upper_module': upper_mod,
                'target': target_name,
                'r': r,
                'p': p,
                'n': len(valid),
                'significant': p < 0.05
            })

            sig_marker = '*' if p < 0.05 else ''
            print(f"{upper_mod:<20} {target_name:<25} {r:>8.3f} {p:>8.3f} {len(valid):>4} {sig_marker}")

module_corr_df = pd.DataFrame(module_correlations)
module_corr_df = module_corr_df.sort_values('p')
module_corr_df.to_csv(f'{output_dir}/module_specific_correlations.csv', index=False)

print("\n[2.2] Strongest module-specific correlations:")
top_corrs = module_corr_df.nsmallest(5, 'p')
for idx, row in top_corrs.iterrows():
    print(f"  {row['upper_module']} → {row['target']}: r={row['r']:.3f}, p={row['p']:.3f}")

print("\n" + "="*80)
print("ANALYSIS 3: Within-Patient Heterogeneity")
print("="*80)

print("\n[3.1] Calculating Upper cell heterogeneity within patients...")
# Hypothesis: Patients with high Upper heterogeneity (stressed subpopulations)
# should have higher VAT1L/Glia stress IF Early Upper drives pathology

patient_heterogeneity = []

for patient_id in df['patient_id'].unique():
    patient_cells = df[df['patient_id'] == patient_id]
    upper_patient = patient_cells[patient_cells['cell_type'].isin(upper_cts)]

    if len(upper_patient) < 20:  # Need sufficient cells
        continue

    # Upper heterogeneity metrics
    upper_stress_sd = upper_patient['stress_total'].std()
    upper_hyperexc_sd = upper_patient['Hyperexcitability'].std()
    upper_ca_sd = upper_patient['Calcium_Signaling'].std()

    # Inter-quantile range (IQR) as robust heterogeneity measure
    upper_stress_iqr = upper_patient['stress_total'].quantile(0.75) - upper_patient['stress_total'].quantile(0.25)
    upper_hyperexc_iqr = upper_patient['Hyperexcitability'].quantile(0.75) - upper_patient['Hyperexcitability'].quantile(0.25)

    # Top 10% vs bottom 90% difference (extreme subpopulation effect)
    top10_stress = upper_patient['stress_total'].quantile(0.90)
    bottom90_stress = upper_patient[upper_patient['stress_total'] <= upper_patient['stress_total'].quantile(0.90)]['stress_total'].mean()
    extreme_stress_diff = top10_stress - bottom90_stress

    top10_hyperexc = upper_patient['Hyperexcitability'].quantile(0.90)
    bottom90_hyperexc = upper_patient[upper_patient['Hyperexcitability'] <= upper_patient['Hyperexcitability'].quantile(0.90)]['Hyperexcitability'].mean()
    extreme_hyperexc_diff = top10_hyperexc - bottom90_hyperexc

    patient_heterogeneity.append({
        'patient_id': patient_id,
        'n_upper': len(upper_patient),
        'upper_stress_sd': upper_stress_sd,
        'upper_hyperexc_sd': upper_hyperexc_sd,
        'upper_ca_sd': upper_ca_sd,
        'upper_stress_iqr': upper_stress_iqr,
        'upper_hyperexc_iqr': upper_hyperexc_iqr,
        'extreme_stress_diff': extreme_stress_diff,
        'extreme_hyperexc_diff': extreme_hyperexc_diff
    })

hetero_df = pd.DataFrame(patient_heterogeneity)
patient_df = patient_df.merge(hetero_df, on='patient_id', how='left')

print(f"Patients with heterogeneity data: {len(hetero_df)}")

print("\n[3.2] Testing heterogeneity → VAT1L/Glia stress...")
heterogeneity_corrs = []

for hetero_metric in ['upper_stress_sd', 'upper_hyperexc_sd', 'upper_hyperexc_iqr', 'extreme_hyperexc_diff']:
    for target, target_name in [('vat1l_stress', 'VAT1L stress'),
                                 ('glia_stress', 'Glia stress'),
                                 ('vat1l_er', 'VAT1L ER stress')]:
        valid = patient_df[[hetero_metric, target]].dropna()
        if len(valid) >= 5:
            r, p = stats.pearsonr(valid[hetero_metric], valid[target])
            heterogeneity_corrs.append({
                'heterogeneity_metric': hetero_metric,
                'target': target_name,
                'r': r,
                'p': p,
                'n': len(valid)
            })
            print(f"{hetero_metric:<25} → {target_name:<20}: r={r:>6.3f}, p={p:>6.3f}")

hetero_corr_df = pd.DataFrame(heterogeneity_corrs)
hetero_corr_df.to_csv(f'{output_dir}/heterogeneity_correlations.csv', index=False)

print("\n" + "="*80)
print("ANALYSIS 4: Subpopulation Refinement")
print("="*80)

print("\n[4.1] Refining Early Upper into extreme vs moderate...")
# Is there a "hyper-extreme" top 3% that drives effects?

early_upper_cells = df[df['is_early_upper']].copy()
print(f"Total Early Upper cells: {len(early_upper_cells):,}")

# Define top 3% hyperexcitable Upper
all_upper = df[df['cell_type'].isin(upper_cts)].copy()
hyperexc_99th = all_upper['Hyperexcitability'].quantile(0.97)  # Top 3%
hyperexc_95th = all_upper['Hyperexcitability'].quantile(0.95)  # Top 5%

df['is_extreme_upper'] = (df['cell_type'].isin(upper_cts)) & (df['Hyperexcitability'] >= hyperexc_99th)
df['is_very_high_upper'] = (df['cell_type'].isin(upper_cts)) & (df['Hyperexcitability'] >= hyperexc_95th)

print(f"Extreme Upper (top 3% hyperexc): {df['is_extreme_upper'].sum():,} cells")
print(f"Very High Upper (top 5% hyperexc): {df['is_very_high_upper'].sum():,} cells")

# Patient-level extreme ratios
patient_extreme = []
for patient_id in df['patient_id'].unique():
    patient_cells = df[df['patient_id'] == patient_id]
    n_upper = patient_cells['cell_type'].isin(upper_cts).sum()
    n_extreme = patient_cells['is_extreme_upper'].sum()
    n_very_high = patient_cells['is_very_high_upper'].sum()

    if n_upper >= 10:
        patient_extreme.append({
            'patient_id': patient_id,
            'extreme_upper_ratio': 100 * n_extreme / n_upper,
            'very_high_upper_ratio': 100 * n_very_high / n_upper
        })

extreme_df = pd.DataFrame(patient_extreme)
patient_df = patient_df.merge(extreme_df, on='patient_id', how='left')

print("\n[4.2] Testing extreme Upper ratio → VAT1L/Glia stress...")
for ratio_var, ratio_name in [('extreme_upper_ratio', 'Extreme (top 3%)'),
                               ('very_high_upper_ratio', 'Very High (top 5%)')]:
    print(f"\n{ratio_name}:")
    for target, target_name in [('vat1l_stress', 'VAT1L stress'),
                                 ('vat1l_er', 'VAT1L ER stress'),
                                 ('glia_stress', 'Glia stress')]:
        valid = patient_df[[ratio_var, target]].dropna()
        if len(valid) >= 5:
            r, p = stats.pearsonr(valid[ratio_var], valid[target])
            print(f"  {ratio_name} → {target_name}: r={r:.3f}, p={p:.3f}")

print("\n[4.3] Comparing Early Upper subgroups...")
# Within Early Upper, are there distinct subgroups?
early_upper_cells['hyperexc_tertile'] = pd.qcut(early_upper_cells['Hyperexcitability'],
                                                  q=3, labels=['Low', 'Mid', 'High'])

print("\nEarly Upper tertiles:")
for tertile in ['Low', 'Mid', 'High']:
    subset = early_upper_cells[early_upper_cells['hyperexc_tertile'] == tertile]
    print(f"\n{tertile} hyperexcitability (n={len(subset)}):")
    print(f"  Hyperexcitability: {subset['Hyperexcitability'].mean():.2f} ± {subset['Hyperexcitability'].std():.2f}")
    print(f"  ER Stress: {subset['ER_Stress'].mean():.2f} ± {subset['ER_Stress'].std():.2f}")
    print(f"  Apoptosis: {subset['Apoptosis'].mean():.2f} ± {subset['Apoptosis'].std():.2f}")
    print(f"  stress_total: {subset['stress_total'].mean():.2f} ± {subset['stress_total'].std():.2f}")

print("\n" + "="*80)
print("ANALYSIS 5: Cross-Cell-Type Module Coupling")
print("="*80)

print("\n[5.1] Cell-level correlation: Early Upper → VAT1L modules...")
# For each patient, correlate Upper cell modules with VAT1L cell modules
# (within-patient cell-level correlation)

patient_cell_corrs = []

for patient_id in df['patient_id'].unique():
    patient_cells = df[df['patient_id'] == patient_id]

    upper_patient = patient_cells[patient_cells['cell_type'].isin(upper_cts)]
    vat1l_patient = patient_cells[patient_cells['cell_type'].isin(vat1l_cts)]

    if len(upper_patient) < 20 or len(vat1l_patient) < 5:
        continue

    # Average Upper modules (patient-level)
    upper_mean_hyperexc = upper_patient['Hyperexcitability'].mean()
    upper_mean_ca = upper_patient['Calcium_Signaling'].mean()
    upper_mean_synaptic = upper_patient['Synaptic'].mean()

    # VAT1L modules (patient-level)
    vat1l_mean_er = vat1l_patient['ER_Stress'].mean()
    vat1l_mean_ca = vat1l_patient['Calcium_Signaling'].mean()
    vat1l_mean_apoptosis = vat1l_patient['Apoptosis'].mean()
    vat1l_mean_stress = vat1l_patient['stress_total'].mean()

    patient_cell_corrs.append({
        'patient_id': patient_id,
        'upper_hyperexc': upper_mean_hyperexc,
        'upper_ca': upper_mean_ca,
        'upper_synaptic': upper_mean_synaptic,
        'vat1l_er': vat1l_mean_er,
        'vat1l_ca': vat1l_mean_ca,
        'vat1l_apoptosis': vat1l_mean_apoptosis,
        'vat1l_stress': vat1l_mean_stress
    })

cell_corr_df = pd.DataFrame(patient_cell_corrs)

print(f"Patients with Upper+VAT1L data: {len(cell_corr_df)}")

# Test specific pathway hypotheses
print("\n[5.2] Pathway-specific correlations:")

pathway_tests = [
    ('upper_hyperexc', 'vat1l_er', 'Upper Hyperexc → VAT1L ER stress'),
    ('upper_ca', 'vat1l_ca', 'Upper Ca²⁺ → VAT1L Ca²⁺'),
    ('upper_ca', 'vat1l_apoptosis', 'Upper Ca²⁺ → VAT1L Apoptosis'),
    ('upper_synaptic', 'vat1l_er', 'Upper Synaptic → VAT1L ER stress'),
    ('upper_synaptic', 'vat1l_apoptosis', 'Upper Synaptic → VAT1L Apoptosis')
]

pathway_results = []

for x_var, y_var, pathway_name in pathway_tests:
    if x_var in cell_corr_df.columns and y_var in cell_corr_df.columns:
        valid = cell_corr_df[[x_var, y_var]].dropna()
        if len(valid) >= 5:
            r, p = stats.pearsonr(valid[x_var], valid[y_var])
            pathway_results.append({
                'pathway': pathway_name,
                'r': r,
                'p': p,
                'n': len(valid)
            })

            sig_marker = '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else ''
            print(f"  {pathway_name:<40}: r={r:>6.3f}, p={p:>6.3f} {sig_marker}")

pathway_df = pd.DataFrame(pathway_results)
pathway_df.to_csv(f'{output_dir}/pathway_specific_correlations.csv', index=False)

print("\n" + "="*80)
print("SUMMARY OF ADVANCED VALIDATION")
print("="*80)

print("\n[Summary] Key findings:")

# 1. Threshold effects
print("\n1. THRESHOLD EFFECTS:")
if len(high_group) >= 3 and len(low_group) >= 3:
    vat1l_diff = high_group['vat1l_stress'].mean() - low_group['vat1l_stress'].mean()
    print(f"   Early Upper >10% vs ≤10%:")
    print(f"   VAT1L stress difference: {vat1l_diff:+.3f} (p={p_vat1l:.3f})")
    if abs(vat1l_diff) > 0.1:
        print(f"   → Possible threshold effect detected")
    else:
        print(f"   → No strong threshold effect")

# 2. Module-specific coupling
print("\n2. MODULE-SPECIFIC COUPLING:")
sig_module_corrs = module_corr_df[module_corr_df['p'] < 0.10].sort_values('p')
if len(sig_module_corrs) > 0:
    print(f"   {len(sig_module_corrs)} correlations with p<0.10:")
    for idx, row in sig_module_corrs.head(3).iterrows():
        print(f"   - {row['upper_module']} → {row['target']}: r={row['r']:.3f}, p={row['p']:.3f}")
else:
    print(f"   No significant module-specific correlations (p<0.10)")

# 3. Heterogeneity
print("\n3. WITHIN-PATIENT HETEROGENEITY:")
hetero_sig = hetero_corr_df[hetero_corr_df['p'] < 0.10].sort_values('p')
if len(hetero_sig) > 0:
    print(f"   {len(hetero_sig)} heterogeneity correlations with p<0.10:")
    for idx, row in hetero_sig.head(3).iterrows():
        print(f"   - {row['heterogeneity_metric']} → {row['target']}: r={row['r']:.3f}, p={row['p']:.3f}")
else:
    print(f"   No significant heterogeneity correlations (p<0.10)")

# 4. Pathway-specific
print("\n4. PATHWAY-SPECIFIC COUPLING:")
pathway_sig = pathway_df[pathway_df['p'] < 0.10].sort_values('p')
if len(pathway_sig) > 0:
    print(f"   {len(pathway_sig)} pathways with p<0.10:")
    for idx, row in pathway_sig.iterrows():
        print(f"   - {row['pathway']}: r={row['r']:.3f}, p={row['p']:.3f}")
else:
    print(f"   No significant pathway coupling (p<0.10)")

print("\n" + "="*80)
print("Analysis complete!")
print("="*80)
print(f"\nResults saved to: {output_dir}/")
print("\nKey outputs:")
print(f"  - patient_level_advanced_metrics.csv")
print(f"  - nonlinear_regression_results.csv")
print(f"  - module_specific_correlations.csv")
print(f"  - heterogeneity_correlations.csv")
print(f"  - pathway_specific_correlations.csv")
