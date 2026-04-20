#!/usr/bin/env python3
"""
Phase 9″: Patient-level Coupling of Vascular Driver-State

Goal:
    Quantify how the abundance of the vascular driver-state
    (identified in Phase 9′) relates to pathology in Upper,
    Glia, and VAT1L at the patient level.

Question:
    Do patients with more vascular driver-state cells
    show stronger Upper hyperexcitability, Glia ER/complement,
    and VAT1L Ca²⁺/ER stress?

Author: Claude Code
Date: 2025-11-24
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr, spearmanr
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import PolynomialFeatures
from sklearn.metrics import r2_score
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# ============================================================
# Configuration
# ============================================================

# Input files
CELL_LEVEL_FILE = 'results/phase9_vascular/cell_level_features_ALL_with_PTdpt.csv'
VDRIVER_CELLS_FILE = 'results/phase9p_vascular_driver/vascular_driver_state_cells.csv'

# Output directory
OUTPUT_DIR = 'results/phase9p_vascular_driver'
PLOT_DIR = f'{OUTPUT_DIR}/patient_coupling_plots'

# Create output directories
Path(OUTPUT_DIR).mkdir(parents=True, exist_ok=True)
Path(PLOT_DIR).mkdir(parents=True, exist_ok=True)

# Cell type definitions
VASCULAR_TYPES = [
    'Vasc.Endo.Capillary',
    'Vasc.Mural.Pericyte',
    'Vasc.Fibro.CLMP.PDGFRA',
    'Vasc.Endo.Venous',
    'Vasc.Mural.SMC',
    'Vasc.Endo.Arterial'
]

UPPER_TYPES = [
    'Ex.L2.L3.CUX2.RASGRF2'  # L2/L3 ONLY (corrected definition)
]

GLIA_TYPES = [
    'Glia.Oligo',
    'Glia.Astro.GFAP-neg',
    'Glia.Micro'
]

VAT1L_TYPES = [
    'Ex.L5.VAT1L.EYA4',
    'Ex.L5.VAT1L.THSD4'
]

print("=" * 80)
print("Phase 9″: Patient-level Vascular Driver Coupling Analysis")
print("=" * 80)
print()

# ============================================================
# Step 1: Load data and attach patient_id to vascular driver cells
# ============================================================

print("[1] Loading data...")

# Load full cell-level data
df_all = pd.read_csv(CELL_LEVEL_FILE)
print(f"  Total cells: {len(df_all):,}")

# Extract patient_id from cell_id (format: BARCODE-PATIENT_ID)
df_all['patient_id'] = df_all['cell_id'].str.split('-').str[1]
print(f"  Patients: {df_all['patient_id'].nunique()}")
print(f"  Patient IDs: {sorted(df_all['patient_id'].unique())}")

# Load vascular driver-state cells
df_vdriver_cells = pd.read_csv(VDRIVER_CELLS_FILE)
print(f"  Vascular driver cells: {len(df_vdriver_cells):,}")

# Merge to get patient_id for driver cells
df_vdriver = df_vdriver_cells.merge(
    df_all[['cell_id', 'patient_id', 'cell_type', 'PT_dpt', 'stress_total']],
    on='cell_id',
    how='left'
)

# Sanity check
missing_patient = df_vdriver['patient_id'].isna().sum()
if missing_patient > 0:
    print(f"  WARNING: {missing_patient} driver cells missing patient_id")
else:
    print(f"  ✓ All driver cells have patient_id")

print()

# ============================================================
# Step 2: Compute per-patient vascular driver abundance
# ============================================================

print("[2] Computing per-patient vascular driver abundance...")

# Get all vascular cells
df_vasc_all = df_all[df_all['cell_type'].isin(VASCULAR_TYPES)].copy()
print(f"  Total vascular cells: {len(df_vasc_all):,}")

# Count total vascular cells per patient
total_vasc_per_patient = df_vasc_all.groupby('patient_id').size().reset_index(name='total_vasc')

# Count driver vascular cells per patient
driver_vasc_per_patient = df_vdriver.groupby('patient_id').size().reset_index(name='driver_vasc')

# Merge
df_patient_vdriver = total_vasc_per_patient.merge(
    driver_vasc_per_patient,
    on='patient_id',
    how='left'
)
df_patient_vdriver['driver_vasc'] = df_patient_vdriver['driver_vasc'].fillna(0).astype(int)
df_patient_vdriver['driver_ratio'] = df_patient_vdriver['driver_vasc'] / df_patient_vdriver['total_vasc']

print(f"  Patients with vascular cells: {len(df_patient_vdriver)}")
print(f"  Driver ratio range: {df_patient_vdriver['driver_ratio'].min():.4f} - {df_patient_vdriver['driver_ratio'].max():.4f}")
print(f"  Mean driver ratio: {df_patient_vdriver['driver_ratio'].mean():.4f}")

# Save
output_file = f'{OUTPUT_DIR}/patient_vascular_driver_abundance.csv'
df_patient_vdriver.to_csv(output_file, index=False)
print(f"  Saved: {output_file}")
print()

# ============================================================
# Step 3: Compute per-patient pathology signatures
# ============================================================

print("[3] Computing per-patient pathology signatures...")

# Module columns
module_cols = [col for col in df_all.columns if col.startswith('module_')]

def compute_patient_signatures(df, cell_types, patient_col='patient_id'):
    """Compute per-patient module signatures for a cell type group."""
    df_group = df[df['cell_type'].isin(cell_types)].copy()

    # Group by patient and compute means
    signatures = {}
    for patient in df[patient_col].unique():
        patient_cells = df_group[df_group[patient_col] == patient]
        if len(patient_cells) == 0:
            continue

        signatures[patient] = {
            'n_cells': len(patient_cells),
            'stress_total': patient_cells['stress_total'].mean()
        }

        # Add module means
        for col in module_cols:
            if col in patient_cells.columns:
                signatures[patient][col.replace('module_', '')] = patient_cells[col].mean()

    return pd.DataFrame.from_dict(signatures, orient='index').reset_index().rename(columns={'index': patient_col})

# Compute signatures for each group
print("  Computing Upper signatures...")
df_upper = compute_patient_signatures(df_all, UPPER_TYPES)

print("  Computing Glia signatures...")
df_glia = compute_patient_signatures(df_all, GLIA_TYPES)

print("  Computing VAT1L signatures...")
df_vat1l = compute_patient_signatures(df_all, VAT1L_TYPES)

# Build pathology signatures
pathology_metrics = []
for _, row in df_patient_vdriver.iterrows():
    patient = row['patient_id']
    metrics = {'patient_id': patient}

    # Upper metrics
    upper_row = df_upper[df_upper['patient_id'] == patient]
    if len(upper_row) > 0:
        upper_row = upper_row.iloc[0]
        # Hyperexcitability composite
        hyperexc_mods = ['Synaptic', 'Calcium_Signaling', 'Ion_Transport']
        hyperexc_vals = [upper_row.get(m, np.nan) for m in hyperexc_mods]
        metrics['Upper_hyperexc'] = np.nanmean(hyperexc_vals)
        metrics['Upper_ER'] = upper_row.get('ER_Stress', np.nan)
        metrics['Upper_Mito'] = upper_row.get('Mitochondria', np.nan)
        metrics['Upper_n_cells'] = upper_row.get('n_cells', 0)
    else:
        metrics['Upper_hyperexc'] = np.nan
        metrics['Upper_ER'] = np.nan
        metrics['Upper_Mito'] = np.nan
        metrics['Upper_n_cells'] = 0

    # Glia metrics
    glia_row = df_glia[df_glia['patient_id'] == patient]
    if len(glia_row) > 0:
        glia_row = glia_row.iloc[0]
        metrics['Glia_ER'] = glia_row.get('ER_Stress', np.nan)
        metrics['Glia_Complement'] = glia_row.get('Complement', np.nan)
        metrics['Glia_Mito'] = glia_row.get('Mitochondria', np.nan)
        metrics['Glia_Inflammation'] = glia_row.get('Inflammation', np.nan)
        metrics['Glia_n_cells'] = glia_row.get('n_cells', 0)
    else:
        metrics['Glia_ER'] = np.nan
        metrics['Glia_Complement'] = np.nan
        metrics['Glia_Mito'] = np.nan
        metrics['Glia_Inflammation'] = np.nan
        metrics['Glia_n_cells'] = 0

    # VAT1L metrics
    vat1l_row = df_vat1l[df_vat1l['patient_id'] == patient]
    if len(vat1l_row) > 0:
        vat1l_row = vat1l_row.iloc[0]
        metrics['VAT1L_Ca'] = vat1l_row.get('Calcium_Signaling', np.nan)
        metrics['VAT1L_ER'] = vat1l_row.get('ER_Stress', np.nan)
        metrics['VAT1L_Stress'] = vat1l_row.get('stress_total', np.nan)
        metrics['VAT1L_n_cells'] = vat1l_row.get('n_cells', 0)
    else:
        metrics['VAT1L_Ca'] = np.nan
        metrics['VAT1L_ER'] = np.nan
        metrics['VAT1L_Stress'] = np.nan
        metrics['VAT1L_n_cells'] = 0

    pathology_metrics.append(metrics)

df_patient_pathology = pd.DataFrame(pathology_metrics)

# Save
output_file = f'{OUTPUT_DIR}/patient_pathology_signatures.csv'
df_patient_pathology.to_csv(output_file, index=False)
print(f"  Saved: {output_file}")
print()

# ============================================================
# Step 4: Merge and run correlation analysis
# ============================================================

print("[4] Running correlation analysis...")

# Merge
df_patient = df_patient_vdriver.merge(
    df_patient_pathology,
    on='patient_id',
    how='inner'
)

print(f"  Patients with complete data: {len(df_patient)}")

# Define predictor-target pairs
predictors = ['driver_ratio', 'driver_vasc']
targets = [
    'Upper_hyperexc', 'Upper_ER', 'Upper_Mito',
    'Glia_ER', 'Glia_Complement', 'Glia_Mito', 'Glia_Inflammation',
    'VAT1L_Ca', 'VAT1L_ER', 'VAT1L_Stress'
]

# Run correlations
correlation_results = []

for predictor in predictors:
    for target in targets:
        # Get valid data (non-NaN)
        valid_mask = ~(df_patient[predictor].isna() | df_patient[target].isna())
        if valid_mask.sum() < 3:
            continue

        x = df_patient.loc[valid_mask, predictor].values
        y = df_patient.loc[valid_mask, target].values
        n = len(x)

        # Pearson correlation
        r_pearson, p_pearson = pearsonr(x, y)

        # Spearman correlation
        rho_spearman, p_spearman = spearmanr(x, y)

        # Linear regression
        X_lin = x.reshape(-1, 1)
        lr = LinearRegression()
        lr.fit(X_lin, y)
        y_pred_lin = lr.predict(X_lin)
        r2_linear = r2_score(y, y_pred_lin)

        # Cubic polynomial regression
        poly = PolynomialFeatures(degree=3)
        X_poly = poly.fit_transform(X_lin)
        lr_poly = LinearRegression()
        lr_poly.fit(X_poly, y)
        y_pred_poly = lr_poly.predict(X_poly)
        r2_cubic = r2_score(y, y_pred_poly)

        correlation_results.append({
            'predictor': predictor,
            'target': target,
            'n': n,
            'pearson_r': r_pearson,
            'pearson_p': p_pearson,
            'spearman_rho': rho_spearman,
            'spearman_p': p_spearman,
            'R2_linear': r2_linear,
            'R2_cubic': r2_cubic,
            'R2_improvement': r2_cubic - r2_linear
        })

df_correlations = pd.DataFrame(correlation_results)

# Sort by absolute Pearson r
df_correlations['abs_r'] = df_correlations['pearson_r'].abs()
df_correlations = df_correlations.sort_values('abs_r', ascending=False)

# Save
output_file = f'{OUTPUT_DIR}/patient_vdriver_correlations.csv'
df_correlations.to_csv(output_file, index=False)
print(f"  Saved: {output_file}")

# Print top correlations
print("\n  Top correlations (|r| > 0.3 or p < 0.10):")
strong_corrs = df_correlations[
    (df_correlations['abs_r'] > 0.3) | (df_correlations['pearson_p'] < 0.10)
]
if len(strong_corrs) > 0:
    for _, row in strong_corrs.iterrows():
        sig = '***' if row['pearson_p'] < 0.01 else '**' if row['pearson_p'] < 0.05 else '*' if row['pearson_p'] < 0.10 else ''
        print(f"    {row['predictor']:15s} vs {row['target']:20s}: r={row['pearson_r']:+.3f} (p={row['pearson_p']:.3f}){sig}, R²={row['R2_linear']:.3f}")
else:
    print("    No strong correlations found")

print()

# ============================================================
# Step 5: Visualizations
# ============================================================

print("[5] Creating visualizations...")

# Plot interesting correlations
plot_threshold_r = 0.3
plot_threshold_r2_improvement = 0.05

interesting = df_correlations[
    (df_correlations['abs_r'] > plot_threshold_r) |
    (df_correlations['R2_improvement'] > plot_threshold_r2_improvement)
]

print(f"  Creating {len(interesting)} scatter plots...")

for idx, row in interesting.iterrows():
    predictor = row['predictor']
    target = row['target']

    # Get valid data
    valid_mask = ~(df_patient[predictor].isna() | df_patient[target].isna())
    x = df_patient.loc[valid_mask, predictor].values
    y = df_patient.loc[valid_mask, target].values

    if len(x) < 3:
        continue

    # Create figure
    fig, ax = plt.subplots(figsize=(8, 6))

    # Scatter plot
    ax.scatter(x, y, s=100, alpha=0.6, edgecolors='black', linewidths=1)

    # Linear fit
    X_lin = x.reshape(-1, 1)
    lr = LinearRegression()
    lr.fit(X_lin, y)
    x_range = np.linspace(x.min(), x.max(), 100).reshape(-1, 1)
    y_pred_lin = lr.predict(x_range)
    ax.plot(x_range, y_pred_lin, 'r--', linewidth=2, label=f'Linear (R²={row["R2_linear"]:.3f})')

    # Cubic fit if improvement is substantial
    if row['R2_improvement'] > 0.05:
        poly = PolynomialFeatures(degree=3)
        X_poly = poly.fit_transform(X_lin)
        lr_poly = LinearRegression()
        lr_poly.fit(X_poly, y)
        X_range_poly = poly.transform(x_range)
        y_pred_poly = lr_poly.predict(X_range_poly)
        ax.plot(x_range, y_pred_poly, 'b-', linewidth=2, label=f'Cubic (R²={row["R2_cubic"]:.3f})')

    # Labels and title
    ax.set_xlabel(predictor.replace('_', ' ').title(), fontsize=12)
    ax.set_ylabel(target.replace('_', ' ').title(), fontsize=12)
    title = f'{predictor} vs {target}\n'
    title += f'r={row["pearson_r"]:.3f} (p={row["pearson_p"]:.3f}), n={row["n"]}'
    ax.set_title(title, fontsize=13, fontweight='bold')
    ax.legend(loc='best')
    ax.grid(True, alpha=0.3)

    # Save
    safe_name = f"{predictor}_vs_{target}".replace('.', '_')
    plt.tight_layout()
    plt.savefig(f'{PLOT_DIR}/plot_{safe_name}.png', dpi=150, bbox_inches='tight')
    plt.close()

print(f"  Saved plots to: {PLOT_DIR}/")
print()

# ============================================================
# Step 6: Summary visualizations
# ============================================================

print("[6] Creating summary visualizations...")

# Histogram of driver_ratio across patients
fig, ax = plt.subplots(figsize=(10, 6))
ax.hist(df_patient_vdriver['driver_ratio'], bins=15, edgecolor='black', alpha=0.7)
ax.axvline(df_patient_vdriver['driver_ratio'].mean(), color='red', linestyle='--', linewidth=2,
           label=f'Mean = {df_patient_vdriver["driver_ratio"].mean():.4f}')
ax.set_xlabel('Vascular Driver Ratio', fontsize=12)
ax.set_ylabel('Number of Patients', fontsize=12)
ax.set_title('Distribution of Vascular Driver-State Abundance Across Patients', fontsize=13, fontweight='bold')
ax.legend()
ax.grid(True, alpha=0.3, axis='y')
plt.tight_layout()
plt.savefig(f'{OUTPUT_DIR}/patient_driver_ratio_histogram.png', dpi=150, bbox_inches='tight')
plt.close()

# Correlation heatmap
fig, ax = plt.subplots(figsize=(14, 10))

# Pivot correlation matrix
heatmap_data = df_correlations.pivot_table(
    index='target',
    columns='predictor',
    values='pearson_r'
)

# Create heatmap
sns.heatmap(heatmap_data, annot=True, fmt='.3f', cmap='RdBu_r', center=0,
            vmin=-1, vmax=1, linewidths=0.5, cbar_kws={'label': 'Pearson r'},
            ax=ax)
ax.set_title('Patient-level Correlations: Vascular Driver vs Pathology Signatures',
             fontsize=14, fontweight='bold', pad=20)
ax.set_xlabel('Predictor', fontsize=12)
ax.set_ylabel('Target (Pathology Signature)', fontsize=12)
plt.tight_layout()
plt.savefig(f'{OUTPUT_DIR}/patient_correlation_heatmap.png', dpi=150, bbox_inches='tight')
plt.close()

print(f"  Saved summary plots")
print()

# ============================================================
# Completion
# ============================================================

print("=" * 80)
print("Phase 9″ Complete!")
print("=" * 80)
print()
print("Key outputs:")
print(f"  - {OUTPUT_DIR}/patient_vascular_driver_abundance.csv")
print(f"  - {OUTPUT_DIR}/patient_pathology_signatures.csv")
print(f"  - {OUTPUT_DIR}/patient_vdriver_correlations.csv")
print(f"  - {OUTPUT_DIR}/patient_driver_ratio_histogram.png")
print(f"  - {OUTPUT_DIR}/patient_correlation_heatmap.png")
print(f"  - {PLOT_DIR}/plot_*.png ({len(interesting)} scatter plots)")
print()
print("Next: Generate comprehensive markdown report")
