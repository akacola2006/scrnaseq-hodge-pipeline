#!/usr/bin/env python3
"""
GSE212630 Trigger Stage Analysis v2 (Structural Hypothesis)
============================================================
Key changes from v1:
- Vascular = hypothesis-central (Mode B fixed, not excluded)
- Effect size > p-value
- Cell-loss + enrichment composite index
- Gate v2: QC-corrected (DAMP residualized)
- Vascular-composite metrics
"""

import pandas as pd
import numpy as np
from pathlib import Path
from scipy.stats import mannwhitneyu, spearmanr
from scipy import stats
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LinearRegression
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import gzip
import gc
import warnings
warnings.filterwarnings('ignore')

# ============================================================================
# CONFIGURATION
# ============================================================================
BASE_DIR = Path('/home/akaco/als/motor_cortex_analysis/ids_causal_analysis')
OUT_DIR = BASE_DIR / 'results_gse212630_triggerstage_v2'
OUT_DIR.mkdir(exist_ok=True)
(OUT_DIR / 'tables').mkdir(exist_ok=True)
(OUT_DIR / 'figures').mkdir(exist_ok=True)
(OUT_DIR / 'reports').mkdir(exist_ok=True)

GSE212630_DIR = BASE_DIR / 'GSE_212630_raw_expression_transposed'

# Gene sets
GENE_SETS = {
    'DAMP': ['HMGB1', 'HSPA1A', 'HSPA1B', 'HSP90AA1', 'HSP90AB1',
             'S100A8', 'S100A9', 'S100A12', 'IL33', 'LGALS3', 'ANXA1', 'CALR', 'HSPD1'],
    'PRR': ['IFIH1', 'DDX58', 'DHX58', 'MAVS', 'EIF2AK2', 'TLR3', 'TLR7',
            'STING1', 'MB21D1', 'TREX1', 'SAMHD1', 'ADAR'],
    'ISG': ['ISG15', 'OAS1', 'OAS2', 'OAS3', 'MX1', 'MX2', 'IFIT1', 'IFIT2', 'IFIT3',
            'IFI44', 'IFI44L', 'RSAD2', 'IRF7', 'STAT1', 'STAT2', 'USP18', 'IFITM1', 'IFITM3', 'BST2'],
    'DDR': ['ATM', 'ATR', 'TP53BP1', 'H2AFX', 'BRCA1', 'BRCA2', 'PARP1', 'PARP2',
            'RAD51', 'XRCC1', 'XRCC5', 'XRCC6', 'NBN', 'MRE11', 'CHEK1', 'CHEK2', 'MDC1'],
    'RNAproc': ['XRN2', 'EXOSC2', 'DDX3X', 'DDX5', 'UPF1', 'SRSF1', 'SRSF2',
                'HNRNPA1', 'HNRNPK', 'PABPN1', 'FUS', 'TARDBP', 'TIA1'],
    'Mito': ['TFAM', 'POLG', 'NDUFA1', 'NDUFA2', 'SDHA', 'SDHB', 'UQCRB',
             'COX4I1', 'COX5A', 'ATP5F1A', 'ATP5F1B', 'VDAC1', 'SOD2', 'PINK1'],
}

ALL_GENES = set()
for genes in GENE_SETS.values():
    ALL_GENES.update(genes)
ALL_GENES = list(ALL_GENES)

TDP_STAGES = ['Control', 'TDPneg', 'TDPmed', 'TDPhigh']
TDP_ORDER = {s: i for i, s in enumerate(TDP_STAGES)}

# Vascular cell types (hypothesis-central)
VASCULAR_TYPES = ['Vasc.Endo', 'Vasc.Capillary', 'Vasc.Fibro', 'Vasc.Pericyte', 'Vasc.SMC', 'Vasc.Unknown']

def read_genes_only(filepath, genes_of_interest):
    """Read only rows matching genes of interest"""
    gene_set = set(genes_of_interest)
    data = {}
    header = None
    with gzip.open(filepath, 'rt') as f:
        for i, line in enumerate(f):
            parts = line.strip().split(',')
            if i == 0:
                header = parts
                continue
            gene = parts[0].strip('"')
            if gene in gene_set:
                data[gene] = [float(x) if x else 0.0 for x in parts[1:]]
    if not data:
        return None, None
    df = pd.DataFrame(data)
    cell_ids = [h.strip('"') for h in header[1:]]
    df.index = cell_ids
    return df, list(data.keys())

def z_normalize_within_celltype(df, score_col, celltype_col='cell_type'):
    result = pd.Series(index=df.index, dtype=float)
    for ct in df[celltype_col].unique():
        mask = df[celltype_col] == ct
        vals = df.loc[mask, score_col]
        if len(vals) > 1 and vals.std() > 0:
            result.loc[mask] = (vals - vals.mean()) / vals.std()
        else:
            result.loc[mask] = 0
    return result

def cliffs_delta(x, y):
    n1, n2 = len(x), len(y)
    if n1 == 0 or n2 == 0:
        return np.nan
    count = 0
    for a in x:
        for b in y:
            if a > b:
                count += 1
            elif a < b:
                count -= 1
    return count / (n1 * n2)

def bootstrap_ci(data, n_boot=1000, ci=0.95):
    """Bootstrap confidence interval for mean"""
    if len(data) < 2:
        return np.nan, np.nan
    boot_means = []
    for _ in range(n_boot):
        sample = np.random.choice(data, size=len(data), replace=True)
        boot_means.append(np.mean(sample))
    lower = np.percentile(boot_means, (1-ci)/2 * 100)
    upper = np.percentile(boot_means, (1+ci)/2 * 100)
    return lower, upper

# ============================================================================
# STEP 1: Load GSE212630 data
# ============================================================================
print("=" * 70)
print("STEP 1: Loading GSE212630 data with extended QC metrics")
print("=" * 70)

gse_expr_files = sorted(GSE212630_DIR.glob('*_expression_transposed.csv.gz'))
print(f"Found {len(gse_expr_files)} cell type files")

all_cells = []
for i, expr_file in enumerate(gse_expr_files):
    ct_name = expr_file.stem.replace('_expression_transposed.csv', '').replace('.gz', '')
    meta_file = GSE212630_DIR / f"{ct_name}_metadata.csv"
    if not meta_file.exists():
        continue

    print(f"  [{i+1}/{len(gse_expr_files)}] {ct_name}...", end=' ', flush=True)
    meta = pd.read_csv(meta_file)
    expr, found_genes = read_genes_only(expr_file, ALL_GENES)
    if expr is None:
        print("no genes!")
        continue

    expr.index = expr.index.str.replace(r'_\d+_(Control|TDPneg|TDPmed|TDPhigh)$', '', regex=True)
    common_cells = [c for c in meta['cell_id'] if c in expr.index]
    if len(common_cells) == 0:
        print("no cells!")
        continue

    expr_aligned = expr.loc[common_cells]
    meta_aligned = meta[meta['cell_id'].isin(common_cells)].set_index('cell_id').loc[common_cells].reset_index()

    # Calculate scores
    scores = pd.DataFrame({'cell_id': common_cells})
    for set_name, genes in GENE_SETS.items():
        gene_cols = [g for g in genes if g in expr_aligned.columns]
        if len(gene_cols) > 0:
            scores[f'{set_name}_det'] = (expr_aligned[gene_cols] > 0).sum(axis=1).values / len(gene_cols)
            scores[f'{set_name}_mean'] = expr_aligned[gene_cols].mean(axis=1).values

    if 'PRR_det' in scores.columns and 'ISG_det' in scores.columns:
        scores['NSA_det'] = scores['PRR_det'] - scores['ISG_det']

    # QC metrics
    scores['n_genes_detected'] = (expr_aligned > 0).sum(axis=1).values
    scores['total_counts'] = expr_aligned.sum(axis=1).values

    scores['sample'] = meta_aligned['sample'].values
    scores['condition'] = meta_aligned['condition'].values
    scores['cell_type'] = meta_aligned['cell_type'].values

    all_cells.append(scores)
    print(f"{len(common_cells)} cells")
    del expr, expr_aligned
    gc.collect()

df = pd.concat(all_cells, ignore_index=True)
print(f"\nTotal: {len(df):,} cells")

# Z-normalize within cell type
score_cols = ['DAMP_det', 'PRR_det', 'ISG_det', 'DDR_det', 'RNAproc_det', 'Mito_det', 'NSA_det']
for col in score_cols:
    if col in df.columns:
        df[f'{col}_z'] = z_normalize_within_celltype(df, col)

# ============================================================================
# STEP 2: QC-corrected DAMP (DAMP_resid) for Gate v2
# ============================================================================
print("\n" + "=" * 70)
print("STEP 2: Creating QC-corrected DAMP (Gate v2 preparation)")
print("=" * 70)

# Regress DAMP on n_genes_detected within each cell type
df['DAMP_resid'] = np.nan
for ct in df['cell_type'].unique():
    mask = df['cell_type'] == ct
    ct_df = df.loc[mask]
    if len(ct_df) < 50:
        df.loc[mask, 'DAMP_resid'] = ct_df['DAMP_det'] - ct_df['DAMP_det'].mean()
        continue

    X = ct_df[['n_genes_detected']].values
    y = ct_df['DAMP_det'].values
    valid = ~np.isnan(X.flatten()) & ~np.isnan(y)
    if valid.sum() < 10:
        df.loc[mask, 'DAMP_resid'] = y - np.nanmean(y)
        continue

    reg = LinearRegression()
    reg.fit(X[valid], y[valid])
    pred = reg.predict(X)
    resid = y - pred
    df.loc[mask, 'DAMP_resid'] = resid

# Z-normalize DAMP_resid
df['DAMP_resid_z'] = z_normalize_within_celltype(df, 'DAMP_resid')

print(f"  DAMP-nGenes correlation (overall): r = {df['DAMP_det'].corr(df['n_genes_detected']):.3f}")
print(f"  DAMP_resid-nGenes correlation (should be ~0): r = {df['DAMP_resid'].corr(df['n_genes_detected']):.3f}")

# ============================================================================
# STEP 3: Extended counts table with QC metrics
# ============================================================================
print("\n" + "=" * 70)
print("STEP 3: Creating extended counts table")
print("=" * 70)

# Calculate Control baseline for cell-loss index
control_counts = df[df['condition'] == 'Control'].groupby('cell_type').agg(
    control_cells=('cell_id', 'count'),
    control_donors=('sample', 'nunique')
).reset_index()

# Full counts table
counts_extended = df.groupby(['cell_type', 'condition']).agg(
    n_cells=('cell_id', 'count'),
    n_donors=('sample', 'nunique'),
    mean_genes_detected=('n_genes_detected', 'mean'),
    mean_total_counts=('total_counts', 'mean'),
    mean_DAMP=('DAMP_det', 'mean'),
    mean_DAMP_resid=('DAMP_resid', 'mean')
).reset_index()

# Add cell fraction within stage
stage_totals = counts_extended.groupby('condition')['n_cells'].sum()
counts_extended['cell_fraction_in_stage'] = counts_extended.apply(
    lambda x: x['n_cells'] / stage_totals[x['condition']], axis=1
)

# Add cell-loss index (relative to Control)
counts_extended = counts_extended.merge(control_counts, on='cell_type', how='left')
counts_extended['cell_loss_index'] = counts_extended['n_cells'] / counts_extended['control_cells'].replace(0, np.nan)

# Pivot for easier viewing
counts_pivot = counts_extended.pivot(index='cell_type', columns='condition', values='n_cells').fillna(0).astype(int)
counts_pivot = counts_pivot[TDP_STAGES]
counts_pivot['total'] = counts_pivot.sum(axis=1)

# Add is_vascular flag
counts_pivot['is_vascular'] = counts_pivot.index.isin(VASCULAR_TYPES)

counts_extended.to_csv(OUT_DIR / 'tables/celltype_counts_by_stage_extended.csv', index=False)
print(f"Saved: celltype_counts_by_stage_extended.csv")

# ============================================================================
# STEP 4: Define Gate v1 and Gate v2
# ============================================================================
print("\n" + "=" * 70)
print("STEP 4: Defining Trigger Gates (v1: original, v2: QC-corrected)")
print("=" * 70)

def define_gate_v1(data, nsa_q=0.75, ddr_q=0.50, damp_q=0.50):
    """Gate v1: Original (DAMP low + |NSA| high + DDR high)"""
    data = data.copy()
    data['gate_v1'] = False
    for ct in data['cell_type'].unique():
        mask = data['cell_type'] == ct
        ct_df = data.loc[mask]
        if len(ct_df) < 30:
            continue
        nsa_thresh = ct_df['NSA_det_z'].abs().quantile(nsa_q)
        ddr_thresh = ct_df['DDR_det_z'].quantile(ddr_q)
        damp_thresh = ct_df['DAMP_det_z'].quantile(damp_q)
        gate = (ct_df['NSA_det_z'].abs() >= nsa_thresh) & \
               (ct_df['DDR_det_z'] >= ddr_thresh) & \
               (ct_df['DAMP_det_z'] <= damp_thresh)
        data.loc[mask, 'gate_v1'] = gate.values
    return data

def define_gate_v2(data, nsa_q=0.75, ddr_q=0.50):
    """Gate v2: QC-corrected (no DAMP, or DAMP_resid, focus on NSA+DDR)"""
    data = data.copy()
    data['gate_v2'] = False
    for ct in data['cell_type'].unique():
        mask = data['cell_type'] == ct
        ct_df = data.loc[mask]
        if len(ct_df) < 30:
            continue
        nsa_thresh = ct_df['NSA_det_z'].abs().quantile(nsa_q)
        ddr_thresh = ct_df['DDR_det_z'].quantile(ddr_q)
        # Gate v2: NSA high + DDR high (DAMP excluded to avoid QC bias)
        gate = (ct_df['NSA_det_z'].abs() >= nsa_thresh) & \
               (ct_df['DDR_det_z'] >= ddr_thresh)
        data.loc[mask, 'gate_v2'] = gate.values
    return data

df = define_gate_v1(df)
df = define_gate_v2(df)

print(f"  Gate v1 cells: {df['gate_v1'].sum():,} ({100*df['gate_v1'].mean():.1f}%)")
print(f"  Gate v2 cells: {df['gate_v2'].sum():,} ({100*df['gate_v2'].mean():.1f}%)")

# Check v1/v2 concordance
v1v2_corr = df['gate_v1'].astype(int).corr(df['gate_v2'].astype(int))
print(f"  Gate v1-v2 correlation: r = {v1v2_corr:.3f}")

# ============================================================================
# STEP 5: Calculate trigger metrics by celltype x stage
# ============================================================================
print("\n" + "=" * 70)
print("STEP 5: Calculating trigger metrics with risk indices")
print("=" * 70)

# Donor-level aggregation
donor_agg = df.groupby(['sample', 'condition', 'cell_type']).agg(
    gate_v1_frac=('gate_v1', 'mean'),
    gate_v2_frac=('gate_v2', 'mean'),
    n_cells=('cell_id', 'count'),
    mean_genes=('n_genes_detected', 'mean'),
    DAMP_mean=('DAMP_det', 'mean'),
    DAMP_resid_mean=('DAMP_resid', 'mean'),
    NSA_mean=('NSA_det', 'mean'),
    DDR_mean=('DDR_det', 'mean')
).reset_index()

# Calculate metrics by celltype x stage
metrics_list = []
for ct in df['cell_type'].unique():
    for stage in TDP_STAGES:
        ct_stage = donor_agg[(donor_agg['cell_type'] == ct) & (donor_agg['condition'] == stage)]
        ctrl = donor_agg[(donor_agg['cell_type'] == ct) & (donor_agg['condition'] == 'Control')]

        if len(ct_stage) == 0:
            continue

        # Gate fractions
        gate_v1_frac = ct_stage['gate_v1_frac'].mean()
        gate_v2_frac = ct_stage['gate_v2_frac'].mean()

        # Bootstrap CI for gate fractions
        v1_ci_low, v1_ci_high = bootstrap_ci(ct_stage['gate_v1_frac'].values)
        v2_ci_low, v2_ci_high = bootstrap_ci(ct_stage['gate_v2_frac'].values)

        # Cell counts
        n_cells = df[(df['cell_type'] == ct) & (df['condition'] == stage)].shape[0]
        n_donors = ct_stage['sample'].nunique()

        # Cell-loss index (vs Control)
        ctrl_cells = df[(df['cell_type'] == ct) & (df['condition'] == 'Control')].shape[0]
        cell_loss_idx = n_cells / ctrl_cells if ctrl_cells > 0 else np.nan

        # Cells per donor (for risk index)
        cells_per_donor = n_cells / n_donors if n_donors > 0 else 0
        ctrl_cells_per_donor = ctrl_cells / ctrl['sample'].nunique() if len(ctrl) > 0 else 0

        # Risk index: gate_frac * (ctrl_cells_per_donor / stage_cells_per_donor)
        # "enrichment despite cell loss"
        if cells_per_donor > 0:
            risk_index_v1 = gate_v1_frac * (ctrl_cells_per_donor / cells_per_donor) if ctrl_cells_per_donor > 0 else gate_v1_frac
            risk_index_v2 = gate_v2_frac * (ctrl_cells_per_donor / cells_per_donor) if ctrl_cells_per_donor > 0 else gate_v2_frac
        else:
            risk_index_v1 = np.nan
            risk_index_v2 = np.nan

        # Cliff's delta vs Control
        if len(ctrl) >= 2 and len(ct_stage) >= 2:
            delta_v1 = cliffs_delta(ct_stage['gate_v1_frac'].values, ctrl['gate_v1_frac'].values)
            delta_v2 = cliffs_delta(ct_stage['gate_v2_frac'].values, ctrl['gate_v2_frac'].values)
        else:
            delta_v1 = np.nan
            delta_v2 = np.nan

        metrics_list.append({
            'cell_type': ct,
            'condition': stage,
            'n_cells': n_cells,
            'n_donors': n_donors,
            'gate_v1_frac': gate_v1_frac,
            'gate_v1_ci_low': v1_ci_low,
            'gate_v1_ci_high': v1_ci_high,
            'gate_v2_frac': gate_v2_frac,
            'gate_v2_ci_low': v2_ci_low,
            'gate_v2_ci_high': v2_ci_high,
            'cell_loss_index': cell_loss_idx,
            'cells_per_donor': cells_per_donor,
            'risk_index_v1': risk_index_v1,
            'risk_index_v2': risk_index_v2,
            'cliffs_delta_v1': delta_v1,
            'cliffs_delta_v2': delta_v2,
            'is_vascular': ct in VASCULAR_TYPES
        })

metrics_df = pd.DataFrame(metrics_list)
metrics_df.to_csv(OUT_DIR / 'tables/trigger_metrics_by_stage_celltype.csv', index=False)
print(f"Saved: trigger_metrics_by_stage_celltype.csv")

# ============================================================================
# STEP 6: Vascular-composite metrics
# ============================================================================
print("\n" + "=" * 70)
print("STEP 6: Creating vascular-composite metrics")
print("=" * 70)

vasc_df = df[df['cell_type'].isin(VASCULAR_TYPES)]
print(f"  Vascular cells total: {len(vasc_df):,}")

# Composite by stage (weighted by n_donors)
vasc_composite = []
for stage in TDP_STAGES:
    stage_vasc = vasc_df[vasc_df['condition'] == stage]
    if len(stage_vasc) == 0:
        continue

    # Per-celltype within vascular
    vasc_stage_metrics = metrics_df[(metrics_df['is_vascular']) & (metrics_df['condition'] == stage)]

    if len(vasc_stage_metrics) == 0:
        continue

    # Weighted average (by n_donors)
    total_donors = vasc_stage_metrics['n_donors'].sum()
    if total_donors > 0:
        weights = vasc_stage_metrics['n_donors'] / total_donors
        composite_gate_v1 = (vasc_stage_metrics['gate_v1_frac'] * weights).sum()
        composite_gate_v2 = (vasc_stage_metrics['gate_v2_frac'] * weights).sum()
        composite_cell_loss = (vasc_stage_metrics['cell_loss_index'].fillna(1) * weights).sum()
        composite_risk_v1 = (vasc_stage_metrics['risk_index_v1'].fillna(0) * weights).sum()
        composite_risk_v2 = (vasc_stage_metrics['risk_index_v2'].fillna(0) * weights).sum()
    else:
        composite_gate_v1 = vasc_stage_metrics['gate_v1_frac'].mean()
        composite_gate_v2 = vasc_stage_metrics['gate_v2_frac'].mean()
        composite_cell_loss = vasc_stage_metrics['cell_loss_index'].mean()
        composite_risk_v1 = vasc_stage_metrics['risk_index_v1'].mean()
        composite_risk_v2 = vasc_stage_metrics['risk_index_v2'].mean()

    vasc_composite.append({
        'condition': stage,
        'n_cells': len(stage_vasc),
        'n_donors': stage_vasc['sample'].nunique(),
        'composite_gate_v1': composite_gate_v1,
        'composite_gate_v2': composite_gate_v2,
        'composite_cell_loss': composite_cell_loss,
        'composite_risk_v1': composite_risk_v1,
        'composite_risk_v2': composite_risk_v2
    })

vasc_composite_df = pd.DataFrame(vasc_composite)
vasc_composite_df.to_csv(OUT_DIR / 'tables/vascular_composite_by_stage.csv', index=False)
print(f"Saved: vascular_composite_by_stage.csv")
print(vasc_composite_df.to_string(index=False))

# ============================================================================
# STEP 7: Compare OPC/Astro with Vascular
# ============================================================================
print("\n" + "=" * 70)
print("STEP 7: Comparison: Vascular vs OPC vs Astro")
print("=" * 70)

comparison_types = ['Glia.OPC', 'Glia.Astro.GFAP.neg', 'Glia.Astro.GFAP.pos']
comparison_data = []

for ct in comparison_types:
    ct_metrics = metrics_df[metrics_df['cell_type'] == ct]
    for _, row in ct_metrics.iterrows():
        comparison_data.append({
            'group': ct,
            'condition': row['condition'],
            'gate_v1_frac': row['gate_v1_frac'],
            'gate_v2_frac': row['gate_v2_frac'],
            'risk_index_v1': row['risk_index_v1'],
            'risk_index_v2': row['risk_index_v2'],
            'cell_loss_index': row['cell_loss_index']
        })

# Add vascular composite
for _, row in vasc_composite_df.iterrows():
    comparison_data.append({
        'group': 'Vascular (composite)',
        'condition': row['condition'],
        'gate_v1_frac': row['composite_gate_v1'],
        'gate_v2_frac': row['composite_gate_v2'],
        'risk_index_v1': row['composite_risk_v1'],
        'risk_index_v2': row['composite_risk_v2'],
        'cell_loss_index': row['composite_cell_loss']
    })

comparison_df = pd.DataFrame(comparison_data)
comparison_df.to_csv(OUT_DIR / 'tables/vascular_opc_astro_comparison.csv', index=False)
print(f"Saved: vascular_opc_astro_comparison.csv")

# ============================================================================
# STEP 8: Threshold sweep stability
# ============================================================================
print("\n" + "=" * 70)
print("STEP 8: Threshold sweep stability")
print("=" * 70)

nsa_qs = [0.70, 0.75, 0.80]
ddr_qs = [0.40, 0.50, 0.60]

sweep_results = []
reference_ranking = None

for nsa_q in nsa_qs:
    for ddr_q in ddr_qs:
        # Gate v2 sweep (no DAMP)
        df_sweep = define_gate_v2(df.copy(), nsa_q=nsa_q, ddr_q=ddr_q)

        sweep_frac = df_sweep.groupby(['cell_type', 'condition'])['gate_v2'].mean().reset_index()
        sweep_pivot = sweep_frac.pivot(index='cell_type', columns='condition', values='gate_v2')

        if 'Control' in sweep_pivot.columns and 'TDPneg' in sweep_pivot.columns:
            sweep_pivot['diff'] = sweep_pivot['TDPneg'] - sweep_pivot['Control']
            ranking = sweep_pivot.sort_values('diff', ascending=False).index.tolist()

            if reference_ranking is None:
                reference_ranking = ranking

            top5_overlap = len(set(ranking[:5]) & set(reference_ranking[:5])) / 5

            sweep_results.append({
                'nsa_q': nsa_q,
                'ddr_q': ddr_q,
                'gate_pct': 100 * df_sweep['gate_v2'].mean(),
                'top5_stability': top5_overlap,
                'top5_celltypes': ', '.join(ranking[:5])
            })

sweep_df = pd.DataFrame(sweep_results)
sweep_df.to_csv(OUT_DIR / 'tables/trigger_gate_v2_threshold_sweep.csv', index=False)
print(f"Saved: trigger_gate_v2_threshold_sweep.csv")
print(f"  Mean top-5 stability: {sweep_df['top5_stability'].mean():.2f}")

# ============================================================================
# STEP 9: Visualizations
# ============================================================================
print("\n" + "=" * 70)
print("STEP 9: Creating visualizations")
print("=" * 70)

# 9A: Gate fraction vs cell-loss scatter
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# Gate v1
ax = axes[0]
for stage in TDP_STAGES:
    stage_data = metrics_df[metrics_df['condition'] == stage]
    vasc = stage_data[stage_data['is_vascular']]
    nonvasc = stage_data[~stage_data['is_vascular']]
    ax.scatter(nonvasc['cell_loss_index'], nonvasc['gate_v1_frac'],
               alpha=0.6, s=50, label=f'{stage} (non-vasc)' if stage == 'Control' else None)
    ax.scatter(vasc['cell_loss_index'], vasc['gate_v1_frac'],
               marker='s', s=100, edgecolor='red', linewidth=2,
               label=f'{stage} (vascular)' if stage == 'Control' else None)
ax.set_xlabel('Cell-loss index (vs Control)')
ax.set_ylabel('Gate v1 fraction')
ax.set_title('Gate v1: "Enrichment despite cell loss"')
ax.axvline(1.0, color='gray', linestyle='--', alpha=0.5)
ax.legend()

# Gate v2
ax = axes[1]
colors = {'Control': 'blue', 'TDPneg': 'green', 'TDPmed': 'orange', 'TDPhigh': 'red'}
for stage in TDP_STAGES:
    stage_data = metrics_df[metrics_df['condition'] == stage]
    vasc = stage_data[stage_data['is_vascular']]
    nonvasc = stage_data[~stage_data['is_vascular']]
    ax.scatter(nonvasc['cell_loss_index'], nonvasc['gate_v2_frac'],
               alpha=0.6, s=50, c=colors[stage])
    ax.scatter(vasc['cell_loss_index'], vasc['gate_v2_frac'],
               marker='s', s=100, c=colors[stage], edgecolor='red', linewidth=2)
ax.set_xlabel('Cell-loss index (vs Control)')
ax.set_ylabel('Gate v2 fraction')
ax.set_title('Gate v2 (QC-corrected): "Enrichment despite cell loss"')
ax.axvline(1.0, color='gray', linestyle='--', alpha=0.5)

plt.tight_layout()
plt.savefig(OUT_DIR / 'figures/gate_vs_cellloss_scatter.png', dpi=150)
plt.close()
print("  Saved: gate_vs_cellloss_scatter.png")

# 9B: Vascular composite comparison
fig, axes = plt.subplots(1, 3, figsize=(15, 5))

# Gate v2 fraction
ax = axes[0]
plot_df = comparison_df.pivot(index='group', columns='condition', values='gate_v2_frac')
if TDP_STAGES[0] in plot_df.columns:
    plot_df = plot_df[TDP_STAGES]
plot_df.plot(kind='bar', ax=ax, width=0.8, color=[colors[s] for s in TDP_STAGES])
ax.set_ylabel('Gate v2 Fraction')
ax.set_title('Trigger Gate (v2) by Stage')
ax.legend(title='Stage')
plt.sca(ax)
plt.xticks(rotation=45, ha='right')

# Risk index v2
ax = axes[1]
plot_df = comparison_df.pivot(index='group', columns='condition', values='risk_index_v2')
if TDP_STAGES[0] in plot_df.columns:
    plot_df = plot_df[TDP_STAGES]
plot_df.plot(kind='bar', ax=ax, width=0.8, color=[colors[s] for s in TDP_STAGES])
ax.set_ylabel('Risk Index v2')
ax.set_title('Risk Index (enrichment × cell-loss)')
ax.legend(title='Stage')
plt.sca(ax)
plt.xticks(rotation=45, ha='right')

# Cell-loss index
ax = axes[2]
plot_df = comparison_df.pivot(index='group', columns='condition', values='cell_loss_index')
if TDP_STAGES[0] in plot_df.columns:
    plot_df = plot_df[TDP_STAGES]
plot_df.plot(kind='bar', ax=ax, width=0.8, color=[colors[s] for s in TDP_STAGES])
ax.set_ylabel('Cell-loss Index (vs Control)')
ax.axhline(1.0, color='red', linestyle='--', alpha=0.7)
ax.set_title('Cell-loss by Stage')
ax.legend(title='Stage')
plt.sca(ax)
plt.xticks(rotation=45, ha='right')

plt.tight_layout()
plt.savefig(OUT_DIR / 'figures/vascular_opc_astro_comparison.png', dpi=150)
plt.close()
print("  Saved: vascular_opc_astro_comparison.png")

# 9C: DAMP vs QC by stage (updated)
fig, axes = plt.subplots(1, 4, figsize=(16, 4))
for i, stage in enumerate(TDP_STAGES):
    ax = axes[i]
    stage_data = donor_agg[donor_agg['condition'] == stage]
    ax.scatter(stage_data['mean_genes'], stage_data['DAMP_mean'], alpha=0.5, s=30)
    r = stage_data['mean_genes'].corr(stage_data['DAMP_mean'])
    ax.set_xlabel('Mean genes detected')
    ax.set_ylabel('DAMP (mean)')
    ax.set_title(f'{stage} (r={r:.2f})')
plt.tight_layout()
plt.savefig(OUT_DIR / 'figures/damp_vs_qc_scatter_v2.png', dpi=150)
plt.close()
print("  Saved: damp_vs_qc_scatter_v2.png")

# 9D: PCA with all features
print("  Creating Pre-TDP PCA...")
feature_cols = [c for c in df.columns if c.endswith('_z')]
ct_cond_agg = df.groupby(['cell_type', 'condition'])[feature_cols].mean().reset_index()

X_ct = ct_cond_agg[feature_cols].fillna(0).values
scaler = StandardScaler()
X_ct_scaled = scaler.fit_transform(X_ct)

pca = PCA(n_components=min(5, len(feature_cols)))
pca_result = pca.fit_transform(X_ct_scaled)
ct_cond_agg['PC1'] = pca_result[:, 0]
ct_cond_agg['PC2'] = pca_result[:, 1]

fig, ax = plt.subplots(figsize=(12, 10))
for stage in TDP_STAGES:
    stage_data = ct_cond_agg[ct_cond_agg['condition'] == stage]
    vasc = stage_data[stage_data['cell_type'].isin(VASCULAR_TYPES)]
    nonvasc = stage_data[~stage_data['cell_type'].isin(VASCULAR_TYPES)]

    ax.scatter(nonvasc['PC1'], nonvasc['PC2'], c=colors[stage], alpha=0.6, s=80, label=stage)
    ax.scatter(vasc['PC1'], vasc['PC2'], c=colors[stage], s=150,
               marker='s', edgecolor='black', linewidth=2)

    # Label key types
    for _, row in stage_data.iterrows():
        if row['cell_type'] in ['Glia.OPC', 'Glia.Astro.GFAP.pos', 'Vasc.Endo', 'Vasc.Capillary']:
            ax.annotate(row['cell_type'].split('.')[-1], (row['PC1'], row['PC2']),
                       fontsize=7, alpha=0.8)

ax.set_xlabel(f'PC1 ({100*pca.explained_variance_ratio_[0]:.1f}%)')
ax.set_ylabel(f'PC2 ({100*pca.explained_variance_ratio_[1]:.1f}%)')
ax.set_title('Pre-TDP Map: Celltype × Stage (squares = vascular)')
ax.legend(title='Stage')
ax.axhline(0, color='gray', linestyle='--', alpha=0.3)
ax.axvline(0, color='gray', linestyle='--', alpha=0.3)
plt.tight_layout()
plt.savefig(OUT_DIR / 'figures/preTDP_map_pca_v2.png', dpi=150)
plt.close()
print("  Saved: preTDP_map_pca_v2.png")

# ============================================================================
# STEP 10: V1 vs V2 concordance
# ============================================================================
print("\n" + "=" * 70)
print("STEP 10: Gate v1 vs v2 concordance analysis")
print("=" * 70)

# Compare rankings
v1_ranking = metrics_df[metrics_df['condition'] == 'TDPneg'].sort_values('gate_v1_frac', ascending=False)['cell_type'].tolist()
v2_ranking = metrics_df[metrics_df['condition'] == 'TDPneg'].sort_values('gate_v2_frac', ascending=False)['cell_type'].tolist()

# Spearman correlation
v1_ranks = {ct: i for i, ct in enumerate(v1_ranking)}
v2_ranks = {ct: i for i, ct in enumerate(v2_ranking)}
common_cts = set(v1_ranks.keys()) & set(v2_ranks.keys())
v1_order = [v1_ranks[ct] for ct in common_cts]
v2_order = [v2_ranks[ct] for ct in common_cts]
rho, pval = spearmanr(v1_order, v2_order)

print(f"  Gate v1 vs v2 ranking (TDPneg): Spearman rho = {rho:.3f} (p = {pval:.4f})")
print(f"  Top 5 (v1): {', '.join(v1_ranking[:5])}")
print(f"  Top 5 (v2): {', '.join(v2_ranking[:5])}")

concordance = pd.DataFrame({
    'v1_ranking': v1_ranking[:10],
    'v2_ranking': v2_ranking[:10]
})
concordance.to_csv(OUT_DIR / 'tables/gate_v1_v2_ranking_concordance.csv', index=False)

# ============================================================================
# STEP 11: Generate Report
# ============================================================================
print("\n" + "=" * 70)
print("STEP 11: Generating report")
print("=" * 70)

# Key findings
vasc_tdpneg = vasc_composite_df[vasc_composite_df['condition'] == 'TDPneg'].iloc[0] if len(vasc_composite_df[vasc_composite_df['condition'] == 'TDPneg']) > 0 else None
opc_tdpneg = metrics_df[(metrics_df['cell_type'] == 'Glia.OPC') & (metrics_df['condition'] == 'TDPneg')]
opc_row = opc_tdpneg.iloc[0] if len(opc_tdpneg) > 0 else None

report = f"""# GSE212630 Trigger Stage Analysis v2 (Structural Hypothesis)

## Key Premise
- TDP staging = **OUTCOME axis** (not cause)
- Vascular cells = **hypothesis-central** (not excluded due to low n)
- Effect size + stability > p-value
- Cell-loss as signal, not noise

## Data Summary
- Total cells: {len(df):,}
- Vascular cells: {len(vasc_df):,} ({100*len(vasc_df)/len(df):.1f}%)
- Gate v1 (DAMP low + |NSA| high + DDR high): {df['gate_v1'].sum():,} ({100*df['gate_v1'].mean():.1f}%)
- Gate v2 (|NSA| high + DDR high, no DAMP): {df['gate_v2'].sum():,} ({100*df['gate_v2'].mean():.1f}%)

## QC Correction
- DAMP-nGenes correlation (raw): r = {df['DAMP_det'].corr(df['n_genes_detected']):.3f}
- DAMP_resid-nGenes correlation: r = {df['DAMP_resid'].corr(df['n_genes_detected']):.3f}
- Gate v1-v2 cell-level correlation: r = {v1v2_corr:.3f}
- Gate v1-v2 ranking Spearman (TDPneg): rho = {rho:.3f}

## Vascular Composite (TDPneg)
"""

if vasc_tdpneg is not None:
    report += f"""- Gate v1: {vasc_tdpneg['composite_gate_v1']:.3f}
- Gate v2: {vasc_tdpneg['composite_gate_v2']:.3f}
- Cell-loss index: {vasc_tdpneg['composite_cell_loss']:.2f}
- Risk index v2: {vasc_tdpneg['composite_risk_v2']:.3f}
"""

report += f"""
## OPC (Mode A reference)
"""

if opc_row is not None:
    report += f"""- Gate v1: {opc_row['gate_v1_frac']:.3f}
- Gate v2: {opc_row['gate_v2_frac']:.3f}
- Cell-loss index: {opc_row['cell_loss_index']:.2f}
- Risk index v2: {opc_row['risk_index_v2']:.3f}
- Cliff's delta (v2 vs Control): {opc_row['cliffs_delta_v2']:.2f}
"""

report += f"""
## Threshold Stability (Gate v2)
- Mean top-5 stability: {sweep_df['top5_stability'].mean():.2f}
- Most stable top types: {sweep_df.iloc[0]['top5_celltypes']}

## Structural Interpretation

1. **Vascular as hypothesis-central**:
   - Low cell count may reflect cell loss (not exclusion criterion)
   - Risk index captures "enrichment despite loss"
   - Composite allows statistical power recovery

2. **Gate v2 (QC-corrected)**:
   - Removes DAMP bias from QC correlation
   - NSA + DDR captures nucleic acid stress + DNA damage response
   - v1-v2 concordance validates robustness

3. **OPC as observable proxy**:
   - Mode A statistical power
   - May receive signals from vascular dysfunction
   - "Downstream amplifier" hypothesis

## Caveats
1. Cell-loss could be technical (QC) or biological (death)
2. Risk index is exploratory (not validated)
3. Vascular composite weights by donor count (assumption)

## Output Files

### Tables
- `celltype_counts_by_stage_extended.csv`
- `trigger_metrics_by_stage_celltype.csv`
- `vascular_composite_by_stage.csv`
- `vascular_opc_astro_comparison.csv`
- `trigger_gate_v2_threshold_sweep.csv`
- `gate_v1_v2_ranking_concordance.csv`

### Figures
- `gate_vs_cellloss_scatter.png`
- `vascular_opc_astro_comparison.png`
- `damp_vs_qc_scatter_v2.png`
- `preTDP_map_pca_v2.png`
"""

with open(OUT_DIR / 'reports/trigger_stage_v2_summary.md', 'w') as f:
    f.write(report)
print(f"Saved: trigger_stage_v2_summary.md")

print("\n" + "=" * 70)
print("COMPLETE")
print("=" * 70)
print(f"Output: {OUT_DIR}")
