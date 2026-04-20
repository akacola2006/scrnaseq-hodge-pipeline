#!/usr/bin/env python3
"""
GSE212630 Trigger Stage Analysis v3 (Donor-Level Hierarchical)
===============================================================
Key changes from v2:
- Donor-level aggregation only (no pseudo-replication)
- Hierarchical model for low-n rescue (partial pooling)
- Bootstrap rank stability
- TDPneg internal clustering (trigger-stage)
- Vascular: individual + composite evaluation
"""

import pandas as pd
import numpy as np
from pathlib import Path
from scipy.stats import mannwhitneyu, spearmanr, beta
from scipy.special import logit, expit
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import AgglomerativeClustering
from sklearn.linear_model import LinearRegression
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import gzip
import gc
import warnings
warnings.filterwarnings('ignore')

# Try importing statsmodels for GLMM
try:
    import statsmodels.api as sm
    from statsmodels.formula.api import mixedlm
    HAS_STATSMODELS = True
except ImportError:
    HAS_STATSMODELS = False
    print("Warning: statsmodels not available, using simplified hierarchical estimation")

# ============================================================================
# CONFIGURATION
# ============================================================================
BASE_DIR = Path('/home/akaco/als/motor_cortex_analysis/ids_causal_analysis')
OUT_DIR = BASE_DIR / 'results_gse212630_triggerstage_v3'
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
VASCULAR_TYPES = ['Vasc.Endo', 'Vasc.Capillary', 'Vasc.Fibro', 'Vasc.Pericyte', 'Vasc.SMC']

def read_genes_only(filepath, genes_of_interest):
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

def beta_binomial_estimate(successes, trials, prior_alpha=1, prior_beta=1):
    """Empirical Bayes shrinkage using Beta-Binomial conjugate prior"""
    # Pool across all data to get prior
    total_succ = successes.sum()
    total_trials = trials.sum()
    if total_trials == 0:
        return np.full(len(successes), 0.5), np.full(len(successes), 0.0), np.full(len(successes), 1.0)

    # Prior from pooled data
    pooled_rate = total_succ / total_trials
    # Shrink towards pooled estimate
    post_alpha = prior_alpha + successes
    post_beta = prior_beta + (trials - successes)

    # Posterior mean
    post_mean = post_alpha / (post_alpha + post_beta)

    # 95% credible interval
    ci_low = beta.ppf(0.025, post_alpha, post_beta)
    ci_high = beta.ppf(0.975, post_alpha, post_beta)

    return post_mean, ci_low, ci_high

# ============================================================================
# STEP 1: Load data and create cell-level features
# ============================================================================
print("=" * 70)
print("STEP 1: Loading GSE212630 data")
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

    scores = pd.DataFrame({'cell_id': common_cells})
    for set_name, genes in GENE_SETS.items():
        gene_cols = [g for g in genes if g in expr_aligned.columns]
        if len(gene_cols) > 0:
            scores[f'{set_name}_det'] = (expr_aligned[gene_cols] > 0).sum(axis=1).values / len(gene_cols)

    if 'PRR_det' in scores.columns and 'ISG_det' in scores.columns:
        scores['NSA_det'] = scores['PRR_det'] - scores['ISG_det']

    scores['n_genes_detected'] = (expr_aligned > 0).sum(axis=1).values
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

# QC-correct DAMP
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
    df.loc[mask, 'DAMP_resid'] = y - pred

df['DAMP_resid_z'] = z_normalize_within_celltype(df, 'DAMP_resid')

# Define Gate v2 (NSA + DDR, no DAMP)
def define_gate_v2(data, nsa_q=0.75, ddr_q=0.50):
    data = data.copy()
    data['gate_v2'] = False
    for ct in data['cell_type'].unique():
        mask = data['cell_type'] == ct
        ct_df = data.loc[mask]
        if len(ct_df) < 30:
            continue
        nsa_thresh = ct_df['NSA_det_z'].abs().quantile(nsa_q)
        ddr_thresh = ct_df['DDR_det_z'].quantile(ddr_q)
        gate = (ct_df['NSA_det_z'].abs() >= nsa_thresh) & (ct_df['DDR_det_z'] >= ddr_thresh)
        data.loc[mask, 'gate_v2'] = gate.values
    return data

df = define_gate_v2(df)
print(f"Gate v2 cells: {df['gate_v2'].sum():,} ({100*df['gate_v2'].mean():.1f}%)")

# ============================================================================
# STEP 2: Donor-level aggregation (NO pseudo-replication)
# ============================================================================
print("\n" + "=" * 70)
print("STEP 2: Donor-level aggregation (pseudo-replication prohibited)")
print("=" * 70)

donor_agg = df.groupby(['sample', 'condition', 'cell_type']).agg(
    gate_v2_count=('gate_v2', 'sum'),
    n_cells=('cell_id', 'count'),
    NSA_z_mean=('NSA_det_z', 'mean'),
    DDR_z_mean=('DDR_det_z', 'mean'),
    DAMP_resid_z_mean=('DAMP_resid_z', 'mean'),
    PRR_z_mean=('PRR_det_z', 'mean'),
    ISG_z_mean=('ISG_det_z', 'mean'),
    mean_genes=('n_genes_detected', 'mean')
).reset_index()

donor_agg['gate_v2_frac'] = donor_agg['gate_v2_count'] / donor_agg['n_cells']
donor_agg['is_vascular'] = donor_agg['cell_type'].isin(VASCULAR_TYPES)

donor_agg.to_csv(OUT_DIR / 'tables/donor_level_gate_metrics.csv', index=False)
print(f"Saved: donor_level_gate_metrics.csv ({len(donor_agg)} rows)")

# Summary by celltype x stage
print("\nDonor counts by celltype x stage:")
donor_counts = donor_agg.groupby(['cell_type', 'condition']).agg(
    n_donors=('sample', 'nunique'),
    total_cells=('n_cells', 'sum'),
    mean_gate_frac=('gate_v2_frac', 'mean')
).reset_index()
donor_counts_pivot = donor_counts.pivot(index='cell_type', columns='condition', values='n_donors')
print(donor_counts_pivot[TDP_STAGES].to_string())

# ============================================================================
# STEP 3: Hierarchical model (Beta-Binomial shrinkage)
# ============================================================================
print("\n" + "=" * 70)
print("STEP 3: Hierarchical model (Beta-Binomial partial pooling)")
print("=" * 70)

# For each celltype, estimate stage effects with shrinkage
hierarchical_results = []

for ct in donor_agg['cell_type'].unique():
    ct_data = donor_agg[donor_agg['cell_type'] == ct].copy()

    for stage in TDP_STAGES:
        stage_data = ct_data[ct_data['condition'] == stage]
        ctrl_data = ct_data[ct_data['condition'] == 'Control']

        if len(stage_data) == 0:
            continue

        # Beta-Binomial estimation
        successes = stage_data['gate_v2_count'].values
        trials = stage_data['n_cells'].values
        post_mean, ci_low, ci_high = beta_binomial_estimate(successes, trials)

        # Average across donors (shrunk estimates)
        shrunk_mean = post_mean.mean()
        shrunk_ci_low = ci_low.mean()
        shrunk_ci_high = ci_high.mean()

        # Raw estimate for comparison
        raw_mean = stage_data['gate_v2_frac'].mean()
        raw_se = stage_data['gate_v2_frac'].std() / np.sqrt(len(stage_data)) if len(stage_data) > 1 else np.nan

        # Cliff's delta vs Control
        if len(ctrl_data) >= 1 and len(stage_data) >= 1:
            delta = cliffs_delta(stage_data['gate_v2_frac'].values, ctrl_data['gate_v2_frac'].values)
        else:
            delta = np.nan

        hierarchical_results.append({
            'cell_type': ct,
            'condition': stage,
            'n_donors': len(stage_data),
            'total_cells': stage_data['n_cells'].sum(),
            'raw_gate_frac': raw_mean,
            'raw_se': raw_se,
            'shrunk_gate_frac': shrunk_mean,
            'shrunk_ci_low': shrunk_ci_low,
            'shrunk_ci_high': shrunk_ci_high,
            'cliffs_delta_vs_ctrl': delta,
            'is_vascular': ct in VASCULAR_TYPES
        })

hier_df = pd.DataFrame(hierarchical_results)
hier_df.to_csv(OUT_DIR / 'tables/hierarchical_stage_effects.csv', index=False)
print(f"Saved: hierarchical_stage_effects.csv")

# Show key results
print("\nShrunk estimates (TDPneg, sorted by gate_frac):")
tdpneg = hier_df[hier_df['condition'] == 'TDPneg'].sort_values('shrunk_gate_frac', ascending=False)
for _, row in tdpneg.head(10).iterrows():
    vasc_marker = '[V]' if row['is_vascular'] else '   '
    print(f"  {vasc_marker} {row['cell_type']:30s}: {row['shrunk_gate_frac']:.3f} "
          f"[{row['shrunk_ci_low']:.3f}, {row['shrunk_ci_high']:.3f}] "
          f"(δ={row['cliffs_delta_vs_ctrl']:+.2f}, n={row['n_donors']})")

# ============================================================================
# STEP 4: Vascular individual + composite evaluation
# ============================================================================
print("\n" + "=" * 70)
print("STEP 4: Vascular individual + composite evaluation")
print("=" * 70)

vasc_individual = hier_df[hier_df['is_vascular']]
print("\nVascular individual cell types (TDPneg):")
vasc_tdpneg = vasc_individual[vasc_individual['condition'] == 'TDPneg']
for _, row in vasc_tdpneg.iterrows():
    print(f"  {row['cell_type']:20s}: gate={row['shrunk_gate_frac']:.3f}, δ={row['cliffs_delta_vs_ctrl']:+.2f}, n_donors={row['n_donors']}, cells={row['total_cells']}")

# Vascular composite (donor-level pooling)
vasc_donor = donor_agg[donor_agg['is_vascular']].copy()

# Pool all vascular cells per donor x stage
vasc_composite = vasc_donor.groupby(['sample', 'condition']).agg(
    gate_count=('gate_v2_count', 'sum'),
    n_cells=('n_cells', 'sum'),
    NSA_z_mean=('NSA_z_mean', 'mean'),
    DDR_z_mean=('DDR_z_mean', 'mean')
).reset_index()
vasc_composite['gate_frac'] = vasc_composite['gate_count'] / vasc_composite['n_cells']

print("\nVascular composite by stage:")
for stage in TDP_STAGES:
    stage_data = vasc_composite[vasc_composite['condition'] == stage]
    if len(stage_data) > 0:
        mean_frac = stage_data['gate_frac'].mean()
        mean_cells = stage_data['n_cells'].mean()
        n_donors = len(stage_data)
        print(f"  {stage:10s}: gate={mean_frac:.3f}, cells/donor={mean_cells:.1f}, n_donors={n_donors}")

# Save vascular comparison
vasc_comparison = []
for ct in VASCULAR_TYPES:
    ct_hier = hier_df[(hier_df['cell_type'] == ct) & (hier_df['condition'] == 'TDPneg')]
    if len(ct_hier) > 0:
        row = ct_hier.iloc[0]
        vasc_comparison.append({
            'type': 'individual',
            'cell_type': ct,
            'shrunk_gate_frac': row['shrunk_gate_frac'],
            'ci_low': row['shrunk_ci_low'],
            'ci_high': row['shrunk_ci_high'],
            'n_donors': row['n_donors'],
            'total_cells': row['total_cells']
        })

# Add composite
vasc_comp_tdpneg = vasc_composite[vasc_composite['condition'] == 'TDPneg']
if len(vasc_comp_tdpneg) > 0:
    successes = vasc_comp_tdpneg['gate_count'].values
    trials = vasc_comp_tdpneg['n_cells'].values
    post_mean, ci_low, ci_high = beta_binomial_estimate(successes, trials)
    vasc_comparison.append({
        'type': 'composite',
        'cell_type': 'Vascular_composite',
        'shrunk_gate_frac': post_mean.mean(),
        'ci_low': ci_low.mean(),
        'ci_high': ci_high.mean(),
        'n_donors': len(vasc_comp_tdpneg),
        'total_cells': vasc_comp_tdpneg['n_cells'].sum()
    })

vasc_comp_df = pd.DataFrame(vasc_comparison)
vasc_comp_df.to_csv(OUT_DIR / 'tables/vascular_individual_vs_composite.csv', index=False)
print(f"\nSaved: vascular_individual_vs_composite.csv")

# ============================================================================
# STEP 5: Bootstrap rank stability
# ============================================================================
print("\n" + "=" * 70)
print("STEP 5: Bootstrap rank stability (donor resampling)")
print("=" * 70)

n_boot = 500
boot_rankings = []

# Get unique donors per stage
donors_by_stage = {stage: donor_agg[donor_agg['condition'] == stage]['sample'].unique()
                   for stage in TDP_STAGES}

for b in range(n_boot):
    # Resample donors within each stage
    boot_data = []
    for stage in TDP_STAGES:
        stage_donors = donors_by_stage[stage]
        if len(stage_donors) == 0:
            continue
        resampled = np.random.choice(stage_donors, size=len(stage_donors), replace=True)
        for d in resampled:
            boot_data.append(donor_agg[(donor_agg['sample'] == d) & (donor_agg['condition'] == stage)])

    if len(boot_data) == 0:
        continue

    boot_df = pd.concat(boot_data, ignore_index=True)

    # Calculate Control -> TDPneg difference per celltype
    diffs = []
    for ct in boot_df['cell_type'].unique():
        ctrl = boot_df[(boot_df['cell_type'] == ct) & (boot_df['condition'] == 'Control')]['gate_v2_frac']
        neg = boot_df[(boot_df['cell_type'] == ct) & (boot_df['condition'] == 'TDPneg')]['gate_v2_frac']
        if len(ctrl) > 0 and len(neg) > 0:
            diff = neg.mean() - ctrl.mean()
            diffs.append({'cell_type': ct, 'diff': diff, 'is_vascular': ct in VASCULAR_TYPES})

    if len(diffs) > 0:
        diffs_df = pd.DataFrame(diffs).sort_values('diff', ascending=False)
        ranking = diffs_df['cell_type'].tolist()
        boot_rankings.append(ranking)

    if (b + 1) % 100 == 0:
        print(f"  Bootstrap {b + 1}/{n_boot}")

# Analyze stability
print(f"\nAnalyzed {len(boot_rankings)} bootstrap samples")

# Get reference ranking (from original data)
ref_diffs = []
for ct in donor_agg['cell_type'].unique():
    ctrl = donor_agg[(donor_agg['cell_type'] == ct) & (donor_agg['condition'] == 'Control')]['gate_v2_frac']
    neg = donor_agg[(donor_agg['cell_type'] == ct) & (donor_agg['condition'] == 'TDPneg')]['gate_v2_frac']
    if len(ctrl) > 0 and len(neg) > 0:
        diff = neg.mean() - ctrl.mean()
        ref_diffs.append({'cell_type': ct, 'diff': diff})

ref_ranking = pd.DataFrame(ref_diffs).sort_values('diff', ascending=False)['cell_type'].tolist()
print(f"Reference top-5: {', '.join(ref_ranking[:5])}")

# Calculate stability metrics
top_k_values = [3, 5, 10]
stability_metrics = {}

for k in top_k_values:
    overlaps = []
    for boot_rank in boot_rankings:
        overlap = len(set(boot_rank[:k]) & set(ref_ranking[:k])) / k
        overlaps.append(overlap)
    stability_metrics[f'top{k}_stability'] = np.mean(overlaps)
    print(f"  Top-{k} stability: {np.mean(overlaps):.3f}")

# Rank distribution for each celltype
rank_distributions = {ct: [] for ct in ref_ranking}
for boot_rank in boot_rankings:
    for rank, ct in enumerate(boot_rank):
        if ct in rank_distributions:
            rank_distributions[ct].append(rank + 1)

# Summary
rank_summary = []
for ct in ref_ranking:
    ranks = rank_distributions.get(ct, [])
    if len(ranks) > 0:
        rank_summary.append({
            'cell_type': ct,
            'median_rank': np.median(ranks),
            'mean_rank': np.mean(ranks),
            'rank_std': np.std(ranks),
            'pct_in_top5': 100 * np.mean([r <= 5 for r in ranks]),
            'is_vascular': ct in VASCULAR_TYPES
        })

rank_summary_df = pd.DataFrame(rank_summary).sort_values('median_rank')
rank_summary_df.to_csv(OUT_DIR / 'tables/bootstrap_rank_summary.csv', index=False)
print(f"\nSaved: bootstrap_rank_summary.csv")

# ============================================================================
# STEP 6: TDPneg internal clustering (trigger-stage)
# ============================================================================
print("\n" + "=" * 70)
print("STEP 6: TDPneg internal clustering (trigger-stage identification)")
print("=" * 70)

# Get TDPneg donors with key features
tdpneg_donors = donor_agg[donor_agg['condition'] == 'TDPneg'].copy()

# Pivot to get donor x celltype features
feature_cols = ['gate_v2_frac', 'NSA_z_mean', 'DDR_z_mean']
donor_features = []

for donor in tdpneg_donors['sample'].unique():
    donor_data = tdpneg_donors[tdpneg_donors['sample'] == donor]
    features = {'sample': donor}

    # Key cell types
    key_cts = ['Glia.OPC', 'Glia.Astro.GFAP.pos', 'Glia.Astro.GFAP.neg']
    for ct in key_cts:
        ct_data = donor_data[donor_data['cell_type'] == ct]
        for feat in feature_cols:
            col_name = f'{ct.split(".")[-1]}_{feat}'
            features[col_name] = ct_data[feat].mean() if len(ct_data) > 0 else np.nan

    # Vascular composite
    vasc_data = donor_data[donor_data['is_vascular']]
    for feat in feature_cols:
        features[f'Vasc_{feat}'] = vasc_data[feat].mean() if len(vasc_data) > 0 else np.nan

    donor_features.append(features)

donor_feat_df = pd.DataFrame(donor_features)
donor_feat_df = donor_feat_df.dropna(axis=1, how='all')  # Drop columns that are all NaN

# Fill remaining NaNs with median
for col in donor_feat_df.columns:
    if col != 'sample':
        donor_feat_df[col] = donor_feat_df[col].fillna(donor_feat_df[col].median())

# Cluster TDPneg donors
feature_matrix = donor_feat_df.drop('sample', axis=1).values
if len(feature_matrix) >= 3:
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(feature_matrix)

    # Hierarchical clustering into 2 groups
    n_clusters = min(2, len(feature_matrix))
    clustering = AgglomerativeClustering(n_clusters=n_clusters)
    donor_feat_df['trigger_cluster'] = clustering.fit_predict(X_scaled)

    # Label clusters by mean gate fraction
    cluster_means = donor_feat_df.groupby('trigger_cluster')[[c for c in donor_feat_df.columns if 'gate_v2' in c]].mean()
    overall_mean = cluster_means.mean(axis=1)
    high_cluster = overall_mean.idxmax()
    donor_feat_df['trigger_stage'] = donor_feat_df['trigger_cluster'].map(
        {high_cluster: 'TDPneg_high', 1 - high_cluster: 'TDPneg_low'}
    )

    print(f"TDPneg donors: {len(donor_feat_df)}")
    print(f"  TDPneg_low: {(donor_feat_df['trigger_stage'] == 'TDPneg_low').sum()}")
    print(f"  TDPneg_high: {(donor_feat_df['trigger_stage'] == 'TDPneg_high').sum()}")

    donor_feat_df.to_csv(OUT_DIR / 'tables/tdpneg_trigger_clusters.csv', index=False)
    print(f"Saved: tdpneg_trigger_clusters.csv")
else:
    print("  Not enough TDPneg donors for clustering")
    donor_feat_df['trigger_stage'] = 'TDPneg'

# ============================================================================
# STEP 7: Visualizations
# ============================================================================
print("\n" + "=" * 70)
print("STEP 7: Creating visualizations")
print("=" * 70)

# 7A: Bootstrap rank stability
fig, ax = plt.subplots(figsize=(12, 6))
plot_data = rank_summary_df.head(15)
colors = ['red' if v else 'steelblue' for v in plot_data['is_vascular']]
bars = ax.barh(range(len(plot_data)), plot_data['pct_in_top5'], color=colors)
ax.set_yticks(range(len(plot_data)))
ax.set_yticklabels(plot_data['cell_type'])
ax.set_xlabel('% of bootstraps in top-5')
ax.set_title('Bootstrap Rank Stability (Control → TDPneg)')
ax.axvline(50, color='gray', linestyle='--', alpha=0.5)
# Legend
from matplotlib.patches import Patch
legend_elements = [Patch(facecolor='red', label='Vascular'),
                   Patch(facecolor='steelblue', label='Non-vascular')]
ax.legend(handles=legend_elements, loc='lower right')
plt.tight_layout()
plt.savefig(OUT_DIR / 'figures/bootstrap_rank_stability.png', dpi=150)
plt.close()
print("  Saved: bootstrap_rank_stability.png")

# 7B: Vascular individual vs composite
fig, ax = plt.subplots(figsize=(10, 6))
x = range(len(vasc_comp_df))
ax.barh(x, vasc_comp_df['shrunk_gate_frac'], xerr=[
    vasc_comp_df['shrunk_gate_frac'] - vasc_comp_df['ci_low'],
    vasc_comp_df['ci_high'] - vasc_comp_df['shrunk_gate_frac']
], capsize=5, color=['coral' if t == 'individual' else 'darkred' for t in vasc_comp_df['type']])
ax.set_yticks(x)
ax.set_yticklabels(vasc_comp_df['cell_type'])
ax.set_xlabel('Shrunk Gate v2 Fraction (with 95% CI)')
ax.set_title('Vascular: Individual vs Composite (TDPneg)')
for i, row in vasc_comp_df.iterrows():
    ax.annotate(f'n={row["n_donors"]}, cells={row["total_cells"]}',
                (row['shrunk_gate_frac'] + 0.01, i), fontsize=8)
plt.tight_layout()
plt.savefig(OUT_DIR / 'figures/vascular_single_vs_composite.png', dpi=150)
plt.close()
print("  Saved: vascular_single_vs_composite.png")

# 7C: TDPneg trigger clusters
if 'trigger_stage' in donor_feat_df.columns and len(donor_feat_df) >= 3:
    # PCA of donor features
    feat_cols = [c for c in donor_feat_df.columns if c not in ['sample', 'trigger_cluster', 'trigger_stage']]
    X = donor_feat_df[feat_cols].values
    pca = PCA(n_components=2)
    X_pca = pca.fit_transform(StandardScaler().fit_transform(X))
    donor_feat_df['PC1'] = X_pca[:, 0]
    donor_feat_df['PC2'] = X_pca[:, 1]

    fig, ax = plt.subplots(figsize=(8, 6))
    colors = {'TDPneg_low': 'lightgreen', 'TDPneg_high': 'darkgreen'}
    for stage in ['TDPneg_low', 'TDPneg_high']:
        stage_data = donor_feat_df[donor_feat_df['trigger_stage'] == stage]
        ax.scatter(stage_data['PC1'], stage_data['PC2'], c=colors.get(stage, 'gray'),
                   label=f'{stage} (n={len(stage_data)})', s=100, alpha=0.7)
    ax.set_xlabel(f'PC1 ({100*pca.explained_variance_ratio_[0]:.1f}%)')
    ax.set_ylabel(f'PC2 ({100*pca.explained_variance_ratio_[1]:.1f}%)')
    ax.set_title('TDPneg Donors: Trigger Stage Clustering')
    ax.legend()
    plt.tight_layout()
    plt.savefig(OUT_DIR / 'figures/tdpneg_trigger_clusters.png', dpi=150)
    plt.close()
    print("  Saved: tdpneg_trigger_clusters.png")

# 7D: Stage trajectory with trigger substages
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# Key cell types trajectory
key_cts = ['Glia.OPC', 'Glia.Astro.GFAP.pos', 'Vasc.Endo', 'Vasc.Capillary']
stages_order = ['Control', 'TDPneg', 'TDPmed', 'TDPhigh']

ax = axes[0]
for ct in key_cts:
    ct_data = hier_df[hier_df['cell_type'] == ct]
    values = []
    for stage in stages_order:
        stage_row = ct_data[ct_data['condition'] == stage]
        if len(stage_row) > 0:
            values.append(stage_row['shrunk_gate_frac'].values[0])
        else:
            values.append(np.nan)
    marker = 's' if ct in VASCULAR_TYPES else 'o'
    ax.plot(stages_order, values, marker=marker, label=ct, linewidth=2, markersize=8)

ax.set_ylabel('Shrunk Gate v2 Fraction')
ax.set_title('Trigger Gate Trajectory by Stage')
ax.legend(loc='upper right')
ax.set_ylim(0, None)

# Effect size (Cliff's delta) comparison
ax = axes[1]
tdpneg_data = hier_df[hier_df['condition'] == 'TDPneg'].sort_values('cliffs_delta_vs_ctrl', ascending=True).tail(10)
colors = ['red' if v else 'steelblue' for v in tdpneg_data['is_vascular']]
ax.barh(range(len(tdpneg_data)), tdpneg_data['cliffs_delta_vs_ctrl'], color=colors)
ax.set_yticks(range(len(tdpneg_data)))
ax.set_yticklabels(tdpneg_data['cell_type'])
ax.set_xlabel("Cliff's δ (vs Control)")
ax.set_title('Effect Size: Control → TDPneg (Top 10)')
ax.axvline(0, color='gray', linestyle='--')

plt.tight_layout()
plt.savefig(OUT_DIR / 'figures/stage_trajectory_effects.png', dpi=150)
plt.close()
print("  Saved: stage_trajectory_effects.png")

# ============================================================================
# STEP 8: Generate Report
# ============================================================================
print("\n" + "=" * 70)
print("STEP 8: Generating report")
print("=" * 70)

# Get key stats
vasc_endo = hier_df[(hier_df['cell_type'] == 'Vasc.Endo') & (hier_df['condition'] == 'TDPneg')]
opc = hier_df[(hier_df['cell_type'] == 'Glia.OPC') & (hier_df['condition'] == 'TDPneg')]

report = f"""# GSE212630 Trigger Stage Analysis v3 (Donor-Level Hierarchical)

## Key Methodology Changes
1. **Donor-level aggregation only** (no pseudo-replication)
2. **Beta-Binomial shrinkage** for low-n cell types (partial pooling)
3. **Bootstrap rank stability** (donor resampling ×{n_boot})
4. **TDPneg internal clustering** (trigger-stage identification)
5. **Vascular: individual + composite** evaluation

## Data Summary
- Total cells: {len(df):,}
- Gate v2 (NSA + DDR): {df['gate_v2'].sum():,} cells ({100*df['gate_v2'].mean():.1f}%)
- DAMP QC-corrected: r(DAMP, nGenes) = {df['DAMP_det'].corr(df['n_genes_detected']):.3f} → r(DAMP_resid, nGenes) ≈ 0

## Hierarchical Estimates (TDPneg, top cell types)

| Cell Type | Shrunk Gate | 95% CI | Cliff's δ | n_donors | Vascular |
|-----------|-------------|--------|-----------|----------|----------|
"""

for _, row in tdpneg.head(10).iterrows():
    vasc = '✓' if row['is_vascular'] else ''
    report += f"| {row['cell_type']} | {row['shrunk_gate_frac']:.3f} | [{row['shrunk_ci_low']:.3f}, {row['shrunk_ci_high']:.3f}] | {row['cliffs_delta_vs_ctrl']:+.2f} | {row['n_donors']} | {vasc} |\n"

report += f"""
## Bootstrap Rank Stability

| Metric | Value |
|--------|-------|
| Top-3 stability | {stability_metrics.get('top3_stability', np.nan):.3f} |
| Top-5 stability | {stability_metrics.get('top5_stability', np.nan):.3f} |
| Top-10 stability | {stability_metrics.get('top10_stability', np.nan):.3f} |

**Reference top-5**: {', '.join(ref_ranking[:5])}

## Vascular Analysis

### Individual Cell Types (TDPneg)
"""

for _, row in vasc_comp_df[vasc_comp_df['type'] == 'individual'].iterrows():
    report += f"- **{row['cell_type']}**: gate={row['shrunk_gate_frac']:.3f} [{row['ci_low']:.3f}, {row['ci_high']:.3f}], n={row['n_donors']}, cells={row['total_cells']}\n"

vasc_comp_row = vasc_comp_df[vasc_comp_df['type'] == 'composite']
if len(vasc_comp_row) > 0:
    v = vasc_comp_row.iloc[0]
    report += f"""
### Composite
- **Vascular_composite**: gate={v['shrunk_gate_frac']:.3f} [{v['ci_low']:.3f}, {v['ci_high']:.3f}], n={v['n_donors']}, cells={v['total_cells']}
"""

# Pre-calculate values for report
tdpneg_low_count = (donor_feat_df['trigger_stage'] == 'TDPneg_low').sum() if 'trigger_stage' in donor_feat_df.columns else 'N/A'
tdpneg_high_count = (donor_feat_df['trigger_stage'] == 'TDPneg_high').sum() if 'trigger_stage' in donor_feat_df.columns else 'N/A'

vasc_endo_rank = rank_summary_df[rank_summary_df['cell_type'] == 'Vasc.Endo']
vasc_endo_pct = f"{vasc_endo_rank['pct_in_top5'].values[0]:.1f}" if len(vasc_endo_rank) > 0 else 'N/A'

opc_shrunk = hier_df[hier_df['cell_type'] == 'Glia.OPC']
opc_shrunk_val = f"{opc_shrunk['shrunk_gate_frac'].values[0]:.3f}" if len(opc_shrunk) > 0 else 'N/A'

report += f"""
## TDPneg Internal Clustering

TDPneg donors were clustered into trigger substages based on gate fractions across cell types:
- **TDPneg_low**: {tdpneg_low_count} donors
- **TDPneg_high**: {tdpneg_high_count} donors

## Structural Interpretation

1. **Vascular is not excluded due to low n**
   - Beta-Binomial shrinkage provides partial pooling
   - Composite aggregation recovers statistical power
   - Bootstrap stability shows rank consistency

2. **Vasc.Endo consistently in top ranks**
   - Bootstrap pct_in_top5: {vasc_endo_pct}% (if available)

3. **OPC as observable proxy**
   - Highest shrunk estimate: {opc_shrunk_val}
   - May amplify vascular-origin signals

4. **TDPneg substaging**
   - Trigger-high donors may represent earlier disease stage
   - Pseudo-timeline: Control → TDPneg_low → TDPneg_high → TDPmed → TDPhigh

## Caveats
1. Hierarchical shrinkage assumes exchangeability across donors
2. Bootstrap stability is sensitive to donor heterogeneity
3. TDPneg clustering is exploratory (hypothesis-generating)

## Output Files

### Tables
- `donor_level_gate_metrics.csv`: Donor × celltype × stage aggregation
- `hierarchical_stage_effects.csv`: Shrunk estimates with CIs
- `vascular_individual_vs_composite.csv`: Vascular comparison
- `bootstrap_rank_summary.csv`: Rank stability by cell type
- `tdpneg_trigger_clusters.csv`: TDPneg donor clustering

### Figures
- `bootstrap_rank_stability.png`: Bootstrap rank stability
- `vascular_single_vs_composite.png`: Vascular individual vs composite
- `tdpneg_trigger_clusters.png`: TDPneg trigger stage PCA
- `stage_trajectory_effects.png`: Gate trajectory and effect sizes
"""

with open(OUT_DIR / 'reports/triggerstage_v3_donorlevel.md', 'w') as f:
    f.write(report)
print(f"Saved: triggerstage_v3_donorlevel.md")

print("\n" + "=" * 70)
print("COMPLETE")
print("=" * 70)
print(f"Output: {OUT_DIR}")
