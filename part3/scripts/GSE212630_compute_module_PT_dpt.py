#!/usr/bin/env python3
"""
GSE_212630: Compute Module-based PT_dpt (Stress-Independent)
=============================================================

This script computes:
1. Module scores from expression data (23 functional modules)
2. PT_dpt_module using Diffusion Pseudotime on module features
3. Validation against TDP-43 stages and comparison with PT_pca

The key difference from the previous PT_dpt (stress-based):
- Uses module scores as features (like the existing project)
- Should be more independent of stress scores
- Provides better temporal ordering

Author: Claude Code + User
Date: 2025-12-07
"""

import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.neighbors import NearestNeighbors
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import eigsh
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import gc
import warnings
warnings.filterwarnings('ignore')

print("=" * 80)
print("GSE_212630: Compute Module-based PT_dpt")
print("=" * 80)

# ============================================================================
# CONFIGURATION
# ============================================================================

DATA_DIR = Path(__file__).parent.parent / 'GSE_212630_raw_expression_transposed'
OUTPUT_DIR = Path(__file__).parent.parent / 'results' / 'GSE212630_ids_analysis'

# Module gene sets (from existing project + extensions)
MODULE_GENE_SETS = {
    # Metabolic Axis
    'Mitochondria': [
        'MT-ND1', 'MT-ND2', 'MT-ND3', 'MT-ND4', 'MT-ND5', 'MT-ND6',
        'MT-CO1', 'MT-CO2', 'MT-CO3', 'MT-ATP6', 'MT-ATP8',
        'NDUFA1', 'NDUFA2', 'NDUFS1', 'NDUFS2', 'SDHA', 'SDHB',
        'COX4I1', 'COX5A', 'COX6A1', 'ATP5F1A', 'ATP5F1B'
    ],
    'ER_Stress': [
        'HSPA5', 'DDIT3', 'ATF4', 'ATF6', 'XBP1', 'ERN1', 'EIF2AK3',
        'DNAJB9', 'DNAJC3', 'CALR', 'CANX', 'PDIA4', 'PDIA6'
    ],
    'Protein_Homeostasis': [
        'HSPA1A', 'HSPA1B', 'HSP90AA1', 'HSP90AB1', 'DNAJA1', 'DNAJB1',
        'UBB', 'UBC', 'PSMA1', 'PSMB1', 'SQSTM1', 'LAMP1', 'LAMP2'
    ],
    'Metabolism': [
        'HK1', 'HK2', 'PFKM', 'PFKP', 'PKM', 'LDHA', 'LDHB',
        'PDK1', 'PDK2', 'CS', 'ACO2', 'IDH1', 'IDH2', 'OGDH'
    ],

    # Hyperexcitability Axis
    'Calcium_Signaling': [
        'CACNA1A', 'CACNA1B', 'CACNA1C', 'CACNA1D', 'CACNA1E',
        'CALM1', 'CALM2', 'CALM3', 'CAMK2A', 'CAMK2B', 'CAMK2G',
        'ATP2A1', 'ATP2A2', 'ATP2B1', 'ATP2B2', 'RYR1', 'RYR2'
    ],
    'Ion_Transport': [
        'KCNA1', 'KCNA2', 'KCNB1', 'KCNQ2', 'KCNQ3',
        'SCN1A', 'SCN2A', 'SCN8A', 'CLCN1', 'CLCN2',
        'SLC12A2', 'SLC12A5', 'ATP1A1', 'ATP1A2', 'ATP1A3'
    ],

    # Inflammatory Axis
    'Oxidative_Stress': [
        'SOD1', 'SOD2', 'CAT', 'GPX1', 'GPX4', 'PRDX1', 'PRDX2',
        'TXN', 'TXNRD1', 'NQO1', 'HMOX1', 'GCLC', 'GCLM', 'GSR'
    ],
    'Inflammation': [
        'IL1B', 'IL6', 'IL18', 'TNF', 'TNFRSF1A', 'TNFRSF1B',
        'NFKB1', 'NFKB2', 'RELA', 'CXCL8', 'CCL2', 'CCL5',
        'TGFB1', 'TGFB2', 'TGFB3'
    ],
    'Complement': [
        'C1QA', 'C1QB', 'C1QC', 'C3', 'C4A', 'C4B',
        'CFB', 'CFH', 'CFI', 'SERPING1', 'C5', 'C6', 'C7'
    ],

    # Structural/ECM
    'Cytoskeleton': [
        'NEFL', 'NEFM', 'NEFH', 'TUBB3', 'TUBB2A', 'TUBA1A',
        'MAP2', 'MAPT', 'ACTB', 'ACTG1', 'VIM', 'DES', 'GFAP'
    ],
    'ECM': [
        'COL1A1', 'COL1A2', 'COL3A1', 'COL4A1', 'FN1', 'LAMA1',
        'LAMB1', 'LAMC1', 'TNC', 'VCAN', 'ACAN'
    ],
    'Myelination': [
        'MBP', 'PLP1', 'MOG', 'MAG', 'MOBP', 'CNP',
        'OLIG1', 'OLIG2', 'SOX10', 'NKX2-2', 'MYRF'
    ],

    # Vascular
    'Angiogenesis': [
        'VEGFA', 'VEGFB', 'VEGFC', 'FLT1', 'KDR', 'FLT4',
        'ANGPT1', 'ANGPT2', 'TEK', 'PECAM1', 'CDH5', 'VWF'
    ],

    # Cellular Processes
    'Apoptosis': [
        'BCL2', 'BCL2L1', 'BAX', 'BAK1', 'BID', 'CASP3',
        'CASP8', 'CASP9', 'CYCS', 'APAF1', 'XIAP', 'BIRC5'
    ],
    'Autophagy': [
        'BECN1', 'ATG5', 'ATG7', 'ATG12', 'MAP1LC3A', 'MAP1LC3B',
        'SQSTM1', 'NBR1', 'OPTN', 'BNIP3', 'BNIP3L'
    ],
    'Cell_Cycle': [
        'CCND1', 'CCND2', 'CCNE1', 'CDK2', 'CDK4', 'CDK6',
        'RB1', 'E2F1', 'TP53', 'CDKN1A', 'CDKN2A'
    ],
    'DNA_Repair': [
        'TP53', 'ATM', 'ATR', 'BRCA1', 'BRCA2', 'RAD51',
        'XRCC1', 'PARP1', 'MSH2', 'MSH6', 'MLH1'
    ],

    # Signaling
    'Growth_Factors': [
        'BDNF', 'NGF', 'GDNF', 'CNTF', 'IGF1', 'IGF2',
        'FGF1', 'FGF2', 'EGF', 'PDGFA', 'PDGFB'
    ],
    'Synaptic': [
        'SYP', 'SYN1', 'SYT1', 'SNAP25', 'STX1A', 'VAMP2',
        'DLG4', 'GRIN1', 'GRIN2A', 'GRIN2B', 'GRIA1', 'GRIA2',
        'GABRA1', 'GABRB2', 'SLC17A7', 'SLC32A1'
    ],

    # Transcriptional
    'Transcription': [
        'FOS', 'JUN', 'MYC', 'EGR1', 'CREB1', 'ATF3',
        'NR4A1', 'NR4A2', 'NPAS4', 'ARC', 'BDNF'
    ],
    'RNA_Processing': [
        'TARDBP', 'FUS', 'HNRNPA1', 'HNRNPA2B1', 'SRSF1', 'SRSF2',
        'ELAVL1', 'ELAVL2', 'ELAVL3', 'ELAVL4', 'RBFOX1', 'RBFOX2'
    ],

    # TDP-43 specific
    'TDP43_targets': [
        'STMN2', 'UNC13A', 'ELAVL4', 'KCNQ2', 'PFKP',
        'TARDBP', 'SORT1', 'KCNQ3', 'ATXN1', 'KCNIP3'
    ],
    'STMN2_pathway': [
        'STMN2', 'STMN1', 'STMN3', 'STMN4',
        'MAPT', 'MAP2', 'MAP1B', 'DPYSL2', 'DPYSL3'
    ]
}

TDP_STAGE_ORDER = ['Control', 'TDPneg', 'TDPmed', 'TDPhigh']
TDP_STAGE_TO_NUMERIC = {'Control': 0, 'TDPneg': 1, 'TDPmed': 2, 'TDPhigh': 3}

MAX_CELLS_PER_TYPE = 500
K_NEIGHBORS = 30
N_COMPONENTS = 10

np.random.seed(42)

# ============================================================================
# STEP 1: LOAD EXPRESSION DATA AND COMPUTE MODULE SCORES
# ============================================================================

print("\n[Step 1] Loading expression data and computing module scores...")

summary_df = pd.read_csv(DATA_DIR / 'extraction_summary.csv')
valid_df = summary_df[(summary_df['n_cells'] >= 50) & (summary_df['n_cells'] <= 10000)].sort_values('n_cells')
valid_cell_types = valid_df['cell_type'].tolist()

print(f"  Cell types to process: {len(valid_cell_types)}")

all_metadata = []

for ct in valid_cell_types:
    ct_safe = ct.replace('/', '_').replace('.', '_')
    expr_file = None
    meta_file = None

    for f in DATA_DIR.glob('*_expression_transposed.csv.gz'):
        fname = f.stem.replace('_expression_transposed.csv', '')
        if ct_safe in fname or fname.replace('_', '.') == ct.replace('/', '_'):
            expr_file = f
            meta_file = DATA_DIR / f.name.replace('_expression_transposed.csv.gz', '_metadata.csv')
            break

    if expr_file is None or not meta_file.exists():
        continue

    print(f"\n  Processing {ct}...")

    try:
        meta_df = pd.read_csv(meta_file)
        meta_df['cell_type'] = ct

        # Sample if too many
        if len(meta_df) > MAX_CELLS_PER_TYPE:
            sampled_idx = []
            for cond in meta_df['condition'].unique():
                cond_idx = meta_df[meta_df['condition'] == cond].index.tolist()
                n_sample = min(len(cond_idx), MAX_CELLS_PER_TYPE // meta_df['condition'].nunique())
                sampled_idx.extend(np.random.choice(cond_idx, size=n_sample, replace=False).tolist())
            meta_df = meta_df.loc[sampled_idx].reset_index(drop=True)

        # Build column names
        meta_df['expr_col'] = meta_df['cell_id'] + '_' + meta_df['sample'].astype(str) + '_' + meta_df['condition']

        # Load expression
        print(f"    Loading expression...")
        expr_df = pd.read_csv(expr_file, index_col=0)

        # Filter columns
        cells_to_load = meta_df['expr_col'].tolist()
        available_cols = [c for c in cells_to_load if c in expr_df.columns]
        if len(available_cols) == 0:
            continue

        expr_df = expr_df[available_cols]
        meta_df = meta_df[meta_df['expr_col'].isin(available_cols)].reset_index(drop=True)

        # Transpose
        expr_T = expr_df.T
        expr_T = expr_T.loc[meta_df['expr_col'].values]

        # Gene names (uppercase for matching)
        gene_names = [g.upper() for g in expr_T.columns.tolist()]
        gene_name_to_idx = {g: i for i, g in enumerate(gene_names)}

        # Compute module scores
        print(f"    Computing module scores...")
        X = expr_T.values

        for module_name, gene_list in MODULE_GENE_SETS.items():
            gene_idx = []
            for g in gene_list:
                g_upper = g.upper()
                if g_upper in gene_name_to_idx:
                    gene_idx.append(gene_name_to_idx[g_upper])

            if len(gene_idx) == 0:
                meta_df[f'module_{module_name}'] = 0.0
            else:
                module_expr = X[:, gene_idx].mean(axis=1)
                # Z-score normalize
                module_mean = module_expr.mean()
                module_std = module_expr.std() + 1e-8
                meta_df[f'module_{module_name}'] = (module_expr - module_mean) / module_std

        all_metadata.append(meta_df)
        print(f"    Done: {len(meta_df)} cells, {len([c for c in meta_df.columns if c.startswith('module_')])} modules")

        del expr_df, expr_T, X
        gc.collect()

    except Exception as e:
        print(f"    [ERROR] {e}")
        gc.collect()
        continue

# Combine
print("\n  Combining all data...")
combined_df = pd.concat(all_metadata, ignore_index=True)
print(f"  Total cells: {len(combined_df):,}")

# ============================================================================
# STEP 2: COMPUTE PT_dpt_module USING MODULE FEATURES
# ============================================================================

print("\n[Step 2] Computing PT_dpt_module using module features...")

# Module columns
module_cols = [col for col in combined_df.columns if col.startswith('module_')]
print(f"  Module columns: {len(module_cols)}")

# Fill NaN
for col in module_cols:
    combined_df[col] = combined_df[col].fillna(0)

# Feature matrix
X = combined_df[module_cols].values
print(f"  Feature matrix: {X.shape}")

# Standardize
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

# Build kNN graph
print(f"\n  Building kNN graph (K={K_NEIGHBORS})...")
nbrs = NearestNeighbors(n_neighbors=K_NEIGHBORS + 1, metric='euclidean', n_jobs=-1)
nbrs.fit(X_scaled)
distances, indices = nbrs.kneighbors(X_scaled)

# Affinity matrix
print("  Building affinity matrix...")
sigma = np.median(distances[:, 1:])
n_cells = len(combined_df)

row_idx = []
col_idx = []
affinity_values = []

for i in range(n_cells):
    for j_idx, j in enumerate(indices[i, 1:]):
        dist = distances[i, j_idx + 1]
        affinity = np.exp(-(dist ** 2) / (2 * sigma ** 2))
        row_idx.append(i)
        col_idx.append(j)
        affinity_values.append(affinity)

W = csr_matrix((affinity_values, (row_idx, col_idx)), shape=(n_cells, n_cells))
W = (W + W.T) / 2

# Diffusion operator
print("  Computing diffusion operator...")
D = np.array(W.sum(axis=1)).flatten()
D_inv_sqrt = 1.0 / np.sqrt(D + 1e-10)
D_inv_sqrt_sparse = csr_matrix((D_inv_sqrt, (range(n_cells), range(n_cells))))
L = D_inv_sqrt_sparse @ W @ D_inv_sqrt_sparse

# Eigendecomposition
print(f"  Computing top {N_COMPONENTS} eigenvectors...")
eigenvalues, eigenvectors = eigsh(L, k=N_COMPONENTS, which='LM')
idx = eigenvalues.argsort()[::-1]
eigenvalues = eigenvalues[idx]
eigenvectors = eigenvectors[:, idx]

# Select root cell: Control cell with lowest stress
print("\n  Selecting root cell...")
control_mask = combined_df['condition'] == 'Control'

# Compute simple stress for root selection
stress_modules = ['module_ER_Stress', 'module_Oxidative_Stress', 'module_Inflammation']
stress_cols_available = [c for c in stress_modules if c in combined_df.columns]
if len(stress_cols_available) > 0:
    combined_df['_temp_stress'] = combined_df[stress_cols_available].mean(axis=1)
else:
    combined_df['_temp_stress'] = 0

control_indices = combined_df[control_mask].index.tolist()
control_stress = combined_df.loc[control_indices, '_temp_stress']
root_idx = control_stress.idxmin()
root_array_idx = combined_df.index.get_loc(root_idx)

root_cell = combined_df.loc[root_idx]
print(f"  Root cell: {root_cell['cell_type']}, condition={root_cell['condition']}")

# DPT distance
dpt_components = eigenvectors[:, 1:6]
root_vector = dpt_components[root_array_idx, :]
dpt_distances = np.linalg.norm(dpt_components - root_vector, axis=1)

# Normalize
PT_dpt_module = (dpt_distances - dpt_distances.min()) / (dpt_distances.max() - dpt_distances.min())
combined_df['PT_dpt_module'] = PT_dpt_module

print(f"  PT_dpt_module range: [{PT_dpt_module.min():.4f}, {PT_dpt_module.max():.4f}]")
print(f"  PT_dpt_module mean: {PT_dpt_module.mean():.4f}")

# Cleanup
combined_df = combined_df.drop(columns=['_temp_stress'])

# ============================================================================
# STEP 3: VALIDATE AGAINST TDP-43 STAGES
# ============================================================================

print("\n[Step 3] Validating against TDP-43 stages...")

combined_df['TDP_numeric'] = combined_df['condition'].map(TDP_STAGE_TO_NUMERIC)

# Spearman correlation
corr, p = stats.spearmanr(combined_df['TDP_numeric'], combined_df['PT_dpt_module'])
print(f"\n  Spearman correlation (TDP stage vs PT_dpt_module):")
print(f"    rho = {corr:.4f}")
print(f"    p-value = {p:.2e}")

# Mean by stage
print(f"\n  Mean PT_dpt_module by TDP-43 stage:")
for stage in TDP_STAGE_ORDER:
    stage_data = combined_df[combined_df['condition'] == stage]['PT_dpt_module']
    print(f"    {stage:10s}: mean={stage_data.mean():.4f} ± {stage_data.std():.4f} (n={len(stage_data)})")

# ============================================================================
# STEP 4: COMPARE WITH OTHER PT METHODS
# ============================================================================

print("\n[Step 4] Comparing with other PT methods...")

# Load previous PT results
prev_file = OUTPUT_DIR / 'combined_metadata_with_all_PTs.csv'
if prev_file.exists():
    prev_df = pd.read_csv(prev_file)

    # Merge PT_pca and PT_dpt
    combined_df = combined_df.merge(
        prev_df[['expr_col', 'PT_pca', 'PT_dpt']].drop_duplicates(),
        on='expr_col',
        how='left'
    )

    # Compare
    pt_methods = ['PT_pca', 'PT_dpt', 'PT_dpt_module']
    results = []

    for pt_col in pt_methods:
        if pt_col not in combined_df.columns or combined_df[pt_col].isna().all():
            continue

        valid_mask = combined_df[pt_col].notna()
        corr, p = stats.spearmanr(
            combined_df.loc[valid_mask, 'TDP_numeric'],
            combined_df.loc[valid_mask, pt_col]
        )

        results.append({
            'method': pt_col,
            'spearman_rho': corr,
            'p_value': p,
            'n_cells': valid_mask.sum()
        })

        print(f"\n  {pt_col}:")
        print(f"    Spearman ρ = {corr:.4f}, p = {p:.2e}")

    results_df = pd.DataFrame(results)
    results_df.to_csv(OUTPUT_DIR / 'PT_comparison_with_module.csv', index=False)
    print(f"\n  Saved: PT_comparison_with_module.csv")

# ============================================================================
# STEP 5: SAVE RESULTS
# ============================================================================

print("\n[Step 5] Saving results...")

output_file = OUTPUT_DIR / 'combined_metadata_with_modules_and_PT.csv'
combined_df.to_csv(output_file, index=False)
print(f"  Saved: {output_file}")

# ============================================================================
# STEP 6: GENERATE FIGURES
# ============================================================================

print("\n[Step 6] Generating figures...")

fig, axes = plt.subplots(2, 3, figsize=(18, 12))

colors = ['#2ecc71', '#f1c40f', '#e67e22', '#e74c3c']

# Panel 1: PT_dpt_module by TDP stage
ax = axes[0, 0]
bp_data = [combined_df[combined_df['condition'] == stage]['PT_dpt_module'].values for stage in TDP_STAGE_ORDER]
bp = ax.boxplot(bp_data, labels=TDP_STAGE_ORDER, patch_artist=True)
for patch, color in zip(bp['boxes'], colors):
    patch.set_facecolor(color)
    patch.set_alpha(0.7)
ax.set_xlabel('TDP-43 Stage', fontsize=12)
ax.set_ylabel('PT_dpt_module', fontsize=12)
ax.set_title(f'PT_dpt_module (Module-based DPT)\nSpearman ρ = {corr:.3f}, p = {p:.2e}', fontsize=12)
ax.grid(True, alpha=0.3)

# Panel 2: Distribution by condition
ax = axes[0, 1]
for stage, color in zip(TDP_STAGE_ORDER, colors):
    data = combined_df[combined_df['condition'] == stage]['PT_dpt_module']
    ax.hist(data, bins=30, alpha=0.5, label=stage, color=color, density=True)
ax.set_xlabel('PT_dpt_module', fontsize=12)
ax.set_ylabel('Density', fontsize=12)
ax.set_title('PT_dpt_module Distribution by Condition', fontsize=12)
ax.legend()

# Panel 3: Mean trend
ax = axes[0, 2]
means = [combined_df[combined_df['condition'] == stage]['PT_dpt_module'].mean() for stage in TDP_STAGE_ORDER]
stds = [combined_df[combined_df['condition'] == stage]['PT_dpt_module'].std() for stage in TDP_STAGE_ORDER]
ax.errorbar(range(4), means, yerr=stds, fmt='o-', capsize=5, capthick=2, linewidth=2, markersize=10, color='purple')
ax.set_xticks(range(4))
ax.set_xticklabels(TDP_STAGE_ORDER)
ax.set_xlabel('TDP-43 Stage', fontsize=12)
ax.set_ylabel('Mean PT_dpt_module ± SD', fontsize=12)
ax.set_title('PT_dpt_module Trend', fontsize=12)
ax.grid(True, alpha=0.3)

# Panel 4: Comparison of all PT methods
ax = axes[1, 0]
if 'PT_pca' in combined_df.columns and 'PT_dpt' in combined_df.columns:
    for pt_col, marker, label, color in [
        ('PT_pca', 'o-', 'PT_pca (key genes)', '#3498db'),
        ('PT_dpt', 's-', 'PT_dpt (stress)', '#e74c3c'),
        ('PT_dpt_module', '^-', 'PT_dpt_module (modules)', '#9b59b6')
    ]:
        if pt_col in combined_df.columns:
            valid = combined_df[pt_col].notna()
            means = [combined_df[valid & (combined_df['condition'] == stage)][pt_col].mean()
                     for stage in TDP_STAGE_ORDER]
            ax.plot(range(4), means, marker, label=label, linewidth=2, markersize=8, color=color)
    ax.set_xticks(range(4))
    ax.set_xticklabels(TDP_STAGE_ORDER)
    ax.set_xlabel('TDP-43 Stage', fontsize=12)
    ax.set_ylabel('Mean PT', fontsize=12)
    ax.set_title('PT Method Comparison', fontsize=12)
    ax.legend()
    ax.grid(True, alpha=0.3)

# Panel 5: Eigenvalue spectrum
ax = axes[1, 1]
ax.plot(range(1, len(eigenvalues) + 1), eigenvalues, 'o-', color='purple', markersize=8)
ax.set_xlabel('Component', fontsize=12)
ax.set_ylabel('Eigenvalue', fontsize=12)
ax.set_title('Diffusion Eigenvalue Spectrum (Module-based)', fontsize=12)
ax.grid(True, alpha=0.3)

# Panel 6: Module heatmap by TDP stage
ax = axes[1, 2]
module_means = combined_df.groupby('condition')[module_cols[:10]].mean()  # Top 10 modules
module_means = module_means.loc[TDP_STAGE_ORDER]
module_labels = [c.replace('module_', '') for c in module_cols[:10]]
sns.heatmap(module_means.T, cmap='RdBu_r', center=0, ax=ax,
            xticklabels=TDP_STAGE_ORDER, yticklabels=module_labels)
ax.set_title('Module Scores by TDP Stage (Top 10)', fontsize=12)

plt.tight_layout()
fig.savefig(OUTPUT_DIR / 'Fig_PT_dpt_module_validation.png', dpi=150, bbox_inches='tight')
plt.close()
print(f"  Saved: Fig_PT_dpt_module_validation.png")

# ============================================================================
# SUMMARY
# ============================================================================

print("\n" + "=" * 80)
print("PT_dpt_module COMPUTATION COMPLETE")
print("=" * 80)

print(f"""
SUMMARY:

1. Module-based PT_dpt computed using {len(module_cols)} functional modules
2. Diffusion Pseudotime with K={K_NEIGHBORS} neighbors
3. Root cell: Control with lowest stress

VALIDATION RESULTS:
  - Spearman correlation (TDP stage vs PT_dpt_module): ρ = {corr:.4f}, p = {p:.2e}

COMPARISON:
""")

if 'results_df' in dir():
    for _, row in results_df.iterrows():
        print(f"  - {row['method']}: ρ = {row['spearman_rho']:.4f}, p = {row['p_value']:.2e}")

print(f"""
INTERPRETATION:
  PT_dpt_module uses functional module scores (not stress scores) as features.
  This provides a more biologically meaningful pseudotime that captures
  transcriptional state changes across multiple cellular processes.

Output: {OUTPUT_DIR}
""")
print("=" * 80)
