#!/usr/bin/env python3
"""
Phase TDP: TDP-43 Pathology Pseudotime and Upstream/Downstream Candidates

Goal:
  1) Define TDP-43 pathology pseudotime PT_TDP using STMN2/TDP43_targets modules
  2) Map cells in 2D space (PT_dpt vs PT_TDP)
  3) Compare L2/3 vs VAT1L trajectories on this phase space
  4) Examine module φ as function of PT_TDP to identify upstream/downstream candidates

Date: 2025-11-25
"""

import os
import sys
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from scipy import stats

# Add config path
sys.path.insert(0, str(Path(__file__).parent.parent))
from config.als_stress_config import ALS_STRESS_COMPONENTS

# =====================================================
# Path Settings
# =====================================================

BASE_DIR = Path("/home/akaco/als/motor_cortex_analysis")
IDS_DIR = BASE_DIR / "ids_causal_analysis"
RESULTS_DIR = IDS_DIR / "results" / "phaseTDP"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

# =====================================================
# Cell Type Definitions
# =====================================================

UPPER_L23 = ['Ex.L2.L3.CUX2.RASGRF2']
VAT1L = ['Ex.L5.VAT1L.EYA4', 'Ex.L5.VAT1L.THSD4']

VASCULAR_TYPES = [
    'Vasc.Endo.Arterial',
    'Vasc.Endo.Capillary',
    'Vasc.Endo.Venous',
    'Vasc.Mural.SMC',
    'Vasc.Fibro.CLMP.PDGFRA'
]

GLIA_TYPES = [
    'Astro',
    'Oligo',
    'OPC',
    'Micro'
]

# =====================================================
# Step 0: Load Cell-Level Features with PT_dpt
# =====================================================

def load_cell_features():
    """Load cell-level features with existing modules and PT_dpt"""
    print("\n" + "="*60)
    print("Step 0: Loading Cell-Level Features")
    print("="*60)

    csv_path = IDS_DIR / "results" / "phase9_vascular" / "cell_level_features_ALL_with_PTdpt.csv"

    print(f"Loading: {csv_path.name}")
    df = pd.read_csv(csv_path)

    print(f"  Loaded {len(df):,} cells")
    print(f"  Columns: {df.shape[1]}")
    print(f"  Conditions: {df['condition'].value_counts().to_dict()}")
    print(f"  Cell types: {df['cell_type'].nunique()}")

    # Standardize column name
    df = df.rename(columns={'condition': 'group'})

    # Check for PT_dpt
    if 'PT_dpt' not in df.columns:
        raise ValueError("PT_dpt column not found!")

    print(f"  PT_dpt range: [{df['PT_dpt'].min():.3f}, {df['PT_dpt'].max():.3f}]")

    return df


# =====================================================
# Step 1: Compute TDP-43 Module Scores
# =====================================================

def load_expression_data(cell_type):
    """Load gene expression data for specific cell type"""
    # Convert cell_type name to file name format
    # Ex.L2.L3.CUX2.RASGRF2 -> Ex.L2_L3.CUX2_RASGRF2
    # Ex.L5.VAT1L.EYA4 -> Ex.L5.VAT1L_EYA4
    ct_filename = cell_type.replace('.L2.L3.', '.L2_L3.').replace('.L3.L5.', '.L3_L5.')
    ct_filename = ct_filename.replace('.L4.L5.', '.L4_L5.').replace('.L4.L6.', '.L4_L6.')
    ct_filename = ct_filename.replace('.L5.L6.', '.L5_L6.').replace('.VAT1L.', '.VAT1L_')
    ct_filename = ct_filename.replace('.THEMIS.', '.THEMIS_').replace('.FEZF2.', '.FEZF2_')
    ct_filename = ct_filename.replace('.5HT3aR.', '.5HT3aR_').replace('.CDH4.', '.CDH4_')

    expr_file_csv = BASE_DIR / f"motor_cortex_{ct_filename}_expression.csv"
    expr_file_gz = BASE_DIR / f"motor_cortex_{ct_filename}_expression.csv.gz"

    if expr_file_csv.exists():
        expr_file = expr_file_csv
    elif expr_file_gz.exists():
        expr_file = expr_file_gz
    else:
        print(f"  Warning: Expression file not found for {cell_type} (tried {ct_filename})")
        return None

    print(f"  Loading: {expr_file.name}")
    df_expr = pd.read_csv(expr_file, index_col=0)

    # Transpose if needed (genes as columns, cells as rows)
    if df_expr.shape[0] > df_expr.shape[1]:
        df_expr = df_expr.T

    print(f"    Shape: {df_expr.shape} (cells × genes)")

    return df_expr


def compute_tdp_module_scores(df_all):
    """Compute TDP43_targets and STMN2_pathway module scores"""
    print("\n" + "="*60)
    print("Step 1: Computing TDP-43 Module Scores")
    print("="*60)

    # Get TDP gene sets from config
    tdp_genes = ALS_STRESS_COMPONENTS['gene_sets']['TDP43_targets']
    stmn2_genes = ALS_STRESS_COMPONENTS['gene_sets']['STMN2_pathway']

    print(f"TDP43_targets genes: {tdp_genes}")
    print(f"STMN2_pathway genes: {stmn2_genes}")

    # Cell types to process (focus on neurons)
    cell_types_to_process = UPPER_L23 + VAT1L + [
        'Ex.L3.L5.CUX2.RORB',
        'Ex.L4.L5.RORB.FOXO1',
        'Ex.L4.L6.RORB.LRRK1',
        'Ex.L5.L6.THEMIS.DCSTAMP',
        'Ex.L6.FEZF2.SCUBE1',
        'In.5HT3aR.CDH4.SCGN'
    ]

    # Initialize TDP module columns
    df_all['TDP43_targets'] = np.nan
    df_all['STMN2_pathway'] = np.nan

    for ct in cell_types_to_process:
        print(f"\nProcessing: {ct}")

        # Load expression data
        df_expr = load_expression_data(ct)
        if df_expr is None:
            continue

        # Match cells
        ct_cells = df_all[df_all['cell_type'] == ct]['cell_id'].values
        common_cells = np.intersect1d(ct_cells, df_expr.index)

        if len(common_cells) == 0:
            print(f"  Warning: No common cells found")
            continue

        print(f"  Common cells: {len(common_cells)}")

        # Get gene names (try both Symbol and ENSG)
        gene_names = [str(g).upper() for g in df_expr.columns]

        # Compute TDP43_targets score
        tdp_idx = []
        for g in tdp_genes:
            g_upper = g.upper()
            if g_upper in gene_names:
                idx = gene_names.index(g_upper)
                tdp_idx.append(df_expr.columns[idx])

        if len(tdp_idx) > 0:
            print(f"  TDP43_targets: {len(tdp_idx)}/{len(tdp_genes)} genes found")
            # For TDP-43 loss, lower expression = worse pathology
            tdp_score = -df_expr.loc[common_cells, tdp_idx].mean(axis=1)
            df_all.loc[df_all['cell_id'].isin(common_cells), 'TDP43_targets'] = tdp_score.values
        else:
            print(f"  TDP43_targets: No genes found")

        # Compute STMN2_pathway score
        stmn2_idx = []
        for g in stmn2_genes:
            g_upper = g.upper()
            if g_upper in gene_names:
                idx = gene_names.index(g_upper)
                stmn2_idx.append(df_expr.columns[idx])

        if len(stmn2_idx) > 0:
            print(f"  STMN2_pathway: {len(stmn2_idx)}/{len(stmn2_genes)} genes found")
            # For STMN2 loss, lower expression = worse pathology
            stmn2_score = -df_expr.loc[common_cells, stmn2_idx].mean(axis=1)
            df_all.loc[df_all['cell_id'].isin(common_cells), 'STMN2_pathway'] = stmn2_score.values
        else:
            print(f"  STMN2_pathway: No genes found")

    # Summary
    n_tdp = df_all['TDP43_targets'].notna().sum()
    n_stmn2 = df_all['STMN2_pathway'].notna().sum()

    print(f"\nSummary:")
    print(f"  TDP43_targets computed for {n_tdp:,} cells")
    print(f"  STMN2_pathway computed for {n_stmn2:,} cells")

    return df_all


# =====================================================
# Step 2: Define PT_TDP (TDP-43 Pathology Pseudotime)
# =====================================================

def construct_PT_TDP(df_all):
    """Construct PT_TDP from TDP module scores"""
    print("\n" + "="*60)
    print("Step 2: Constructing PT_TDP")
    print("="*60)

    # Combine TDP43_targets and STMN2_pathway
    # Use equal weights if both available, else use whichever exists
    df_all['TDP_module'] = np.nan

    has_both = df_all['TDP43_targets'].notna() & df_all['STMN2_pathway'].notna()
    has_tdp_only = df_all['TDP43_targets'].notna() & df_all['STMN2_pathway'].isna()
    has_stmn2_only = df_all['TDP43_targets'].isna() & df_all['STMN2_pathway'].notna()

    # Combine with equal weights where both available
    df_all.loc[has_both, 'TDP_module'] = (
        0.5 * df_all.loc[has_both, 'TDP43_targets'] +
        0.5 * df_all.loc[has_both, 'STMN2_pathway']
    )
    df_all.loc[has_tdp_only, 'TDP_module'] = df_all.loc[has_tdp_only, 'TDP43_targets']
    df_all.loc[has_stmn2_only, 'TDP_module'] = df_all.loc[has_stmn2_only, 'STMN2_pathway']

    print(f"TDP_module computed for {df_all['TDP_module'].notna().sum():,} cells")

    # Normalize per cell_type using Control cells
    TDP_norm_params = {}

    df_ctrl = df_all[df_all['group'] == 'Control']

    for ct in df_all['cell_type'].unique():
        df_ct_ctrl = df_ctrl[df_ctrl['cell_type'] == ct]

        if len(df_ct_ctrl) >= 30 and df_ct_ctrl['TDP_module'].notna().sum() >= 20:
            mu_ct = df_ct_ctrl['TDP_module'].mean()
            sigma_ct = df_ct_ctrl['TDP_module'].std()
            TDP_norm_params[ct] = (mu_ct, sigma_ct)
        else:
            # Use global Control mean/std
            mu_global = df_ctrl['TDP_module'].mean()
            sigma_global = df_ctrl['TDP_module'].std()
            TDP_norm_params[ct] = (mu_global, sigma_global)

    # Compute Z-score
    df_all['TDP_Z'] = np.nan

    for ct in df_all['cell_type'].unique():
        if ct in TDP_norm_params:
            mu_ct, sigma_ct = TDP_norm_params[ct]
            mask = (df_all['cell_type'] == ct) & df_all['TDP_module'].notna()
            df_all.loc[mask, 'TDP_Z'] = (
                (df_all.loc[mask, 'TDP_module'] - mu_ct) / (sigma_ct + 1e-8)
            )

    # Map Z to PT_TDP [0, 1]
    # Higher TDP_Z (worse pathology) → higher PT_TDP
    df_als = df_all[df_all['group'] == 'ALS']

    # Use only ALS cells to define scaling
    raw_TDP_depth = -df_als['TDP_Z']  # Negate so higher = more TDP loss

    min_depth = raw_TDP_depth.min()
    max_depth = raw_TDP_depth.max()

    # Apply to all cells
    df_all['PT_TDP'] = (-df_all['TDP_Z'] - min_depth) / (max_depth - min_depth + 1e-8)
    df_all['PT_TDP'] = df_all['PT_TDP'].clip(0, 1)

    print(f"\nPT_TDP Statistics:")
    print(f"  Control: mean={df_all[df_all['group']=='Control']['PT_TDP'].mean():.3f}, "
          f"std={df_all[df_all['group']=='Control']['PT_TDP'].std():.3f}")
    print(f"  ALS: mean={df_all[df_all['group']=='ALS']['PT_TDP'].mean():.3f}, "
          f"std={df_all[df_all['group']=='ALS']['PT_TDP'].std():.3f}")

    # QC plot
    fig, ax = plt.subplots(figsize=(8, 5))

    for grp, color in [('Control', 'blue'), ('ALS', 'red')]:
        data = df_all[df_all['group'] == grp]['PT_TDP'].dropna()
        ax.hist(data, bins=50, alpha=0.5, label=f'{grp} (n={len(data):,})', color=color)

    ax.set_xlabel('PT_TDP', fontsize=12)
    ax.set_ylabel('Cell count', fontsize=12)
    ax.set_title('PT_TDP Distribution: Control vs ALS', fontsize=14)
    ax.legend()

    plt.tight_layout()
    plt.savefig(RESULTS_DIR / 'PT_TDP_histogram_Control_vs_ALS.png', dpi=300, bbox_inches='tight')
    plt.close()

    print(f"  Saved: PT_TDP_histogram_Control_vs_ALS.png")

    return df_all


# =====================================================
# Step 3: 2D Phase Space (PT_dpt vs PT_TDP)
# =====================================================

def plot_phase_space_L23_vs_VAT1L(df_all):
    """Create 2D phase space plots for L2/3 vs VAT1L"""
    print("\n" + "="*60)
    print("Step 3: Creating 2D Phase Space (PT_dpt vs PT_TDP)")
    print("="*60)

    df_als = df_all[df_all['group'] == 'ALS'].copy()
    df_L23 = df_als[df_als['cell_type'].isin(UPPER_L23)]
    df_VAT1L = df_als[df_als['cell_type'].isin(VAT1L)]

    # Remove rows with missing PT_dpt or PT_TDP
    df_L23 = df_L23[df_L23['PT_dpt'].notna() & df_L23['PT_TDP'].notna()]
    df_VAT1L = df_VAT1L[df_VAT1L['PT_dpt'].notna() & df_VAT1L['PT_TDP'].notna()]

    print(f"L2/3 cells: {len(df_L23):,}")
    print(f"VAT1L cells: {len(df_VAT1L):,}")

    # 3a. Scatter plot
    fig, ax = plt.subplots(figsize=(10, 8))

    ax.scatter(df_L23['PT_dpt'], df_L23['PT_TDP'],
              s=10, alpha=0.3, c='blue', label=f'L2/3 (n={len(df_L23):,})')
    ax.scatter(df_VAT1L['PT_dpt'], df_VAT1L['PT_TDP'],
              s=10, alpha=0.3, c='red', label=f'VAT1L (n={len(df_VAT1L):,})')

    ax.set_xlabel('PT_dpt (Overall Disease Space)', fontsize=14)
    ax.set_ylabel('PT_TDP (TDP-43 Pathology Depth)', fontsize=14)
    ax.set_title('2D Phase Space: PT_dpt vs PT_TDP (L2/3 vs VAT1L)', fontsize=16)
    ax.legend(fontsize=12)
    ax.grid(alpha=0.3)

    plt.tight_layout()
    plt.savefig(RESULTS_DIR / 'Fig_PT_dpt_vs_PT_TDP_L23_vs_VAT1L_scatter.png',
                dpi=300, bbox_inches='tight')
    plt.close()

    print("  Saved: Fig_PT_dpt_vs_PT_TDP_L23_vs_VAT1L_scatter.png")

    # 3b. Density plots (hexbin)
    for name, df_sub, color in [('L23', df_L23, 'Blues'), ('VAT1L', df_VAT1L, 'Reds')]:
        fig, ax = plt.subplots(figsize=(10, 8))

        hb = ax.hexbin(df_sub['PT_dpt'], df_sub['PT_TDP'],
                      gridsize=30, cmap=color, mincnt=1)

        ax.set_xlabel('PT_dpt (Overall Disease Space)', fontsize=14)
        ax.set_ylabel('PT_TDP (TDP-43 Pathology Depth)', fontsize=14)
        ax.set_title(f'2D Phase Space Density: {name} (n={len(df_sub):,})', fontsize=16)

        cb = plt.colorbar(hb, ax=ax)
        cb.set_label('Cell count', fontsize=12)

        plt.tight_layout()
        plt.savefig(RESULTS_DIR / f'Fig_PT_dpt_vs_PT_TDP_{name}_density.png',
                   dpi=300, bbox_inches='tight')
        plt.close()

    print("  Saved: Fig_PT_dpt_vs_PT_TDP_L23_density.png")
    print("  Saved: Fig_PT_dpt_vs_PT_TDP_VAT1L_density.png")

    # 3c. Quantitative comparison
    stats_list = []

    for name, df_sub in [('L23', df_L23), ('VAT1L', df_VAT1L)]:
        if len(df_sub) < 2:
            print(f"  Warning: {name} has fewer than 2 cells, skipping stats")
            continue

        mean_pt_dpt = df_sub['PT_dpt'].mean()
        std_pt_dpt = df_sub['PT_dpt'].std()
        mean_pt_tdp = df_sub['PT_TDP'].mean()
        std_pt_tdp = df_sub['PT_TDP'].std()

        r_dpt_tdp, p_dpt_tdp = stats.pearsonr(df_sub['PT_dpt'], df_sub['PT_TDP'])

        stats_list.append({
            'cell_group': name,
            'n_cells': len(df_sub),
            'mean_PT_dpt': mean_pt_dpt,
            'std_PT_dpt': std_pt_dpt,
            'mean_PT_TDP': mean_pt_tdp,
            'std_PT_TDP': std_pt_tdp,
            'corr_PT_dpt_PT_TDP': r_dpt_tdp,
            'corr_pval': p_dpt_tdp
        })

    df_stats = pd.DataFrame(stats_list)

    if len(df_stats) > 0:
        df_stats.to_csv(RESULTS_DIR / 'L23_vs_VAT1L_phase_space_stats.csv', index=False)
        print("\nPhase Space Statistics:")
        print(df_stats.to_string(index=False))
        print(f"  Saved: L23_vs_VAT1L_phase_space_stats.csv")
    else:
        print("\n  Warning: No phase space statistics computed (insufficient cells)")

    return df_all


# =====================================================
# Step 4: Module φ vs PT_TDP
# =====================================================

def compute_phi_vs_PT_TDP(df_all):
    """Compute module φ as function of PT_TDP"""
    print("\n" + "="*60)
    print("Step 4: Computing Module φ vs PT_TDP Profiles")
    print("="*60)

    # Get module columns (exclude TDP modules we just created)
    module_cols = [c for c in df_all.columns if c.startswith('module_')]
    print(f"Found {len(module_cols)} modules")

    # Define cell groups
    cell_groups = {
        'L23': UPPER_L23,
        'VAT1L': VAT1L,
        'Glia': GLIA_TYPES,
        'Vascular': VASCULAR_TYPES
    }

    # Get Control mean/std for each module and cell_type
    df_ctrl = df_all[df_all['group'] == 'Control']

    ctrl_stats = {}
    for module in module_cols:
        ctrl_stats[module] = {}
        for ct in df_all['cell_type'].unique():
            df_ct_ctrl = df_ctrl[df_ctrl['cell_type'] == ct]
            if len(df_ct_ctrl) >= 10:
                mu = df_ct_ctrl[module].mean()
                sigma = df_ct_ctrl[module].std()
                ctrl_stats[module][ct] = (mu, sigma)

    # Compute φ (Z²) for each cell
    print("\nComputing φ for all cells...")
    for module in module_cols:
        phi_col = f'phi_{module.replace("module_", "")}'
        df_all[phi_col] = np.nan

        for ct in df_all['cell_type'].unique():
            if ct in ctrl_stats[module]:
                mu, sigma = ctrl_stats[module][ct]
                mask = df_all['cell_type'] == ct
                Z = (df_all.loc[mask, module] - mu) / (sigma + 1e-8)
                df_all.loc[mask, phi_col] = Z ** 2

    phi_cols = [c for c in df_all.columns if c.startswith('phi_')]
    print(f"  Computed {len(phi_cols)} φ columns")

    # Bin PT_TDP and compute mean φ
    df_als = df_all[df_all['group'] == 'ALS'].copy()
    df_als = df_als[df_als['PT_TDP'].notna()]

    bins_TDP = np.linspace(0, 1, 21)
    bin_centers = (bins_TDP[:-1] + bins_TDP[1:]) / 2

    results = []

    for module_base in [m.replace('module_', '') for m in module_cols]:
        phi_col = f'phi_{module_base}'

        if phi_col not in df_als.columns:
            continue

        for group_name, cell_types in cell_groups.items():
            df_group = df_als[df_als['cell_type'].isin(cell_types)]

            if len(df_group) < 50:
                continue

            for i in range(len(bins_TDP) - 1):
                mask = (df_group['PT_TDP'] >= bins_TDP[i]) & (df_group['PT_TDP'] < bins_TDP[i+1])
                df_bin = df_group[mask]

                if len(df_bin) >= 5:
                    phi_mean = df_bin[phi_col].mean()
                    phi_std = df_bin[phi_col].std()
                    n_cells = len(df_bin)

                    results.append({
                        'module': module_base,
                        'cell_group': group_name,
                        'PT_TDP_bin_center': bin_centers[i],
                        'phi_mean': phi_mean,
                        'phi_std': phi_std,
                        'n_cells': n_cells
                    })

    df_profiles = pd.DataFrame(results)
    df_profiles.to_csv(RESULTS_DIR / 'module_phi_vs_PT_TDP_profiles.csv', index=False)

    print(f"\nSaved: module_phi_vs_PT_TDP_profiles.csv ({len(df_profiles)} rows)")

    # Visualization: φ vs PT_TDP for each module (L23 vs VAT1L)
    print("\nGenerating φ vs PT_TDP plots...")

    for module_base in df_profiles['module'].unique():
        df_mod = df_profiles[df_profiles['module'] == module_base]

        fig, ax = plt.subplots(figsize=(10, 6))

        for group_name, color, ls in [('L23', 'blue', '-'), ('VAT1L', 'red', '--')]:
            df_grp = df_mod[df_mod['cell_group'] == group_name]

            if len(df_grp) > 0:
                ax.plot(df_grp['PT_TDP_bin_center'], df_grp['phi_mean'],
                       color=color, linestyle=ls, linewidth=2,
                       marker='o', markersize=4, label=group_name)

        ax.set_xlabel('PT_TDP (TDP-43 Pathology Depth)', fontsize=12)
        ax.set_ylabel(f'φ ({module_base})', fontsize=12)
        ax.set_title(f'Module φ vs PT_TDP: {module_base}', fontsize=14)
        ax.legend()
        ax.grid(alpha=0.3)

        plt.tight_layout()
        plt.savefig(RESULTS_DIR / f'Fig_phi_vs_PT_TDP_{module_base}_L23_vs_VAT1L.png',
                   dpi=200, bbox_inches='tight')
        plt.close()

    print(f"  Saved {len(df_profiles['module'].unique())} φ vs PT_TDP plots")

    return df_all, df_profiles


# =====================================================
# Step 5: Upstream/Downstream Candidate Ranking
# =====================================================

def rank_upstream_downstream_candidates(df_profiles):
    """Rank modules by onset in PT_TDP space"""
    print("\n" + "="*60)
    print("Step 5: Ranking Upstream/Downstream Candidates")
    print("="*60)

    # For each module and cell_group, find onset_PT_TDP and peak_PT_TDP
    phi_threshold = 10.0  # Threshold for "onset"

    results = []

    for module in df_profiles['module'].unique():
        for group in ['L23', 'VAT1L']:
            df_mg = df_profiles[(df_profiles['module'] == module) &
                               (df_profiles['cell_group'] == group)]

            if len(df_mg) == 0:
                continue

            # Onset: first bin where phi_mean > threshold
            onset_bins = df_mg[df_mg['phi_mean'] > phi_threshold]
            if len(onset_bins) > 0:
                onset_PT_TDP = onset_bins['PT_TDP_bin_center'].min()
            else:
                onset_PT_TDP = np.nan

            # Peak: bin with maximum phi_mean
            peak_idx = df_mg['phi_mean'].idxmax()
            peak_PT_TDP = df_mg.loc[peak_idx, 'PT_TDP_bin_center']
            peak_phi = df_mg.loc[peak_idx, 'phi_mean']

            results.append({
                'module': module,
                'cell_group': group,
                'onset_PT_TDP': onset_PT_TDP,
                'peak_PT_TDP': peak_PT_TDP,
                'peak_phi': peak_phi
            })

    df_flow = pd.DataFrame(results)

    # Compute onset ranks (L23 vs VAT1L)
    for module in df_flow['module'].unique():
        df_mod = df_flow[df_flow['module'] == module]

        if len(df_mod) == 2:
            l23_onset = df_mod[df_mod['cell_group'] == 'L23']['onset_PT_TDP'].values[0]
            vat_onset = df_mod[df_mod['cell_group'] == 'VAT1L']['onset_PT_TDP'].values[0]

            if pd.notna(l23_onset) and pd.notna(vat_onset):
                if l23_onset < vat_onset:
                    df_flow.loc[(df_flow['module'] == module) & (df_flow['cell_group'] == 'L23'),
                               'onset_rank_L23_vs_VAT1L'] = 1
                    df_flow.loc[(df_flow['module'] == module) & (df_flow['cell_group'] == 'VAT1L'),
                               'onset_rank_L23_vs_VAT1L'] = 2
                elif vat_onset < l23_onset:
                    df_flow.loc[(df_flow['module'] == module) & (df_flow['cell_group'] == 'L23'),
                               'onset_rank_L23_vs_VAT1L'] = 2
                    df_flow.loc[(df_flow['module'] == module) & (df_flow['cell_group'] == 'VAT1L'),
                               'onset_rank_L23_vs_VAT1L'] = 1

    df_flow.to_csv(RESULTS_DIR / 'module_TDP_flow_summary.csv', index=False)

    print(f"\nSaved: module_TDP_flow_summary.csv ({len(df_flow)} rows)")

    # Summary: Upstream-like (L23 earlier) vs Downstream-like (VAT1L earlier)
    print("\n--- Upstream Candidates (L23 earlier than VAT1L) ---")
    upstream = df_flow[(df_flow['cell_group'] == 'L23') & (df_flow['onset_rank_L23_vs_VAT1L'] == 1)]
    print(upstream[['module', 'onset_PT_TDP', 'peak_phi']].sort_values('onset_PT_TDP'))

    print("\n--- Downstream Candidates (VAT1L earlier than L23) ---")
    downstream = df_flow[(df_flow['cell_group'] == 'VAT1L') & (df_flow['onset_rank_L23_vs_VAT1L'] == 1)]
    print(downstream[['module', 'onset_PT_TDP', 'peak_phi']].sort_values('onset_PT_TDP'))

    return df_flow


# =====================================================
# Main Execution
# =====================================================

def main():
    """Main execution pipeline"""
    print("\n" + "="*60)
    print("Phase TDP: TDP-43 Pathology Pseudotime Analysis")
    print("="*60)

    # Step 0: Load data
    df_all = load_cell_features()

    # Step 1: Compute TDP module scores
    df_all = compute_tdp_module_scores(df_all)

    # Step 2: Construct PT_TDP
    df_all = construct_PT_TDP(df_all)

    # Step 3: 2D Phase space
    df_all = plot_phase_space_L23_vs_VAT1L(df_all)

    # Step 4: Module φ vs PT_TDP
    df_all, df_profiles = compute_phi_vs_PT_TDP(df_all)

    # Step 5: Rank upstream/downstream
    df_flow = rank_upstream_downstream_candidates(df_profiles)

    # Save updated cell-level features with PT_TDP
    output_path = RESULTS_DIR / 'cell_level_features_with_PT_TDP.csv'
    df_all.to_csv(output_path, index=False)
    print(f"\nSaved: {output_path.name}")

    print("\n" + "="*60)
    print("Phase TDP Analysis Complete!")
    print(f"Results saved to: {RESULTS_DIR}")
    print("="*60)


if __name__ == '__main__':
    main()
