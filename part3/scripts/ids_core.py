#!/usr/bin/env python3
"""
IDS Core Library (v1.2)
=======================
Common core functions for Integrated Disease Space analysis.

This library provides reusable functions for:
- φ (phi) energy computation from Control-normalized scores
- Pseudotime (PT) binning and aggregation
- onset/peak detection along PT trajectories
- Hierarchical aggregation across cell types/modules
- Axis Discovery: Automatic detection of independent disease axes

Originally extracted from ALS motor cortex project (Phase 13, TDP-43/STMN2 analyses).

Author: Claude Code + User
Date: 2025-11-27 (Phase IDS-Core v1)
Date: 2025-11-30 (Phase IDS-Core v1.2 - Axis Discovery Engine)
"""

import numpy as np
import pandas as pd
from typing import Optional, Union, Tuple, List, Dict
import warnings

# Axis Discovery dependencies (optional, imported only when needed)
try:
    from sklearn.cluster import KMeans, AgglomerativeClustering
    from sklearn.decomposition import PCA, NMF
    from sklearn.preprocessing import StandardScaler
    from sklearn.metrics import silhouette_score
    SKLEARN_AVAILABLE = True
except ImportError:
    SKLEARN_AVAILABLE = False


# ============================================================================
# Core Function 1: φ (phi) Energy Computation
# ============================================================================

def compute_phi_from_scores(
    scores: pd.Series,
    control_mask: pd.Series,
    eps: float = 1e-8,
    return_z: bool = False,
) -> Union[pd.Series, Tuple[pd.Series, pd.Series]]:
    """
    Compute φ (phi) energy = Z² from Control-normalized scores.

    φ represents the squared deviation from Control group statistics,
    quantifying the "energy" or "abnormality" of each cell.

    Parameters
    ----------
    scores : pd.Series
        Module/gene scores for all cells (index = cell_id or similar)
    control_mask : pd.Series
        Boolean mask indicating Control cells (same index as scores)
    eps : float, default=1e-8
        Small epsilon to avoid division by zero
    return_z : bool, default=False
        If True, also return Z-scores (not just Z²)

    Returns
    -------
    phi : pd.Series
        φ = Z² values (same index as scores)
    z : pd.Series (optional, if return_z=True)
        Z-scores

    Examples
    --------
    >>> scores = pd.Series([1.0, 2.0, 3.0, 4.0, 5.0], index=['c1', 'c2', 'c3', 'c4', 'c5'])
    >>> control_mask = pd.Series([True, True, False, False, False], index=['c1', 'c2', 'c3', 'c4', 'c5'])
    >>> phi = compute_phi_from_scores(scores, control_mask)
    >>> print(phi)

    Notes
    -----
    - Based on Phase_STMN2_single_gene_analysis_v2.py lines 165-167
    - Formula: Z = (score - μ_Control) / σ_Control, φ = Z²
    - φ loses sign information (both up- and down-regulation increase φ)
    """
    # Validate inputs
    if not scores.index.equals(control_mask.index):
        raise ValueError("scores and control_mask must have the same index")

    # Convert to float
    scores = scores.astype(float)

    # Extract Control scores
    ctrl_scores = scores[control_mask]

    if len(ctrl_scores) == 0:
        warnings.warn("No Control cells found. Returning NaN for phi.", UserWarning)
        return pd.Series(np.nan, index=scores.index)

    # Compute Control statistics
    mu = ctrl_scores.mean()
    sigma = ctrl_scores.std(ddof=1)  # Use sample std (ddof=1) to match pandas default

    if sigma < eps:
        warnings.warn(f"Control std is near zero (σ={sigma:.2e}). Adding eps.", UserWarning)

    # Compute Z-score
    z = (scores - mu) / (sigma + eps)

    # Compute φ = Z²
    phi = z ** 2

    if return_z:
        return phi, z
    else:
        return phi


# ============================================================================
# Core Function 2: Pseudotime (PT) Binning
# ============================================================================

def bin_by_pt(
    pt: pd.Series,
    n_bins: int = 25,
) -> Tuple[pd.Series, np.ndarray, np.ndarray]:
    """
    Bin cells by pseudotime (PT) into equal-width bins.

    Parameters
    ----------
    pt : pd.Series
        Pseudotime values (index = cell_id or similar)
    n_bins : int, default=25
        Number of bins

    Returns
    -------
    bin_index : pd.Series
        Bin assignment for each cell (0 to n_bins-1), same index as pt
    bin_edges : np.ndarray
        Bin edges (length = n_bins + 1)
    bin_centers : np.ndarray
        Bin center values (length = n_bins)

    Examples
    --------
    >>> pt = pd.Series([0.1, 0.2, 0.3, 0.4, 0.5], index=['c1', 'c2', 'c3', 'c4', 'c5'])
    >>> bin_index, bin_edges, bin_centers = bin_by_pt(pt, n_bins=5)
    >>> print(bin_index)

    Notes
    -----
    - Based on Phase13c_NVU_energy_complete.py lines 146-147
    - Uses np.linspace for equal-width binning
    - Bins are 0-indexed (0 to n_bins-1)
    """
    pt = pt.astype(float).dropna()

    if len(pt) == 0:
        raise ValueError("PT series is empty after dropna()")

    pt_min = pt.min()
    pt_max = pt.max()

    # Create bin edges
    edges = np.linspace(pt_min, pt_max, n_bins + 1)

    # Assign bins using pd.cut (returns categorical, convert to int)
    bin_cat = pd.cut(pt, bins=edges, labels=False, include_lowest=True)
    bin_index = pd.Series(bin_cat, index=pt.index, dtype='Int64')  # nullable int

    # Compute bin centers
    centers = (edges[:-1] + edges[1:]) / 2

    return bin_index, edges, centers


# ============================================================================
# Core Function 3: Aggregate φ by Group and Bin
# ============================================================================

def summarize_phi_by_group_and_bin(
    phi: pd.Series,
    group: pd.Series,
    bin_index: pd.Series,
    group_order: Optional[List[str]] = None,
) -> pd.DataFrame:
    """
    Aggregate φ values by [group × bin].

    Typically used to compute mean φ for each (cell_type, PT_bin) combination.

    Parameters
    ----------
    phi : pd.Series
        φ values (index = cell_id or similar)
    group : pd.Series
        Group labels (e.g., cell_type, NVU node name), same index as phi
    bin_index : pd.Series
        Bin indices (0 to n_bins-1), same index as phi
    group_order : list of str, optional
        Desired order of groups in output (for categorical sorting)

    Returns
    -------
    summary : pd.DataFrame
        Columns: ['group', 'bin', 'phi_mean', 'phi_std', 'n_cells']
        One row per (group, bin) combination

    Examples
    --------
    >>> phi = pd.Series([1.0, 2.0, 3.0, 4.0], index=['c1', 'c2', 'c3', 'c4'])
    >>> group = pd.Series(['A', 'A', 'B', 'B'], index=['c1', 'c2', 'c3', 'c4'])
    >>> bin_index = pd.Series([0, 1, 0, 1], index=['c1', 'c2', 'c3', 'c4'])
    >>> summary = summarize_phi_by_group_and_bin(phi, group, bin_index)
    >>> print(summary)

    Notes
    -----
    - Based on Phase13c_NVU_energy_complete.py lines 154
    - Uses groupby(['group', 'bin']).agg(['mean', 'std', 'count'])
    """
    # Validate inputs
    if not (phi.index.equals(group.index) and phi.index.equals(bin_index.index)):
        raise ValueError("phi, group, and bin_index must have the same index")

    # Create DataFrame
    df = pd.DataFrame({
        'phi': phi,
        'group': group,
        'bin': bin_index,
    }).dropna()

    # Group and aggregate
    gb = df.groupby(['group', 'bin'])
    summary = gb['phi'].agg(['mean', 'std', 'count']).reset_index()
    summary = summary.rename(columns={'mean': 'phi_mean', 'std': 'phi_std', 'count': 'n_cells'})

    # Optional: enforce group order
    if group_order is not None:
        summary['group'] = pd.Categorical(summary['group'], categories=group_order, ordered=True)
        summary = summary.sort_values(['group', 'bin']).reset_index(drop=True)

    return summary


# ============================================================================
# Core Function 4: Detect Onset and Peak
# ============================================================================

def detect_onset_and_peak(
    phi_mean_df: pd.DataFrame,
    bin_centers: Optional[np.ndarray] = None,
    threshold: float = 0.5,
    min_bins: int = 1,
) -> pd.DataFrame:
    """
    Detect onset_PT and peak_PT for each group from aggregated φ profiles.

    Onset: First bin where φ_mean > threshold
    Peak: Bin with maximum φ_mean

    Parameters
    ----------
    phi_mean_df : pd.DataFrame
        Output of summarize_phi_by_group_and_bin()
        Must have columns: ['group', 'bin', 'phi_mean']
    bin_centers : np.ndarray, optional
        Bin center values (length = n_bins) for mapping bin index to PT values.
        If None, only bin indices are returned (not PT values).
    threshold : float, default=0.5
        φ threshold for onset detection
    min_bins : int, default=1
        Minimum number of bins required for onset detection (optional constraint)

    Returns
    -------
    flow_summary : pd.DataFrame
        Columns: ['group', 'onset_bin', 'onset_PT', 'peak_bin', 'peak_PT', 'peak_phi']
        One row per group

    Examples
    --------
    >>> phi_mean_df = pd.DataFrame({
    ...     'group': ['A', 'A', 'B', 'B'],
    ...     'bin': [0, 1, 0, 1],
    ...     'phi_mean': [0.3, 0.6, 0.7, 0.5]
    ... })
    >>> bin_centers = np.array([0.25, 0.75])
    >>> flow_summary = detect_onset_and_peak(phi_mean_df, bin_centers, threshold=0.5)
    >>> print(flow_summary)

    Notes
    -----
    - Based on Phase13c_NVU_energy_complete.py lines 188-197
    - Onset: first bin where φ_mean > threshold
    - Peak: bin with maximum φ_mean (argmax)
    """
    results = []

    for group, sub in phi_mean_df.groupby('group'):
        sub = sub.sort_values('bin').reset_index(drop=True)

        # Peak: bin with max φ_mean
        idx_peak = sub['phi_mean'].idxmax()
        peak_bin = int(sub.loc[idx_peak, 'bin'])
        peak_phi = float(sub.loc[idx_peak, 'phi_mean'])

        # Map bin to PT (if bin_centers provided)
        if bin_centers is not None:
            peak_PT = float(bin_centers[peak_bin])
        else:
            peak_PT = None

        # Onset: first bin where φ_mean > threshold
        onset_candidates = sub[sub['phi_mean'] > threshold]

        if len(onset_candidates) >= min_bins:
            onset_bin = int(onset_candidates.iloc[0]['bin'])
            if bin_centers is not None:
                onset_PT = float(bin_centers[onset_bin])
            else:
                onset_PT = None
        else:
            onset_bin = None
            onset_PT = None

        results.append({
            'group': group,
            'onset_bin': onset_bin,
            'onset_PT': onset_PT,
            'peak_bin': peak_bin,
            'peak_PT': peak_PT,
            'peak_phi': peak_phi,
        })

    return pd.DataFrame(results)


# ============================================================================
# Utility Function: Map Cell Types to Hierarchical Groups
# ============================================================================

def map_cell_types_to_groups(
    cell_types: pd.Series,
    group_mapping: dict,
    default_group: str = 'Other'
) -> pd.Series:
    """
    Map cell types to higher-level groups (e.g., NVU nodes, broad cell classes).

    Parameters
    ----------
    cell_types : pd.Series
        Cell type labels (index = cell_id or similar)
    group_mapping : dict
        Mapping from group name to list of cell type patterns
        Example: {'Vascular': ['Vasc.Endo', 'Vasc.Fibro'], 'Glia': ['Glia.Astro', 'Glia.Oligo']}
    default_group : str, default='Other'
        Group label for cells that don't match any pattern

    Returns
    -------
    groups : pd.Series
        Group labels (same index as cell_types)

    Examples
    --------
    >>> cell_types = pd.Series(['Vasc.Endo', 'Glia.Astro', 'Unknown'], index=['c1', 'c2', 'c3'])
    >>> group_mapping = {'Vascular': ['Vasc.Endo'], 'Glia': ['Glia.Astro']}
    >>> groups = map_cell_types_to_groups(cell_types, group_mapping)
    >>> print(groups)

    Notes
    -----
    - Based on Phase13_complete_module_analysis.py lines 95-102
    - Uses substring matching (any pattern in cell_type string)
    """
    def map_to_group(cell_type):
        for group_name, patterns in group_mapping.items():
            if any(pattern in cell_type for pattern in patterns):
                return group_name
        return default_group

    groups = cell_types.apply(map_to_group)
    return groups


# ============================================================================
# Convenience Wrapper: One-Call φ-Flow Pipeline
# ============================================================================

def compute_phi_flow(
    scores: pd.Series,
    control_mask: pd.Series,
    pt: pd.Series,
    group: pd.Series,
    n_bins: int = 25,
    threshold: float = 0.5,
    group_order: Optional[List[str]] = None,
    eps: float = 1e-8,
) -> Tuple[pd.DataFrame, pd.DataFrame, np.ndarray]:
    """
    Convenience wrapper: Compute full φ-flow analysis in one call.

    This function combines all core IDS operations:
      1. Compute φ = Z² from Control-normalized scores
      2. Bin cells by pseudotime (PT)
      3. Aggregate φ by [group × bin]
      4. Detect onset and peak for each group

    Parameters
    ----------
    scores : pd.Series
        Module/gene scores for all cells (index = cell_id)
    control_mask : pd.Series
        Boolean mask indicating Control cells (same index as scores)
    pt : pd.Series
        Pseudotime values (same index as scores)
    group : pd.Series
        Group labels (e.g., cell_type, NVU node), same index as scores
    n_bins : int, default=25
        Number of PT bins
    threshold : float, default=0.5
        φ threshold for onset detection
    group_order : list of str, optional
        Desired order of groups in output
    eps : float, default=1e-8
        Small epsilon for φ computation (avoid division by zero)

    Returns
    -------
    profiles : pd.DataFrame
        Aggregated φ profiles by [group × bin]
        Columns: ['group', 'bin', 'phi_mean', 'phi_std', 'n_cells']
    flow : pd.DataFrame
        Flow summary for each group
        Columns: ['group', 'onset_bin', 'onset_PT', 'peak_bin', 'peak_PT', 'peak_phi']
    bin_centers : np.ndarray
        PT bin center values (length = n_bins)

    Examples
    --------
    >>> # Compute φ-flow for a single module across NVU nodes
    >>> profiles, flow, centers = compute_phi_flow(
    ...     scores=df['module_Metabolism'],
    ...     control_mask=(df['condition'] == 'Control'),
    ...     pt=df['PT_dpt'],
    ...     group=df['nvu_node'],
    ...     n_bins=30,
    ...     threshold=0.5,
    ...     group_order=['Vascular', 'Glia', 'Upper', 'VAT1L']
    ... )

    Notes
    -----
    - This is a high-level convenience function for common φ-flow workflows
    - For more control over individual steps, use the core functions directly:
      compute_phi_from_scores() → bin_by_pt() → summarize_phi_by_group_and_bin() → detect_onset_and_peak()
    - All input Series must have the same index (typically cell_id)
    """
    # Validate inputs
    if not (scores.index.equals(control_mask.index) and
            scores.index.equals(pt.index) and
            scores.index.equals(group.index)):
        raise ValueError("All input Series must have the same index")

    # Step 1: Compute φ from Control-normalized scores
    phi = compute_phi_from_scores(scores, control_mask, eps=eps)

    # Step 2: Bin by pseudotime
    bin_index, bin_edges, bin_centers = bin_by_pt(pt, n_bins=n_bins)

    # Step 3: Aggregate φ by [group × bin]
    profiles = summarize_phi_by_group_and_bin(
        phi=phi,
        group=group,
        bin_index=bin_index,
        group_order=group_order
    )

    # Step 4: Detect onset and peak
    phi_mean_df = profiles[['group', 'bin', 'phi_mean']].copy()
    flow = detect_onset_and_peak(
        phi_mean_df=phi_mean_df,
        bin_centers=bin_centers,
        threshold=threshold,
        min_bins=1
    )

    return profiles, flow, bin_centers


# ============================================================================
# Axis Discovery Engine (v1.2)
# ============================================================================

def build_module_feature_matrix(
    phi_profiles: pd.DataFrame,
    group_col: str = 'node',
    module_col: str = 'module',
    bin_col: str = 'pt_bin',
    value_col: str = 'phi_mean',
    normalize: bool = True,
) -> Tuple[np.ndarray, pd.Index, List[str]]:
    """
    Build module × feature matrix for axis discovery.

    Converts φ-flow profiles (long format) into a matrix where:
    - Rows = modules
    - Columns = (group, bin) combinations

    This matrix can be used for clustering or factorization to discover
    independent disease axes.

    Parameters
    ----------
    phi_profiles : pd.DataFrame
        φ-flow profiles in long format (output of summarize_phi_by_group_and_bin)
        Must have columns: [group_col, module_col, bin_col, value_col]
    group_col : str, default='node'
        Column name for groups (e.g., 'node', 'cell_type')
    module_col : str, default='module'
        Column name for modules
    bin_col : str, default='pt_bin'
        Column name for PT bins
    value_col : str, default='phi_mean'
        Column name for φ values
    normalize : bool, default=True
        If True, standardize features (mean=0, std=1) across modules

    Returns
    -------
    X : np.ndarray
        Feature matrix, shape [n_modules, n_features]
        where n_features = n_groups × n_bins
    module_names : pd.Index
        Module names corresponding to rows
    feature_names : List[str]
        Feature names (e.g., ['Vascular_bin0', 'Vascular_bin1', ...])

    Examples
    --------
    >>> # Load φ-flow profiles
    >>> profiles = pd.read_csv('nvu_phi_profiles_final_by_node_and_bin.csv')
    >>> X, modules, features = build_module_feature_matrix(profiles)
    >>> print(X.shape)  # (n_modules, n_groups * n_bins)

    Notes
    -----
    - Missing (group, bin) combinations are filled with 0
    - If normalize=True, StandardScaler is applied (requires sklearn)
    """
    # Pivot to wide format: modules × (group, bin)
    df_pivot = phi_profiles.pivot_table(
        index=module_col,
        columns=[group_col, bin_col],
        values=value_col,
        fill_value=0.0  # Fill missing combinations with 0
    )

    # Extract module names and feature names
    module_names = df_pivot.index
    feature_names = [f"{grp}_bin{bin}" for grp, bin in df_pivot.columns]

    # Convert to numpy array
    X = df_pivot.values

    # Normalize if requested
    if normalize:
        if not SKLEARN_AVAILABLE:
            warnings.warn("sklearn not available. Skipping normalization.", UserWarning)
        else:
            scaler = StandardScaler()
            X = scaler.fit_transform(X)

    return X, module_names, feature_names


def discover_axes_clustering(
    X: np.ndarray,
    module_names: pd.Index,
    n_axes: int = 3,
    method: str = 'kmeans',
    random_state: int = 42,
) -> pd.DataFrame:
    """
    Discover disease axes using clustering methods.

    Groups modules with similar φ-flow patterns into clusters (axes).

    Parameters
    ----------
    X : np.ndarray
        Feature matrix from build_module_feature_matrix()
        Shape: [n_modules, n_features]
    module_names : pd.Index
        Module names corresponding to rows of X
    n_axes : int, default=3
        Number of axes (clusters) to discover
    method : str, default='kmeans'
        Clustering method: 'kmeans' or 'hierarchical'
    random_state : int, default=42
        Random seed for reproducibility

    Returns
    -------
    axis_assignment : pd.DataFrame
        Columns: ['module', 'axis', 'axis_name']
        One row per module

    Examples
    --------
    >>> X, modules, features = build_module_feature_matrix(profiles)
    >>> axes = discover_axes_clustering(X, modules, n_axes=3)
    >>> print(axes)

    Notes
    -----
    - Axis names are automatically assigned as 'Axis_1', 'Axis_2', etc.
    - User should rename based on biological interpretation
    - Requires sklearn
    """
    if not SKLEARN_AVAILABLE:
        raise ImportError("sklearn is required for axis discovery. Install with: pip install scikit-learn")

    # Perform clustering
    if method == 'kmeans':
        clusterer = KMeans(n_clusters=n_axes, random_state=random_state, n_init=10)
    elif method == 'hierarchical':
        clusterer = AgglomerativeClustering(n_clusters=n_axes)
    else:
        raise ValueError(f"Unknown method: {method}. Use 'kmeans' or 'hierarchical'.")

    # Fit and predict
    cluster_labels = clusterer.fit_predict(X)

    # Compute silhouette score
    if len(np.unique(cluster_labels)) > 1:
        sil_score = silhouette_score(X, cluster_labels)
    else:
        sil_score = np.nan

    # Create result DataFrame
    result = pd.DataFrame({
        'module': module_names,
        'axis': cluster_labels,
        'axis_name': [f'Axis_{i+1}' for i in cluster_labels]
    })

    # Add silhouette score as metadata (store in a separate attribute if needed)
    result.attrs['silhouette_score'] = sil_score
    result.attrs['n_axes'] = n_axes
    result.attrs['method'] = method

    return result


def discover_axes_factorization(
    X: np.ndarray,
    module_names: pd.Index,
    feature_names: List[str],
    n_axes: int = 3,
    method: str = 'nmf',
    random_state: int = 42,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Discover disease axes using matrix factorization methods.

    Decomposes module × feature matrix into:
    - W: module × axis loadings
    - H: axis × feature patterns

    Parameters
    ----------
    X : np.ndarray
        Feature matrix from build_module_feature_matrix()
        Shape: [n_modules, n_features]
    module_names : pd.Index
        Module names corresponding to rows of X
    feature_names : List[str]
        Feature names corresponding to columns of X
    n_axes : int, default=3
        Number of latent axes to extract
    method : str, default='nmf'
        Factorization method: 'nmf' (Non-negative Matrix Factorization) or 'pca'
    random_state : int, default=42
        Random seed for reproducibility

    Returns
    -------
    W : pd.DataFrame
        Module loadings on each axis
        Columns: ['module', 'Axis_1', 'Axis_2', ...]
        Higher values = module strongly associated with that axis
    H : pd.DataFrame
        Axis patterns (φ-flow signatures)
        Columns: feature_names
        Index: ['Axis_1', 'Axis_2', ...]

    Examples
    --------
    >>> X, modules, features = build_module_feature_matrix(profiles)
    >>> W, H = discover_axes_factorization(X, modules, features, n_axes=3, method='nmf')
    >>> # Find top modules for Axis_1
    >>> print(W.nlargest(5, 'Axis_1'))

    Notes
    -----
    - NMF requires non-negative input (φ is always ≥ 0, so this is fine)
    - PCA gives orthogonal components but may have negative values
    - For biological interpretation, NMF is often preferred
    - Requires sklearn
    """
    if not SKLEARN_AVAILABLE:
        raise ImportError("sklearn is required for axis discovery. Install with: pip install scikit-learn")

    # Perform factorization
    if method == 'nmf':
        # Ensure non-negative (φ should already be non-negative, but just in case)
        if np.any(X < 0):
            warnings.warn("NMF requires non-negative values. Shifting X by min value.", UserWarning)
            X = X - X.min()

        model = NMF(n_components=n_axes, random_state=random_state, init='nndsvda', max_iter=500)
    elif method == 'pca':
        model = PCA(n_components=n_axes, random_state=random_state)
    else:
        raise ValueError(f"Unknown method: {method}. Use 'nmf' or 'pca'.")

    # Fit and transform
    W_matrix = model.fit_transform(X)  # modules × axes
    H_matrix = model.components_  # axes × features

    # Convert to DataFrames
    axis_names = [f'Axis_{i+1}' for i in range(n_axes)]

    W = pd.DataFrame(W_matrix, columns=axis_names)
    W.insert(0, 'module', module_names)

    H = pd.DataFrame(H_matrix, columns=feature_names, index=axis_names)

    # Add explained variance (for PCA) or reconstruction error (for NMF)
    if method == 'pca':
        W.attrs['explained_variance_ratio'] = model.explained_variance_ratio_
    elif method == 'nmf':
        W.attrs['reconstruction_error'] = model.reconstruction_err_

    return W, H


def reconstruct_axis_profile(
    H: pd.DataFrame,
    axis_name: str,
    group_order: Optional[List[str]] = None,
) -> pd.DataFrame:
    """
    Reconstruct φ-flow profile for a discovered axis.

    Converts axis × feature vector back to (group, bin, phi) format
    for visualization.

    Parameters
    ----------
    H : pd.DataFrame
        Axis patterns from discover_axes_factorization()
        Columns: feature_names (e.g., 'Vascular_bin0', 'Glia_bin5')
        Index: axis names
    axis_name : str
        Name of axis to reconstruct (e.g., 'Axis_1')
    group_order : List[str], optional
        Desired order of groups for plotting

    Returns
    -------
    profile : pd.DataFrame
        Columns: ['group', 'bin', 'phi_pattern']
        Ready for plotting with standard φ-flow visualization tools

    Examples
    --------
    >>> W, H = discover_axes_factorization(X, modules, features, n_axes=3)
    >>> axis1_profile = reconstruct_axis_profile(H, 'Axis_1')
    >>> # Plot this profile to see the "canonical" φ-flow for Axis 1

    Notes
    -----
    - The reconstructed profile represents the "canonical pattern" for this axis
    - Modules with high loading on this axis will resemble this pattern
    """
    # Extract axis pattern
    axis_pattern = H.loc[axis_name]

    # Parse feature names back to (group, bin)
    records = []
    for feature_name, value in axis_pattern.items():
        # feature_name format: 'GroupName_binN'
        parts = feature_name.rsplit('_bin', 1)
        if len(parts) == 2:
            group = parts[0]
            bin_num = int(parts[1])
            records.append({
                'group': group,
                'bin': bin_num,
                'phi_pattern': value
            })

    profile = pd.DataFrame(records)

    # Optional: enforce group order
    if group_order is not None:
        profile['group'] = pd.Categorical(profile['group'], categories=group_order, ordered=True)
        profile = profile.sort_values(['group', 'bin']).reset_index(drop=True)

    return profile


# ============================================================================
# End of ids_core.py v1.2
# ============================================================================
