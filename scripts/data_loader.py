"""
scRNAseq Hodge Pipeline — Data Loading & Filtering
====================================================
Loads h5ad files, filters cells by cell type,
and provides per-donor cell count tables and sample metadata.
"""
import gc
import logging
from pathlib import Path
from typing import Dict, Optional, Tuple

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

from . import config

logger = logging.getLogger(__name__)


# ── Gene annotation (mito gene identification) ────────────────
_GENE_ANNOTATION_CACHE: Optional[pd.DataFrame] = None


def _load_gene_annotation() -> Optional[pd.DataFrame]:
    """Load and cache gene annotation with chromosome info."""
    global _GENE_ANNOTATION_CACHE
    if _GENE_ANNOTATION_CACHE is not None:
        return _GENE_ANNOTATION_CACHE

    if config.GENE_ANNOTATION_PATH and config.GENE_ANNOTATION_PATH.exists():
        annot = pd.read_csv(config.GENE_ANNOTATION_PATH)
        annot = annot.set_index("gene_id")
        _GENE_ANNOTATION_CACHE = annot
        return annot
    return None


def get_mito_gene_ids(var_names: pd.Index) -> np.ndarray:
    """Return boolean mask for mitochondrial genes among var_names.

    Uses gene_annotation.csv if available (chrM chromosome),
    otherwise falls back to gene names starting with 'MT-'.
    """
    annot = _load_gene_annotation()

    if annot is not None and "chromosome" in annot.columns:
        mt_ids = set(annot.index[annot["chromosome"] == "chrM"])
        bare = [v.split(".")[0] for v in var_names]
        return np.array([b in mt_ids for b in bare], dtype=bool)
    else:
        # Fallback: gene name starts with "MT-" or "mt-"
        return np.array(
            [str(v).upper().startswith("MT-") for v in var_names],
            dtype=bool,
        )


# ── Sample metadata ───────────────────────────────────────────

def get_sample_info() -> pd.DataFrame:
    """Load sample_info.csv with donor_id as string."""
    path = config.SAMPLE_INFO_PATH
    if not path or not path.exists():
        raise FileNotFoundError(
            f"sample_info.csv not found at {path}. "
            "Please create data/metadata/sample_info.csv with columns: "
            "donor_id, file, condition"
        )
    df = pd.read_csv(path)
    if "donor_id" not in df.columns:
        raise ValueError("sample_info.csv must have a 'donor_id' column")
    if "file" not in df.columns:
        raise ValueError("sample_info.csv must have a 'file' column")
    df["donor_id"] = df["donor_id"].astype(str)
    return df


# ── Per-sample loading ─────────────────────────────────────────

def load_sample(filename: str, backed: Optional[str] = None) -> ad.AnnData:
    """Load a single h5ad file."""
    path = config.H5AD_DIR / filename
    if not path.exists():
        raise FileNotFoundError(f"h5ad file not found: {path}")
    return ad.read_h5ad(path, backed=backed)


def load_all_samples(backed: Optional[str] = None) -> Dict[str, ad.AnnData]:
    """Load all h5ad files into a dict keyed by donor_id (string)."""
    info = get_sample_info()
    result = {}
    for _, row in info.iterrows():
        donor_id = str(row["donor_id"])
        fname = row["file"]
        logger.info("Loading %s (donor %s)", fname, donor_id)
        try:
            result[donor_id] = load_sample(fname, backed=backed)
        except Exception as e:
            logger.error("Failed to load %s: %s", fname, e)
    logger.info("Loaded %d / %d samples", len(result), len(info))
    return result


# ── Cell-type filtering ───────────────────────────────────────

def get_cells_for_celltype(
    adata_dict: Dict[str, ad.AnnData],
    celltype: str,
    min_cells_per_donor: int = None,
) -> ad.AnnData:
    """Extract and combine cells of a given cell type from all donors.

    For each donor, selects cells where obs[CELLTYPE_COLUMN] == celltype.
    Donors with fewer than min_cells_per_donor cells are excluded.
    Adds obs columns: 'donor_id', 'nUMI', 'pct_mito'.

    Parameters
    ----------
    adata_dict : dict
        {donor_id_str: AnnData} as returned by load_all_samples().
    celltype : str
        Cell type label matching obs[CELLTYPE_COLUMN].
    min_cells_per_donor : int, optional
        Minimum cells required per donor. Defaults to config.MIN_CELLS_PER_DONOR.

    Returns
    -------
    AnnData
        Combined AnnData with cells from all qualifying donors.
    """
    if min_cells_per_donor is None:
        min_cells_per_donor = config.MIN_CELLS_PER_DONOR

    ct_col = config.CELLTYPE_COLUMN
    pieces = []
    included_donors = []
    excluded_donors = []

    for donor_id, adata in sorted(adata_dict.items()):
        if ct_col not in adata.obs.columns:
            logger.warning(
                "Donor %s: column '%s' not found in obs. Available: %s",
                donor_id, ct_col, list(adata.obs.columns),
            )
            continue

        mask = adata.obs[ct_col] == celltype
        n_cells = int(mask.sum())

        if n_cells < min_cells_per_donor:
            excluded_donors.append((donor_id, n_cells))
            continue

        sub = adata[mask].copy()

        # Densify X
        if sp.issparse(sub.X):
            sub.X = np.asarray(sub.X.todense(), dtype=np.float32)
        else:
            sub.X = np.asarray(sub.X, dtype=np.float32)

        # Compute QC metrics
        nUMI = sub.X.sum(axis=1)
        sub.obs["nUMI"] = nUMI.astype(np.float32)

        # pct_mito
        mito_mask = get_mito_gene_ids(sub.var_names)
        if mito_mask.any():
            mito_counts = sub.X[:, mito_mask].sum(axis=1)
            safe_total = np.where(nUMI > 0, nUMI, 1.0)
            sub.obs["pct_mito"] = (mito_counts / safe_total).astype(np.float32)
        else:
            sub.obs["pct_mito"] = np.float32(0.0)

        sub.obs["donor_id"] = str(donor_id)

        # Carry over condition metadata
        cond_col = config.CONDITION_COLUMN
        if cond_col in sub.obs.columns:
            sub.obs["condition"] = sub.obs[cond_col].values
        elif "Condition" in sub.obs.columns:
            sub.obs["condition"] = sub.obs["Condition"].values

        pieces.append(sub)
        included_donors.append(donor_id)

    if not pieces:
        raise ValueError(
            f"No donors have >= {min_cells_per_donor} cells for celltype '{celltype}'. "
            f"Excluded: {excluded_donors}"
        )

    logger.info(
        "Celltype %s: %d donors included (%d excluded), %d total cells",
        celltype,
        len(included_donors),
        len(excluded_donors),
        sum(p.n_obs for p in pieces),
    )

    combined = ad.concat(pieces, merge="same")
    combined.obs["donor_id"] = combined.obs["donor_id"].astype(str)
    combined.obs["nUMI"] = combined.obs["nUMI"].astype(np.float32)
    combined.obs["pct_mito"] = combined.obs["pct_mito"].astype(np.float32)

    if sp.issparse(combined.X):
        combined.X = np.asarray(combined.X.todense(), dtype=np.float32)
    else:
        combined.X = np.asarray(combined.X, dtype=np.float32)

    return combined
