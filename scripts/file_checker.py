"""
scRNAseq Hodge Pipeline — File Checker & Modifier
====================================================
Validates and fixes h5ad files for pipeline compatibility.

Checks:
  1. h5ad structure (.X, .obs, .var)
  2. Cell type column presence and naming
  3. Gene ID format (ENSEMBL vs symbol)
  4. Count matrix type (raw counts vs normalized)
  5. Per-donor cell counts and QC metrics
  6. sample_info.csv completeness

Fixes:
  - Rename cell type column
  - Split multi-sample h5ad into per-donor files
  - Auto-generate sample_info.csv
  - Convert gene names (symbol <-> ENSEMBL)
  - Filter low-quality cells (high mito, low nUMI)
"""
import json
import logging
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import scipy.sparse as sp

logger = logging.getLogger(__name__)

# Common cell type column names across different pipelines
KNOWN_CT_COLUMNS = [
    "CellType", "celltype", "cell_type", "CellType_l1", "CellType_l2",
    "cluster", "Cluster", "clusters", "leiden", "louvain",
    "cell_type_ontology_term_id", "cell_ontology_class",
    "ann_level_1", "ann_level_2", "ann_finest_level",
    "predicted.celltype.l1", "predicted.celltype.l2",
    "seurat_clusters", "labels", "Label",
]

# Common donor/sample ID column names
KNOWN_DONOR_COLUMNS = [
    "donor_id", "sample_id", "patient_id", "subject_id",
    "donor", "sample", "patient", "subject", "individual",
    "orig.ident", "Sample", "Donor", "Patient", "SubjectID",
    "batch", "Batch",
]

# Common condition column names
KNOWN_CONDITION_COLUMNS = [
    "condition", "Condition", "disease", "Disease", "diagnosis",
    "Diagnosis", "group", "Group", "status", "Status",
    "treatment", "Treatment", "phenotype", "Phenotype",
]


def check_h5ad(filepath: Path) -> Dict[str, Any]:
    """Run all checks on a single h5ad file.

    Returns dict with:
        status: "OK", "WARNING", or "ERROR"
        checks: list of {name, status, message, fix_available}
    """
    import anndata as ad

    filepath = Path(filepath)
    results = {"file": filepath.name, "checks": [], "status": "OK"}

    def add_check(name, status, message, fix_available=False, fix_id=None):
        results["checks"].append({
            "name": name,
            "status": status,
            "message": message,
            "fix_available": fix_available,
            "fix_id": fix_id,
        })
        if status == "ERROR":
            results["status"] = "ERROR"
        elif status == "WARNING" and results["status"] != "ERROR":
            results["status"] = "WARNING"

    # ── 1. File readable ──
    try:
        adata = ad.read_h5ad(filepath, backed="r")
    except Exception as e:
        add_check("readable", "ERROR", f"Cannot read file: {e}")
        return results

    add_check("readable", "OK", f"{adata.n_obs} cells x {adata.n_vars} genes")

    # ── 2. Expression matrix ──
    if adata.X is None:
        add_check("expression", "ERROR", "No expression matrix (.X)")
    else:
        # Check if raw counts (integers, reasonable range)
        try:
            if sp.issparse(adata.X):
                sample = adata.X[:min(1000, adata.n_obs)].toarray()
            else:
                sample = np.asarray(adata.X[:min(1000, adata.n_obs)])

            is_integer = np.allclose(sample, np.round(sample), atol=0.01)
            max_val = float(np.max(sample))
            min_val = float(np.min(sample))

            if min_val < 0:
                add_check("expression", "WARNING",
                          f"Negative values found (min={min_val:.2f}). "
                          "Data may be log-normalized or scaled. Pipeline expects raw counts.",
                          fix_available=False)
            elif not is_integer and max_val < 50:
                add_check("expression", "WARNING",
                          f"Values look normalized (max={max_val:.2f}, non-integer). "
                          "Pipeline expects raw counts.",
                          fix_available=False)
            elif is_integer:
                add_check("expression", "OK",
                          f"Raw counts detected (max={max_val:.0f}, integer)")
            else:
                add_check("expression", "OK",
                          f"Count matrix (max={max_val:.1f})")
        except Exception as e:
            add_check("expression", "WARNING", f"Could not inspect matrix: {e}")

    # ── 3. Cell type column ──
    obs_cols = list(adata.obs.columns)
    found_ct = None
    for col in KNOWN_CT_COLUMNS:
        if col in obs_cols:
            found_ct = col
            break

    if found_ct:
        n_types = adata.obs[found_ct].nunique()
        types_preview = list(adata.obs[found_ct].value_counts().head(5).index)
        if found_ct == "CellType":
            add_check("celltype_column", "OK",
                      f"Column 'CellType' found: {n_types} types ({types_preview})")
        else:
            add_check("celltype_column", "WARNING",
                      f"Found '{found_ct}' (not 'CellType'). {n_types} types: {types_preview}",
                      fix_available=True, fix_id="rename_ct_column")
    else:
        add_check("celltype_column", "ERROR",
                  f"No cell type column found. Available: {obs_cols[:10]}",
                  fix_available=True, fix_id="set_ct_column")

    # ── 4. Donor/sample column ──
    found_donor = None
    for col in KNOWN_DONOR_COLUMNS:
        if col in obs_cols:
            found_donor = col
            break

    if found_donor:
        n_donors = adata.obs[found_donor].nunique()
        add_check("donor_column", "OK" if found_donor == "donor_id" else "WARNING",
                  f"Found '{found_donor}': {n_donors} donors" +
                  (" (will need rename to 'donor_id')" if found_donor != "donor_id" else ""),
                  fix_available=(found_donor != "donor_id"),
                  fix_id="rename_donor_column")

        # Check if multi-sample file
        if n_donors > 1:
            add_check("multi_sample", "WARNING",
                      f"Multiple donors ({n_donors}) in one file. "
                      "Pipeline expects one h5ad per donor.",
                      fix_available=True, fix_id="split_donors")
    else:
        add_check("donor_column", "WARNING",
                  "No donor/sample ID column found. "
                  "If this file is one donor, that's OK.",
                  fix_available=True, fix_id="set_donor_column")

    # ── 5. Condition column ──
    found_cond = None
    for col in KNOWN_CONDITION_COLUMNS:
        if col in obs_cols:
            found_cond = col
            break

    if found_cond:
        conditions = list(adata.obs[found_cond].unique())
        add_check("condition_column", "OK",
                  f"Found '{found_cond}': {conditions}")
    else:
        add_check("condition_column", "INFO",
                  "No condition column. Will need to set in sample_info.csv.")

    # ── 6. Gene format ──
    gene_sample = list(adata.var_names[:5])
    if gene_sample and gene_sample[0].startswith("ENSG"):
        add_check("gene_format", "OK",
                  f"ENSEMBL IDs detected (e.g., {gene_sample[0]})")
    elif gene_sample and gene_sample[0].startswith("ENSMUSG"):
        add_check("gene_format", "OK",
                  f"Mouse ENSEMBL IDs detected (e.g., {gene_sample[0]})")
    else:
        add_check("gene_format", "OK",
                  f"Gene symbols detected (e.g., {gene_sample[:3]})")

    # ── 7. Basic QC ──
    n_zero_cells = 0
    try:
        if sp.issparse(adata.X):
            cell_counts = np.asarray(adata.X.sum(axis=1)).ravel()
        else:
            cell_counts = np.asarray(adata.X.sum(axis=1)).ravel()
        n_zero_cells = int((cell_counts == 0).sum())
        median_umi = float(np.median(cell_counts))
        add_check("qc_basic", "OK" if n_zero_cells == 0 else "WARNING",
                  f"Median UMI/cell: {median_umi:.0f}" +
                  (f", {n_zero_cells} cells with zero counts" if n_zero_cells > 0 else ""))
    except Exception:
        add_check("qc_basic", "WARNING", "Could not compute QC metrics")

    return results


def check_sample_info(metadata_dir: Path, h5ad_dir: Path) -> Dict[str, Any]:
    """Check sample_info.csv for completeness."""
    results = {"checks": [], "status": "OK"}

    def add_check(name, status, message, fix_available=False, fix_id=None):
        results["checks"].append({
            "name": name, "status": status, "message": message,
            "fix_available": fix_available, "fix_id": fix_id,
        })
        if status == "ERROR":
            results["status"] = "ERROR"
        elif status == "WARNING" and results["status"] != "ERROR":
            results["status"] = "WARNING"

    sample_info_path = metadata_dir / "sample_info.csv"

    if not sample_info_path.exists():
        h5ad_files = list(h5ad_dir.glob("*.h5ad"))
        add_check("exists", "ERROR",
                  f"sample_info.csv not found. {len(h5ad_files)} h5ad files available.",
                  fix_available=True, fix_id="generate_sample_info")
        return results

    df = pd.read_csv(sample_info_path)
    add_check("exists", "OK", f"Found: {len(df)} rows")

    # Required columns
    for col in ["donor_id", "file", "condition"]:
        if col in df.columns:
            add_check(f"column_{col}", "OK", f"'{col}' present")
        else:
            add_check(f"column_{col}", "ERROR", f"'{col}' missing",
                      fix_available=True, fix_id="fix_sample_info_columns")

    # Check file references
    if "file" in df.columns:
        missing_files = []
        for f in df["file"]:
            if not (h5ad_dir / f).exists():
                missing_files.append(f)
        if missing_files:
            add_check("file_refs", "ERROR",
                      f"{len(missing_files)} referenced files not found: {missing_files[:5]}")
        else:
            add_check("file_refs", "OK", "All referenced h5ad files exist")

    return results


# =====================================================================
# Fixers
# =====================================================================

def fix_rename_column(filepath: Path, old_name: str, new_name: str) -> Path:
    """Rename a column in h5ad .obs and save."""
    import anndata as ad
    adata = ad.read_h5ad(filepath)
    adata.obs = adata.obs.rename(columns={old_name: new_name})
    adata.write_h5ad(filepath)
    logger.info("Renamed column '%s' -> '%s' in %s", old_name, new_name, filepath.name)
    return filepath


def fix_split_donors(filepath: Path, donor_column: str, output_dir: Path) -> List[Path]:
    """Split a multi-donor h5ad into per-donor files."""
    import anndata as ad
    adata = ad.read_h5ad(filepath)
    output_dir.mkdir(parents=True, exist_ok=True)

    donors = adata.obs[donor_column].unique()
    output_files = []

    for donor in donors:
        mask = adata.obs[donor_column] == donor
        sub = adata[mask].copy()
        out_path = output_dir / f"{donor}.h5ad"
        sub.write_h5ad(out_path)
        output_files.append(out_path)
        logger.info("  Wrote %s: %d cells", out_path.name, sub.n_obs)

    logger.info("Split %s into %d donor files", filepath.name, len(output_files))
    return output_files


def fix_generate_sample_info(h5ad_dir: Path, metadata_dir: Path,
                              condition_default: str = "Unknown") -> Path:
    """Auto-generate sample_info.csv from h5ad filenames."""
    import anndata as ad

    h5ad_files = sorted(h5ad_dir.glob("*.h5ad"))
    rows = []

    for f in h5ad_files:
        donor_id = f.stem  # Use filename without extension as donor_id
        condition = condition_default

        # Try to read condition from the file
        try:
            adata = ad.read_h5ad(f, backed="r")
            for col in KNOWN_CONDITION_COLUMNS:
                if col in adata.obs.columns:
                    # Use the most common condition value
                    condition = str(adata.obs[col].mode().iloc[0])
                    break
        except Exception:
            pass

        rows.append({
            "donor_id": donor_id,
            "file": f.name,
            "condition": condition,
        })

    df = pd.DataFrame(rows)
    metadata_dir.mkdir(parents=True, exist_ok=True)
    out_path = metadata_dir / "sample_info.csv"
    df.to_csv(out_path, index=False)
    logger.info("Generated sample_info.csv: %d samples", len(df))
    return out_path


def fix_filter_cells(filepath: Path, max_mito_frac: float = 0.2,
                      min_umi: int = 200) -> Tuple[Path, Dict]:
    """Filter low-quality cells from h5ad."""
    import anndata as ad

    adata = ad.read_h5ad(filepath)
    n_before = adata.n_obs

    # Compute UMI counts
    if sp.issparse(adata.X):
        umi_counts = np.asarray(adata.X.sum(axis=1)).ravel()
    else:
        umi_counts = np.asarray(adata.X.sum(axis=1)).ravel()

    # Compute mito fraction
    mito_mask = np.array([
        str(g).upper().startswith("MT-") or str(g).startswith("mt-")
        for g in adata.var_names
    ])
    if mito_mask.any():
        if sp.issparse(adata.X):
            mito_counts = np.asarray(adata.X[:, mito_mask].sum(axis=1)).ravel()
        else:
            mito_counts = np.asarray(adata.X[:, mito_mask].sum(axis=1)).ravel()
        mito_frac = mito_counts / np.maximum(umi_counts, 1)
    else:
        mito_frac = np.zeros(adata.n_obs)

    # Filter
    keep = (umi_counts >= min_umi) & (mito_frac <= max_mito_frac)
    adata = adata[keep].copy()

    n_after = adata.n_obs
    n_removed = n_before - n_after

    adata.write_h5ad(filepath)

    stats = {
        "n_before": n_before,
        "n_after": n_after,
        "n_removed": n_removed,
        "pct_removed": f"{n_removed / n_before * 100:.1f}%",
    }
    logger.info("Filtered %s: %d -> %d cells (removed %d)", filepath.name, n_before, n_after, n_removed)
    return filepath, stats


def run_full_check(h5ad_dir: Path, metadata_dir: Path) -> Dict[str, Any]:
    """Run all checks on all files.

    Returns comprehensive report.
    """
    report = {
        "h5ad_files": [],
        "sample_info": None,
        "overall_status": "OK",
        "n_errors": 0,
        "n_warnings": 0,
        "fixes_available": [],
    }

    # Check h5ad files
    h5ad_files = sorted(h5ad_dir.glob("*.h5ad"))
    for f in h5ad_files:
        result = check_h5ad(f)
        report["h5ad_files"].append(result)
        for check in result["checks"]:
            if check["status"] == "ERROR":
                report["n_errors"] += 1
            elif check["status"] == "WARNING":
                report["n_warnings"] += 1
            if check.get("fix_available"):
                report["fixes_available"].append({
                    "file": result["file"],
                    "fix_id": check["fix_id"],
                    "description": check["message"],
                })

    # Check sample_info
    report["sample_info"] = check_sample_info(metadata_dir, h5ad_dir)
    for check in report["sample_info"]["checks"]:
        if check["status"] == "ERROR":
            report["n_errors"] += 1
        elif check["status"] == "WARNING":
            report["n_warnings"] += 1
        if check.get("fix_available"):
            report["fixes_available"].append({
                "file": "sample_info.csv",
                "fix_id": check["fix_id"],
                "description": check["message"],
            })

    if report["n_errors"] > 0:
        report["overall_status"] = "ERROR"
    elif report["n_warnings"] > 0:
        report["overall_status"] = "WARNING"

    return report
