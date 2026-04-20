"""
3φ Residual Framework (Section 4.1-4.2, 6.15, 8.9 of the paper).

Computes three Hodge potential (φ) variants per cell type:
  - φ_static    : from pathologically normal (PN) donors only
                  → captures healthy-state co-expression network topology
  - φ_disease   : from disease donors only, pseudotime-based multi-window
                  → captures disease-state co-expression structure
  - φ_condition : from pooled disease + PN donors, direct ALS vs PN Δ
                  → captures the binary case-control co-expression difference

Then regresses φ_disease on φ_static via a cubic polynomial fit and computes
per-gene residuals. Standardised residual z-scores identify:
  - "rewiring" genes  (z > +2): pathological co-expression beyond topology
  - "collapse" genes  (z < -2): disruption exceeding topology prediction

This framework is the basis for:
  - Table 1 R² manifold preservation gradient (Section 4.1)
  - Rewiring/collapse gene counts (Section 4.2)
  - NEMF universal disease-specific downshift (Section 8.9.1)
  - Pathway-level matched null z-score tests (Section 4.2, Figure 5)

Usage:
    from scripts.three_phi_residual import run_three_phi, compute_residual_z

    results = run_three_phi(
        adata=per_donor_adata_dict,
        donor_condition=donor_condition_dict,
        base_gene_indices=base_idx,
        pseudotime_df=pt_df,
        cell_types=["Oligo", "OPC", "Astro", ...],
    )

    # results[ct] = {
    #     "phi_static": array(n_genes),
    #     "phi_disease": array(n_genes),
    #     "phi_condition": array(n_genes),
    #     "gf_static", "gf_disease", "gf_condition": float,
    # }

    residual_z = compute_residual_z(
        phi_disease=results["Oligo"]["phi_disease"],
        phi_static=results["Oligo"]["phi_static"],
        degree=3,  # cubic polynomial
    )

Migrated from: sals_analysis_frozen_20260211/scripts/run_condition_phi.py
              sals_analysis_frozen_20260211/scripts/run_3phi_verification.py
"""

from __future__ import annotations

import gc
import logging
from collections import Counter
from pathlib import Path
from typing import Dict, List, Optional, Sequence

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from scripts.gene_hodge import (
    precompute_donor_log_corr,
    run_single_transition_hodge,
    hodge_gradient_kn,
)

logger = logging.getLogger(__name__)


def _edge_weight_flow(delta: np.ndarray) -> np.ndarray:
    """Construct antisymmetric edge-weight flow f(i,j) = |Δ_ij| · sign(d_i - d_j).

    This is the axiom-consistent edge flow from Methods 6.4 / Supplementary Note 3.
    """
    n = delta.shape[0]
    d = np.linalg.norm(delta, axis=1)
    n_edges = n * (n - 1) // 2
    flow = np.empty(n_edges, dtype=np.float64)
    idx = 0
    for i in range(n):
        n_pairs = n - i - 1
        if n_pairs == 0:
            break
        js = np.arange(i + 1, n)
        flow[idx : idx + n_pairs] = np.abs(delta[i, js]) * np.sign(d[i] - d[js])
        idx += n_pairs

    std_f = flow.std()
    if std_f > 0:
        flow /= std_f
    return flow, d


def compute_phi_from_log_corr(
    log_corr: np.ndarray,
) -> Dict[str, np.ndarray | float]:
    """Compute Hodge φ directly from a log-correlation matrix (no pseudotime delta).

    Used for φ_static: treats the PN-only mean log-correlation matrix as a
    static "distance" object whose row norms define an engagement scalar d_i.
    """
    flow, d = _edge_weight_flow(log_corr)
    n_genes = log_corr.shape[0]
    result = hodge_gradient_kn(flow, n_genes)
    result["d"] = d
    return result


def compute_phi_condition(
    mean_log_als: np.ndarray,
    mean_log_pn: np.ndarray,
) -> Dict[str, np.ndarray | float]:
    """Compute Hodge φ_condition from direct ALS vs PN delta (no pseudotime windows).

    Δ_condition = mean(log_corr across ALS donors) - mean(log_corr across PN donors)
    """
    delta = mean_log_als - mean_log_pn
    flow, d = _edge_weight_flow(delta)
    n_genes = delta.shape[0]
    result = hodge_gradient_kn(flow, n_genes)
    result["d"] = d
    result["delta"] = delta
    return result


def run_three_phi_for_celltype(
    residuals: np.ndarray,
    donor_ids: np.ndarray,
    donor_condition: Dict[str, str],
    base_gene_indices: np.ndarray,
    pseudotime_df: pd.DataFrame,
    flow_mode: str = "edge_weight",
    min_donors_per_window: int = 2,
    min_donors_per_condition: int = 3,
    ct_label: str = "",
) -> Optional[Dict]:
    """Run the three-φ framework for a single cell type.

    Parameters
    ----------
    residuals : (n_cells, n_genes) ndarray of regression residuals
    donor_ids : (n_cells,) donor labels
    donor_condition : dict mapping donor_id -> "ALS" or "PN" (or other)
    base_gene_indices : indices into residuals columns defining the base gene set
    pseudotime_df : DataFrame with columns [donor_id, pt, window]
    control_label : label in donor_condition for pathological normal donors
    disease_label : label in donor_condition for disease donors

    Returns
    -------
    dict with keys:
        phi_static, phi_disease, phi_condition : arrays of shape (n_base_genes,)
        gf_static, gf_disease, gf_condition : gradient fraction floats
        d_static, d_condition : engagement scalars (row L2 norms)
        n_als, n_pn : donor counts
        w_star : pseudotime window index used for phi_disease (or None)
    """
    n_genes = len(base_gene_indices)

    # Step 1: per-donor log(corr) on base set
    dlc, dw = precompute_donor_log_corr(residuals, donor_ids, base_gene_indices, pseudotime_df)
    dids = sorted(dlc.keys())

    als_donors = [d for d in dids if donor_condition.get(d) == "ALS"]
    pn_donors = [d for d in dids if donor_condition.get(d) == "PN"]
    logger.info(f"[{ct_label}] donors: {len(als_donors)} ALS, {len(pn_donors)} PN")

    if len(als_donors) < min_donors_per_condition or len(pn_donors) < min_donors_per_condition:
        logger.warning(f"[{ct_label}] Insufficient donors per condition, skipping")
        return None

    # Step 2: φ_condition = Hodge on Δ(ALS − PN)
    mean_log_als = np.mean([dlc[d] for d in als_donors], axis=0)
    mean_log_pn = np.mean([dlc[d] for d in pn_donors], axis=0)
    cond = compute_phi_condition(mean_log_als, mean_log_pn)
    logger.info(f"[{ct_label}] GF_condition = {cond['gradient_fraction']:.4f}")

    # Step 3: φ_static = Hodge on PN-only mean log(corr), no delta
    static = compute_phi_from_log_corr(mean_log_pn)
    logger.info(f"[{ct_label}] GF_static = {static['gradient_fraction']:.4f}")

    # Step 4: φ_disease = pseudotime-based single-transition Hodge (disease cohort only)
    # Selects the window transition w_star with maximum d_corr (most informative Δ).
    phi_disease = None
    gf_disease = None
    w_star = None

    disease_dids = als_donors
    wc = Counter(dw[d] for d in disease_dids)
    valid_transitions = [w for w in sorted(wc) if wc.get(w, 0) >= min_donors_per_window and wc.get(w + 1, 0) >= min_donors_per_window]

    if valid_transitions:
        # Default: use the first valid transition
        w_star = valid_transitions[0]
        disease_dlc = {d: dlc[d] for d in disease_dids}
        disease_dw = {d: dw[d] for d in disease_dids}
        baseline = run_single_transition_hodge(
            disease_dlc, disease_dw, disease_dids, w_star, n_genes, flow_mode=flow_mode
        )
        if baseline is not None:
            phi_disease = baseline["phi"]
            gf_disease = baseline["gradient_fraction"]
            logger.info(f"[{ct_label}] GF_disease = {gf_disease:.4f} (window {w_star}→{w_star+1})")

    # Step 5: Concordance diagnostics
    comparisons = {
        "condition_vs_static": spearmanr(cond["phi"], static["phi"]).correlation,
    }
    if phi_disease is not None:
        comparisons["disease_vs_condition"] = spearmanr(phi_disease, cond["phi"]).correlation
        comparisons["disease_vs_static"] = spearmanr(phi_disease, static["phi"]).correlation

    del dlc
    gc.collect()

    return {
        "phi_static": static["phi"],
        "phi_disease": phi_disease,
        "phi_condition": cond["phi"],
        "gf_static": static["gradient_fraction"],
        "gf_disease": gf_disease,
        "gf_condition": cond["gradient_fraction"],
        "d_static": static["d"],
        "d_condition": cond["d"],
        "delta_condition": cond["delta"],
        "n_als": len(als_donors),
        "n_pn": len(pn_donors),
        "n_genes": n_genes,
        "w_star": w_star,
        "comparisons": comparisons,
    }


def compute_residual_z(
    phi_disease: np.ndarray,
    phi_static: np.ndarray,
    degree: int = 3,
) -> Dict[str, np.ndarray | float]:
    """Compute poly3 residual z-scores (Section 4.2 / 6.15).

    Regresses φ_disease on φ_static via a polynomial of given degree (default cubic),
    then standardises residuals per cell type.

    Returns
    -------
    dict:
        z         : standardised residuals, shape (n_genes,)
        residuals : raw residuals before standardisation
        fit_r2    : R² of the polynomial fit (manifold preservation metric)
        coeffs    : polynomial coefficients (high order first, numpy convention)
    """
    if phi_disease is None or phi_static is None:
        return {"z": None, "residuals": None, "fit_r2": None, "coeffs": None}

    coeffs = np.polyfit(phi_static, phi_disease, deg=degree)
    predicted = np.polyval(coeffs, phi_static)
    residuals = phi_disease - predicted

    ss_res = np.sum(residuals ** 2)
    ss_tot = np.sum((phi_disease - phi_disease.mean()) ** 2)
    fit_r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0

    sigma = residuals.std()
    if sigma > 0:
        z = residuals / sigma
    else:
        z = np.zeros_like(residuals)

    return {
        "z": z,
        "residuals": residuals,
        "fit_r2": float(fit_r2),
        "coeffs": coeffs,
    }


def matched_null_z_score(
    residuals: np.ndarray,
    donor_ids: np.ndarray,
    donor_condition: Dict[str, str],
    base_gene_indices: np.ndarray,
    pseudotime_df: pd.DataFrame,
    gene_indices_of_interest: Sequence[int],
    n_permutations: int = 1000,
    seed: int = 42,
) -> np.ndarray:
    """Matched null z-score for selected genes (Section 6.15, 8.9.1).

    For each permutation:
      1. Shuffle donor ALS/PN labels (preserving cell type assignment)
      2. Re-run three_phi + poly3 regression
      3. Record residual z under the permutation

    The matched null z-score is (z_observed − z_null_mean) / z_null_std,
    corrected for gene-specific variance structure.

    Returns
    -------
    matched_null_z : array shape (n_genes_of_interest,)

    Notes
    -----
    - Used for the genome-wide NEMF screen (single-gene z < -2 in ≥7/10 CTs).
    - Used for pathway-level upstream-ness tests (Translation 9-10/10 CTs,
      chaperone 4/10, ATP synthase 3/10).
    - Donor labels are permuted while preserving the cell-type assignment
      (matching the experimental unit).
    """
    rng = np.random.RandomState(seed)
    all_dids = sorted(donor_condition.keys())

    null_z_matrix = np.full((n_permutations, len(gene_indices_of_interest)), np.nan)

    # Compute observed z once (no permutation)
    obs = run_three_phi_for_celltype(
        residuals, donor_ids, donor_condition,
        base_gene_indices, pseudotime_df,
    )
    if obs is None or obs["phi_disease"] is None:
        return np.full(len(gene_indices_of_interest), np.nan)

    obs_resid = compute_residual_z(obs["phi_disease"], obs["phi_static"])
    if obs_resid["z"] is None:
        return np.full(len(gene_indices_of_interest), np.nan)
    z_observed = obs_resid["z"][gene_indices_of_interest]

    # Permutation loop
    for p in range(n_permutations):
        perm_labels = rng.permutation([donor_condition[d] for d in all_dids])
        perm_cond = dict(zip(all_dids, perm_labels))

        try:
            perm_obs = run_three_phi_for_celltype(
                residuals, donor_ids, perm_cond,
                base_gene_indices, pseudotime_df,
            )
            if perm_obs is None or perm_obs["phi_disease"] is None:
                continue
            perm_resid = compute_residual_z(perm_obs["phi_disease"], perm_obs["phi_static"])
            if perm_resid["z"] is not None:
                null_z_matrix[p, :] = perm_resid["z"][gene_indices_of_interest]
        except Exception as e:
            logger.warning(f"Permutation {p} failed: {e}")
            continue

    null_mean = np.nanmean(null_z_matrix, axis=0)
    null_std = np.nanstd(null_z_matrix, axis=0, ddof=1)
    null_std = np.where(null_std > 0, null_std, 1.0)

    matched_z = (z_observed - null_mean) / null_std
    return matched_z


def classify_rewiring_collapse(
    z: np.ndarray,
    threshold: float = 2.0,
) -> Dict[str, np.ndarray]:
    """Classify genes as rewiring (z > +T) or collapse (z < -T).

    Section 4.2: T = 2 for the main classification; T = 1 is used for
    Metascape enrichment to preserve statistical power.
    """
    return {
        "rewiring_idx": np.where(z > threshold)[0],
        "collapse_idx": np.where(z < -threshold)[0],
        "n_rewiring": int((z > threshold).sum()),
        "n_collapse": int((z < -threshold).sum()),
        "ratio": float((z > threshold).sum()) / max(1, int((z < -threshold).sum())),
    }


# --------------------------------------------------------------------------
# CLI entry point (runs 3φ framework for all configured cell types)
# --------------------------------------------------------------------------
def main():
    """Execute the 3φ residual framework for all configured cell types.

    Reads from project_config.yaml and writes results to
    results/three_phi_residual/{ct}_3phi.csv, plus a summary JSON.
    """
    import json
    from scripts import config
    from scripts.data_loader import load_all_samples, get_cells_for_celltype, get_sample_info
    from scripts.residuals import compute_residuals

    logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")

    out_dir = Path("results/three_phi_residual")
    out_dir.mkdir(parents=True, exist_ok=True)

    sample_info = get_sample_info()
    donor_condition = dict(zip(sample_info["donor_id"].astype(str), sample_info["condition"]))

    pt_df_path = Path("results/pseudotime/pt_b_with_windows.csv")
    if not pt_df_path.exists():
        raise FileNotFoundError(
            f"Pseudotime file not found: {pt_df_path}. "
            f"Run `python run_pipeline.py --step pseudotime` first."
        )
    pt_df = pd.read_csv(pt_df_path)

    logger.info("Loading all samples...")
    all_adata = load_all_samples()

    summary = []
    for ct in config.CELL_TYPES:
        logger.info(f"\n=== {ct} ===")
        try:
            adata = get_cells_for_celltype(all_adata, ct, min_cells_per_donor=config.MIN_CELLS_PER_DONOR)
            donor_ids = adata.obs["donor_id"].values.astype(str)
            residuals, _, res_genes = compute_residuals(adata, hvg_n=None)
            del adata
            gc.collect()

            n_base = getattr(config, "N_BASE_GENES", 3500)
            gene_var = np.var(residuals, axis=0)
            base_idx = np.sort(np.argsort(gene_var)[::-1][:n_base])

            result = run_three_phi_for_celltype(
                residuals, donor_ids, donor_condition,
                base_idx, pt_df, ct_label=ct,
            )
            if result is None:
                continue

            # Poly3 residual z
            resid = compute_residual_z(result["phi_disease"], result["phi_static"])

            # Per-gene output
            gene_symbols = [res_genes[i] for i in base_idx]
            rows = []
            for i in range(result["n_genes"]):
                rows.append({
                    "gene": gene_symbols[i],
                    "phi_static": float(result["phi_static"][i]),
                    "phi_disease": float(result["phi_disease"][i]) if result["phi_disease"] is not None else np.nan,
                    "phi_condition": float(result["phi_condition"][i]),
                    "residual_z": float(resid["z"][i]) if resid["z"] is not None else np.nan,
                })
            df = pd.DataFrame(rows)
            df.to_csv(out_dir / f"{ct.lower()}_3phi.csv", index=False)

            classification = classify_rewiring_collapse(resid["z"]) if resid["z"] is not None else None
            summary.append({
                "ct": ct,
                "n_genes": result["n_genes"],
                "n_als": result["n_als"],
                "n_pn": result["n_pn"],
                "gf_static": result["gf_static"],
                "gf_disease": result["gf_disease"],
                "gf_condition": result["gf_condition"],
                "manifold_R2": resid["fit_r2"],
                "n_rewiring_z2": classification["n_rewiring"] if classification else None,
                "n_collapse_z2": classification["n_collapse"] if classification else None,
                "ratio_rew_col": classification["ratio"] if classification else None,
                "comparisons": result["comparisons"],
            })

            del residuals
            gc.collect()

        except Exception as e:
            logger.error(f"{ct} FAILED: {e}")
            import traceback
            traceback.print_exc()
            continue

    del all_adata
    gc.collect()

    with open(out_dir / "three_phi_summary.json", "w") as f:
        json.dump(summary, f, indent=2, default=str)

    logger.info(f"\nSaved 3φ results to {out_dir}")
    logger.info("\nManifold preservation R² (disease ~ static):")
    for s in sorted(summary, key=lambda x: -(x["manifold_R2"] or 0)):
        r2 = s["manifold_R2"]
        if r2 is not None:
            logger.info(f"  {s['ct']:>8s} : {r2:.3f}")


if __name__ == "__main__":
    main()
