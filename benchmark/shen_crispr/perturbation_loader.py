"""
P4 Validation — Perturbation Data Loader
=========================================
Loads Spatial Perturb-Seq scRNA-seq data (10X filtered matrices),
classifies cells by perturbation status, applies QC, and normalizes.

Reuses IDS-SALS core modules (spd.py) for downstream analysis.
"""
import gzip
import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import scipy.io
import scipy.sparse as sp

logger = logging.getLogger(__name__)

# ── Perturbation gene metadata ────────────────────────────────────
PERTURBATION_GENES = {
    # gRNA_name_in_matrix: (short_name, disease_category)
    "mSafe_gRNA": ("mSafe", "Control"),
    "Trem2_gRNA": ("Trem2", "AD"),
    "Rraga_gRNA": ("Rraga", "Control"),
    "Fasn_gRNA": ("Fasn", "Control"),
    "Clu_gRNA": ("Clu", "AD"),
    "Dpp6_gRNA": ("Dpp6", "ALS"),
    "Tbk1_gRNA": ("Tbk1", "ALS"),
    "Flcn_gRNA": ("Flcn", "Control"),
    "Gfap_gRNA": ("Gfap", "Control"),
    "C9orf72_gRNA": ("C9orf72", "ALS"),
    "Cfap410_gRNA": ("Cfap410", "ALS"),
    "Stk39_gRNA": ("Stk39", "PD"),
    "Lrrk2_gRNA": ("Lrrk2", "PD"),
    "Ndufaf_gRNA": ("Ndufaf2", "AD"),
    "Sh3gl2_gRNA": ("Sh3gl2", "PD"),
    "Srf_gRNA": ("Srf", "Control"),
    "Rbfox_gRNA": ("Rbfox1", "AD"),
    "Olig2_gRNA": ("Olig2", "Control"),
}

ALS_GENES = ["C9orf72", "Tbk1", "Dpp6", "Cfap410"]
PD_GENES = ["Lrrk2", "Stk39", "Sh3gl2"]
AD_GENES = ["Trem2", "Clu", "Ndufaf2", "Rbfox1"]
CONTROL_GENES = ["mSafe", "Fasn", "Flcn", "Gfap", "Olig2", "Rraga", "Srf"]


class PerturbationDataset:
    """Merged scRNA-seq data from Spatial Perturb-Seq experiment."""

    def __init__(self, matrix_dirs: List[str], min_genes: int = 200,
                 max_genes: int = 7000, max_pct_mito: float = 0.15):
        """
        Load and merge 10X filtered matrices from multiple libraries.

        Parameters
        ----------
        matrix_dirs : list of str
            Paths to directories containing matrix.mtx.gz, features.tsv.gz, barcodes.tsv.gz
        min_genes, max_genes : int
            QC thresholds for genes detected per cell
        max_pct_mito : float
            Maximum mitochondrial fraction
        """
        self.library_ids = []
        self.gene_names = None
        self.gene_ids = None
        self.n_libraries = len(matrix_dirs)

        # ── Load and merge ──
        all_matrices = []
        all_cell_barcodes = []
        all_library_labels = []

        for i, mdir in enumerate(matrix_dirs):
            mdir = Path(mdir)
            lib_id = mdir.name
            self.library_ids.append(lib_id)

            logger.info(f"Loading library {lib_id}...")

            # Features
            with gzip.open(mdir / "features.tsv.gz", "rt") as f:
                features = [line.strip().split("\t") for line in f]
            gnames = [ft[1] for ft in features]
            gids = [ft[0] for ft in features]

            if self.gene_names is None:
                self.gene_names = np.array(gnames)
                self.gene_ids = np.array(gids)
            else:
                assert list(self.gene_names) == gnames, f"Feature mismatch in {lib_id}"

            # Barcodes
            with gzip.open(mdir / "barcodes.tsv.gz", "rt") as f:
                barcodes = [line.strip() for line in f]

            # Matrix (genes × cells) → transpose to (cells × genes)
            mat = scipy.io.mmread(mdir / "matrix.mtx.gz").tocsc()
            mat_t = mat.T.tocsr()  # (cells × genes)

            all_matrices.append(mat_t)
            all_cell_barcodes.extend([f"{lib_id}_{bc}" for bc in barcodes])
            all_library_labels.extend([lib_id] * len(barcodes))

            logger.info(f"  {lib_id}: {mat_t.shape[0]} cells, {mat_t.shape[1]} genes")

        # Stack
        self.X = sp.vstack(all_matrices, format="csr")  # (total_cells × genes)
        self.cell_barcodes = np.array(all_cell_barcodes)
        self.library_labels = np.array(all_library_labels)
        n_total = self.X.shape[0]
        logger.info(f"Merged: {n_total} cells × {self.X.shape[1]} genes")

        # ── Identify gRNA features ──
        self.grna_indices = {}
        for grna_name in PERTURBATION_GENES:
            idx = np.where(self.gene_names == grna_name)[0]
            if len(idx) == 1:
                self.grna_indices[grna_name] = idx[0]

        # Exclude Myrf (all zeros)
        self.grna_indices = {k: v for k, v in self.grna_indices.items()
                            if self.X[:, v].nnz > 0}
        logger.info(f"Active gRNA features: {len(self.grna_indices)}")

        # ── Classify cells by perturbation ──
        self._classify_cells()

        # ── QC filtering ──
        self._qc_filter(min_genes, max_genes, max_pct_mito)

        # ── Identify endogenous gene indices (exclude gRNA features) ──
        all_grna_idx = set(self.grna_indices.values())
        myrf_idx = np.where(self.gene_names == "Myrf_gRNA")[0]
        if len(myrf_idx) > 0:
            all_grna_idx.add(myrf_idx[0])
        self.endogenous_mask = np.array([i not in all_grna_idx
                                          for i in range(len(self.gene_names))])
        self.n_endogenous_genes = self.endogenous_mask.sum()
        logger.info(f"Endogenous genes: {self.n_endogenous_genes}")

    def _classify_cells(self):
        """Classify each cell as perturbation-positive or negative."""
        n_cells = self.X.shape[0]
        self.perturbation_label = np.full(n_cells, "Negative", dtype=object)
        self.n_barcodes_per_cell = np.zeros(n_cells, dtype=int)

        for grna_name, gidx in self.grna_indices.items():
            col = np.asarray(self.X[:, gidx].todense()).ravel()
            positive = col > 0
            self.n_barcodes_per_cell[positive] += 1

        # Assign labels for single-barcode cells
        for grna_name, gidx in self.grna_indices.items():
            col = np.asarray(self.X[:, gidx].todense()).ravel()
            positive = col > 0
            single_bc = positive & (self.n_barcodes_per_cell == 1)
            short_name = PERTURBATION_GENES[grna_name][0]
            self.perturbation_label[single_bc] = short_name

        # Mark multi-barcode cells
        multi = self.n_barcodes_per_cell > 1
        self.perturbation_label[multi] = "Multiple"

        # Summary
        unique, counts = np.unique(self.perturbation_label, return_counts=True)
        logger.info("Cell classification:")
        for u, c in sorted(zip(unique, counts), key=lambda x: -x[1]):
            logger.info(f"  {u}: {c}")

    def _qc_filter(self, min_genes, max_genes, max_pct_mito):
        """Apply QC filters."""
        n_before = self.X.shape[0]

        # Genes per cell
        genes_per_cell = np.diff(self.X.tocsr().indptr)

        # Detect mito genes (mouse: mt- prefix)
        mito_mask = np.array([g.startswith("mt-") for g in self.gene_names])
        if mito_mask.any():
            mito_counts = np.asarray(self.X[:, mito_mask].sum(axis=1)).ravel()
            total_counts = np.asarray(self.X.sum(axis=1)).ravel()
            pct_mito = np.where(total_counts > 0, mito_counts / total_counts, 0)
        else:
            pct_mito = np.zeros(n_before)

        self.pct_mito = pct_mito

        # Filter
        keep = (genes_per_cell >= min_genes) & (genes_per_cell <= max_genes) & (pct_mito <= max_pct_mito)

        self.X = self.X[keep]
        self.cell_barcodes = self.cell_barcodes[keep]
        self.library_labels = self.library_labels[keep]
        self.perturbation_label = self.perturbation_label[keep]
        self.n_barcodes_per_cell = self.n_barcodes_per_cell[keep]
        self.pct_mito = self.pct_mito[keep]

        n_after = self.X.shape[0]
        logger.info(f"QC filter: {n_before} → {n_after} cells ({n_before - n_after} removed)")

    def normalize_log1p_cpm(self):
        """Normalize to log1p(CPM) on endogenous genes only. Returns dense matrix."""
        X_endo = self.X[:, self.endogenous_mask].astype(np.float64)
        if sp.issparse(X_endo):
            X_endo = X_endo.toarray()
        else:
            X_endo = np.asarray(X_endo)

        # CPM normalization
        total = X_endo.sum(axis=1, keepdims=True)
        total = np.where(total > 0, total, 1.0)
        X_cpm = X_endo / total * 1e6

        # log1p
        X_norm = np.log1p(X_cpm).astype(np.float32)
        self.X_norm = X_norm
        self.gene_names_endogenous = self.gene_names[self.endogenous_mask]
        logger.info(f"Normalized: {X_norm.shape[0]} cells × {X_norm.shape[1]} genes")
        return X_norm

    def get_perturbation_cells(self, gene_short_name: str) -> np.ndarray:
        """Return boolean mask for cells with a specific perturbation."""
        return self.perturbation_label == gene_short_name

    def get_control_cells(self) -> np.ndarray:
        """Return boolean mask for barcode-negative cells."""
        return self.perturbation_label == "Negative"

    def get_cell_indices(self, gene_short_name: str,
                         control_ratio: int = 10,
                         seed: int = 42) -> Tuple[np.ndarray, np.ndarray]:
        """Return (KO_indices, Control_indices) for a perturbation.

        Control cells are subsampled to control_ratio × n_KO.
        """
        ko_mask = self.get_perturbation_cells(gene_short_name)
        ctrl_mask = self.get_control_cells()

        ko_idx = np.where(ko_mask)[0]
        ctrl_all = np.where(ctrl_mask)[0]

        n_ctrl = min(len(ctrl_all), control_ratio * len(ko_idx))
        rng = np.random.default_rng(seed)
        ctrl_idx = rng.choice(ctrl_all, size=n_ctrl, replace=False)

        return ko_idx, ctrl_idx

    def summary(self) -> dict:
        """Return summary statistics."""
        unique, counts = np.unique(self.perturbation_label, return_counts=True)
        perturb_counts = dict(zip(unique, counts.astype(int)))
        return {
            "total_cells": self.X.shape[0],
            "total_genes": self.X.shape[1],
            "endogenous_genes": int(self.n_endogenous_genes),
            "n_libraries": self.n_libraries,
            "perturbation_counts": perturb_counts,
        }
