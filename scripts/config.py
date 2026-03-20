"""
scRNAseq Hodge Decomposition Pipeline — Master Configuration
=============================================================
Loads project_config.yaml and provides all parameters as module-level constants.
This file is the single source of truth for the entire pipeline.
"""
import hashlib
import os
from pathlib import Path
from dataclasses import dataclass, field
from typing import List, Tuple, Optional, Dict, Any

import yaml

# ── Load project config ──────────────────────────────────────────
PIPELINE_ROOT = Path(__file__).resolve().parent.parent
_CONFIG_PATH = PIPELINE_ROOT / "project_config.yaml"

def _load_yaml_config() -> Dict[str, Any]:
    """Load project_config.yaml. Returns empty dict if not found."""
    if _CONFIG_PATH.exists():
        with open(_CONFIG_PATH, "r", encoding="utf-8") as f:
            return yaml.safe_load(f) or {}
    return {}

_cfg = _load_yaml_config()

# ── Project Info ─────────────────────────────────────────────────
PROJECT_NAME = _cfg.get("project_name", "scrnaseq_project")

# ── Paths ────────────────────────────────────────────────────────
_data_cfg = _cfg.get("data", {})

def _resolve_path(p: Optional[str]) -> Optional[Path]:
    if p is None:
        return None
    path = Path(p)
    if not path.is_absolute():
        path = PIPELINE_ROOT / path
    return path

H5AD_DIR = _resolve_path(_data_cfg.get("h5ad_dir", "data/h5ad"))
SAMPLE_INFO_PATH = _resolve_path(_data_cfg.get("sample_info", "data/metadata/sample_info.csv"))
GENE_ANNOTATION_PATH = _resolve_path(_data_cfg.get("gene_annotation"))
FUNCTIONAL_MODULES_PATH = _resolve_path(_data_cfg.get("functional_modules"))
RESULTS_DIR = PIPELINE_ROOT / "results"

# ── Hardware ─────────────────────────────────────────────────────
_hw_cfg = _cfg.get("hardware", {})
N_CPU_CORES = _hw_cfg.get("n_cpu_cores", 8)
GPU_DEVICE = _hw_cfg.get("gpu_device", "cuda:0")
USE_GPU = _hw_cfg.get("use_gpu", True)

# ── Cell Types ───────────────────────────────────────────────────
_ct_cfg = _cfg.get("cell_types", {})
CELLTYPE_COLUMN = _ct_cfg.get("column_name", "CellType")
CELLTYPE_LIST = _ct_cfg.get("include", [])
PT_CELLTYPE_SET = _ct_cfg.get("pseudotime_set", CELLTYPE_LIST[:4] if CELLTYPE_LIST else [])

# ── Conditions ───────────────────────────────────────────────────
_cond_cfg = _cfg.get("conditions", {})
CONDITION_COLUMN = _cond_cfg.get("column_name", "condition")
CONTROL_LABEL = _cond_cfg.get("control_label", "Control")
DISEASE_LABELS = _cond_cfg.get("disease_labels", [])

# ── QC & Preprocessing ──────────────────────────────────────────
_prep_cfg = _cfg.get("preprocessing", {})
MIN_CELLS_PER_DONOR = _prep_cfg.get("min_cells_per_donor", 100)
GENE_DETECTION_MIN_FRAC = _prep_cfg.get("gene_detection_min_frac", 0.01)
NORMALIZATION = _prep_cfg.get("normalization", "log1p_cpm")

# ── PCA ──────────────────────────────────────────────────────────
_pca_cfg = _cfg.get("pca", {})
PCA_K = _pca_cfg.get("k", 100)
PCA_N_OVERSAMPLES = _pca_cfg.get("n_oversamples", 10)
PCA_N_ITER = _pca_cfg.get("n_iter", 5)
HVG_N = _pca_cfg.get("hvg_n", None)

# ── Pseudotime ───────────────────────────────────────────────────
_pt_cfg = _cfg.get("pseudotime", {})
PT_METHOD = _pt_cfg.get("method", "diffusion_map")
PT_EIGENVECTOR = _pt_cfg.get("eigenvector_index", 1)
PT_SIGN_CONVENTION = _pt_cfg.get("sign_convention", "control_low")

# ── Lane A ───────────────────────────────────────────────────────
_la_cfg = _cfg.get("lane_a", {})
N_WINDOWS = _la_cfg.get("n_windows", 8)
MIN_DONORS_PER_WINDOW = _la_cfg.get("min_donors_per_window", 3)
SPD_REPRESENTATIVE = "log_euclidean_arithmetic_mean"

# ── Hodge Quality ────────────────────────────────────────────────
_hodge_cfg = _cfg.get("hodge", {})
GRADIENT_FRACTION_THRESHOLD = _hodge_cfg.get("gradient_fraction_threshold", 0.5)
CURL_FRACTION_WARNING = _hodge_cfg.get("curl_fraction_warning", 0.3)

# ── Bootstrap ────────────────────────────────────────────────────
_boot_cfg = _cfg.get("bootstrap", {})
BOOTSTRAP_N = _boot_cfg.get("n_iterations", 100)
BOOTSTRAP_TOP_K = _boot_cfg.get("top_k", 3)
BOOTSTRAP_STABLE_THRESHOLD = _boot_cfg.get("stable_threshold", 0.80)
BOOTSTRAP_SUGGESTIVE_THRESHOLD = _boot_cfg.get("suggestive_threshold", 0.50)
BOOTSTRAP_N_WORKERS = _boot_cfg.get("n_workers", 1)
BOOTSTRAP_BASE_SEED = _boot_cfg.get("base_seed", 42)

# ── Lane B ───────────────────────────────────────────────────────
_lb_cfg = _cfg.get("lane_b", {})
UPSTREAM_PC_TOP_N = _lb_cfg.get("upstream_pc_top_n", 10)
GENE_LOADING_TOP_N = _lb_cfg.get("gene_loading_top_n", 200)

# ── Gene Hodge ───────────────────────────────────────────────────
_gh_cfg = _cfg.get("gene_hodge", {})
GENE_HODGE_PERMUTATION_N = _gh_cfg.get("permutation_n", 1000)
GENE_HODGE_HIGH_SIGMA = _gh_cfg.get("high_sigma", 1.0)
GENE_HODGE_MEDIUM_SIGMA = _gh_cfg.get("medium_sigma", 0.0)
GENE_HODGE_MIN_GENES = _gh_cfg.get("min_genes", 10)
GENE_HODGE_FLOW_MODE = _gh_cfg.get("flow_mode", "sign")
GENE_HODGE_TIE_EPSILON_FACTOR = 0.01
GENE_HODGE_VALID_FLOW_MODES = ("sign", "weighted", "edge_product", "edge_weight")

# ── Rank Shift ───────────────────────────────────────────────────
_rs_cfg = _cfg.get("rank_shift", {})
RANK_SHIFT_MIN_DONORS_PER_WINDOW = _rs_cfg.get("min_donors_per_window", 1)
RANK_SHIFT_PERMUTATION_N = _rs_cfg.get("permutation_n", 1000)
RANK_SHIFT_ALPHA = _rs_cfg.get("alpha", 0.05)
RANK_SHIFT_MIN_VALID_TRANSITIONS = 2

# ── Permutation Test ─────────────────────────────────────────────
_perm_cfg = _cfg.get("permutation", {})
PERMUTATION_N = _perm_cfg.get("n", 1000)
PERMUTATION_ALPHA = _perm_cfg.get("alpha", 0.05)

# ── Global Seed ──────────────────────────────────────────────────
GLOBAL_SEED = _cfg.get("seed", 42)

# ── Seed Management ──────────────────────────────────────────────
def compute_seed(run_id: str) -> int:
    """Deterministic seed from run_id."""
    return int(hashlib.sha256(run_id.encode()).hexdigest(), 16) % (2**32)

def make_run_id() -> str:
    """Generate run_id with timestamp."""
    from datetime import datetime
    ts = datetime.now().strftime("%Y%m%d_%H%M")
    return f"run_{ts}"

# ── Output Directory Builder ─────────────────────────────────────
def setup_run_dirs(run_id: str) -> dict:
    """Create the full output directory structure for a run."""
    base = RESULTS_DIR / run_id
    dirs = {
        "root": base,
        "logs": base / "logs",
        "pt": base / "pt",
        "laneA": base / "laneA",
        "hodge_quality": base / "laneA" / "hodge_quality",
        "laneB": base / "laneB",
        "gene_hodge": base / "laneB" / "gene_hodge",
        "rank_shift": base / "laneB" / "rank_shift",
        "bootstrap": base / "bootstrap",
        "qc": base / "qc",
        "figures": base / "figures",
    }
    for d in dirs.values():
        d.mkdir(parents=True, exist_ok=True)
    return dirs

# ── Validation ───────────────────────────────────────────────────
def validate_config() -> List[str]:
    """Validate configuration and return list of warnings/errors."""
    issues = []
    if not H5AD_DIR or not H5AD_DIR.exists():
        issues.append(f"h5ad directory not found: {H5AD_DIR}")
    if not SAMPLE_INFO_PATH or not SAMPLE_INFO_PATH.exists():
        issues.append(f"sample_info.csv not found: {SAMPLE_INFO_PATH}")
    if not CELLTYPE_LIST:
        issues.append("No cell types defined in project_config.yaml")
    if not CONTROL_LABEL:
        issues.append("No control_label defined in project_config.yaml")
    if not PT_CELLTYPE_SET:
        issues.append("No pseudotime cell type set defined")
    return issues
