"""Seed management for reproducibility."""

import hashlib
import random

import numpy as np


def set_global_seed(seed: int) -> None:
    """Set random seeds for Python, NumPy, and PyTorch (CPU + CUDA)."""
    random.seed(seed)
    np.random.seed(seed)
    try:
        import torch
        torch.manual_seed(seed)
        if torch.cuda.is_available():
            torch.cuda.manual_seed(seed)
            torch.cuda.manual_seed_all(seed)
            torch.backends.cudnn.deterministic = True
            torch.backends.cudnn.benchmark = False
    except ImportError:
        pass


def get_bootstrap_seed(base_seed: int, iteration: int) -> int:
    """Derive a deterministic seed for a specific bootstrap iteration."""
    key = f"{base_seed}:{iteration}".encode()
    h = hashlib.sha256(key).hexdigest()
    return int(h, 16) % (2**32)
