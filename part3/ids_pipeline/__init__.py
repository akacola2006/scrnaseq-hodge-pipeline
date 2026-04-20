"""
IDS Pipeline - ALS Motor Cortex Causal Analysis Package
========================================================

Unified pipeline for φ-flow energy analysis using IDS Core methodology.

Usage:
    from ids_pipeline import IDSPipeline

    pipeline = IDSPipeline()  # Uses default config
    # or
    pipeline = IDSPipeline(config_path='my_config.yaml')

    # Run full analysis
    results = pipeline.run(
        expression_file='expression.csv',
        metadata_file='metadata.csv',
        output_dir='results/'
    )

Author: Claude Code + User
Date: 2025-12-11
"""

from .pipeline import IDSPipeline
from .config import IDSConfig
from .gene_utils import GeneMapper
from .module_utils import ModuleScorer

__version__ = '1.0.0'
__all__ = ['IDSPipeline', 'IDSConfig', 'GeneMapper', 'ModuleScorer']
