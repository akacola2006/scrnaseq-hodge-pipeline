"""
IDS Pipeline - Main Pipeline Runner
====================================

Unified pipeline for φ-flow energy analysis.
"""

import sys
import numpy as np
import pandas as pd
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Optional, Tuple, Union
import warnings

# Import local modules
from .config import IDSConfig
from .gene_utils import GeneMapper
from .module_utils import ModuleScorer

# Import ids_core functions (if available)
try:
    sys.path.insert(0, str(Path(__file__).parent.parent / 'scripts'))
    from ids_core import (
        compute_phi_from_scores,
        bin_by_pt,
        summarize_phi_by_group_and_bin,
        detect_onset_and_peak,
        compute_phi_flow,
        map_cell_types_to_groups,
        build_module_feature_matrix,
        discover_axes_clustering,
        discover_axes_factorization
    )
    IDS_CORE_AVAILABLE = True
except ImportError:
    IDS_CORE_AVAILABLE = False
    warnings.warn("ids_core.py not found. Some functions will be limited.")


class IDSPipeline:
    """
    Unified IDS Analysis Pipeline.

    Provides a simple interface to run the complete φ-flow analysis:
    1. Load expression data
    2. Compute module scores
    3. Compute φ energy
    4. Bin by pseudotime
    5. Detect onset/peak timing
    6. Discover disease axes

    Examples
    --------
    >>> from ids_pipeline import IDSPipeline

    >>> # Initialize with default config
    >>> pipeline = IDSPipeline()

    >>> # Or with custom config
    >>> pipeline = IDSPipeline(config_path='my_config.yaml')

    >>> # Run full analysis
    >>> results = pipeline.run(
    ...     expression_file='expression.csv',
    ...     metadata_file='metadata.csv',
    ...     output_dir='results/'
    ... )

    >>> # Or run step by step
    >>> pipeline.load_data('expression.csv', 'metadata.csv')
    >>> pipeline.compute_modules()
    >>> pipeline.compute_phi()
    >>> pipeline.analyze_flow()
    >>> pipeline.save_results('results/')
    """

    def __init__(
        self,
        config: Optional[IDSConfig] = None,
        config_path: Optional[str] = None
    ):
        """
        Initialize IDS Pipeline.

        Parameters
        ----------
        config : IDSConfig, optional
            Configuration object
        config_path : str, optional
            Path to YAML config file
        """
        # Load configuration
        if config is not None:
            self.config = config
        elif config_path is not None:
            self.config = IDSConfig.from_yaml(config_path)
        else:
            self.config = IDSConfig()

        # Initialize utilities
        self.gene_mapper = GeneMapper(self.config.gene_mapping_file)
        self.module_scorer = ModuleScorer(
            self.config.module_file,
            gene_mapper=self.gene_mapper
        )

        # Data containers
        self.expression: Optional[pd.DataFrame] = None
        self.metadata: Optional[pd.DataFrame] = None
        self.module_scores: Optional[pd.DataFrame] = None
        self.phi_values: Optional[pd.DataFrame] = None
        self.flow_profiles: Optional[pd.DataFrame] = None
        self.flow_summary: Optional[pd.DataFrame] = None

        # Validate configuration
        issues = self.config.validate()
        if issues:
            print("⚠️  Configuration issues:")
            for issue in issues:
                print(f"  - {issue}")

    def load_data(
        self,
        expression_file: Union[str, Path],
        metadata_file: Union[str, Path],
        cell_column: str = 'cell_id',
        condition_column: str = 'condition',
        celltype_column: str = 'cell_type',
        pt_column: Optional[str] = None
    ) -> 'IDSPipeline':
        """
        Load expression and metadata files.

        Parameters
        ----------
        expression_file : str or Path
            Gene expression matrix (genes × cells or cells × genes)
        metadata_file : str or Path
            Cell metadata with condition, cell type, etc.
        cell_column : str
            Column name for cell IDs
        condition_column : str
            Column name for condition (Control/ALS)
        celltype_column : str
            Column name for cell types
        pt_column : str, optional
            Column name for pseudotime (if pre-computed)

        Returns
        -------
        self
            For method chaining
        """
        expression_file = Path(expression_file)
        metadata_file = Path(metadata_file)

        print(f"\n{'='*60}")
        print("Loading Data")
        print(f"{'='*60}")

        # Load expression
        print(f"\nExpression: {expression_file}")
        if expression_file.suffix == '.gz':
            import gzip
            with gzip.open(expression_file, 'rt') as f:
                self.expression = pd.read_csv(f, index_col=0)
        else:
            self.expression = pd.read_csv(expression_file, index_col=0)

        # Auto-detect orientation (genes × cells vs cells × genes)
        if self.expression.shape[0] > self.expression.shape[1]:
            # Likely genes × cells (more genes than cells is unusual)
            # Check if index looks like genes
            if not self.expression.index[0].startswith('ENSG'):
                # Might need transpose
                print("  Note: Auto-detecting matrix orientation...")

        print(f"  Shape: {self.expression.shape}")

        # Load metadata
        print(f"\nMetadata: {metadata_file}")
        self.metadata = pd.read_csv(metadata_file)
        print(f"  Shape: {self.metadata.shape}")

        # Standardize column names
        self._cell_column = cell_column
        self._condition_column = condition_column
        self._celltype_column = celltype_column
        self._pt_column = pt_column

        # Verify columns exist
        required_cols = [cell_column, condition_column, celltype_column]
        missing = [c for c in required_cols if c not in self.metadata.columns]
        if missing:
            print(f"\n⚠️  Missing columns: {missing}")
            print(f"   Available: {list(self.metadata.columns)}")

        # Add NVU group mapping
        if celltype_column in self.metadata.columns:
            self.metadata['nvu_group'] = self.metadata[celltype_column].apply(
                self.config.get_nvu_group
            )
            nvu_counts = self.metadata['nvu_group'].value_counts()
            print(f"\nNVU Groups:")
            for group, count in nvu_counts.items():
                print(f"  {group}: {count:,}")

        print(f"\n✓ Data loaded successfully")

        return self

    def compute_modules(
        self,
        modules: Optional[List[str]] = None,
        method: str = 'mean'
    ) -> 'IDSPipeline':
        """
        Compute module scores from expression data.

        Parameters
        ----------
        modules : list of str, optional
            Modules to score. If None, score all modules.
        method : str
            Scoring method: 'mean' or 'sum'

        Returns
        -------
        self
            For method chaining
        """
        if self.expression is None:
            raise ValueError("Expression data not loaded. Call load_data() first.")

        print(f"\n{'='*60}")
        print("Computing Module Scores")
        print(f"{'='*60}")

        self.module_scores = self.module_scorer.score_expression(
            self.expression,
            modules=modules,
            method=method
        )

        # Merge with metadata
        if self.metadata is not None:
            # Align indices
            common_cells = self.module_scores.index.intersection(
                self.metadata[self._cell_column] if self._cell_column in self.metadata.columns
                else self.metadata.index
            )
            print(f"\n✓ {len(common_cells)} cells with module scores and metadata")

        return self

    def compute_phi(
        self,
        control_value: str = 'Control',
        modules: Optional[List[str]] = None
    ) -> 'IDSPipeline':
        """
        Compute φ (phi) energy for each module.

        Parameters
        ----------
        control_value : str
            Value in condition column indicating Control cells
        modules : list of str, optional
            Modules to analyze. If None, use all scored modules.

        Returns
        -------
        self
            For method chaining
        """
        if self.module_scores is None:
            raise ValueError("Module scores not computed. Call compute_modules() first.")

        if not IDS_CORE_AVAILABLE:
            raise ImportError("ids_core.py required for φ computation")

        print(f"\n{'='*60}")
        print("Computing φ Energy")
        print(f"{'='*60}")

        # Get module columns
        if modules is None:
            module_cols = [c for c in self.module_scores.columns if c.startswith('module_')]
        else:
            module_cols = [f'module_{m}' for m in modules if f'module_{m}' in self.module_scores.columns]

        print(f"\nModules: {len(module_cols)}")

        # Create control mask
        control_mask = self.metadata[self._condition_column] == control_value
        print(f"Control cells: {control_mask.sum():,}")

        # Compute φ for each module
        phi_results = {}
        for col in module_cols:
            phi = compute_phi_from_scores(
                scores=self.module_scores[col],
                control_mask=control_mask,
                eps=self.config.eps
            )
            phi_results[col.replace('module_', 'phi_')] = phi

        self.phi_values = pd.DataFrame(phi_results)

        print(f"\n✓ Computed φ for {len(phi_results)} modules")

        return self

    def analyze_flow(
        self,
        group_column: str = 'nvu_group',
        pt_column: Optional[str] = None
    ) -> 'IDSPipeline':
        """
        Analyze φ-flow across pseudotime.

        Parameters
        ----------
        group_column : str
            Column for grouping (e.g., 'nvu_group', 'cell_type')
        pt_column : str, optional
            Column for pseudotime. If None, uses config setting.

        Returns
        -------
        self
            For method chaining
        """
        if self.phi_values is None:
            raise ValueError("φ values not computed. Call compute_phi() first.")

        if not IDS_CORE_AVAILABLE:
            raise ImportError("ids_core.py required for flow analysis")

        print(f"\n{'='*60}")
        print("Analyzing φ-Flow")
        print(f"{'='*60}")

        # Get pseudotime column
        if pt_column is None:
            pt_column = self._pt_column

        if pt_column is None or pt_column not in self.metadata.columns:
            print("⚠️  Pseudotime column not found. Skipping flow analysis.")
            print(f"   Available columns: {list(self.metadata.columns)}")
            return self

        pt = self.metadata[pt_column]
        group = self.metadata[group_column]

        # Bin by pseudotime (once for all modules)
        bin_index, bin_edges, bin_centers = bin_by_pt(pt, n_bins=self.config.n_bins)

        # Analyze each module
        all_profiles = []
        all_flows = []

        phi_cols = [c for c in self.phi_values.columns if c.startswith('phi_')]

        for phi_col in phi_cols:
            module_name = phi_col.replace('phi_', '')

            # Use pre-computed φ values directly (NOT compute_phi_flow which would re-compute)
            # This avoids the φ² bug
            profiles = summarize_phi_by_group_and_bin(
                phi=self.phi_values[phi_col],
                group=group,
                bin_index=bin_index
            )

            # Detect onset and peak
            phi_mean_df = profiles[['group', 'bin', 'phi_mean']].copy()
            flow = detect_onset_and_peak(
                phi_mean_df=phi_mean_df,
                bin_centers=bin_centers,
                threshold=self.config.phi_threshold
            )

            profiles['module'] = module_name
            flow['module'] = module_name

            all_profiles.append(profiles)
            all_flows.append(flow)

        self.flow_profiles = pd.concat(all_profiles, ignore_index=True)
        self.flow_summary = pd.concat(all_flows, ignore_index=True)

        print(f"\n✓ Analyzed {len(phi_cols)} modules across {len(group.unique())} groups")

        # Print summary
        print(f"\nOnset Summary (first bins):")
        onset_summary = self.flow_summary.groupby('group')['onset_PT'].mean()
        for grp, onset in onset_summary.sort_values().items():
            print(f"  {grp}: PT = {onset:.3f}")

        return self

    def save_results(
        self,
        output_dir: Union[str, Path],
        prefix: str = ''
    ) -> 'IDSPipeline':
        """
        Save all results to output directory.

        Parameters
        ----------
        output_dir : str or Path
            Output directory
        prefix : str
            Prefix for output filenames

        Returns
        -------
        self
            For method chaining
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        print(f"\n{'='*60}")
        print(f"Saving Results to {output_dir}")
        print(f"{'='*60}")

        timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')

        if self.module_scores is not None:
            path = output_dir / f'{prefix}module_scores.csv'
            self.module_scores.to_csv(path)
            print(f"  ✓ {path}")

        if self.phi_values is not None:
            path = output_dir / f'{prefix}phi_values.csv'
            self.phi_values.to_csv(path)
            print(f"  ✓ {path}")

        if self.flow_profiles is not None:
            path = output_dir / f'{prefix}flow_profiles.csv'
            self.flow_profiles.to_csv(path, index=False)
            print(f"  ✓ {path}")

        if self.flow_summary is not None:
            path = output_dir / f'{prefix}flow_summary.csv'
            self.flow_summary.to_csv(path, index=False)
            print(f"  ✓ {path}")

        # Save config
        config_path = output_dir / f'{prefix}config.yaml'
        self.config.to_yaml(str(config_path))
        print(f"  ✓ {config_path}")

        return self

    def run(
        self,
        expression_file: Union[str, Path],
        metadata_file: Union[str, Path],
        output_dir: Union[str, Path],
        **kwargs
    ) -> Dict:
        """
        Run complete pipeline.

        Parameters
        ----------
        expression_file : str or Path
            Gene expression matrix
        metadata_file : str or Path
            Cell metadata
        output_dir : str or Path
            Output directory
        **kwargs
            Additional arguments passed to load_data()

        Returns
        -------
        dict
            Dictionary with results DataFrames
        """
        print(f"\n{'#'*60}")
        print("IDS Pipeline - Full Analysis")
        print(f"{'#'*60}")
        print(f"\nStarted: {datetime.now()}")

        # Run steps
        self.load_data(expression_file, metadata_file, **kwargs)
        self.compute_modules()
        self.compute_phi()
        self.analyze_flow()
        self.save_results(output_dir)

        print(f"\n{'#'*60}")
        print(f"Completed: {datetime.now()}")
        print(f"{'#'*60}")

        return {
            'module_scores': self.module_scores,
            'phi_values': self.phi_values,
            'flow_profiles': self.flow_profiles,
            'flow_summary': self.flow_summary
        }

    def get_results(self) -> Dict:
        """Get all results as dictionary."""
        return {
            'module_scores': self.module_scores,
            'phi_values': self.phi_values,
            'flow_profiles': self.flow_profiles,
            'flow_summary': self.flow_summary
        }

    def print_summary(self):
        """Print analysis summary."""
        print(f"\n{'='*60}")
        print("Analysis Summary")
        print(f"{'='*60}")

        if self.expression is not None:
            print(f"\nExpression: {self.expression.shape[0]} genes × {self.expression.shape[1]} cells")

        if self.metadata is not None:
            print(f"Metadata: {len(self.metadata)} cells")

        if self.module_scores is not None:
            n_modules = len([c for c in self.module_scores.columns if c.startswith('module_')])
            print(f"Module Scores: {n_modules} modules")

        if self.phi_values is not None:
            n_phi = len([c for c in self.phi_values.columns if c.startswith('phi_')])
            print(f"φ Values: {n_phi} modules")

        if self.flow_summary is not None:
            print(f"Flow Summary: {len(self.flow_summary)} entries")

        print(f"{'='*60}")
