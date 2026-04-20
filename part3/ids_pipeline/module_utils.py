"""
Module Scoring Utilities
========================

Handles functional module score computation from gene expression data.
"""

import json
import numpy as np
import pandas as pd
from pathlib import Path
from typing import Dict, List, Optional, Union

from .gene_utils import GeneMapper


class ModuleScorer:
    """
    Compute functional module scores from gene expression data.

    Modules are defined as groups of functionally related genes.
    Score = mean expression of module genes (Z-score normalized).

    Examples
    --------
    >>> scorer = ModuleScorer()  # Auto-loads default modules
    >>> scorer = ModuleScorer('/path/to/modules.json')

    >>> # Score single cell expression
    >>> scores = scorer.score_expression(expression_df)

    >>> # Get genes for a module
    >>> genes = scorer.get_module_genes('Mitochondria')

    >>> # List all modules
    >>> modules = scorer.list_modules()
    """

    def __init__(
        self,
        module_file: Optional[Union[str, Path]] = None,
        gene_mapper: Optional[GeneMapper] = None
    ):
        """
        Initialize ModuleScorer.

        Parameters
        ----------
        module_file : str or Path, optional
            Path to module definition JSON file.
            If None, attempts to auto-detect from common locations.
        gene_mapper : GeneMapper, optional
            Gene mapper for ID conversion.
            If None, creates new instance.
        """
        self.modules: Dict[str, List[str]] = {}
        self.gene_mapper = gene_mapper or GeneMapper()

        if module_file is not None:
            self._load_modules(Path(module_file))
        else:
            self._auto_load_modules()

    def _auto_load_modules(self):
        """Auto-detect and load module file."""
        search_paths = [
            Path(__file__).parent.parent.parent.parent / '24_functional_modules_fixed.json',
            Path(__file__).parent.parent / 'config' / '24_functional_modules_fixed.json',
            Path(__file__).parent.parent / 'data' / '24_functional_modules_fixed.json',
            Path.home() / 'als' / 'motor_cortex_analysis' / '24_functional_modules_fixed.json',
        ]

        for path in search_paths:
            if path.exists():
                self._load_modules(path)
                return

        print("Warning: Module definition file not found. Using empty modules.")
        print("Searched locations:")
        for path in search_paths:
            print(f"  - {path}")

    def _load_modules(self, module_file: Path):
        """Load module definitions from JSON file."""
        with open(module_file, 'r') as f:
            data = json.load(f)

        # Handle nested format (with 'modules' key)
        if 'modules' in data:
            self.modules = data['modules']
        else:
            # Assume flat format
            self.modules = data

        # Filter out non-module keys (like 'summary')
        self.modules = {k: v for k, v in self.modules.items()
                       if isinstance(v, list)}

        total_genes = sum(len(genes) for genes in self.modules.values())
        print(f"Loaded {len(self.modules)} modules ({total_genes:,} total genes) from {module_file}")

    def list_modules(self) -> List[str]:
        """Return list of module names."""
        return list(self.modules.keys())

    def get_module_genes(
        self,
        module_name: str,
        as_ensg: bool = True
    ) -> List[str]:
        """
        Get genes for a module.

        Parameters
        ----------
        module_name : str
            Module name
        as_ensg : bool
            If True, return Ensembl IDs. If False, return symbols.

        Returns
        -------
        list of str
            Gene identifiers
        """
        if module_name not in self.modules:
            raise ValueError(f"Unknown module: {module_name}. "
                           f"Available: {self.list_modules()}")

        genes = self.modules[module_name]

        if not as_ensg:
            # Convert to symbols
            genes = self.gene_mapper.get_symbol_list(genes)

        return genes

    def score_expression(
        self,
        expression: pd.DataFrame,
        modules: Optional[List[str]] = None,
        method: str = 'mean',
        center: bool = True,
        scale: bool = True,
        min_genes: int = 3
    ) -> pd.DataFrame:
        """
        Compute module scores from expression matrix.

        Parameters
        ----------
        expression : pd.DataFrame
            Gene expression matrix.
            - Rows = genes (index can be ENSG or Symbol)
            - Columns = cells
        modules : list of str, optional
            Modules to score. If None, score all modules.
        method : str
            Scoring method: 'mean' (default) or 'sum'
        center : bool
            If True, center scores (subtract mean)
        scale : bool
            If True, scale scores (divide by std)
        min_genes : int
            Minimum genes required for a module to be scored

        Returns
        -------
        pd.DataFrame
            Module scores.
            - Rows = cells
            - Columns = module names (prefixed with 'module_')
        """
        if modules is None:
            modules = self.list_modules()

        # Detect gene ID type
        sample_gene = expression.index[0]
        is_ensg = sample_gene.startswith('ENSG')

        results = {}

        for module_name in modules:
            # Get module genes in appropriate format
            module_genes = self.get_module_genes(module_name, as_ensg=is_ensg)

            # Find genes present in expression data
            present_genes = [g for g in module_genes if g in expression.index]

            if len(present_genes) < min_genes:
                print(f"  Warning: {module_name} has only {len(present_genes)}/{len(module_genes)} "
                      f"genes present (min={min_genes}). Skipping.")
                continue

            # Extract module expression
            module_expr = expression.loc[present_genes]

            # Compute score
            if method == 'mean':
                score = module_expr.mean(axis=0)
            elif method == 'sum':
                score = module_expr.sum(axis=0)
            else:
                raise ValueError(f"Unknown method: {method}")

            # Normalize
            if center:
                score = score - score.mean()
            if scale and score.std() > 0:
                score = score / score.std()

            results[f'module_{module_name}'] = score

        scores_df = pd.DataFrame(results)

        print(f"Computed {len(results)} module scores for {len(scores_df)} cells")

        return scores_df

    def score_from_file(
        self,
        expression_file: Union[str, Path],
        output_file: Optional[Union[str, Path]] = None,
        **kwargs
    ) -> pd.DataFrame:
        """
        Score expression from file.

        Parameters
        ----------
        expression_file : str or Path
            Path to expression CSV/TSV file (genes × cells)
        output_file : str or Path, optional
            Path to save scores
        **kwargs
            Additional arguments passed to score_expression()

        Returns
        -------
        pd.DataFrame
            Module scores
        """
        expression_file = Path(expression_file)

        # Load expression data
        print(f"Loading expression from {expression_file}...")

        if expression_file.suffix == '.gz':
            import gzip
            with gzip.open(expression_file, 'rt') as f:
                expression = pd.read_csv(f, index_col=0)
        else:
            expression = pd.read_csv(expression_file, index_col=0)

        print(f"  Loaded: {expression.shape[0]} genes × {expression.shape[1]} cells")

        # Score
        scores = self.score_expression(expression, **kwargs)

        # Save if requested
        if output_file:
            scores.to_csv(output_file)
            print(f"  Saved scores to {output_file}")

        return scores

    def get_module_stats(self, expression: pd.DataFrame) -> pd.DataFrame:
        """
        Get statistics about module gene coverage in expression data.

        Parameters
        ----------
        expression : pd.DataFrame
            Expression matrix

        Returns
        -------
        pd.DataFrame
            Module statistics (genes_defined, genes_present, coverage)
        """
        sample_gene = expression.index[0]
        is_ensg = sample_gene.startswith('ENSG')

        stats = []
        for module_name in self.list_modules():
            module_genes = self.get_module_genes(module_name, as_ensg=is_ensg)
            present_genes = [g for g in module_genes if g in expression.index]

            stats.append({
                'module': module_name,
                'genes_defined': len(module_genes),
                'genes_present': len(present_genes),
                'coverage': len(present_genes) / len(module_genes) if module_genes else 0
            })

        return pd.DataFrame(stats)

    def __len__(self) -> int:
        """Return number of modules."""
        return len(self.modules)

    def __contains__(self, module: str) -> bool:
        """Check if module exists."""
        return module in self.modules


# Singleton instance
_default_scorer = None


def get_scorer() -> ModuleScorer:
    """Get default ModuleScorer instance (singleton)."""
    global _default_scorer
    if _default_scorer is None:
        _default_scorer = ModuleScorer()
    return _default_scorer
