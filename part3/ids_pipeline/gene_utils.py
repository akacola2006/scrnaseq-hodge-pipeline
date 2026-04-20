"""
Gene Mapping Utilities
======================

Handles conversion between gene identifiers (ENSG <-> Symbol).
"""

import json
from pathlib import Path
from typing import Dict, List, Optional, Union
import pandas as pd


class GeneMapper:
    """
    Gene identifier conversion utility.

    Provides bidirectional mapping between:
    - Ensembl Gene IDs (ENSG00000000003)
    - Gene Symbols (TSPAN6)

    Examples
    --------
    >>> mapper = GeneMapper()  # Auto-loads default mapping
    >>> mapper = GeneMapper('/path/to/gene_mapping.json')

    >>> mapper.symbol_to_ensg('STMN2')
    'ENSG00000104435'

    >>> mapper.ensg_to_symbol('ENSG00000104435')
    'STMN2'

    >>> # Convert list of genes
    >>> mapper.convert_genes(['STMN2', 'TARDBP'], to='ensg')
    ['ENSG00000104435', 'ENSG00000120948']

    >>> # Convert DataFrame columns
    >>> df = mapper.convert_dataframe_genes(df, gene_column='gene_id', to='symbol')
    """

    def __init__(self, mapping_file: Optional[Union[str, Path]] = None):
        """
        Initialize GeneMapper.

        Parameters
        ----------
        mapping_file : str or Path, optional
            Path to gene mapping JSON file.
            If None, attempts to auto-detect from common locations.
        """
        self.ensg_to_sym: Dict[str, str] = {}
        self.sym_to_ensg: Dict[str, str] = {}

        if mapping_file is not None:
            self._load_mapping(Path(mapping_file))
        else:
            self._auto_load_mapping()

    def _auto_load_mapping(self):
        """Auto-detect and load gene mapping file."""
        # Search common locations
        search_paths = [
            Path(__file__).parent.parent.parent.parent / 'gene_mapping_enhanced.json',
            Path(__file__).parent.parent / 'config' / 'gene_mapping_enhanced.json',
            Path(__file__).parent.parent / 'data' / 'gene_mapping_enhanced.json',
            Path.home() / 'als' / 'motor_cortex_analysis' / 'gene_mapping_enhanced.json',
        ]

        for path in search_paths:
            if path.exists():
                self._load_mapping(path)
                return

        print("Warning: Gene mapping file not found. Using empty mapping.")
        print("Searched locations:")
        for path in search_paths:
            print(f"  - {path}")

    def _load_mapping(self, mapping_file: Path):
        """Load gene mapping from JSON file."""
        with open(mapping_file, 'r') as f:
            data = json.load(f)

        # Handle different JSON formats
        if 'ensg_to_symbol' in data:
            self.ensg_to_sym = data['ensg_to_symbol']
        else:
            # Assume flat dict format (ENSG -> Symbol)
            self.ensg_to_sym = data

        # Build reverse mapping
        self.sym_to_ensg = {v: k for k, v in self.ensg_to_sym.items()}

        print(f"Loaded {len(self.ensg_to_sym):,} gene mappings from {mapping_file}")

    def symbol_to_ensg(self, symbol: str) -> Optional[str]:
        """
        Convert gene symbol to Ensembl ID.

        Parameters
        ----------
        symbol : str
            Gene symbol (e.g., 'STMN2')

        Returns
        -------
        str or None
            Ensembl ID or None if not found
        """
        return self.sym_to_ensg.get(symbol.upper(), None)

    def ensg_to_symbol(self, ensg: str) -> Optional[str]:
        """
        Convert Ensembl ID to gene symbol.

        Parameters
        ----------
        ensg : str
            Ensembl ID (e.g., 'ENSG00000104435')

        Returns
        -------
        str or None
            Gene symbol or None if not found
        """
        return self.ensg_to_sym.get(ensg, None)

    def convert_genes(
        self,
        genes: List[str],
        to: str = 'symbol',
        keep_missing: bool = True
    ) -> List[Optional[str]]:
        """
        Convert list of gene identifiers.

        Parameters
        ----------
        genes : list of str
            List of gene identifiers
        to : str
            Target format: 'symbol' or 'ensg'
        keep_missing : bool
            If True, keep original ID for missing mappings.
            If False, return None for missing.

        Returns
        -------
        list
            Converted gene identifiers
        """
        if to == 'symbol':
            converter = self.ensg_to_symbol
        elif to == 'ensg':
            converter = self.symbol_to_ensg
        else:
            raise ValueError(f"Unknown target format: {to}. Use 'symbol' or 'ensg'.")

        results = []
        for gene in genes:
            converted = converter(gene)
            if converted is None and keep_missing:
                converted = gene
            results.append(converted)

        return results

    def convert_dataframe_genes(
        self,
        df: pd.DataFrame,
        gene_column: str = 'gene_id',
        to: str = 'symbol',
        inplace: bool = False
    ) -> pd.DataFrame:
        """
        Convert gene column in DataFrame.

        Parameters
        ----------
        df : pd.DataFrame
            Input DataFrame
        gene_column : str
            Column containing gene identifiers
        to : str
            Target format: 'symbol' or 'ensg'
        inplace : bool
            If True, modify DataFrame in place

        Returns
        -------
        pd.DataFrame
            DataFrame with converted gene column
        """
        if not inplace:
            df = df.copy()

        df[gene_column] = self.convert_genes(
            df[gene_column].tolist(),
            to=to,
            keep_missing=True
        )

        return df

    def convert_expression_matrix_index(
        self,
        df: pd.DataFrame,
        to: str = 'symbol'
    ) -> pd.DataFrame:
        """
        Convert gene index of expression matrix.

        Parameters
        ----------
        df : pd.DataFrame
            Expression matrix with genes as index
        to : str
            Target format: 'symbol' or 'ensg'

        Returns
        -------
        pd.DataFrame
            Expression matrix with converted index
        """
        df = df.copy()
        df.index = self.convert_genes(df.index.tolist(), to=to, keep_missing=True)
        return df

    def get_ensg_list(self, symbols: List[str]) -> List[str]:
        """
        Get list of Ensembl IDs for gene symbols.

        Skips genes not found in mapping.

        Parameters
        ----------
        symbols : list of str
            Gene symbols

        Returns
        -------
        list of str
            Ensembl IDs (only found genes)
        """
        result = []
        for sym in symbols:
            ensg = self.symbol_to_ensg(sym.upper())
            if ensg:
                result.append(ensg)
        return result

    def get_symbol_list(self, ensgs: List[str]) -> List[str]:
        """
        Get list of gene symbols for Ensembl IDs.

        Skips genes not found in mapping.

        Parameters
        ----------
        ensgs : list of str
            Ensembl IDs

        Returns
        -------
        list of str
            Gene symbols (only found genes)
        """
        result = []
        for ensg in ensgs:
            sym = self.ensg_to_symbol(ensg)
            if sym:
                result.append(sym)
        return result

    def detect_id_type(self, gene_id: str) -> str:
        """
        Detect the type of gene identifier.

        Parameters
        ----------
        gene_id : str
            Gene identifier

        Returns
        -------
        str
            'ensg', 'symbol', or 'unknown'
        """
        if gene_id.startswith('ENSG'):
            return 'ensg'
        elif gene_id in self.sym_to_ensg:
            return 'symbol'
        elif gene_id in self.ensg_to_sym:
            return 'ensg'
        else:
            return 'unknown'

    def __len__(self) -> int:
        """Return number of gene mappings."""
        return len(self.ensg_to_sym)

    def __contains__(self, gene: str) -> bool:
        """Check if gene is in mapping."""
        return gene in self.ensg_to_sym or gene in self.sym_to_ensg


# Singleton instance for convenience
_default_mapper = None


def get_mapper() -> GeneMapper:
    """Get default GeneMapper instance (singleton)."""
    global _default_mapper
    if _default_mapper is None:
        _default_mapper = GeneMapper()
    return _default_mapper
