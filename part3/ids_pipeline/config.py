"""
IDS Pipeline Configuration
==========================

Central configuration management for the IDS analysis pipeline.
Handles file paths, gene mappings, module definitions, and analysis parameters.
"""

import os
import json
import yaml
from pathlib import Path
from typing import Dict, List, Optional, Any
from dataclasses import dataclass, field


@dataclass
class IDSConfig:
    """
    Central configuration for IDS analysis pipeline.

    Can be initialized from:
    - Default values (auto-detects project structure)
    - YAML configuration file
    - Dictionary

    Examples
    --------
    >>> config = IDSConfig()  # Use defaults
    >>> config = IDSConfig.from_yaml('config.yaml')
    >>> config = IDSConfig.from_dict({'n_bins': 30})
    """

    # ==========================================================================
    # Directory Paths
    # ==========================================================================

    # Base project directory (auto-detected)
    base_dir: Path = field(default_factory=lambda: Path(__file__).parent.parent.parent)

    # Subdirectories
    data_dir: Optional[Path] = None
    results_dir: Optional[Path] = None
    config_dir: Optional[Path] = None

    # ==========================================================================
    # Data Files
    # ==========================================================================

    # Gene mapping file (ENSG <-> Symbol)
    gene_mapping_file: Optional[Path] = None

    # Module definition file (JSON)
    module_file: Optional[Path] = None

    # ==========================================================================
    # NVU Group Mapping
    # ==========================================================================

    nvu_mapping: Dict[str, List[str]] = field(default_factory=lambda: {
        'Vascular': [
            'Vasc.Endo', 'Vasc.Fibro', 'Vasc.Mural', 'Vasc.Pericyte', 'Vasc.SMC',
            'Endo', 'VLMC', 'Peri', 'Immune'
        ],
        'Glia': [
            'Glia.Astro', 'Glia.Oligo', 'Glia.Micro', 'Glia.OPC',
            'Astro', 'Oligo', 'OPC', 'Micro'
        ],
        'Excitatory': [
            'Ex.L2', 'Ex.L3', 'Ex.L4', 'Ex.L5', 'Ex.L6',
            'Ex.RORB', 'Ex.THEMIS', 'Ex.OPRK1'
        ],
        'Inhibitory': [
            'In.', 'Inh.', 'VIP', 'PVALB', 'SST', 'LAMP5', 'PAX6'
        ],
        'VAT1L': [
            'VAT1L', 'Ex.L5.VAT1L'
        ]
    })

    # ==========================================================================
    # Analysis Parameters
    # ==========================================================================

    # Pseudotime binning
    n_bins: int = 25

    # φ threshold for onset detection
    phi_threshold: float = 0.5

    # Minimum cells per group/bin for analysis
    min_cells: int = 10

    # Statistical parameters
    eps: float = 1e-8

    # ==========================================================================
    # ALS Stress Components (weights must sum to 1.0)
    # ==========================================================================

    stress_components: Dict[str, Dict[str, Any]] = field(default_factory=lambda: {
        'TDP43_targets': {
            'weight': 0.30,
            'genes': ['STMN2', 'UNC13A', 'ELAVL4', 'KCNQ2', 'PFKP',
                     'KALRN', 'SYT1', 'GRIN2B', 'CACNA1A', 'SCN2A']
        },
        'STMN2_pathway': {
            'weight': 0.20,
            'genes': ['STMN1', 'STMN2', 'STMN3', 'STMN4', 'MAPT',
                     'MAP2', 'MAP1B', 'DPYSL2', 'DPYSL3']
        },
        'stress_granule': {
            'weight': 0.15,
            'genes': ['G3BP1', 'G3BP2', 'TIA1', 'TIAL1', 'FUS',
                     'TARDBP', 'HNRNPA1', 'HNRNPA2B1', 'ATXN2']
        },
        'ER_stress': {
            'weight': 0.15,
            'genes': ['HSPA5', 'DDIT3', 'ATF4', 'XBP1', 'ATF6',
                     'ERN1', 'EIF2AK3', 'PDIA4', 'CALR']
        },
        'cytoskeleton': {
            'weight': 0.10,
            'genes': ['NEFL', 'NEFM', 'NEFH', 'TUBB3', 'TUBA1A',
                     'ACTB', 'ACTG1', 'MAP1A', 'MAP4', 'MAPT', 'VIM']
        },
        'inflammation': {
            'weight': 0.10,
            'genes': ['IL1B', 'IL6', 'TNF', 'C1QA', 'C1QB', 'C1QC',
                     'C3', 'CD68', 'TREM2', 'TYROBP', 'AIF1', 'CX3CR1', 'P2RY12']
        }
    })

    # ==========================================================================
    # Output Settings
    # ==========================================================================

    # Figure DPI
    figure_dpi: int = 150

    # Save intermediate files
    save_intermediate: bool = True

    # Verbose output
    verbose: bool = True

    def __post_init__(self):
        """Auto-detect paths after initialization."""
        self.base_dir = Path(self.base_dir)

        # Auto-detect directories
        if self.data_dir is None:
            # Try multiple locations
            candidates = [
                self.base_dir.parent,  # motor_cortex_analysis
                self.base_dir / 'data',
                self.base_dir
            ]
            for cand in candidates:
                if cand.exists():
                    self.data_dir = cand
                    break
        else:
            self.data_dir = Path(self.data_dir)

        if self.results_dir is None:
            self.results_dir = self.base_dir / 'results'
        else:
            self.results_dir = Path(self.results_dir)

        if self.config_dir is None:
            self.config_dir = self.base_dir / 'config'
        else:
            self.config_dir = Path(self.config_dir)

        # Auto-detect gene mapping file
        if self.gene_mapping_file is None:
            candidates = [
                self.base_dir / 'gene_mapping_enhanced.json',
                self.base_dir.parent / 'gene_mapping_enhanced.json',
                self.data_dir / 'gene_mapping_enhanced.json' if self.data_dir else None,
                self.config_dir / 'gene_mapping_enhanced.json' if self.config_dir else None,
            ]
            for cand in candidates:
                if cand and cand.exists():
                    self.gene_mapping_file = cand
                    break
        else:
            self.gene_mapping_file = Path(self.gene_mapping_file)

        # Auto-detect module file
        if self.module_file is None:
            candidates = [
                self.base_dir / '24_functional_modules_fixed.json',
                self.base_dir.parent / '24_functional_modules_fixed.json',
                self.data_dir / '24_functional_modules_fixed.json' if self.data_dir else None,
                self.config_dir / '24_functional_modules_fixed.json' if self.config_dir else None,
            ]
            for cand in candidates:
                if cand and cand.exists():
                    self.module_file = cand
                    break
        else:
            self.module_file = Path(self.module_file)

    @classmethod
    def from_yaml(cls, yaml_path: str) -> 'IDSConfig':
        """Load configuration from YAML file."""
        with open(yaml_path, 'r') as f:
            config_dict = yaml.safe_load(f)
        return cls.from_dict(config_dict)

    @classmethod
    def from_dict(cls, config_dict: Dict[str, Any]) -> 'IDSConfig':
        """Create configuration from dictionary."""
        # Convert path strings to Path objects
        path_fields = ['base_dir', 'data_dir', 'results_dir', 'config_dir',
                       'gene_mapping_file', 'module_file']
        for field in path_fields:
            if field in config_dict and config_dict[field] is not None:
                config_dict[field] = Path(config_dict[field])

        return cls(**config_dict)

    def to_dict(self) -> Dict[str, Any]:
        """Convert configuration to dictionary."""
        result = {}
        for key, value in self.__dict__.items():
            if isinstance(value, Path):
                result[key] = str(value)
            else:
                result[key] = value
        return result

    def to_yaml(self, yaml_path: str):
        """Save configuration to YAML file."""
        with open(yaml_path, 'w') as f:
            yaml.dump(self.to_dict(), f, default_flow_style=False)

    def get_nvu_group(self, cell_type: str) -> str:
        """
        Map a cell type to its NVU group.

        Parameters
        ----------
        cell_type : str
            Cell type name (e.g., 'Glia.Astro.GFAP-pos')

        Returns
        -------
        str
            NVU group name (e.g., 'Glia') or 'Other' if not matched
        """
        for group, patterns in self.nvu_mapping.items():
            for pattern in patterns:
                if pattern in cell_type:
                    return group
        return 'Other'

    def validate(self) -> List[str]:
        """
        Validate configuration and return list of issues.

        Returns
        -------
        List[str]
            List of validation issues (empty if valid)
        """
        issues = []

        # Check directories
        if self.data_dir and not self.data_dir.exists():
            issues.append(f"Data directory not found: {self.data_dir}")

        # Check gene mapping
        if self.gene_mapping_file is None:
            issues.append("Gene mapping file not found (searched multiple locations)")
        elif not self.gene_mapping_file.exists():
            issues.append(f"Gene mapping file not found: {self.gene_mapping_file}")

        # Check module file
        if self.module_file is None:
            issues.append("Module file not found (searched multiple locations)")
        elif not self.module_file.exists():
            issues.append(f"Module file not found: {self.module_file}")

        # Check stress component weights
        total_weight = sum(c['weight'] for c in self.stress_components.values())
        if abs(total_weight - 1.0) > 0.01:
            issues.append(f"Stress component weights sum to {total_weight}, should be 1.0")

        return issues

    def print_summary(self):
        """Print configuration summary."""
        print("=" * 60)
        print("IDS Pipeline Configuration")
        print("=" * 60)
        print(f"\nDirectories:")
        print(f"  Base:    {self.base_dir}")
        print(f"  Data:    {self.data_dir}")
        print(f"  Results: {self.results_dir}")
        print(f"  Config:  {self.config_dir}")

        print(f"\nData Files:")
        print(f"  Gene Mapping: {self.gene_mapping_file}")
        print(f"  Module File:  {self.module_file}")

        print(f"\nAnalysis Parameters:")
        print(f"  PT bins:       {self.n_bins}")
        print(f"  φ threshold:   {self.phi_threshold}")
        print(f"  Min cells:     {self.min_cells}")

        print(f"\nNVU Groups:")
        for group, patterns in self.nvu_mapping.items():
            print(f"  {group}: {len(patterns)} patterns")

        print(f"\nStress Components:")
        for name, comp in self.stress_components.items():
            print(f"  {name}: {comp['weight']*100:.0f}% ({len(comp['genes'])} genes)")

        # Validation
        issues = self.validate()
        if issues:
            print(f"\n⚠️  Validation Issues:")
            for issue in issues:
                print(f"  - {issue}")
        else:
            print(f"\n✓ Configuration valid")

        print("=" * 60)


# Default configuration instance
DEFAULT_CONFIG = IDSConfig()
