#!/usr/bin/env python3
"""
IDS Pipeline Runner
===================

Command-line interface for the IDS analysis pipeline.

Usage:
------
# Run with default config
python run_pipeline.py expression.csv metadata.csv -o results/

# Run with custom config
python run_pipeline.py expression.csv metadata.csv -o results/ -c my_config.yaml

# Show configuration
python run_pipeline.py --show-config

# Generate config template
python run_pipeline.py --generate-config output_config.yaml

Examples:
---------
# Full analysis on ALS motor cortex data
python run_pipeline.py \\
    /path/to/expression.csv \\
    /path/to/metadata.csv \\
    -o results/my_analysis/ \\
    --condition-column condition \\
    --celltype-column cell_type \\
    --pt-column PT_dpt

Author: Claude Code + User
Date: 2025-12-11
"""

import argparse
import sys
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from ids_pipeline import IDSPipeline, IDSConfig


def main():
    parser = argparse.ArgumentParser(
        description='IDS Pipeline - φ-Flow Energy Analysis',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Run analysis
  python run_pipeline.py expression.csv metadata.csv -o results/

  # With custom config
  python run_pipeline.py expression.csv metadata.csv -o results/ -c config.yaml

  # Show configuration
  python run_pipeline.py --show-config

  # Generate config template
  python run_pipeline.py --generate-config my_config.yaml
"""
    )

    # Input files (positional, optional)
    parser.add_argument(
        'expression_file',
        nargs='?',
        help='Gene expression matrix (genes × cells)'
    )
    parser.add_argument(
        'metadata_file',
        nargs='?',
        help='Cell metadata CSV'
    )

    # Output
    parser.add_argument(
        '-o', '--output-dir',
        default='ids_results',
        help='Output directory (default: ids_results)'
    )

    # Config
    parser.add_argument(
        '-c', '--config',
        help='Configuration YAML file'
    )

    # Column names
    parser.add_argument(
        '--cell-column',
        default='cell_id',
        help='Cell ID column name (default: cell_id)'
    )
    parser.add_argument(
        '--condition-column',
        default='condition',
        help='Condition column name (default: condition)'
    )
    parser.add_argument(
        '--celltype-column',
        default='cell_type',
        help='Cell type column name (default: cell_type)'
    )
    parser.add_argument(
        '--pt-column',
        help='Pseudotime column name (optional)'
    )

    # Analysis options
    parser.add_argument(
        '--n-bins',
        type=int,
        default=25,
        help='Number of PT bins (default: 25)'
    )
    parser.add_argument(
        '--phi-threshold',
        type=float,
        default=0.5,
        help='φ threshold for onset detection (default: 0.5)'
    )

    # Utility options
    parser.add_argument(
        '--show-config',
        action='store_true',
        help='Show current configuration and exit'
    )
    parser.add_argument(
        '--generate-config',
        metavar='FILE',
        help='Generate config template and exit'
    )
    parser.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='Verbose output'
    )

    args = parser.parse_args()

    # Handle utility options
    if args.generate_config:
        print(f"Generating config template: {args.generate_config}")
        config = IDSConfig()
        config.to_yaml(args.generate_config)
        print("✓ Done")
        return

    if args.show_config:
        config = IDSConfig.from_yaml(args.config) if args.config else IDSConfig()
        config.print_summary()
        return

    # Require input files for analysis
    if not args.expression_file or not args.metadata_file:
        parser.error("Expression and metadata files are required for analysis")

    # Load config
    if args.config:
        config = IDSConfig.from_yaml(args.config)
    else:
        config = IDSConfig()

    # Override config with command-line arguments
    if args.n_bins:
        config.n_bins = args.n_bins
    if args.phi_threshold:
        config.phi_threshold = args.phi_threshold
    config.verbose = args.verbose

    # Print configuration
    if args.verbose:
        config.print_summary()

    # Initialize pipeline
    pipeline = IDSPipeline(config=config)

    # Run analysis
    try:
        results = pipeline.run(
            expression_file=args.expression_file,
            metadata_file=args.metadata_file,
            output_dir=args.output_dir,
            cell_column=args.cell_column,
            condition_column=args.condition_column,
            celltype_column=args.celltype_column,
            pt_column=args.pt_column
        )

        print(f"\n✓ Analysis complete!")
        print(f"  Results saved to: {args.output_dir}/")

    except Exception as e:
        print(f"\n✗ Error: {e}")
        if args.verbose:
            import traceback
            traceback.print_exc()
        sys.exit(1)


if __name__ == '__main__':
    main()
