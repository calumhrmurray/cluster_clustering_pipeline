#!/usr/bin/env python3
"""
Create a template configuration file for the clustering pipeline.
"""

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from src.config import create_template_config


def main():
    parser = argparse.ArgumentParser(
        description='Create a template configuration file')
    parser.add_argument('--output', '-o', type=str,
                       default='config_template.yaml',
                       help='Output path for template config')

    args = parser.parse_args()

    create_template_config(args.output)

    print(f"\nTemplate configuration created at: {args.output}")
    print("\nNext steps:")
    print("1. Edit the configuration file with your catalogue paths and parameters")
    print("2. Run: python scripts/run_analysis.py --config config.yaml --test")
    print("3. Generate SLURM jobs: python scripts/generate_jobs.py --config config.yaml")
    print("4. Submit jobs using the generated submission script")


if __name__ == '__main__':
    main()
