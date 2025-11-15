#!/usr/bin/env python3
"""
Create n(z) diagnostic plots showing galaxy redshift distributions in tomographic bins.
"""

import argparse
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.table import Table
from pathlib import Path
import yaml


def load_galaxy_catalogue(filepath, z_col='phz_median'):
    """
    Load galaxy catalogue and extract redshifts.

    Parameters
    ----------
    filepath : str
        Path to galaxy catalogue FITS file
    z_col : str
        Name of redshift column

    Returns
    -------
    array
        Array of redshifts
    """
    print(f"Loading galaxy catalogue from {filepath}")
    galaxies = Table(fits.open(filepath)[1].data)

    # Get redshifts and remove NaNs
    redshifts = galaxies[z_col]
    valid = ~np.isnan(redshifts)
    redshifts = redshifts[valid]

    print(f"Loaded {len(redshifts)} galaxies with valid redshifts")
    return redshifts


def create_nz_plot(redshifts, bin_edges, output_path, normalized=True):
    """
    Create n(z) plot with all tomographic bins on same plot.

    Parameters
    ----------
    redshifts : array
        All galaxy redshifts
    bin_edges : list
        Redshift bin edges for tomographic bins
    output_path : str or Path
        Output path for figure
    normalized : bool
        Whether to normalize histograms
    """
    fig, ax = plt.subplots(figsize=(11, 9))

    # Color palette for different bins
    colors = plt.cm.tab10(np.linspace(0, 1, len(bin_edges) - 1))

    # Create histogram bins for plotting (higher resolution)
    z_bins = np.linspace(0, max(bin_edges), 100)

    # Plot each tomographic bin
    for i in range(len(bin_edges) - 1):
        z_min, z_max = bin_edges[i], bin_edges[i + 1]

        # Select galaxies in this bin
        in_bin = (redshifts >= z_min) & (redshifts < z_max)
        z_in_bin = redshifts[in_bin]

        # Create histogram
        counts, edges = np.histogram(z_in_bin, bins=z_bins)

        # Normalize if requested
        if normalized:
            # Normalize to peak = 1 for each bin
            if counts.max() > 0:
                counts = counts / counts.max()

        # Plot as histogram
        bin_centers = 0.5 * (edges[1:] + edges[:-1])
        ax.bar(bin_centers, counts, width=np.diff(edges),
               alpha=0.7, color=colors[i],
               label=f'{z_min:.1f} ≤ z < {z_max:.1f} (N={len(z_in_bin):,})',
               align='center')

    # Formatting
    ax.set_xlabel('Redshift', fontsize=14)
    if normalized:
        ax.set_ylabel('Normalised counts', fontsize=14)
    else:
        ax.set_ylabel('Counts', fontsize=14)

    ax.set_xlim(0, max(bin_edges) * 1.1)
    ax.set_ylim(0, None)

    # Add grid
    ax.grid(True, alpha=0.3, linestyle='--')

    # Legend
    ax.legend(loc='upper right', fontsize=12, framealpha=0.9)

    # Title
    title = f'Galaxy Redshift Distribution - Tomographic Bins'
    ax.set_title(title, fontsize=16, pad=15)

    # Tight layout
    plt.tight_layout()

    # Save figure
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    print(f"\nFigure saved to {output_path}")

    plt.close()


def main():
    parser = argparse.ArgumentParser(
        description='Create n(z) diagnostic plots for tomographic galaxy bins'
    )
    parser.add_argument('--catalogue', type=str, required=True,
                       help='Path to galaxy catalogue FITS file')
    parser.add_argument('--config', type=str,
                       help='Path to config file to read bin edges from')
    parser.add_argument('--bins', type=float, nargs='+',
                       help='Redshift bin edges (alternative to --config)')
    parser.add_argument('--z-col', type=str, default='phz_median',
                       help='Name of redshift column (default: phz_median)')
    parser.add_argument('--output', type=str,
                       default='outputs/figures/galaxy_nz_all_bins.png',
                       help='Output path for figure')
    parser.add_argument('--no-normalize', action='store_true',
                       help='Do not normalize histograms')

    args = parser.parse_args()

    # Get bin edges
    if args.config:
        print(f"Reading bin edges from config: {args.config}")
        with open(args.config, 'r') as f:
            config = yaml.safe_load(f)

        galaxy_bins = config.get('galaxy_bins', {})
        redshift_config = galaxy_bins.get('redshift', {})
        bin_edges = redshift_config.get('bins')

        if not bin_edges:
            raise ValueError("No redshift bins found in config")

    elif args.bins:
        bin_edges = args.bins
    else:
        raise ValueError("Must provide either --config or --bins")

    print(f"\nTomographic bin edges: {bin_edges}")
    print(f"Number of bins: {len(bin_edges) - 1}")

    # Load galaxy redshifts
    redshifts = load_galaxy_catalogue(args.catalogue, z_col=args.z_col)

    # Create plot
    create_nz_plot(redshifts, bin_edges, args.output,
                   normalized=not args.no_normalize)

    # Print summary statistics
    print("\n" + "="*70)
    print("SUMMARY STATISTICS")
    print("="*70)

    for i in range(len(bin_edges) - 1):
        z_min, z_max = bin_edges[i], bin_edges[i + 1]
        in_bin = (redshifts >= z_min) & (redshifts < z_max)
        n_in_bin = np.sum(in_bin)
        fraction = 100 * n_in_bin / len(redshifts)

        print(f"Bin {i+1}: {z_min:.1f} ≤ z < {z_max:.1f}")
        print(f"  N = {n_in_bin:,} ({fraction:.1f}%)")

        if n_in_bin > 0:
            z_in_bin = redshifts[in_bin]
            print(f"  Mean z = {np.mean(z_in_bin):.3f}")
            print(f"  Median z = {np.median(z_in_bin):.3f}")
            print(f"  Std z = {np.std(z_in_bin):.3f}")
        print()


if __name__ == '__main__':
    main()
