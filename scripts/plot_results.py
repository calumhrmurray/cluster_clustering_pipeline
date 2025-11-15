#!/usr/bin/env python3
"""
Plot clustering results from the angular correlation analysis.

This script creates various diagnostic and publication-quality plots:
- Correlation functions vs separation for different bins
- Tomographic redshift comparison
- Cluster richness/redshift dependencies
"""

import argparse
import pickle
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import re
from collections import defaultdict


def load_result(filepath):
    """Load a single result file."""
    with open(filepath, 'rb') as f:
        return pickle.load(f)


def parse_filename(filename):
    """
    Parse result filename to extract bin information.

    Example: clustering_c_RICHNESS_CLUSTER_40.00_80.00_Z_CLUSTER_0.20_0.38_g_redshift_0.00_0.30.pkl

    Returns
    -------
    dict
        Dictionary with cluster and galaxy bin information
    """
    stem = Path(filename).stem

    info = {}

    # Extract cluster richness
    rich_match = re.search(r'RICHNESS_CLUSTER_([\d.]+)_([\d.]+)', stem)
    if rich_match:
        info['rich_min'] = float(rich_match.group(1))
        info['rich_max'] = float(rich_match.group(2))

    # Extract cluster redshift
    cz_match = re.search(r'Z_CLUSTER_([\d.]+)_([\d.]+)', stem)
    if cz_match:
        info['cluster_z_min'] = float(cz_match.group(1))
        info['cluster_z_max'] = float(cz_match.group(2))

    # Extract galaxy redshift
    gz_match = re.search(r'redshift_([\d.]+)_([\d.]+)', stem)
    if gz_match:
        info['galaxy_z_min'] = float(gz_match.group(1))
        info['galaxy_z_max'] = float(gz_match.group(2))

    return info


def plot_tomographic_comparison(results_dir, output_path, cluster_bin=None):
    """
    Plot correlation functions for all galaxy redshift bins for a given cluster bin.

    Parameters
    ----------
    results_dir : Path
        Directory containing result files
    output_path : Path
        Output path for figure
    cluster_bin : dict, optional
        Cluster bin to plot (richness and redshift ranges)
        If None, uses the first available cluster bin
    """
    # Load all results
    result_files = sorted(Path(results_dir).glob('clustering_*.pkl'))

    if len(result_files) == 0:
        print(f"No result files found in {results_dir}")
        return

    # Group by cluster bin
    cluster_bins = defaultdict(list)

    for filepath in result_files:
        info = parse_filename(filepath.name)
        result = load_result(filepath)

        # Create cluster bin key
        cluster_key = (info.get('rich_min'), info.get('rich_max'),
                      info.get('cluster_z_min'), info.get('cluster_z_max'))

        cluster_bins[cluster_key].append({
            'info': info,
            'result': result,
            'filepath': filepath
        })

    # Select cluster bin to plot
    if cluster_bin is None:
        cluster_key = sorted(cluster_bins.keys())[0]
    else:
        cluster_key = (cluster_bin['rich_min'], cluster_bin['rich_max'],
                      cluster_bin['cluster_z_min'], cluster_bin['cluster_z_max'])

    if cluster_key not in cluster_bins:
        print(f"Cluster bin {cluster_key} not found")
        return

    results_list = cluster_bins[cluster_key]

    # Sort by galaxy redshift
    results_list.sort(key=lambda x: x['info']['galaxy_z_min'])

    # Create figure
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(11, 12), sharex=True)

    # Color palette
    colors = plt.cm.viridis(np.linspace(0, 0.9, len(results_list)))

    for i, item in enumerate(results_list):
        info = item['info']
        result = item['result']

        # Get correlation function data
        results_data = result.get('results', result)  # Handle both formats
        meanr = results_data.get('r', results_data.get('meanr', results_data.get('rnom')))
        xi = results_data.get('xi', results_data.get('correlation'))

        # Calculate error (from var_xi or varxi if available)
        varxi = results_data.get('var_xi', results_data.get('varxi'))
        if varxi is not None and not np.all(varxi == 0):
            xi_err = np.sqrt(varxi)
        else:
            xi_err = results_data.get('sigma_xi')  # Try direct sigma if available

        # Galaxy redshift bin label
        gz_min = info['galaxy_z_min']
        gz_max = info['galaxy_z_max']
        label = f'{gz_min:.1f} ≤ z_gal < {gz_max:.1f}'

        # Plot correlation function
        ax1.errorbar(meanr, xi, yerr=xi_err, marker='o', ms=6,
                    label=label, color=colors[i], alpha=0.8, capsize=3)

        # Plot log version
        valid = xi > 0
        if np.any(valid):
            ax2.errorbar(meanr[valid], xi[valid], yerr=xi_err[valid] if xi_err is not None else None,
                        marker='o', ms=6, label=label, color=colors[i], alpha=0.8, capsize=3)

    # Formatting
    rich_min, rich_max, cz_min, cz_max = cluster_key
    if rich_min is not None and rich_max is not None:
        title = (f'Cluster-Galaxy Clustering (Tomographic)\n'
                f'Clusters: {rich_min:.0f} ≤ N_rich < {rich_max:.0f}, '
                f'{cz_min:.2f} ≤ z_cl < {cz_max:.2f}')
    else:
        title = (f'Cluster-Galaxy Clustering (Tomographic)\n'
                f'Clusters: {cz_min:.2f} ≤ z_cl < {cz_max:.2f}')

    ax1.set_ylabel(r'$\omega(\theta)$', fontsize=14)
    ax1.axhline(0, color='k', linestyle='--', alpha=0.3)
    ax1.legend(loc='upper right', fontsize=11, framealpha=0.9)
    ax1.grid(True, alpha=0.3)
    ax1.set_title(title, fontsize=14, pad=15)

    ax2.set_xlabel('Angular separation [arcmin]', fontsize=14)
    ax2.set_ylabel(r'$\omega(\theta)$', fontsize=14)
    ax2.set_xscale('log')
    ax2.set_yscale('log')
    ax2.legend(loc='lower left', fontsize=11, framealpha=0.9)
    ax2.grid(True, alpha=0.3, which='both')

    plt.tight_layout()

    # Save
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    print(f"Saved tomographic plot to {output_path}")
    plt.close()


def plot_all_cluster_bins(results_dir, output_dir):
    """
    Create tomographic plots for each cluster bin.

    Parameters
    ----------
    results_dir : Path
        Directory containing result files
    output_dir : Path
        Output directory for figures
    """
    # Load all results
    result_files = sorted(Path(results_dir).glob('clustering_*.pkl'))

    if len(result_files) == 0:
        print(f"No result files found in {results_dir}")
        return

    # Group by cluster bin
    cluster_bins = defaultdict(list)

    for filepath in result_files:
        info = parse_filename(filepath.name)
        result = load_result(filepath)

        # Create cluster bin key
        cluster_key = (info.get('rich_min'), info.get('rich_max'),
                      info.get('cluster_z_min'), info.get('cluster_z_max'))

        cluster_bins[cluster_key].append({
            'info': info,
            'result': result,
            'filepath': filepath
        })

    print(f"\nFound {len(cluster_bins)} unique cluster bins")

    # Create plot for each cluster bin
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    for cluster_key in sorted(cluster_bins.keys()):
        rich_min, rich_max, cz_min, cz_max = cluster_key

        # Handle case where richness bins don't exist
        if rich_min is not None and rich_max is not None:
            output_filename = (f"tomographic_rich{rich_min:.0f}_{rich_max:.0f}_"
                              f"z{cz_min:.2f}_{cz_max:.2f}.png")
        else:
            output_filename = f"tomographic_z{cz_min:.2f}_{cz_max:.2f}.png"

        output_path = output_dir / output_filename

        cluster_bin = {
            'rich_min': rich_min,
            'rich_max': rich_max,
            'cluster_z_min': cz_min,
            'cluster_z_max': cz_max
        }

        plot_tomographic_comparison(results_dir, output_path, cluster_bin)


def plot_summary_grid(results_dir, output_path):
    """
    Create a summary grid showing all bins.

    Parameters
    ----------
    results_dir : Path
        Directory containing result files
    output_path : Path
        Output path for figure
    """
    # Load all results
    result_files = sorted(Path(results_dir).glob('clustering_*.pkl'))

    if len(result_files) == 0:
        print(f"No result files found in {results_dir}")
        return

    # Parse all results
    all_results = []
    for filepath in result_files:
        info = parse_filename(filepath.name)
        result = load_result(filepath)
        all_results.append({'info': info, 'result': result})

    # Get unique bins
    cluster_z_bins = sorted(set((r['info']['cluster_z_min'], r['info']['cluster_z_max'])
                                for r in all_results))

    # Check if richness bins exist
    if 'rich_min' in all_results[0]['info']:
        rich_bins = sorted(set((r['info']['rich_min'], r['info']['rich_max'])
                               for r in all_results))
        n_rich = len(rich_bins)
    else:
        rich_bins = [(None, None)]  # Dummy richness bin
        n_rich = 1

    galaxy_z_bins = sorted(set((r['info']['galaxy_z_min'], r['info']['galaxy_z_max'])
                               for r in all_results))

    n_cluster_z = len(cluster_z_bins)
    n_galaxy_z = len(galaxy_z_bins)

    print(f"\nCreating summary grid:")
    print(f"  Cluster redshift bins: {n_cluster_z}")
    print(f"  Richness bins: {n_rich}")
    print(f"  Galaxy redshift bins: {n_galaxy_z}")
    print(f"  Total combinations: {n_cluster_z * n_rich * n_galaxy_z}")

    # Create grid (cluster bins × galaxy bins)
    n_cluster_bins = n_cluster_z * n_rich
    fig, axes = plt.subplots(n_cluster_bins, 1, figsize=(10, 4 * n_cluster_bins),
                            sharex=True)

    if n_cluster_bins == 1:
        axes = [axes]

    # Colors for galaxy bins
    colors = plt.cm.viridis(np.linspace(0, 0.9, n_galaxy_z))

    # Plot each cluster bin
    idx = 0
    for cz_bin in cluster_z_bins:
        for rich_bin in rich_bins:
            ax = axes[idx]

            # Find results for this cluster bin
            for gal_idx, gz_bin in enumerate(galaxy_z_bins):
                # Find matching result
                if rich_bin[0] is not None:
                    # With richness bins
                    matches = [r for r in all_results
                              if (r['info']['cluster_z_min'], r['info']['cluster_z_max']) == cz_bin
                              and (r['info'].get('rich_min'), r['info'].get('rich_max')) == rich_bin
                              and (r['info']['galaxy_z_min'], r['info']['galaxy_z_max']) == gz_bin]
                else:
                    # Without richness bins
                    matches = [r for r in all_results
                              if (r['info']['cluster_z_min'], r['info']['cluster_z_max']) == cz_bin
                              and (r['info']['galaxy_z_min'], r['info']['galaxy_z_max']) == gz_bin]

                if len(matches) > 0:
                    result_data = matches[0]['result']
                    results = result_data.get('results', result_data)  # Handle both formats
                    meanr = results.get('r', results.get('meanr', results.get('rnom')))
                    xi = results.get('xi', results.get('correlation'))

                    varxi = results.get('var_xi', results.get('varxi'))
                    if varxi is not None and not np.all(varxi == 0):
                        xi_err = np.sqrt(varxi)
                    else:
                        xi_err = results.get('sigma_xi')

                    label = f'{gz_bin[0]:.1f}-{gz_bin[1]:.1f}'
                    ax.errorbar(meanr, xi, yerr=xi_err, marker='o', ms=4,
                               label=label, color=colors[gal_idx], alpha=0.7)

            # Formatting
            ax.axhline(0, color='k', linestyle='--', alpha=0.3, lw=0.5)
            ax.set_ylabel(r'$\omega(\theta)$', fontsize=10)
            ax.legend(loc='best', fontsize=8, ncol=5, title=r'$z_{gal}$')
            ax.grid(True, alpha=0.3)

            # Handle case where richness bins don't exist
            if rich_bin[0] is not None and rich_bin[1] is not None:
                title = f'Rich: {rich_bin[0]:.0f}-{rich_bin[1]:.0f}, Z_cl: {cz_bin[0]:.2f}-{cz_bin[1]:.2f}'
            else:
                title = f'Z_cl: {cz_bin[0]:.2f}-{cz_bin[1]:.2f}'
            ax.set_title(title, fontsize=10)

            idx += 1

    axes[-1].set_xlabel('Angular separation [arcmin]', fontsize=12)
    axes[-1].set_xscale('log')

    plt.tight_layout()

    # Save
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    print(f"\nSaved summary grid to {output_path}")
    plt.close()


def main():
    parser = argparse.ArgumentParser(
        description='Plot clustering results from angular correlation analysis'
    )
    parser.add_argument('--results-dir', type=str,
                       default='outputs_rr2_angular',
                       help='Directory containing result .pkl files')
    parser.add_argument('--output-dir', type=str,
                       default='outputs_rr2_angular/figures',
                       help='Output directory for figures')
    parser.add_argument('--plot-type', type=str, default='all',
                       choices=['all', 'tomographic', 'summary'],
                       help='Type of plots to create')

    args = parser.parse_args()

    results_dir = Path(args.results_dir)
    output_dir = Path(args.output_dir)

    if not results_dir.exists():
        print(f"Error: Results directory not found: {results_dir}")
        return

    # Count available results
    n_results = len(list(results_dir.glob('clustering_*.pkl')))
    print(f"Found {n_results} result files in {results_dir}")

    if n_results == 0:
        print("No results to plot!")
        return

    # Create plots
    if args.plot_type in ['all', 'tomographic']:
        print("\n" + "="*70)
        print("Creating tomographic comparison plots for each cluster bin...")
        print("="*70)
        plot_all_cluster_bins(results_dir, output_dir / 'tomographic')

    if args.plot_type in ['all', 'summary']:
        print("\n" + "="*70)
        print("Creating summary grid plot...")
        print("="*70)
        plot_summary_grid(results_dir, output_dir / 'clustering_summary_grid.png')

    print("\n" + "="*70)
    print("PLOTTING COMPLETE")
    print("="*70)
    print(f"Figures saved to {output_dir}")


if __name__ == '__main__':
    main()
