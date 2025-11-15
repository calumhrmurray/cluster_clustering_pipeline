#!/usr/bin/env python3
"""
Plot clustering results from output files.
"""

import argparse
import sys
from pathlib import Path
import pickle
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent))


def load_result(result_file):
    """Load a result file."""
    with open(result_file, 'rb') as f:
        return pickle.load(f)


def plot_clustering_results(result_files, output_path=None, log_scale=True):
    """
    Plot clustering results from multiple files.

    Parameters
    ----------
    result_files : list of Path
        List of result file paths
    output_path : Path, optional
        Where to save the plot
    log_scale : bool
        Use log scale for axes
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

    colors = plt.cm.viridis(np.linspace(0, 1, len(result_files)))

    for i, result_file in enumerate(result_files):
        data = load_result(result_file)
        results = data['results']
        metadata = data.get('metadata', {})

        # Create label from metadata
        if 'cluster_filters' in metadata and 'galaxy_filters' in metadata:
            c_filters = metadata['cluster_filters']
            g_filters = metadata['galaxy_filters']

            # Extract redshift info for label
            label_parts = []
            for col, val in c_filters.items():
                if 'redshift' in col.lower() or 'z_cluster' in col.lower():
                    if isinstance(val, tuple):
                        label_parts.append(f"z_c={val[0]:.1f}-{val[1]:.1f}")
            for col, val in g_filters.items():
                if 'redshift' in col.lower():
                    if isinstance(val, tuple):
                        label_parts.append(f"z_g={val[0]:.1f}-{val[1]:.1f}")

            label = ", ".join(label_parts) if label_parts else result_file.stem
        else:
            label = result_file.stem

        # Plot correlation function
        r = results['r']
        xi = results['xi']
        sigma_xi = results['sigma_xi']

        # Filter out zero or negative values for log plot
        valid = (results['npairs'] > 0) & np.isfinite(xi)

        ax1.errorbar(r[valid], xi[valid], yerr=sigma_xi[valid],
                    fmt='o-', label=label, color=colors[i], capsize=3)

        # Plot S/N
        sn = np.abs(xi) / sigma_xi
        ax2.plot(r[valid], sn[valid], 'o-', label=label, color=colors[i])

    # Format correlation function plot
    ax1.set_xlabel(r'$r_{\perp}$ [Mpc/h]', fontsize=12)
    ax1.set_ylabel(r'$w(\theta)$', fontsize=12)
    ax1.legend(fontsize=10)
    ax1.grid(True, alpha=0.3)

    if log_scale:
        ax1.set_xscale('log')
        ax1.set_yscale('log')
        ax1.set_title('Correlation Function (log-log)')
    else:
        ax1.set_title('Correlation Function')

    # Format S/N plot
    ax2.set_xlabel(r'$r_{\perp}$ [Mpc/h]', fontsize=12)
    ax2.set_ylabel('Signal-to-Noise', fontsize=12)
    ax2.axhline(y=1, color='k', linestyle='--', alpha=0.3)
    ax2.legend(fontsize=10)
    ax2.grid(True, alpha=0.3)

    if log_scale:
        ax2.set_xscale('log')
        ax2.set_yscale('log')
        ax2.set_title('Signal-to-Noise (log-log)')
    else:
        ax2.set_title('Signal-to-Noise')

    plt.tight_layout()

    if output_path:
        plt.savefig(output_path, dpi=150, bbox_inches='tight')
        print(f"Plot saved to {output_path}")
    else:
        plt.show()


def plot_comparison(result_files, labels=None, output_path=None):
    """
    Plot comparison of different clustering measurements.

    Parameters
    ----------
    result_files : list of Path
        Result files to compare
    labels : list of str, optional
        Custom labels for each file
    output_path : Path, optional
        Where to save plot
    """
    fig, ax = plt.subplots(figsize=(10, 7))

    colors = plt.cm.tab10(np.arange(len(result_files)))

    for i, result_file in enumerate(result_files):
        data = load_result(result_file)
        results = data['results']

        label = labels[i] if labels else result_file.stem

        r = results['r']
        xi = results['xi']
        sigma_xi = results['sigma_xi']
        valid = (results['npairs'] > 0) & np.isfinite(xi)

        ax.errorbar(r[valid], xi[valid], yerr=sigma_xi[valid],
                   fmt='o-', label=label, color=colors[i],
                   capsize=3, markersize=6)

    ax.set_xlabel(r'$r_{\perp}$ [Mpc/h]', fontsize=14)
    ax.set_ylabel(r'$w(\theta)$', fontsize=14)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)
    ax.set_title('Clustering Comparison', fontsize=15)

    plt.tight_layout()

    if output_path:
        plt.savefig(output_path, dpi=150, bbox_inches='tight')
        print(f"Comparison plot saved to {output_path}")
    else:
        plt.show()


def main():
    parser = argparse.ArgumentParser(
        description='Plot clustering results')
    parser.add_argument('results', nargs='+', type=str,
                       help='Result file(s) to plot')
    parser.add_argument('--output', '-o', type=str,
                       help='Output path for plot')
    parser.add_argument('--labels', nargs='+', type=str,
                       help='Custom labels for plots')
    parser.add_argument('--linear', action='store_true',
                       help='Use linear scale instead of log')
    parser.add_argument('--compare', action='store_true',
                       help='Create comparison plot')

    args = parser.parse_args()

    # Convert paths
    result_files = [Path(f) for f in args.results]

    # Check files exist
    for f in result_files:
        if not f.exists():
            print(f"Error: File not found: {f}")
            return 1

    output_path = Path(args.output) if args.output else None

    if args.compare:
        plot_comparison(result_files, labels=args.labels,
                       output_path=output_path)
    else:
        plot_clustering_results(result_files, output_path=output_path,
                               log_scale=not args.linear)

    return 0


if __name__ == '__main__':
    sys.exit(main())
