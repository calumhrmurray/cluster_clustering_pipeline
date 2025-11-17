#!/usr/bin/env python3
"""
Compare weighted vs unweighted clustering results.

This script creates comparison plots showing:
- Weighted and unweighted correlation functions on same axes
- Ratio plots (weighted / unweighted)
- Difference plots (weighted - unweighted)
- Summary grids showing impact across all bins
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
    """Extract bin information from filename."""
    stem = Path(filename).stem
    info = {}

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


def find_matching_files(unweighted_dir, weighted_dir):
    """
    Find pairs of weighted and unweighted result files.

    Returns
    -------
    dict
        Mapping of bin identifier to (unweighted_path, weighted_path) tuples
    """
    unweighted_files = list(Path(unweighted_dir).glob('clustering_*.pkl'))
    weighted_files = list(Path(weighted_dir).glob('clustering_*.pkl'))

    # Create mapping by filename
    pairs = {}

    for uw_file in unweighted_files:
        # Look for matching weighted file
        w_file = Path(weighted_dir) / uw_file.name

        if w_file.exists():
            info = parse_filename(uw_file.name)
            key = (info['cluster_z_min'], info['cluster_z_max'],
                   info['galaxy_z_min'], info['galaxy_z_max'])
            pairs[key] = (uw_file, w_file)

    print(f"Found {len(pairs)} matching pairs of weighted/unweighted results")
    return pairs


def plot_single_comparison(unweighted_path, weighted_path, output_path):
    """
    Create comparison plot for a single bin.

    Shows:
    - Top panel: Weighted and unweighted ω(θ) on same axes
    - Bottom panel: Ratio (weighted/unweighted)
    """
    # Load results
    uw_data = load_result(unweighted_path)
    w_data = load_result(weighted_path)

    # Extract data
    uw_results = uw_data.get('results', uw_data)
    w_results = w_data.get('results', w_data)

    r_uw = uw_results.get('r', uw_results.get('meanr'))
    xi_uw = uw_results.get('xi')
    sig_uw = uw_results.get('sigma_xi', np.sqrt(uw_results.get('var_xi', 0)))

    r_w = w_results.get('r', w_results.get('meanr'))
    xi_w = w_results.get('xi')
    sig_w = w_results.get('sigma_xi', np.sqrt(w_results.get('var_xi', 0)))

    # Extract bin info
    info = parse_filename(unweighted_path.name)
    cz_min, cz_max = info['cluster_z_min'], info['cluster_z_max']
    gz_min, gz_max = info['galaxy_z_min'], info['galaxy_z_max']

    # Create figure with 2 subplots
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10), sharex=True,
                                    gridspec_kw={'height_ratios': [2, 1]})

    # Top panel: Both correlation functions
    ax1.errorbar(r_uw, xi_uw, yerr=sig_uw, marker='o', ms=6,
                label='Unweighted', color='C0', alpha=0.8, capsize=3)
    ax1.errorbar(r_w, xi_w, yerr=sig_w, marker='s', ms=6,
                label='Weighted', color='C1', alpha=0.8, capsize=3)

    ax1.axhline(0, color='k', linestyle='--', alpha=0.3, lw=1)
    ax1.set_ylabel(r'$\omega(\theta)$', fontsize=14)
    ax1.legend(loc='best', fontsize=12, framealpha=0.9)
    ax1.grid(True, alpha=0.3)
    ax1.set_xscale('log')

    title = (f'Cluster-Galaxy Clustering: Weighted vs Unweighted\n'
            f'Clusters: {cz_min:.2f} ≤ z < {cz_max:.2f}, '
            f'Galaxies: {gz_min:.2f} ≤ z < {gz_max:.2f}')
    ax1.set_title(title, fontsize=13, pad=10)

    # Bottom panel: Ratio
    # Only plot where both are non-zero and positive
    valid = (xi_uw != 0) & (xi_w != 0) & (xi_uw > 0) & (xi_w > 0)

    if np.any(valid):
        ratio = xi_w[valid] / xi_uw[valid]

        # Propagate errors for ratio
        if sig_uw is not None and sig_w is not None:
            # σ(w/uw) = (w/uw) * sqrt((σ_w/w)² + (σ_uw/uw)²)
            rel_err_uw = sig_uw[valid] / np.abs(xi_uw[valid])
            rel_err_w = sig_w[valid] / np.abs(xi_w[valid])
            ratio_err = np.abs(ratio) * np.sqrt(rel_err_uw**2 + rel_err_w**2)
        else:
            ratio_err = None

        ax2.errorbar(r_w[valid], ratio, yerr=ratio_err, marker='o', ms=6,
                    color='C2', alpha=0.8, capsize=3)

        ax2.axhline(1, color='k', linestyle='--', alpha=0.5, lw=1.5,
                   label='No change')
        ax2.set_xlabel('Angular separation [arcmin]', fontsize=14)
        ax2.set_ylabel('Weighted / Unweighted', fontsize=12)
        ax2.legend(loc='best', fontsize=10)
        ax2.grid(True, alpha=0.3)
        ax2.set_xscale('log')
        ax2.set_ylim([0.5, 1.5])  # Reasonable range for ratio

    plt.tight_layout()

    # Save
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    print(f"Saved comparison to {output_path}")
    plt.close()


def plot_all_comparisons(unweighted_dir, weighted_dir, output_dir):
    """Create individual comparison plots for all matching bins."""
    pairs = find_matching_files(unweighted_dir, weighted_dir)

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    for key, (uw_path, w_path) in sorted(pairs.items()):
        cz_min, cz_max, gz_min, gz_max = key
        output_filename = f"comparison_cl_{cz_min:.2f}_{cz_max:.2f}_gal_{gz_min:.2f}_{gz_max:.2f}.png"
        output_path = output_dir / output_filename

        plot_single_comparison(uw_path, w_path, output_path)


def plot_summary_comparison_grid(unweighted_dir, weighted_dir, output_path):
    """
    Create summary grid showing weighted vs unweighted for all bins.

    Each panel shows one cluster z-bin with all galaxy z-bins overlaid.
    Solid lines = unweighted, dashed lines = weighted.
    """
    pairs = find_matching_files(unweighted_dir, weighted_dir)

    # Group by cluster redshift
    cluster_bins = defaultdict(list)
    for key, (uw_path, w_path) in pairs.items():
        cz_min, cz_max, gz_min, gz_max = key
        cluster_bins[(cz_min, cz_max)].append({
            'galaxy_z': (gz_min, gz_max),
            'unweighted': uw_path,
            'weighted': w_path
        })

    # Sort cluster bins
    cluster_keys = sorted(cluster_bins.keys())
    n_cluster_bins = len(cluster_keys)

    # Create figure
    fig, axes = plt.subplots(n_cluster_bins, 1, figsize=(12, 4 * n_cluster_bins),
                            sharex=True)
    if n_cluster_bins == 1:
        axes = [axes]

    # Colors for galaxy bins
    galaxy_z_bins = sorted(set(item['galaxy_z'] for items in cluster_bins.values() for item in items))
    colors = plt.cm.viridis(np.linspace(0, 0.9, len(galaxy_z_bins)))
    color_map = {gz: colors[i] for i, gz in enumerate(galaxy_z_bins)}

    # Plot each cluster bin
    for idx, cz_key in enumerate(cluster_keys):
        ax = axes[idx]
        items = sorted(cluster_bins[cz_key], key=lambda x: x['galaxy_z'])

        for item in items:
            gz_min, gz_max = item['galaxy_z']
            color = color_map[item['galaxy_z']]

            # Load unweighted
            uw_data = load_result(item['unweighted'])
            uw_results = uw_data.get('results', uw_data)
            r_uw = uw_results.get('r', uw_results.get('meanr'))
            xi_uw = uw_results.get('xi')

            # Load weighted
            w_data = load_result(item['weighted'])
            w_results = w_data.get('results', w_data)
            r_w = w_results.get('r', w_results.get('meanr'))
            xi_w = w_results.get('xi')

            # Plot
            label = f'{gz_min:.1f}-{gz_max:.1f}'
            ax.plot(r_uw, xi_uw, marker='o', ms=4, linestyle='-', color=color,
                   alpha=0.7, label=f'{label} (unweighted)')
            ax.plot(r_w, xi_w, marker='s', ms=4, linestyle='--', color=color,
                   alpha=0.7, label=f'{label} (weighted)')

        # Formatting
        cz_min, cz_max = cz_key
        ax.axhline(0, color='k', linestyle=':', alpha=0.3, lw=1)
        ax.set_ylabel(r'$\omega(\theta)$', fontsize=11)
        ax.set_xscale('log')
        ax.grid(True, alpha=0.3)
        ax.legend(loc='best', fontsize=8, ncol=2)
        ax.set_title(f'Clusters: {cz_min:.2f} ≤ z < {cz_max:.2f}', fontsize=11)

    axes[-1].set_xlabel('Angular separation [arcmin]', fontsize=12)

    fig.suptitle('Clustering: Weighted (dashed) vs Unweighted (solid)',
                fontsize=14, y=0.995)
    plt.tight_layout()

    # Save
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    print(f"\nSaved summary grid to {output_path}")
    plt.close()


def plot_ratio_heatmap(unweighted_dir, weighted_dir, output_path, scale_index=5):
    """
    Create 2D heatmap showing weighted/unweighted ratio across all bins.

    Parameters
    ----------
    scale_index : int
        Which angular scale bin to show (default: 5, typically ~1 arcmin)
    """
    pairs = find_matching_files(unweighted_dir, weighted_dir)

    # Get all cluster and galaxy z bins
    cluster_z_bins = sorted(set((key[0], key[1]) for key in pairs.keys()))
    galaxy_z_bins = sorted(set((key[2], key[3]) for key in pairs.keys()))

    # Create ratio matrix
    ratio_matrix = np.full((len(cluster_z_bins), len(galaxy_z_bins)), np.nan)

    for key, (uw_path, w_path) in pairs.items():
        cz_min, cz_max, gz_min, gz_max = key

        # Find indices
        cz_idx = cluster_z_bins.index((cz_min, cz_max))
        gz_idx = galaxy_z_bins.index((gz_min, gz_max))

        # Load results
        uw_data = load_result(uw_path)
        w_data = load_result(w_path)

        uw_results = uw_data.get('results', uw_data)
        w_results = w_data.get('results', w_data)

        xi_uw = uw_results.get('xi')
        xi_w = w_results.get('xi')

        # Compute ratio at specified scale
        if len(xi_uw) > scale_index and len(xi_w) > scale_index:
            if xi_uw[scale_index] != 0 and xi_uw[scale_index] > 0:
                ratio = xi_w[scale_index] / xi_uw[scale_index]
                ratio_matrix[cz_idx, gz_idx] = ratio

    # Create figure
    fig, ax = plt.subplots(figsize=(10, 8))

    # Plot heatmap
    im = ax.imshow(ratio_matrix, cmap='RdBu_r', vmin=0.7, vmax=1.3,
                  aspect='auto', origin='lower')

    # Add colorbar
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label('Weighted / Unweighted', fontsize=12)

    # Set ticks
    ax.set_xticks(range(len(galaxy_z_bins)))
    ax.set_yticks(range(len(cluster_z_bins)))

    ax.set_xticklabels([f'{gz[0]:.1f}-{gz[1]:.1f}' for gz in galaxy_z_bins], rotation=45)
    ax.set_yticklabels([f'{cz[0]:.1f}-{cz[1]:.1f}' for cz in cluster_z_bins])

    ax.set_xlabel('Galaxy Redshift Bin', fontsize=12)
    ax.set_ylabel('Cluster Redshift Bin', fontsize=12)
    ax.set_title(f'Impact of Systematic Weights (ratio at scale bin {scale_index})',
                fontsize=13, pad=15)

    # Add text annotations
    for i in range(len(cluster_z_bins)):
        for j in range(len(galaxy_z_bins)):
            if not np.isnan(ratio_matrix[i, j]):
                text_color = 'white' if abs(ratio_matrix[i, j] - 1.0) > 0.15 else 'black'
                ax.text(j, i, f'{ratio_matrix[i, j]:.2f}',
                       ha='center', va='center', fontsize=9, color=text_color)

    plt.tight_layout()

    # Save
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    print(f"Saved ratio heatmap to {output_path}")
    plt.close()


def main():
    parser = argparse.ArgumentParser(
        description='Compare weighted vs unweighted clustering results'
    )
    parser.add_argument('--unweighted-dir', type=str,
                       default='outputs_rr2_angular',
                       help='Directory with unweighted results')
    parser.add_argument('--weighted-dir', type=str,
                       default='outputs_rr2_angular_weighted',
                       help='Directory with weighted results')
    parser.add_argument('--output-dir', type=str,
                       default='outputs_rr2_angular_weighted/figures/comparison',
                       help='Output directory for comparison plots')
    parser.add_argument('--plot-type', type=str, default='all',
                       choices=['all', 'individual', 'summary', 'heatmap'],
                       help='Type of plots to create')

    args = parser.parse_args()

    unweighted_dir = Path(args.unweighted_dir)
    weighted_dir = Path(args.weighted_dir)
    output_dir = Path(args.output_dir)

    # Check directories exist
    if not unweighted_dir.exists():
        print(f"Error: Unweighted directory not found: {unweighted_dir}")
        return

    if not weighted_dir.exists():
        print(f"Error: Weighted directory not found: {weighted_dir}")
        return

    # Count results
    n_unweighted = len(list(unweighted_dir.glob('clustering_*.pkl')))
    n_weighted = len(list(weighted_dir.glob('clustering_*.pkl')))

    print(f"Found {n_unweighted} unweighted results")
    print(f"Found {n_weighted} weighted results")

    if n_unweighted == 0 or n_weighted == 0:
        print("Need both weighted and unweighted results!")
        return

    # Create plots
    if args.plot_type in ['all', 'individual']:
        print("\n" + "="*70)
        print("Creating individual comparison plots...")
        print("="*70)
        plot_all_comparisons(unweighted_dir, weighted_dir,
                           output_dir / 'individual')

    if args.plot_type in ['all', 'summary']:
        print("\n" + "="*70)
        print("Creating summary comparison grid...")
        print("="*70)
        plot_summary_comparison_grid(unweighted_dir, weighted_dir,
                                     output_dir / 'comparison_summary_grid.png')

    if args.plot_type in ['all', 'heatmap']:
        print("\n" + "="*70)
        print("Creating ratio heatmap...")
        print("="*70)
        plot_ratio_heatmap(unweighted_dir, weighted_dir,
                          output_dir / 'ratio_heatmap.png')

    print("\n" + "="*70)
    print("COMPARISON PLOTS COMPLETE")
    print("="*70)
    print(f"Figures saved to {output_dir}")


if __name__ == '__main__':
    main()
