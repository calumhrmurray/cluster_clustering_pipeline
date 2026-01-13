"""
Plotting utilities for cluster-galaxy clustering analysis.
Provides helper functions for creating publication-quality plots.
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import pickle


def setup_plot_style(font_size_minor=15, font_size_major=20):
    """
    Set up matplotlib style for publication-quality plots.

    Parameters
    ----------
    font_size_minor : int
        Font size for minor labels (tick labels, etc.)
    font_size_major : int
        Font size for major labels (axis labels, titles)
    """
    plt.rcParams['font.family'] = 'serif'
    plt.rcParams['font.size'] = font_size_minor
    plt.rcParams['axes.labelsize'] = font_size_major
    plt.rcParams['axes.titlesize'] = font_size_major
    plt.rcParams['xtick.labelsize'] = font_size_minor
    plt.rcParams['ytick.labelsize'] = font_size_minor
    plt.rcParams['legend.fontsize'] = font_size_minor
    plt.rcParams['figure.titlesize'] = font_size_major
    plt.rcParams['xtick.major.size'] = 8
    plt.rcParams['xtick.minor.size'] = 4
    plt.rcParams['ytick.major.size'] = 8
    plt.rcParams['ytick.minor.size'] = 4
    plt.rcParams['xtick.major.width'] = 1.5
    plt.rcParams['xtick.minor.width'] = 1.0
    plt.rcParams['ytick.major.width'] = 1.5
    plt.rcParams['ytick.minor.width'] = 1.0
    plt.rcParams['axes.linewidth'] = 1.5


def load_clustering_result(result_file):
    """
    Load clustering results from a pickle file.

    Parameters
    ----------
    result_file : str or Path
        Path to the results pickle file

    Returns
    -------
    dict
        Dictionary containing clustering results with keys:
        - 'r': separation bins (arcmin)
        - 'xi': correlation function
        - 'sigma_xi': error bars (jackknife)
        - 'npairs': number of pairs
    """
    result_file = Path(result_file)

    if not result_file.exists():
        raise FileNotFoundError(f"Result file not found: {result_file}")

    with open(result_file, 'rb') as f:
        data = pickle.load(f)

    # Extract results from nested structure
    if 'results' in data:
        results = data['results']
    else:
        results = data

    return results


def get_redshift_bin_label(z_min, z_max):
    """
    Create a formatted redshift bin label.

    Parameters
    ----------
    z_min : float
        Minimum redshift
    z_max : float
        Maximum redshift

    Returns
    -------
    str
        Formatted label like "0.1 < z < 0.4"
    """
    return f"{z_min:.1f} < z < {z_max:.1f}"


def get_richness_bin_label(lambda_min, lambda_max):
    """
    Create a formatted richness bin label.

    Parameters
    ----------
    lambda_min : float
        Minimum richness
    lambda_max : float
        Maximum richness

    Returns
    -------
    str
        Formatted label like "40 < λ < 80"
    """
    if lambda_max >= 450:
        return f"{int(lambda_min)} < λ"
    return f"{int(lambda_min)} < λ < {int(lambda_max)}"


def plot_correlation_function(ax, r, xi, sigma_xi, label=None, color=None,
                               marker='o', markersize=6, show_errors=True,
                               linestyle='none', alpha=1.0):
    """
    Plot a correlation function with error bars on a given axis.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axis to plot on
    r : array-like
        Separation bins in arcmin
    xi : array-like
        Correlation function values
    sigma_xi : array-like
        Error bars (jackknife)
    label : str, optional
        Label for legend
    color : str, optional
        Color for the plot
    marker : str
        Marker style
    markersize : float
        Size of markers
    show_errors : bool
        Whether to show error bars
    linestyle : str
        Line style (default: 'none' for no lines)
    alpha : float
        Transparency
    """
    # Filter out invalid values
    valid = (r > 0) & np.isfinite(xi) & np.isfinite(sigma_xi) & (sigma_xi > 0)

    if show_errors:
        ax.errorbar(r[valid], xi[valid], yerr=sigma_xi[valid],
                   label=label, color=color, marker=marker,
                   markersize=markersize, linestyle=linestyle,
                   capsize=3, capthick=1.5, alpha=alpha)
    else:
        ax.plot(r[valid], xi[valid], label=label, color=color,
               marker=marker, markersize=markersize,
               linestyle=linestyle, alpha=alpha)


def format_clustering_axis(ax, xlabel=r'$\theta$ [arcmin]',
                           ylabel=r'$\omega_{cg}(\theta)$',
                           xscale='log', yscale='log',
                           xlim=None, ylim=None, grid=False):
    """
    Format an axis for clustering plots.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axis to format
    xlabel : str
        Label for x-axis
    ylabel : str
        Label for y-axis
    xscale : str
        Scale for x-axis ('log' or 'linear')
    yscale : str
        Scale for y-axis ('log' or 'linear')
    xlim : tuple, optional
        X-axis limits
    ylim : tuple, optional
        Y-axis limits
    grid : bool
        Whether to show grid (default: False)
    """
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_xscale(xscale)
    ax.set_yscale(yscale)

    if xlim is not None:
        ax.set_xlim(xlim)
    if ylim is not None:
        ax.set_ylim(ylim)

    if grid:
        ax.grid(True, alpha=0.3, which='both', linestyle='--', linewidth=0.5)

    ax.legend(frameon=True, framealpha=0.9)
    ax.minorticks_on()


def find_result_files(output_dir, galaxy_type='all', z_cluster_bin=None,
                      z_galaxy_bin=None, richness_bin=None):
    """
    Find result files matching specific bin criteria.

    Parameters
    ----------
    output_dir : str or Path
        Output directory to search
    galaxy_type : str
        Galaxy type filter: 'all', 'early', or 'late'
    z_cluster_bin : tuple, optional
        Cluster redshift bin (z_min, z_max)
    z_galaxy_bin : tuple, optional
        Galaxy redshift bin (z_min, z_max)
    richness_bin : tuple, optional
        Richness bin (lambda_min, lambda_max)

    Returns
    -------
    list
        List of Path objects for matching result files
    """
    output_dir = Path(output_dir)

    if not output_dir.exists():
        return []

    # Build pattern based on criteria
    # File format: clustering_c_RICHNESS_CLUSTER_40.00_80.00_Z_CLUSTER_0.10_0.40_g_redshift_0.10_0.40*.pkl
    pattern_parts = ["clustering_c_"]

    if richness_bin is not None:
        lambda_min, lambda_max = richness_bin
        pattern_parts.append(f"RICHNESS_CLUSTER_{lambda_min:.2f}_{lambda_max:.2f}_")
    else:
        pattern_parts.append("*")

    if z_cluster_bin is not None:
        z_min, z_max = z_cluster_bin
        pattern_parts.append(f"Z_CLUSTER_{z_min:.2f}_{z_max:.2f}_")
    else:
        pattern_parts.append("Z_CLUSTER_*_")

    if z_galaxy_bin is not None:
        z_min, z_max = z_galaxy_bin
        pattern_parts.append(f"g_redshift_{z_min:.2f}_{z_max:.2f}")
    else:
        pattern_parts.append("g_redshift_*")

    # Search for files (directly in output_dir, not in results/)
    pattern = "".join(pattern_parts) + "*.pkl"
    files = list(output_dir.glob(pattern))

    return sorted(files)


def get_galaxy_type_color(galaxy_type):
    """
    Get standard color for galaxy type.

    Parameters
    ----------
    galaxy_type : str
        'all', 'early', or 'late'

    Returns
    -------
    str
        Color code
    """
    colors = {
        'all': 'black',
        'early': '#e81010',  # Red for early-type
        'late': '#008fd6'    # Blue for late-type
    }
    return colors.get(galaxy_type, 'gray')


def get_galaxy_type_label(galaxy_type):
    """
    Get formatted label for galaxy type.

    Parameters
    ----------
    galaxy_type : str
        'all', 'early', or 'late'

    Returns
    -------
    str
        Formatted label
    """
    labels = {
        'all': 'All galaxies',
        'early': 'Early-type (n ≥ 2.0)',
        'late': 'Late-type (n < 2.0)'
    }
    return labels.get(galaxy_type, galaxy_type)
