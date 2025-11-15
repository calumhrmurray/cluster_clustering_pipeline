"""
Systematic weight computation for galaxy catalogues.

This module provides infrastructure for applying systematic weights
to correct for observational effects like depth variations, seeing, etc.
"""

import numpy as np
from astropy.table import Table
import pickle
from pathlib import Path


class WeightCalculator:
    """Class for computing and applying systematic weights."""

    def __init__(self):
        """Initialize weight calculator."""
        self.weight_maps = {}
        self.weight_functions = {}

    def compute_uniform_weights(self, catalogue, weight_value=1.0):
        """
        Compute uniform weights for all objects.

        Parameters
        ----------
        catalogue : astropy.table.Table
            Input catalogue
        weight_value : float
            Weight value to assign

        Returns
        -------
        np.ndarray
            Array of weights
        """
        return np.full(len(catalogue), weight_value)

    def compute_depth_weights(self, catalogue, depth_column, reference_depth=None):
        """
        Compute weights based on survey depth variations.

        Weight = N_reference / N_local

        Parameters
        ----------
        catalogue : astropy.table.Table
            Galaxy catalogue
        depth_column : str
            Name of column containing depth information
        reference_depth : float, optional
            Reference depth. If None, use median depth.

        Returns
        -------
        np.ndarray
            Depth-based weights
        """
        depths = catalogue[depth_column]

        if reference_depth is None:
            reference_depth = np.median(depths[~np.isnan(depths)])

        # Simple linear weighting: deeper observations get lower weight
        # This assumes number density scales linearly with depth
        weights = reference_depth / depths

        # Handle NaNs and infinities
        weights = np.where(np.isfinite(weights), weights, 1.0)

        # Normalize to mean of 1
        weights /= np.mean(weights)

        return weights

    def compute_property_weights(self, catalogue, property_column,
                                 target_distribution=None, nbins=20):
        """
        Compute weights to match a target distribution in some property.

        This can be used to reweight galaxies to match a specific
        redshift distribution, color distribution, etc.

        Parameters
        ----------
        catalogue : astropy.table.Table
            Input catalogue
        property_column : str
            Column name for the property to reweight
        target_distribution : tuple of (bin_centers, target_counts), optional
            Target distribution. If None, returns uniform weights.
        nbins : int
            Number of bins for histogramming

        Returns
        -------
        np.ndarray
            Property-based weights
        """
        if target_distribution is None:
            return self.compute_uniform_weights(catalogue)

        prop_values = catalogue[property_column]
        target_bins, target_counts = target_distribution

        # Compute histogram of current distribution
        observed_counts, bin_edges = np.histogram(prop_values, bins=nbins,
                                                  range=(target_bins.min(),
                                                        target_bins.max()))

        # Avoid division by zero
        observed_counts = np.where(observed_counts > 0, observed_counts, 1)

        # Compute weight for each bin
        bin_weights = target_counts / observed_counts

        # Assign weights based on which bin each object falls into
        bin_indices = np.digitize(prop_values, bin_edges) - 1
        bin_indices = np.clip(bin_indices, 0, len(bin_weights) - 1)

        weights = bin_weights[bin_indices]

        # Handle NaNs
        weights = np.where(np.isfinite(weights), weights, 1.0)

        return weights

    def compute_spatial_weights(self, catalogue, ra_col='ra', dec_col='dec',
                               nside=64):
        """
        Compute spatial weights based on local density variations.

        This uses HEALPix pixelization to estimate local overdensity
        and assigns weights to flatten the distribution.

        Parameters
        ----------
        catalogue : astropy.table.Table
            Input catalogue
        ra_col : str
            RA column name
        dec_col : str
            Dec column name
        nside : int
            HEALPix nside parameter

        Returns
        -------
        np.ndarray
            Spatial weights

        Notes
        -----
        Requires healpy to be installed.
        """
        try:
            import healpy as hp
        except ImportError:
            print("Warning: healpy not installed. Returning uniform weights.")
            return self.compute_uniform_weights(catalogue)

        # Convert to HEALPix pixels
        theta = np.deg2rad(90.0 - catalogue[dec_col])
        phi = np.deg2rad(catalogue[ra_col])
        pixels = hp.ang2pix(nside, theta, phi)

        # Count objects per pixel
        unique_pixels, pixel_counts = np.unique(pixels, return_counts=True)
        pixel_map = dict(zip(unique_pixels, pixel_counts))

        # Mean density
        mean_density = len(catalogue) / len(unique_pixels)

        # Assign weights inversely proportional to local density
        object_pixels = np.array([pixel_map.get(p, 1) for p in pixels])
        weights = mean_density / object_pixels

        # Normalize
        weights /= np.mean(weights)

        return weights

    def combine_weights(self, *weight_arrays):
        """
        Combine multiple weight arrays multiplicatively.

        Parameters
        ----------
        *weight_arrays : np.ndarray
            Variable number of weight arrays

        Returns
        -------
        np.ndarray
            Combined weights
        """
        combined = np.ones(len(weight_arrays[0]))
        for weights in weight_arrays:
            combined *= weights

        # Normalize to mean of 1
        combined /= np.mean(combined)

        return combined

    def save_weights(self, weights, output_path, metadata=None):
        """
        Save weights to file.

        Parameters
        ----------
        weights : np.ndarray
            Weight array
        output_path : str or Path
            Output file path
        metadata : dict, optional
            Metadata about weight calculation
        """
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        save_dict = {
            'weights': weights,
            'metadata': metadata or {}
        }

        if output_path.suffix == '.npy':
            np.save(output_path, weights)
        else:
            with open(output_path, 'wb') as f:
                pickle.dump(save_dict, f)

        print(f"Weights saved to {output_path}")

    @staticmethod
    def load_weights(input_path):
        """
        Load weights from file.

        Parameters
        ----------
        input_path : str or Path
            Input file path

        Returns
        -------
        np.ndarray or dict
            Loaded weights (array if .npy, dict if pickle)
        """
        input_path = Path(input_path)

        if input_path.suffix == '.npy':
            return np.load(input_path)
        else:
            with open(input_path, 'rb') as f:
                return pickle.load(f)


def apply_systematic_weights(galaxy_catalogue, weight_config):
    """
    Apply systematic weights to a galaxy catalogue based on configuration.

    Parameters
    ----------
    galaxy_catalogue : astropy.table.Table
        Galaxy catalogue
    weight_config : dict
        Configuration dictionary specifying weight calculations.
        Example:
        {
            'depth_weights': {'column': 'limiting_magnitude'},
            'spatial_weights': {'nside': 64},
            'property_weights': {'column': 'color', 'target': (bins, counts)}
        }

    Returns
    -------
    np.ndarray
        Systematic weights array
    """
    calculator = WeightCalculator()
    weight_arrays = []

    # Compute different types of weights
    if 'depth_weights' in weight_config:
        params = weight_config['depth_weights']
        weights = calculator.compute_depth_weights(
            galaxy_catalogue, **params)
        weight_arrays.append(weights)
        print(f"Depth weights computed. Mean={np.mean(weights):.3f}, "
              f"Std={np.std(weights):.3f}")

    if 'spatial_weights' in weight_config:
        params = weight_config['spatial_weights']
        weights = calculator.compute_spatial_weights(
            galaxy_catalogue, **params)
        weight_arrays.append(weights)
        print(f"Spatial weights computed. Mean={np.mean(weights):.3f}, "
              f"Std={np.std(weights):.3f}")

    if 'property_weights' in weight_config:
        params = weight_config['property_weights']
        weights = calculator.compute_property_weights(
            galaxy_catalogue, **params)
        weight_arrays.append(weights)
        print(f"Property weights computed. Mean={np.mean(weights):.3f}, "
              f"Std={np.std(weights):.3f}")

    # Combine all weights
    if len(weight_arrays) == 0:
        return calculator.compute_uniform_weights(galaxy_catalogue)
    elif len(weight_arrays) == 1:
        return weight_arrays[0]
    else:
        combined = calculator.combine_weights(*weight_arrays)
        print(f"Combined weights: Mean={np.mean(combined):.3f}, "
              f"Std={np.std(combined):.3f}")
        return combined
