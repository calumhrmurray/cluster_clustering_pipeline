"""
Shape-position correlation measurements for weak lensing and intrinsic alignments.
"""

import numpy as np
import treecorr
import pickle
from pathlib import Path


class ShapeCorrelationAnalysis:
    """Class for computing shape-position correlations."""

    def __init__(self, min_sep, max_sep, nbins, min_rpar=None, max_rpar=None,
                 bin_type='Log', sep_units=None, metric='Euclidean'):
        """
        Initialize shape correlation analysis parameters.

        Parameters
        ----------
        min_sep : float
            Minimum separation
        max_sep : float
            Maximum separation
        nbins : int
            Number of bins
        min_rpar : float, optional
            Minimum separation parallel to line of sight
        max_rpar : float, optional
            Maximum separation parallel to line of sight
        bin_type : str
            'Log' or 'Linear' binning
        sep_units : str, optional
            Units for separation (e.g., 'deg', 'arcmin')
        metric : str
            Distance metric ('Euclidean', 'Rperp', etc.)
            Note: For weak lensing projections, use 'Euclidean' with angular coordinates
        """
        self.min_sep = min_sep
        self.max_sep = max_sep
        self.nbins = nbins
        self.min_rpar = min_rpar
        self.max_rpar = max_rpar
        self.bin_type = bin_type
        self.sep_units = sep_units
        self.metric = metric

    def compute_ng_correlation(self, position_cat, shape_cat, random_cat):
        """
        Compute the position-shape (NG) cross-correlation.

        This can be used for:
        - Weak lensing measurements (lenses at low-z, sources at high-z)
        - Intrinsic alignment measurements (same galaxy population)

        Parameters
        ----------
        position_cat : treecorr.Catalog
            Catalogue of positions (e.g., clusters or lens galaxies)
        shape_cat : treecorr.Catalog
            Catalogue of galaxies with shape information
        random_cat : treecorr.Catalog
            Random catalogue for positions

        Returns
        -------
        dict
            Dictionary containing:
            - 'r': separation bins
            - 'xi_plus': tangential shear (gamma_t)
            - 'xi_cross': cross component (gamma_x)
            - 'var_xi': variance
            - 'npairs': number of pairs
        """
        print(f"Computing NG correlation with {self.metric} metric")
        print(f"  Separation range: {self.min_sep} - {self.max_sep}")

        # Create correlation objects
        corr_kwargs = {
            'min_sep': self.min_sep,
            'max_sep': self.max_sep,
            'nbins': self.nbins,
            'bin_type': self.bin_type,
        }

        if self.sep_units is not None:
            corr_kwargs['sep_units'] = self.sep_units

        if self.min_rpar is not None:
            corr_kwargs['min_rpar'] = self.min_rpar
            corr_kwargs['max_rpar'] = self.max_rpar

        # Create NG correlation objects
        ng = treecorr.NGCorrelation(**corr_kwargs)
        rg = treecorr.NGCorrelation(**corr_kwargs)

        # Process correlations
        print("Processing NG (data-shapes)...")
        ng.process(position_cat, shape_cat, metric=self.metric)

        print("Processing RG (random-shapes)...")
        rg.process(random_cat, shape_cat, metric=self.metric)

        # Calculate correlation with random subtraction
        print("Calculating shape correlation...")
        xi_plus, xi_cross, var_xi = ng.calculateXi(rg=rg)

        # Get mean separation
        r = np.exp(ng.meanlogr)

        # Get number of pairs
        npairs = ng.npairs

        print(f"Shape correlation complete. Non-zero bins: {np.sum(npairs > 0)}/{self.nbins}")

        return {
            'r': r,
            'xi_plus': xi_plus,      # Tangential component (gamma_t)
            'xi_cross': xi_cross,    # Cross component (gamma_x, should be ~0)
            'var_xi': var_xi,
            'npairs': npairs,
            'sigma_xi': np.sqrt(var_xi)
        }

    def compute_weak_lensing(self, lens_cat, source_cat, random_cat):
        """
        Compute weak lensing signal around lenses.

        Parameters
        ----------
        lens_cat : treecorr.Catalog
            Lens catalogue (should use angular coordinates)
        source_cat : treecorr.Catalog
            Source catalogue with shapes (should use angular coordinates)
        random_cat : treecorr.Catalog
            Random catalogue

        Returns
        -------
        dict
            Lensing results with tangential and cross shear profiles
        """
        print("Computing weak lensing signal")
        return self.compute_ng_correlation(lens_cat, source_cat, random_cat)

    def compute_intrinsic_alignment(self, position_cat, shape_cat, random_cat):
        """
        Compute intrinsic alignment signal.

        For IA, both position and shape catalogues are typically from
        the same galaxy population.

        Parameters
        ----------
        position_cat : treecorr.Catalog
            Galaxy positions
        shape_cat : treecorr.Catalog
            Galaxy shapes (can be the same as position_cat)
        random_cat : treecorr.Catalog
            Random catalogue

        Returns
        -------
        dict
            IA results
        """
        print("Computing intrinsic alignment signal")
        return self.compute_ng_correlation(position_cat, shape_cat, random_cat)

    def save_results(self, results, output_path, metadata=None):
        """
        Save shape correlation results to file.

        Parameters
        ----------
        results : dict
            Results dictionary
        output_path : str or Path
            Path to output file
        metadata : dict, optional
            Additional metadata
        """
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        save_dict = {
            'results': results,
            'parameters': {
                'min_sep': self.min_sep,
                'max_sep': self.max_sep,
                'nbins': self.nbins,
                'min_rpar': self.min_rpar,
                'max_rpar': self.max_rpar,
                'bin_type': self.bin_type,
                'sep_units': self.sep_units,
                'metric': self.metric
            }
        }

        if metadata is not None:
            save_dict['metadata'] = metadata

        with open(output_path, 'wb') as f:
            pickle.dump(save_dict, f)

        print(f"Shape correlation results saved to {output_path}")

    @staticmethod
    def load_results(input_path):
        """
        Load shape correlation results from file.

        Parameters
        ----------
        input_path : str or Path
            Path to results file

        Returns
        -------
        dict
            Loaded results
        """
        with open(input_path, 'rb') as f:
            return pickle.load(f)


def compute_cluster_lensing(cluster_filters, source_filters,
                           cluster_catalogue, galaxy_catalogue,
                           random_catalogue, catalogue_manager,
                           analysis_params, output_dir):
    """
    Convenience function to compute weak lensing around clusters.

    Parameters
    ----------
    cluster_filters : dict
        Filters for cluster selection
    source_filters : dict
        Filters for source galaxy selection
    cluster_catalogue : astropy.table.Table
        Cluster catalogue
    galaxy_catalogue : astropy.table.Table
        Galaxy catalogue with shapes
    random_catalogue : astropy.table.Table
        Random catalogue
    catalogue_manager : CatalogueManager
        Catalogue manager instance
    analysis_params : dict
        Analysis parameters
    output_dir : str or Path
        Output directory

    Returns
    -------
    dict
        Lensing results
    """
    from .filters import CatalogueFilter

    print("\n" + "="*60)
    print("CLUSTER FILTERS:")
    filtered_clusters = CatalogueFilter.apply_filters(cluster_catalogue,
                                                      cluster_filters)

    print("\nSOURCE GALAXY FILTERS:")
    filtered_sources = CatalogueFilter.apply_filters(galaxy_catalogue,
                                                     source_filters)

    if len(filtered_clusters) == 0 or len(filtered_sources) == 0:
        print("Warning: No objects after filtering. Skipping.")
        return None

    # For weak lensing, use angular coordinates
    print("\nCreating angular catalogues for weak lensing...")

    source_cat = catalogue_manager.create_treecorr_catalogue(
        filtered_sources, use_shapes=True, use_weights=True, mode='2d')

    cluster_cat = catalogue_manager.create_treecorr_catalogue(
        filtered_clusters, patch_centers=source_cat.patch_centers, mode='2d')

    random_cat = catalogue_manager.create_treecorr_catalogue(
        random_catalogue, patch_centers=source_cat.patch_centers, mode='2d')

    # Create analysis object
    analysis = ShapeCorrelationAnalysis(**analysis_params)

    # Compute lensing
    results = analysis.compute_weak_lensing(cluster_cat, source_cat, random_cat)

    # Save results
    cluster_label = CatalogueFilter.get_bin_label(cluster_filters)
    source_label = CatalogueFilter.get_bin_label(source_filters)
    output_filename = f"lensing_c_{cluster_label}_s_{source_label}.pkl"
    output_path = Path(output_dir) / output_filename

    metadata = {
        'cluster_filters': cluster_filters,
        'source_filters': source_filters,
        'n_clusters': len(filtered_clusters),
        'n_sources': len(filtered_sources),
        'analysis_type': 'weak_lensing'
    }

    analysis.save_results(results, output_path, metadata=metadata)

    return results
