"""
Weighted clustering measurements using TreeCorr.
"""

import numpy as np
import treecorr
import pickle
from pathlib import Path


class ClusteringAnalysis:
    """Class for computing weighted clustering around clusters."""

    def __init__(self, min_sep, max_sep, nbins, min_rpar=None, max_rpar=None,
                 bin_type='Log', sep_units=None, metric='Rperp'):
        """
        Initialize clustering analysis parameters.

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
            Distance metric ('Rperp', 'Euclidean', etc.)
        """
        self.min_sep = min_sep
        self.max_sep = max_sep
        self.nbins = nbins
        self.min_rpar = min_rpar
        self.max_rpar = max_rpar
        self.bin_type = bin_type
        self.sep_units = sep_units
        self.metric = metric

    def compute_clustering(self, cluster_cat, galaxy_cat, random_cat,
                          use_weights=False):
        """
        Compute the cluster-galaxy cross-correlation function.

        Uses the Landy-Szalay estimator:
        xi = (DD - DR - RD + RR) / RR

        Parameters
        ----------
        cluster_cat : treecorr.Catalog
            Cluster catalogue
        galaxy_cat : treecorr.Catalog
            Galaxy catalogue (with optional weights)
        random_cat : treecorr.Catalog
            Random catalogue
        use_weights : bool
            Whether to use weights in the correlation

        Returns
        -------
        dict
            Dictionary containing:
            - 'r': separation bins
            - 'xi': correlation function
            - 'var_xi': variance of correlation function
            - 'npairs': number of pairs in each bin
        """
        print(f"Computing clustering with {self.metric} metric")
        print(f"  Separation range: {self.min_sep} - {self.max_sep}")
        if self.min_rpar is not None:
            print(f"  R_parallel range: {self.min_rpar} - {self.max_rpar}")

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

        # Create four correlation objects for Landy-Szalay
        nn = treecorr.NNCorrelation(**corr_kwargs)  # Data-Data
        nr = treecorr.NNCorrelation(**corr_kwargs)  # Data-Random
        rn = treecorr.NNCorrelation(**corr_kwargs)  # Random-Data
        rr = treecorr.NNCorrelation(**corr_kwargs)  # Random-Random

        # Process correlations
        print("Processing DD...")
        nn.process(cluster_cat, galaxy_cat, metric=self.metric)

        print("Processing DR...")
        nr.process(cluster_cat, random_cat, metric=self.metric)

        print("Processing RD...")
        rn.process(random_cat, galaxy_cat, metric=self.metric)

        print("Processing RR...")
        rr.process(random_cat, random_cat, metric=self.metric)

        # Calculate Landy-Szalay estimator
        print("Calculating correlation function...")
        xi, var_xi = nn.calculateXi(rr=rr, dr=nr, rd=rn)

        # Get mean separation
        r = np.exp(nn.meanlogr)

        # Get number of pairs
        npairs = nn.npairs

        print(f"Clustering measurement complete. Non-zero bins: {np.sum(npairs > 0)}/{self.nbins}")

        return {
            'r': r,
            'xi': xi,
            'var_xi': var_xi,
            'npairs': npairs,
            'sigma_xi': np.sqrt(var_xi)
        }

    def compute_autocorrelation(self, galaxy_cat, random_cat, use_weights=False):
        """
        Compute the galaxy auto-correlation function.

        Parameters
        ----------
        galaxy_cat : treecorr.Catalog
            Galaxy catalogue
        random_cat : treecorr.Catalog
            Random catalogue
        use_weights : bool
            Whether to use weights

        Returns
        -------
        dict
            Dictionary with correlation function results
        """
        print("Computing galaxy auto-correlation")
        return self.compute_clustering(galaxy_cat, galaxy_cat, random_cat,
                                       use_weights=use_weights)

    def save_results(self, results, output_path, metadata=None):
        """
        Save clustering results to file.

        Parameters
        ----------
        results : dict
            Results dictionary from compute_clustering
        output_path : str or Path
            Path to output file
        metadata : dict, optional
            Additional metadata to save
        """
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        # Prepare data to save
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

        # Save as pickle for Python or numpy format
        if output_path.suffix == '.pkl':
            with open(output_path, 'wb') as f:
                pickle.dump(save_dict, f)
        elif output_path.suffix == '.npz':
            np.savez(output_path, **save_dict)
        else:
            # Default to pickle
            with open(output_path.with_suffix('.pkl'), 'wb') as f:
                pickle.dump(save_dict, f)

        print(f"Results saved to {output_path}")

    @staticmethod
    def load_results(input_path):
        """
        Load clustering results from file.

        Parameters
        ----------
        input_path : str or Path
            Path to results file

        Returns
        -------
        dict
            Loaded results dictionary
        """
        input_path = Path(input_path)

        if input_path.suffix == '.pkl':
            with open(input_path, 'rb') as f:
                return pickle.load(f)
        elif input_path.suffix == '.npz':
            data = np.load(input_path, allow_pickle=True)
            return dict(data)
        else:
            raise ValueError(f"Unknown file format: {input_path.suffix}")


def compute_weighted_clustering(cluster_filters, galaxy_filters,
                                cluster_catalogue, galaxy_catalogue,
                                random_catalogue, catalogue_manager,
                                analysis_params, output_dir,
                                galaxy_weights=None):
    """
    Convenience function to compute weighted clustering for specific bins.

    Parameters
    ----------
    cluster_filters : dict
        Filters for cluster selection
    galaxy_filters : dict
        Filters for galaxy selection
    cluster_catalogue : astropy.table.Table
        Full cluster catalogue
    galaxy_catalogue : astropy.table.Table
        Full galaxy catalogue
    random_catalogue : astropy.table.Table
        Random catalogue
    catalogue_manager : CatalogueManager
        Catalogue manager instance
    analysis_params : dict
        Dictionary of analysis parameters (min_sep, max_sep, etc.)
    output_dir : str or Path
        Directory for output files
    galaxy_weights : array-like, optional
        Systematic weights for galaxies

    Returns
    -------
    dict
        Clustering results
    """
    from .filters import CatalogueFilter

    # Apply filters
    print("\n" + "="*60)
    print("CLUSTER FILTERS:")
    filtered_clusters = CatalogueFilter.apply_filters(cluster_catalogue,
                                                      cluster_filters)

    print("\nGALAXY FILTERS:")
    filtered_galaxies = CatalogueFilter.apply_filters(galaxy_catalogue,
                                                      galaxy_filters)

    if len(filtered_clusters) == 0 or len(filtered_galaxies) == 0:
        print("Warning: No objects after filtering. Skipping this bin.")
        return None

    # Apply systematic weights if provided
    if galaxy_weights is not None:
        if 'w' not in filtered_galaxies.colnames:
            filtered_galaxies['w'] = np.ones(len(filtered_galaxies))
        filtered_galaxies['w'] *= galaxy_weights

    # Get mode from analysis parameters
    mode = analysis_params.get('mode', '3d')

    # Validate mode and parameters
    if mode == '2d':
        # For 2D mode, check that we're not using 3D-only parameters
        if 'min_rpar' in analysis_params or 'max_rpar' in analysis_params:
            print("WARNING: min_rpar/max_rpar are ignored in 2D mode")
        if analysis_params.get('metric') not in ['Euclidean', 'Arc', None]:
            print(f"WARNING: metric='{analysis_params.get('metric')}' unusual for 2D mode")

    # Add coordinates (only Cartesian for 3D mode)
    if mode == '3d':
        print("\nAdding Cartesian coordinates for 3D mode...")
        filtered_clusters = catalogue_manager.add_cartesian_coordinates(filtered_clusters)
        filtered_galaxies = catalogue_manager.add_cartesian_coordinates(filtered_galaxies)
        random_catalogue = catalogue_manager.add_cartesian_coordinates(random_catalogue)
    else:
        print("\nUsing angular coordinates for 2D mode...")

    # Filter randoms by galaxy redshift bin for tomographic analysis
    if 'redshift' in galaxy_filters:
        z_range = galaxy_filters['redshift']
        if isinstance(z_range, tuple) and len(z_range) == 2:
            z_min, z_max = z_range
            n_randoms_before = len(random_catalogue)
            random_catalogue = random_catalogue[
                (random_catalogue['redshift'] >= z_min) &
                (random_catalogue['redshift'] < z_max)
            ]
            n_randoms_after = len(random_catalogue)
            print(f"Filtered randoms by galaxy redshift bin [{z_min:.2f}, {z_max:.2f}): "
                  f"{n_randoms_before} → {n_randoms_after}")

    print(f"\nCATALOGUE SIZES:")
    print(f"  Filtered clusters: {len(filtered_clusters)}")
    print(f"  Filtered galaxies: {len(filtered_galaxies)}")
    print(f"  Randoms: {len(random_catalogue)}")

    # Create TreeCorr catalogues
    use_shapes = 'e1' in filtered_galaxies.colnames
    use_weights = 'w' in filtered_galaxies.colnames

    print("\nCreating TreeCorr catalogues...")
    galaxy_cat = catalogue_manager.create_treecorr_catalogue(
        filtered_galaxies, use_shapes=use_shapes, use_weights=use_weights,
        mode=analysis_params.get('mode', '3d'))

    cluster_cat = catalogue_manager.create_treecorr_catalogue(
        filtered_clusters, patch_centers=galaxy_cat.patch_centers,
        mode=analysis_params.get('mode', '3d'))

    random_cat = catalogue_manager.create_treecorr_catalogue(
        random_catalogue, patch_centers=galaxy_cat.patch_centers,
        mode=analysis_params.get('mode', '3d'))

    # Create analysis object (filter out mode parameter)
    clustering_params = {k: v for k, v in analysis_params.items() if k != 'mode'}
    analysis = ClusteringAnalysis(**clustering_params)

    # Compute clustering
    results = analysis.compute_clustering(cluster_cat, galaxy_cat, random_cat,
                                         use_weights=use_weights)

    # Save results
    cluster_label = CatalogueFilter.get_bin_label(cluster_filters)
    galaxy_label = CatalogueFilter.get_bin_label(galaxy_filters)
    output_filename = f"clustering_c_{cluster_label}_g_{galaxy_label}.pkl"
    output_path = Path(output_dir) / output_filename

    metadata = {
        'cluster_filters': cluster_filters,
        'galaxy_filters': galaxy_filters,
        'n_clusters': len(filtered_clusters),
        'n_galaxies': len(filtered_galaxies),
        'n_randoms': len(random_catalogue),
        'use_weights': use_weights
    }

    analysis.save_results(results, output_path, metadata=metadata)

    return results
