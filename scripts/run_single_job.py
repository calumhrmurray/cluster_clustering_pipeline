#!/usr/bin/env python3
"""
Run a single clustering analysis job.

This script is called by SLURM job scripts to execute a single
cluster-galaxy bin combination.
"""

import argparse
import sys
import yaml
from pathlib import Path
from astropy.cosmology import FlatLambdaCDM

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from src.config import PipelineConfig
from src.catalogue import CatalogueManager
from src.clustering import compute_weighted_clustering
from src.weights import apply_systematic_weights


def main():
    parser = argparse.ArgumentParser(
        description='Run clustering analysis for a single job')
    parser.add_argument('--config', type=str, required=True,
                       help='Path to pipeline configuration file')
    parser.add_argument('--job-spec', type=str, required=True,
                       help='Path to job specification YAML file')

    args = parser.parse_args()

    print("="*70)
    print("CLUSTERING ANALYSIS - SINGLE JOB")
    print("="*70)

    # Load configuration
    print("\nLoading configuration...")
    config = PipelineConfig(args.config)

    # Load job specification
    print(f"Loading job specification from {args.job_spec}")
    with open(args.job_spec, 'r') as f:
        job_spec = yaml.safe_load(f)

    job_id = job_spec['job_id']
    cluster_filters = job_spec['cluster_filters']
    galaxy_filters = job_spec['galaxy_filters']

    print(f"\nJob ID: {job_id}")
    print(f"Cluster filters: {cluster_filters}")
    print(f"Galaxy filters: {galaxy_filters}")

    # Setup cosmology
    cosmo_params = config.get_cosmology_params()
    cosmology = FlatLambdaCDM(H0=cosmo_params['H0'], Om0=cosmo_params['Om0'])
    print(f"\nCosmology: H0={cosmo_params['H0']}, Om0={cosmo_params['Om0']}")

    # Initialize catalogue manager
    cat_manager = CatalogueManager(cosmology=cosmology)

    # Load catalogues
    catalogues = config.get_catalogues()
    col_mappings = config.get_column_mappings()
    analysis_params = config.get_analysis_params()
    mode = analysis_params.get('mode', '3d')

    print("\n" + "-"*70)
    print("LOADING CATALOGUES")
    print("-"*70)

    # Determine which columns to load for memory efficiency
    cluster_col = col_mappings.get('cluster', {})
    cluster_ra = cluster_col.get('ra', 'RIGHT_ASCENSION_CLUSTER_pzwav')
    cluster_dec = cluster_col.get('dec', 'DECLINATION_CLUSTER_pzwav')
    cluster_z = cluster_col.get('redshift', 'Z_CLUSTER_pzwav')

    # For clusters, also need filter columns
    cluster_cols_to_load = [cluster_ra, cluster_dec, cluster_z]
    for filter_col in cluster_filters.keys():
        if filter_col not in cluster_cols_to_load:
            cluster_cols_to_load.append(filter_col)

    clusters = cat_manager.load_cluster_catalogue(
        catalogues['clusters'],
        ra_col=cluster_ra,
        dec_col=cluster_dec,
        z_col=cluster_z
    )

    galaxy_col = col_mappings.get('galaxy', {})
    galaxy_ra = galaxy_col.get('ra', 'right_ascension')
    galaxy_dec = galaxy_col.get('dec', 'declination')
    galaxy_z = galaxy_col.get('redshift', 'phz_median')

    # For 2D mode, only load RA, Dec, redshift (and filter columns)
    # This saves significant memory (62M rows × fewer columns)
    if mode == '2d':
        gal_cols_to_load = [galaxy_ra, galaxy_dec, galaxy_z]
        # Add any filter columns
        for filter_col in galaxy_filters.keys():
            if filter_col not in gal_cols_to_load and filter_col != 'redshift':
                gal_cols_to_load.append(filter_col)
        print(f"2D mode: Loading only {len(gal_cols_to_load)} galaxy columns to save memory")
    else:
        gal_cols_to_load = None  # Load all columns for 3D mode

    galaxies = cat_manager.load_galaxy_catalogue(
        catalogues['galaxies'],
        ra_col=galaxy_ra,
        dec_col=galaxy_dec,
        z_col=galaxy_z,
        columns=gal_cols_to_load
    )

    # Add shapes if available and needed
    if mode == '3d' and 'ellipticity' in galaxy_col:
        if galaxy_col['ellipticity'] in galaxies.colnames:
            print("Adding galaxy shape information...")
            galaxies = cat_manager.add_galaxy_shapes(
                galaxies,
                ellipticity_col=galaxy_col['ellipticity'],
                pa_col=galaxy_col.get('position_angle', 'position_angle')
            )

    random_col = col_mappings.get('random', {})
    random_ra = random_col.get('ra', 'right_ascension')
    random_dec = random_col.get('dec', 'declination')
    random_z = random_col.get('redshift', 'z')

    # For randoms, only load RA, Dec, redshift in 2D mode
    if mode == '2d':
        random_cols_to_load = [random_ra, random_dec, random_z]
        print(f"2D mode: Loading only {len(random_cols_to_load)} random columns to save memory")
    else:
        random_cols_to_load = None

    randoms = cat_manager.load_random_catalogue(
        catalogues['randoms'],
        ra_col=random_ra,
        dec_col=random_dec,
        z_col=random_z,
        columns=random_cols_to_load
    )

    # Apply systematic weights if configured
    galaxy_weights = None
    weight_config = config.get_weight_config()
    if weight_config.get('enabled', False):
        print("\n" + "-"*70)
        print("COMPUTING SYSTEMATIC WEIGHTS")
        print("-"*70)

        if 'weight_file' in weight_config and weight_config['weight_file']:
            # Load precomputed weights
            print(f"Loading weights from {weight_config['weight_file']}")
            from src.weights import WeightCalculator
            galaxy_weights = WeightCalculator.load_weights(weight_config['weight_file'])
        else:
            # Compute weights on the fly
            print("Computing systematic weights...")
            galaxy_weights = apply_systematic_weights(galaxies, weight_config)

    # Get analysis parameters
    analysis_params = config.get_analysis_params()

    # Create output directory
    output_dir = config.get_output_dir()
    output_dir.mkdir(parents=True, exist_ok=True)

    print("\n" + "-"*70)
    print("RUNNING CLUSTERING ANALYSIS")
    print("-"*70)

    # Run the analysis
    results = compute_weighted_clustering(
        cluster_filters=cluster_filters,
        galaxy_filters=galaxy_filters,
        cluster_catalogue=clusters,
        galaxy_catalogue=galaxies,
        random_catalogue=randoms,
        catalogue_manager=cat_manager,
        analysis_params=analysis_params,
        output_dir=output_dir,
        galaxy_weights=galaxy_weights
    )

    if results is not None:
        print("\n" + "="*70)
        print("JOB COMPLETED SUCCESSFULLY")
        print("="*70)
        return 0
    else:
        print("\n" + "="*70)
        print("JOB COMPLETED WITH NO RESULTS (empty bin)")
        print("="*70)
        return 0


if __name__ == '__main__':
    sys.exit(main())
