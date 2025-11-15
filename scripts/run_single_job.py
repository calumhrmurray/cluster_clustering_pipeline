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

    print("\n" + "-"*70)
    print("LOADING CATALOGUES")
    print("-"*70)

    cluster_col = col_mappings.get('cluster', {})
    clusters = cat_manager.load_cluster_catalogue(
        catalogues['clusters'],
        ra_col=cluster_col.get('ra', 'RIGHT_ASCENSION_CLUSTER_pzwav'),
        dec_col=cluster_col.get('dec', 'DECLINATION_CLUSTER_pzwav'),
        z_col=cluster_col.get('redshift', 'Z_CLUSTER_pzwav')
    )

    galaxy_col = col_mappings.get('galaxy', {})
    galaxies = cat_manager.load_galaxy_catalogue(
        catalogues['galaxies'],
        ra_col=galaxy_col.get('ra', 'right_ascension'),
        dec_col=galaxy_col.get('dec', 'declination'),
        z_col=galaxy_col.get('redshift', 'phz_median')
    )

    # Add shapes if available
    if 'ellipticity' in galaxy_col and galaxy_col['ellipticity'] in galaxies.colnames:
        print("Adding galaxy shape information...")
        galaxies = cat_manager.add_galaxy_shapes(
            galaxies,
            ellipticity_col=galaxy_col['ellipticity'],
            pa_col=galaxy_col.get('position_angle', 'position_angle')
        )

    random_col = col_mappings.get('random', {})
    randoms = cat_manager.load_random_catalogue(
        catalogues['randoms'],
        ra_col=random_col.get('ra', 'right_ascension'),
        dec_col=random_col.get('dec', 'declination'),
        z_col=random_col.get('redshift', 'z')
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
