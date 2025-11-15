#!/usr/bin/env python3
"""
Run clustering analysis for all bins defined in configuration.

This can be used to run the full analysis locally (not recommended for large jobs)
or to test the pipeline before submitting to SLURM.
"""

import argparse
import sys
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
        description='Run full clustering analysis')
    parser.add_argument('--config', type=str, required=True,
                       help='Path to configuration file')
    parser.add_argument('--dry-run', action='store_true',
                       help='Print what would be done without running')
    parser.add_argument('--test', action='store_true',
                       help='Run only the first job as a test')

    args = parser.parse_args()

    # Load configuration
    print("Loading configuration...")
    config = PipelineConfig(args.config)

    # Validate configuration
    print("Validating configuration...")
    config.validate()

    # Get all job specifications
    job_specs = config.get_all_job_specs()

    print(f"\nTotal jobs to run: {len(job_specs)}")

    if args.dry_run:
        print("\nDry run - jobs that would be executed:")
        for job in job_specs:
            print(f"  Job {job['job_id']}: "
                  f"Clusters={job['cluster_filters']}, "
                  f"Galaxies={job['galaxy_filters']}")
        return 0

    if args.test:
        print("\nTest mode - running only first job")
        job_specs = job_specs[:1]

    # Setup cosmology
    cosmo_params = config.get_cosmology_params()
    cosmology = FlatLambdaCDM(H0=cosmo_params['H0'], Om0=cosmo_params['Om0'])

    # Initialize catalogue manager
    cat_manager = CatalogueManager(cosmology=cosmology)

    # Load catalogues
    print("\nLoading catalogues...")
    catalogues = config.get_catalogues()
    col_mappings = config.get_column_mappings()

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
        print("\nComputing systematic weights...")
        if 'weight_file' in weight_config and weight_config['weight_file']:
            from src.weights import WeightCalculator
            galaxy_weights = WeightCalculator.load_weights(weight_config['weight_file'])
        else:
            galaxy_weights = apply_systematic_weights(galaxies, weight_config)

    # Get analysis parameters
    analysis_params = config.get_analysis_params()

    # Create output directory
    output_dir = config.get_output_dir()
    output_dir.mkdir(parents=True, exist_ok=True)

    # Run all jobs
    print("\n" + "="*70)
    print(f"RUNNING {len(job_specs)} JOBS")
    print("="*70)

    successful_jobs = 0
    failed_jobs = 0

    for i, job_spec in enumerate(job_specs, 1):
        print(f"\n{'='*70}")
        print(f"JOB {i}/{len(job_specs)} - ID: {job_spec['job_id']}")
        print(f"{'='*70}")

        try:
            results = compute_weighted_clustering(
                cluster_filters=job_spec['cluster_filters'],
                galaxy_filters=job_spec['galaxy_filters'],
                cluster_catalogue=clusters,
                galaxy_catalogue=galaxies,
                random_catalogue=randoms,
                catalogue_manager=cat_manager,
                analysis_params=analysis_params,
                output_dir=output_dir,
                galaxy_weights=galaxy_weights
            )

            if results is not None:
                successful_jobs += 1
                print(f"Job {job_spec['job_id']} completed successfully")
            else:
                print(f"Job {job_spec['job_id']} completed with no results (empty bin)")
                successful_jobs += 1

        except Exception as e:
            failed_jobs += 1
            print(f"Job {job_spec['job_id']} FAILED with error:")
            print(f"  {type(e).__name__}: {e}")
            import traceback
            traceback.print_exc()

    # Summary
    print("\n" + "="*70)
    print("ANALYSIS COMPLETE")
    print("="*70)
    print(f"Successful jobs: {successful_jobs}/{len(job_specs)}")
    print(f"Failed jobs: {failed_jobs}/{len(job_specs)}")
    print(f"Results saved to: {output_dir}")

    return 0 if failed_jobs == 0 else 1


if __name__ == '__main__':
    sys.exit(main())
