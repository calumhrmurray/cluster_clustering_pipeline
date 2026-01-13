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

    galaxy_filter_cols = set()
    for job_spec in job_specs:
        galaxy_filter_cols.update(job_spec['galaxy_filters'].keys())

    # Setup cosmology
    cosmo_params = config.get_cosmology_params()
    cosmology = FlatLambdaCDM(H0=cosmo_params['H0'], Om0=cosmo_params['Om0'])

    # Initialize catalogue manager
    cat_manager = CatalogueManager(cosmology=cosmology)

    # Load catalogues
    print("\nLoading catalogues...")
    catalogues = config.get_catalogues()
    col_mappings = config.get_column_mappings()
    analysis_params = config.get_analysis_params()
    mode = analysis_params.get('mode', '3d')

    cluster_col = col_mappings.get('cluster', {})
    clusters = cat_manager.load_cluster_catalogue(
        catalogues['clusters'],
        ra_col=cluster_col.get('ra', 'RIGHT_ASCENSION_CLUSTER_pzwav'),
        dec_col=cluster_col.get('dec', 'DECLINATION_CLUSTER_pzwav'),
        z_col=cluster_col.get('redshift', 'Z_CLUSTER_pzwav')
    )

    galaxy_col = col_mappings.get('galaxy', {})
    galaxy_ra = galaxy_col.get('ra', 'right_ascension')
    galaxy_dec = galaxy_col.get('dec', 'declination')
    galaxy_z = galaxy_col.get('redshift', 'phz_median')
    if mode == '2d':
        gal_cols_to_load = [galaxy_ra, galaxy_dec, galaxy_z]
        for filter_col in galaxy_filter_cols:
            if filter_col not in gal_cols_to_load and filter_col != 'redshift':
                gal_cols_to_load.append(filter_col)
        print(f"2D mode: Loading only {len(gal_cols_to_load)} galaxy columns to save memory")
    else:
        gal_cols_to_load = None

    galaxies = cat_manager.load_galaxy_catalogue(
        catalogues['galaxies'],
        ra_col=galaxy_ra,
        dec_col=galaxy_dec,
        z_col=galaxy_z,
        columns=gal_cols_to_load
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
    random_ra = random_col.get('ra', 'right_ascension')
    random_dec = random_col.get('dec', 'declination')
    random_z = random_col.get('redshift', 'z')
    if random_z in (None, ""):
        random_z = None
    if mode == '2d':
        random_cols_to_load = [random_ra, random_dec]
        if random_z:
            random_cols_to_load.append(random_z)
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

    # Apply optional HEALPix mask filtering
    mask_cfg = config.get_mask_filter_config()
    if mask_cfg and mask_cfg.get('enabled', True):
        from src.filters import CatalogueFilter

        mask_file = mask_cfg.get('map_file') or mask_cfg.get('mask_file')
        min_value = mask_cfg.get('min_value', mask_cfg.get('threshold'))
        max_value = mask_cfg.get('max_value')
        nest = mask_cfg.get('nest', False)
        nside = mask_cfg.get('nside')
        value_field = mask_cfg.get('value_field')
        targets = mask_cfg.get('targets', mask_cfg.get('apply_to', ['galaxies']))
        if isinstance(targets, str):
            targets = [targets]
        elif isinstance(targets, dict):
            targets = [name for name, enabled in targets.items() if enabled]

        print("\nApplying HEALPix mask filter...")
        if 'clusters' in targets:
            clusters = CatalogueFilter.apply_healpix_mask(
                clusters,
                mask_file,
                min_value=min_value,
                max_value=max_value,
                nest=nest,
                nside=nside,
                value_field=value_field,
                label='clusters'
            )
        if 'galaxies' in targets:
            galaxies = CatalogueFilter.apply_healpix_mask(
                galaxies,
                mask_file,
                min_value=min_value,
                max_value=max_value,
                nest=nest,
                nside=nside,
                value_field=value_field,
                label='galaxies'
            )
        if 'randoms' in targets:
            randoms = CatalogueFilter.apply_healpix_mask(
                randoms,
                mask_file,
                min_value=min_value,
                max_value=max_value,
                nest=nest,
                nside=nside,
                value_field=value_field,
                label='randoms'
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
