"""
Configuration management for clustering pipeline.
"""

import yaml
from pathlib import Path
from typing import Dict, List, Any
import numpy as np


class PipelineConfig:
    """Configuration manager for clustering pipeline."""

    def __init__(self, config_path=None):
        """
        Initialize configuration from YAML file.

        Parameters
        ----------
        config_path : str or Path, optional
            Path to YAML configuration file
        """
        self.config = {}
        if config_path is not None:
            self.load_config(config_path)

    def load_config(self, config_path):
        """
        Load configuration from YAML file.

        Parameters
        ----------
        config_path : str or Path
            Path to configuration file
        """
        config_path = Path(config_path)
        with open(config_path, 'r') as f:
            self.config = yaml.safe_load(f)

        print(f"Configuration loaded from {config_path}")

    def save_config(self, output_path):
        """
        Save configuration to YAML file.

        Parameters
        ----------
        output_path : str or Path
            Output path for configuration file
        """
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        with open(output_path, 'w') as f:
            yaml.dump(self.config, f, default_flow_style=False, sort_keys=False)

        print(f"Configuration saved to {output_path}")

    def get_catalogues(self):
        """Get catalogue file paths."""
        return self.config.get('catalogues', {})

    def get_cosmology_params(self):
        """Get cosmology parameters."""
        cosmo_params = self.config.get('cosmology', {})
        return {
            'H0': cosmo_params.get('H0', 70),
            'Om0': cosmo_params.get('Om0', 0.3)
        }

    def get_analysis_params(self):
        """Get analysis parameters for correlation functions."""
        return self.config.get('analysis_parameters', {})

    def get_cluster_bins(self) -> List[Dict[str, Any]]:
        """
        Get cluster binning specifications.

        Returns
        -------
        list of dict
            List of filter dictionaries for cluster bins
        """
        bins_config = self.config.get('cluster_bins', {})
        return self._parse_bin_config(bins_config)

    def get_galaxy_bins(self) -> List[Dict[str, Any]]:
        """
        Get galaxy binning specifications.

        Returns
        -------
        list of dict
            List of filter dictionaries for galaxy bins
        """
        bins_config = self.config.get('galaxy_bins', {})
        return self._parse_bin_config(bins_config)

    def _parse_bin_config(self, bins_config):
        """
        Parse bin configuration into filter specifications.

        Parameters
        ----------
        bins_config : dict
            Binning configuration from YAML

        Returns
        -------
        list of dict
            List of filter dictionaries

        Examples
        --------
        Input config:
        {
            'redshift': {
                'bins': [0.0, 0.5, 1.0],
                'column': 'redshift'
            },
            'richness': {
                'bins': [10, 20, 50],
                'column': 'RICHNESS_CLUSTER_pzwav'
            }
        }

        Output: List of filter dicts for all combinations
        """
        if not bins_config:
            return [{}]  # No filtering

        # Separate properties and their bins
        properties = {}
        for prop_name, prop_config in bins_config.items():
            if isinstance(prop_config, dict) and 'bins' in prop_config:
                column = prop_config.get('column', prop_name)
                bins = prop_config['bins']
                properties[column] = bins
            elif isinstance(prop_config, dict) and 'values' in prop_config:
                # Categorical values
                column = prop_config.get('column', prop_name)
                values = prop_config['values']
                properties[column] = values

        # Generate all combinations of bins
        filter_list = self._generate_bin_combinations(properties)

        return filter_list

    def _generate_bin_combinations(self, properties):
        """
        Generate all combinations of bin filters.

        Parameters
        ----------
        properties : dict
            Dictionary mapping column names to bin edges

        Returns
        -------
        list of dict
            List of filter specifications
        """
        if not properties:
            return [{}]

        # Start with a single empty filter
        filters = [{}]

        for column, bins in properties.items():
            new_filters = []

            # Create bins from edges
            if isinstance(bins, list) and len(bins) > 1:
                # Check if it's bin edges (numeric) or categorical values
                if all(isinstance(b, (int, float)) for b in bins):
                    # Numeric bins - create ranges
                    for i in range(len(bins) - 1):
                        for existing_filter in filters:
                            new_filter = existing_filter.copy()
                            new_filter[column] = (bins[i], bins[i + 1])
                            new_filters.append(new_filter)
                else:
                    # Categorical values
                    for value in bins:
                        for existing_filter in filters:
                            new_filter = existing_filter.copy()
                            new_filter[column] = [value]
                            new_filters.append(new_filter)
            else:
                # Single value or range
                for existing_filter in filters:
                    new_filter = existing_filter.copy()
                    new_filter[column] = bins
                    new_filters.append(new_filter)

            filters = new_filters

        return filters

    def get_weight_config(self):
        """Get systematic weight configuration."""
        return self.config.get('systematic_weights', {})

    def get_output_dir(self):
        """Get output directory path."""
        return Path(self.config.get('output_directory', 'outputs'))

    def get_slurm_config(self):
        """Get SLURM job configuration."""
        return self.config.get('slurm', {})

    def get_column_mappings(self):
        """Get column name mappings for catalogues."""
        return self.config.get('column_mappings', {})

    def validate(self):
        """
        Validate the configuration.

        Returns
        -------
        bool
            True if configuration is valid

        Raises
        ------
        ValueError
            If configuration is invalid
        """
        required_sections = ['catalogues', 'analysis_parameters']

        for section in required_sections:
            if section not in self.config:
                raise ValueError(f"Missing required section: {section}")

        # Check catalogue paths exist
        catalogues = self.get_catalogues()
        for cat_type, cat_path in catalogues.items():
            if not Path(cat_path).exists():
                raise ValueError(f"{cat_type} catalogue not found: {cat_path}")

        # Check analysis parameters
        analysis = self.get_analysis_params()
        required_params = ['min_sep', 'max_sep', 'nbins']
        for param in required_params:
            if param not in analysis:
                raise ValueError(f"Missing required analysis parameter: {param}")

        # Validate mode-specific parameters
        mode = analysis.get('mode', '3d')
        if mode == '2d':
            # 2D mode should use angular metrics
            metric = analysis.get('metric', 'Euclidean')
            if metric not in ['Euclidean', 'Arc']:
                raise ValueError(f"2D mode requires metric='Euclidean' or 'Arc', got '{metric}'")

            # 2D mode shouldn't have line-of-sight parameters
            if 'min_rpar' in analysis or 'max_rpar' in analysis:
                print("WARNING: min_rpar/max_rpar are ignored in 2D mode")

            # Check sep_units if provided
            sep_units = analysis.get('sep_units')
            if sep_units and sep_units not in ['deg', 'arcmin', 'arcsec', 'radians']:
                raise ValueError(f"Invalid sep_units for 2D mode: {sep_units}")

        elif mode == '3d':
            # 3D mode should use Rperp metric
            metric = analysis.get('metric', 'Rperp')
            if metric != 'Rperp':
                print(f"WARNING: 3D mode typically uses metric='Rperp', got '{metric}'")

        print("Configuration validation passed")
        return True

    def get_all_job_specs(self):
        """
        Generate all job specifications for batch processing.

        Returns
        -------
        list of dict
            List of job specifications, each containing cluster and galaxy filters
        """
        cluster_bins = self.get_cluster_bins()
        galaxy_bins = self.get_galaxy_bins()

        job_specs = []
        job_id = 0

        for cluster_filter in cluster_bins:
            for galaxy_filter in galaxy_bins:
                job_specs.append({
                    'job_id': job_id,
                    'cluster_filters': cluster_filter,
                    'galaxy_filters': galaxy_filter
                })
                job_id += 1

        print(f"Generated {len(job_specs)} job specifications")
        return job_specs


def create_template_config(output_path='config_template.yaml'):
    """
    Create a template configuration file.

    Parameters
    ----------
    output_path : str
        Path for output template file
    """
    template = {
        'catalogues': {
            'clusters': '/path/to/cluster_catalogue.fits',
            'galaxies': '/path/to/galaxy_catalogue.fits',
            'randoms': '/path/to/random_catalogue.fits'
        },
        'column_mappings': {
            'cluster': {
                'ra': 'RIGHT_ASCENSION_CLUSTER_pzwav',
                'dec': 'DECLINATION_CLUSTER_pzwav',
                'redshift': 'Z_CLUSTER_pzwav'
            },
            'galaxy': {
                'ra': 'right_ascension',
                'dec': 'declination',
                'redshift': 'phz_median',
                'ellipticity': 'ellipticity',
                'position_angle': 'position_angle'
            },
            'random': {
                'ra': 'right_ascension',
                'dec': 'declination',
                'redshift': 'z'
            }
        },
        'cosmology': {
            'H0': 70,
            'Om0': 0.3
        },
        'analysis_parameters': {
            'min_sep': 0.1,
            'max_sep': 10.0,
            'nbins': 16,
            'min_rpar': -60,
            'max_rpar': 60,
            'bin_type': 'Log',
            'metric': 'Rperp',
            'mode': '3d'
        },
        'cluster_bins': {
            'redshift': {
                'column': 'Z_CLUSTER_pzwav',
                'bins': [0.2, 0.5, 0.8, 1.2]
            },
            'richness': {
                'column': 'RICHNESS_CLUSTER_pzwav',
                'bins': [10, 20, 50, 100]
            }
        },
        'galaxy_bins': {
            'redshift': {
                'column': 'redshift',
                'bins': [0.0, 0.5, 1.0, 1.5, 2.0]
            },
            'type': {
                'column': 'sersic_sersic_vis_index',
                'values': ['early', 'late'],
                'early': lambda x: x > 2.0,
                'late': lambda x: x < 2.0
            }
        },
        'systematic_weights': {
            'enabled': False,
            'weight_file': None,
            'depth_weights': {
                'column': 'limiting_magnitude'
            }
        },
        'output_directory': './outputs',
        'slurm': {
            'time': '02:00:00',
            'mem': '16G',
            'cpus': 4,
            'partition': 'batch',
            'account': None
        }
    }

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, 'w') as f:
        yaml.dump(template, f, default_flow_style=False, sort_keys=False)

    print(f"Template configuration created at {output_path}")
