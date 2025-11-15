"""
Catalogue management module for loading and preparing cluster and galaxy catalogues.
"""

import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy.cosmology import FlatLambdaCDM
import treecorr


class CatalogueManager:
    """Manager for loading and preparing astronomical catalogues."""

    def __init__(self, cosmology=None):
        """
        Initialize the catalogue manager.

        Parameters
        ----------
        cosmology : astropy.cosmology object, optional
            Cosmology to use for distance calculations.
            Default: FlatLambdaCDM(H0=70, Om0=0.3)
        """
        if cosmology is None:
            self.cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
        else:
            self.cosmo = cosmology

    def load_cluster_catalogue(self, filepath, ra_col='RIGHT_ASCENSION_CLUSTER_pzwav',
                               dec_col='DECLINATION_CLUSTER_pzwav',
                               z_col='Z_CLUSTER_pzwav'):
        """
        Load a cluster catalogue from FITS file.

        Parameters
        ----------
        filepath : str
            Path to the FITS file
        ra_col : str
            Name of the RA column
        dec_col : str
            Name of the declination column
        z_col : str
            Name of the redshift column

        Returns
        -------
        astropy.table.Table
            The loaded catalogue
        """
        print(f"Loading cluster catalogue from {filepath}")
        clusters = Table(fits.open(filepath)[1].data)

        # Add coordinate columns
        clusters['ra'] = clusters[ra_col]
        clusters['dec'] = clusters[dec_col]
        clusters['redshift'] = clusters[z_col]

        print(f"Loaded {len(clusters)} clusters")
        return clusters

    def load_galaxy_catalogue(self, filepath, ra_col='right_ascension',
                             dec_col='declination', z_col='phz_median',
                             columns=None):
        """
        Load a galaxy catalogue from FITS file.

        Parameters
        ----------
        filepath : str
            Path to the FITS file
        ra_col : str
            Name of the RA column
        dec_col : str
            Name of the declination column
        z_col : str
            Name of the redshift column
        columns : list, optional
            List of column names to load. If None, loads all columns.
            For memory efficiency with large catalogues, specify only needed columns.

        Returns
        -------
        astropy.table.Table
            The loaded catalogue
        """
        print(f"Loading galaxy catalogue from {filepath}")

        with fits.open(filepath) as hdul:
            if columns is not None:
                # Only load specified columns for memory efficiency
                print(f"  Loading only {len(columns)} columns to save memory")
                galaxies = Table()
                for col in columns:
                    galaxies[col] = hdul[1].data[col]
            else:
                galaxies = Table(hdul[1].data)

        # Standardize column names
        galaxies['ra'] = galaxies[ra_col]
        galaxies['dec'] = galaxies[dec_col]
        galaxies['redshift'] = galaxies[z_col]

        print(f"Loaded {len(galaxies)} galaxies")
        return galaxies

    def load_random_catalogue(self, filepath, ra_col='right_ascension',
                             dec_col='declination', z_col='z', columns=None):
        """
        Load a random catalogue from FITS file.

        Parameters
        ----------
        filepath : str
            Path to the FITS file
        ra_col : str
            Name of the RA column
        dec_col : str
            Name of the declination column
        z_col : str
            Name of the redshift column
        columns : list, optional
            List of column names to load. If None, loads all columns.

        Returns
        -------
        astropy.table.Table
            The loaded catalogue
        """
        print(f"Loading random catalogue from {filepath}")

        with fits.open(filepath) as hdul:
            if columns is not None:
                # Only load specified columns for memory efficiency
                print(f"  Loading only {len(columns)} columns to save memory")
                randoms = Table()
                for col in columns:
                    randoms[col] = hdul[1].data[col]
            else:
                randoms = Table(hdul[1].data)

        # Standardize column names
        randoms['ra'] = randoms[ra_col]
        randoms['dec'] = randoms[dec_col]
        randoms['redshift'] = randoms[z_col]

        print(f"Loaded {len(randoms)} randoms")
        return randoms

    def add_cartesian_coordinates(self, catalogue):
        """
        Add Cartesian coordinates (x, y, z) to a catalogue.

        Parameters
        ----------
        catalogue : astropy.table.Table
            Catalogue with 'ra', 'dec', and 'redshift' columns

        Returns
        -------
        astropy.table.Table
            Catalogue with added x, y, z columns
        """
        # Handle NaN redshifts
        valid_z = ~np.isnan(catalogue['redshift'])
        d = np.zeros(len(catalogue))
        d[valid_z] = self.cosmo.comoving_distance(catalogue['redshift'][valid_z]).value

        catalogue['comoving_distance'] = d

        # Convert to radians
        ra_rad = np.deg2rad(catalogue['ra'])
        dec_rad = np.deg2rad(catalogue['dec'])

        # Calculate Cartesian coordinates
        catalogue['x'] = d * np.cos(dec_rad) * np.cos(ra_rad)
        catalogue['y'] = d * np.cos(dec_rad) * np.sin(ra_rad)
        catalogue['z'] = d * np.sin(dec_rad)

        return catalogue

    def add_galaxy_shapes(self, galaxies, ellipticity_col='ellipticity',
                         pa_col='position_angle', shape_error=0.3):
        """
        Add shape information (e1, e2, weights) to galaxy catalogue.

        Parameters
        ----------
        galaxies : astropy.table.Table
            Galaxy catalogue
        ellipticity_col : str
            Name of the ellipticity column
        pa_col : str
            Name of the position angle column (in degrees)
        shape_error : float
            Assumed shape measurement error for weights

        Returns
        -------
        astropy.table.Table
            Catalogue with added e1, e2, and w columns
        """
        # Convert sextractor ellipticity to weak lensing convention
        e_sex = galaxies[ellipticity_col]
        e = -(1 - (1 - e_sex)**2) / (1 + (1 - e_sex)**2)

        # Position angle in radians
        pa_rad = np.deg2rad(galaxies[pa_col])

        # Complex ellipticity components
        galaxies['e1'] = e * np.cos(2 * pa_rad)
        galaxies['e2'] = e * np.sin(2 * pa_rad)

        # Weights (inverse variance)
        galaxies['w'] = np.ones(len(galaxies)) / shape_error**2

        return galaxies

    def create_treecorr_catalogue(self, catalogue, npatch=50,
                                  use_shapes=False, use_weights=False,
                                  patch_centers=None, mode='3d'):
        """
        Create a TreeCorr catalogue object.

        Parameters
        ----------
        catalogue : astropy.table.Table
            Input catalogue with appropriate columns
        npatch : int
            Number of patches for jackknife resampling
        use_shapes : bool
            Whether to include shape information (g1, g2)
        use_weights : bool
            Whether to include weights
        patch_centers : array-like, optional
            Pre-computed patch centers to use
        mode : str
            '3d' for Cartesian coordinates (x,y,z) or '2d' for angles (ra,dec)

        Returns
        -------
        treecorr.Catalog
            TreeCorr catalogue object
        """
        kwargs = {}

        if mode == '3d':
            kwargs['x'] = catalogue['x']
            kwargs['y'] = catalogue['y']
            kwargs['z'] = catalogue['z']
        elif mode == '2d':
            kwargs['ra'] = catalogue['ra']
            kwargs['dec'] = catalogue['dec']
            kwargs['ra_units'] = 'deg'
            kwargs['dec_units'] = 'deg'
        else:
            raise ValueError("mode must be '3d' or '2d'")

        if use_shapes:
            kwargs['g1'] = catalogue['e1']
            kwargs['g2'] = catalogue['e2']

        if use_weights and 'w' in catalogue.colnames:
            kwargs['w'] = catalogue['w']

        if patch_centers is not None:
            kwargs['patch_centers'] = patch_centers
        else:
            kwargs['npatch'] = npatch

        return treecorr.Catalog(**kwargs)
