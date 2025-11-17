"""
Systematic weight loading and assignment from HEALPix maps.

This module provides functions to:
1. Load pre-computed systematic weights from FITS files
2. Assign weights to catalogue objects based on sky position
3. Combine systematic weights with other weights

Systematic weights correct for spatial variations in survey properties
(Galactic extinction, stellar density, survey depth, etc.) that can bias
clustering measurements.
"""

import numpy as np
import healpy as hp
from astropy.io import fits
from astropy.table import Table
from pathlib import Path


def load_weight_map(weight_file_path):
    """
    Load systematic weights from HEALPix FITS file.

    Parameters
    ----------
    weight_file_path : str or Path
        Path to FITS file with PIXEL and WEIGHT columns

    Returns
    -------
    weights : np.ndarray
        Weight values for each HEALPix pixel
    nside : int
        HEALPix NSIDE parameter
    nest : bool
        True if NESTED ordering, False if RING

    Examples
    --------
    >>> weights, nside, nest = load_weight_map('AMICO_SNR8_weights.fits')
    >>> print(f"Loaded {len(weights)} pixels at nside={nside}")
    """
    weight_file_path = Path(weight_file_path)

    if not weight_file_path.exists():
        raise FileNotFoundError(f"Weight file not found: {weight_file_path}")

    print(f"Loading systematic weight map from {weight_file_path}")

    with fits.open(weight_file_path) as hdul:
        data = hdul[1].data
        weights = np.array(data['WEIGHT'], dtype=float)

        # Determine nside from number of pixels
        npix = len(weights)
        nside = hp.npix2nside(npix)

        # Check header for ordering (default: RING)
        header = hdul[1].header
        ordering = header.get('ORDERING', 'RING').upper()
        nest = (ordering == 'NEST')

    print(f"  NSIDE={nside}, npix={npix}, ordering={'NEST' if nest else 'RING'}")
    print(f"  Weight statistics: min={weights.min():.3f}, max={weights.max():.3f}, mean={weights.mean():.3f}")

    return weights, nside, nest


def assign_weights_to_catalogue(catalogue, weight_map_info,
                                  ra_col='ra', dec_col='dec',
                                  output_col='systematic_weight'):
    """
    Assign systematic weights to catalogue objects based on position.

    Uses HEALPix to determine which pixel each object falls in, then
    assigns the corresponding weight.

    Parameters
    ----------
    catalogue : astropy.table.Table
        Catalogue with RA/Dec columns
    weight_map_info : tuple
        (weights, nside, nest) from load_weight_map()
    ra_col : str
        Column name for right ascension (degrees)
    dec_col : str
        Column name for declination (degrees)
    output_col : str
        Name for output weight column

    Returns
    -------
    catalogue : astropy.table.Table
        Input catalogue with added weight column

    Examples
    --------
    >>> weight_info = load_weight_map('weights.fits')
    >>> catalogue = assign_weights_to_catalogue(galaxies, weight_info)
    >>> print(f"Mean weight: {catalogue['systematic_weight'].mean():.3f}")
    """
    weights, nside, nest = weight_map_info

    # Extract coordinates
    ra = np.array(catalogue[ra_col], dtype=float)
    dec = np.array(catalogue[dec_col], dtype=float)

    # Initialize output weights
    object_weights = np.ones(len(catalogue), dtype=float)

    # Only process valid coordinates
    valid = np.isfinite(ra) & np.isfinite(dec)
    n_valid = valid.sum()

    if n_valid == 0:
        print("WARNING: No valid coordinates found!")
        catalogue[output_col] = object_weights
        return catalogue

    # Convert to HEALPix pixel indices
    # theta = co-latitude (0 at North pole, 180 at South pole)
    # phi = longitude (0-360 degrees)
    theta = np.radians(90.0 - dec[valid])  # Co-latitude
    phi = np.radians(ra[valid] % 360.0)     # Longitude [0, 360)
    pixels = hp.ang2pix(nside, theta, phi, nest=nest)

    # Safety: clip to valid pixel range
    pixels = np.clip(pixels, 0, len(weights) - 1)

    # Look up weights for each pixel
    pixel_weights = weights[pixels]

    # Handle any NaN/inf in weight map
    pixel_weights[~np.isfinite(pixel_weights)] = 1.0

    # Assign to output array
    object_weights[valid] = pixel_weights

    # Add column to catalogue
    catalogue[output_col] = object_weights

    print(f"  Assigned systematic weights to {len(catalogue)} objects ({n_valid} valid)")
    print(f"  Weight stats: min={object_weights[valid].min():.3f}, "
          f"max={object_weights[valid].max():.3f}, mean={object_weights[valid].mean():.3f}")

    return catalogue


def combine_with_existing_weights(catalogue, systematic_col='systematic_weight',
                                    existing_col='w', output_col='total_weight'):
    """
    Multiply systematic weights with existing weights.

    This is useful when you have other weights (e.g., lensing shape weights)
    that need to be combined with systematic weights.

    Parameters
    ----------
    catalogue : astropy.table.Table
        Catalogue with weight columns
    systematic_col : str
        Column name for systematic weights
    existing_col : str
        Column name for existing weights (e.g., lensing weights)
    output_col : str
        Column name for combined output

    Returns
    -------
    catalogue : astropy.table.Table
        Catalogue with combined weight column
    """
    if existing_col not in catalogue.colnames:
        print(f"  No existing '{existing_col}' column found, using systematic weights only")
        catalogue[output_col] = catalogue[systematic_col]
    else:
        existing = np.array(catalogue[existing_col])
        systematic = np.array(catalogue[systematic_col])
        catalogue[output_col] = existing * systematic
        print(f"  Combined {existing_col} and {systematic_col} → {output_col}")
        print(f"  Mean combined weight: {catalogue[output_col].mean():.3f}")

    return catalogue


def apply_systematic_weights(catalogue, weight_file,
                              ra_col='ra', dec_col='dec',
                              existing_weight_col=None,
                              output_col='w'):
    """
    Convenience function: load weight map and assign to catalogue in one step.

    This is the main function you'll typically use. It handles loading the
    weight map, assigning weights to objects, and optionally combining with
    existing weights.

    Parameters
    ----------
    catalogue : astropy.table.Table
        Catalogue to add weights to
    weight_file : str or Path
        Path to systematic weight FITS file
    ra_col, dec_col : str
        Coordinate column names
    existing_weight_col : str, optional
        If specified, multiply with these existing weights
    output_col : str
        Name for final weight column (default: 'w' for TreeCorr)

    Returns
    -------
    catalogue : astropy.table.Table
        Catalogue with assigned weights

    Examples
    --------
    >>> # Simple case: just add systematic weights
    >>> galaxies = apply_systematic_weights(
    ...     galaxies, 'RR2_GALAXIES_weights.fits', output_col='w')
    >>>
    >>> # Combine with existing weights
    >>> galaxies = apply_systematic_weights(
    ...     galaxies, 'RR2_GALAXIES_weights.fits',
    ...     existing_weight_col='lensing_weight', output_col='w')
    """
    # Load weight map
    weight_info = load_weight_map(weight_file)

    # Assign to catalogue
    catalogue = assign_weights_to_catalogue(
        catalogue, weight_info, ra_col=ra_col, dec_col=dec_col,
        output_col='_systematic_temp')

    # Combine with existing weights if requested
    if existing_weight_col is not None and existing_weight_col in catalogue.colnames:
        catalogue = combine_with_existing_weights(
            catalogue, systematic_col='_systematic_temp',
            existing_col=existing_weight_col, output_col=output_col)
        # Clean up temp column
        catalogue.remove_column('_systematic_temp')
    else:
        # Just rename to output
        catalogue.rename_column('_systematic_temp', output_col)

    return catalogue
