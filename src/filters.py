"""
Filtering utilities for selecting subsets of clusters and galaxies.
"""

import numpy as np
import healpy as hp
from astropy.io import fits
from pathlib import Path


class CatalogueFilter:
    """Class for applying filters to astronomical catalogues."""

    @staticmethod
    def _load_healpix_map(map_file, nest=False, nside=None, value_field=None):
        map_path = Path(map_file)
        if value_field:
            with fits.open(map_path) as hdul:
                data = np.asarray(hdul[1].data[value_field], dtype=float).reshape(-1)
                header = hdul[1].header
                map_nside = int(header.get("NSIDE", nside or 0))
        else:
            data = hp.read_map(map_path, nest=nest, dtype=np.float64, verbose=False)
            map_nside = hp.npix2nside(len(data))

        if nside and nside != map_nside:
            order_in = "NEST" if nest else "RING"
            order_out = "NEST" if nest else "RING"
            data = hp.ud_grade(data, nside_out=nside, order_in=order_in, order_out=order_out)
            map_nside = nside

        return data, map_nside

    @staticmethod
    def apply_healpix_mask(catalogue, map_file, min_value=None, max_value=None,
                           nest=False, nside=None, value_field=None, label=None):
        """
        Filter a catalogue using a HEALPix mask map.

        Parameters
        ----------
        catalogue : astropy.table.Table
            Input catalogue with 'ra' and 'dec' columns.
        map_file : str or Path
            HEALPix FITS map file.
        min_value : float, optional
            Minimum mask value to keep (e.g. 0.8).
        max_value : float, optional
            Maximum mask value to keep.
        nest : bool, optional
            Whether the map uses NESTED ordering.
        nside : int, optional
            Override map NSIDE (map is upgraded/downgraded if needed).
        value_field : str, optional
            Column name if the FITS file is a table with values.
        label : str, optional
            Label for logging.
        """
        if 'ra' not in catalogue.colnames or 'dec' not in catalogue.colnames:
            print("Warning: Catalogue missing 'ra'/'dec' columns; skipping mask filter.")
            return catalogue

        map_values, map_nside = CatalogueFilter._load_healpix_map(
            map_file, nest=nest, nside=nside, value_field=value_field
        )

        ra = np.asarray(catalogue['ra'], dtype=float)
        dec = np.asarray(catalogue['dec'], dtype=float)
        valid = np.isfinite(ra) & np.isfinite(dec)

        keep = np.zeros(len(catalogue), dtype=bool)
        if not np.any(valid):
            print("Warning: No valid coordinates for mask filtering.")
            return catalogue[keep]

        theta = np.radians(90.0 - dec[valid])
        phi = np.radians(np.mod(ra[valid], 360.0))
        pix = hp.ang2pix(map_nside, theta, phi, nest=nest)

        values = map_values[pix]
        good = np.isfinite(values) & (values != hp.UNSEEN) & (values > -1e20)
        if min_value is not None:
            good &= values >= min_value
        if max_value is not None:
            good &= values <= max_value

        keep[valid] = good
        tag = f" ({label})" if label else ""
        print(f"Mask filter{tag}: {np.sum(keep)} / {len(catalogue)} objects pass")
        return catalogue[keep]

    @staticmethod
    def apply_filters(catalogue, filters):
        """
        Apply a set of filters to a catalogue.

        Parameters
        ----------
        catalogue : astropy.table.Table
            Input catalogue
        filters : dict
            Dictionary of filter specifications. Each key is a column name,
            and the value is either:
            - A tuple (min, max) for range filters
            - A list of allowed values for categorical filters
            - A callable that takes the column and returns a boolean mask

        Returns
        -------
        astropy.table.Table
            Filtered catalogue

        Examples
        --------
        >>> filters = {
        ...     'redshift': (0.2, 0.5),
        ...     'RICHNESS_CLUSTER_pzwav': (10, None),  # min only
        ...     'sersic_sersic_vis_index': lambda x: x > 2  # custom function
        ... }
        >>> filtered = CatalogueFilter.apply_filters(catalogue, filters)
        """
        mask = np.ones(len(catalogue), dtype=bool)

        for col_name, filter_spec in filters.items():
            if col_name not in catalogue.colnames:
                print(f"Warning: Column '{col_name}' not found in catalogue. Skipping filter.")
                continue

            col_data = catalogue[col_name]

            # Handle different filter types
            if callable(filter_spec):
                # Custom function filter
                col_mask = filter_spec(col_data)
            elif isinstance(filter_spec, (tuple, list)) and len(filter_spec) == 2:
                # Range filter (min, max)
                min_val, max_val = filter_spec
                col_mask = np.ones(len(catalogue), dtype=bool)

                if min_val is not None:
                    col_mask &= (col_data >= min_val)
                if max_val is not None:
                    col_mask &= (col_data < max_val)
            else:
                # List of allowed values
                col_mask = np.isin(col_data, filter_spec)

            # Handle NaN values
            col_mask &= ~np.isnan(col_data)

            mask &= col_mask
            print(f"Filter '{col_name}': {np.sum(col_mask)} / {len(catalogue)} objects pass")

        print(f"Total: {np.sum(mask)} / {len(catalogue)} objects pass all filters")
        return catalogue[mask]

    @staticmethod
    def create_redshift_bins(catalogue, z_col='redshift', z_bins=None):
        """
        Create redshift bin filters.

        Parameters
        ----------
        catalogue : astropy.table.Table
            Input catalogue
        z_col : str
            Name of the redshift column
        z_bins : array-like
            Bin edges for redshift. E.g., [0.0, 0.5, 1.0, 1.5]

        Returns
        -------
        list of dict
            List of filter dictionaries, one per bin
        """
        if z_bins is None:
            # Default bins
            z_bins = [0.0, 0.5, 1.0, 1.5, 2.0]

        filters_list = []
        for i in range(len(z_bins) - 1):
            z_min, z_max = z_bins[i], z_bins[i + 1]
            filters_list.append({z_col: (z_min, z_max)})

        return filters_list

    @staticmethod
    def create_property_bins(catalogue, prop_col, bins):
        """
        Create binned filters for any property.

        Parameters
        ----------
        catalogue : astropy.table.Table
            Input catalogue
        prop_col : str
            Name of the property column
        bins : array-like
            Bin edges

        Returns
        -------
        list of dict
            List of filter dictionaries, one per bin
        """
        filters_list = []
        for i in range(len(bins) - 1):
            min_val, max_val = bins[i], bins[i + 1]
            filters_list.append({prop_col: (min_val, max_val)})

        return filters_list

    @staticmethod
    def select_early_type(catalogue, sersic_col='sersic_sersic_vis_index',
                         threshold=2.0):
        """
        Select early-type galaxies using Sersic index.

        Parameters
        ----------
        catalogue : astropy.table.Table
            Galaxy catalogue
        sersic_col : str
            Name of the Sersic index column
        threshold : float
            Sersic index threshold (early-type if > threshold)

        Returns
        -------
        astropy.table.Table
            Early-type galaxy catalogue
        """
        filters = {sersic_col: lambda x: x > threshold}
        return CatalogueFilter.apply_filters(catalogue, filters)

    @staticmethod
    def select_late_type(catalogue, sersic_col='sersic_sersic_vis_index',
                        threshold=2.0):
        """
        Select late-type galaxies using Sersic index.

        Parameters
        ----------
        catalogue : astropy.table.Table
            Galaxy catalogue
        sersic_col : str
            Name of the Sersic index column
        threshold : float
            Sersic index threshold (late-type if < threshold)

        Returns
        -------
        astropy.table.Table
            Late-type galaxy catalogue
        """
        filters = {sersic_col: lambda x: x < threshold}
        return CatalogueFilter.apply_filters(catalogue, filters)

    @staticmethod
    def create_lens_source_split(galaxy_catalogue, z_lens_range, z_source_range,
                                z_col='redshift'):
        """
        Create lens and source catalogues for weak lensing analysis.

        Parameters
        ----------
        galaxy_catalogue : astropy.table.Table
            Full galaxy catalogue
        z_lens_range : tuple
            (z_min, z_max) for lens selection
        z_source_range : tuple
            (z_min, z_max) for source selection
        z_col : str
            Name of the redshift column

        Returns
        -------
        tuple of astropy.table.Table
            (lenses, sources) catalogues
        """
        lens_filters = {z_col: z_lens_range}
        source_filters = {z_col: z_source_range}

        lenses = CatalogueFilter.apply_filters(galaxy_catalogue, lens_filters)
        sources = CatalogueFilter.apply_filters(galaxy_catalogue, source_filters)

        return lenses, sources

    @staticmethod
    def get_bin_label(filters):
        """
        Create a human-readable label from filter specification.

        Parameters
        ----------
        filters : dict
            Filter dictionary

        Returns
        -------
        str
            Label describing the filters
        """
        labels = []
        for col_name, filter_spec in filters.items():
            if isinstance(filter_spec, (tuple, list)) and len(filter_spec) == 2:
                min_val, max_val = filter_spec
                if min_val is not None and max_val is not None:
                    labels.append(f"{col_name}_{min_val:.2f}_{max_val:.2f}")
                elif min_val is not None:
                    labels.append(f"{col_name}_gt_{min_val:.2f}")
                elif max_val is not None:
                    labels.append(f"{col_name}_lt_{max_val:.2f}")
            else:
                labels.append(f"{col_name}_custom")

        return "_".join(labels) if labels else "all"
