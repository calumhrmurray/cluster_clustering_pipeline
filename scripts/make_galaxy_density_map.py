#!/usr/bin/env python3
"""Build a HEALPix density map from a parquet catalogue and plot a TAN projection.

The TAN projection reuses the cluster validation pipeline's systematics map
visualisation code in:
  /sps/euclid/Users/cmurray/euclid_clusters_tr1_validation/cluster_validation/pipeline
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path
import sys

os.environ.setdefault("MPLBACKEND", "Agg")

import healpy as hp
import numpy as np
import pyarrow.parquet as pq
from astropy.io import fits


DEFAULT_CATALOGUE = Path(
    "/sps/euclid/OU-LE3/CL/COMB/TR1/shear_catalogues/clustering_catalogue_130125_calum.pq"
)
DEFAULT_VALIDATION_ROOT = Path(
    "/sps/euclid/Users/cmurray/euclid_clusters_tr1_validation"
)


def parse_args() -> argparse.Namespace:
    script_root = Path(__file__).resolve().parents[1]
    default_output = script_root / "outputs" / "galaxy_density_maps"
    parser = argparse.ArgumentParser(description="Create a HEALPix density map from a parquet catalogue.")
    parser.add_argument("--input", type=Path, default=DEFAULT_CATALOGUE, help="Path to the parquet catalogue.")
    parser.add_argument("--output-dir", type=Path, default=default_output, help="Directory for map + plots.")
    parser.add_argument("--nside", type=int, default=1024, help="HEALPix NSIDE for the density map.")
    parser.add_argument("--batch-size", type=int, default=1_000_000, help="Rows per parquet batch read.")
    parser.add_argument("--ra-col", type=str, default=None, help="RA column name override.")
    parser.add_argument("--dec-col", type=str, default=None, help="Dec column name override.")
    parser.add_argument("--name", type=str, default=None, help="Base name for outputs.")
    parser.add_argument("--center-ra", type=float, default=55.0, help="TAN projection center RA (deg).")
    parser.add_argument("--center-dec", type=float, default=-46.0, help="TAN projection center Dec (deg).")
    parser.add_argument("--npix-x", type=int, default=1000, help="TAN projection x pixels.")
    parser.add_argument("--npix-y", type=int, default=900, help="TAN projection y pixels.")
    parser.add_argument("--pixel-scale", type=float, default=0.05, help="TAN projection pixel scale (deg).")
    parser.add_argument("--validation-root", type=Path, default=DEFAULT_VALIDATION_ROOT,
                        help="Path to the cluster validation repo root.")
    return parser.parse_args()


def detect_coord_columns(path: Path, ra_override: str | None, dec_override: str | None) -> tuple[str, str]:
    schema = pq.ParquetFile(path).schema_arrow
    names = set(schema.names)

    if ra_override or dec_override:
        if not (ra_override and dec_override):
            raise ValueError("Provide both --ra-col and --dec-col when overriding.")
        if ra_override not in names or dec_override not in names:
            raise ValueError(f"Columns not found in parquet: {ra_override}, {dec_override}")
        return ra_override, dec_override

    candidates = [
        ("she_metacal_ra", "she_metacal_dec"),
        ("she_lensmc_ra", "she_lensmc_dec"),
        ("right_ascension", "declination"),
        ("ra", "dec"),
        ("RA", "DEC"),
    ]
    for ra_col, dec_col in candidates:
        if ra_col in names and dec_col in names:
            return ra_col, dec_col

    raise ValueError(
        "Could not auto-detect RA/Dec columns. "
        "Use --ra-col and --dec-col to specify explicitly."
    )


def build_density_map(
    path: Path,
    ra_col: str,
    dec_col: str,
    nside: int,
    batch_size: int,
) -> tuple[np.ndarray, int, int]:
    pf = pq.ParquetFile(path)
    npix = hp.nside2npix(nside)
    counts = np.zeros(npix, dtype=np.int64)

    total_rows = 0
    used_rows = 0

    for batch in pf.iter_batches(columns=[ra_col, dec_col], batch_size=batch_size):
        ra = batch.column(ra_col).to_numpy(zero_copy_only=False)
        dec = batch.column(dec_col).to_numpy(zero_copy_only=False)

        total_rows += len(ra)
        mask = np.isfinite(ra) & np.isfinite(dec)
        if not np.any(mask):
            continue

        ra = np.mod(ra[mask].astype(float, copy=False), 360.0)
        dec = np.clip(dec[mask].astype(float, copy=False), -90.0, 90.0)
        used_rows += ra.size

        theta = np.radians(90.0 - dec)
        phi = np.radians(ra)
        pixels = hp.ang2pix(nside, theta, phi, nest=False)
        counts += np.bincount(pixels, minlength=npix)

    return counts.astype(float), total_rows, used_rows


def run_tan_projection(
    map_path: Path,
    output_dir: Path,
    name: str,
    center_ra: float,
    center_dec: float,
    npix_x: int,
    npix_y: int,
    pixel_scale: float,
    validation_root: Path,
) -> None:
    if not validation_root.exists():
        raise FileNotFoundError(f"Validation repo not found: {validation_root}")

    sys.path.insert(0, str(validation_root))
    from cluster_validation.pipeline.config import SystematicConfig  # noqa: E402
    from cluster_validation.pipeline.io import SystematicData  # noqa: E402
    from cluster_validation.pipeline.systematics.tan_projection import run_tan_projections  # noqa: E402

    cfg = SystematicConfig(type="density", name=name, path=map_path)
    cfg.colorbar_label = "Galaxy surface density (counts per pixel)"
    systematic = SystematicData(config=cfg, hdu=fits.open(map_path))
    try:
        run_tan_projections(
            [systematic],
            output_dir,
            center_ra=center_ra,
            center_dec=center_dec,
            npix_x=npix_x,
            npix_y=npix_y,
            pixel_scale=pixel_scale,
        )
    finally:
        systematic.hdu.close()


def main() -> int:
    args = parse_args()
    input_path = args.input.resolve()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    ra_col, dec_col = detect_coord_columns(input_path, args.ra_col, args.dec_col)
    name = args.name or input_path.stem

    print(f"Input: {input_path}")
    print(f"Columns: RA={ra_col}, Dec={dec_col}")
    print(f"NSIDE: {args.nside}")

    density_map, total_rows, used_rows = build_density_map(
        input_path,
        ra_col,
        dec_col,
        args.nside,
        args.batch_size,
    )

    map_dir = output_dir / "healpix_maps"
    map_dir.mkdir(parents=True, exist_ok=True)
    map_path = map_dir / f"{name}_density_nside{args.nside}.fits"
    hp.write_map(map_path, density_map, nest=False, dtype=float, overwrite=True)

    print(f"Rows processed: {used_rows}/{total_rows}")
    print(f"Saved HEALPix map: {map_path}")

    tan_dir = output_dir / "tan_projections"
    tan_dir.mkdir(parents=True, exist_ok=True)
    run_tan_projection(
        map_path,
        tan_dir,
        name,
        center_ra=args.center_ra,
        center_dec=args.center_dec,
        npix_x=args.npix_x,
        npix_y=args.npix_y,
        pixel_scale=args.pixel_scale,
        validation_root=args.validation_root,
    )
    print(f"TAN projection saved under: {tan_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
