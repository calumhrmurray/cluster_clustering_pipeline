#!/usr/bin/env python3
"""Generate a uniform random catalogue within a HEALPix footprint mask.

Pixels are sampled directly (optionally weighted), then points are drawn
uniformly inside each pixel for speed.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import warnings

import healpy as hp
from healpy.utils.deprecation import HealpyDeprecationWarning
import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq
from astropy.io import fits


DEFAULT_MASK = Path(
    "/sps/euclid/Users/cmurray/euclid_clusters_tr1_validation/systematic_maps/maps/TR1/tr1_south_mask_nside2048.fits"
)


def parse_args() -> argparse.Namespace:
    script_root = Path(__file__).resolve().parents[1]
    default_output = script_root / "outputs" / "random_catalogues"
    parser = argparse.ArgumentParser(
        description="Create a random catalogue uniformly inside a HEALPix footprint."
    )
    parser.add_argument("--mask", type=Path, default=DEFAULT_MASK, help="HEALPix mask FITS path.")
    parser.add_argument("--n-randoms", type=int, required=True, help="Number of random points.")
    parser.add_argument(
        "--min-mask",
        type=float,
        default=0.8,
        help="Minimum mask value (exclusive) to keep.",
    )
    parser.add_argument(
        "--use-mask-weights",
        action="store_true",
        help="Sample pixels with probability proportional to the mask value.",
    )
    parser.add_argument("--seed", type=int, default=12345, help="Random seed.")
    parser.add_argument("--chunk-size", type=int, default=1_000_000,
                        help="Number of random points per chunk.")
    parser.add_argument("--ra-name", type=str, default="right_ascension", help="RA column name.")
    parser.add_argument("--dec-name", type=str, default="declination", help="Dec column name.")
    parser.add_argument("--output", type=Path, default=default_output,
                        help="Output directory or parquet file path.")
    return parser.parse_args()


def _read_mask(path: Path) -> tuple[np.ndarray, bool]:
    """Read a HEALPix mask and return (map, nest_flag)."""
    with fits.open(path, memmap=False) as hdul:
        header = None
        for hdu in hdul:
            if getattr(hdu, "data", None) is None:
                continue
            header = hdu.header
            break
        ordering = (header.get("ORDERING") or header.get("ORDER") or "RING").upper()
        nest = ordering.startswith("NEST")

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=HealpyDeprecationWarning)
        mask = hp.read_map(path, nest=nest, verbose=False)
    return np.asarray(mask, dtype=float), nest


def _draw_vectors(
    rng: np.random.Generator,
    cos_max: float,
    u: np.ndarray,
    v: np.ndarray,
    c: np.ndarray,
) -> np.ndarray:
    """Draw unit vectors uniformly within a cap around each center vector."""
    n = len(c)
    u_rand = rng.random(n)
    v_rand = rng.random(n)
    cos_a = 1.0 - u_rand * (1.0 - cos_max)
    sin_a = np.sqrt(1.0 - cos_a * cos_a)
    phi = 2.0 * np.pi * v_rand
    cos_phi = np.cos(phi)
    sin_phi = np.sin(phi)
    return (
        u * (sin_a * cos_phi)[:, None]
        + v * (sin_a * sin_phi)[:, None]
        + c * cos_a[:, None]
    )


def _sample_in_pixels(
    nside: int,
    pixels: np.ndarray,
    nest: bool,
    rng: np.random.Generator,
    max_iter: int = 100,
) -> tuple[np.ndarray, np.ndarray]:
    """Sample RA/Dec uniformly within each HEALPix pixel."""
    max_rad = hp.max_pixrad(nside)
    cos_max = np.cos(max_rad)

    cx, cy, cz = hp.pix2vec(nside, pixels, nest=nest)
    c = np.stack([cx, cy, cz], axis=1)

    ref = np.zeros_like(c)
    ref[:, 2] = 1.0
    near_pole = np.abs(cz) > 0.99
    ref[near_pole] = np.array([0.0, 1.0, 0.0])

    u = np.cross(ref, c)
    u_norm = np.linalg.norm(u, axis=1)
    zero = u_norm == 0
    if np.any(zero):
        u[zero] = np.array([1.0, 0.0, 0.0])
        u_norm[zero] = 1.0
    u /= u_norm[:, None]
    v = np.cross(c, u)

    vec = _draw_vectors(rng, cos_max, u, v, c)
    pix_check = hp.vec2pix(nside, vec[:, 0], vec[:, 1], vec[:, 2], nest=nest)
    bad = pix_check != pixels

    it = 0
    while np.any(bad):
        it += 1
        if it > max_iter:
            raise RuntimeError(
                f"Sampling did not converge after {max_iter} iterations; "
                "reduce --chunk-size or increase max_iter."
            )
        vec_bad = _draw_vectors(rng, cos_max, u[bad], v[bad], c[bad])
        vec[bad] = vec_bad
        pix_check = hp.vec2pix(nside, vec[bad, 0], vec[bad, 1], vec[bad, 2], nest=nest)
        bad_idx = np.flatnonzero(bad)
        bad[bad_idx] = pix_check != pixels[bad_idx]

    theta, phi = hp.vec2ang(vec)
    ra = np.degrees(phi)
    dec = 90.0 - np.degrees(theta)
    return ra, dec


def main() -> int:
    args = parse_args()
    mask_path = args.mask.resolve()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True) if output.suffix != ".parquet" else output.parent.mkdir(parents=True, exist_ok=True)

    mask_map, nest = _read_mask(mask_path)
    nside = hp.npix2nside(len(mask_map))

    valid = np.isfinite(mask_map) & (mask_map != hp.UNSEEN) & (mask_map > args.min_mask)
    if not np.any(valid):
        raise ValueError("Mask has no valid pixels after thresholding.")

    area_frac = valid.sum() / len(mask_map)
    print(f"Mask: {mask_path}")
    print(f"NSIDE: {nside} | valid frac: {area_frac:.3f}")

    valid_pixels = np.flatnonzero(valid)
    if args.use_mask_weights:
        weights = np.clip(mask_map[valid_pixels], 0.0, 1.0)
        weights_sum = weights.sum()
        if weights_sum <= 0:
            raise ValueError("Mask weights are all zero; cannot weight-sample.")
        weights /= weights_sum
    else:
        weights = None

    out_path = output
    if out_path.suffix != ".parquet":
        out_path = output / f"randoms_{nside}_n{args.n_randoms}.parquet"

    schema = pa.schema(
        [
            (args.ra_name, pa.float64()),
            (args.dec_name, pa.float64()),
        ]
    )
    writer = pq.ParquetWriter(out_path, schema)

    rng = np.random.default_rng(args.seed)
    remaining = args.n_randoms
    generated = 0

    try:
        while remaining > 0:
            n_chunk = min(args.chunk_size, remaining)
            if weights is None:
                idx = rng.integers(0, valid_pixels.size, size=n_chunk)
                pix = valid_pixels[idx]
            else:
                pix = rng.choice(valid_pixels, size=n_chunk, replace=True, p=weights)

            ra, dec = _sample_in_pixels(nside, pix, nest, rng)
            table = pa.Table.from_arrays([ra, dec], schema=schema)
            writer.write_table(table)

            remaining -= n_chunk
            generated += n_chunk
            print(f"Generated {generated}/{args.n_randoms} randoms...")
    finally:
        writer.close()

    print(f"Wrote random catalogue: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
