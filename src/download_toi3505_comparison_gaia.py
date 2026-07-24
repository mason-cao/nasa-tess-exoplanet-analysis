"""Download small Gaia DR3 neighborhoods around the measured ground stars.

The comparison stars lie outside the 2.5-arcminute target catalog used for the
TESS crowding check.  This script converts the first-frame AstroImageJ
positions to sky coordinates with the existing plate solution, then asks Gaia
for sources within 12 arcseconds of each measured position.  The small cones
are enough to identify each star and count cataloged blends inside the
25-pixel (about 9-arcsecond) ground aperture without downloading the full
camera field.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from datetime import datetime, timezone
from pathlib import Path

import certifi
import pandas as pd
from astropy.io import fits
from astropy.wcs import WCS
from astroquery.gaia import Gaia

from analyze_toi3505_photometry import COMPARISON_STARS, load_table


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_TABLE = (
    ROOT
    / "outputs"
    / "toi3505_aperture_check"
    / "TOI_3505.01_2022-07-22_R_measurements_25px_70-139_AIJ.tbl"
)
DEFAULT_WCS_IMAGE = (
    ROOT
    / "data"
    / "ground"
    / "toi3505"
    / "plate_solved"
    / "TOI_3505.01_50.000s_R-0001_wcs.fits"
)
DEFAULT_OUTPUT_DIR = ROOT / "data" / "catalogs" / "toi3505"
SEARCH_RADIUS_ARCSEC = 12.0
STAR_NAMES = ("T1", *COMPARISON_STARS)

os.environ.setdefault("SSL_CERT_FILE", certifi.where())


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--table", type=Path, default=DEFAULT_TABLE)
    parser.add_argument("--wcs-image", type=Path, default=DEFAULT_WCS_IMAGE)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    return parser.parse_args()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(1024 * 1024):
            digest.update(block)
    return digest.hexdigest()


def measured_sky_positions(table: pd.DataFrame, wcs: WCS) -> pd.DataFrame:
    """Convert first-frame FITS positions to ICRS coordinates."""
    first = table.iloc[0]
    rows: list[dict[str, float | str]] = []
    for star in STAR_NAMES:
        # AIJ FITS coordinates are one-indexed; Astropy array coordinates are
        # zero-indexed.
        x = float(first[f"X(FITS)_{star}"]) - 1.0
        y = float(first[f"Y(FITS)_{star}"]) - 1.0
        coordinate = wcs.pixel_to_world(x, y)
        rows.append(
            {
                "star": star,
                "x_zero_indexed": x,
                "y_zero_indexed": y,
                "wcs_ra_deg": float(coordinate.ra.deg),
                "wcs_dec_deg": float(coordinate.dec.deg),
            }
        )
    return pd.DataFrame(rows)


def gaia_query(positions: pd.DataFrame) -> str:
    radius_deg = SEARCH_RADIUS_ARCSEC / 3600.0
    circles = []
    for row in positions.itertuples(index=False):
        circles.append(
            "1 = CONTAINS(POINT('ICRS', ra, dec), "
            f"CIRCLE('ICRS', {row.wcs_ra_deg:.9f}, {row.wcs_dec_deg:.9f}, "
            f"{radius_deg:.10f}))"
        )
    where = "\n    OR ".join(circles)
    return f"""SELECT
    source_id,
    ra,
    dec,
    phot_g_mean_mag,
    phot_bp_mean_mag,
    phot_rp_mean_mag,
    bp_rp,
    ruwe,
    parallax,
    parallax_error,
    pmra,
    pmdec,
    phot_variable_flag,
    non_single_star
FROM gaiadr3.gaia_source
WHERE {where}
"""


def attach_nearest_measured_star(
    sources: pd.DataFrame, positions: pd.DataFrame
) -> pd.DataFrame:
    """Label every returned Gaia source by its nearest measured aperture."""
    from astropy import units as u
    from astropy.coordinates import SkyCoord

    measured = SkyCoord(
        positions["wcs_ra_deg"].to_numpy() * u.deg,
        positions["wcs_dec_deg"].to_numpy() * u.deg,
    )
    catalog = SkyCoord(
        sources["ra"].to_numpy() * u.deg,
        sources["dec"].to_numpy() * u.deg,
    )
    nearest, separation, _ = catalog.match_to_catalog_sky(measured)
    labeled = sources.copy()
    labeled.insert(0, "measured_star", positions.iloc[nearest]["star"].to_numpy())
    labeled.insert(1, "separation_from_measured_position_arcsec", separation.arcsec)
    return labeled.sort_values(
        ["measured_star", "separation_from_measured_position_arcsec"]
    ).reset_index(drop=True)


def main() -> None:
    args = parse_args()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    table = load_table(
        args.table.resolve(), expected_outer_radius=139, expected_source_radius=25
    )
    wcs = WCS(fits.getheader(args.wcs_image.resolve()))
    positions = measured_sky_positions(table, wcs)
    positions_path = output_dir / "ground_star_wcs_positions.csv"
    positions.to_csv(positions_path, index=False, float_format="%.10f")

    query = gaia_query(positions)
    query_path = output_dir / "ground_star_gaia_dr3_query.adql"
    query_path.write_text(query, encoding="utf-8")
    print("Querying Gaia DR3 around the 11 measured ground stars...")
    result = Gaia.launch_job_async(query, dump_to_file=False).get_results()
    raw = result.to_pandas()
    labeled = attach_nearest_measured_star(raw, positions)
    catalog_path = output_dir / "ground_star_gaia_dr3_12arcsec.csv"
    labeled.to_csv(catalog_path, index=False, float_format="%.10f")

    counts = labeled.groupby("measured_star").size().reindex(STAR_NAMES, fill_value=0)
    missing = [star for star, count in counts.items() if count == 0]
    summary = {
        "target": "TOI-3505.01",
        "catalog": "Gaia DR3 gaiadr3.gaia_source",
        "retrieved_utc": datetime.now(timezone.utc).isoformat(),
        "search_radius_arcsec": SEARCH_RADIUS_ARCSEC,
        "measured_positions": len(positions),
        "returned_sources": len(labeled),
        "sources_per_measured_position": {
            star: int(count) for star, count in counts.items()
        },
        "positions_without_a_gaia_source": missing,
        "plate_solution": str(args.wcs_image.resolve()),
        "files": {
            positions_path.name: sha256_file(positions_path),
            query_path.name: sha256_file(query_path),
            catalog_path.name: sha256_file(catalog_path),
        },
    }
    (output_dir / "ground_star_catalog_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(
        f"Saved {len(labeled)} Gaia sources near {len(positions)} measured stars "
        f"to {output_dir}"
    )


if __name__ == "__main__":
    main()
