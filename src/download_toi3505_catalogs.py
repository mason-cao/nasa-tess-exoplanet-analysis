"""Download the public Gaia DR3 and TIC neighborhood of TOI-3505.01.

The 2.5-arcminute radius covers the ground-based nearby-star check and most of
the 15x15-pixel TESS cutout.  The exact queries are saved beside the returned
tables so the catalog census can be repeated later.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from datetime import datetime, timezone
from pathlib import Path

import astropy.units as u
import certifi
from astropy.coordinates import SkyCoord
from astroquery.gaia import Gaia
from astroquery.mast import Catalogs


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT_DIR = ROOT / "data" / "catalogs" / "toi3505"
RA_DEG = 297.043476
DEC_DEG = 18.698914
RADIUS_ARCMIN = 2.5

os.environ.setdefault("SSL_CERT_FILE", certifi.where())


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    return parser.parse_args()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(1024 * 1024):
            digest.update(block)
    return digest.hexdigest()


def gaia_query() -> str:
    radius_deg = RADIUS_ARCMIN / 60.0
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
    non_single_star,
    DISTANCE(
        POINT('ICRS', ra, dec),
        POINT('ICRS', {RA_DEG:.6f}, {DEC_DEG:.6f})
    ) * 3600.0 AS separation_arcsec
FROM gaiadr3.gaia_source
WHERE 1 = CONTAINS(
    POINT('ICRS', ra, dec),
    CIRCLE('ICRS', {RA_DEG:.6f}, {DEC_DEG:.6f}, {radius_deg:.10f})
)
ORDER BY separation_arcsec
"""


def main() -> None:
    args = parse_args()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    query = gaia_query()
    (output_dir / "gaia_dr3_query.adql").write_text(query, encoding="utf-8")
    print("Querying the Gaia DR3 archive...")
    gaia = Gaia.launch_job_async(query, dump_to_file=False).get_results()
    gaia_path = output_dir / "gaia_dr3_2p5arcmin.csv"
    gaia.write(gaia_path, format="ascii.csv", overwrite=True)

    print("Querying the TESS Input Catalog v8...")
    center = SkyCoord(RA_DEG * u.deg, DEC_DEG * u.deg, frame="icrs")
    tic = Catalogs.query_region(
        center,
        radius=RADIUS_ARCMIN * u.arcmin,
        catalog="TIC",
        version=8,
        pagesize=5000,
    )
    useful_columns = [
        name
        for name in (
            "ID",
            "ra",
            "dec",
            "Tmag",
            "e_Tmag",
            "GAIA",
            "GAIAmag",
            "gaiabp",
            "gaiarp",
            "pmRA",
            "pmDEC",
            "plx",
            "Teff",
            "rad",
            "mass",
            "d",
            "numcont",
            "contratio",
            "disposition",
            "duplicate_id",
            "dstArcSec",
        )
        if name in tic.colnames
    ]
    tic_path = output_dir / "tic_v8_2p5arcmin.csv"
    tic[useful_columns].write(tic_path, format="ascii.csv", overwrite=True)
    (output_dir / "tic_query.txt").write_text(
        "astroquery.mast.Catalogs.query_region\n"
        f"center_icrs_deg = {RA_DEG:.6f}, {DEC_DEG:.6f}\n"
        f"radius_arcmin = {RADIUS_ARCMIN}\n"
        "catalog = TIC\nversion = 8\n",
        encoding="utf-8",
    )

    target_rows = tic[[str(value) == "390988385" for value in tic["ID"]]]
    if len(target_rows) != 1:
        raise RuntimeError(f"Expected one TIC 390988385 row; found {len(target_rows)}")
    target = target_rows[0]
    summary = {
        "target": "TOI-3505.01",
        "tic_id": 390988385,
        "center_icrs_degrees": [RA_DEG, DEC_DEG],
        "radius_arcmin": RADIUS_ARCMIN,
        "retrieved_utc": datetime.now(timezone.utc).isoformat(),
        "gaia_dr3_rows": len(gaia),
        "tic_v8_rows": len(tic),
        "target_tmag": float(target["Tmag"]),
        "target_tic_contamination_ratio": float(target["contratio"]),
        "files": {
            gaia_path.name: sha256_file(gaia_path),
            tic_path.name: sha256_file(tic_path),
        },
        "sources": {
            "gaia": "ESA Gaia Archive, gaiadr3.gaia_source",
            "tic": "MAST TESS Input Catalog version 8 cone service",
        },
    }
    (output_dir / "catalog_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(f"Saved {len(gaia)} Gaia sources and {len(tic)} TIC rows to {output_dir}")


if __name__ == "__main__":
    main()
