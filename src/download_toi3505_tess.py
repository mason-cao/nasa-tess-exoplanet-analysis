"""Download and record the public TESS data for TOI-3505.01.

The four-sector analysis uses the QLP light curves as its common data set.
SPOC and TESS-SPOC light curves are also saved for Sectors 54 and 81 so that
the QLP result can be checked against another reduction of the same pixels.
Small TESScut images are downloaded for all four sectors.  The 2-minute SPOC
target-pixel files are saved for Sectors 54 and 81, where they exist.
The corresponding SPOC Data Validation report, mini-report, summary, XML, and
transit-data FITS products are also saved for those two sectors.  A later SPOC
multi-sector search combined the Sector 54 and 81 observations at 10-minute
cadence; its complete report set and the target-specific Exo.MAST records are
saved separately so the search range is not confused with actual coverage.

Run from the project directory with::

    .venv/bin/python src/download_toi3505_tess.py

The downloaded products are public MAST files and are intentionally ignored
by git. Compact search tables and a checksum manifest are written to
``outputs/toi3505_tess_download``.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import platform
import shutil
import ssl
import sys
import tempfile
import urllib.request
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

os.environ.setdefault(
    "MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "nasa-tess-matplotlib-cache")
)
os.environ.setdefault(
    "XDG_CACHE_HOME", str(Path(tempfile.gettempdir()) / "nasa-tess-download-cache")
)

import astropy
import certifi
import lightkurve as lk
import numpy as np
from astropy.io import fits
from astroquery.mast import Observations


ROOT = Path(__file__).resolve().parents[1]
TIC_ID = 390988385
TARGET = f"TIC {TIC_ID}"
SECTORS = (14, 41, 54, 81)
PIXEL_SECTORS = (54, 81)
DATA_VALIDATION_SUBGROUPS = ("DVR", "DVM", "DVS", "DVT")
MULTISECTOR_OBSID = 262995968
MULTISECTOR_SCOPE = "s0014-s0086"
MULTISECTOR_OBSERVATION_ID = (
    "tess2019199201929-s0014-s0086-0000000390988385"
)
MULTISECTOR_CONTRIBUTING_SECTORS = (54, 81)
CUTOUT_SIZE = (15, 15)
LIGHT_CURVE_AUTHORS = ("QLP", "SPOC", "TESS-SPOC")
DEFAULT_DATA_DIR = ROOT / "data" / "tess" / "toi3505"
DEFAULT_OUTPUT_DIR = ROOT / "outputs" / "toi3505_tess_download"
EXOMAST_BASE = f"https://exo.mast.stsci.edu/api/v0.1/dvdata/tess/{TIC_ID}"
EXOMAST_VETTING_URL = (
    "https://exo.mast.stsci.edu/api/v0.1/exoplanets/vetting/"
    f"?obs={MULTISECTOR_OBSID}"
)
DR122_RELEASE_NOTE_URL = (
    "https://archive.stsci.edu/missions/tess/doc/tess_drn/"
    "tess_multisector_14_86_drn122_v01.pdf"
)
DR122_TARGET_INFO_URL = (
    "https://archive.stsci.edu/missions/tess/catalogs/targetinfo/"
    "tess_multisector_14_86_drn122_targetinfo_v01.txt"
)

SEARCH_COLUMNS = (
    "sequence_number",
    "mission",
    "author",
    "exptime",
    "target_name",
    "productFilename",
    "dataURI",
    "size",
    "t_min",
    "t_max",
    "obsid",
    "distance",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=DEFAULT_DATA_DIR)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument(
        "--inventory-only",
        action="store_true",
        help="query MAST and save search tables without downloading FITS files",
    )
    parser.add_argument(
        "--skip-pixels",
        action="store_true",
        help="download light curves but skip TESScut and target-pixel files",
    )
    return parser.parse_args()


def plain_value(value: Any) -> str | int | float | bool | None:
    """Convert an Astropy table value into a CSV- and JSON-safe value."""
    if value is None or np.ma.is_masked(value):
        return None
    if isinstance(value, np.generic):
        value = value.item()
    if isinstance(value, bytes):
        return value.decode(errors="replace")
    if isinstance(value, (str, int, float, bool)):
        if isinstance(value, float) and not np.isfinite(value):
            return None
        return value
    return str(value)


def search_rows(result: lk.SearchResult, product_kind: str) -> list[dict[str, Any]]:
    """Return the useful part of a Lightkurve search result as plain rows."""
    rows: list[dict[str, Any]] = []
    available = set(result.table.colnames)
    for row in result.table:
        record: dict[str, Any] = {"product_kind": product_kind}
        for column in SEARCH_COLUMNS:
            record[column] = (
                plain_value(row[column]) if column in available else None
            )
        rows.append(record)
    return rows


def write_csv(path: Path, rows: Iterable[dict[str, Any]]) -> None:
    records = list(rows)
    path.parent.mkdir(parents=True, exist_ok=True)
    if not records:
        path.write_text("", encoding="utf-8")
        return
    fieldnames = list(records[0])
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(records)


def result_for_authors(
    result: lk.SearchResult, authors: tuple[str, ...]
) -> lk.SearchResult:
    """Select named pipelines from a Lightkurve search result."""
    author_values = np.asarray(result.table["author"], dtype=str)
    return result[np.isin(author_values, authors)]


def check_primary_coverage(result: lk.SearchResult) -> None:
    """Stop if QLP does not provide exactly one result in every sector."""
    qlp = result_for_authors(result, ("QLP",))
    found = [int(value) for value in qlp.table["sequence_number"]]
    if sorted(found) != list(SECTORS):
        raise RuntimeError(
            f"Expected one QLP light curve in Sectors {SECTORS}; found {found}"
        )


def sha256_file(path: Path, block_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(block_size):
            digest.update(block)
    return digest.hexdigest()


def download_public_file(url: str, destination: Path) -> None:
    """Download one public archive file atomically with certificate checking."""
    destination.parent.mkdir(parents=True, exist_ok=True)
    request = urllib.request.Request(
        url,
        headers={"User-Agent": "TOI-3505-research-download/1.0"},
    )
    context = ssl.create_default_context(cafile=certifi.where())
    temporary: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            dir=destination.parent,
            prefix=f".{destination.name}.",
            suffix=".part",
            delete=False,
        ) as handle:
            temporary = Path(handle.name)
            with urllib.request.urlopen(
                request, context=context, timeout=180
            ) as response:
                shutil.copyfileobj(response, handle, length=1024 * 1024)
        if temporary.stat().st_size == 0:
            raise RuntimeError(f"Archive returned an empty file for {url}")
        temporary.replace(destination)
    finally:
        if temporary is not None and temporary.exists():
            temporary.unlink()


def download_multisector_supporting_records(data_dir: Path) -> None:
    """Save target-specific Exo.MAST records and the DR122 documentation."""
    exomast_dir = data_dir / "exomast"
    release_dir = data_dir / "release_notes"
    records = (
        (f"{EXOMAST_BASE}/tces/", exomast_dir / "tce_list.json"),
        (
            f"{EXOMAST_BASE}/info/?tce=TCE_1&sector={MULTISECTOR_SCOPE}",
            exomast_dir / f"{MULTISECTOR_SCOPE}_tce1_info.json",
        ),
        (
            f"{EXOMAST_BASE}/table/?tce=TCE_1&sector={MULTISECTOR_SCOPE}",
            exomast_dir / f"{MULTISECTOR_SCOPE}_tce1_table.json",
        ),
        (
            EXOMAST_VETTING_URL,
            exomast_dir / f"{MULTISECTOR_SCOPE}_vetting_products.json",
        ),
        (
            DR122_RELEASE_NOTE_URL,
            release_dir / "tess_multisector_14_86_drn122_v01.pdf",
        ),
        (
            DR122_TARGET_INFO_URL,
            release_dir
            / "tess_multisector_14_86_drn122_targetinfo_v01.txt",
        ),
    )
    for url, destination in records:
        print(f"Downloading {destination.name}...")
        download_public_file(url, destination)


def download_mast_products_flat(products: Any, destination: Path) -> None:
    """Download a compact product table without adding another directory tree."""
    for product in products:
        filename = str(plain_value(product["productFilename"]))
        uri = str(plain_value(product["dataURI"]))
        url = "https://mast.stsci.edu/api/v0.1/Download/file/?uri=" + uri
        print(f"Downloading {filename}...")
        download_public_file(url, destination / filename)


def first_header_value(headers: list[fits.Header], *names: str) -> Any:
    for name in names:
        for header in headers:
            if name in header:
                return plain_value(header[name])
    return None


def fits_record(path: Path, data_dir: Path) -> dict[str, Any]:
    """Read compact provenance fields from one downloaded FITS file."""
    with fits.open(path, memmap=True) as hdul:
        headers = [hdu.header for hdu in hdul[:3]]
        extension_names = ",".join(
            str(hdu.header.get("EXTNAME", "PRIMARY")) for hdu in hdul
        )
    cadence = first_header_value(headers, "TIMEDEL", "FRAMETIM", "EXPOSURE")
    if isinstance(cadence, (int, float)) and 0 < float(cadence) < 2:
        cadence_seconds: float | None = float(cadence) * 86400.0
    elif isinstance(cadence, (int, float)):
        cadence_seconds = float(cadence)
    else:
        cadence_seconds = None
    return {
        "relative_path": path.relative_to(data_dir).as_posix(),
        "file_name": path.name,
        "bytes": path.stat().st_size,
        "sha256": sha256_file(path),
        "target": first_header_value(headers, "OBJECT", "TICID", "TIC_ID"),
        "sector": first_header_value(headers, "SECTOR"),
        "pipeline": first_header_value(headers, "AUTHOR", "ORIGIN", "PROCVER"),
        "data_release": first_header_value(headers, "DATA_REL", "DR_NUM"),
        "cadence_seconds": cadence_seconds,
        "camera": first_header_value(headers, "CAMERA"),
        "ccd": first_header_value(headers, "CCD"),
        "time_start_btjd": first_header_value(headers, "TSTART"),
        "time_stop_btjd": first_header_value(headers, "TSTOP"),
        "extensions": extension_names,
    }


def build_file_manifest(data_dir: Path) -> list[dict[str, Any]]:
    files = sorted(
        path
        for path in data_dir.rglob("*")
        if path.is_file()
        and path.suffix.lower()
        in {".fits", ".fit", ".fz", ".pdf", ".xml", ".json", ".txt"}
    )
    rows = []
    for path in files:
        if path.suffix.lower() in {".fits", ".fit", ".fz"}:
            rows.append(fits_record(path, data_dir))
        else:
            relative_path = path.relative_to(data_dir).as_posix()
            if "data_validation_multi_sector" in relative_path:
                sector: str | int | None = "54,81"
                pipeline = "SPOC Data Validation, combined sectors"
            elif "data_validation" in relative_path:
                sector = next(
                    (
                        value
                        for value in PIXEL_SECTORS
                        if f"s{value:04d}" in path.name
                    ),
                    None,
                )
                pipeline = "SPOC Data Validation"
            elif "exomast" in relative_path:
                sector = "54,81" if MULTISECTOR_SCOPE in path.name else None
                pipeline = "Exo.MAST"
            else:
                sector = None
                pipeline = "TESS archive documentation"
            rows.append(
                {
                    "relative_path": relative_path,
                    "file_name": path.name,
                    "bytes": path.stat().st_size,
                    "sha256": sha256_file(path),
                    "target": TIC_ID,
                    "sector": sector,
                    "pipeline": pipeline,
                    "data_release": None,
                    "cadence_seconds": None,
                    "camera": None,
                    "ccd": None,
                    "time_start_btjd": None,
                    "time_stop_btjd": None,
                    "extensions": path.suffix.lower().lstrip("."),
                }
            )
    return rows


def data_validation_products(
    target_pixel_rows: list[dict[str, Any]],
) -> tuple[Any, list[dict[str, Any]]]:
    """Return the two sector reports and the later combined report set."""
    tables = []
    rows: list[dict[str, Any]] = []
    observations = [
        (
            str(search_row["obsid"]),
            f"s{int(search_row['sequence_number']):04d}-"
            f"s{int(search_row['sequence_number']):04d}",
            str(int(search_row["sequence_number"])),
            "per-sector",
        )
        for search_row in target_pixel_rows
    ]
    observations.append(
        (
            str(MULTISECTOR_OBSID),
            MULTISECTOR_SCOPE,
            ";".join(str(value) for value in MULTISECTOR_CONTRIBUTING_SECTORS),
            "combined",
        )
    )
    for obsid, sector_scope, contributing_sectors, report_kind in observations:
        products = Observations.get_product_list(obsid)
        subgroup = np.asarray(products["productSubGroupDescription"], dtype=str)
        selected = products[np.isin(subgroup, DATA_VALIDATION_SUBGROUPS)]
        if report_kind == "combined":
            obs_ids = {str(value) for value in selected["obs_id"]}
            if obs_ids != {MULTISECTOR_OBSERVATION_ID}:
                raise RuntimeError(
                    "The expected combined SPOC observation was not returned: "
                    f"{sorted(obs_ids)}"
                )
        tables.append(selected)
        for product in selected:
            rows.append(
                {
                    "report_kind": report_kind,
                    "sector_scope": sector_scope,
                    "contributing_sectors": contributing_sectors,
                    "obsid": obsid,
                    "product_subgroup": plain_value(
                        product["productSubGroupDescription"]
                    ),
                    "product_filename": plain_value(product["productFilename"]),
                    "data_uri": plain_value(product["dataURI"]),
                    "bytes": plain_value(product["size"]),
                }
            )
    if not tables:
        raise RuntimeError("No SPOC observations were available for Data Validation")
    from astropy.table import vstack

    return vstack(tables), rows


def print_search_summary(rows: list[dict[str, Any]]) -> None:
    for row in rows:
        sector = row.get("sequence_number")
        author = row.get("author") or row.get("product_kind")
        exposure = row.get("exptime")
        name = row.get("productFilename") or "TESScut product"
        print(f"Sector {sector}: {author}, {exposure} s, {name}")


def write_download_readme(output_dir: Path) -> None:
    text = """# TOI-3505.01 TESS download record

This folder records the public TESS products used by the project.

- QLP light curves provide one common reduction for Sectors 14, 41, 54, and 81.
- SPOC and TESS-SPOC light curves provide same-observation pipeline checks in
  Sectors 54 and 81.
- TESScut files support the same custom-aperture test in all four sectors.
- SPOC target-pixel and per-sector Data Validation products are available for
  Sectors 54 and 81. The search found no corresponding SPOC target observation
  for Sectors 14 or 41.
- SPOC Data Release 122 also contains a combined multi-sector report. Its run
  is named `s0014-s0086`, but the target-information table and FITS sector
  vector show that only Sectors 54 and 81 contributed to TOI-3505.01. The
  combined time series uses 10-minute bins.

The checksum manifest includes FITS, PDF, XML, JSON, and text products. The
official SPOC reports are comparisons and diagnostic records, not independent
observations.

Official references:

- [QLP high-level science products](https://archive.stsci.edu/hlsp/qlp)
- [TESS data products](https://heasarc.gsfc.nasa.gov/docs/tess/data-products.html)
- [Exo.MAST light-curve services](https://exo.mast.stsci.edu/docs/dvdata_ws.html)
- [TESS Data Release 122](https://archive.stsci.edu/missions/tess/doc/tess_drn/tess_multisector_14_86_drn122_v01.pdf)
- [TESS target-pixel tutorial](https://heasarc.gsfc.nasa.gov/docs/tess/Target-Pixel-File-Tutorial.html)

The QLP archive page lists a CC BY 4.0 license and asks users to cite the QLP
papers and the TESS mission when publishing results.
"""
    (output_dir / "README.md").write_text(text, encoding="utf-8")


def main() -> None:
    args = parse_args()
    data_dir = args.data_dir.resolve()
    output_dir = args.output_dir.resolve()
    light_curve_dir = data_dir / "light_curves"
    tesscut_dir = data_dir / "tesscut"
    target_pixel_dir = data_dir / "target_pixels"
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"Querying MAST for {TARGET}...")
    light_curve_search = lk.search_lightcurve(
        TARGET, mission="TESS", sector=list(SECTORS)
    )
    check_primary_coverage(light_curve_search)
    light_curve_rows = search_rows(light_curve_search, "light curve")
    write_csv(output_dir / "light_curve_search.csv", light_curve_rows)
    print_search_summary(light_curve_rows)

    target_pixel_search = lk.search_targetpixelfile(
        TARGET,
        mission="TESS",
        sector=list(PIXEL_SECTORS),
        author="SPOC",
        exptime=120,
    )
    target_pixel_rows = search_rows(target_pixel_search, "target pixel file")
    write_csv(output_dir / "target_pixel_search.csv", target_pixel_rows)
    validation_products, validation_rows = data_validation_products(
        target_pixel_rows
    )
    write_csv(output_dir / "data_validation_search.csv", validation_rows)

    tesscut_searches: list[lk.SearchResult] = []
    tesscut_rows: list[dict[str, Any]] = []
    for sector in SECTORS:
        result = lk.search_tesscut(TARGET, sector=sector)
        if len(result) != 1:
            raise RuntimeError(
                f"Expected one TESScut result for Sector {sector}; found {len(result)}"
            )
        tesscut_searches.append(result)
        tesscut_rows.extend(search_rows(result, "TESScut"))
    write_csv(output_dir / "tesscut_search.csv", tesscut_rows)

    if not args.inventory_only:
        selected_light_curves = result_for_authors(
            light_curve_search, LIGHT_CURVE_AUTHORS
        )
        print(f"Downloading {len(selected_light_curves)} light curves...")
        selected_light_curves.download_all(download_dir=str(light_curve_dir))

        validation_dir = data_dir / "data_validation"
        multisector_validation_dir = data_dir / "data_validation_multi_sector"
        combined_mask = np.asarray(
            validation_products["obsID"], dtype=int
        ) == MULTISECTOR_OBSID
        per_sector_products = validation_products[~combined_mask]
        multisector_products = validation_products[combined_mask]
        print(
            f"Downloading {len(per_sector_products)} per-sector SPOC "
            "Data Validation products..."
        )
        Observations.download_products(
            per_sector_products,
            download_dir=str(validation_dir),
            mrp_only=False,
        )
        print(
            f"Downloading {len(multisector_products)} combined-sector SPOC "
            "Data Validation products..."
        )
        download_mast_products_flat(
            multisector_products, multisector_validation_dir
        )
        download_multisector_supporting_records(data_dir)

        if not args.skip_pixels:
            print(f"Downloading {len(target_pixel_search)} SPOC target-pixel files...")
            target_pixel_search.download_all(download_dir=str(target_pixel_dir))
            for sector, result in zip(SECTORS, tesscut_searches):
                print(f"Downloading {CUTOUT_SIZE[0]}x{CUTOUT_SIZE[1]} TESScut Sector {sector}...")
                result.download_all(
                    download_dir=str(tesscut_dir), cutout_size=CUTOUT_SIZE
                )

    manifest = build_file_manifest(data_dir) if data_dir.exists() else []
    write_csv(output_dir / "file_manifest.csv", manifest)
    retrieved_at = datetime.now(timezone.utc).isoformat()
    summary = {
        "target": "TOI-3505.01",
        "tic_id": TIC_ID,
        "sectors": list(SECTORS),
        "primary_light_curve_pipeline": "QLP",
        "cross_check_pipelines": ["SPOC", "TESS-SPOC"],
        "official_spoc_report_sets": [
            "Sector 54",
            "Sector 81",
            "combined Sectors 54 and 81 from the s0014-s0086 search",
        ],
        "pixel_cutout_size": list(CUTOUT_SIZE),
        "retrieved_utc": retrieved_at,
        "mast_query": f"{TARGET}; mission=TESS; sectors={list(SECTORS)}",
        "downloaded_file_count": len(manifest),
        "downloaded_bytes": int(sum(row["bytes"] for row in manifest)),
        "software": {
            "python": platform.python_version(),
            "operating_system": platform.platform(),
            "lightkurve": lk.__version__,
            "astropy": astropy.__version__,
            "script": Path(__file__).name,
        },
        "notes": [
            "QLP is the uniform four-sector light-curve set.",
            "SPOC and TESS-SPOC are checks of the same TESS observations, not independent observations.",
            "TESScut files support a common custom-aperture check in all four sectors.",
            "Official SPOC Data Validation products exist for Sectors 54 and 81; Sectors 14 and 41 have no SPOC target observation in this search.",
            "DR122's s0014-s0086 label is the pipeline search range; only Sectors 54 and 81 contributed for this target, at 10-minute cadence.",
        ],
    }
    (output_dir / "download_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    write_download_readme(output_dir)
    print(f"Saved search tables and manifest to {output_dir}")
    print(
        f"Manifested {len(manifest)} downloaded files "
        f"({summary['downloaded_bytes']:,} bytes)."
    )


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise
