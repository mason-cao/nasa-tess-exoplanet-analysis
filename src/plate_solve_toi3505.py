"""Plate-solve the first usable TOI-3505.01 image with Astrometry.net.

The source positions are measured locally, so the full science image does not
need to be uploaded. The Astrometry.net key is read without displaying it and
is never written to a file. The search settings follow the Schar tutorial:
0.350 +/- 0.050 arcseconds per pixel, a 20-arcminute position radius, and a
second-order sky model.
"""

from __future__ import annotations

import argparse
import getpass
import json
import math
import os
import time
from datetime import UTC, datetime
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS
from astropy.wcs.utils import skycoord_to_pixel
from astroquery.astrometry_net import AstrometryNet
from scipy.ndimage import maximum_filter

from align_toi3505 import measure_star, robust_center_and_noise


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = (
    ROOT
    / "data"
    / "ground"
    / "toi3505"
    / "reduced"
    / "TOI_3505.01_50.000s_R-0001_out.fits"
)
DEFAULT_OUTPUT = (
    ROOT
    / "data"
    / "ground"
    / "toi3505"
    / "plate_solved"
    / "TOI_3505.01_50.000s_R-0001_wcs.fits"
)
DEFAULT_REPORT_DIR = ROOT / "outputs" / "toi3505_plate_solve"
PLATE_SOLVE_VERSION = "1.0.0"
EXPECTED_SHAPE = (4096, 4096)


def find_sources(image: np.ndarray, *, maximum_sources: int = 120) -> pd.DataFrame:
    center, noise = robust_center_and_noise(image[::8, ::8])
    local_maximum = maximum_filter(image, size=15, mode="nearest")
    yy, xx = np.nonzero(
        (image == local_maximum)
        & (image > center + 10.0 * noise)
        & (image < 45_000.0)
    )
    rows: list[dict[str, float]] = []
    for x_peak, y_peak in zip(xx, yy):
        star = measure_star(image, int(x_peak), int(y_peak))
        if star is None:
            continue
        if not (6.0 <= star.fwhm <= 14.0):
            continue
        if star.flux < 20_000.0 or star.peak_fraction >= 0.08:
            continue
        if math.hypot(star.x - x_peak, star.y - y_peak) >= 5.0:
            continue
        rows.append(
            {
                "x_pixel": star.x,
                "y_pixel": star.y,
                "brightness_adu": star.flux,
                "width_pixels": star.fwhm,
            }
        )
    table = pd.DataFrame(rows).sort_values(
        "brightness_adu", ascending=False, ignore_index=True
    )
    table = table.head(maximum_sources).copy()
    table.insert(0, "source", np.arange(1, len(table) + 1))
    if len(table) < 30:
        raise RuntimeError(f"Only {len(table)} suitable stars were found")
    return table


def target_coordinates(header: fits.Header) -> SkyCoord:
    ra_hours = float(header["RAOBJ2K"])
    dec_degrees = float(header["DECOBJ2K"])
    return SkyCoord(ra=ra_hours, dec=dec_degrees, unit=("hourangle", "degree"))


def merge_wcs_header(
    source_header: fits.Header,
    wcs_header: fits.Header,
    *,
    submission_id: int | str,
) -> fits.Header:
    merged = source_header.copy()
    for keyword in ("CHECKSUM", "DATASUM"):
        if keyword in merged:
            del merged[keyword]
    for card in wcs_header.cards:
        if card.keyword in {"SIMPLE", "BITPIX", "NAXIS", "NAXIS1", "NAXIS2"}:
            continue
        if card.keyword == "COMMENT":
            continue
        if card.keyword == "HISTORY":
            merged.add_history(str(card.value))
            continue
        merged[card.keyword] = (card.value, card.comment)
    merged["PLTSOLVD"] = (True, "Plate solution added with Astrometry.net")
    merged["PLTSUBID"] = (str(submission_id), "Astrometry.net submission number")
    merged["PLTVER"] = (PLATE_SOLVE_VERSION, "Plate-solving script version")
    merged["PLTDATE"] = (
        datetime.now(UTC).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "UTC plate-solving time",
    )
    merged.add_history("WCS solution found from a local list of measured star positions.")
    return merged


def write_fits_atomic(
    path: Path,
    image: np.ndarray,
    header: fits.Header,
    *,
    overwrite: bool,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists() and not overwrite:
        raise FileExistsError(path)
    temporary = path.with_name(f".{path.name}.partial")
    if temporary.exists():
        temporary.unlink()
    fits.PrimaryHDU(data=image, header=header).writeto(
        temporary,
        overwrite=True,
        checksum=True,
        output_verify="fix+warn",
    )
    os.replace(temporary, path)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--report-dir", type=Path, default=DEFAULT_REPORT_DIR)
    parser.add_argument("--maximum-sources", type=int, default=120)
    parser.add_argument("--solve-timeout", type=int, default=300)
    parser.add_argument(
        "--detect-only",
        action="store_true",
        help="Measure the source list without contacting Astrometry.net",
    )
    parser.add_argument(
        "--finish-existing",
        action="store_true",
        help="Finish the report from an existing plate-solved FITS file",
    )
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    started = time.perf_counter()
    input_path = args.input.resolve()
    output_path = args.output.resolve()
    report_dir = args.report_dir.resolve()
    report_dir.mkdir(parents=True, exist_ok=True)

    with fits.open(input_path, memmap=True, checksum=False) as hdul:
        image = np.array(hdul[0].data, dtype=np.float32, copy=True)
        header = hdul[0].header.copy()
    if image.shape != EXPECTED_SHAPE:
        raise RuntimeError(f"Unexpected image shape: {image.shape}")
    if not bool(header.get("REDUCED", False)):
        raise RuntimeError("The input image is not marked as calibrated")

    sources = find_sources(image, maximum_sources=args.maximum_sources)
    source_path = report_dir / "source_list.csv"
    sources.to_csv(source_path, index=False)
    print(f"Found {len(sources)} suitable stars", flush=True)
    if args.detect_only:
        return

    target = target_coordinates(header)
    if args.finish_existing:
        if not output_path.exists():
            raise FileNotFoundError(output_path)
        solved_header = fits.getheader(output_path)
        if not bool(solved_header.get("PLTSOLVD", False)):
            raise RuntimeError("The existing output is not marked as plate-solved")
        submission_id = str(solved_header["PLTSUBID"])
    else:
        api_key = os.environ.get("ASTROMETRY_NET_API_KEY")
        if not api_key:
            api_key = getpass.getpass("Astrometry.net API key: ")
        if not api_key.strip():
            raise RuntimeError("An Astrometry.net API key is required")

        solver = AstrometryNet()
        solver.api_key = api_key.strip()
        print("Submitting only the measured star positions", flush=True)
        result = solver.solve_from_source_list(
            sources["x_pixel"].to_numpy() + 1.0,
            sources["y_pixel"].to_numpy() + 1.0,
            image_width=image.shape[1],
            image_height=image.shape[0],
            solve_timeout=args.solve_timeout,
            verbose=True,
            return_submission_id=True,
            scale_units="arcsecperpix",
            scale_type="ul",
            scale_lower=0.300,
            scale_upper=0.400,
            center_ra=float(target.ra.degree),
            center_dec=float(target.dec.degree),
            radius=20.0 / 60.0,
            tweak_order=2,
            publicly_visible="n",
        )
        api_key = ""
        solver.api_key = ""
        wcs_header, submission_id = result
        if not isinstance(wcs_header, fits.Header) or not wcs_header:
            raise RuntimeError("Astrometry.net did not find a plate solution")

        solved_header = merge_wcs_header(
            header,
            wcs_header,
            submission_id=submission_id,
        )
    wcs = WCS(solved_header)
    target_x_value, target_y_value = skycoord_to_pixel(target, wcs, origin=0)
    target_x = float(np.asarray(target_x_value).squeeze())
    target_y = float(np.asarray(target_y_value).squeeze())
    if not (0 <= target_x < image.shape[1] and 0 <= target_y < image.shape[0]):
        raise RuntimeError(
            f"The solved target position is outside the image: ({target_x}, {target_y})"
        )
    measured_target = measure_star(image, target_x, target_y)
    if measured_target is None:
        raise RuntimeError("No star was found at the solved target position")
    target_distance = math.hypot(
        measured_target.x - target_x,
        measured_target.y - target_y,
    )
    if target_distance > 8.0:
        raise RuntimeError(
            f"The nearest measured star is {target_distance:.2f} pixels from the target"
        )

    solved_header["TARG_X"] = (measured_target.x + 1.0, "Target x position, FITS pixels")
    solved_header["TARG_Y"] = (measured_target.y + 1.0, "Target y position, FITS pixels")
    if not args.finish_existing:
        write_fits_atomic(
            output_path,
            image,
            solved_header,
            overwrite=args.overwrite,
        )
        wcs_header.totextfile(report_dir / "wcs_header.txt", overwrite=True)

    scales_degrees = np.asarray(
        [scale.to_value("deg") for scale in wcs.proj_plane_pixel_scales()]
    )
    pixel_scale = float(np.mean(np.abs(scales_degrees)) * 3600.0)
    report = {
        "status": "pass",
        "plate_solve_version": PLATE_SOLVE_VERSION,
        "created_utc": datetime.now(UTC).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "input_file": input_path.name,
        "output_file": output_path.name,
        "submission_id": str(submission_id),
        "sources_submitted": len(sources),
        "full_image_uploaded": False,
        "image_width_pixels": image.shape[1],
        "image_height_pixels": image.shape[0],
        "pixel_scale_arcsec": pixel_scale,
        "search_scale_arcsec_per_pixel": [0.300, 0.400],
        "search_radius_degrees": 20.0 / 60.0,
        "target_ra_degrees": float(target.ra.degree),
        "target_dec_degrees": float(target.dec.degree),
        "plate_solution_target_x_pixel_zero_based": target_x,
        "plate_solution_target_y_pixel_zero_based": target_y,
        "plate_solution_target_x_pixel_fits": target_x + 1.0,
        "plate_solution_target_y_pixel_fits": target_y + 1.0,
        "target_x_pixel_zero_based": measured_target.x,
        "target_y_pixel_zero_based": measured_target.y,
        "target_x_pixel_fits": measured_target.x + 1.0,
        "target_y_pixel_fits": measured_target.y + 1.0,
        "wcs_to_star_distance_pixels": target_distance,
        "target_width_pixels": measured_target.fwhm,
        "api_key_stored": False,
        "elapsed_seconds": float(time.perf_counter() - started),
    }
    (report_dir / "solution.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
