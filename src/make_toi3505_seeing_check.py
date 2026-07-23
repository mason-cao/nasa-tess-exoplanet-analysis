"""Make the target-position and aperture checks for TOI-3505.01.

The aperture values are measurements from AstroImageJ's Seeing Profile on the
first plate-solved science image. This script checks the measured center
against the plate solution and makes two figures that are easy to review.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import tempfile
from datetime import UTC, datetime
from pathlib import Path

os.environ.setdefault(
    "MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "nasa-tess-matplotlib-cache")
)

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Patch
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.visualization import AsinhStretch, ImageNormalize, ZScaleInterval
from astropy.wcs import WCS
from astropy.wcs.utils import skycoord_to_pixel


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = (
    ROOT
    / "data"
    / "ground"
    / "toi3505"
    / "plate_solved"
    / "TOI_3505.01_50.000s_R-0001_wcs.fits"
)
DEFAULT_SOLUTION = ROOT / "outputs" / "toi3505_plate_solve" / "solution.json"
DEFAULT_SOURCES = ROOT / "outputs" / "toi3505_plate_solve" / "source_list.csv"
DEFAULT_OUTPUT_DIR = ROOT / "outputs" / "toi3505_seeing"

# Values shown by AstroImageJ's Seeing Profile window for the target.
AIJ_CENTER_X_FITS = 1850.58
AIJ_CENTER_Y_FITS = 1742.82
AIJ_FWHM_PIXELS = 13.64
SOURCE_RADIUS_PIXELS = 35.0
BACKGROUND_INNER_PIXELS = 70.0
BACKGROUND_OUTER_PIXELS = 139.0
PROFILE_CUTOFF = 0.010
POSITION_LIMIT_PIXELS = 2.0
SEEING_CHECK_VERSION = "1.0.0"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--solution", type=Path, default=DEFAULT_SOLUTION)
    parser.add_argument("--sources", type=Path, default=DEFAULT_SOURCES)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--center-x-fits", type=float, default=AIJ_CENTER_X_FITS)
    parser.add_argument("--center-y-fits", type=float, default=AIJ_CENTER_Y_FITS)
    parser.add_argument("--fwhm", type=float, default=AIJ_FWHM_PIXELS)
    parser.add_argument("--source-radius", type=float, default=SOURCE_RADIUS_PIXELS)
    parser.add_argument(
        "--background-inner", type=float, default=BACKGROUND_INNER_PIXELS
    )
    parser.add_argument(
        "--background-outer", type=float, default=BACKGROUND_OUTER_PIXELS
    )
    return parser.parse_args()


def crop_with_fits_extent(
    image: np.ndarray,
    center_x_fits: float,
    center_y_fits: float,
    half_width: int,
) -> tuple[np.ndarray, tuple[float, float, float, float]]:
    center_x = int(round(center_x_fits - 1.0))
    center_y = int(round(center_y_fits - 1.0))
    x_start = max(0, center_x - half_width)
    x_end = min(image.shape[1], center_x + half_width + 1)
    y_start = max(0, center_y - half_width)
    y_end = min(image.shape[0], center_y + half_width + 1)
    crop = np.asarray(image[y_start:y_end, x_start:x_end], dtype=float)
    extent = (
        x_start + 0.5,
        x_end + 0.5,
        y_start + 0.5,
        y_end + 0.5,
    )
    return crop, extent


def image_stretch(crop: np.ndarray) -> ImageNormalize:
    return ImageNormalize(crop, interval=ZScaleInterval(), stretch=AsinhStretch())


def nearest_other_source(
    source_path: Path,
    *,
    center_x_fits: float,
    center_y_fits: float,
) -> dict[str, float]:
    nearby: list[dict[str, float]] = []
    with source_path.open(newline="", encoding="utf-8") as file_obj:
        for row in csv.DictReader(file_obj):
            x_fits = float(row["x_pixel"]) + 1.0
            y_fits = float(row["y_pixel"]) + 1.0
            distance = math.hypot(
                x_fits - center_x_fits,
                y_fits - center_y_fits,
            )
            if distance > 5.0:
                nearby.append(
                    {
                        "x_fits": x_fits,
                        "y_fits": y_fits,
                        "distance_pixels": distance,
                        "width_pixels": float(row["width_pixels"]),
                    }
                )
    if not nearby:
        raise RuntimeError("No other detected stars were available for the aperture check")
    return min(nearby, key=lambda row: row["distance_pixels"])


def make_aperture_figure(
    image: np.ndarray,
    output_path: Path,
    *,
    center_x_fits: float,
    center_y_fits: float,
    source_radius: float,
    background_inner: float,
    background_outer: float,
    nearby_star_x_fits: float,
    nearby_star_y_fits: float,
) -> None:
    crop, extent = crop_with_fits_extent(
        image,
        center_x_fits,
        center_y_fits,
        half_width=175,
    )
    fig, ax = plt.subplots(figsize=(7.2, 7.2))
    ax.imshow(
        crop,
        origin="lower",
        extent=extent,
        cmap="gray_r",
        norm=image_stretch(crop),
        interpolation="nearest",
    )

    source_color = "#d24a3a"
    background_color = "#18807a"
    for radius, color, width in (
        (source_radius, source_color, 2.0),
        (background_inner, background_color, 1.8),
        (background_outer, background_color, 1.8),
    ):
        ax.add_patch(
            Circle(
                (center_x_fits, center_y_fits),
                radius,
                fill=False,
                edgecolor=color,
                linewidth=width,
            )
        )
    ax.scatter(
        [center_x_fits],
        [center_y_fits],
        marker="+",
        s=90,
        linewidths=1.8,
        color="#e2ad21",
        label="TOI-3505.01",
        zorder=5,
    )
    ax.scatter(
        [nearby_star_x_fits],
        [nearby_star_y_fits],
        marker="x",
        s=80,
        linewidths=1.7,
        color="#7b4b91",
        label="Nearby star",
        zorder=5,
    )
    handles, labels = ax.get_legend_handles_labels()
    handles.extend(
        [
            Patch(facecolor="none", edgecolor=source_color, label="Source aperture"),
            Patch(
                facecolor="none",
                edgecolor=background_color,
                label="Background area",
            ),
        ]
    )
    labels.extend(["Source aperture", "Background area"])
    ax.legend(handles, labels, loc="upper right", frameon=True)
    ax.set_title("TOI-3505.01 aperture on the first image")
    ax.set_xlabel("FITS x pixel")
    ax.set_ylabel("FITS y pixel")
    ax.text(
        0.02,
        0.02,
        (
            f"Source radius: {source_radius:.0f} pixels\n"
            f"Background: {background_inner:.0f}-{background_outer:.0f} pixels"
        ),
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=10,
        bbox={"facecolor": "white", "edgecolor": "#bdbdbd", "alpha": 0.9},
    )
    fig.tight_layout()
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def make_position_figure(
    image: np.ndarray,
    output_path: Path,
    *,
    predicted_x_fits: float,
    predicted_y_fits: float,
    measured_x_fits: float,
    measured_y_fits: float,
) -> None:
    crop, extent = crop_with_fits_extent(
        image,
        measured_x_fits,
        measured_y_fits,
        half_width=28,
    )
    fig, ax = plt.subplots(figsize=(6.5, 6.1))
    ax.imshow(
        crop,
        origin="lower",
        extent=extent,
        cmap="gray_r",
        norm=image_stretch(crop),
        interpolation="nearest",
    )
    ax.scatter(
        [predicted_x_fits],
        [predicted_y_fits],
        marker="+",
        s=170,
        linewidths=2.2,
        color="#cc3d3d",
        label="Position from the plate solution",
        zorder=6,
    )
    ax.scatter(
        [measured_x_fits],
        [measured_y_fits],
        marker="x",
        s=90,
        linewidths=1.8,
        color="#2474a6",
        label="Center measured by AstroImageJ",
        zorder=7,
    )
    ax.legend(loc="upper right", frameon=True, fontsize=9)
    ax.set_title("TOI-3505.01 position check")
    ax.set_xlabel("FITS x pixel")
    ax.set_ylabel("FITS y pixel")
    fig.tight_layout()
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def write_notes(
    output_dir: Path,
    *,
    image_name: str,
    pixel_scale: float,
    predicted_x: float,
    predicted_y: float,
    center_x: float,
    center_y: float,
    position_difference: float,
    fwhm: float,
    source_radius: float,
    background_inner: float,
    background_outer: float,
    nearby_star_distance: float,
    clean_background_outer: float,
) -> None:
    notes = f"""# TOI-3505.01 target and aperture check

I plate-solved the first reduced image and used the target coordinates in its
header to locate TOI-3505.01. The file I used is `{image_name}`. The plate
solution puts the target at FITS pixel
({predicted_x:.2f}, {predicted_y:.2f}). AstroImageJ centered the star at
({center_x:.2f}, {center_y:.2f}), a difference of {position_difference:.2f}
pixel. This confirms that the Seeing Profile was measured on the target.

The AstroImageJ Seeing Profile measured a FWHM of {fwhm:.2f} pixels. Its
starting photometry settings are:

- Source radius: {source_radius:.0f} pixels
- Background inner radius: {background_inner:.0f} pixels
- Background outer radius: {background_outer:.0f} pixels

The plate scale is {pixel_scale:.3f} arcseconds per pixel. A nearby star is
{nearby_star_distance:.2f} pixels from the target, so it falls inside the
{background_outer:.0f}-pixel outer background ring. I will also test a
{background_inner:.0f}-{clean_background_outer:.0f} pixel background area,
which stays one source radius away from that star. These are starting values
for Multi-Aperture. I still need to choose the final aperture by checking the
seeing and comparison-star trends across the full set of usable images.
"""
    (output_dir / "README.md").write_text(notes, encoding="utf-8")

    discord = f"""I plate-solved the first reduced TOI-3505.01 image at {pixel_scale:.3f} arcsec/pixel. The target position from the plate solution and the center measured by AstroImageJ are only {position_difference:.2f} pixel apart. The Seeing Profile gave a FWHM of {fwhm:.2f} pixels and starting radii of {source_radius:.0f}, {background_inner:.0f}, and {background_outer:.0f} pixels. A nearby star is {nearby_star_distance:.1f} pixels away, so it falls inside the outer background ring. I am going to compare that with a {background_inner:.0f}-{clean_background_outer:.0f} pixel background area while I check the seeing and comparison-star trends across all 281 usable images.
"""
    (output_dir / "discord_update.txt").write_text(discord, encoding="utf-8")


def main() -> None:
    args = parse_args()
    input_path = args.input.resolve()
    solution_path = args.solution.resolve()
    source_path = args.sources.resolve()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    solution = json.loads(solution_path.read_text(encoding="utf-8"))
    if solution.get("status") != "pass":
        raise RuntimeError("The plate-solution report did not pass")
    nearby_star = nearest_other_source(
        source_path,
        center_x_fits=args.center_x_fits,
        center_y_fits=args.center_y_fits,
    )
    background_contains_star = (
        args.background_inner
        <= nearby_star["distance_pixels"]
        <= args.background_outer
    )
    clear_limit = nearby_star["distance_pixels"] - args.source_radius
    clean_background_outer = 5.0 * math.floor(clear_limit / 5.0)
    if clean_background_outer <= args.background_inner:
        raise RuntimeError("There is not enough clear space for a background area")

    with fits.open(input_path, memmap=True, checksum=True) as hdul:
        image = hdul[0].data
        header = hdul[0].header.copy()
        checksum_ok = bool(hdul[0].verify_checksum())
        datasum_ok = bool(hdul[0].verify_datasum())
        if not bool(header.get("PLTSOLVD", False)):
            raise RuntimeError("The input image is not marked as plate-solved")

        target = SkyCoord(
            ra=float(solution["target_ra_degrees"]),
            dec=float(solution["target_dec_degrees"]),
            unit="deg",
        )
        predicted_x_value, predicted_y_value = skycoord_to_pixel(
            target,
            WCS(header),
            origin=1,
        )
        predicted_x = float(np.asarray(predicted_x_value).squeeze())
        predicted_y = float(np.asarray(predicted_y_value).squeeze())
        position_difference = math.hypot(
            predicted_x - args.center_x_fits,
            predicted_y - args.center_y_fits,
        )
        if position_difference > POSITION_LIMIT_PIXELS:
            raise RuntimeError(
                "The AstroImageJ center does not match the plate-solved target position"
            )

        make_aperture_figure(
            image,
            output_dir / "01_target_and_aperture.png",
            center_x_fits=args.center_x_fits,
            center_y_fits=args.center_y_fits,
            source_radius=args.source_radius,
            background_inner=args.background_inner,
            background_outer=args.background_outer,
            nearby_star_x_fits=nearby_star["x_fits"],
            nearby_star_y_fits=nearby_star["y_fits"],
        )
        make_position_figure(
            image,
            output_dir / "03_position_check.png",
            predicted_x_fits=predicted_x,
            predicted_y_fits=predicted_y,
            measured_x_fits=args.center_x_fits,
            measured_y_fits=args.center_y_fits,
        )

    pixel_scale = float(solution["pixel_scale_arcsec"])
    report = {
        "status": "pass",
        "seeing_check_version": SEEING_CHECK_VERSION,
        "created_utc": datetime.now(UTC).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "image_file": input_path.name,
        "target": "TOI-3505.01",
        "filter": "Johnson-Cousins R",
        "fits_checksum_valid": checksum_ok,
        "fits_datasum_valid": datasum_ok,
        "pixel_scale_arcsec": pixel_scale,
        "plate_solution_x_fits": predicted_x,
        "plate_solution_y_fits": predicted_y,
        "astroimagej_center_x_fits": args.center_x_fits,
        "astroimagej_center_y_fits": args.center_y_fits,
        "position_difference_pixels": position_difference,
        "position_limit_pixels": POSITION_LIMIT_PIXELS,
        "fwhm_pixels": args.fwhm,
        "fwhm_arcsec": args.fwhm * pixel_scale,
        "profile_cutoff_fraction": PROFILE_CUTOFF,
        "source_radius_pixels": args.source_radius,
        "background_inner_radius_pixels": args.background_inner,
        "background_outer_radius_pixels": args.background_outer,
        "nearest_other_star_distance_pixels": nearby_star["distance_pixels"],
        "nearest_other_star_width_pixels": nearby_star["width_pixels"],
        "starting_background_contains_detected_star": background_contains_star,
        "clear_background_inner_radius_pixels": args.background_inner,
        "clear_background_outer_radius_pixels": clean_background_outer,
        "clear_background_rule": "Keep one source-radius between the background edge and the nearby star center",
        "source": "AstroImageJ Seeing Profile on the first plate-solved image",
    }
    (output_dir / "aperture_settings.json").write_text(
        json.dumps(report, indent=2) + "\n",
        encoding="utf-8",
    )
    write_notes(
        output_dir,
        image_name=input_path.name,
        pixel_scale=pixel_scale,
        predicted_x=predicted_x,
        predicted_y=predicted_y,
        center_x=args.center_x_fits,
        center_y=args.center_y_fits,
        position_difference=position_difference,
        fwhm=args.fwhm,
        source_radius=args.source_radius,
        background_inner=args.background_inner,
        background_outer=args.background_outer,
        nearby_star_distance=nearby_star["distance_pixels"],
        clean_background_outer=clean_background_outer,
    )
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
