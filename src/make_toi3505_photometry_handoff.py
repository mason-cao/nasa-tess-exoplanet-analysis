"""Prepare a starting comparison-star list for TOI-3505.01 photometry."""

from __future__ import annotations

import argparse
import json
import os
import tempfile
from pathlib import Path

os.environ.setdefault(
    "MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "nasa-tess-matplotlib-cache")
)

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.visualization import AsinhStretch, ImageNormalize, ZScaleInterval


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_IMAGE = (
    ROOT
    / "data"
    / "ground"
    / "toi3505"
    / "plate_solved"
    / "TOI_3505.01_50.000s_R-0001_wcs.fits"
)
DEFAULT_REFERENCES = ROOT / "outputs" / "toi3505_alignment" / "reference_stars.csv"
DEFAULT_SETTINGS = ROOT / "outputs" / "toi3505_seeing" / "aperture_settings.json"
DEFAULT_OUTPUT_DIR = ROOT / "outputs" / "toi3505_seeing"
COMPARISON_COUNT = 10
EDGE_MARGIN_PIXELS = 500.0
MINIMUM_TARGET_DISTANCE_PIXELS = 300.0
MINIMUM_FLUX_RATIO = 0.50
MAXIMUM_FLUX_RATIO = 1.50
MAXIMUM_WIDTH_DIFFERENCE_PIXELS = 1.0
COURSE_CIRCLE_RADIUS_PIXELS = 416.667


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--image", type=Path, default=DEFAULT_IMAGE)
    parser.add_argument("--references", type=Path, default=DEFAULT_REFERENCES)
    parser.add_argument("--settings", type=Path, default=DEFAULT_SETTINGS)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--comparison-count", type=int, default=COMPARISON_COUNT)
    return parser.parse_args()


def choose_comparison_stars(
    references: pd.DataFrame,
    *,
    target_x_zero_based: float,
    target_y_zero_based: float,
    image_width: int,
    image_height: int,
    comparison_count: int = COMPARISON_COUNT,
) -> tuple[pd.Series, pd.DataFrame]:
    table = references.copy()
    table["target_distance_pixels"] = np.hypot(
        table["x_pixel"] - target_x_zero_based,
        table["y_pixel"] - target_y_zero_based,
    )
    target_index = table["target_distance_pixels"].idxmin()
    target = table.loc[target_index].copy()
    if float(target["target_distance_pixels"]) > 5.0:
        raise RuntimeError("The target was not found in the aligned reference-star list")

    target_flux = float(target["brightness_adu"])
    target_width = float(target["width_pixels"])
    table["brightness_ratio"] = table["brightness_adu"] / target_flux
    table["width_difference_pixels"] = np.abs(table["width_pixels"] - target_width)
    table["edge_distance_pixels"] = np.minimum.reduce(
        [
            table["x_pixel"],
            image_width - 1.0 - table["x_pixel"],
            table["y_pixel"],
            image_height - 1.0 - table["y_pixel"],
        ]
    )
    table["brightness_score"] = np.abs(np.log(table["brightness_ratio"]))
    choices = table.loc[
        (table.index != target_index)
        & (table["target_distance_pixels"] >= MINIMUM_TARGET_DISTANCE_PIXELS)
        & (table["edge_distance_pixels"] >= EDGE_MARGIN_PIXELS)
        & (table["brightness_ratio"] >= MINIMUM_FLUX_RATIO)
        & (table["brightness_ratio"] <= MAXIMUM_FLUX_RATIO)
        & (table["width_difference_pixels"] <= MAXIMUM_WIDTH_DIFFERENCE_PIXELS)
    ].copy()
    choices = choices.sort_values(
        ["brightness_score", "width_difference_pixels", "edge_distance_pixels"],
        ascending=[True, True, False],
    ).head(comparison_count)
    if len(choices) < comparison_count:
        raise RuntimeError(
            f"Only {len(choices)} suitable comparison stars were found; "
            f"{comparison_count} were requested"
        )
    choices = choices.reset_index(drop=True)
    choices.insert(0, "aperture", [f"C{number}" for number in range(2, 2 + len(choices))])
    return target, choices


def make_field_figure(
    image: np.ndarray,
    output_path: Path,
    *,
    target_x_fits: float,
    target_y_fits: float,
    comparisons: pd.DataFrame,
) -> None:
    sampled = np.asarray(image[::4, ::4], dtype=float)
    norm = ImageNormalize(
        sampled,
        interval=ZScaleInterval(),
        stretch=AsinhStretch(),
    )
    fig, ax = plt.subplots(figsize=(8.2, 8.2))
    ax.imshow(
        sampled,
        origin="lower",
        extent=(0.5, image.shape[1] + 0.5, 0.5, image.shape[0] + 0.5),
        cmap="gray_r",
        norm=norm,
        interpolation="nearest",
    )
    ax.add_patch(
        Circle(
            (target_x_fits, target_y_fits),
            COURSE_CIRCLE_RADIUS_PIXELS,
            fill=False,
            edgecolor="#d69a22",
            linestyle="--",
            linewidth=1.5,
            label="2.5 arcminute circle",
        )
    )
    ax.scatter(
        [target_x_fits],
        [target_y_fits],
        marker="o",
        s=75,
        facecolors="none",
        edgecolors="#24844f",
        linewidths=2.0,
        label="Target",
        zorder=5,
    )
    ax.text(
        target_x_fits + 35,
        target_y_fits + 35,
        "T1",
        color="#17643a",
        weight="bold",
        fontsize=10,
        zorder=6,
    )
    x_values = comparisons["x_fits"].to_numpy()
    y_values = comparisons["y_fits"].to_numpy()
    ax.scatter(
        x_values,
        y_values,
        marker="o",
        s=60,
        facecolors="none",
        edgecolors="#c4423b",
        linewidths=1.7,
        label="Comparison stars",
        zorder=5,
    )
    for row in comparisons.itertuples(index=False):
        ax.text(
            row.x_fits + 30,
            row.y_fits + 30,
            row.aperture,
            color="#a92f2a",
            weight="bold",
            fontsize=9,
            zorder=6,
        )
    ax.legend(loc="upper right", frameon=True)
    ax.set_title("TOI-3505.01 comparison-star starting list")
    ax.set_xlabel("FITS x pixel")
    ax.set_ylabel("FITS y pixel")
    ax.set_xlim(0.5, image.shape[1] + 0.5)
    ax.set_ylim(0.5, image.shape[0] + 0.5)
    fig.tight_layout()
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def write_handoff(
    output_dir: Path,
    *,
    settings: dict,
    comparison_count: int,
) -> None:
    source_radius = float(settings["source_radius_pixels"])
    background_inner = float(settings["background_inner_radius_pixels"])
    background_outer = float(settings["background_outer_radius_pixels"])
    clear_outer = float(settings["clear_background_outer_radius_pixels"])
    handoff = f"""# Next step: Multi-Aperture photometry

The next step in the Schar light-curve tutorial is to make the measurement
table in AstroImageJ.

## Images

Use the 281 aligned images ending in `_aligned.fits`. Frames 0282 and 0283
are not usable and should stay out of the stack. Only the first image has a
plate solution, so use the aligned-image settings:

- Uncheck **Use RA/Dec to locate aperture positions**.
- Check **Halt processing on WCS or centroid error**.

## Starting aperture settings

- Source radius: {source_radius:.0f} pixels
- Background inner radius: {background_inner:.0f} pixels
- AstroImageJ starting outer radius: {background_outer:.0f} pixels
- Clear outer-radius test: {clear_outer:.0f} pixels
- CCD readout noise: 1.414214
- CCD dark current per second: 0.012283

The {background_outer:.0f}-pixel outer ring crosses a nearby star, so compare
it with the {background_inner:.0f}-{clear_outer:.0f} pixel background area
before choosing the final setting.

## Apertures

Click TOI-3505.01 first so it is T1. Then use the {comparison_count} stars in
`comparison_star_candidates.csv` as a starting list. They have similar sizes
and brightnesses to the target and are away from the image edges. They are not
automatically the final comparison set. After the table is made, plot each
comparison star's relative flux and remove stars with large scatter or clear
trends.

Save both the aperture file and the measurement table. The measurement table
will be used with `SCHAR_Plot_config.plotcfg` for the light curve, seeing,
airmass, total comparison counts, and target-position plots.
"""
    (output_dir / "photometry_handoff.md").write_text(handoff, encoding="utf-8")


def main() -> None:
    args = parse_args()
    image_path = args.image.resolve()
    reference_path = args.references.resolve()
    settings_path = args.settings.resolve()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    settings = json.loads(settings_path.read_text(encoding="utf-8"))
    if settings.get("status") != "pass":
        raise RuntimeError("The target and aperture check did not pass")

    references = pd.read_csv(reference_path)
    target_x_zero = float(settings["astroimagej_center_x_fits"]) - 1.0
    target_y_zero = float(settings["astroimagej_center_y_fits"]) - 1.0

    with fits.open(image_path, memmap=True, checksum=True) as hdul:
        image = hdul[0].data
        target, comparisons = choose_comparison_stars(
            references,
            target_x_zero_based=target_x_zero,
            target_y_zero_based=target_y_zero,
            image_width=image.shape[1],
            image_height=image.shape[0],
            comparison_count=args.comparison_count,
        )
        comparisons["x_fits"] = comparisons["x_pixel"] + 1.0
        comparisons["y_fits"] = comparisons["y_pixel"] + 1.0
        make_field_figure(
            image,
            output_dir / "04_comparison_star_candidates.png",
            target_x_fits=float(settings["astroimagej_center_x_fits"]),
            target_y_fits=float(settings["astroimagej_center_y_fits"]),
            comparisons=comparisons,
        )

    columns = [
        "aperture",
        "x_fits",
        "y_fits",
        "brightness_ratio",
        "width_pixels",
        "target_distance_pixels",
        "edge_distance_pixels",
    ]
    comparisons[columns].to_csv(
        output_dir / "comparison_star_candidates.csv",
        index=False,
        float_format="%.6f",
    )
    write_handoff(
        output_dir,
        settings=settings,
        comparison_count=len(comparisons),
    )
    summary = {
        "status": "ready",
        "target_reference_star": int(target["star"]),
        "comparison_stars": len(comparisons),
        "minimum_brightness_ratio": float(comparisons["brightness_ratio"].min()),
        "maximum_brightness_ratio": float(comparisons["brightness_ratio"].max()),
        "maximum_width_difference_pixels": float(
            comparisons["width_difference_pixels"].max()
        ),
        "minimum_edge_distance_pixels": float(
            comparisons["edge_distance_pixels"].min()
        ),
    }
    (output_dir / "photometry_handoff.json").write_text(
        json.dumps(summary, indent=2) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
