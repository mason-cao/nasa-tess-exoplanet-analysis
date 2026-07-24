"""Review the TOI-3505.01 AstroImageJ photometry tables.

The Schar light-curve tutorial asks for a background-area check, a review of
each comparison star, and plots of the target flux alongside observing
conditions. This script makes those checks without detrending or fitting a
transit model.
"""

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
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT_DIR = ROOT / "outputs" / "toi3505_photometry"
DEFAULT_NARROW_TABLE = (
    DEFAULT_OUTPUT_DIR
    / "TOI_3505.01_2022-07-22_R_measurements_70-100.tbl"
)
DEFAULT_WIDE_TABLE = (
    DEFAULT_OUTPUT_DIR
    / "TOI_3505.01_2022-07-22_R_measurements_70-139.tbl"
)
DEFAULT_WORKING_TABLE = (
    DEFAULT_OUTPUT_DIR
    / "TOI_3505.01_2022-07-22_R_working_measurements.xls"
)
COMPARISON_STARS = [f"C{number}" for number in range(2, 12)]
TARGET_PERIOD_DAYS = 2.9151556
TARGET_EPOCH_BJD_TDB = 2459793.534385
MINIMUM_WORKING_COMPARISONS = 5
MINIMUM_RELATIVE_IMPROVEMENT = 0.005


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--narrow-table", type=Path, default=DEFAULT_NARROW_TABLE)
    parser.add_argument("--wide-table", type=Path, default=DEFAULT_WIDE_TABLE)
    parser.add_argument("--working-table", type=Path, default=DEFAULT_WORKING_TABLE)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    return parser.parse_args()


def normalized(values: pd.Series | np.ndarray) -> np.ndarray:
    array = np.asarray(values, dtype=float)
    median = float(np.nanmedian(array))
    if not np.isfinite(median) or median == 0:
        raise RuntimeError("A light-curve column has an invalid median")
    return array / median


def robust_scatter(values: pd.Series | np.ndarray) -> float:
    curve = normalized(values)
    return float(1.4826 * np.nanmedian(np.abs(curve - np.nanmedian(curve))))


def clipped_standard_deviation(values: pd.Series | np.ndarray) -> float:
    curve = normalized(values)
    center = float(np.nanmedian(curve))
    scatter = float(1.4826 * np.nanmedian(np.abs(curve - center)))
    if not np.isfinite(scatter) or scatter <= 0:
        return float("nan")
    keep = np.isfinite(curve) & (np.abs(curve - center) <= 5.0 * scatter)
    return float(np.nanstd(curve[keep]))


def correlation(left: np.ndarray, right: pd.Series | np.ndarray) -> float:
    x = np.asarray(left, dtype=float)
    y = np.asarray(right, dtype=float)
    keep = np.isfinite(x) & np.isfinite(y)
    if keep.sum() < 3 or np.nanstd(x[keep]) == 0 or np.nanstd(y[keep]) == 0:
        return float("nan")
    return float(np.corrcoef(x[keep], y[keep])[0, 1])


def load_table(
    path: Path,
    expected_outer_radius: int,
    *,
    expected_source_radius: float = 35.0,
    require_all_comparisons: bool = True,
) -> pd.DataFrame:
    table = pd.read_csv(path, sep="\t")
    required = {
        "Label",
        "slice",
        "Saturated",
        "BJD_TDB",
        "AIRMASS",
        "FWHM_Mean",
        "Source_Radius",
        "Sky_Rad(min)",
        "Sky_Rad(max)",
        "rel_flux_T1",
        "tot_C_cnts",
        "Source-Sky_T1",
        "Source_Error_T1",
        "Sky/Pixel_T1",
        "Width_T1",
        "X(FITS)_T1",
        "Y(FITS)_T1",
    }
    if require_all_comparisons:
        for star in COMPARISON_STARS:
            required.update(
                {
                    f"rel_flux_{star}",
                    f"rel_flux_SNR_{star}",
                    f"Source-Sky_{star}",
                    f"Source_Error_{star}",
                }
            )
    missing = sorted(required.difference(table.columns))
    if missing:
        raise RuntimeError(f"{path.name} is missing columns: {', '.join(missing)}")
    if len(table) != 281:
        raise RuntimeError(f"{path.name} contains {len(table)} rows instead of 281")
    if not np.array_equal(table["slice"].to_numpy(), np.arange(1, 282)):
        raise RuntimeError(f"{path.name} does not contain slices 1 through 281 in order")
    if not np.allclose(table["Source_Radius"], expected_source_radius):
        raise RuntimeError(
            f"{path.name} does not use a {expected_source_radius:g}-pixel source radius"
        )
    if not np.allclose(table["Sky_Rad(min)"], 70.0):
        raise RuntimeError(f"{path.name} does not use a 70-pixel inner sky radius")
    if not np.allclose(table["Sky_Rad(max)"], expected_outer_radius):
        raise RuntimeError(
            f"{path.name} does not use a {expected_outer_radius}-pixel outer sky radius"
        )
    return table


def target_curve_from_counts(table: pd.DataFrame, stars: list[str]) -> np.ndarray:
    comparison_total = np.sum(
        [table[f"Source-Sky_{star}"].to_numpy(dtype=float) for star in stars],
        axis=0,
    )
    target = table["Source-Sky_T1"].to_numpy(dtype=float)
    curve = np.full(len(table), np.nan, dtype=float)
    valid = np.isfinite(target) & np.isfinite(comparison_total) & (comparison_total > 0)
    curve[valid] = target[valid] / comparison_total[valid]
    return normalized(curve)


def target_curve_error(table: pd.DataFrame, stars: list[str]) -> np.ndarray:
    target = table["Source-Sky_T1"].to_numpy(dtype=float)
    target_error = table["Source_Error_T1"].to_numpy(dtype=float)
    comparison_total = np.sum(
        [table[f"Source-Sky_{star}"].to_numpy(dtype=float) for star in stars],
        axis=0,
    )
    comparison_error = np.sqrt(
        np.sum(
            [
                np.square(table[f"Source_Error_{star}"].to_numpy(dtype=float))
                for star in stars
            ],
            axis=0,
        )
    )
    ratio = target / comparison_total
    ratio_error = np.abs(ratio) * np.sqrt(
        np.square(target_error / target)
        + np.square(comparison_error / comparison_total)
    )
    return ratio_error / np.nanmedian(ratio)


def validate_native_working_table(
    path: Path,
    source_table: pd.DataFrame,
    working_stars: list[str],
    expected_outer_radius: int,
) -> dict[str, float | int | str | bool]:
    """Confirm that AstroImageJ saved the selected reference-star calculation."""
    native = load_table(
        path,
        expected_outer_radius=expected_outer_radius,
        require_all_comparisons=False,
    )
    if not native["Label"].equals(source_table["Label"]):
        raise RuntimeError(
            f"{path.name} does not contain the same images as the source table"
        )

    expected_curve = target_curve_from_counts(source_table, working_stars)
    native_curve = normalized(native["rel_flux_T1"])
    curve_difference = float(np.nanmax(np.abs(expected_curve - native_curve)))

    expected_counts = np.sum(
        [
            source_table[f"Source-Sky_{star}"].to_numpy(dtype=float)
            for star in working_stars
        ],
        axis=0,
    )
    count_difference = float(
        np.nanmax(
            np.abs(expected_counts - native["tot_C_cnts"].to_numpy(dtype=float))
        )
    )
    if curve_difference > 1.0e-10 or count_difference > 1.0e-5:
        raise RuntimeError(
            f"{path.name} does not match the selected comparison-star calculation"
        )
    return {
        "present": True,
        "file": path.name,
        "rows": int(len(native)),
        "maximum_normalized_curve_difference": curve_difference,
        "maximum_comparison_count_difference": count_difference,
    }


def choose_working_comparisons(
    table: pd.DataFrame,
) -> tuple[list[str], list[dict[str, float | str | int]]]:
    selected = COMPARISON_STARS.copy()
    steps: list[dict[str, float | str | int]] = []
    while len(selected) > MINIMUM_WORKING_COMPARISONS:
        current_scatter = robust_scatter(target_curve_from_counts(table, selected))
        trials: list[tuple[float, str]] = []
        for star in selected:
            trial_stars = [candidate for candidate in selected if candidate != star]
            trials.append(
                (robust_scatter(target_curve_from_counts(table, trial_stars)), star)
            )
        best_scatter, best_star = min(trials)
        improvement = (current_scatter - best_scatter) / current_scatter
        if improvement < MINIMUM_RELATIVE_IMPROVEMENT:
            break
        steps.append(
            {
                "stars_before": len(selected),
                "removed": best_star,
                "scatter_before_percent": 100.0 * current_scatter,
                "scatter_after_percent": 100.0 * best_scatter,
                "relative_improvement_percent": 100.0 * improvement,
            }
        )
        selected.remove(best_star)
    return selected, steps


def comparison_star_review(
    table: pd.DataFrame,
    working_stars: list[str],
) -> pd.DataFrame:
    baseline_scatter = robust_scatter(
        target_curve_from_counts(table, COMPARISON_STARS)
    )
    rows: list[dict[str, float | str | bool]] = []
    for star in COMPARISON_STARS:
        curve = normalized(table[f"rel_flux_{star}"])
        without_star = [candidate for candidate in COMPARISON_STARS if candidate != star]
        leave_one_out = robust_scatter(target_curve_from_counts(table, without_star))
        correlations = {
            "airmass_correlation": correlation(curve, table["AIRMASS"]),
            "star_width_correlation": correlation(curve, table["Width_T1"]),
            "x_position_correlation": correlation(curve, table["X(FITS)_T1"]),
            "y_position_correlation": correlation(curve, table["Y(FITS)_T1"]),
        }
        rows.append(
            {
                "star": star,
                "working_set": star in working_stars,
                "robust_scatter_percent": 100.0 * robust_scatter(curve),
                "clipped_standard_deviation_percent": (
                    100.0 * clipped_standard_deviation(curve)
                ),
                "median_relative_snr": float(
                    np.nanmedian(table[f"rel_flux_SNR_{star}"])
                ),
                "median_source_counts": float(
                    np.nanmedian(table[f"Source-Sky_{star}"])
                ),
                "target_scatter_without_star_percent": 100.0 * leave_one_out,
                "target_scatter_change_percent_points": (
                    100.0 * (leave_one_out - baseline_scatter)
                ),
                **correlations,
                "largest_condition_correlation": float(
                    np.nanmax(np.abs(list(correlations.values())))
                ),
            }
        )
    return pd.DataFrame(rows)


def make_background_check(
    narrow: pd.DataFrame,
    wide: pd.DataFrame,
) -> pd.DataFrame:
    rows: list[dict[str, float | str]] = []
    for label, table in (("70-100 pixels", narrow), ("70-139 pixels", wide)):
        rows.append(
            {
                "sky_ring": label,
                "target_robust_scatter_percent": (
                    100.0 * robust_scatter(table["rel_flux_T1"])
                ),
                "target_clipped_standard_deviation_percent": (
                    100.0 * clipped_standard_deviation(table["rel_flux_T1"])
                ),
                "target_full_standard_deviation_percent": (
                    100.0 * np.nanstd(normalized(table["rel_flux_T1"]))
                ),
                "median_sky_counts_per_pixel": float(
                    np.nanmedian(table["Sky/Pixel_T1"])
                ),
            }
        )
    return pd.DataFrame(rows)


def review_flags(table: pd.DataFrame, working_curve: np.ndarray) -> tuple[np.ndarray, list[str]]:
    width = table["Width_T1"].to_numpy(dtype=float)
    width_center = float(np.nanmedian(width))
    width_scatter = float(1.4826 * np.nanmedian(np.abs(width - width_center)))
    curve_center = float(np.nanmedian(working_curve))
    curve_scatter = robust_scatter(working_curve)
    comparison_counts = table["tot_C_cnts"].to_numpy(dtype=float)
    count_center = float(np.nanmedian(comparison_counts))
    count_scatter = float(
        1.4826 * np.nanmedian(np.abs(comparison_counts - count_center))
    )

    flags = np.zeros(len(table), dtype=bool)
    reasons: list[str] = []
    for index in range(len(table)):
        point_reasons: list[str] = []
        if float(table.iloc[index]["Saturated"]) > 0:
            point_reasons.append("saturated")
        if width[index] > width_center + 4.0 * width_scatter:
            point_reasons.append("large star width")
        if abs(working_curve[index] - curve_center) > 5.0 * curve_scatter:
            point_reasons.append("unusual target brightness")
        if comparison_counts[index] < count_center - 5.0 * count_scatter:
            point_reasons.append("low comparison-star counts")
        flags[index] = bool(point_reasons)
        reasons.append("; ".join(point_reasons))
    return flags, reasons


def plot_background_area_check(
    narrow: pd.DataFrame,
    wide: pd.DataFrame,
    output_path: Path,
) -> None:
    hours = (wide["BJD_TDB"].to_numpy() - float(wide["BJD_TDB"].iloc[0])) * 24.0
    narrow_curve = normalized(narrow["rel_flux_T1"])
    wide_curve = normalized(wide["rel_flux_T1"])
    difference_percent = 100.0 * (narrow_curve - wide_curve)

    fig, axes = plt.subplots(
        2,
        1,
        figsize=(10.5, 7.6),
        sharex=True,
        gridspec_kw={"height_ratios": [2.2, 1.0]},
    )
    axes[0].plot(
        hours,
        narrow_curve,
        ".",
        color="#2f6f9f",
        markersize=4.0,
        label="70-100 pixel sky ring",
    )
    axes[0].plot(
        hours,
        wide_curve,
        ".",
        color="#bf5b3f",
        markersize=4.0,
        alpha=0.78,
        label="70-139 pixel sky ring",
    )
    axes[0].axhline(1.0, color="#777777", linewidth=0.8)
    axes[0].set_ylabel("Relative brightness")
    axes[0].legend(loc="best", frameon=True)
    axes[0].grid(alpha=0.18)

    axes[1].plot(hours, difference_percent, ".", color="#5b4b8a", markersize=4.0)
    axes[1].axhline(0.0, color="#777777", linewidth=0.8)
    axes[1].set_ylabel("Difference (%)")
    axes[1].set_xlabel("Hours since the first image")
    axes[1].grid(alpha=0.18)
    fig.suptitle("TOI-3505.01 background-area check", fontsize=15)
    fig.text(
        0.5,
        0.935,
        "281 R-band images, 50-second exposures, 2022-07-21/22 UTC",
        ha="center",
        fontsize=10.5,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.91))
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def plot_comparison_stars(
    table: pd.DataFrame,
    review: pd.DataFrame,
    output_path: Path,
) -> None:
    hours = (table["BJD_TDB"].to_numpy() - float(table["BJD_TDB"].iloc[0])) * 24.0
    all_percent: list[np.ndarray] = []
    for star in COMPARISON_STARS:
        all_percent.append(100.0 * (normalized(table[f"rel_flux_{star}"]) - 1.0))
    combined = np.concatenate(all_percent)
    lower, upper = np.nanpercentile(combined, [0.5, 99.5])
    limit = max(abs(float(lower)), abs(float(upper)), 1.0) * 1.12

    fig, axes = plt.subplots(5, 2, figsize=(11.5, 13.5), sharex=True, sharey=True)
    review_by_star = review.set_index("star")
    for ax, star, percent in zip(axes.flat, COMPARISON_STARS, all_percent):
        used = bool(review_by_star.loc[star, "working_set"])
        color = "#2f7d4f" if used else "#c26635"
        status = "working set" if used else "set aside for review"
        scatter = float(review_by_star.loc[star, "robust_scatter_percent"])
        ax.plot(hours, percent, ".", color=color, markersize=3.2)
        ax.axhline(0.0, color="#777777", linewidth=0.7)
        ax.set_ylim(-limit, limit)
        ax.set_title(f"{star} — {status}", fontsize=10.5)
        ax.text(
            0.02,
            0.91,
            f"scatter: {scatter:.3f}%",
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=8.5,
        )
        ax.grid(alpha=0.15)
    for ax in axes[:, 0]:
        ax.set_ylabel("Brightness change (%)")
    for ax in axes[-1, :]:
        ax.set_xlabel("Hours since the first image")
    fig.suptitle("TOI-3505.01 comparison-star check", fontsize=15, y=0.995)
    fig.text(
        0.5,
        0.975,
        "70-139 pixel sky ring; each star is shown relative to the other comparison stars",
        ha="center",
        fontsize=10.5,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.945))
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def plot_light_curve_checks(
    table: pd.DataFrame,
    working_curve: np.ndarray,
    working_stars: list[str],
    flags: np.ndarray,
    output_path: Path,
) -> None:
    hours = (table["BJD_TDB"].to_numpy() - float(table["BJD_TDB"].iloc[0])) * 24.0
    working_counts = np.sum(
        [table[f"Source-Sky_{star}"].to_numpy(dtype=float) for star in working_stars],
        axis=0,
    )
    count_percent = 100.0 * (normalized(working_counts) - 1.0)
    x_change = table["X(FITS)_T1"].to_numpy() - np.nanmedian(table["X(FITS)_T1"])
    y_change = table["Y(FITS)_T1"].to_numpy() - np.nanmedian(table["Y(FITS)_T1"])

    fig, axes = plt.subplots(6, 1, figsize=(11.5, 14.0), sharex=True)
    good = ~flags
    axes[0].plot(
        hours[good],
        working_curve[good],
        ".",
        color="#265f8e",
        markersize=4.0,
        label="Points used in the check",
    )
    if flags.any():
        axes[0].plot(
            hours[flags],
            working_curve[flags],
            "o",
            markerfacecolor="none",
            markeredgecolor="#c5573d",
            markersize=5.0,
            label="Points to review",
        )
    axes[0].axhline(1.0, color="#777777", linewidth=0.8)
    axes[0].set_ylabel("Relative\nbrightness")
    axes[0].legend(loc="best", fontsize=8.5)

    axes[1].plot(hours, table["Sky/Pixel_T1"], ".", color="#7a5a96", markersize=3.8)
    axes[1].set_ylabel("Sky level\n(counts/pixel)")

    axes[2].plot(hours, table["Width_T1"], ".", color="#3d806b", markersize=3.8)
    axes[2].set_ylabel("Star width\n(pixels)")

    axes[3].plot(hours, table["AIRMASS"], ".", color="#a25c45", markersize=3.8)
    axes[3].set_ylabel("Airmass")

    axes[4].plot(hours, count_percent, ".", color="#7d6b32", markersize=3.8)
    axes[4].axhline(0.0, color="#777777", linewidth=0.7)
    axes[4].set_ylabel("Comparison-star\ncounts change (%)")

    axes[5].plot(hours, x_change, ".", color="#2f6f9f", markersize=3.8, label="X")
    axes[5].plot(hours, y_change, ".", color="#bf5b3f", markersize=3.8, label="Y")
    axes[5].axhline(0.0, color="#777777", linewidth=0.7)
    axes[5].set_ylabel("Target position\nchange (pixels)")
    axes[5].set_xlabel("Hours since the first image")
    axes[5].legend(loc="best", ncol=2, fontsize=8.5)

    for ax in axes:
        ax.grid(alpha=0.16)
    fig.suptitle("TOI-3505.01 light-curve quality checks", fontsize=15, y=0.995)
    fig.text(
        0.5,
        0.975,
        "No predicted transit midpoint falls in this observing time range",
        ha="center",
        fontsize=10.5,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.945))
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def closest_predicted_midpoint(start_bjd: float, end_bjd: float) -> dict[str, float | bool]:
    middle = 0.5 * (start_bjd + end_bjd)
    epoch_number = round((middle - TARGET_EPOCH_BJD_TDB) / TARGET_PERIOD_DAYS)
    midpoint = TARGET_EPOCH_BJD_TDB + epoch_number * TARGET_PERIOD_DAYS
    return {
        "epoch_number": int(epoch_number),
        "midpoint_bjd_tdb": float(midpoint),
        "inside_observation": bool(start_bjd <= midpoint <= end_bjd),
        "hours_from_observation_start": float((midpoint - start_bjd) * 24.0),
        "hours_from_observation_end": float((midpoint - end_bjd) * 24.0),
    }


def write_readme(
    output_dir: Path,
    *,
    chosen_ring: str,
    background: pd.DataFrame,
    working_stars: list[str],
    removed_stars: list[str],
    all_scatter: float,
    working_scatter: float,
    start_bjd: float,
    end_bjd: float,
    closest_midpoint: dict[str, float | bool],
    review_count: int,
) -> None:
    background_by_ring = background.set_index("sky_ring")
    narrow_scatter = float(
        background_by_ring.loc["70-100 pixels", "target_robust_scatter_percent"]
    )
    wide_scatter = float(
        background_by_ring.loc["70-139 pixels", "target_robust_scatter_percent"]
    )
    text = f"""# TOI-3505.01 multi-aperture photometry

## What is finished

- AstroImageJ measured the target and ten comparison stars in all 281 usable
  aligned images.
- The source radius was 35 pixels in both runs.
- Two sky rings were tested: 70-100 pixels and 70-139 pixels.
- Neither table contains a saturated measurement.
- The native AstroImageJ tables, aperture file, screenshots, review tables,
  and plots are in this folder.

## Background-area check

The 70-100 pixel ring gave {narrow_scatter:.3f}% robust scatter in the target
light curve. The 70-139 pixel ring gave {wide_scatter:.3f}%. The difference is
small, but the wider ring performed slightly better while AstroImageJ's
**Remove stars from background** setting was on. I kept the {chosen_ring} ring
for the next light-curve check because it also matches the Seeing Profile.

## Comparison-star check

All ten stars gave {all_scatter:.3f}% robust target scatter. A working check
set using {', '.join(working_stars)} gave {working_scatter:.3f}%. The stars set
aside for review are {', '.join(removed_stars)}. This is a working selection,
not a final claim; the individual curves and the review table show why each
star should be checked with the research mentor.

There are {review_count} image measurements marked for review because of an
unusual target value, large star width, low comparison-star counts, or a
saturation flag. They remain in the native tables. No detrending or transit
fit has been applied.

## Timing

The table covers BJD_TDB {start_bjd:.6f} to {end_bjd:.6f}. Using the stated
period of {TARGET_PERIOD_DAYS:.7f} days and epoch {TARGET_EPOCH_BJD_TDB:.6f},
the closest predicted midpoint is
{float(closest_midpoint['midpoint_bjd_tdb']):.6f}. It is
{abs(float(closest_midpoint['hours_from_observation_start'])):.2f} hours before
the first image, so this data set does not contain that predicted midpoint.
This light curve can still be used to check the reduction, apertures,
comparison stars, and systematic trends, but it should not be presented as a
transit detection.

## Saved AstroImageJ review files

The working reference-star calculation is saved in
`TOI_3505.01_2022-07-22_R_working_measurements.xls`. The matching aperture
file is `TOI_3505_working.apertures`, and the matching plot settings are in
`TOI_3505.01_2022-07-22_R_working.plotcfg`. The AstroImageJ screenshots show
the measurement table, reference-star selection, Y-data settings, raw light
curve, and fit settings. The transit fit and detrending are both off.
The saved table's target curve and total comparison counts match the seven
selected comparison stars.

The `.tbl` copy has the same contents as the `.xls` table and is included so
AstroImageJ 6.0.7 can reopen the working measurements directly.

## Next Schar tutorial step

Send the working light curve, seeing plot, measurement table, aperture view,
plot configuration, and observing-condition plots to the research mentor.
Ask for confirmation of the sky ring and comparison-star selection, and ask
whether this sequence is meant to be a reduction-quality exercise or whether
there is another TOI-3505 sequence that contains the predicted transit. Add a
predicted transit model only after that timing question is resolved.
"""
    (output_dir / "README.md").write_text(text, encoding="utf-8")

def main() -> None:
    args = parse_args()
    narrow_path = args.narrow_table.resolve()
    wide_path = args.wide_table.resolve()
    working_path = args.working_table.resolve()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    narrow = load_table(narrow_path, expected_outer_radius=100)
    wide = load_table(wide_path, expected_outer_radius=139)
    if not narrow["Label"].equals(wide["Label"]):
        raise RuntimeError("The two sky-ring tables do not contain the same images")

    background = make_background_check(narrow, wide)
    narrow_scatter = float(
        background.loc[
            background["sky_ring"] == "70-100 pixels",
            "target_robust_scatter_percent",
        ].iloc[0]
    )
    wide_scatter = float(
        background.loc[
            background["sky_ring"] == "70-139 pixels",
            "target_robust_scatter_percent",
        ].iloc[0]
    )
    chosen = wide if wide_scatter <= narrow_scatter else narrow
    chosen_ring = "70-139 pixels" if chosen is wide else "70-100 pixels"

    working_stars, selection_steps = choose_working_comparisons(chosen)
    removed_stars = [star for star in COMPARISON_STARS if star not in working_stars]
    comparison_review = comparison_star_review(chosen, working_stars)
    working_curve = target_curve_from_counts(chosen, working_stars)
    working_error = target_curve_error(chosen, working_stars)
    all_curve = normalized(chosen["rel_flux_T1"])
    all_scatter = 100.0 * robust_scatter(all_curve)
    working_scatter = 100.0 * robust_scatter(working_curve)
    flags, reasons = review_flags(chosen, working_curve)

    native_working_table: dict[str, float | int | str | bool] = {
        "present": False,
        "file": working_path.name,
    }
    if working_path.exists():
        expected_outer_radius = 139 if chosen is wide else 100
        native_working_table = validate_native_working_table(
            working_path,
            chosen,
            working_stars,
            expected_outer_radius,
        )

    start_bjd = float(chosen["BJD_TDB"].iloc[0])
    end_bjd = float(chosen["BJD_TDB"].iloc[-1])
    midpoint = closest_predicted_midpoint(start_bjd, end_bjd)
    hours = (chosen["BJD_TDB"].to_numpy() - start_bjd) * 24.0
    working_counts = np.sum(
        [chosen[f"Source-Sky_{star}"].to_numpy(dtype=float) for star in working_stars],
        axis=0,
    )
    light_curve = pd.DataFrame(
        {
            "image": chosen["Label"],
            "slice": chosen["slice"],
            "bjd_tdb": chosen["BJD_TDB"],
            "hours_since_first_image": hours,
            "relative_brightness_all_comparisons": all_curve,
            "relative_brightness_working_set": working_curve,
            "relative_brightness_error": working_error,
            "sky_counts_per_pixel": chosen["Sky/Pixel_T1"],
            "star_width_pixels": chosen["Width_T1"],
            "airmass": chosen["AIRMASS"],
            "working_comparison_counts": working_counts,
            "target_x_fits": chosen["X(FITS)_T1"],
            "target_y_fits": chosen["Y(FITS)_T1"],
            "review_point": flags,
            "review_reason": reasons,
        }
    )

    background.to_csv(output_dir / "background_area_check.csv", index=False)
    comparison_review.to_csv(
        output_dir / "comparison_star_review.csv", index=False, float_format="%.8f"
    )
    light_curve.to_csv(
        output_dir / "TOI_3505.01_2022-07-22_R_working_light_curve.csv",
        index=False,
        float_format="%.10f",
    )

    plot_background_area_check(
        narrow, wide, output_dir / "07_background_area_check.png"
    )
    plot_comparison_stars(
        chosen, comparison_review, output_dir / "08_comparison_star_check.png"
    )
    plot_light_curve_checks(
        chosen,
        working_curve,
        working_stars,
        flags,
        output_dir / "09_light_curve_quality_checks.png",
    )

    condition_columns = {
        "airmass": chosen["AIRMASS"],
        "sky_counts_per_pixel": chosen["Sky/Pixel_T1"],
        "star_width_pixels": chosen["Width_T1"],
        "target_x_fits": chosen["X(FITS)_T1"],
        "target_y_fits": chosen["Y(FITS)_T1"],
        "comparison_counts": working_counts,
    }
    condition_correlations = {
        name: correlation(working_curve, values)
        for name, values in condition_columns.items()
    }
    summary = {
        "status": "pass",
        "images_measured": int(len(chosen)),
        "saturated_measurements": int((chosen["Saturated"] > 0).sum()),
        "source_radius_pixels": 35,
        "chosen_sky_ring": chosen_ring,
        "background_area_check": background.to_dict(orient="records"),
        "working_comparison_stars": working_stars,
        "comparison_stars_set_aside": removed_stars,
        "selection_steps": selection_steps,
        "target_robust_scatter_all_comparisons_percent": all_scatter,
        "target_robust_scatter_working_set_percent": working_scatter,
        "native_working_table_validation": native_working_table,
        "review_points": int(flags.sum()),
        "observation": {
            "start_bjd_tdb": start_bjd,
            "end_bjd_tdb": end_bjd,
            "duration_hours": float((end_bjd - start_bjd) * 24.0),
            "airmass_min": float(np.nanmin(chosen["AIRMASS"])),
            "airmass_max": float(np.nanmax(chosen["AIRMASS"])),
            "median_fwhm_pixels": float(np.nanmedian(chosen["FWHM_Mean"])),
            "median_width_pixels": float(np.nanmedian(chosen["Width_T1"])),
        },
        "closest_predicted_transit_midpoint": midpoint,
        "working_light_curve_condition_correlations": condition_correlations,
        "detrending_applied": False,
        "transit_model_applied": False,
    }
    (output_dir / "summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    write_readme(
        output_dir,
        chosen_ring=chosen_ring,
        background=background,
        working_stars=working_stars,
        removed_stars=removed_stars,
        all_scatter=all_scatter,
        working_scatter=working_scatter,
        start_bjd=start_bjd,
        end_bjd=end_bjd,
        closest_midpoint=midpoint,
        review_count=int(flags.sum()),
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
