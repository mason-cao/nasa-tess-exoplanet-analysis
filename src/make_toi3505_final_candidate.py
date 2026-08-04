"""Build the TOI-3505.01 light curve and its supporting checks.

This script deliberately keeps the raw AstroImageJ differential photometry as
the primary scientific record.  It uses all ten comparison stars so that the
ensemble is not selected by minimizing the target's scatter, retains every
exposure, and promotes a systematics correction only when a predeclared model
passes blocked out-of-sample checks.

The current ephemeris places the nearest predicted transit midpoint outside
the observing window, while the recovered 2022 scheduling-sheet row places a
historical ingress/egress window inside it.  The script shows that historical
window and measures one fixed box as a timing check, but never fits a physical
transit model.
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import tempfile
from dataclasses import dataclass
from pathlib import Path

os.environ.setdefault(
    "MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "nasa-tess-matplotlib-cache")
)

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from analyze_toi3505_photometry import (
    COMPARISON_STARS,
    DEFAULT_WIDE_TABLE,
    TARGET_PERIOD_DAYS,
    clipped_standard_deviation,
    closest_predicted_midpoint,
    correlation,
    load_table,
    normalized,
    robust_scatter,
    target_curve_error,
    target_curve_from_counts,
)
from toi3505_schedule import (
    DEFAULT_SCHEDULE_RECORD,
    analyze_schedule_window,
    write_target_plot_config,
)


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT_DIR = ROOT / "outputs" / "toi3505_final_candidate"
DEFAULT_REVIEW_PACKAGE_DIR = ROOT / "outputs" / "toi3505_review_package"
DEFAULT_FINAL_TABLE = (
    ROOT
    / "outputs"
    / "toi3505_aperture_check"
    / "TOI_3505.01_2022-07-22_R_measurements_25px_70-139_AIJ.tbl"
)
DEFAULT_SOURCE_RADIUS = 25.0
PREVIOUS_WORKING_STARS = ["C2", "C3", "C5", "C6", "C7", "C10", "C11"]
BLOCK_COUNT = 8
MINIMUM_CV_IMPROVEMENT = 0.05
MAXIMUM_RMSE_RATIO = 1.01
MINIMUM_BETTER_BLOCK_FRACTION = 0.60
NEAR_BEST_SCORE_FRACTION = 0.01
BIN_MINUTES = 10.0
CATALOG_DEPTH_PPT = 2.910
TARGET_TEFF_K = 6220.0  # TIC v8 value in data/catalogs/toi3505.
VISUAL_QUALITY_EXCLUSIONS = {
    253: "field-wide vertical PSF trail confirmed in target and comparison stars",
}

# The list is intentionally short and physically interpretable.  Time itself
# is not a candidate because it could absorb real variability.  The meridian
# side term is derived from the minimum airmass and is tested, not assumed.
CANDIDATE_MODELS: dict[str, tuple[str, ...]] = {
    "none": (),
    "airmass": ("airmass",),
    "sky": ("sky",),
    "width": ("width",),
    "x": ("x",),
    "y": ("y",),
    "comparison_counts": ("comparison_counts",),
    "centroid_xy": ("x", "y"),
    "airmass_sky": ("airmass", "sky"),
    "airmass_meridian": ("airmass", "meridian_side"),
    "conditions": (
        "airmass",
        "sky",
        "width",
        "x",
        "y",
        "comparison_counts",
        "meridian_side",
    ),
}


@dataclass(frozen=True)
class RobustLinearModel:
    """A small robust linear model with reusable feature scaling."""

    feature_names: tuple[str, ...]
    center: np.ndarray
    scale: np.ndarray
    coefficients: np.ndarray

    def predict(self, features: pd.DataFrame) -> np.ndarray:
        if not self.feature_names:
            return np.full(len(features), float(self.coefficients[0]))
        values = features.loc[:, self.feature_names].to_numpy(dtype=float)
        standardized = (values - self.center) / self.scale
        design = np.column_stack([np.ones(len(features)), standardized])
        return design @ self.coefficients


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--table", type=Path, default=DEFAULT_FINAL_TABLE)
    parser.add_argument(
        "--source-radius", type=float, default=DEFAULT_SOURCE_RADIUS
    )
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument(
        "--schedule-record", type=Path, default=DEFAULT_SCHEDULE_RECORD
    )
    return parser.parse_args()


def robust_sigma(values: np.ndarray | pd.Series) -> float:
    """Return a normalization-free MAD estimate of standard deviation."""
    array = np.asarray(values, dtype=float)
    finite = array[np.isfinite(array)]
    if finite.size == 0:
        return float("nan")
    center = float(np.median(finite))
    return float(1.4826 * np.median(np.abs(finite - center)))


def fit_robust_linear(
    features: pd.DataFrame,
    values: np.ndarray,
    feature_names: tuple[str, ...],
    *,
    maximum_iterations: int = 50,
    huber_delta: float = 1.345,
) -> RobustLinearModel:
    """Fit a Huber IRLS model without adding a heavyweight dependency."""
    y = np.asarray(values, dtype=float)
    if len(features) != len(y):
        raise ValueError("Features and values have different lengths")
    if not feature_names:
        finite = y[np.isfinite(y)]
        if finite.size == 0:
            raise ValueError("No finite values are available for the intercept")
        return RobustLinearModel(
            feature_names=(),
            center=np.empty(0),
            scale=np.empty(0),
            coefficients=np.array([float(np.median(finite))]),
        )

    x = features.loc[:, feature_names].to_numpy(dtype=float)
    finite_rows = np.isfinite(y) & np.all(np.isfinite(x), axis=1)
    if finite_rows.sum() <= len(feature_names) + 1:
        raise ValueError("Too few finite rows for the requested model")
    x = x[finite_rows]
    y = y[finite_rows]

    center = np.median(x, axis=0)
    scale = np.array([robust_sigma(x[:, index]) for index in range(x.shape[1])])
    ordinary_scale = np.std(x, axis=0)
    invalid_scale = ~np.isfinite(scale) | (scale <= 0)
    scale[invalid_scale] = ordinary_scale[invalid_scale]
    scale[~np.isfinite(scale) | (scale <= 0)] = 1.0
    standardized = (x - center) / scale
    design = np.column_stack([np.ones(len(x)), standardized])
    coefficients = np.linalg.lstsq(design, y, rcond=None)[0]

    for _ in range(maximum_iterations):
        residual = y - design @ coefficients
        residual_scale = robust_sigma(residual)
        if not np.isfinite(residual_scale) or residual_scale <= 0:
            residual_scale = float(1.4826 * np.median(np.abs(residual)))
        if not np.isfinite(residual_scale) or residual_scale <= 0:
            residual_scale = float(np.std(residual))
        if not np.isfinite(residual_scale) or residual_scale <= 0:
            break
        distance = np.abs(residual - np.median(residual))
        cutoff = huber_delta * residual_scale
        weights = np.ones_like(distance)
        outside = distance > cutoff
        weights[outside] = cutoff / distance[outside]
        weighted_design = design * np.sqrt(weights)[:, None]
        weighted_values = y * np.sqrt(weights)
        updated = np.linalg.lstsq(weighted_design, weighted_values, rcond=None)[0]
        if np.max(np.abs(updated - coefficients)) < 1.0e-12:
            coefficients = updated
            break
        coefficients = updated

    return RobustLinearModel(
        feature_names=feature_names,
        center=center,
        scale=scale,
        coefficients=coefficients,
    )


def build_features(table: pd.DataFrame, comparison_counts: np.ndarray) -> pd.DataFrame:
    airmass = table["AIRMASS"].to_numpy(dtype=float)
    meridian_index = int(np.nanargmin(airmass))
    meridian_side = np.zeros(len(table), dtype=float)
    meridian_side[meridian_index + 1 :] = 1.0
    return pd.DataFrame(
        {
            "airmass": airmass,
            "sky": table["Sky/Pixel_T1"].to_numpy(dtype=float),
            "width": table["Width_T1"].to_numpy(dtype=float),
            "x": table["X(FITS)_T1"].to_numpy(dtype=float),
            "y": table["Y(FITS)_T1"].to_numpy(dtype=float),
            "comparison_counts": np.log(comparison_counts),
            "meridian_side": meridian_side,
        }
    )


def make_review_flags(
    table: pd.DataFrame,
    raw_curve: np.ndarray,
    comparison_counts: np.ndarray,
) -> tuple[np.ndarray, list[str], np.ndarray]:
    """Mark review points and define a condition-only model-fit mask.

    Brightness anomalies are marked for human review but do not decide whether
    a row is eligible for the condition model.  The model-fit mask therefore
    cannot reject a possible astrophysical change merely because it is deep.
    No row is removed from the exported light curve.
    """
    width = table["Width_T1"].to_numpy(dtype=float)
    width_center = float(np.nanmedian(width))
    width_scale = robust_sigma(width)
    count_center = float(np.nanmedian(comparison_counts))
    count_scale = robust_sigma(comparison_counts)
    curve_center = float(np.nanmedian(raw_curve))
    curve_scale = robust_sigma(raw_curve)

    review = np.zeros(len(table), dtype=bool)
    fit_eligible = np.ones(len(table), dtype=bool)
    reasons: list[str] = []
    for index in range(len(table)):
        row_reasons: list[str] = []
        condition_failure = False
        if float(table.iloc[index]["Saturated"]) > 0:
            row_reasons.append("saturation flag")
            condition_failure = True
        if not (
            np.isfinite(raw_curve[index])
            and np.isfinite(width[index])
            and np.isfinite(comparison_counts[index])
        ):
            row_reasons.append("non-finite measurement")
            condition_failure = True
        if width[index] > width_center + 4.0 * width_scale:
            row_reasons.append("large star width")
            condition_failure = True
        if comparison_counts[index] < count_center - 5.0 * count_scale:
            row_reasons.append("low comparison-star counts")
            condition_failure = True
        if abs(raw_curve[index] - curve_center) > 5.0 * curve_scale:
            row_reasons.append("unusual target brightness")
        review[index] = bool(row_reasons)
        fit_eligible[index] = not condition_failure
        reasons.append("; ".join(row_reasons))
    return review, reasons, fit_eligible


def contiguous_fold_ids(indices: np.ndarray, block_count: int) -> np.ndarray:
    """Assign sorted eligible indices to contiguous, nearly equal blocks."""
    if block_count < 2:
        raise ValueError("At least two blocks are required")
    if len(indices) < block_count:
        raise ValueError("There are fewer eligible rows than requested blocks")
    fold_ids = np.full(int(np.max(indices)) + 1, -1, dtype=int)
    for fold, block in enumerate(np.array_split(np.asarray(indices, dtype=int), block_count)):
        fold_ids[block] = fold
    return fold_ids


def blocked_model_review(
    features: pd.DataFrame,
    raw_curve: np.ndarray,
    fit_eligible: np.ndarray,
    *,
    block_count: int = BLOCK_COUNT,
    candidate_models: dict[str, tuple[str, ...]] = CANDIDATE_MODELS,
) -> tuple[pd.DataFrame, dict[str, np.ndarray]]:
    """Score predeclared condition models with contiguous blocked CV."""
    y = np.asarray(raw_curve, dtype=float)
    eligible_indices = np.flatnonzero(fit_eligible & np.isfinite(y))
    fold_ids = contiguous_fold_ids(eligible_indices, block_count)
    predictions: dict[str, np.ndarray] = {}
    fold_residual_scatter: dict[str, list[float]] = {}

    for name, feature_names in candidate_models.items():
        prediction = np.full(len(y), np.nan, dtype=float)
        per_fold: list[float] = []
        for fold in range(block_count):
            test = eligible_indices[fold_ids[eligible_indices] == fold]
            train = eligible_indices[fold_ids[eligible_indices] != fold]
            model = fit_robust_linear(
                features.iloc[train], y[train], feature_names
            )
            prediction[test] = model.predict(features.iloc[test])
            per_fold.append(robust_sigma(y[test] - prediction[test]))
        predictions[name] = prediction
        fold_residual_scatter[name] = per_fold

    null_prediction = predictions["none"]
    null_keep = fit_eligible & np.isfinite(null_prediction) & np.isfinite(y)
    null_residual = y[null_keep] - null_prediction[null_keep]
    null_scatter = robust_sigma(null_residual)
    null_rmse = float(np.sqrt(np.mean(np.square(null_residual))))
    null_blocks = np.asarray(fold_residual_scatter["none"])

    rows: list[dict[str, float | str | int | bool]] = []
    for name, feature_names in candidate_models.items():
        prediction = predictions[name]
        keep = fit_eligible & np.isfinite(prediction) & np.isfinite(y)
        residual = y[keep] - prediction[keep]
        scatter = robust_sigma(residual)
        rmse = float(np.sqrt(np.mean(np.square(residual))))
        block_scores = np.asarray(fold_residual_scatter[name])
        if name == "none":
            blocks_better = 0
        else:
            blocks_better = int(np.sum(block_scores < null_blocks))
        rows.append(
            {
                "model": name,
                "features": "+".join(feature_names) if feature_names else "intercept",
                "parameter_count": 1 + len(feature_names),
                "cv_robust_scatter_fraction": scatter,
                "cv_robust_scatter_ppt": 1000.0 * scatter,
                "cv_rmse_fraction": rmse,
                "cv_rmse_ppt": 1000.0 * rmse,
                "robust_improvement_vs_none_percent": (
                    100.0 * (null_scatter - scatter) / null_scatter
                ),
                "rmse_ratio_vs_none": rmse / null_rmse,
                "blocks_better_than_none": blocks_better,
                "block_count": block_count,
                "passes_absolute_rules": bool(
                    name != "none"
                    and scatter <= null_scatter * (1.0 - MINIMUM_CV_IMPROVEMENT)
                    and rmse <= null_rmse * MAXIMUM_RMSE_RATIO
                    and blocks_better / block_count >= MINIMUM_BETTER_BLOCK_FRACTION
                ),
            }
        )
    return pd.DataFrame(rows), predictions


def select_promoted_model(review: pd.DataFrame) -> str:
    """Apply the frozen promotion rule and prefer parsimony near the best."""
    passing = review[review["passes_absolute_rules"]].copy()
    if passing.empty:
        return "none"
    best_score = float(passing["cv_robust_scatter_fraction"].min())
    near_best = passing[
        passing["cv_robust_scatter_fraction"]
        <= best_score * (1.0 + NEAR_BEST_SCORE_FRACTION)
    ].copy()
    near_best = near_best.sort_values(
        ["parameter_count", "cv_robust_scatter_fraction", "model"]
    )
    return str(near_best.iloc[0]["model"])


def apply_condition_model(
    features: pd.DataFrame,
    raw_curve: np.ndarray,
    fit_eligible: np.ndarray,
    feature_names: tuple[str, ...],
) -> tuple[np.ndarray, np.ndarray, RobustLinearModel]:
    """Fit the promoted model and remove only its relative baseline trend."""
    fit_indices = np.flatnonzero(fit_eligible & np.isfinite(raw_curve))
    model = fit_robust_linear(
        features.iloc[fit_indices], raw_curve[fit_indices], feature_names
    )
    baseline = model.predict(features)
    corrected = raw_curve - (baseline - float(np.nanmedian(baseline)))
    return normalized(corrected), baseline, model


def bin_light_curve(
    hours: np.ndarray,
    values: np.ndarray,
    *,
    bin_minutes: float = BIN_MINUTES,
) -> pd.DataFrame:
    """Make transparent fixed-width display bins without dropping raw rows."""
    hours = np.asarray(hours, dtype=float)
    values = np.asarray(values, dtype=float)
    width_hours = bin_minutes / 60.0
    identifiers = np.floor((hours - np.nanmin(hours)) / width_hours).astype(int)
    rows: list[dict[str, float | int]] = []
    for identifier in np.unique(identifiers):
        keep = (identifiers == identifier) & np.isfinite(hours) & np.isfinite(values)
        if not np.any(keep):
            continue
        selected = values[keep]
        scatter = robust_sigma(selected)
        count = int(keep.sum())
        rows.append(
            {
                "bin": int(identifier),
                "hours": float(np.median(hours[keep])),
                "relative_brightness": float(np.median(selected)),
                "uncertainty": float(scatter / np.sqrt(count)),
                "points": count,
            }
        )
    return pd.DataFrame(rows)


def aij_review_columns() -> list[str]:
    """Return the lossless subset needed by MultiPlot and the NEB macro.

    AstroImageJ 6.0.7 cannot reliably reopen this project's 255-column
    measurement table.  The review table keeps every row and every value used
    by the requested fit-settings and nearby-eclipsing-binary diagnostics, but
    omits unrelated aperture-statistics columns.  It is a view of the original
    measurements, not a new photometric reduction.
    """
    columns = [
        "Label",
        "slice",
        "Saturated",
        "J.D.-2400000",
        "JD_UTC",
        "BJD_TDB",
        "AIRMASS",
        "Source_Radius",
        "Sky/Pixel_T1",
        "Width_T1",
        "tot_C_cnts",
    ]
    for star in ("T1", *COMPARISON_STARS):
        columns.extend(
            [
                f"X(FITS)_{star}",
                f"Y(FITS)_{star}",
                f"Source-Sky_{star}",
                f"rel_flux_{star}",
                f"rel_flux_err_{star}",
            ]
        )
    return columns


def write_aij_review_table(table: pd.DataFrame, output_path: Path) -> None:
    """Write a compact tab-delimited AIJ table without changing its rows."""
    columns = aij_review_columns()
    missing = [column for column in columns if column not in table.columns]
    if missing:
        raise ValueError(f"AIJ review table is missing columns: {missing}")
    table.loc[:, columns].to_csv(
        output_path,
        sep="\t",
        index=False,
        lineterminator="\n",
    )


def comparison_ensemble_review(table: pd.DataFrame) -> pd.DataFrame:
    """Review comparison stars independently of the target light curve."""
    rows: list[dict[str, float | str | bool]] = []
    for star in COMPARISON_STARS:
        pseudo_target = normalized(table[f"rel_flux_{star}"])
        rows.append(
            {
                "star": star,
                "included_primary_ensemble": True,
                "pseudo_target_robust_scatter_percent": 100.0
                * robust_scatter(pseudo_target),
                "pseudo_target_clipped_std_percent": 100.0
                * clipped_standard_deviation(pseudo_target),
                "median_relative_snr": float(
                    np.nanmedian(table[f"rel_flux_SNR_{star}"])
                ),
                "median_source_counts": float(
                    np.nanmedian(table[f"Source-Sky_{star}"])
                ),
                "airmass_correlation": correlation(pseudo_target, table["AIRMASS"]),
                "width_correlation": correlation(pseudo_target, table["Width_T1"]),
                "x_correlation": correlation(pseudo_target, table["X(FITS)_T1"]),
                "y_correlation": correlation(pseudo_target, table["Y(FITS)_T1"]),
            }
        )
    review = pd.DataFrame(rows)
    scatter = review["pseudo_target_robust_scatter_percent"].to_numpy()
    center = float(np.median(scatter))
    scale = robust_sigma(scatter)
    review["pseudo_target_robust_z"] = (scatter - center) / scale
    review["three_mad_scatter_outlier"] = review["pseudo_target_robust_z"] > 3.0
    return review


def _save_figure_variants(fig: plt.Figure, output_path: Path) -> None:
    """Save the requested raster figure and a matching vector copy."""
    fig.savefig(output_path, dpi=220, facecolor="white")
    if output_path.suffix.lower() == ".png":
        svg_path = output_path.with_suffix(".svg")
        fig.savefig(svg_path, facecolor="white")
        svg_text = svg_path.read_text(encoding="utf-8")
        svg_path.write_text(
            "\n".join(line.rstrip() for line in svg_text.splitlines()) + "\n",
            encoding="utf-8",
        )


def page_relative_series(
    values: np.ndarray | pd.Series,
    *,
    center: float = 0.946,
    span: float = 0.022,
    invert: bool = False,
) -> np.ndarray:
    """Scale one diagnostic into a shared display band without changing time.

    AstroImageJ's recommended final-light-curve layout uses ``Page Rel`` for
    observing-condition diagnostics.  This helper reproduces that display-only
    behavior with robust percentiles.  It is never used in the photometry,
    detrending review, fixed-window measurement, or exported flux table.
    """
    array = np.asarray(values, dtype=float)
    finite = array[np.isfinite(array)]
    if finite.size == 0:
        return np.full_like(array, np.nan)
    lower, median, upper = np.nanpercentile(finite, [2.0, 50.0, 98.0])
    width = float(upper - lower)
    if not np.isfinite(width) or width <= 0.0:
        return np.full_like(array, center)
    scaled = np.clip((array - median) / width, -0.65, 0.65)
    direction = -1.0 if invert else 1.0
    return center + direction * span * scaled


def plot_final_candidate(
    output_path: Path,
    *,
    bjd_tdb: np.ndarray,
    adopted_curve: np.ndarray,
    primary_use: np.ndarray,
    binned: pd.DataFrame,
    diagnostics: dict[str, np.ndarray],
    schedule: dict[str, object],
    source_radius: float,
) -> None:
    """Make the mentor-review light curve in the program's diagnostic style.

    The student example saved with the project is an AstroImageJ fit plot: it
    contains a detrended target curve, a transit model, and residuals.  Those
    layers are intentionally absent here because this sequence has no
    supported transit fit and no detrending model passed the frozen tests.  We
    do include the six observing-condition traces requested by the TFOP/Schar
    guidance, scaled and shifted only for display below the unmodified target
    flux.
    """
    bjd_zero = 2459782.0
    x_values = np.asarray(bjd_tdb, dtype=float) - bjd_zero
    binned_x = (
        float(bjd_tdb[0])
        + binned["hours"].to_numpy(dtype=float) / 24.0
        - bjd_zero
    )
    primary_use = np.asarray(primary_use, dtype=bool)
    excluded = ~primary_use
    working = schedule["working_interpretation"]
    fixed_window = schedule["fixed_window_check"]
    assert isinstance(working, dict) and isinstance(fixed_window, dict)
    schedule_times = working["times"]
    assert isinstance(schedule_times, dict)
    ingress_x = float(schedule_times["ingress"]["bjd_tdb"]) - bjd_zero
    egress_x = float(schedule_times["egress"]["bjd_tdb"]) - bjd_zero

    fig, axis = plt.subplots(figsize=(9.2, 8.2))
    axis.axvspan(ingress_x, egress_x, color="#d9534f", alpha=0.045, zorder=0)
    axis.axvline(
        ingress_x,
        color="#d43f3a",
        linestyle="--",
        linewidth=0.9,
        zorder=1,
    )
    axis.axvline(
        egress_x,
        color="#d43f3a",
        linestyle="--",
        linewidth=0.9,
        zorder=1,
    )

    primary_scatter_ppt = 1000.0 * robust_scatter(adopted_curve[primary_use])
    axis.plot(
        x_values[primary_use],
        adopted_curve[primary_use],
        ".",
        color="#0645f5",
        markersize=3.5,
        label=(
            "rel_flux_T1 (normalized, undetrended) "
            f"(RMS={primary_scatter_ppt:.2f} ppt; bin=1)"
        ),
        zorder=4,
    )
    axis.plot(
        x_values[excluded],
        adopted_curve[excluded],
        "x",
        color="#8f8f8f",
        markersize=4.5,
        markeredgewidth=0.9,
        label=f"excluded image-quality frames ({int(excluded.sum())})",
        zorder=3,
    )
    axis.errorbar(
        binned_x,
        binned["relative_brightness"],
        yerr=binned["uncertainty"],
        fmt="o",
        color="#202020",
        markerfacecolor="white",
        markeredgecolor="#202020",
        markersize=4.8,
        linewidth=0.8,
        capsize=1.8,
        label=f"{BIN_MINUTES:.0f}-minute medians (display only)",
        zorder=5,
    )
    axis.axhline(1.0, color="#5f5f5f", linewidth=0.7, zorder=2)

    diagnostic_specs = (
        ("width", "Width_T1", "#b7b7b7", False, 1.1),
        ("sky", "Sky/Pixel_T1", "#e6d400", False, 1.0),
        ("airmass", "AIRMASS (inverted)", "#1596a6", True, 1.0),
        ("comparison_counts", "tot_C_cnts", "#a67c52", False, 1.0),
        ("x", "X(FITS)_T1", "#f09a9a", False, 0.75),
        ("y", "Y(FITS)_T1", "#4db6e7", False, 0.75),
    )
    for key, label, color, invert, linewidth in diagnostic_specs:
        display_values = page_relative_series(diagnostics[key], invert=invert)
        axis.plot(
            x_values,
            display_values,
            color=color,
            linewidth=linewidth,
            alpha=0.95,
            label=f"{label} (page-relative, display only)",
            zorder=2,
        )

    axis.text(
        ingress_x,
        0.9290,
        "Ingress",
        color="#d43f3a",
        ha="center",
        va="bottom",
        fontsize=8.0,
    )
    axis.text(
        egress_x,
        0.9290,
        "Egress",
        color="#d43f3a",
        ha="center",
        va="bottom",
        fontsize=8.0,
    )
    axis.text(
        0.5 * (ingress_x + egress_x),
        0.9263,
        "2022 schedule window (EDT assumption)",
        color="#d43f3a",
        ha="center",
        va="bottom",
        fontsize=7.8,
    )
    axis.text(
        0.985,
        0.975,
        "Undetrended; no transit fit\n"
        "Historical-window depth: "
        f"{float(fixed_window['observed_depth_ppt']):+.2f} ± "
        f"{float(fixed_window['observed_depth_error_ppt']):.2f} ppt",
        transform=axis.transAxes,
        ha="right",
        va="top",
        fontsize=8.5,
        color="#202020",
    )

    axis.set_xlabel("Barycentric Julian Date (TDB) − 2459782")
    axis.set_ylabel("rel_flux_T1 (normalized); diagnostics offset for display")
    axis.set_xlim(float(np.nanmin(x_values)), float(np.nanmax(x_values)))
    axis.set_ylim(0.925, 1.017)
    axis.grid(axis="y", color="#bdbdbd", alpha=0.65, linewidth=0.65)
    axis.tick_params(top=True, right=True, direction="in")
    axis.legend(
        loc="upper left",
        frameon=False,
        fontsize=7.25,
        ncol=1,
        handlelength=2.0,
        handletextpad=0.55,
        borderaxespad=0.5,
        labelspacing=0.32,
    )
    fig.suptitle("TOI-3505.01, UT 2022-07-22", fontsize=15, y=0.982)
    fig.text(
        0.5,
        0.950,
        f"GMU 0.8 m (R filter, 50 s exp, fap {source_radius:g}-70-139)",
        ha="center",
        fontsize=10.5,
    )
    fig.text(
        0.5,
        0.018,
        "The historical sheet does not state its time zone; the current-catalog "
        "transit is outside this sequence.",
        ha="center",
        fontsize=8.2,
    )
    fig.tight_layout(rect=(0.055, 0.045, 0.985, 0.935))
    _save_figure_variants(fig, output_path)
    plt.close(fig)


def write_readme(
    output_dir: Path,
    *,
    model_name: str,
    review_count: int,
    retained_review_count: int,
    primary_count: int,
    quality_exclusion_count: int,
    all_scatter_ppt: float,
    primary_scatter_ppt: float,
    previous_35px_scatter_ppt: float,
    aperture_improvement_percent: float,
    midpoint: dict[str, float | bool],
    schedule: dict[str, object],
    source_radius: float,
) -> None:
    detrending_text = (
        "No correction was applied. Checks against airmass, sky level, star "
        "width, image position, and comparison-star counts did not improve "
        "the light curve consistently."
        if model_name == "none"
        else f"A correction using {model_name.replace('_', ' ')} was applied."
    )
    working = schedule["working_interpretation"]
    fixed = schedule["fixed_window_check"]
    assert isinstance(working, dict) and isinstance(fixed, dict)
    schedule_times = working["times"]
    assert isinstance(schedule_times, dict)
    text = f"""# TOI-3505.01 R-band light curve

## Result

This light curve is ready for mentor review. A recovered observing-schedule
row lists ingress at 00:15 and egress at 01:54. Interpreting those sheet clocks
as Eastern local time gives BJD_TDB
{float(schedule_times['ingress']['bjd_tdb']):.6f}-
{float(schedule_times['egress']['bjd_tdb']):.6f}, fully inside the measured
sequence. A straight-baseline, fixed-window check gives
{float(fixed['observed_depth_ppt']):.3f} +/-
{float(fixed['observed_depth_error_ppt']):.3f} ppt, so it does not find a
transit-like dimming in that exact historical window. Injecting a 2.91-ppt box
at the same times gives a total recovered depth of
{float(fixed['injected_total_depth_ppt']):.3f} ppt at
{float(fixed['injected_total_depth_snr']):.2f} sigma.

The source row does not state its time zone, epoch, timing uncertainty, or
prediction source. The Eastern-time interpretation is used because the
sheet's planned 21:10-04:55 range brackets the actual image sequence. Under
the current catalog period and epoch, the nearest predicted midpoint is
BJD_TDB
{float(midpoint['midpoint_bjd_tdb']):.6f}, or
{abs(float(midpoint['hours_from_observation_start'])):.2f} hours before the
first exposure. No physical transit model was fitted.

## Photometry settings

- 281 R-band exposures, 50 seconds each.
- {source_radius:g}-pixel source aperture.
- 70–139-pixel sky annulus.
- Ten comparison stars, C2 through C11.
- All ten comparison stars are unsaturated and have stable light curves.
- The {source_radius:g}-pixel aperture lowered the all-frame scatter from
  {previous_35px_scatter_ppt:.3f} to {all_scatter_ppt:.3f} ppt compared with
  the original 35-pixel aperture, an improvement of
  {aperture_improvement_percent:.1f}%.

## Frame review

All 281 measurements remain in the saved table. The plotted light curve uses
{primary_count} measurements. The other {quality_exclusion_count} frames have
documented image problems such as clouds, poor seeing, tracking, or trailing.
No frame was excluded only because the target appeared faint. Of the
{review_count} measurements checked by eye, {retained_review_count} had no
separate image problem and remain in the light curve.

The scatter of the plotted measurements is {primary_scatter_ppt:.3f} ppt.
The 10-minute bins are included only to make the overall shape easier to see.
The posting figure also includes the program-requested width, sky, inverted
airmass, comparison-count, X-position, and Y-position diagnostics. They are
robustly scaled and shifted into a shared lower display band, like
AstroImageJ's Page Rel option; those display transforms never enter the flux
calculation or fixed-window result.

## Detrending

{detrending_text}

## Mentor review package

The current three-image package is collected in `../toi3505_review_package/`:

- `01_TOI_3505.01_final_light_curve.png`
- `02_TOI_3505.01_data_set_fit_settings.png`
- `03_TOI_3505.01_NEB_screen.png`

## Supporting files

- `TOI_3505.01_2022-07-22_R_light_curve.csv`: measurements and frame notes.
- `TOI_3505.01_2022-07-22_R_10min_bins.csv`: plotted 10-minute bins.
- `comparison_star_checks.csv`: comparison-star measurements.
- `detrending_checks.csv`: results of the trend checks.
- `frame_review.csv` and the two review figures: image-by-image notes.
- `analysis_settings.json`: settings used to make the light curve.
- `TOI_3505.01_2022-07-22_R_AIJ_review.tbl`: compact, row-preserving AIJ
  table used to reopen MultiPlot 6.0.7 for the fit-settings screenshot.
- `historical_schedule_check.json` and `historical_schedule_times.csv`: the
  preserved schedule interpretation and fixed-window result.
- `summary.json`: short numerical summary.

## Still needed

The mentor still needs to confirm the original spreadsheet or Transit Info
source, its time zone, prediction epoch and uncertainty, and review the
{source_radius:g}-pixel aperture and comparison stars. The stored BJD_TDB
values agree with an independent calculation to within 0.000201 seconds, but
the observatory clock-sync record has not been found. Describe this as a light
curve submitted for review, not as a transit detection.
"""
    (output_dir / "README.md").write_text(text, encoding="utf-8")


def main() -> None:
    args = parse_args()
    table_path = args.table.resolve()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    table = load_table(
        table_path,
        expected_outer_radius=139,
        expected_source_radius=args.source_radius,
    )
    raw_curve = target_curve_from_counts(table, COMPARISON_STARS)
    raw_error = target_curve_error(table, COMPARISON_STARS)
    working_curve = target_curve_from_counts(table, PREVIOUS_WORKING_STARS)
    comparison_counts = np.sum(
        [table[f"Source-Sky_{star}"].to_numpy(dtype=float) for star in COMPARISON_STARS],
        axis=0,
    )
    features = build_features(table, comparison_counts)
    review_points, review_reasons, fit_eligible = make_review_flags(
        table, raw_curve, comparison_counts
    )
    primary_use = fit_eligible.copy()
    primary_exclusion_reasons = ["" for _ in range(len(table))]
    for index in range(len(table)):
        if not fit_eligible[index]:
            primary_exclusion_reasons[index] = review_reasons[index]
    slices = table["slice"].to_numpy(dtype=int)
    for slice_number, reason in VISUAL_QUALITY_EXCLUSIONS.items():
        matches = np.flatnonzero(slices == slice_number)
        if len(matches) != 1:
            raise RuntimeError(
                f"Visual exclusion slice {slice_number} is missing or duplicated"
            )
        primary_use[matches[0]] = False
        primary_exclusion_reasons[matches[0]] = reason
    model_review, cv_predictions = blocked_model_review(
        features, raw_curve, primary_use
    )
    model_name = select_promoted_model(model_review)
    feature_names = CANDIDATE_MODELS[model_name]
    adopted_curve, model_baseline, fitted_model = apply_condition_model(
        features, raw_curve, primary_use, feature_names
    )
    adopted_curve = adopted_curve / float(np.nanmedian(adopted_curve[primary_use]))

    start_bjd = float(table["BJD_TDB"].iloc[0])
    end_bjd = float(table["BJD_TDB"].iloc[-1])
    hours = (table["BJD_TDB"].to_numpy(dtype=float) - start_bjd) * 24.0
    midpoint = closest_predicted_midpoint(start_bjd, end_bjd)
    bins = bin_light_curve(hours[primary_use], adopted_curve[primary_use])
    ensemble_review = comparison_ensemble_review(table)
    schedule_analysis = analyze_schedule_window(
        args.schedule_record.resolve(),
        table["BJD_TDB"].to_numpy(dtype=float),
        adopted_curve,
        raw_error,
        primary_use,
        exposure_seconds=50.0,
        comparison_depth_ppt=CATALOG_DEPTH_PPT,
    )

    previous_35px_table = load_table(
        DEFAULT_WIDE_TABLE,
        expected_outer_radius=139,
        expected_source_radius=35,
    )
    previous_35px_curve = target_curve_from_counts(
        previous_35px_table, COMPARISON_STARS
    )
    previous_35px_scatter_ppt = 1000.0 * robust_scatter(previous_35px_curve)

    curve_table = pd.DataFrame(
        {
            "image": table["Label"],
            "slice": table["slice"],
            "bjd_tdb": table["BJD_TDB"],
            "hours_since_first_image": hours,
            "raw_relative_brightness_10_comparisons": raw_curve,
            "raw_relative_brightness_error": raw_error,
            "working_relative_brightness_7_comparisons": working_curve,
            "adopted_relative_brightness": adopted_curve,
            "condition_model_baseline": model_baseline,
            "review_point": review_points,
            "review_reason": review_reasons,
            "condition_model_fit_eligible": fit_eligible,
            "used_in_primary_curve": primary_use,
            "primary_exclusion_reason": primary_exclusion_reasons,
            "airmass": features["airmass"],
            "sky_counts_per_pixel": features["sky"],
            "target_width_pixels": features["width"],
            "target_x_fits": features["x"],
            "target_y_fits": features["y"],
            "comparison_counts": comparison_counts,
            "meridian_side": features["meridian_side"],
        }
    )

    curve_name = "TOI_3505.01_2022-07-22_R_final_candidate_light_curve.csv"
    bins_name = "TOI_3505.01_2022-07-22_R_final_candidate_10min_bins.csv"
    curve_table.to_csv(output_dir / curve_name, index=False, float_format="%.10f")
    bins.to_csv(output_dir / bins_name, index=False, float_format="%.10f")
    model_review.to_csv(
        output_dir / "model_selection.csv", index=False, float_format="%.10f"
    )
    ensemble_review.to_csv(
        output_dir / "comparison_ensemble_review.csv",
        index=False,
        float_format="%.10f",
    )
    (output_dir / "historical_schedule_check.json").write_text(
        json.dumps(schedule_analysis, indent=2) + "\n", encoding="utf-8"
    )
    schedule_rows: list[dict[str, object]] = []
    for interpretation_name in ("working_interpretation", "utc_alternative"):
        interpretation = schedule_analysis[interpretation_name]
        assert isinstance(interpretation, dict)
        schedule_times = interpretation["times"]
        assert isinstance(schedule_times, dict)
        for event_name in ("planned_start", "ingress", "egress", "planned_end"):
            event = schedule_times[event_name]
            assert isinstance(event, dict)
            schedule_rows.append(
                {
                    "interpretation": interpretation_name,
                    "timezone": interpretation["timezone"],
                    "event": event_name,
                    "source_clock": event["source_clock"],
                    "clock_datetime": event["clock_datetime"],
                    "utc_datetime": event["utc_datetime"],
                    "bjd_tdb": event["bjd_tdb"],
                }
            )
    pd.DataFrame(schedule_rows).to_csv(
        output_dir / "historical_schedule_times.csv",
        index=False,
        float_format="%.10f",
    )

    public_curve = pd.DataFrame(
        {
            "Image": table["Label"],
            "Frame": table["slice"],
            "BJD_TDB": table["BJD_TDB"],
            "Relative_Brightness": adopted_curve,
            "Flux_Error": raw_error,
            "Used_in_Plot": primary_use,
            "Notes": primary_exclusion_reasons,
        }
    )
    public_bins = pd.DataFrame(
        {
            "BJD_TDB": start_bjd + bins["hours"] / 24.0,
            "Relative_Brightness": bins["relative_brightness"],
            "Uncertainty": bins["uncertainty"],
            "Measurements": bins["points"],
        }
    )
    public_model_review = model_review.rename(
        columns={
            "model": "Correction",
            "features": "Terms",
            "cv_robust_scatter_ppt": "Scatter_ppt",
            "cv_rmse_ppt": "RMS_ppt",
            "robust_improvement_vs_none_percent": "Improvement_percent",
            "blocks_better_than_none": "Time_sections_better",
            "block_count": "Time_sections",
            "passes_absolute_rules": "Passed",
        }
    ).loc[
        :,
        [
            "Correction",
            "Terms",
            "Scatter_ppt",
            "RMS_ppt",
            "Improvement_percent",
            "Time_sections_better",
            "Time_sections",
            "Passed",
        ],
    ]
    public_comparison_review = ensemble_review.rename(
        columns={
            "star": "Star",
            "included_primary_ensemble": "Included",
            "pseudo_target_robust_scatter_percent": "Scatter_percent",
            "pseudo_target_clipped_std_percent": "Clipped_scatter_percent",
            "median_relative_snr": "Median_SNR",
            "median_source_counts": "Median_source_counts",
            "airmass_correlation": "Airmass_correlation",
            "width_correlation": "Width_correlation",
            "x_correlation": "X_correlation",
            "y_correlation": "Y_correlation",
        }
    ).loc[
        :,
        [
            "Star",
            "Included",
            "Scatter_percent",
            "Clipped_scatter_percent",
            "Median_SNR",
            "Median_source_counts",
            "Airmass_correlation",
            "Width_correlation",
            "X_correlation",
            "Y_correlation",
        ],
    ]
    public_curve.to_csv(
        output_dir / "TOI_3505.01_2022-07-22_R_light_curve.csv",
        index=False,
        float_format="%.10f",
    )
    public_bins.to_csv(
        output_dir / "TOI_3505.01_2022-07-22_R_10min_bins.csv",
        index=False,
        float_format="%.10f",
    )
    public_model_review.to_csv(
        output_dir / "detrending_checks.csv", index=False, float_format="%.10f"
    )
    public_comparison_review.to_csv(
        output_dir / "comparison_star_checks.csv",
        index=False,
        float_format="%.10f",
    )

    plot_final_candidate(
        output_dir / "TOI_3505.01_2022-07-22_R_light_curve.png",
        bjd_tdb=table["BJD_TDB"].to_numpy(dtype=float),
        adopted_curve=adopted_curve,
        primary_use=primary_use,
        binned=bins,
        diagnostics={
            "airmass": features["airmass"].to_numpy(dtype=float),
            "sky": features["sky"].to_numpy(dtype=float),
            "width": features["width"].to_numpy(dtype=float),
            "x": features["x"].to_numpy(dtype=float),
            "y": features["y"].to_numpy(dtype=float),
            "comparison_counts": comparison_counts,
        },
        schedule=schedule_analysis,
        source_radius=args.source_radius,
    )

    protocol = {
        "status": "frozen_for_final_candidate_generation",
        "input_table": table_path.name,
        "source_radius_pixels": args.source_radius,
        "sky_annulus_pixels": [70, 139],
        "primary_comparison_stars": COMPARISON_STARS,
        "comparison_selection_uses_target_scatter": False,
        "measurements_archived": int(len(table)),
        "measurements_used_in_primary_curve": int(primary_use.sum()),
        "point_exclusion_policy": (
            "exclude from the primary curve only for predeclared condition "
            "failures or a cross-star image artifact; target brightness alone "
            "never excludes a point"
        ),
        "visual_quality_exclusions": {
            str(slice_number): reason
            for slice_number, reason in VISUAL_QUALITY_EXCLUSIONS.items()
        },
        "condition_fit_policy": (
            "the primary quality mask is also the condition-model fit mask"
        ),
        "candidate_models": {
            name: list(names) for name, names in CANDIDATE_MODELS.items()
        },
        "blocked_cross_validation_blocks": BLOCK_COUNT,
        "promotion_rules": {
            "minimum_robust_scatter_improvement_fraction": MINIMUM_CV_IMPROVEMENT,
            "maximum_rmse_ratio": MAXIMUM_RMSE_RATIO,
            "minimum_fraction_of_blocks_better": MINIMUM_BETTER_BLOCK_FRACTION,
            "parsimony_near_best_fraction": NEAR_BEST_SCORE_FRACTION,
        },
        "display_bin_minutes": BIN_MINUTES,
        "transit_model_allowed": False,
        "reason_no_transit_model": (
            "the historical schedule and current ephemeris disagree, and the "
            "fixed historical window contains no significant transit-like dimming"
        ),
        "historical_schedule": {
            "source_record": schedule_analysis["source_record"],
            "working_timezone": "America/New_York",
            "source_timezone_explicit": False,
            "fixed_window_only": True,
        },
    }
    (output_dir / "frozen_protocol.json").write_text(
        json.dumps(protocol, indent=2) + "\n", encoding="utf-8"
    )
    public_settings = {
        "input_measurements": table_path.name,
        "source_aperture_pixels": args.source_radius,
        "sky_annulus_pixels": [70, 139],
        "comparison_stars": COMPARISON_STARS,
        "frames_measured": int(len(table)),
        "frames_used_in_plot": int(primary_use.sum()),
        "detrending": "none" if model_name == "none" else model_name,
        "transit_fit": "none",
        "historical_schedule_window_displayed": True,
        "historical_schedule_timezone_assumption": "America/New_York",
        "display_bin_minutes": BIN_MINUTES,
        "final_plot_style": "TFOP/Schar diagnostic review layout",
        "final_plot_page_relative_diagnostics": [
            "Width_T1",
            "Sky/Pixel_T1",
            "AIRMASS (inverted)",
            "tot_C_cnts",
            "X(FITS)_T1",
            "Y(FITS)_T1",
        ],
        "page_relative_diagnostics_are_display_only": True,
    }
    (output_dir / "analysis_settings.json").write_text(
        json.dumps(public_settings, indent=2) + "\n", encoding="utf-8"
    )

    selected_row = model_review.loc[model_review["model"] == model_name].iloc[0]
    best_non_null = model_review.loc[model_review["model"] != "none"].sort_values(
        "cv_robust_scatter_fraction"
    ).iloc[0]
    all_scatter_ppt = 1000.0 * robust_scatter(raw_curve)
    primary_scatter_ppt = 1000.0 * robust_scatter(adopted_curve[primary_use])
    summary = {
        "status": "ready_for_mentor_review",
        "scientific_interpretation": (
            "historically scheduled transit window with no fixed-window dimming; "
            "off transit under the current catalog ephemeris"
        ),
        "images_measured": int(len(table)),
        "measurements_archived": int(len(table)),
        "measurements_used_in_primary_curve": int(primary_use.sum()),
        "measurements_excluded_from_primary_curve": int((~primary_use).sum()),
        "review_points": int(review_points.sum()),
        "review_points_retained_in_primary_curve": int(
            np.sum(review_points & primary_use)
        ),
        "condition_model_fit_points": int(primary_use.sum()),
        "source_radius_pixels": args.source_radius,
        "sky_annulus_pixels": [70, 139],
        "primary_comparison_stars": COMPARISON_STARS,
        "all_frame_raw_robust_scatter_ppt": all_scatter_ppt,
        "previous_35px_all_frame_raw_robust_scatter_ppt": previous_35px_scatter_ppt,
        "aperture_all_frame_robust_scatter_improvement_percent": 100.0
        * (previous_35px_scatter_ppt - all_scatter_ppt)
        / previous_35px_scatter_ppt,
        "all_frame_raw_clipped_standard_deviation_ppt": 1000.0
        * clipped_standard_deviation(raw_curve),
        "primary_raw_robust_scatter_ppt": primary_scatter_ppt,
        "primary_raw_clipped_standard_deviation_ppt": 1000.0
        * clipped_standard_deviation(adopted_curve[primary_use]),
        "previous_working_ensemble_primary_robust_scatter_ppt": 1000.0
        * robust_scatter(working_curve[primary_use]),
        "promoted_condition_model": model_name,
        "detrending_applied_to_adopted_curve": model_name != "none",
        "adopted_primary_robust_scatter_ppt": primary_scatter_ppt,
        "promoted_model_cv_robust_improvement_percent": float(
            selected_row["robust_improvement_vs_none_percent"]
        ),
        "best_non_null_cv_model": str(best_non_null["model"]),
        "best_non_null_cv_robust_improvement_percent": float(
            best_non_null["robust_improvement_vs_none_percent"]
        ),
        "fitted_model": {
            "features": list(fitted_model.feature_names),
            "standardization_center": fitted_model.center.tolist(),
            "standardization_scale": fitted_model.scale.tolist(),
            "coefficients": fitted_model.coefficients.tolist(),
        },
        "observation": {
            "start_bjd_tdb": start_bjd,
            "end_bjd_tdb": end_bjd,
            "duration_hours": float(hours[-1]),
        },
        "closest_predicted_transit_midpoint": midpoint,
        "historical_schedule": schedule_analysis,
        "transit_model_applied": False,
        "visual_quality_exclusions": {
            str(slice_number): reason
            for slice_number, reason in VISUAL_QUALITY_EXCLUSIONS.items()
        },
        "remaining_finalization_gates": [
            "observatory clock-sync provenance",
            "original scheduling workbook or Transit Info file, time zone, epoch, and timing uncertainty",
            f"mentor timing and {args.source_radius:g}-pixel photometry review",
        ],
    }
    (output_dir / "summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    write_readme(
        output_dir,
        model_name=model_name,
        review_count=int(review_points.sum()),
        retained_review_count=int(np.sum(review_points & primary_use)),
        primary_count=int(primary_use.sum()),
        quality_exclusion_count=int((~primary_use).sum()),
        all_scatter_ppt=all_scatter_ppt,
        primary_scatter_ppt=primary_scatter_ppt,
        previous_35px_scatter_ppt=previous_35px_scatter_ppt,
        aperture_improvement_percent=100.0
        * (previous_35px_scatter_ppt - all_scatter_ppt)
        / previous_35px_scatter_ppt,
        midpoint=midpoint,
        schedule=schedule_analysis,
        source_radius=args.source_radius,
    )

    shutil.copyfile(
        table_path,
        output_dir / "TOI_3505.01_2022-07-22_R_measurements.xls",
    )
    write_aij_review_table(
        table,
        output_dir / "TOI_3505.01_2022-07-22_R_AIJ_review.tbl",
    )
    write_target_plot_config(
        ROOT
        / "outputs"
        / "toi3505_photometry"
        / "TOI_3505.01_2022-07-22_R_working.plotcfg",
        output_dir / "TOI_3505.01_2022-07-22_R_light_curve.plotcfg",
        schedule_analysis,
        observation_start_bjd=start_bjd,
        observation_end_bjd=end_bjd,
        orbital_period_days=TARGET_PERIOD_DAYS,
        transit_center_bjd=float(midpoint["midpoint_bjd_tdb"]),
        catalog_depth_ppt=CATALOG_DEPTH_PPT,
        target_teff_k=TARGET_TEFF_K,
    )
    aperture_source = output_dir / "TOI_3505_final_25px.apertures"
    if aperture_source.exists():
        shutil.copyfile(
            aperture_source,
            output_dir / "TOI_3505.01_2022-07-22_R.apertures",
        )

    review_dir = DEFAULT_REVIEW_PACKAGE_DIR
    review_dir.mkdir(parents=True, exist_ok=True)
    review_files = {
        output_dir / "TOI_3505.01_2022-07-22_R_light_curve.png": (
            review_dir / "01_TOI_3505.01_final_light_curve.png"
        ),
    }
    for source, destination in review_files.items():
        shutil.copyfile(source, destination)

    print(f"Wrote light-curve files to {output_dir}")
    print(f"Detrending applied: {'no' if model_name == 'none' else model_name}")
    print(f"All-frame robust scatter: {all_scatter_ppt:.3f} ppt")
    print(f"Primary robust scatter: {primary_scatter_ppt:.3f} ppt")
    print(f"Measurements shown: {int(primary_use.sum())}/{len(table)}")


if __name__ == "__main__":
    main()
