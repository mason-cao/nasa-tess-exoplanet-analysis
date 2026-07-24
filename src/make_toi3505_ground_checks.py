"""Run the remaining ground-based checks for TOI-3505.01.

The current ephemeris does not put a transit in the sequence, but a recovered
2022 scheduling row does.  This script therefore keeps the current-ephemeris
and historical-window interpretations separate.  It compares three reasonable
comparison-star ensembles, checks binned scatter, injects simple box dips,
attaches Gaia metadata, and screens nearby TIC stars during the fixed
historical window without claiming a physical transit fit.

Every check uses simple aperture photometry or linear least squares so the
calculation can be followed directly from the saved tables.
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

import astropy.units as u
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS

from analyze_toi3505_photometry import (
    COMPARISON_STARS,
    TARGET_EPOCH_BJD_TDB,
    TARGET_PERIOD_DAYS,
    load_table,
    normalized,
    target_curve_from_counts,
)
from check_toi3505_aperture_radii import measure_frame
from make_toi3505_final_candidate import robust_sigma
from toi3505_schedule import (
    DEFAULT_SCHEDULE_RECORD,
    fixed_window_fraction,
    schedule_context,
)


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_TABLE = (
    ROOT
    / "outputs"
    / "toi3505_aperture_check"
    / "TOI_3505.01_2022-07-22_R_measurements_25px_70-139_AIJ.tbl"
)
DEFAULT_FINAL_CURVE = (
    ROOT
    / "outputs"
    / "toi3505_final_candidate"
    / "TOI_3505.01_2022-07-22_R_final_candidate_light_curve.csv"
)
DEFAULT_WCS_IMAGE = (
    ROOT
    / "data"
    / "ground"
    / "toi3505"
    / "plate_solved"
    / "TOI_3505.01_50.000s_R-0001_wcs.fits"
)
DEFAULT_ALIGNED_DIR = ROOT / "data" / "ground" / "toi3505" / "aligned"
DEFAULT_CATALOG_DIR = ROOT / "data" / "catalogs" / "toi3505"
DEFAULT_OUTPUT_DIR = ROOT / "outputs" / "toi3505_ground_checks"

STAR_NAMES = ("T1", *COMPARISON_STARS)
PIXEL_SCALE_ARCSEC = 0.36212365
GROUND_SOURCE_RADIUS_PIXELS = 25.0
GROUND_SOURCE_RADIUS_ARCSEC = GROUND_SOURCE_RADIUS_PIXELS * PIXEL_SCALE_ARCSEC
INJECTION_DURATION_HOURS = 2.0
INJECTION_DEPTHS_PPT = (1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0)
INJECTION_MIDPOINT_COUNT = 13
BIN_MINUTES = (2.0, 5.0, 10.0, 15.0, 20.0, 30.0, 45.0, 60.0)
NEARBY_RADIUS_ARCSEC = 60.0
NEARBY_SOURCE_RADIUS_PIXELS = 15.0
NEARBY_SKY_INNER_PIXELS = 30.0
NEARBY_SKY_OUTER_PIXELS = 50.0
CATALOG_DEPTH = 0.002910


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--table", type=Path, default=DEFAULT_TABLE)
    parser.add_argument("--final-curve", type=Path, default=DEFAULT_FINAL_CURVE)
    parser.add_argument("--wcs-image", type=Path, default=DEFAULT_WCS_IMAGE)
    parser.add_argument("--aligned-dir", type=Path, default=DEFAULT_ALIGNED_DIR)
    parser.add_argument("--catalog-dir", type=Path, default=DEFAULT_CATALOG_DIR)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument(
        "--schedule-record", type=Path, default=DEFAULT_SCHEDULE_RECORD
    )
    parser.add_argument(
        "--skip-nearby-images",
        action="store_true",
        help="Skip the slow image-level nearby-star measurement.",
    )
    return parser.parse_args()


def curve_from_ensemble(target: np.ndarray, ensemble: np.ndarray) -> np.ndarray:
    """Return a median-normalized target/ensemble light curve."""
    target = np.asarray(target, dtype=float)
    ensemble = np.asarray(ensemble, dtype=float)
    curve = np.full(len(target), np.nan, dtype=float)
    valid = (
        np.isfinite(target)
        & np.isfinite(ensemble)
        & (target > 0.0)
        & (ensemble > 0.0)
    )
    curve[valid] = target[valid] / ensemble[valid]
    return normalized(curve)


def build_ensemble_curves(table: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Compare flux-sum, equal-star, and inverse-error ensembles."""
    target = table["Source-Sky_T1"].to_numpy(dtype=float)
    comparison_flux = np.column_stack(
        [table[f"Source-Sky_{star}"].to_numpy(dtype=float) for star in COMPARISON_STARS]
    )
    comparison_error = np.column_stack(
        [table[f"Source_Error_{star}"].to_numpy(dtype=float) for star in COMPARISON_STARS]
    )
    medians = np.nanmedian(comparison_flux, axis=0)
    normalized_flux = comparison_flux / medians
    normalized_error = comparison_error / medians

    flux_sum_curve = target_curve_from_counts(table, COMPARISON_STARS)
    equal_ensemble = np.nanmean(normalized_flux, axis=1)
    equal_curve = curve_from_ensemble(target / np.nanmedian(target), equal_ensemble)

    weights = np.zeros_like(normalized_error)
    valid_error = np.isfinite(normalized_error) & (normalized_error > 0)
    weights[valid_error] = 1.0 / np.square(normalized_error[valid_error])
    weighted_ensemble = np.sum(weights * normalized_flux, axis=1) / np.sum(
        weights, axis=1
    )
    inverse_error_curve = curve_from_ensemble(
        target / np.nanmedian(target), weighted_ensemble
    )

    curves = pd.DataFrame(
        {
            "slice": table["slice"].to_numpy(dtype=int),
            "bjd_tdb": table["BJD_TDB"].to_numpy(dtype=float),
            "flux_sum_ensemble": flux_sum_curve,
            "equal_star_ensemble": equal_curve,
            "inverse_error_ensemble": inverse_error_curve,
        }
    )
    rows = []
    for name, description in (
        (
            "flux_sum_ensemble",
            "target divided by the sum of all ten comparison-star counts",
        ),
        (
            "equal_star_ensemble",
            "target divided by the mean of ten individually normalized stars",
        ),
        (
            "inverse_error_ensemble",
            "target divided by an exposure-by-exposure inverse-error weighted mean",
        ),
    ):
        values = curves[name].to_numpy(dtype=float)
        rows.append(
            {
                "ensemble": name,
                "description": description,
                "all_frame_robust_scatter_ppt": robust_sigma(values) * 1000.0,
            }
        )
    return curves, pd.DataFrame(rows)


def add_primary_scatter(
    metrics: pd.DataFrame, curves: pd.DataFrame, primary_use: np.ndarray
) -> pd.DataFrame:
    """Add same-mask statistics without using them to change the adopted curve."""
    updated = metrics.copy()
    primary_scatter = []
    standard_deviation = []
    for name in updated["ensemble"]:
        values = curves[str(name)].to_numpy(dtype=float)[primary_use]
        primary_scatter.append(robust_sigma(values) * 1000.0)
        standard_deviation.append(float(np.nanstd(values, ddof=1)) * 1000.0)
    updated["primary_robust_scatter_ppt"] = primary_scatter
    updated["primary_standard_deviation_ppt"] = standard_deviation
    updated["adopted_for_final_light_curve"] = (
        updated["ensemble"] == "flux_sum_ensemble"
    )
    updated["selection_note"] = (
        "The adopted all-ten flux sum remains fixed; this table is a robustness comparison."
    )
    return updated


def noise_vs_bin_size(
    bjd_tdb: np.ndarray, values: np.ndarray, primary_use: np.ndarray
) -> pd.DataFrame:
    """Measure how robust scatter changes in fixed-width time bins."""
    time = np.asarray(bjd_tdb, dtype=float)[primary_use]
    flux = np.asarray(values, dtype=float)[primary_use]
    finite = np.isfinite(time) & np.isfinite(flux)
    time = time[finite]
    flux = flux[finite]
    point_scatter = robust_sigma(flux)
    rows: list[dict[str, float | int]] = []
    for minutes in BIN_MINUTES:
        width_days = minutes / 1440.0
        identifiers = np.floor((time - np.min(time)) / width_days).astype(int)
        centers: list[float] = []
        counts: list[int] = []
        for identifier in np.unique(identifiers):
            selected = flux[identifiers == identifier]
            if len(selected) < 2:
                continue
            centers.append(float(np.mean(selected)))
            counts.append(len(selected))
        if len(centers) < 3:
            continue
        typical_count = float(np.median(counts))
        measured = robust_sigma(np.asarray(centers))
        white_expectation = point_scatter / np.sqrt(typical_count)
        rows.append(
            {
                "bin_minutes": minutes,
                "bins": len(centers),
                "median_points_per_bin": typical_count,
                "measured_robust_scatter_ppt": measured * 1000.0,
                "white_noise_expectation_ppt": white_expectation * 1000.0,
                "beta": measured / white_expectation,
            }
        )
    return pd.DataFrame(rows)


def integrated_box_fraction(
    time_hours: np.ndarray,
    midpoint_hours: float,
    duration_hours: float,
    exposure_hours: float,
) -> np.ndarray:
    """Fraction of each exposure inside a box-shaped event."""
    time = np.asarray(time_hours, dtype=float)
    exposure_start = time - exposure_hours / 2.0
    exposure_end = time + exposure_hours / 2.0
    event_start = midpoint_hours - duration_hours / 2.0
    event_end = midpoint_hours + duration_hours / 2.0
    overlap = np.maximum(
        0.0, np.minimum(exposure_end, event_end) - np.maximum(exposure_start, event_start)
    )
    return overlap / exposure_hours


def fit_linear_box(
    time_hours: np.ndarray,
    flux: np.ndarray,
    flux_error: np.ndarray,
    box_fraction: np.ndarray,
    use: np.ndarray,
) -> dict[str, float]:
    """Fit a straight baseline plus one fixed box depth."""
    time = np.asarray(time_hours, dtype=float)
    values = np.asarray(flux, dtype=float)
    errors = np.asarray(flux_error, dtype=float)
    box = np.asarray(box_fraction, dtype=float)
    keep = (
        np.asarray(use, dtype=bool)
        & np.isfinite(time)
        & np.isfinite(values)
        & np.isfinite(errors)
        & (errors > 0)
        & np.isfinite(box)
    )
    if keep.sum() < 10 or np.sum(box[keep] > 0) < 3:
        raise ValueError("Too few points for the requested box fit")
    centered_time = time[keep] - float(np.mean(time[keep]))
    design = np.column_stack(
        [np.ones(keep.sum()), centered_time, -box[keep]]
    )
    weights = 1.0 / np.square(errors[keep])
    normal = design.T @ (weights[:, None] * design)
    right = design.T @ (weights * values[keep])
    coefficients = np.linalg.solve(normal, right)
    residual = values[keep] - design @ coefficients
    degrees_of_freedom = keep.sum() - design.shape[1]
    reduced_chi_square = float(
        np.sum(np.square(residual / errors[keep])) / degrees_of_freedom
    )
    covariance = np.linalg.inv(normal) * reduced_chi_square
    return {
        "depth": float(coefficients[2]),
        "depth_error": float(np.sqrt(covariance[2, 2])),
        "baseline": float(coefficients[0]),
        "slope_per_hour": float(coefficients[1]),
        "residual_scatter": robust_sigma(residual),
        "reduced_chi_square": reduced_chi_square,
        "points": int(keep.sum()),
        "in_event_points": int(np.sum(box[keep] > 0)),
    }


def run_ground_injections(
    light_curve: pd.DataFrame, exposure_seconds: float = 50.0
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Inject two-hour dips at interior times and refit the same simple box."""
    time = light_curve["hours_since_first_image"].to_numpy(dtype=float)
    flux = light_curve["adopted_relative_brightness"].to_numpy(dtype=float)
    error = light_curve["raw_relative_brightness_error"].to_numpy(dtype=float)
    use = light_curve["used_in_primary_curve"].to_numpy(dtype=bool)
    margin = INJECTION_DURATION_HOURS / 2.0 + 0.25
    start = float(np.nanmin(time[use])) + margin
    end = float(np.nanmax(time[use])) - margin
    if end <= start:
        raise RuntimeError("The light curve is too short for the injection grid")
    midpoints = np.linspace(start, end, INJECTION_MIDPOINT_COUNT)
    rows: list[dict[str, float | int | bool | str]] = []
    for trial, midpoint in enumerate(midpoints, start=1):
        box = integrated_box_fraction(
            time,
            midpoint,
            INJECTION_DURATION_HOURS,
            exposure_seconds / 3600.0,
        )
        control = fit_linear_box(time, flux, error, box, use)
        for depth_ppt in INJECTION_DEPTHS_PPT:
            depth = depth_ppt / 1000.0
            injected = flux * (1.0 - depth * box)
            recovered = fit_linear_box(time, injected, error, box, use)
            increment = recovered["depth"] - control["depth"]
            rows.append(
                {
                    "trial": trial,
                    "midpoint_hours_since_first_image": midpoint,
                    "duration_hours": INJECTION_DURATION_HOURS,
                    "injected_depth_ppt": depth_ppt,
                    "control_depth_ppt": control["depth"] * 1000.0,
                    "recovered_total_depth_ppt": recovered["depth"] * 1000.0,
                    "recovered_increment_ppt": increment * 1000.0,
                    "total_bias_ppt": (recovered["depth"] - depth) * 1000.0,
                    "increment_bias_ppt": (increment - depth) * 1000.0,
                    "total_recovery_fraction": recovered["depth"] / depth,
                    "increment_recovery_fraction": increment / depth,
                    "formal_depth_error_ppt": recovered["depth_error"] * 1000.0,
                    "total_depth_snr": recovered["depth"] / recovered["depth_error"],
                    "recovered_above_3_sigma": bool(
                        recovered["depth"] >= 3.0 * recovered["depth_error"]
                    ),
                    "fit_residual_scatter_ppt": recovered["residual_scatter"]
                    * 1000.0,
                    "in_event_points": recovered["in_event_points"],
                    "scope": (
                        "two-hour box placed inside the observed sequence; the "
                        "trial grid is not tied to a verified event window"
                    ),
                }
            )
    table = pd.DataFrame(rows)
    summary = (
        table.groupby("injected_depth_ppt", as_index=False)
        .agg(
            trials=("trial", "size"),
            median_total_recovered_ppt=("recovered_total_depth_ppt", "median"),
            total_recovered_std_ppt=("recovered_total_depth_ppt", "std"),
            median_total_bias_ppt=("total_bias_ppt", "median"),
            median_increment_bias_ppt=("increment_bias_ppt", "median"),
            median_formal_error_ppt=("formal_depth_error_ppt", "median"),
            recovery_above_3_sigma_fraction=("recovered_above_3_sigma", "mean"),
        )
        .reset_index(drop=True)
    )
    return table, summary


def catalog_match_measured_stars(
    table: pd.DataFrame, wcs_image: Path, catalog_dir: Path
) -> pd.DataFrame:
    """Match T1 and C2-C11 to the targeted Gaia download when available."""
    wcs = WCS(fits.getheader(wcs_image))
    first = table.iloc[0]
    position_rows: list[dict[str, float | str]] = []
    for star in STAR_NAMES:
        x = float(first[f"X(FITS)_{star}"]) - 1.0
        y = float(first[f"Y(FITS)_{star}"]) - 1.0
        coordinate = wcs.pixel_to_world(x, y)
        position_rows.append(
            {
                "star": star,
                "x_zero_indexed": x,
                "y_zero_indexed": y,
                "wcs_ra_deg": float(coordinate.ra.deg),
                "wcs_dec_deg": float(coordinate.dec.deg),
            }
        )
    positions = pd.DataFrame(position_rows)
    catalog_path = catalog_dir / "ground_star_gaia_dr3_12arcsec.csv"
    if not catalog_path.exists():
        positions["catalog_match_status"] = (
            "targeted Gaia query not downloaded; run download_toi3505_comparison_gaia.py"
        )
        return positions

    catalog = pd.read_csv(catalog_path)
    rows: list[dict[str, object]] = []
    for position in positions.itertuples(index=False):
        neighborhood = catalog[catalog["measured_star"] == position.star].copy()
        if neighborhood.empty:
            row = position._asdict()
            row.update(
                {
                    "catalog_match_status": "no Gaia source within 12 arcseconds",
                    "gaia_source_id": np.nan,
                    "gaia_separation_arcsec": np.nan,
                    "gaia_g_mag": np.nan,
                    "gaia_bp_rp_mag": np.nan,
                    "gaia_ruwe": np.nan,
                    "gaia_variable_flag": "",
                    "gaia_non_single_star": np.nan,
                    "gaia_sources_within_ground_aperture": 0,
                    "additional_gaia_sources_within_ground_aperture": 0,
                }
            )
            rows.append(row)
            continue
        neighborhood = neighborhood.sort_values(
            "separation_from_measured_position_arcsec"
        )
        nearest = neighborhood.iloc[0]
        in_aperture = neighborhood[
            neighborhood["separation_from_measured_position_arcsec"]
            <= GROUND_SOURCE_RADIUS_ARCSEC
        ]
        row = position._asdict()
        row.update(
            {
                "catalog_match_status": (
                    "matched"
                    if float(nearest["separation_from_measured_position_arcsec"]) <= 3.0
                    else "nearest Gaia source is more than 3 arcseconds away"
                ),
                "gaia_source_id": str(int(nearest["source_id"])),
                "gaia_separation_arcsec": float(
                    nearest["separation_from_measured_position_arcsec"]
                ),
                "gaia_g_mag": float(nearest["phot_g_mean_mag"]),
                "gaia_bp_rp_mag": float(nearest["bp_rp"]),
                "gaia_ruwe": float(nearest["ruwe"]),
                "gaia_variable_flag": str(nearest["phot_variable_flag"]),
                "gaia_non_single_star": int(nearest["non_single_star"]),
                "gaia_sources_within_ground_aperture": int(len(in_aperture)),
                "additional_gaia_sources_within_ground_aperture": max(
                    0, int(len(in_aperture)) - 1
                ),
            }
        )
        rows.append(row)
    return pd.DataFrame(rows)


def prepare_nearby_candidates(catalog_dir: Path, wcs: WCS) -> pd.DataFrame:
    """Select TIC sources bright enough to mimic the catalog depth."""
    tic = pd.read_csv(catalog_dir / "tic_v8_2p5arcmin.csv")
    target_rows = tic[tic["ID"] == 390988385]
    if len(target_rows) != 1:
        raise RuntimeError("TIC 390988385 is missing or duplicated")
    target_tmag = float(target_rows.iloc[0]["Tmag"])
    candidates = tic[
        (tic["ID"] != 390988385)
        & np.isfinite(tic["Tmag"])
        & np.isfinite(tic["dstArcSec"])
        & (tic["dstArcSec"] <= NEARBY_RADIUS_ARCSEC)
    ].copy()
    candidates["delta_tmag"] = candidates["Tmag"] - target_tmag
    candidates["flux_ratio_to_target"] = np.power(
        10.0, -0.4 * candidates["delta_tmag"]
    )
    candidates["eclipse_fraction_needed_simple"] = (
        CATALOG_DEPTH / candidates["flux_ratio_to_target"]
    )
    candidates = candidates[candidates["eclipse_fraction_needed_simple"] <= 1.0]

    # A few TIC rows can describe the same Gaia source.  Keep the brightest
    # entry for each Gaia identifier, or the TIC identifier when Gaia is blank.
    gaia_key = candidates["GAIA"].astype("string")
    usable_gaia = gaia_key.notna() & (gaia_key != "0") & (gaia_key != "nan")
    candidates["deduplication_key"] = np.where(
        usable_gaia, "gaia_" + gaia_key.fillna(""), "tic_" + candidates["ID"].astype(str)
    )
    candidates = (
        candidates.sort_values("Tmag")
        .drop_duplicates("deduplication_key", keep="first")
        .sort_values("dstArcSec")
        .reset_index(drop=True)
    )
    coordinates = SkyCoord(
        candidates["ra"].to_numpy(dtype=float) * u.deg,
        candidates["dec"].to_numpy(dtype=float) * u.deg,
    )
    x, y = wcs.world_to_pixel(coordinates)
    candidates["x_zero_indexed"] = x
    candidates["y_zero_indexed"] = y
    candidates["catalog_screen"] = (
        "bright enough to mimic 2.91 ppt only under a total-throughput simple eclipse calculation"
    )
    return candidates


def measure_nearby_stars(
    candidates: pd.DataFrame,
    table: pd.DataFrame,
    final_curve: pd.DataFrame,
    aligned_dir: Path,
    schedule: dict[str, object],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Measure catalog candidates on every aligned frame."""
    positions = candidates[["x_zero_indexed", "y_zero_indexed"]].to_numpy(
        dtype=float
    )
    flux = np.full((len(table), len(candidates)), np.nan, dtype=float)
    background = np.full_like(flux, np.nan)
    for row_index, row in table.iterrows():
        image_path = aligned_dir / str(row["Label"])
        if not image_path.exists():
            raise FileNotFoundError(image_path)
        with fits.open(image_path, memmap=True) as hdul:
            measured, sky = measure_frame(
                hdul[0].data,
                positions,
                radii=(NEARBY_SOURCE_RADIUS_PIXELS,),
                annulus_inner=NEARBY_SKY_INNER_PIXELS,
                annulus_outer=NEARBY_SKY_OUTER_PIXELS,
            )
        flux[row_index] = measured[0]
        background[row_index] = sky
        if (row_index + 1) % 25 == 0 or row_index + 1 == len(table):
            print(
                f"Measured nearby stars in {row_index + 1}/{len(table)} frames",
                flush=True,
            )

    comparison_counts = table["tot_C_cnts"].to_numpy(dtype=float)
    primary = final_curve["used_in_primary_curve"].to_numpy(dtype=bool)
    bjd = table["BJD_TDB"].to_numpy(dtype=float)
    phase = ((bjd - TARGET_EPOCH_BJD_TDB) / TARGET_PERIOD_DAYS) % 1.0
    working = schedule["working_interpretation"]
    assert isinstance(working, dict)
    schedule_times = working["times"]
    assert isinstance(schedule_times, dict)
    ingress = float(schedule_times["ingress"]["bjd_tdb"])
    egress = float(schedule_times["egress"]["bjd_tdb"])
    schedule_fraction = fixed_window_fraction(bjd, ingress, egress, 50.0)
    curves: dict[str, np.ndarray | pd.Series] = {
        "slice": table["slice"],
        "bjd_tdb": bjd,
        "orbital_phase": phase,
        "historical_schedule_window_fraction": schedule_fraction,
        "used_in_primary_ground_curve": primary,
    }
    rows: list[dict[str, object]] = []
    for index, candidate in candidates.iterrows():
        differential = curve_from_ensemble(flux[:, index], comparison_counts)
        key = f"tic_{int(candidate['ID'])}_relative_brightness"
        curves[key] = differential
        values = differential[primary & np.isfinite(differential)]
        center = float(np.nanmedian(values)) if len(values) else float("nan")
        scatter = robust_sigma(values)
        deviation = np.abs(values - center)
        finite_fraction = float(np.mean(np.isfinite(differential)))
        valid = primary & np.isfinite(differential)
        inside = valid & (schedule_fraction > 0.0)
        outside = valid & (schedule_fraction == 0.0)
        if inside.sum() >= 10 and outside.sum() >= 20:
            inside_values = differential[inside]
            outside_values = differential[outside]
            inside_scatter = robust_sigma(inside_values)
            outside_scatter = robust_sigma(outside_values)
            schedule_depth = float(
                np.median(outside_values) - np.median(inside_values)
            )
            schedule_error = float(
                np.sqrt(
                    inside_scatter**2 / inside.sum()
                    + outside_scatter**2 / outside.sum()
                )
            )
        else:
            schedule_depth = float("nan")
            schedule_error = float("nan")
        required_eclipse = float(candidate["eclipse_fraction_needed_simple"])
        target_aperture_overlap = bool(
            float(candidate["dstArcSec"])
            < (GROUND_SOURCE_RADIUS_PIXELS + NEARBY_SOURCE_RADIUS_PIXELS)
            * PIXEL_SCALE_ARCSEC
        )
        sufficient_measurements = bool(
            finite_fraction >= 0.80
            and inside.sum() >= 60
            and outside.sum() >= 120
            and np.isfinite(schedule_error)
            and schedule_error > 0.0
        )
        three_sigma_upper = (
            schedule_depth + 3.0 * schedule_error
            if np.isfinite(schedule_depth) and np.isfinite(schedule_error)
            else float("nan")
        )
        cleared = bool(
            sufficient_measurements
            and not target_aperture_overlap
            and three_sigma_upper < required_eclipse
        )
        rows.append(
            {
                "tic_id": int(candidate["ID"]),
                "gaia_id": candidate["GAIA"],
                "separation_arcsec": float(candidate["dstArcSec"]),
                "tmag": float(candidate["Tmag"]),
                "delta_tmag": float(candidate["delta_tmag"]),
                "flux_ratio_to_target": float(candidate["flux_ratio_to_target"]),
                "eclipse_fraction_needed_simple": float(
                    candidate["eclipse_fraction_needed_simple"]
                ),
                "valid_measurement_fraction": finite_fraction,
                "night_robust_scatter_ppt": scatter * 1000.0,
                "night_standard_deviation_ppt": float(np.nanstd(values, ddof=1))
                * 1000.0,
                "maximum_absolute_deviation_ppt": float(np.nanmax(deviation))
                * 1000.0,
                "points_beyond_five_robust_sigma": int(
                    np.sum(deviation > 5.0 * scatter)
                )
                if np.isfinite(scatter) and scatter > 0
                else 0,
                "observed_phase_min": float(np.nanmin(phase[primary])),
                "observed_phase_max": float(np.nanmax(phase[primary])),
                "historical_schedule_timezone_assumption": "America/New_York",
                "historical_window_points": int(inside.sum()),
                "outside_window_points": int(outside.sum()),
                "historical_window_depth_ppt": schedule_depth * 1000.0,
                "historical_window_depth_error_ppt": schedule_error * 1000.0,
                "historical_window_three_sigma_upper_ppt": three_sigma_upper
                * 1000.0,
                "required_eclipse_depth_ppt_simple": required_eclipse * 1000.0,
                "target_aperture_overlap": target_aperture_overlap,
                "sufficient_schedule_window_measurements": sufficient_measurements,
                "predicted_transit_covered": True,
                "transit_relevant_clearance": cleared,
                "interpretation": (
                    "inconsistent with the required eclipse at the three-sigma screening level"
                    if cleared
                    else "not cleared because of aperture overlap, incomplete/noisy measurements, or an insufficient depth limit"
                ),
            }
        )
    return pd.DataFrame(rows), pd.DataFrame(curves)


def plot_ground_robustness(
    output_dir: Path,
    curves: pd.DataFrame,
    metrics: pd.DataFrame,
    noise: pd.DataFrame,
    injections: pd.DataFrame,
    injection_summary: pd.DataFrame,
    primary_use: np.ndarray,
) -> None:
    """Plot the ensemble, binning, and injection checks."""
    figure, axes = plt.subplots(2, 2, figsize=(12.2, 9.0))
    hours = (curves["bjd_tdb"] - float(curves["bjd_tdb"].iloc[0])) * 24.0
    ensemble_names = (
        ("flux_sum_ensemble", "Flux sum", "o"),
        ("equal_star_ensemble", "Equal stars", "s"),
        ("inverse_error_ensemble", "Inverse error", "^"),
    )
    for offset, (name, label, marker) in enumerate(ensemble_names):
        display = 1000.0 * (curves[name].to_numpy(dtype=float) - 1.0) + 9.0 * offset
        axes[0, 0].plot(
            hours[primary_use],
            display[primary_use],
            marker,
            markersize=2.8,
            linestyle="none",
            alpha=0.68,
            label=label,
        )
    axes[0, 0].set_xlabel("Hours since first exposure")
    axes[0, 0].set_ylabel("Brightness change + display offset (ppt)")
    axes[0, 0].set_title("Comparison-star ensemble check")
    axes[0, 0].legend(fontsize=8.5, ncol=3)

    short_labels = ["Flux sum", "Equal stars", "Inverse error"]
    axes[0, 1].bar(
        short_labels,
        metrics["primary_robust_scatter_ppt"],
        color=["#31688e", "#6a7f58", "#9b6b43"],
    )
    axes[0, 1].set_ylabel("Robust scatter (ppt)")
    axes[0, 1].set_title("Same-mask scatter")
    axes[0, 1].tick_params(axis="x", rotation=15)
    for position, value in enumerate(metrics["primary_robust_scatter_ppt"]):
        axes[0, 1].text(position, value, f"{value:.2f}", ha="center", va="bottom")

    axes[1, 0].plot(
        noise["bin_minutes"],
        noise["measured_robust_scatter_ppt"],
        "o-",
        label="Measured",
    )
    axes[1, 0].plot(
        noise["bin_minutes"],
        noise["white_noise_expectation_ppt"],
        "--",
        color="0.35",
        label="White-noise expectation",
    )
    axes[1, 0].set_xlabel("Bin width (minutes)")
    axes[1, 0].set_ylabel("Robust scatter of bin means (ppt)")
    axes[1, 0].set_title("Noise versus bin width")
    axes[1, 0].legend(fontsize=8.5)

    grouped = injections.groupby("injected_depth_ppt")
    medians = grouped["recovered_total_depth_ppt"].median()
    low = grouped["recovered_total_depth_ppt"].quantile(0.16)
    high = grouped["recovered_total_depth_ppt"].quantile(0.84)
    x = medians.index.to_numpy(dtype=float)
    axes[1, 1].errorbar(
        x,
        medians.to_numpy(),
        yerr=np.vstack([medians.to_numpy() - low.to_numpy(), high.to_numpy() - medians.to_numpy()]),
        fmt="o",
        capsize=3,
        label="13 trial midpoints",
    )
    limits = [0.0, max(INJECTION_DEPTHS_PPT) + 1.0]
    axes[1, 1].plot(limits, limits, "--", color="0.35", label="Exact recovery")
    axes[1, 1].set_xlim(limits)
    axes[1, 1].set_ylim(limits)
    axes[1, 1].set_xlabel("Injected depth (ppt)")
    axes[1, 1].set_ylabel("Recovered total depth (ppt)")
    axes[1, 1].set_title("Two-hour box recovery")
    axes[1, 1].legend(fontsize=8.5)
    threshold_rows = injection_summary[
        injection_summary["recovery_above_3_sigma_fraction"] >= 0.9
    ]
    threshold_text = (
        f"At least 90% above 3σ by {threshold_rows.iloc[0]['injected_depth_ppt']:.0f} ppt"
        if len(threshold_rows)
        else "No tested depth reached 90% above 3σ"
    )
    axes[1, 1].text(
        0.04,
        0.94,
        threshold_text,
        transform=axes[1, 1].transAxes,
        ha="left",
        va="top",
        fontsize=9,
    )

    for axis in axes.flat:
        axis.grid(alpha=0.16)
    figure.suptitle("TOI-3505.01 ground light-curve checks", fontsize=15.5)
    figure.tight_layout(rect=(0, 0, 1, 0.96))
    figure.savefig(output_dir / "01_ground_light_curve_checks.png", dpi=230)
    figure.savefig(output_dir / "01_ground_light_curve_checks.svg")
    plt.close(figure)


def plot_nearby_screen(
    output_dir: Path,
    image_path: Path,
    wcs_image: Path,
    candidates: pd.DataFrame,
    measurements: pd.DataFrame,
) -> None:
    """Show candidate locations and the conditional schedule-window screen."""
    with fits.open(image_path, memmap=True) as hdul:
        image = np.asarray(hdul[0].data, dtype=float)
    wcs = WCS(fits.getheader(wcs_image))
    target = SkyCoord(297.043476 * u.deg, 18.698914 * u.deg)
    target_x, target_y = wcs.world_to_pixel(target)
    radius_pixels = NEARBY_RADIUS_ARCSEC / PIXEL_SCALE_ARCSEC
    x_min = max(0, int(np.floor(target_x - radius_pixels - 20)))
    x_max = min(image.shape[1], int(np.ceil(target_x + radius_pixels + 20)))
    y_min = max(0, int(np.floor(target_y - radius_pixels - 20)))
    y_max = min(image.shape[0], int(np.ceil(target_y + radius_pixels + 20)))
    crop = image[y_min:y_max, x_min:x_max]
    low, high = np.nanpercentile(crop, [20, 99.7])

    figure, axes = plt.subplots(1, 2, figsize=(12.0, 5.7))
    axes[0].imshow(crop, origin="lower", cmap="gray", vmin=low, vmax=high)
    axes[0].scatter(
        candidates["x_zero_indexed"] - x_min,
        candidates["y_zero_indexed"] - y_min,
        facecolors="none",
        edgecolors="#d18435",
        s=45,
        linewidths=0.9,
        label="Catalog candidates",
    )
    axes[0].scatter(
        target_x - x_min,
        target_y - y_min,
        marker="+",
        s=120,
        linewidths=2,
        color="#2d6f9e",
        label="TOI-3505.01",
    )
    axes[0].set_xlabel("Image column (pixels)")
    axes[0].set_ylabel("Image row (pixels)")
    axes[0].set_title("Sources within 60 arcseconds")
    axes[0].legend(fontsize=8.5)

    cleared = measurements["transit_relevant_clearance"].astype(bool)
    axes[1].scatter(
        100.0 * measurements.loc[~cleared, "eclipse_fraction_needed_simple"],
        measurements.loc[~cleared, "night_robust_scatter_ppt"],
        s=34,
        alpha=0.75,
        color="#a45c4a",
        label="Not cleared",
    )
    axes[1].scatter(
        100.0 * measurements.loc[cleared, "eclipse_fraction_needed_simple"],
        measurements.loc[cleared, "night_robust_scatter_ppt"],
        s=34,
        alpha=0.75,
        color="#35745c",
        label="Cleared by this screen",
    )
    axes[1].set_xscale("log")
    axes[1].set_yscale("log")
    axes[1].set_xlabel("Eclipse needed to mimic 2.91 ppt (%)")
    axes[1].set_ylabel("Nightly robust scatter (ppt)")
    axes[1].set_title("Fixed schedule-window source screen")
    axes[1].legend(fontsize=8.5)
    axes[1].grid(alpha=0.16, which="both")
    figure.suptitle(
        "Nearby-star screen for the 2022 schedule window (EDT assumed)",
        fontsize=15,
    )
    figure.tight_layout(rect=(0, 0, 1, 0.94))
    figure.savefig(output_dir / "02_nearby_star_screen.png", dpi=230)
    plt.close(figure)


def measured_star_photometry(table: pd.DataFrame) -> pd.DataFrame:
    """Summarize brightness and scatter for the measured target and stars."""
    rows: list[dict[str, object]] = []
    for star in STAR_NAMES:
        relative = normalized(table[f"rel_flux_{star}"].to_numpy(dtype=float))
        rows.append(
            {
                "star": star,
                "role": "target" if star == "T1" else "comparison",
                "median_source_counts": float(
                    np.nanmedian(table[f"Source-Sky_{star}"])
                ),
                "median_source_snr": float(
                    np.nanmedian(table[f"Source_SNR_{star}"])
                ),
                "maximum_peak_counts": float(np.nanmax(table[f"Peak_{star}"])),
                "relative_flux_robust_scatter_ppt": robust_sigma(relative)
                * 1000.0,
                "frames_with_global_saturation_flag": int(
                    np.sum(table["Saturated"].to_numpy(dtype=float) > 0)
                ),
            }
        )
    return pd.DataFrame(rows)


def write_readme(
    output_dir: Path,
    ensemble_metrics: pd.DataFrame,
    noise: pd.DataFrame,
    injection_summary: pd.DataFrame,
    catalog_matches: pd.DataFrame,
    nearby_measurements: pd.DataFrame | None,
    schedule: dict[str, object],
) -> None:
    adopted = ensemble_metrics[
        ensemble_metrics["adopted_for_final_light_curve"]
    ].iloc[0]
    beta_30 = noise.iloc[(noise["bin_minutes"] - 30.0).abs().argsort()[:1]].iloc[0]
    threshold_rows = injection_summary[
        injection_summary["recovery_above_3_sigma_fraction"] >= 0.9
    ]
    threshold_text = (
        f"The first tested depth with at least 90% of placements above three sigma was {threshold_rows.iloc[0]['injected_depth_ppt']:.0f} ppt."
        if len(threshold_rows)
        else "None of the tested depths reached three sigma in at least 90% of placements."
    )
    matched = int(np.sum(catalog_matches["catalog_match_status"] == "matched"))
    variable_count = (
        int(
            np.sum(
                catalog_matches.get(
                    "gaia_variable_flag", pd.Series(dtype=str)
                ).astype(str)
                == "VARIABLE"
            )
        )
        if "gaia_variable_flag" in catalog_matches
        else 0
    )
    if nearby_measurements is None:
        nearby_text = (
            "The image-level nearby-star measurement was skipped. Run this script "
            "without `--skip-nearby-images` to complete it."
        )
    else:
        cleared = int(nearby_measurements["transit_relevant_clearance"].sum())
        overlap = int(nearby_measurements["target_aperture_overlap"].sum())
        nearby_text = (
            f"{len(nearby_measurements)} deduplicated TIC sources within 60 arcseconds "
            "were bright enough to mimic a 2.91-ppt event in the simple total-eclipse "
            "screen and were measured on all 281 aligned images. Under the documented "
            "Eastern-time interpretation of the 2022 schedule, "
            f"{cleared} sources are inconsistent with the required eclipse at this "
            f"screen's three-sigma level. {overlap} source apertures overlap the "
            "target aperture and are not cleared."
        )
    working = schedule["working_interpretation"]
    assert isinstance(working, dict)
    schedule_times = working["times"]
    assert isinstance(schedule_times, dict)
    text = f"""# TOI-3505.01 ground-data checks

## Comparison ensemble

The final light curve still uses the predeclared sum of all ten comparison
stars. Its robust scatter on the final frame mask is
{adopted['primary_robust_scatter_ppt']:.3f} ppt. Equal-star and inverse-error
ensembles are saved as robustness comparisons; neither was used to retune the
published curve after looking at the target.

At a 30-minute bin width, the measured robust scatter is
{beta_30['measured_robust_scatter_ppt']:.3f} ppt, compared with a
{beta_30['white_noise_expectation_ppt']:.3f}-ppt white-noise expectation
(beta = {beta_30['beta']:.2f}). This is a descriptive time-correlation check,
not a replacement for a transit fit.

## Injection check

Two-hour box dips from 1 to 10 ppt were placed at 13 interior times. Each was
fit with a straight baseline and a fixed box. The control depth at the same
time is retained, so the saved total recovery shows the real phase-dependent
structure rather than an artificially perfect control-subtracted result.
{threshold_text} These trials test this light curve and fitting method; they do
not simulate a target-only point-spread function or establish a corrected TESS
dilution.

## Star catalog check

{matched} of {len(catalog_matches)} measured positions have a Gaia source
within 3 arcseconds in the targeted query. Gaia marks {variable_count} matched
sources as `VARIABLE`. A Gaia value of `NOT_AVAILABLE` is not proof that a star
is constant; the ground pseudo-target curves remain the direct stability
check. The table also counts Gaia sources inside the 9.05-arcsecond ground
aperture so blends are visible instead of hidden.

## Nearby-star scope

{nearby_text}

The schedule-window result uses ingress BJD_TDB
{float(schedule_times['ingress']['bjd_tdb']):.6f} and egress
{float(schedule_times['egress']['bjd_tdb']):.6f}. The source row does not state
its time zone, epoch, uncertainty, or prediction source. Clearance here is a
conservative image-level screen, not the program's formal AstroImageJ NEB
procedure or planet validation. It cannot resolve or clear the known
0.517-arcsecond companion. The current ephemeris still places its nearest
event about 17.4 hours before this sequence.

## Files

- `01_ground_light_curve_checks.png`: ensemble, binning, and injection checks.
- `comparison_ensemble_light_curves.csv` and `comparison_ensemble_metrics.csv`.
- `noise_vs_bin_size.csv`.
- `ground_light_curve_injections.csv` and `ground_injection_summary.csv`.
- `comparison_star_catalog_matches.csv`.
- `nearby_star_catalog_candidates.csv`, plus image measurements when run.
- `summary.json`: key settings, counts, and limitations.
"""
    (output_dir / "README.md").write_text(text, encoding="utf-8")


def main() -> None:
    args = parse_args()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    table = load_table(
        args.table.resolve(), expected_outer_radius=139, expected_source_radius=25
    )
    final_curve = pd.read_csv(args.final_curve.resolve())
    if len(final_curve) != len(table):
        raise RuntimeError("The final light curve and AIJ table have different lengths")
    if not np.array_equal(
        final_curve["slice"].to_numpy(dtype=int), table["slice"].to_numpy(dtype=int)
    ):
        raise RuntimeError("The final light curve and AIJ table slices do not align")
    primary_use = final_curve["used_in_primary_curve"].to_numpy(dtype=bool)
    schedule = schedule_context(
        args.schedule_record.resolve(),
        observation_start_bjd=float(table["BJD_TDB"].min()),
        observation_end_bjd=float(table["BJD_TDB"].max()),
    )

    ensemble_curves, ensemble_metrics = build_ensemble_curves(table)
    ensemble_metrics = add_primary_scatter(
        ensemble_metrics, ensemble_curves, primary_use
    )
    noise = noise_vs_bin_size(
        final_curve["bjd_tdb"].to_numpy(dtype=float),
        final_curve["adopted_relative_brightness"].to_numpy(dtype=float),
        primary_use,
    )
    injections, injection_summary = run_ground_injections(final_curve)

    catalog_matches = catalog_match_measured_stars(
        table, args.wcs_image.resolve(), args.catalog_dir.resolve()
    )
    photometry = measured_star_photometry(table)
    catalog_matches = catalog_matches.merge(photometry, on="star", how="left")

    wcs = WCS(fits.getheader(args.wcs_image.resolve()))
    candidates = prepare_nearby_candidates(args.catalog_dir.resolve(), wcs)
    nearby_measurements: pd.DataFrame | None = None
    nearby_curves: pd.DataFrame | None = None
    if not args.skip_nearby_images:
        nearby_measurements, nearby_curves = measure_nearby_stars(
            candidates,
            table,
            final_curve,
            args.aligned_dir.resolve(),
            schedule,
        )

    ensemble_curves.to_csv(
        output_dir / "comparison_ensemble_light_curves.csv",
        index=False,
        float_format="%.10f",
    )
    ensemble_metrics.to_csv(
        output_dir / "comparison_ensemble_metrics.csv",
        index=False,
        float_format="%.10f",
    )
    noise.to_csv(
        output_dir / "noise_vs_bin_size.csv", index=False, float_format="%.10f"
    )
    injections.to_csv(
        output_dir / "ground_light_curve_injections.csv",
        index=False,
        float_format="%.10f",
    )
    injection_summary.to_csv(
        output_dir / "ground_injection_summary.csv",
        index=False,
        float_format="%.10f",
    )
    catalog_matches.to_csv(
        output_dir / "comparison_star_catalog_matches.csv",
        index=False,
        float_format="%.10f",
    )
    candidates.to_csv(
        output_dir / "nearby_star_catalog_candidates.csv",
        index=False,
        float_format="%.10f",
    )
    if nearby_measurements is not None and nearby_curves is not None:
        nearby_measurements.to_csv(
            output_dir / "nearby_star_image_measurements.csv",
            index=False,
            float_format="%.10f",
        )
        nearby_curves.to_csv(
            output_dir / "nearby_star_light_curves.csv",
            index=False,
            float_format="%.10f",
        )

    plot_ground_robustness(
        output_dir,
        ensemble_curves,
        ensemble_metrics,
        noise,
        injections,
        injection_summary,
        primary_use,
    )
    if nearby_measurements is not None:
        first_image = args.aligned_dir.resolve() / str(table.iloc[0]["Label"])
        plot_nearby_screen(
            output_dir,
            first_image,
            args.wcs_image.resolve(),
            candidates,
            nearby_measurements,
        )
    write_readme(
        output_dir,
        ensemble_metrics,
        noise,
        injection_summary,
        catalog_matches,
        nearby_measurements,
        schedule,
    )

    summary = {
        "target": "TOI-3505.01",
        "observation_interpretation": (
            "historical scheduled window under an EDT assumption; off transit "
            "under the current ephemeris"
        ),
        "historical_schedule": schedule,
        "adopted_ensemble": "sum of C2-C11 source counts",
        "alternative_ensembles_change_primary_curve": False,
        "ensemble_metrics": ensemble_metrics.to_dict(orient="records"),
        "noise_vs_bin_size": noise.to_dict(orient="records"),
        "injection": {
            "duration_hours": INJECTION_DURATION_HOURS,
            "depths_ppt": list(INJECTION_DEPTHS_PPT),
            "midpoint_trials": INJECTION_MIDPOINT_COUNT,
            "summary": injection_summary.to_dict(orient="records"),
        },
        "catalog": {
            "measured_stars": len(catalog_matches),
            "matched_within_3_arcsec": int(
                np.sum(catalog_matches["catalog_match_status"] == "matched")
            ),
            "ground_aperture_radius_arcsec": GROUND_SOURCE_RADIUS_ARCSEC,
        },
        "nearby_star_screen": {
            "radius_arcsec": NEARBY_RADIUS_ARCSEC,
            "bright_enough_catalog_candidates": len(candidates),
            "image_measurement_completed": nearby_measurements is not None,
            "historical_schedule_window_covered": True,
            "working_timezone": "America/New_York",
            "sources_cleared_by_conditional_screen": int(
                nearby_measurements["transit_relevant_clearance"].sum()
            )
            if nearby_measurements is not None
            else None,
            "source_apertures_overlapping_target": int(
                nearby_measurements["target_aperture_overlap"].sum()
            )
            if nearby_measurements is not None
            else None,
            "can_resolve_known_close_companion": False,
        },
    }
    (output_dir / "summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(f"Saved ground-data checks to {output_dir}")


if __name__ == "__main__":
    main()
