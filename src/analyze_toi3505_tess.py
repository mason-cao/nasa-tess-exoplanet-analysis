"""Measure the TOI-3505.01 signal in all four public TESS sectors.

This is the main space-based analysis for the project.  It uses QLP as the
uniform four-sector light-curve set, checks the signal with SPOC products where
they exist, runs a period search in each sector, measures individual events,
and checks whether Sector 54 contains data during the GMU observation.

The fitted shape is an exposure-integrated box, not a physical planet model.
That keeps the first-pass measurements transparent and prevents the 30-minute
Sector 14 data from appearing more precise than they are.
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
from astropy.timeseries import LombScargle

from toi3505_tess_tools import (
    LightCurveData,
    bls_period_search,
    event_cycles,
    fit_one_event,
    fit_repeating_box,
    grid_box_fit,
    independent_peaks,
    integrated_box_fraction,
    load_light_curve,
    phase_offset,
    robust_scatter,
    weighted_linear_ephemeris,
)
from toi3505_schedule import DEFAULT_SCHEDULE_RECORD, schedule_context


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DATA_DIR = ROOT / "data" / "tess" / "toi3505" / "light_curves"
DEFAULT_OUTPUT_DIR = ROOT / "outputs" / "toi3505_tess_analysis"
DEFAULT_GROUND_CURVE = (
    ROOT
    / "outputs"
    / "toi3505_final_candidate"
    / "TOI_3505.01_2022-07-22_R_light_curve.csv"
)

SECTORS = (14, 41, 54, 81)
PERIOD_DAYS = 2.9151556
PERIOD_ERROR_DAYS = 0.0000117
EPOCH_BJD = 2459793.534385
EPOCH_ERROR_DAYS = 0.0020787
CATALOG_DEPTH = 0.002910
CATALOG_DEPTH_ERROR = 0.000196
CATALOG_DURATION_DAYS = 2.004 / 24.0
FIT_WINDOW_DAYS = 0.22
GROUND_START_BJD = 2459782.598234811
GROUND_END_BJD = 2459782.809458706
RANDOM_SEED = 3505


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=DEFAULT_DATA_DIR)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--ground-curve", type=Path, default=DEFAULT_GROUND_CURVE)
    parser.add_argument(
        "--schedule-record", type=Path, default=DEFAULT_SCHEDULE_RECORD
    )
    return parser.parse_args()


def light_curve_files(data_dir: Path) -> list[Path]:
    files = sorted(data_dir.rglob("*lc.fits"))
    if len(files) < 8:
        raise FileNotFoundError(
            f"Expected at least eight QLP/SPOC light curves under {data_dir}; found {len(files)}"
        )
    return files


def load_primary_curves(paths: list[Path]) -> dict[int, LightCurveData]:
    curves: dict[int, LightCurveData] = {}
    for path in paths:
        if "hlsp_qlp" not in path.name.lower():
            continue
        curve = load_light_curve(path)
        if curve.sector in curves:
            raise RuntimeError(f"More than one QLP light curve for Sector {curve.sector}")
        curves[curve.sector] = curve
    if set(curves) != set(SECTORS):
        raise RuntimeError(f"QLP sectors found: {sorted(curves)}; expected {SECTORS}")
    return curves


def phase_bins(
    curve: LightCurveData,
    epoch_bjd: float,
    period_days: float,
    half_width_days: float = 0.16,
) -> pd.DataFrame:
    phase = phase_offset(curve.time_bjd, period_days, epoch_bjd)
    width = max(curve.cadence_days, 10.0 / 1440.0)
    keep = curve.good & (np.abs(phase) <= half_width_days)
    bin_number = np.floor((phase[keep] + half_width_days) / width).astype(int)
    rows = []
    for number in np.unique(bin_number):
        values = curve.flux[keep][bin_number == number]
        phases = phase[keep][bin_number == number]
        if len(values) == 0:
            continue
        scatter = robust_scatter(values)
        uncertainty = (
            scatter / np.sqrt(len(values))
            if len(values) > 1 and np.isfinite(scatter)
            else float(np.nanmedian(curve.flux_error[keep][bin_number == number]))
        )
        rows.append(
            {
                "phase_days": float(np.median(phases)),
                "flux": float(np.median(values)),
                "uncertainty": uncertainty,
                "points": len(values),
            }
        )
    return pd.DataFrame(rows)


def save_clean_curves(
    curves: dict[int, LightCurveData], output_dir: Path
) -> None:
    clean_dir = output_dir / "clean_light_curves"
    clean_dir.mkdir(parents=True, exist_ok=True)
    for sector, curve in curves.items():
        table = pd.DataFrame(
            {
                "time_bjd_tdb": curve.time_bjd,
                "time_btjd": curve.time_bjd - 2457000.0,
                "normalized_flux": curve.flux,
                "normalized_flux_error": curve.flux_error,
                "quality": curve.quality,
                "quality_zero": curve.good,
                "phase_days_current_ephemeris": phase_offset(
                    curve.time_bjd, PERIOD_DAYS, EPOCH_BJD
                ),
                "cycle_current_ephemeris": event_cycles(
                    curve.time_bjd, PERIOD_DAYS, EPOCH_BJD
                ),
            }
        )
        table.to_csv(clean_dir / f"sector_{sector}_qlp.csv", index=False)


def measure_sectors(
    curves: dict[int, LightCurveData], output_dir: Path
) -> tuple[pd.DataFrame, dict[int, dict[str, object]]]:
    durations = np.linspace(1.3 / 24.0, 2.8 / 24.0, 31)
    offsets = np.linspace(-0.035, 0.035, 71)
    rows: list[dict[str, object]] = []
    details: dict[int, dict[str, object]] = {}
    grid_dir = output_dir / "fit_grids"
    grid_dir.mkdir(parents=True, exist_ok=True)
    for sector in SECTORS:
        curve = curves[sector]
        best, chi2_grid, bounds = grid_box_fit(
            curve,
            period_days=PERIOD_DAYS,
            epoch_bjd=EPOCH_BJD,
            durations_days=durations,
            offsets_days=offsets,
            window_days=FIT_WINDOW_DAYS,
        )
        np.savez_compressed(
            grid_dir / f"sector_{sector}_box_grid.npz",
            durations_days=durations,
            offsets_days=offsets,
            chi2=chi2_grid,
        )
        odd_cycles = set(
            int(value)
            for value in np.unique(
                event_cycles(curve.time_bjd[curve.good], PERIOD_DAYS, EPOCH_BJD)
            )
            if int(value) % 2 != 0
        )
        even_cycles = set(
            int(value)
            for value in np.unique(
                event_cycles(curve.time_bjd[curve.good], PERIOD_DAYS, EPOCH_BJD)
            )
            if int(value) % 2 == 0
        )
        odd = fit_repeating_box(
            curve,
            PERIOD_DAYS,
            EPOCH_BJD,
            best.duration_days,
            best.time_offset_days,
            allowed_cycles=odd_cycles,
        )
        even = fit_repeating_box(
            curve,
            PERIOD_DAYS,
            EPOCH_BJD,
            best.duration_days,
            best.time_offset_days,
            allowed_cycles=even_cycles,
        )
        secondary = fit_repeating_box(
            curve,
            PERIOD_DAYS,
            EPOCH_BJD + PERIOD_DAYS / 2.0,
            best.duration_days,
            0.0,
        )
        odd_even_difference = odd.depth - even.depth
        odd_even_error = float(
            np.hypot(odd.depth_error, even.depth_error)
        )
        row = {
            "sector": sector,
            "pipeline": curve.pipeline,
            "cadence_minutes": curve.cadence_days * 1440.0,
            "quality_zero_points": int(curve.good.sum()),
            "events_in_sector_fit": best.events,
            "depth_ppt": best.depth * 1000.0,
            "depth_error_ppt": best.depth_error * 1000.0,
            "depth_snr": best.depth / best.depth_error,
            "duration_hours": best.duration_days * 24.0,
            "duration_low_hours": bounds["duration_low_days"] * 24.0,
            "duration_high_hours": bounds["duration_high_days"] * 24.0,
            "midpoint_offset_minutes": best.time_offset_days * 1440.0,
            "midpoint_offset_low_minutes": bounds["offset_low_days"] * 1440.0,
            "midpoint_offset_high_minutes": bounds["offset_high_days"] * 1440.0,
            "residual_scatter_ppt": best.residual_scatter * 1000.0,
            "reduced_chi2": best.reduced_chi2,
            "odd_depth_ppt": odd.depth * 1000.0,
            "odd_depth_error_ppt": odd.depth_error * 1000.0,
            "even_depth_ppt": even.depth * 1000.0,
            "even_depth_error_ppt": even.depth_error * 1000.0,
            "odd_even_difference_ppt": odd_even_difference * 1000.0,
            "odd_even_difference_sigma": (
                odd_even_difference / odd_even_error if odd_even_error > 0 else np.nan
            ),
            "phase_0_5_depth_ppt": secondary.depth * 1000.0,
            "phase_0_5_depth_error_ppt": secondary.depth_error * 1000.0,
            "phase_0_5_snr": secondary.depth / secondary.depth_error,
            "qlp_crowding_header": curve.crowdsap,
            "depth_interpretation": "observed depth in normalized QLP aperture flux; no CROWDSAP keyword present",
        }
        rows.append(row)
        details[sector] = {
            "curve": curve,
            "fit": best,
            "bounds": bounds,
            "odd": odd,
            "even": even,
            "secondary": secondary,
        }
        print(
            f"Sector {sector}: depth {row['depth_ppt']:.3f} +/- "
            f"{row['depth_error_ppt']:.3f} ppt; duration {row['duration_hours']:.2f} h"
        )
    table = pd.DataFrame(rows)
    table.to_csv(output_dir / "sector_measurements.csv", index=False)
    return table, details


def run_period_searches(
    curves: dict[int, LightCurveData], output_dir: Path
) -> tuple[pd.DataFrame, dict[int, dict[str, object]]]:
    rows: list[dict[str, object]] = []
    details: dict[int, dict[str, object]] = {}
    for sector in SECTORS:
        print(f"Running Sector {sector} box-least-squares period search...")
        result = bls_period_search(curves[sector])
        periods = result["period_grid_days"]
        power = result["power"]
        expected_index = int(np.argmin(np.abs(periods - PERIOD_DAYS)))
        peaks = independent_peaks(periods, power)
        row: dict[str, object] = {
            "sector": sector,
            "best_period_days": result["best_period_days"],
            "best_duration_hours": float(result["best_duration_days"]) * 24.0,
            "best_depth_ppt": float(result["best_depth"]) * 1000.0,
            "best_depth_error_ppt": float(result["best_depth_error"]) * 1000.0,
            "best_snr": result["best_power"],
            "catalog_period_days": PERIOD_DAYS,
            "snr_at_catalog_period": float(power[expected_index]),
            "best_period_difference_percent": 100.0
            * (float(result["best_period_days"]) - PERIOD_DAYS)
            / PERIOD_DAYS,
        }
        for index, (period, peak_power) in enumerate(peaks, start=1):
            row[f"peak_{index}_period_days"] = period
            row[f"peak_{index}_snr"] = peak_power
        rows.append(row)
        details[sector] = result
    table = pd.DataFrame(rows)
    table.to_csv(output_dir / "period_search.csv", index=False)
    return table, details


def measure_event_times(
    curves: dict[int, LightCurveData],
    sector_details: dict[int, dict[str, object]],
    output_dir: Path,
    schedule_record: Path,
) -> tuple[pd.DataFrame, dict[str, object]]:
    offsets = np.linspace(-0.07, 0.07, 281)
    rows: list[dict[str, object]] = []
    for sector in SECTORS:
        curve = curves[sector]
        duration = sector_details[sector]["fit"].duration_days
        good_times = curve.time_bjd[curve.good]
        first_cycle = int(np.floor((good_times.min() - EPOCH_BJD) / PERIOD_DAYS)) - 1
        last_cycle = int(np.ceil((good_times.max() - EPOCH_BJD) / PERIOD_DAYS)) + 1
        for cycle in range(first_cycle, last_cycle + 1):
            predicted = EPOCH_BJD + cycle * PERIOD_DAYS
            if predicted < good_times.min() - 0.1 or predicted > good_times.max() + 0.1:
                continue
            try:
                row = fit_one_event(
                    curve,
                    PERIOD_DAYS,
                    EPOCH_BJD,
                    cycle,
                    duration,
                    offsets,
                )
            except RuntimeError:
                continue
            row["used_in_ephemeris"] = bool(
                row["full_local_window"]
                and row["depth_snr"] >= 2.5
                and row["time_error_days"] <= 0.03
            )
            rows.append(row)
    events = pd.DataFrame(rows).sort_values(["measured_bjd", "sector"])
    events.to_csv(output_dir / "event_times.csv", index=False)

    selected = events[events["used_in_ephemeris"]]
    all_fit = weighted_linear_ephemeris(
        selected["cycle"].to_numpy(),
        selected["measured_bjd"].to_numpy(),
        selected["time_error_days"].to_numpy(),
    )
    sector14 = selected[selected["sector"] == 14]
    sector14_fit = weighted_linear_ephemeris(
        sector14["cycle"].to_numpy(),
        sector14["measured_bjd"].to_numpy(),
        sector14["time_error_days"].to_numpy(),
    )
    ground_midpoint = 0.5 * (GROUND_START_BJD + GROUND_END_BJD)
    ground_cycle = int(round((ground_midpoint - EPOCH_BJD) / PERIOD_DAYS))

    def prediction(fit: dict[str, float], cycle: int) -> dict[str, float]:
        time = fit["epoch_bjd"] + cycle * fit["period_days"]
        variance = (
            fit["epoch_error_days"] ** 2
            + cycle**2 * fit["period_error_days"] ** 2
            + 2.0 * cycle * fit["epoch_period_covariance_days2"]
        )
        return {
            "cycle": cycle,
            "bjd": time,
            "uncertainty_minutes": np.sqrt(max(0.0, variance)) * 1440.0,
        }

    current_time = EPOCH_BJD + ground_cycle * PERIOD_DAYS
    current_error = np.sqrt(
        EPOCH_ERROR_DAYS**2 + ground_cycle**2 * PERIOD_ERROR_DAYS**2
    )
    sector14_prediction = prediction(sector14_fit, ground_cycle)
    all_sector_prediction = prediction(all_fit, ground_cycle)
    historical_schedule = schedule_context(
        schedule_record,
        observation_start_bjd=GROUND_START_BJD,
        observation_end_bjd=GROUND_END_BJD,
    )
    historical_working = historical_schedule["working_interpretation"]
    assert isinstance(historical_working, dict)
    historical_midpoint = float(historical_working["event_midpoint_bjd_tdb"])
    schedule_row = historical_schedule["row"]
    assert isinstance(schedule_row, dict)
    schedule_period = float(schedule_row["Orbital Period"])
    historical_schedule["comparison_with_ephemerides"] = {
        "schedule_period_days": schedule_period,
        "schedule_minus_current_period_seconds": (
            schedule_period - PERIOD_DAYS
        )
        * 86400.0,
        "schedule_midpoint_minus_current_prediction_hours": (
            historical_midpoint - current_time
        )
        * 24.0,
        "schedule_midpoint_minus_sector14_prediction_hours": (
            historical_midpoint - float(sector14_prediction["bjd"])
        )
        * 24.0,
        "schedule_midpoint_minus_four_sector_prediction_hours": (
            historical_midpoint - float(all_sector_prediction["bjd"])
        )
        * 24.0,
        "period_interpretation": (
            "The schedule period is close to the later periods; without the "
            "schedule epoch, the large event-time difference cannot be assigned "
            "to period drift alone."
        ),
    }
    ephemeris = {
        "selection_rule": "full local window, box-depth SNR at least 2.5, timing error at most 0.03 day",
        "selected_event_count": int(len(selected)),
        "catalog": {
            "epoch_bjd": EPOCH_BJD,
            "epoch_error_days": EPOCH_ERROR_DAYS,
            "period_days": PERIOD_DAYS,
            "period_error_days": PERIOD_ERROR_DAYS,
            "epoch_period_covariance": "not available in the local catalog snapshot; zero used only for this displayed propagation",
            "ground_cycle": ground_cycle,
            "ground_prediction_bjd": current_time,
            "ground_prediction_uncertainty_minutes": current_error * 1440.0,
        },
        "sector_14_box_timing_fit": sector14_fit,
        "sector_14_prediction_at_ground_cycle": sector14_prediction,
        "four_sector_box_timing_fit": all_fit,
        "four_sector_prediction_at_ground_cycle": all_sector_prediction,
        "historical_schedule": historical_schedule,
        "warning": (
            "The scheduling row recovers a historical ingress/egress window but "
            "not its epoch, uncertainty, time-zone field, or prediction source. "
            "The TESS box timings remain preliminary checks."
        ),
    }
    (output_dir / "ephemeris_checks.json").write_text(
        json.dumps(ephemeris, indent=2) + "\n", encoding="utf-8"
    )
    return events, ephemeris


def pipeline_comparison(
    paths: list[Path], output_dir: Path
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    offsets = np.linspace(-0.035, 0.035, 71)
    for path in paths:
        for series in ("corrected", "raw"):
            curve = load_light_curve(path, series=series)
            best, _, bounds = grid_box_fit(
                curve,
                PERIOD_DAYS,
                EPOCH_BJD,
                durations_days=np.array([CATALOG_DURATION_DAYS]),
                offsets_days=offsets,
            )
            rows.append(
                {
                    "sector": curve.sector,
                    "pipeline": curve.pipeline,
                    "series": series,
                    "flux_column": curve.flux_name,
                    "cadence_minutes": curve.cadence_days * 1440.0,
                    "depth_ppt_fixed_2_004h": best.depth * 1000.0,
                    "depth_error_ppt": best.depth_error * 1000.0,
                    "midpoint_offset_minutes": best.time_offset_days * 1440.0,
                    "offset_low_minutes": bounds["offset_low_days"] * 1440.0,
                    "offset_high_minutes": bounds["offset_high_days"] * 1440.0,
                    "crowdsap": curve.crowdsap,
                    "flfrcsap": curve.flfrcsap,
                    "same_observation_note": "pipeline/extraction check; not an independent astronomical observation",
                    "file": path.name,
                }
            )
    table = pd.DataFrame(rows).sort_values(
        ["sector", "pipeline", "series", "cadence_minutes"]
    )
    table.to_csv(output_dir / "pipeline_comparison.csv", index=False)
    return table


def segment_detrend(
    curve: LightCurveData, transit_duration_days: float
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Remove a quadratic in each continuous section for a variability screen."""
    good_indices = np.flatnonzero(curve.good)
    time = curve.time_bjd[good_indices]
    flux = curve.flux[good_indices]
    phase = phase_offset(time, PERIOD_DAYS, EPOCH_BJD)
    transit_free = np.abs(phase) > 1.5 * transit_duration_days
    breaks = np.flatnonzero(np.diff(time) > 0.3) + 1
    sections = np.split(np.arange(len(time)), breaks)
    corrected = np.full(len(time), np.nan)
    for section in sections:
        fit_indices = section[transit_free[section]]
        if len(fit_indices) < 20:
            continue
        center = np.median(time[section])
        x_fit = time[fit_indices] - center
        coefficients = np.polyfit(x_fit, flux[fit_indices], deg=2)
        baseline = np.polyval(coefficients, time[section] - center)
        corrected[section] = flux[section] / baseline
    keep = transit_free & np.isfinite(corrected)
    return time[keep], corrected[keep], curve.flux_error[good_indices][keep]


def variability_screen(
    primary_paths: dict[int, Path],
    sector_details: dict[int, dict[str, object]],
    output_dir: Path,
) -> tuple[pd.DataFrame, dict[int, dict[str, np.ndarray | float]]]:
    rows = []
    details: dict[int, dict[str, np.ndarray | float]] = {}
    for sector in SECTORS:
        raw_curve = load_light_curve(primary_paths[sector], series="raw")
        duration = sector_details[sector]["fit"].duration_days
        time, flux, error = segment_detrend(raw_curve, duration)
        frequency, power = LombScargle(time, flux, error).autopower(
            minimum_frequency=1.0 / 13.0,
            maximum_frequency=10.0,
            samples_per_peak=8,
        )
        best = int(np.nanargmax(power))
        best_period = 1.0 / float(frequency[best])
        false_alarm = float(
            LombScargle(time, flux, error).false_alarm_probability(power[best])
        )
        rows.append(
            {
                "sector": sector,
                "strongest_residual_period_days": best_period,
                "lomb_scargle_power": float(power[best]),
                "formal_false_alarm_probability": false_alarm,
                "points": len(time),
                "interpretation": "variability screen only; QLP/systematics processing can alter long periods",
            }
        )
        details[sector] = {
            "period_days": 1.0 / frequency,
            "power": power,
            "best_period_days": best_period,
        }
    table = pd.DataFrame(rows)
    table.to_csv(output_dir / "variability_screen.csv", index=False)
    return table, details


def make_full_sector_plot(
    curves: dict[int, LightCurveData], output_dir: Path
) -> None:
    figure, axes = plt.subplots(4, 1, figsize=(11, 10), sharey=True)
    for ax, sector in zip(axes, SECTORS):
        curve = curves[sector]
        good = curve.good
        time_btjd = curve.time_bjd - 2457000.0
        ax.scatter(time_btjd[good], curve.flux[good], s=3, alpha=0.5, color="#244f78")
        first_cycle = int(np.floor((curve.time_bjd[good].min() - EPOCH_BJD) / PERIOD_DAYS))
        last_cycle = int(np.ceil((curve.time_bjd[good].max() - EPOCH_BJD) / PERIOD_DAYS))
        for cycle in range(first_cycle, last_cycle + 1):
            midpoint = EPOCH_BJD + cycle * PERIOD_DAYS - 2457000.0
            ax.axvspan(
                midpoint - CATALOG_DURATION_DAYS / 2.0,
                midpoint + CATALOG_DURATION_DAYS / 2.0,
                color="#c47a36",
                alpha=0.12,
                linewidth=0,
            )
        ax.set_ylabel(f"S{sector}\nrelative flux")
        ax.text(
            0.99,
            0.90,
            f"{curve.cadence_days * 1440:.0f}-minute cadence",
            transform=ax.transAxes,
            ha="right",
            va="top",
            fontsize=9,
        )
        ax.grid(alpha=0.16)
    axes[-1].set_xlabel("Time (BJD TDB − 2457000)")
    figure.suptitle("TOI-3505.01 QLP light curves", fontsize=15)
    figure.tight_layout(rect=(0, 0, 1, 0.97))
    figure.savefig(output_dir / "01_full_sector_light_curves.png", dpi=220)
    figure.savefig(output_dir / "01_full_sector_light_curves.svg")
    plt.close(figure)


def make_phase_plot(
    curves: dict[int, LightCurveData],
    details: dict[int, dict[str, object]],
    output_dir: Path,
) -> None:
    figure, axes = plt.subplots(2, 2, figsize=(11, 8.2), sharex=True, sharey=True)
    model_phase = np.linspace(-0.16, 0.16, 1200)
    for ax, sector in zip(axes.flat, SECTORS):
        curve = curves[sector]
        fit = details[sector]["fit"]
        phase = phase_offset(curve.time_bjd, PERIOD_DAYS, EPOCH_BJD)
        keep = curve.good & (np.abs(phase) <= 0.16)
        ax.scatter(
            phase[keep] * 24.0,
            curve.flux[keep],
            s=5,
            alpha=0.18,
            color="#6d7780",
            label="individual measurements",
        )
        bins = phase_bins(curve, EPOCH_BJD, PERIOD_DAYS)
        ax.errorbar(
            bins["phase_days"] * 24.0,
            bins["flux"],
            yerr=bins["uncertainty"],
            fmt="o",
            markersize=3.2,
            linewidth=0.8,
            color="#244f78",
            label="time bins",
        )
        model_fraction = integrated_box_fraction(
            model_phase - fit.time_offset_days,
            fit.duration_days,
            curve.cadence_days,
        )
        ax.plot(
            model_phase * 24.0,
            1.0 - fit.depth * model_fraction,
            color="#b35a2a",
            linewidth=1.8,
            label="box estimate",
        )
        ax.axhline(1.0, color="0.35", linewidth=0.8)
        ax.set_title(
            f"Sector {sector}: {fit.depth * 1000:.2f} +/- {fit.depth_error * 1000:.2f} ppt"
        )
        ax.grid(alpha=0.16)
    for ax in axes[:, 0]:
        ax.set_ylabel("Relative flux")
    for ax in axes[-1, :]:
        ax.set_xlabel("Hours from catalog midpoint")
    axes[0, 0].legend(loc="lower left", fontsize=8)
    figure.suptitle("TOI-3505.01 in each TESS sector", fontsize=15)
    figure.tight_layout(rect=(0, 0, 1, 0.96))
    figure.savefig(output_dir / "02_phase_folded_sectors.png", dpi=240)
    figure.savefig(output_dir / "02_phase_folded_sectors.svg")
    plt.close(figure)


def make_period_plot(
    details: dict[int, dict[str, object]], output_dir: Path
) -> None:
    figure, axes = plt.subplots(2, 2, figsize=(11, 7.6), sharex=True)
    for ax, sector in zip(axes.flat, SECTORS):
        result = details[sector]
        ax.plot(
            result["period_grid_days"], result["power"], color="#244f78", linewidth=1.0
        )
        ax.axvline(PERIOD_DAYS, color="#b35a2a", linestyle="--", linewidth=1.2)
        ax.axvline(PERIOD_DAYS / 2.0, color="0.45", linestyle=":", linewidth=1.0)
        ax.set_title(
            f"Sector {sector}: strongest period {result['best_period_days']:.4f} d"
        )
        ax.set_ylabel("Box-search S/N")
        ax.grid(alpha=0.16)
    for ax in axes[-1, :]:
        ax.set_xlabel("Trial period (days)")
    figure.suptitle("Independent period search in each sector", fontsize=15)
    figure.tight_layout(rect=(0, 0, 1, 0.96))
    figure.savefig(output_dir / "03_period_search.png", dpi=230)
    figure.savefig(output_dir / "03_period_search.svg")
    plt.close(figure)


def make_depth_plot(measurements: pd.DataFrame, output_dir: Path) -> None:
    figure, ax = plt.subplots(figsize=(8.6, 5.2))
    x = np.arange(len(measurements))
    ax.errorbar(
        x,
        measurements["depth_ppt"],
        yerr=measurements["depth_error_ppt"],
        fmt="o",
        markersize=7,
        capsize=4,
        color="#244f78",
        label="QLP box estimates",
    )
    ax.axhspan(
        (CATALOG_DEPTH - CATALOG_DEPTH_ERROR) * 1000.0,
        (CATALOG_DEPTH + CATALOG_DEPTH_ERROR) * 1000.0,
        color="#c47a36",
        alpha=0.18,
        label="catalog depth range",
    )
    ax.axhline(CATALOG_DEPTH * 1000.0, color="#b35a2a", linewidth=1.2)
    ax.set_xticks(x, [f"Sector {sector}" for sector in measurements["sector"]])
    ax.set_ylabel("Observed depth (ppt)")
    ax.set_title("TOI-3505.01 depth by sector")
    ax.legend()
    ax.grid(axis="y", alpha=0.18)
    figure.tight_layout()
    figure.savefig(output_dir / "04_depth_by_sector.png", dpi=230)
    figure.savefig(output_dir / "04_depth_by_sector.svg")
    plt.close(figure)


def make_event_time_plot(events: pd.DataFrame, output_dir: Path) -> None:
    figure, ax = plt.subplots(figsize=(10.5, 5.5))
    colors = {14: "#355f8a", 41: "#3f7d64", 54: "#b16a2f", 81: "#7b5a92"}
    for sector in SECTORS:
        rows = events[(events["sector"] == sector) & events["used_in_ephemeris"]]
        ax.errorbar(
            rows["measured_bjd"] - 2457000.0,
            rows["observed_minus_calculated_minutes"],
            yerr=rows["time_error_days"] * 1440.0,
            fmt="o",
            markersize=5,
            capsize=2,
            color=colors[sector],
            label=f"Sector {sector}",
        )
    ax.axhline(0.0, color="0.3", linewidth=1.0)
    ax.set_xlabel("Measured midpoint (BJD TDB − 2457000)")
    ax.set_ylabel("Observed minus catalog prediction (minutes)")
    ax.set_title("TOI-3505.01 box timing check")
    ax.legend(ncol=4)
    ax.grid(alpha=0.18)
    figure.tight_layout()
    figure.savefig(output_dir / "05_event_times.png", dpi=230)
    figure.savefig(output_dir / "05_event_times.svg")
    plt.close(figure)


def make_timing_window_plot(
    ephemeris: dict[str, object], output_dir: Path
) -> None:
    ground_cycle = ephemeris["catalog"]["ground_cycle"]
    current_times = [
        EPOCH_BJD + cycle * PERIOD_DAYS
        for cycle in (ground_cycle, ground_cycle + 1)
    ]
    sector14_prediction = ephemeris["sector_14_prediction_at_ground_cycle"]
    historical = ephemeris["historical_schedule"]
    assert isinstance(historical, dict)
    historical_working = historical["working_interpretation"]
    assert isinstance(historical_working, dict)
    historical_times = historical_working["times"]
    assert isinstance(historical_times, dict)
    historical_ingress = float(historical_times["ingress"]["bjd_tdb"])
    historical_egress = float(historical_times["egress"]["bjd_tdb"])
    display_zero = 2459780.0
    display_current_times = [value - display_zero for value in current_times]
    display_ground_start = GROUND_START_BJD - display_zero
    display_ground_end = GROUND_END_BJD - display_zero
    display_historical_ingress = historical_ingress - display_zero
    display_historical_egress = historical_egress - display_zero
    minimum = current_times[0] - 0.20
    maximum = current_times[1] + 0.20
    figure, ax = plt.subplots(figsize=(11, 3.8))
    ax.hlines(
        1.0,
        display_ground_start,
        display_ground_end,
        color="#244f78",
        linewidth=9,
    )
    ax.text(
        (display_ground_start + display_ground_end) / 2.0,
        1.07,
        "GMU usable images",
        ha="center",
        color="#244f78",
    )
    for index, midpoint in enumerate(display_current_times):
        ax.axvspan(
            midpoint - CATALOG_DURATION_DAYS / 2.0,
            midpoint + CATALOG_DURATION_DAYS / 2.0,
            color="#c47a36",
            alpha=0.24,
        )
        ax.axvline(midpoint, color="#b35a2a", linewidth=1.2)
        ax.text(
            midpoint,
            0.88,
            "catalog prediction" if index == 0 else "next catalog prediction",
            ha="center",
            va="top",
            fontsize=9,
        )
    s14_time = sector14_prediction["bjd"]
    s14_error = sector14_prediction["uncertainty_minutes"] / 1440.0
    ax.errorbar(
        s14_time - display_zero,
        0.55,
        xerr=s14_error,
        fmt="D",
        capsize=4,
        color="#4e7658",
        label="Sector 14-only box ephemeris",
    )
    ax.hlines(
        0.76,
        display_historical_ingress,
        display_historical_egress,
        color="#8a5a9b",
        linewidth=9,
        label="2022 schedule window (EDT assumed)",
    )
    ax.vlines(
        [display_historical_ingress, display_historical_egress],
        0.71,
        0.81,
        color="#6e407f",
        linewidth=1.0,
    )
    ax.set_xlim(minimum - display_zero, maximum - display_zero)
    ax.set_ylim(0.35, 1.22)
    ax.set_yticks([])
    ax.set_xlabel("Barycentric Julian Date (TDB) − 2459780")
    ax.set_title("GMU observation window and timing references")
    ax.legend(loc="lower right")
    ax.grid(axis="x", alpha=0.18)
    figure.tight_layout()
    figure.savefig(output_dir / "06_gmu_timing_window.png", dpi=240)
    figure.savefig(output_dir / "06_gmu_timing_window.svg")
    plt.close(figure)


def make_ground_gap_plot(
    sector54: LightCurveData, ground_path: Path, output_dir: Path
) -> None:
    if not ground_path.exists():
        return
    ground = pd.read_csv(ground_path)
    used = ground["Used_in_Plot"].astype(bool).to_numpy()
    ground_time = ground["BJD_TDB"].to_numpy(dtype=float)
    ground_flux = ground["Relative_Brightness"].to_numpy(dtype=float)
    ground_error = ground["Flux_Error"].to_numpy(dtype=float)
    ground_scale = np.nanmedian(ground_flux[used])
    ground_flux /= ground_scale
    ground_error /= ground_scale
    tess_keep = (
        sector54.good
        & (sector54.time_bjd >= ground_time.min())
        & (sector54.time_bjd <= ground_time.max())
    )
    before = sector54.time_bjd[sector54.good & (sector54.time_bjd < ground_time.min())]
    after = sector54.time_bjd[sector54.good & (sector54.time_bjd > ground_time.max())]
    nearest_before = float(np.max(before)) if len(before) else None
    nearest_after = float(np.min(after)) if len(after) else None
    nearby = (
        sector54.good
        & (sector54.time_bjd >= ground_time.min() - 0.5)
        & (sector54.time_bjd <= ground_time.max() + 1.0)
    )
    tess_scale = np.nanmedian(sector54.flux[nearby])
    tess_flux = sector54.flux / tess_scale
    hours = (ground_time - ground_time.min()) * 24.0
    figure, axes = plt.subplots(2, 1, figsize=(10.5, 7.2))
    axes[0].errorbar(
        hours[used],
        ground_flux[used],
        yerr=ground_error[used],
        fmt=".",
        markersize=4,
        linewidth=0.5,
        alpha=0.75,
        color="#244f78",
    )
    axes[0].scatter(
        hours[~used],
        ground_flux[~used],
        facecolors="none",
        edgecolors="#a34d3f",
        s=22,
        label="image-quality exclusions",
    )
    axes[0].set_ylabel("GMU relative flux")
    axes[0].legend(loc="best", fontsize=8)
    axes[1].scatter(
        sector54.time_bjd[nearby],
        tess_flux[nearby],
        s=12,
        alpha=0.7,
        color="#b16a2f",
    )
    axes[1].axvspan(
        ground_time.min(),
        ground_time.max(),
        color="#244f78",
        alpha=0.17,
        label="GMU observing window",
    )
    axes[1].set_xlim(ground_time.min() - 0.5, ground_time.max() + 1.0)
    axes[1].set_ylabel("TESS relative flux")
    axes[1].set_xlabel("Time (BJD TDB)")
    axes[1].legend(loc="best", fontsize=8)
    for ax in axes:
        ax.axhline(1.0, color="0.35", linewidth=0.8)
        ax.grid(alpha=0.18)
    if tess_keep.sum() == 0:
        title = "The GMU observation falls in a TESS Sector 54 data gap"
    else:
        title = "GMU and TESS Sector 54 observing overlap"
    figure.suptitle(title, fontsize=15)
    figure.tight_layout(rect=(0, 0, 1, 0.96))
    figure.savefig(output_dir / "07_sector54_ground_data_gap.png", dpi=230)
    figure.savefig(output_dir / "07_sector54_ground_data_gap.svg")
    plt.close(figure)

    overlap = {
        "ground_start_bjd_tdb": float(ground_time.min()),
        "ground_end_bjd_tdb": float(ground_time.max()),
        "quality_zero_tess_points_during_ground_window": int(tess_keep.sum()),
        "nearest_quality_zero_tess_before_bjd_tdb": nearest_before,
        "hours_from_last_tess_point_to_ground_start": (
            (float(ground_time.min()) - nearest_before) * 24.0
            if nearest_before is not None
            else None
        ),
        "nearest_quality_zero_tess_after_bjd_tdb": nearest_after,
        "hours_from_ground_end_to_next_tess_point": (
            (nearest_after - float(ground_time.max())) * 24.0
            if nearest_after is not None
            else None
        ),
        "conclusion": (
            "No direct GMU/TESS flux comparison is possible because the GMU sequence is inside the Sector 54 data gap."
            if tess_keep.sum() == 0
            else "A direct overlap comparison is available."
        ),
    }
    (output_dir / "sector54_ground_overlap.json").write_text(
        json.dumps(overlap, indent=2) + "\n", encoding="utf-8"
    )


def make_variability_plot(
    details: dict[int, dict[str, np.ndarray | float]], output_dir: Path
) -> None:
    figure, axes = plt.subplots(2, 2, figsize=(11, 7.5), sharex=True)
    for ax, sector in zip(axes.flat, SECTORS):
        result = details[sector]
        order = np.argsort(result["period_days"])
        ax.plot(
            np.asarray(result["period_days"])[order],
            np.asarray(result["power"])[order],
            color="#244f78",
            linewidth=0.9,
        )
        ax.axvline(
            result["best_period_days"], color="#b35a2a", linestyle="--", linewidth=1.0
        )
        ax.set_xscale("log")
        ax.set_title(f"Sector {sector}: {result['best_period_days']:.3f} d peak")
        ax.set_ylabel("Lomb-Scargle power")
        ax.grid(alpha=0.16)
    for ax in axes[-1, :]:
        ax.set_xlabel("Trial period (days)")
    figure.suptitle("Out-of-transit variability screen", fontsize=15)
    figure.tight_layout(rect=(0, 0, 1, 0.96))
    figure.savefig(output_dir / "08_variability_screen.png", dpi=230)
    figure.savefig(output_dir / "08_variability_screen.svg")
    plt.close(figure)


def write_readme(
    output_dir: Path,
    measurements: pd.DataFrame,
    period_search: pd.DataFrame,
    events: pd.DataFrame,
    ephemeris: dict[str, object],
    pipeline_table: pd.DataFrame,
    variability: pd.DataFrame,
) -> None:
    lines = [
        "# TOI-3505.01 four-sector TESS analysis",
        "",
        "This folder contains a reproducible first-pass measurement of the public TESS signal in Sectors 14, 41, 54, and 81.",
        "",
        "## Method",
        "",
        "- QLP is the common four-sector data set.",
        "- Only finite measurements with QUALITY = 0 are used.",
        "- The fitted shape is a box integrated across each sector's exposure time.",
        "- Every event has its own straight local baseline.",
        "- SPOC and TESS-SPOC products are pipeline checks of the same observations, not extra observations.",
        "- QLP files do not contain a CROWDSAP keyword, so the QLP depths here are observed aperture depths and are not labeled dilution-corrected.",
        "",
        "## Sector measurements",
        "",
        "| Sector | Cadence (min) | Depth (ppt) | Duration (h) | Midpoint offset (min) | BLS period (d) |",
        "|---:|---:|---:|---:|---:|---:|",
    ]
    for _, row in measurements.iterrows():
        bls = period_search.loc[period_search["sector"] == row["sector"]].iloc[0]
        lines.append(
            f"| {int(row['sector'])} | {row['cadence_minutes']:.1f} | "
            f"{row['depth_ppt']:.3f} +/- {row['depth_error_ppt']:.3f} | "
            f"{row['duration_hours']:.2f} | {row['midpoint_offset_minutes']:.1f} | "
            f"{bls['best_period_days']:.5f} |"
        )
    selected = int(events["used_in_ephemeris"].sum())
    s14_prediction = ephemeris["sector_14_prediction_at_ground_cycle"]
    historical = ephemeris["historical_schedule"]
    assert isinstance(historical, dict)
    historical_working = historical["working_interpretation"]
    historical_comparison = historical["comparison_with_ephemerides"]
    assert isinstance(historical_working, dict)
    assert isinstance(historical_comparison, dict)
    historical_times = historical_working["times"]
    assert isinstance(historical_times, dict)
    lines.extend(
        [
            "",
            "## Timing interpretation",
            "",
            f"{selected} individual events pass the stated timing-quality rule. The Sector 14-only box ephemeris predicts the catalog cycle nearest the GMU night at BJD_TDB {s14_prediction['bjd']:.6f} +/- {s14_prediction['uncertainty_minutes']:.1f} minutes. This is a modern reconstruction from the public Sector 14 light curve.",
            "",
            f"The recovered 2022 schedule row lists ingress at 00:15 and egress at 01:54. Under the documented Eastern-time interpretation, those become BJD_TDB {float(historical_times['ingress']['bjd_tdb']):.6f}-{float(historical_times['egress']['bjd_tdb']):.6f}, inside the GMU images. Its midpoint is {float(historical_comparison['schedule_midpoint_minus_current_prediction_hours']):.2f} hours later than the nearest current-catalog prediction.",
            "",
            "The schedule period is close to the later measured periods, but the source row contains no epoch or timing uncertainty. The disagreement therefore cannot be attributed to period drift alone. The historical window is a recovered scheduling prediction, not proof that a transit occurred.",
            "",
            "The GMU sequence also falls entirely inside a Sector 54 TESS data gap. There are no quality-zero TESS points during the ground window, so a simultaneous flux comparison is not possible.",
            "",
            "## Limits",
            "",
            "- A box estimate does not measure limb darkening, impact parameter, or a physical planet radius.",
            "- The reported formal errors include the local fit scatter but do not include a complete astrophysical dilution model.",
            "- The variability periodograms are screens; QLP detrending and the simple quadratic section correction can alter long periods.",
            "- Difference imaging and custom-aperture tests are kept in the separate pixel-analysis output.",
            "- Official SPOC per-sector and combined-report values are kept in outputs/toi3505_data_validation. They are same-observation comparisons.",
            "",
            "## Output tables",
            "",
            "- `sector_measurements.csv`: depth, duration, timing, odd/even, and phase-0.5 checks.",
            "- `period_search.csv`: independent BLS result in each sector.",
            "- `event_times.csv`: one row per measurable transit window.",
            "- `ephemeris_checks.json`: Sector 14-only and four-sector box timing fits.",
            "- `pipeline_comparison.csv`: QLP, SPOC, and TESS-SPOC extraction checks.",
            "- `variability_screen.csv`: out-of-transit residual-period search.",
            "- `clean_light_curves/`: normalized QLP data with quality, phase, and cycle columns.",
            "- `../toi3505_data_validation/`: official Sector 54, Sector 81, and combined-sector SPOC report comparison.",
            "",
            f"Pipeline-comparison rows: {len(pipeline_table)}. Variability-screen rows: {len(variability)}.",
        ]
    )
    (output_dir / "README.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    args = parse_args()
    data_dir = args.data_dir.resolve()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    paths = light_curve_files(data_dir)
    curves = load_primary_curves(paths)
    primary_paths = {sector: curve.path for sector, curve in curves.items()}

    save_clean_curves(curves, output_dir)
    measurements, sector_details = measure_sectors(curves, output_dir)
    period_table, period_details = run_period_searches(curves, output_dir)
    events, ephemeris = measure_event_times(
        curves, sector_details, output_dir, args.schedule_record.resolve()
    )
    pipeline_table = pipeline_comparison(paths, output_dir)
    variability_table, variability_details = variability_screen(
        primary_paths, sector_details, output_dir
    )

    make_full_sector_plot(curves, output_dir)
    make_phase_plot(curves, sector_details, output_dir)
    make_period_plot(period_details, output_dir)
    make_depth_plot(measurements, output_dir)
    make_event_time_plot(events, output_dir)
    make_timing_window_plot(ephemeris, output_dir)
    make_ground_gap_plot(curves[54], args.ground_curve.resolve(), output_dir)
    make_variability_plot(variability_details, output_dir)

    summary = {
        "target": "TOI-3505.01",
        "tic_id": 390988385,
        "sectors": list(SECTORS),
        "catalog_ephemeris": {
            "period_days": PERIOD_DAYS,
            "period_error_days": PERIOD_ERROR_DAYS,
            "epoch_bjd_tdb": EPOCH_BJD,
            "epoch_error_days": EPOCH_ERROR_DAYS,
        },
        "catalog_depth_ppt": CATALOG_DEPTH * 1000.0,
        "catalog_duration_hours": CATALOG_DURATION_DAYS * 24.0,
        "quality_rule": "finite flux/error and QUALITY == 0",
        "model": "finite-exposure box with one local straight baseline per event",
        "random_seed_reserved_for_injection_tests": RANDOM_SEED,
        "sector_measurements": measurements.to_dict(orient="records"),
        "period_search": period_table.to_dict(orient="records"),
        "ephemeris_checks": ephemeris,
    }
    (output_dir / "analysis_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    write_readme(
        output_dir,
        measurements,
        period_table,
        events,
        ephemeris,
        pipeline_table,
        variability_table,
    )
    print(f"Saved the four-sector TESS analysis to {output_dir}")


if __name__ == "__main__":
    main()
