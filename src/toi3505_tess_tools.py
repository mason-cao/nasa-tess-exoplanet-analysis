"""Small, testable tools used by the TOI-3505.01 TESS analysis.

The measurements in this project deliberately use a transparent box model.
It is not a physical transit model: it estimates the depth, duration, and time
of a repeatable dip while accounting for the finite TESS exposure time and a
separate straight-line baseline around every event.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
from astropy.io import fits
from astropy.timeseries import BoxLeastSquares


@dataclass
class LightCurveData:
    path: Path
    sector: int
    pipeline: str
    flux_name: str
    time_bjd: np.ndarray
    flux: np.ndarray
    flux_error: np.ndarray
    quality: np.ndarray
    cadence_days: float
    crowdsap: float | None
    flfrcsap: float | None

    @property
    def good(self) -> np.ndarray:
        return (
            np.isfinite(self.time_bjd)
            & np.isfinite(self.flux)
            & np.isfinite(self.flux_error)
            & (self.flux_error > 0)
            & (self.quality == 0)
        )


@dataclass
class BoxFit:
    depth: float
    depth_error: float
    duration_days: float
    time_offset_days: float
    chi2: float
    reduced_chi2: float
    points: int
    events: int
    residual_scatter: float
    coefficients: np.ndarray
    model: np.ndarray
    used_indices: np.ndarray


def _header_value(headers: Iterable[fits.Header], name: str, default=None):
    for header in headers:
        if name in header:
            return header[name]
    return default


def pipeline_name(path: Path) -> str:
    name = path.name.lower()
    if "hlsp_qlp" in name:
        return "QLP"
    if "hlsp_tess-spoc" in name:
        return "TESS-SPOC"
    if name.endswith("_lc.fits") and name.startswith("tess"):
        return "SPOC"
    raise ValueError(f"Cannot identify the light-curve pipeline for {path.name}")


def load_light_curve(path: Path, series: str = "corrected") -> LightCurveData:
    """Load one QLP, TESS-SPOC, or SPOC light curve and normalize its flux."""
    path = Path(path)
    pipeline = pipeline_name(path)
    with fits.open(path, memmap=True) as hdul:
        primary = hdul[0].header.copy()
        table_header = hdul[1].header.copy()
        table = hdul[1].data
        names = set(table.columns.names)
        time = np.asarray(table["TIME"], dtype=float)
        quality = np.asarray(table["QUALITY"], dtype=np.int64)

        if pipeline == "QLP":
            if series == "raw":
                flux_name = "SAP_FLUX"
                error_name = (
                    "DET_FLUX_ERR" if "DET_FLUX_ERR" in names else "KSPSAP_FLUX_ERR"
                )
            else:
                flux_name = "DET_FLUX" if "DET_FLUX" in names else "KSPSAP_FLUX"
                error_name = (
                    "DET_FLUX_ERR" if "DET_FLUX_ERR" in names else "KSPSAP_FLUX_ERR"
                )
        else:
            if series == "raw":
                flux_name, error_name = "SAP_FLUX", "SAP_FLUX_ERR"
            else:
                flux_name, error_name = "PDCSAP_FLUX", "PDCSAP_FLUX_ERR"

        flux = np.asarray(table[flux_name], dtype=float)
        flux_error = np.asarray(table[error_name], dtype=float)

    headers = [primary, table_header]
    reference = float(_header_value(headers, "BJDREFI", 2457000.0)) + float(
        _header_value(headers, "BJDREFF", 0.0)
    )
    time_bjd = time + reference
    initial_good = (
        np.isfinite(time_bjd)
        & np.isfinite(flux)
        & np.isfinite(flux_error)
        & (flux_error > 0)
        & (quality == 0)
    )
    if initial_good.sum() < 20:
        raise RuntimeError(f"Too few good measurements in {path.name}")
    scale = float(np.nanmedian(flux[initial_good]))
    if not np.isfinite(scale) or scale <= 0:
        raise RuntimeError(f"Invalid median flux in {path.name}")

    return LightCurveData(
        path=path,
        sector=int(_header_value(headers, "SECTOR")),
        pipeline=pipeline,
        flux_name=flux_name,
        time_bjd=time_bjd,
        flux=flux / scale,
        flux_error=flux_error / scale,
        quality=quality,
        cadence_days=float(_header_value(headers, "TIMEDEL")),
        crowdsap=(
            float(_header_value(headers, "CROWDSAP"))
            if _header_value(headers, "CROWDSAP") is not None
            else None
        ),
        flfrcsap=(
            float(_header_value(headers, "FLFRCSAP"))
            if _header_value(headers, "FLFRCSAP") is not None
            else None
        ),
    )


def phase_offset(time_bjd: np.ndarray, period_days: float, epoch_bjd: float) -> np.ndarray:
    """Signed time from the nearest predicted event, in days."""
    time = np.asarray(time_bjd, dtype=float)
    return (time - epoch_bjd + period_days / 2.0) % period_days - period_days / 2.0


def event_cycles(time_bjd: np.ndarray, period_days: float, epoch_bjd: float) -> np.ndarray:
    time = np.asarray(time_bjd, dtype=float)
    cycles = np.zeros(time.shape, dtype=int)
    finite = np.isfinite(time)
    cycles[finite] = np.rint((time[finite] - epoch_bjd) / period_days).astype(int)
    return cycles


def integrated_box_fraction(
    phase_days: np.ndarray, duration_days: float, cadence_days: float
) -> np.ndarray:
    """Fraction of each exposure that overlaps a centered box-shaped event."""
    phase = np.asarray(phase_days, dtype=float)
    if duration_days <= 0 or cadence_days <= 0:
        raise ValueError("Duration and cadence must be positive")
    exposure_left = phase - cadence_days / 2.0
    exposure_right = phase + cadence_days / 2.0
    event_left = -duration_days / 2.0
    event_right = duration_days / 2.0
    overlap = np.minimum(exposure_right, event_right) - np.maximum(
        exposure_left, event_left
    )
    return np.clip(overlap / cadence_days, 0.0, 1.0)


def robust_scatter(values: np.ndarray) -> float:
    values = np.asarray(values, dtype=float)
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        return float("nan")
    center = float(np.median(finite))
    return float(1.4826 * np.median(np.abs(finite - center)))


def _solve_weighted_model(
    design: np.ndarray, flux: np.ndarray, error: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray, float, float]:
    weight_root = 1.0 / error
    weighted_design = design * weight_root[:, None]
    weighted_flux = flux * weight_root
    coefficients, _, _, _ = np.linalg.lstsq(
        weighted_design, weighted_flux, rcond=None
    )
    model = design @ coefficients
    residual = flux - model
    chi2 = float(np.sum(np.square(residual / error)))
    degrees_of_freedom = max(1, len(flux) - design.shape[1])
    reduced_chi2 = chi2 / degrees_of_freedom
    covariance = np.linalg.pinv(weighted_design.T @ weighted_design)
    covariance *= max(1.0, reduced_chi2)
    return coefficients, covariance, model, chi2, reduced_chi2


def fit_repeating_box(
    curve: LightCurveData,
    period_days: float,
    epoch_bjd: float,
    duration_days: float,
    time_offset_days: float = 0.0,
    window_days: float = 0.22,
    allowed_cycles: set[int] | None = None,
    extra_mask: np.ndarray | None = None,
) -> BoxFit:
    """Fit one depth plus a local line around each predicted event."""
    reference_phase = phase_offset(curve.time_bjd, period_days, epoch_bjd)
    cycles = event_cycles(curve.time_bjd, period_days, epoch_bjd)
    keep = curve.good & (np.abs(reference_phase) <= window_days)
    if allowed_cycles is not None:
        keep &= np.isin(cycles, list(allowed_cycles))
    if extra_mask is not None:
        keep &= np.asarray(extra_mask, dtype=bool)
    indices = np.flatnonzero(keep)
    if len(indices) < 10:
        raise RuntimeError("Too few points near the predicted events")

    selected_cycles = cycles[indices]
    useful_cycles = [
        int(cycle)
        for cycle in np.unique(selected_cycles)
        if np.count_nonzero(selected_cycles == cycle) >= 4
    ]
    keep_useful = np.isin(selected_cycles, useful_cycles)
    indices = indices[keep_useful]
    selected_cycles = cycles[indices]
    if not useful_cycles:
        raise RuntimeError("No event has enough points for a local baseline")

    local_phase = reference_phase[indices]
    columns: list[np.ndarray] = []
    for cycle in useful_cycles:
        event = (selected_cycles == cycle).astype(float)
        columns.append(event)
        columns.append(event * local_phase / window_days)
    shifted_phase = phase_offset(
        curve.time_bjd[indices], period_days, epoch_bjd + time_offset_days
    )
    box = integrated_box_fraction(
        shifted_phase, duration_days=duration_days, cadence_days=curve.cadence_days
    )
    columns.append(-box)
    design = np.column_stack(columns)
    coefficients, covariance, model, chi2, reduced_chi2 = _solve_weighted_model(
        design, curve.flux[indices], curve.flux_error[indices]
    )
    depth = float(coefficients[-1])
    depth_error = float(np.sqrt(max(0.0, covariance[-1, -1])))
    return BoxFit(
        depth=depth,
        depth_error=depth_error,
        duration_days=duration_days,
        time_offset_days=time_offset_days,
        chi2=chi2,
        reduced_chi2=reduced_chi2,
        points=len(indices),
        events=len(useful_cycles),
        residual_scatter=robust_scatter(curve.flux[indices] - model),
        coefficients=coefficients,
        model=model,
        used_indices=indices,
    )


def grid_box_fit(
    curve: LightCurveData,
    period_days: float,
    epoch_bjd: float,
    durations_days: np.ndarray,
    offsets_days: np.ndarray,
    window_days: float = 0.22,
    allowed_cycles: set[int] | None = None,
) -> tuple[BoxFit, np.ndarray, dict[str, float]]:
    """Search a small duration/time grid and return joint 68% profile bounds."""
    durations = np.asarray(durations_days, dtype=float)
    offsets = np.asarray(offsets_days, dtype=float)
    chi2_grid = np.full((len(durations), len(offsets)), np.nan)
    fits: dict[tuple[int, int], BoxFit] = {}
    for duration_index, duration in enumerate(durations):
        for offset_index, offset in enumerate(offsets):
            fit = fit_repeating_box(
                curve,
                period_days=period_days,
                epoch_bjd=epoch_bjd,
                duration_days=float(duration),
                time_offset_days=float(offset),
                window_days=window_days,
                allowed_cycles=allowed_cycles,
            )
            chi2_grid[duration_index, offset_index] = fit.chi2
            fits[(duration_index, offset_index)] = fit

    best_flat = int(np.nanargmin(chi2_grid))
    best_index = np.unravel_index(best_flat, chi2_grid.shape)
    best = fits[best_index]
    scale = max(1.0, best.reduced_chi2)
    joint_region = (chi2_grid - best.chi2) / scale <= 2.30
    duration_values = np.broadcast_to(durations[:, None], chi2_grid.shape)
    offset_values = np.broadcast_to(offsets[None, :], chi2_grid.shape)
    bounds = {
        "duration_low_days": float(np.min(duration_values[joint_region])),
        "duration_high_days": float(np.max(duration_values[joint_region])),
        "offset_low_days": float(np.min(offset_values[joint_region])),
        "offset_high_days": float(np.max(offset_values[joint_region])),
    }
    return best, chi2_grid, bounds


def fit_one_event(
    curve: LightCurveData,
    period_days: float,
    epoch_bjd: float,
    cycle: int,
    duration_days: float,
    offsets_days: np.ndarray,
    window_days: float = 0.20,
) -> dict[str, float | int | bool]:
    """Measure one event time using a fixed duration and a local line."""
    predicted = epoch_bjd + cycle * period_days
    dt = curve.time_bjd - predicted
    keep = curve.good & (np.abs(dt) <= window_days)
    indices = np.flatnonzero(keep)
    if len(indices) < 6:
        raise RuntimeError(f"Cycle {cycle} has too few nearby points")
    full_window = bool(
        np.min(dt[indices]) <= -0.75 * window_days
        and np.max(dt[indices]) >= 0.75 * window_days
    )
    x = dt[indices] / window_days
    y = curve.flux[indices]
    error = curve.flux_error[indices]
    chi2_values = np.full(len(offsets_days), np.nan)
    saved: list[tuple[np.ndarray, np.ndarray, np.ndarray, float, float]] = []
    for index, offset in enumerate(offsets_days):
        box = integrated_box_fraction(
            dt[indices] - offset,
            duration_days=duration_days,
            cadence_days=curve.cadence_days,
        )
        design = np.column_stack([np.ones(len(indices)), x, -box])
        result = _solve_weighted_model(design, y, error)
        saved.append(result)
        chi2_values[index] = result[3]
    best_index = int(np.nanargmin(chi2_values))
    coefficients, covariance, model, chi2, reduced_chi2 = saved[best_index]
    scale = max(1.0, reduced_chi2)
    interval = offsets_days[(chi2_values - chi2) / scale <= 1.0]
    timing_error = (
        0.5 * float(np.max(interval) - np.min(interval))
        if len(interval) >= 2
        else float(np.median(np.diff(offsets_days)))
    )
    depth = float(coefficients[-1])
    depth_error = float(np.sqrt(max(0.0, covariance[-1, -1])))
    measured_time = predicted + float(offsets_days[best_index])
    return {
        "sector": curve.sector,
        "pipeline": curve.pipeline,
        "cycle": cycle,
        "predicted_bjd": predicted,
        "measured_bjd": measured_time,
        "time_error_days": timing_error,
        "observed_minus_calculated_minutes": (
            measured_time - predicted
        )
        * 1440.0,
        "depth_ppt": depth * 1000.0,
        "depth_error_ppt": depth_error * 1000.0,
        "depth_snr": depth / depth_error if depth_error > 0 else float("nan"),
        "points": len(indices),
        "full_local_window": full_window,
        "reduced_chi2": reduced_chi2,
        "residual_scatter_ppt": robust_scatter(y - model) * 1000.0,
    }


def weighted_linear_ephemeris(
    cycles: np.ndarray, times_bjd: np.ndarray, errors_days: np.ndarray
) -> dict[str, float]:
    """Fit T(E) = T0 + P E with a centered weighted straight line."""
    cycle = np.asarray(cycles, dtype=float)
    time = np.asarray(times_bjd, dtype=float)
    error = np.asarray(errors_days, dtype=float)
    good = np.isfinite(cycle) & np.isfinite(time) & np.isfinite(error) & (error > 0)
    cycle, time, error = cycle[good], time[good], error[good]
    if len(cycle) < 3:
        raise RuntimeError("At least three event times are required for an ephemeris")
    weights = 1.0 / np.square(error)
    reference_cycle = int(np.rint(np.average(cycle, weights=weights)))
    centered_cycle = cycle - reference_cycle
    design = np.column_stack([np.ones(len(cycle)), centered_cycle])
    coefficients, covariance, model, chi2, reduced_chi2 = _solve_weighted_model(
        design, time, error
    )
    center_time, period = coefficients
    transform = np.array([[1.0, -reference_cycle], [0.0, 1.0]])
    uncentered_covariance = transform @ covariance @ transform.T
    epoch_zero = center_time - reference_cycle * period
    return {
        "epoch_bjd": float(epoch_zero),
        "epoch_error_days": float(np.sqrt(uncentered_covariance[0, 0])),
        "period_days": float(period),
        "period_error_days": float(np.sqrt(uncentered_covariance[1, 1])),
        "epoch_period_covariance_days2": float(uncentered_covariance[0, 1]),
        "reference_cycle": reference_cycle,
        "reference_time_bjd": float(center_time),
        "reference_time_error_days": float(np.sqrt(covariance[0, 0])),
        "events": len(cycle),
        "chi2": chi2,
        "reduced_chi2": reduced_chi2,
        "residual_rms_minutes": float(np.sqrt(np.mean(np.square(time - model))) * 1440.0),
    }


def bls_period_search(
    curve: LightCurveData,
    minimum_period_days: float = 1.0,
    maximum_period_days: float = 5.0,
) -> dict[str, np.ndarray | float]:
    """Run an independent box-least-squares search on one sector."""
    good = curve.good
    time = curve.time_bjd[good]
    flux = curve.flux[good]
    error = curve.flux_error[good]
    model = BoxLeastSquares(time, flux, dy=error)
    durations = np.array([0.055, 0.075, 0.095, 0.115])
    result = model.autopower(
        durations,
        minimum_period=minimum_period_days,
        maximum_period=maximum_period_days,
        frequency_factor=2.0,
        objective="snr",
    )
    best = int(np.nanargmax(result.power))
    return {
        "period_grid_days": np.asarray(result.period, dtype=float),
        "power": np.asarray(result.power, dtype=float),
        "best_period_days": float(result.period[best]),
        "best_duration_days": float(result.duration[best]),
        "best_time_bjd": float(result.transit_time[best]),
        "best_depth": float(result.depth[best]),
        "best_depth_error": float(result.depth_err[best]),
        "best_power": float(result.power[best]),
    }


def independent_peaks(
    periods: np.ndarray,
    power: np.ndarray,
    count: int = 5,
    separation_days: float = 0.03,
) -> list[tuple[float, float]]:
    """Return the strongest separated periods from a periodogram."""
    order = np.argsort(np.asarray(power))[::-1]
    selected: list[tuple[float, float]] = []
    for index in order:
        period = float(periods[index])
        if all(abs(period - old_period) >= separation_days for old_period, _ in selected):
            selected.append((period, float(power[index])))
        if len(selected) == count:
            break
    return selected
