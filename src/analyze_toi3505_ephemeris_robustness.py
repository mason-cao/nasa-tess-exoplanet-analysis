"""Quantify robustness of the adopted TOI-3505.01 linear ephemeris.

The canonical fit uses 27 measured TESS mid-transit times and their event-level
uncertainties. This script does not remeasure the light curves. It reads the
accepted event table and asks how the period changes under sector aggregation,
delete-one-sector and delete-one-event resampling, stricter event selection, a
quadratic timing model, and one public MuSCAT2 timing point. The TESS-only
linear fit remains primary; every other result is an explicitly named control.
"""

from __future__ import annotations

import csv
import json
import os
import tempfile
from dataclasses import asdict, dataclass
from pathlib import Path

os.environ.setdefault(
    "MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "nasa-tess-matplotlib-cache")
)

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
EVENT_PATH = (
    ROOT
    / "outputs"
    / "toi3505_ephemeris_refined"
    / "event_times_best_per_sector.csv"
)
EPHEMERIS_PATH = (
    ROOT / "outputs" / "toi3505_ephemeris_refined" / "ephemeris_refined.json"
)
EXOFOP_PATH = ROOT / "data" / "catalogs" / "toi3505" / "exofop_ground_followup.json"
OUTPUT_DIR = ROOT / "outputs" / "toi3505_ephemeris_robustness"

SECTORS = (14, 41, 54, 81)
CATALOG_EPOCH_BJD = 2459793.534385
CATALOG_PERIOD_DAYS = 2.9151556
CATALOG_PERIOD_ERROR_DAYS = 0.0000117
SPOC_PERIOD_DAYS = 2.915145579641331
SPOC_PERIOD_ERROR_DAYS = 6.9513185e-06


@dataclass(frozen=True)
class LinearFit:
    """One weighted linear ephemeris and its error treatment."""

    label: str
    events: int
    epoch_bjd_tdb: float
    epoch_error_days: float
    period_days: float
    period_error_days: float
    covariance_days2: float
    chi_square: float
    degrees_of_freedom: int
    reduced_chi2: float
    residual_rms_minutes: float
    error_scale: float

    def as_dict(self) -> dict[str, object]:
        return asdict(self)


def load_json(path: Path) -> dict[str, object]:
    """Read and validate a JSON object."""
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError(f"Expected a JSON object in {path}")
    return payload


def load_events(path: Path = EVENT_PATH, *, accepted_only: bool = True) -> list[dict[str, object]]:
    """Load the canonical event table without changing its acceptance decisions."""
    events: list[dict[str, object]] = []
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            accepted = row["used_in_ephemeris"] == "True"
            if accepted_only and not accepted:
                continue
            events.append(
                {
                    "sector": int(row["sector"]),
                    "pipeline": row["pipeline"],
                    "cycle": int(row["cycle"]),
                    "measured_bjd": float(row["measured_bjd"]),
                    "time_error_days": float(row["time_error_days"]),
                    "depth_ppt": float(row["depth_ppt"]),
                    "depth_snr": float(row["depth_snr"]),
                    "window_coverage": float(row["window_coverage"]),
                    "used_in_ephemeris": accepted,
                }
            )
    if not events:
        raise RuntimeError(f"No events loaded from {path}")
    return events


def event_arrays(
    events: list[dict[str, object]],
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return cycle, time, and one-sigma uncertainty arrays."""
    cycle = np.array([event["cycle"] for event in events], dtype=float)
    time = np.array([event["measured_bjd"] for event in events], dtype=float)
    error = np.array([event["time_error_days"] for event in events], dtype=float)
    if len(cycle) < 3:
        raise ValueError("At least three timing points are required")
    if np.any(~np.isfinite(cycle)) or np.any(~np.isfinite(time)):
        raise ValueError("Cycle and time arrays must be finite")
    if np.any(~np.isfinite(error)) or np.any(error <= 0.0):
        raise ValueError("Timing uncertainties must be finite and positive")
    if np.unique(cycle).size < 2:
        raise ValueError("At least two distinct cycle values are required")
    return cycle, time, error


def weighted_linear_fit(
    events: list[dict[str, object]], label: str, *, inflate_errors: bool = True
) -> LinearFit:
    """Fit time = epoch + cycle * period with a weighted-cycle pivot.

    The covariance is scaled by sqrt(reduced chi-square) only when the reduced
    chi-square exceeds one, matching the canonical ephemeris implementation.
    """
    cycle, time, error = event_arrays(events)
    weight = 1.0 / error**2
    pivot = float(np.sum(weight * cycle) / np.sum(weight))
    design = np.column_stack((np.ones_like(cycle), cycle - pivot))
    normal = design.T @ (design * weight[:, None])
    covariance = np.linalg.inv(normal)
    coefficients = covariance @ (design.T @ (weight * time))
    residual = time - design @ coefficients
    chi_square = float(np.sum(weight * residual**2))
    dof = len(cycle) - 2
    reduced_chi2 = chi_square / dof
    scale = (
        float(np.sqrt(reduced_chi2))
        if inflate_errors and reduced_chi2 > 1.0
        else 1.0
    )
    covariance = covariance * scale**2

    period = float(coefficients[1])
    epoch = float(coefficients[0] - pivot * period)
    epoch_jacobian = np.array([1.0, -pivot])
    return LinearFit(
        label=label,
        events=len(cycle),
        epoch_bjd_tdb=epoch,
        epoch_error_days=float(
            np.sqrt(epoch_jacobian @ covariance @ epoch_jacobian)
        ),
        period_days=period,
        period_error_days=float(np.sqrt(covariance[1, 1])),
        covariance_days2=float(
            epoch_jacobian @ covariance @ np.array([0.0, 1.0])
        ),
        chi_square=chi_square,
        degrees_of_freedom=dof,
        reduced_chi2=reduced_chi2,
        residual_rms_minutes=float(np.sqrt(np.mean(residual**2)) * 1440.0),
        error_scale=scale,
    )


def quadratic_fit(
    events: list[dict[str, object]], linear: LinearFit
) -> dict[str, object]:
    """Fit a quadratic timing model and compare it to the linear model by BIC."""
    cycle, time, error = event_arrays(events)
    weight = 1.0 / error**2
    pivot = float(np.sum(weight * cycle) / np.sum(weight))
    centered = cycle - pivot
    design = np.column_stack((np.ones_like(cycle), centered, centered**2))
    covariance = np.linalg.inv(design.T @ (design * weight[:, None]))
    coefficients = covariance @ (design.T @ (weight * time))
    residual = time - design @ coefficients
    chi_square = float(np.sum(weight * residual**2))
    dof = len(cycle) - 3
    reduced_chi2 = chi_square / dof
    scale = float(np.sqrt(reduced_chi2)) if reduced_chi2 > 1.0 else 1.0

    time_pivot, period_pivot, quadratic = (float(value) for value in coefficients)
    epoch_zero = time_pivot - period_pivot * pivot + quadratic * pivot**2
    period_zero = period_pivot - 2.0 * quadratic * pivot
    linear_bic = linear.chi_square + 2.0 * np.log(len(cycle))
    quadratic_bic = chi_square + 3.0 * np.log(len(cycle))
    quadratic_error = float(np.sqrt(covariance[2, 2]) * scale)
    return {
        "model": "time = a + b(E-E_pivot) + q(E-E_pivot)^2",
        "events": len(cycle),
        "pivot_cycle": pivot,
        "time_at_pivot_bjd_tdb": time_pivot,
        "period_at_pivot_days": period_pivot,
        "quadratic_term_days_per_cycle2": quadratic,
        "quadratic_term_error_days_per_cycle2": quadratic_error,
        "period_derivative_days_per_cycle": 2.0 * quadratic,
        "epoch_at_cycle_zero_bjd_tdb": epoch_zero,
        "period_at_cycle_zero_days": period_zero,
        "chi_square": chi_square,
        "degrees_of_freedom": dof,
        "reduced_chi2": reduced_chi2,
        "linear_bic": float(linear_bic),
        "quadratic_bic": float(quadratic_bic),
        "delta_chi_square_linear_minus_quadratic": linear.chi_square - chi_square,
        "delta_bic_quadratic_minus_linear": float(quadratic_bic - linear_bic),
        "quadratic_term_significance_sigma": quadratic / quadratic_error,
        "assessment": (
            "The quadratic model is not favored over the linear ephemeris."
            if quadratic_bic >= linear_bic
            else "The quadratic model has a lower BIC than the linear ephemeris."
        ),
    }


def quadratic_prediction(qfit: dict[str, object], cycle: np.ndarray) -> np.ndarray:
    """Evaluate the pivoted quadratic timing model."""
    pivot = float(qfit["pivot_cycle"])
    centered = np.asarray(cycle, dtype=float) - pivot
    return (
        float(qfit["time_at_pivot_bjd_tdb"])
        + float(qfit["period_at_pivot_days"]) * centered
        + float(qfit["quadratic_term_days_per_cycle2"]) * centered**2
    )


def sector_residual_rows(
    events: list[dict[str, object]], linear: LinearFit
) -> list[dict[str, object]]:
    """Summarize timing residuals by observing sector."""
    rows: list[dict[str, object]] = []
    for sector in SECTORS:
        selected = [event for event in events if event["sector"] == sector]
        cycle, time, error = event_arrays(selected)
        residual_minutes = (
            time - (linear.epoch_bjd_tdb + cycle * linear.period_days)
        ) * 1440.0
        weight = 1.0 / error**2
        weighted_mean = float(np.sum(weight * residual_minutes) / np.sum(weight))
        formal_sem = float(1440.0 / np.sqrt(np.sum(weight)))
        lag1 = (
            float(np.corrcoef(residual_minutes[:-1], residual_minutes[1:])[0, 1])
            if len(residual_minutes) > 2
            else float("nan")
        )
        rows.append(
            {
                "sector": sector,
                "events": len(selected),
                "weighted_mean_residual_minutes": weighted_mean,
                "formal_sem_minutes": formal_sem,
                "rms_residual_minutes": float(np.sqrt(np.mean(residual_minutes**2))),
                "lag1_residual_correlation": lag1,
            }
        )
    return rows


def sector_anchor_rows(
    events: list[dict[str, object]], reference_period_days: float
) -> list[dict[str, object]]:
    """Compress each sector to one timing anchor with a conservative error.

    The empirical error is the sample scatter of within-sector residuals about
    the reference period divided by sqrt(N). The adopted anchor error is the
    larger of that value and the formal weighted-mean error.
    """
    rows: list[dict[str, object]] = []
    for sector in SECTORS:
        selected = [event for event in events if event["sector"] == sector]
        cycle, time, error = event_arrays(selected)
        weight = 1.0 / error**2
        anchor_cycle = float(np.sum(weight * cycle) / np.sum(weight))
        reduced_time = time - reference_period_days * (cycle - anchor_cycle)
        anchor_time = float(np.sum(weight * reduced_time) / np.sum(weight))
        formal_error = float(1.0 / np.sqrt(np.sum(weight)))
        empirical_error = float(np.std(reduced_time, ddof=1) / np.sqrt(len(time)))
        rows.append(
            {
                "sector": sector,
                "events": len(selected),
                "anchor_cycle": anchor_cycle,
                "anchor_bjd_tdb": anchor_time,
                "formal_error_days": formal_error,
                "empirical_error_days": empirical_error,
                "adopted_error_days": max(formal_error, empirical_error),
                "error_source": (
                    "empirical within-sector scatter"
                    if empirical_error > formal_error
                    else "formal weighted mean"
                ),
            }
        )
    return rows


def anchor_events(
    rows: list[dict[str, object]], error_field: str
) -> list[dict[str, object]]:
    """Convert sector anchor rows to the common timing-event interface."""
    return [
        {
            "sector": int(row["sector"]),
            "cycle": float(row["anchor_cycle"]),
            "measured_bjd": float(row["anchor_bjd_tdb"]),
            "time_error_days": float(row[error_field]),
        }
        for row in rows
    ]


def jackknife_standard_error(values: list[float]) -> float:
    """Delete-one jackknife standard error for a list of replicate estimates."""
    array = np.asarray(values, dtype=float)
    if array.size < 2 or np.any(~np.isfinite(array)):
        raise ValueError("At least two finite jackknife replicates are required")
    return float(
        np.sqrt((array.size - 1.0) / array.size * np.sum((array - array.mean()) ** 2))
    )


def leave_one_sector_out(
    events: list[dict[str, object]], primary: LinearFit
) -> tuple[list[dict[str, object]], float]:
    """Refit after dropping each sector and return the four-cluster jackknife SE."""
    rows: list[dict[str, object]] = []
    periods: list[float] = []
    for sector in SECTORS:
        selected = [event for event in events if event["sector"] != sector]
        fit = weighted_linear_fit(selected, f"drop Sector {sector}")
        periods.append(fit.period_days)
        rows.append(
            {
                "omitted_sector": sector,
                "events": fit.events,
                "period_days": fit.period_days,
                "period_error_days": fit.period_error_days,
                "period_shift_days": fit.period_days - primary.period_days,
                "period_shift_seconds_per_orbit": (
                    fit.period_days - primary.period_days
                )
                * 86400.0,
                "shift_in_primary_formal_sigma": (
                    fit.period_days - primary.period_days
                )
                / primary.period_error_days,
                "reduced_chi2": fit.reduced_chi2,
            }
        )
    return rows, jackknife_standard_error(periods)


def leave_one_event_out(
    events: list[dict[str, object]], primary: LinearFit
) -> dict[str, object]:
    """Summarize event-level delete-one sensitivity without listing 27 fits."""
    periods = []
    shifts = []
    for index in range(len(events)):
        selected = events[:index] + events[index + 1 :]
        fit = weighted_linear_fit(selected, f"drop event {index}")
        periods.append(fit.period_days)
        shifts.append(fit.period_days - primary.period_days)
    return {
        "replicates": len(periods),
        "jackknife_standard_error_days": jackknife_standard_error(periods),
        "maximum_absolute_period_shift_days": max(abs(value) for value in shifts),
        "maximum_absolute_period_shift_seconds_per_orbit": max(
            abs(value) for value in shifts
        )
        * 86400.0,
    }


def selection_sensitivity(
    events: list[dict[str, object]], primary: LinearFit
) -> list[dict[str, object]]:
    """Repeat the fit at stricter depth-S/N thresholds."""
    rows: list[dict[str, object]] = []
    for threshold in (2.5, 3.0, 4.0):
        selected = [event for event in events if float(event["depth_snr"]) >= threshold]
        fit = weighted_linear_fit(selected, f"depth S/N >= {threshold:g}")
        rows.append(
            {
                "minimum_depth_snr": threshold,
                "events": fit.events,
                "period_days": fit.period_days,
                "period_error_days": fit.period_error_days,
                "period_shift_days": fit.period_days - primary.period_days,
                "period_shift_seconds_per_orbit": (
                    fit.period_days - primary.period_days
                )
                * 86400.0,
                "reduced_chi2": fit.reduced_chi2,
            }
        )
    return rows


def external_muscat2_control(
    events: list[dict[str, object]], exofop: dict[str, object]
) -> tuple[LinearFit, dict[str, object]]:
    """Append the one public ground timing with a reported uncertainty."""
    report_values = exofop["report_values"]
    if not isinstance(report_values, dict):
        raise ValueError("report_values must be an object")
    muscat2 = report_values["muscat2_2023_07_14"]
    if not isinstance(muscat2, dict):
        raise ValueError("MuSCAT2 report values must be an object")
    midpoint = float(muscat2["midpoint_bjd_tdb"])
    error = float(muscat2["midpoint_error_days"])
    cycle = int(round((midpoint - CATALOG_EPOCH_BJD) / CATALOG_PERIOD_DAYS))
    external = {
        "sector": "MuSCAT2",
        "cycle": cycle,
        "measured_bjd": midpoint,
        "time_error_days": error,
    }
    fit = weighted_linear_fit(events + [external], "TESS + MuSCAT2 timing control")
    return fit, {
        "facility": "TCS 1.52 m / MuSCAT2",
        "cycle": cycle,
        "midpoint_bjd_tdb": midpoint,
        "midpoint_error_days": error,
        "source_url": muscat2["url"],
        "report_language": muscat2["report_language"],
        "role": "external control only; not part of the primary TESS-only fit",
    }


def model_comparison_rows(
    canonical: dict[str, object],
    primary: LinearFit,
    conservative: LinearFit,
    selections: list[dict[str, object]],
    external: LinearFit,
) -> list[dict[str, object]]:
    """Assemble the paper-facing period comparison table."""
    ephemerides = canonical["ephemeris"]
    if not isinstance(ephemerides, dict):
        raise ValueError("canonical ephemeris block must be an object")

    def row(
        model: str,
        role: str,
        period: float,
        error: float,
        events: int | None,
        reduced_chi2: float | None,
        note: str,
    ) -> dict[str, object]:
        return {
            "model": model,
            "role": role,
            "events": events,
            "period_days": period,
            "period_error_days": error,
            "difference_from_primary_seconds_per_orbit": (
                period - primary.period_days
            )
            * 86400.0,
            "reduced_chi2": reduced_chi2,
            "note": note,
        }

    qlp = ephemerides["qlp_only"]
    spoc = ephemerides["spoc_only"]
    assert isinstance(qlp, dict) and isinstance(spoc, dict)
    selection_by_threshold = {
        float(item["minimum_depth_snr"]): item for item in selections
    }
    rows = [
        row(
            "Adopted TESS linear",
            "primary",
            primary.period_days,
            primary.period_error_days,
            primary.events,
            primary.reduced_chi2,
            "Best available pipeline per sector; trapezoid event model.",
        ),
        row(
            "Four-sector conservative anchors",
            "robustness",
            conservative.period_days,
            conservative.period_error_days,
            conservative.events,
            conservative.reduced_chi2,
            "One anchor per sector; error is max(formal, empirical sector SEM).",
        ),
        row(
            "Depth S/N >= 3",
            "selection control",
            float(selection_by_threshold[3.0]["period_days"]),
            float(selection_by_threshold[3.0]["period_error_days"]),
            int(selection_by_threshold[3.0]["events"]),
            float(selection_by_threshold[3.0]["reduced_chi2"]),
            "Stricter accepted-event depth threshold.",
        ),
        row(
            "Depth S/N >= 4",
            "selection control",
            float(selection_by_threshold[4.0]["period_days"]),
            float(selection_by_threshold[4.0]["period_error_days"]),
            int(selection_by_threshold[4.0]["events"]),
            float(selection_by_threshold[4.0]["reduced_chi2"]),
            "Stricter accepted-event depth threshold.",
        ),
        row(
            "TESS + MuSCAT2",
            "external control",
            external.period_days,
            external.period_error_days,
            external.events,
            external.reduced_chi2,
            "Adds one public multiband ground midpoint with reported uncertainty.",
        ),
        row(
            "QLP all sectors",
            "pipeline control",
            float(qlp["period_days"]),
            float(qlp["period_error_days"]),
            int(qlp["events"]),
            float(qlp["reduced_chi2"]),
            "Same TESS observations; not an independent detection.",
        ),
        row(
            "SPOC Sectors 54 and 81",
            "pipeline subset control",
            float(spoc["period_days"]),
            float(spoc["period_error_days"]),
            int(spoc["events"]),
            float(spoc["reduced_chi2"]),
            "Two-minute subset only; same TESS observations.",
        ),
        row(
            "ExoFOP TOI catalog",
            "catalog comparison",
            CATALOG_PERIOD_DAYS,
            CATALOG_PERIOD_ERROR_DAYS,
            None,
            None,
            "Public catalog value defining cycle zero.",
        ),
        row(
            "Official SPOC multi-sector",
            "pipeline comparison",
            SPOC_PERIOD_DAYS,
            SPOC_PERIOD_ERROR_DAYS,
            None,
            None,
            "Official SPOC fit; overlaps the same TESS observations.",
        ),
    ]
    return rows


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    """Write a list of consistently shaped dictionaries."""
    if not rows:
        raise ValueError(f"Cannot write empty table {path}")
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def plot_robustness(
    events: list[dict[str, object]],
    primary: LinearFit,
    quadratic: dict[str, object],
    sector_rows: list[dict[str, object]],
    leave_sector_rows: list[dict[str, object]],
    conservative: LinearFit,
    selections: list[dict[str, object]],
    external: LinearFit,
    sector_jackknife_error: float,
    path: Path,
) -> None:
    """Plot timing residuals and the main period-sensitivity estimates."""
    colors = {14: "#7b61a8", 41: "#3478a8", 54: "#d18b2c", 81: "#07845a"}
    ink = "#202a35"
    muted = "#697886"
    fig, (timing, periods) = plt.subplots(
        1, 2, figsize=(12.4, 5.3), dpi=240, gridspec_kw={"width_ratios": [1.45, 1.0]}
    )

    for sector in SECTORS:
        selected = [event for event in events if event["sector"] == sector]
        cycle, time, error = event_arrays(selected)
        residual = (
            time - (primary.epoch_bjd_tdb + cycle * primary.period_days)
        ) * 1440.0
        timing.errorbar(
            cycle,
            residual,
            yerr=error * 1440.0,
            fmt="o",
            markersize=4.2,
            color=colors[sector],
            ecolor=colors[sector],
            elinewidth=0.85,
            capsize=1.5,
            alpha=0.88,
            label=f"Sector {sector}",
            zorder=3,
        )
        summary = next(row for row in sector_rows if row["sector"] == sector)
        timing.errorbar(
            np.mean(cycle),
            float(summary["weighted_mean_residual_minutes"]),
            yerr=float(summary["formal_sem_minutes"]),
            fmt="s",
            markersize=7,
            markeredgecolor="white",
            markeredgewidth=0.8,
            color=colors[sector],
            capsize=3,
            zorder=5,
        )

    grid = np.linspace(
        min(float(event["cycle"]) for event in events),
        max(float(event["cycle"]) for event in events),
        500,
    )
    quadratic_residual = (
        quadratic_prediction(quadratic, grid)
        - (primary.epoch_bjd_tdb + grid * primary.period_days)
    ) * 1440.0
    timing.plot(
        grid,
        quadratic_residual,
        color=muted,
        linestyle="--",
        linewidth=1.5,
        label="quadratic control",
        zorder=2,
    )
    timing.axhline(0.0, color=ink, linewidth=0.9)
    timing.set_xlabel("Transit cycle relative to catalog epoch")
    timing.set_ylabel("Observed minus adopted linear time (minutes)")
    timing.set_title("Twenty-seven TESS timings over 5.05 years")
    timing.grid(color="#dce2e7", linewidth=0.7)
    timing.set_axisbelow(True)
    timing.legend(ncol=3, fontsize=8.3, frameon=False, loc="upper right")

    estimates: list[tuple[str, float, float, str]] = [
        ("Primary", primary.period_days, primary.period_error_days, ink),
        (
            "4-sector anchors",
            conservative.period_days,
            conservative.period_error_days,
            "#006633",
        ),
        (
            "S/N >= 3",
            float(selections[1]["period_days"]),
            float(selections[1]["period_error_days"]),
            "#3478a8",
        ),
        (
            "S/N >= 4",
            float(selections[2]["period_days"]),
            float(selections[2]["period_error_days"]),
            "#3478a8",
        ),
        (
            "+ MuSCAT2",
            external.period_days,
            external.period_error_days,
            "#d18b2c",
        ),
    ]
    for row in leave_sector_rows:
        estimates.append(
            (
                f"drop S{int(row['omitted_sector'])}",
                float(row["period_days"]),
                float(row["period_error_days"]),
                colors[int(row["omitted_sector"])],
            )
        )

    y = np.arange(len(estimates))[::-1]
    shifts = np.array(
        [(period - primary.period_days) * 86400.0 for _, period, _, _ in estimates]
    )
    errors = np.array([error * 86400.0 for _, _, error, _ in estimates])
    for y_value, shift, error, (_, _, _, color) in zip(
        y, shifts, errors, estimates, strict=True
    ):
        periods.errorbar(
            shift,
            y_value,
            xerr=error,
            fmt="o",
            color=color,
            ecolor=color,
            markersize=5.5,
            capsize=2.5,
            zorder=3,
        )
    periods.axvspan(
        -primary.period_error_days * 86400.0,
        primary.period_error_days * 86400.0,
        color="#cfd6dd",
        alpha=0.45,
        label="primary formal 1 sigma",
    )
    periods.axvspan(
        -sector_jackknife_error * 86400.0,
        sector_jackknife_error * 86400.0,
        color="#ffcc33",
        alpha=0.13,
        label="4-sector jackknife scale",
    )
    periods.axvline(0.0, color=ink, linewidth=0.9)
    periods.set_yticks(y, [item[0] for item in estimates])
    periods.set_xlabel("Period difference from primary (seconds per orbit)")
    periods.set_title("Period is stable across controls")
    periods.grid(axis="x", color="#dce2e7", linewidth=0.7)
    periods.set_axisbelow(True)
    periods.legend(loc="upper right", fontsize=8.2, frameon=False)

    for axis in (timing, periods):
        axis.spines[["top", "right"]].set_visible(False)
        axis.tick_params(colors=ink)
    fig.tight_layout()
    fig.savefig(path, facecolor="white", bbox_inches="tight")
    fig.savefig(path.with_suffix(".svg"), facecolor="white", bbox_inches="tight")
    plt.close(fig)


def write_readme(summary: dict[str, object], path: Path) -> None:
    """Document which uncertainty is suitable for which scientific statement."""
    primary = summary["primary_tess_linear"]
    conservative = summary["four_sector_anchor_fit"]
    jackknife = summary["delete_one_sector_jackknife"]
    quadratic = summary["quadratic_model_control"]
    external = summary["external_muscat2_control"]
    assert all(
        isinstance(item, dict)
        for item in (primary, conservative, jackknife, quadratic, external)
    )
    path.write_text(
        f"""# TOI-3505.01 ephemeris robustness

## Primary result

The primary result remains the TESS-only, 27-event linear fit:

`P = {float(primary['period_days']):.10f} +/- {float(primary['period_error_days']):.10f} days`.

The quoted uncertainty is the event-level weighted-fit uncertainty after one
sqrt(reduced chi-square) inflation. It is not replaced by any control below.

## Robustness controls

- Four conservative sector anchors give
  `{float(conservative['period_days']):.10f} +/- {float(conservative['period_error_days']):.10f}` days.
  Each sector error is the larger of its formal weighted-mean error and the
  empirical within-sector standard error.
- Deleting one sector at a time gives a four-cluster jackknife scale of
  `{float(jackknife['jackknife_standard_error_days']):.10f}` days. With only
  four clusters this is a sensitivity diagnostic, not a calibrated replacement
  for the primary uncertainty.
- A quadratic timing model changes chi-square by only
  `{float(quadratic['delta_chi_square_linear_minus_quadratic']):.2f}` and has
  `Delta BIC = {float(quadratic['delta_bic_quadratic_minus_linear']):+.2f}`
  relative to the linear model, so curvature is not favored.
- Adding the public MuSCAT2 timing only as an external control gives
  `P = {float(external['fit']['period_days']):.10f} +/- {float(external['fit']['period_error_days']):.10f}` days.

The stricter depth-S/N selections, delete-one-event fits, pipeline controls,
and model table are recorded in the accompanying JSON and CSV files.

## Reporting rule

Use the primary TESS-only period in the abstract and headline result. Report
the conservative four-sector fit and delete-one-sector scale as robustness
checks. Do not describe MuSCAT2, QLP, and SPOC as mutually independent
detections: MuSCAT2 is external but tentative, while QLP and SPOC reuse TESS
observations.

Regenerate with:

```bash
.venv/bin/python src/analyze_toi3505_ephemeris_robustness.py
```
""",
        encoding="utf-8",
    )


def build_analysis() -> tuple[
    dict[str, object],
    dict[str, list[dict[str, object]]],
    list[dict[str, object]],
]:
    """Compute every robustness product without writing files."""
    events = load_events()
    canonical = load_json(EPHEMERIS_PATH)
    exofop = load_json(EXOFOP_PATH)
    primary = weighted_linear_fit(events, "adopted TESS linear ephemeris")
    quadratic = quadratic_fit(events, primary)
    sector_residuals = sector_residual_rows(events, primary)
    anchors = sector_anchor_rows(events, primary.period_days)
    conservative = weighted_linear_fit(
        anchor_events(anchors, "adopted_error_days"),
        "four conservative sector anchors",
    )
    formal_anchors = weighted_linear_fit(
        anchor_events(anchors, "formal_error_days"),
        "four formal sector anchors",
    )
    leave_sector, sector_jackknife_error = leave_one_sector_out(events, primary)
    leave_event = leave_one_event_out(events, primary)
    selections = selection_sensitivity(events, primary)
    external_fit, external_point = external_muscat2_control(events, exofop)
    models = model_comparison_rows(
        canonical, primary, conservative, selections, external_fit
    )

    summary = {
        "target": "TOI-3505.01",
        "analysis_scope": (
            "Controls derived from the accepted TESS event table plus one named "
            "public MuSCAT2 external timing; no light curves are remeasured here."
        ),
        "primary_tess_linear": primary.as_dict(),
        "sector_residuals": sector_residuals,
        "four_sector_anchors": anchors,
        "four_sector_anchor_fit": conservative.as_dict(),
        "four_sector_formal_anchor_fit": formal_anchors.as_dict(),
        "delete_one_sector_jackknife": {
            "clusters": len(SECTORS),
            "jackknife_standard_error_days": sector_jackknife_error,
            "ratio_to_primary_formal_error": (
                sector_jackknife_error / primary.period_error_days
            ),
            "replicates": leave_sector,
            "interpretation": (
                "Sensitivity diagnostic only: four sectors are too few for a "
                "well-calibrated cluster-resampling uncertainty."
            ),
        },
        "delete_one_event_jackknife": leave_event,
        "selection_sensitivity": selections,
        "quadratic_model_control": quadratic,
        "external_muscat2_control": {
            "point": external_point,
            "fit": external_fit.as_dict(),
            "period_precision_change_percent": (
                100.0
                * (primary.period_error_days - external_fit.period_error_days)
                / primary.period_error_days
            ),
        },
        "reporting_decision": {
            "headline": "primary TESS-only linear fit",
            "robustness": (
                "Report the conservative sector-anchor fit and four-sector "
                "delete-one scale separately."
            ),
            "curvature": "Do not claim a period derivative or transit-timing variation.",
            "external": "MuSCAT2 remains a named external control, not primary input.",
        },
    }
    tables = {
        "sector_residuals.csv": sector_residuals,
        "sector_anchors.csv": anchors,
        "leave_one_sector_out.csv": leave_sector,
        "selection_sensitivity.csv": selections,
        "model_comparison.csv": models,
    }
    return summary, tables, events


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    summary, tables, events = build_analysis()
    (OUTPUT_DIR / "ephemeris_robustness.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    for name, rows in tables.items():
        write_csv(OUTPUT_DIR / name, rows)

    primary = LinearFit(**summary["primary_tess_linear"])  # type: ignore[arg-type]
    conservative = LinearFit(**summary["four_sector_anchor_fit"])  # type: ignore[arg-type]
    external_block = summary["external_muscat2_control"]
    jackknife = summary["delete_one_sector_jackknife"]
    assert isinstance(external_block, dict) and isinstance(jackknife, dict)
    external = LinearFit(**external_block["fit"])  # type: ignore[arg-type]
    plot_robustness(
        events,
        primary,
        summary["quadratic_model_control"],  # type: ignore[arg-type]
        tables["sector_residuals.csv"],
        tables["leave_one_sector_out.csv"],
        conservative,
        tables["selection_sensitivity.csv"],
        external,
        float(jackknife["jackknife_standard_error_days"]),
        OUTPUT_DIR / "01_ephemeris_robustness.png",
    )
    write_readme(summary, OUTPUT_DIR / "README.md")

    quadratic = summary["quadratic_model_control"]
    assert isinstance(quadratic, dict)
    print("TOI-3505.01 ephemeris robustness")
    print(
        f"  primary P = {primary.period_days:.10f} +/- "
        f"{primary.period_error_days:.2e} d"
    )
    print(
        f"  sector anchors P = {conservative.period_days:.10f} +/- "
        f"{conservative.period_error_days:.2e} d"
    )
    print(
        "  delete-one-sector jackknife SE = "
        f"{float(jackknife['jackknife_standard_error_days']):.2e} d"
    )
    print(
        "  quadratic Delta BIC = "
        f"{float(quadratic['delta_bic_quadratic_minus_linear']):+.2f}"
    )
    print(f"Wrote {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
