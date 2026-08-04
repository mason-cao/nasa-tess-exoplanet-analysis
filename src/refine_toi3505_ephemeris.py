"""Refit the TOI-3505.01 linear ephemeris using the best pipeline per sector.

The frozen four-sector analysis in ``analyze_toi3505_tess.py`` deliberately uses
QLP everywhere so that all four sectors are treated identically.  That is the
right choice for comparing sectors against each other, but it throws away real
timing precision: SPOC published two-minute light curves for Sectors 54 and 81,
and a two-minute sample constrains a mid-transit time far better than a
ten-minute Full Frame Image sample of the same transit.

This script therefore builds one more event list on a different rule:

* Sectors 14 and 41 keep QLP, because no SPOC light curve exists for them.
* Sectors 54 and 81 switch to the SPOC two-minute PDCSAP series.

QLP and SPOC measure the *same photons* for Sectors 54 and 81, so their event
times are not independent and must never both enter one fit.  Choosing the
better-sampled pipeline per sector keeps every transit represented exactly once.

Everything else matches the frozen pipeline except the event shape: the same
local linear baseline, the same core acceptance rules, and the same
square-root-of-reduced-chi-square error inflation are retained, while an
exposure-integrated trapezoid represents the near-grazing transit. Four fits
are reported so the pipeline and shape choices remain auditable:

1. ``qlp_only`` reproduces the frozen four-sector result.
2. ``best_per_sector`` is the adopted fit.
3. ``spoc_only`` uses Sectors 54 and 81 alone, as a consistency check.
4. ``box_shape_control`` repeats the adopted pipeline choice with a box.
"""

from __future__ import annotations

import json
import os
import tempfile
from dataclasses import dataclass
from pathlib import Path

os.environ.setdefault(
    "MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "nasa-tess-matplotlib-cache")
)

import numpy as np
import pandas as pd

from toi3505_tess_tools import (
    LightCurveData,
    fit_one_event,
    flat_fraction_from_geometry,
    grid_box_fit,
    load_light_curve,
)

ROOT = Path(__file__).resolve().parents[1]
LIGHT_CURVES = ROOT / "data" / "tess" / "toi3505" / "light_curves" / "mastDownload"
OUTPUT_DIR = ROOT / "outputs" / "toi3505_ephemeris_refined"

# Catalog ephemeris that defines cycle zero.  Every cycle number in this file
# and in the poster is counted from this epoch so the numbers stay comparable
# with the ExoFOP row.
PERIOD_DAYS = 2.9151556
PERIOD_ERROR_DAYS = 0.0000117
EPOCH_BJD = 2459793.534385
EPOCH_ERROR_DAYS = 0.0020787

# Official SPOC multi-sector fit, for the honest side-by-side comparison, and
# the geometry that sets the transit shape.  TOI-3505.01 is nearly grazing, so
# ingress and egress take up most of the event and a box is a poor stand-in.
SPOC_PERIOD_DAYS = 2.915145579641331
SPOC_PERIOD_ERROR_DAYS = 6.9513185e-06
SPOC_RADIUS_RATIO = 0.061769367661928656
SPOC_IMPACT_PARAMETER = 0.9159642085456406
FLAT_FRACTION = flat_fraction_from_geometry(SPOC_RADIUS_RATIO, SPOC_IMPACT_PARAMETER)

SECTORS = (14, 41, 54, 81)
FIT_WINDOW_DAYS = 0.22
EVENT_WINDOW_DAYS = 0.20

# Acceptance rules.  The first two are copied from the frozen pipeline so the
# event lists are selected the same way.
MINIMUM_DEPTH_SNR = 2.5
MAXIMUM_TIME_ERROR_DAYS = 0.03

# The frozen pipeline only asked for points near both ends of the window, which
# lets an event through even when the middle is missing.  Sector 81 cycle 247
# sits on the edge of a downlink gap: it keeps points at both ends but covers
# less than half the window, the local baseline and the transit depth stop being
# separable, and the fit returns a 24 ppt "transit" where the real one is 3 ppt.
# Requiring most of the window to be filled removes it.  The accepted events
# cover 75 to 100 per cent and the two rejected ones cover 34 and 46 per cent,
# so any threshold between those two groups gives the same answer.
MINIMUM_WINDOW_COVERAGE = 0.70

# Forward-propagation dates, as barycentric Julian dates at 00:00 UT.
FORWARD_DATES = {
    "2026-08-01": 2461253.5,
    "2027-08-01": 2461618.5,
    "2030-08-01": 2462714.5,
}


@dataclass(frozen=True)
class Ephemeris:
    label: str
    events: int
    sectors: tuple[int, ...]
    epoch_bjd: float
    epoch_error_days: float
    period_days: float
    period_error_days: float
    covariance_days2: float
    reduced_chi2: float
    residual_rms_minutes: float
    error_scale: float

    def as_dict(self) -> dict[str, object]:
        return {
            "label": self.label,
            "events": self.events,
            "sectors": list(self.sectors),
            "epoch_bjd_tdb": self.epoch_bjd,
            "epoch_error_days": self.epoch_error_days,
            "period_days": self.period_days,
            "period_error_days": self.period_error_days,
            "covariance_days2": self.covariance_days2,
            "reduced_chi2": self.reduced_chi2,
            "residual_rms_minutes": self.residual_rms_minutes,
            "error_scale": self.error_scale,
        }


def curve_paths() -> dict[int, Path]:
    """Best available light curve per sector: SPOC where it exists, else QLP."""
    chosen: dict[int, Path] = {}
    for sector in SECTORS:
        spoc = sorted(
            (LIGHT_CURVES / "TESS").glob(f"*-s{sector:04d}-*/*_lc.fits")
        )
        if spoc:
            if len(spoc) > 1:
                raise RuntimeError(f"More than one SPOC curve for Sector {sector}")
            chosen[sector] = spoc[0]
            continue
        qlp = sorted(
            (LIGHT_CURVES / "HLSP").glob(f"hlsp_qlp_*s{sector:04d}*/*.fits")
        )
        if len(qlp) != 1:
            raise RuntimeError(f"Expected exactly one QLP curve for Sector {sector}")
        chosen[sector] = qlp[0]
    return chosen


def qlp_paths() -> dict[int, Path]:
    """QLP everywhere, to reproduce the frozen four-sector fit."""
    chosen: dict[int, Path] = {}
    for sector in SECTORS:
        qlp = sorted(
            (LIGHT_CURVES / "HLSP").glob(f"hlsp_qlp_*s{sector:04d}*/*.fits")
        )
        if len(qlp) != 1:
            raise RuntimeError(f"Expected exactly one QLP curve for Sector {sector}")
        chosen[sector] = qlp[0]
    return chosen


def offset_grid(curve: LightCurveData) -> np.ndarray:
    """Timing search grid, refined so the step never limits the measurement.

    The frozen pipeline uses a 43-second step, which is well below the timing
    precision a ten-minute sample can reach but is coarse next to two-minute
    data.  The step is tied to the cadence here so the grid is never the thing
    that sets the error bar.
    """
    step_days = min(0.0005, curve.cadence_days / 8.0)
    half_width = 0.07
    count = int(round(2 * half_width / step_days)) + 1
    return np.linspace(-half_width, half_width, count)


def measure_events(
    curves: dict[int, LightCurveData], flat_fraction: float
) -> pd.DataFrame:
    """Fit one mid-transit time per observed event, sector by sector."""
    durations = np.linspace(1.3 / 24.0, 3.4 / 24.0, 43)
    duration_offsets = np.linspace(-0.035, 0.035, 71)
    rows: list[dict[str, object]] = []
    for sector in SECTORS:
        curve = curves[sector]
        best, _, _ = grid_box_fit(
            curve,
            period_days=PERIOD_DAYS,
            epoch_bjd=EPOCH_BJD,
            durations_days=durations,
            offsets_days=duration_offsets,
            window_days=FIT_WINDOW_DAYS,
            flat_fraction=flat_fraction,
        )
        offsets = offset_grid(curve)
        good_times = curve.time_bjd[curve.good]
        first = int(np.floor((good_times.min() - EPOCH_BJD) / PERIOD_DAYS)) - 1
        last = int(np.ceil((good_times.max() - EPOCH_BJD) / PERIOD_DAYS)) + 1
        for cycle in range(first, last + 1):
            predicted = EPOCH_BJD + cycle * PERIOD_DAYS
            if predicted < good_times.min() - 0.1 or predicted > good_times.max() + 0.1:
                continue
            try:
                row = fit_one_event(
                    curve,
                    PERIOD_DAYS,
                    EPOCH_BJD,
                    cycle,
                    best.duration_days,
                    offsets,
                    window_days=EVENT_WINDOW_DAYS,
                    flat_fraction=flat_fraction,
                )
            except RuntimeError:
                continue
            expected_points = 2.0 * EVENT_WINDOW_DAYS / curve.cadence_days
            coverage = float(row["points"]) / expected_points
            row["cadence_minutes"] = curve.cadence_days * 1440.0
            row["fitted_duration_hours"] = best.duration_days * 24.0
            row["flat_fraction"] = flat_fraction
            row["window_coverage"] = coverage
            row["used_in_ephemeris"] = bool(
                row["full_local_window"]
                and coverage >= MINIMUM_WINDOW_COVERAGE
                and row["depth_snr"] >= MINIMUM_DEPTH_SNR
                and row["time_error_days"] <= MAXIMUM_TIME_ERROR_DAYS
            )
            rows.append(row)
    return pd.DataFrame(rows).sort_values(["measured_bjd", "sector"]).reset_index(
        drop=True
    )


def fit_ephemeris(events: pd.DataFrame, label: str) -> Ephemeris:
    """Weighted straight line through the accepted times, with error inflation.

    The fit is pivoted at the weighted-mean cycle so the epoch and the period
    come out very nearly uncorrelated, then transformed back to catalog cycle
    zero.  The residuals are wider than the formal error bars, so the covariance
    is scaled by the reduced chi-square before anything is quoted.

    The covariance is built here rather than through
    ``toi3505_tess_tools.weighted_linear_ephemeris`` on purpose: that helper
    routes through ``_solve_weighted_model``, which already applies the same
    reduced-chi-square inflation internally, and scaling its output again would
    count the penalty twice.
    """
    selected = events[events["used_in_ephemeris"]]
    cycles = selected["cycle"].to_numpy(dtype=float)
    times = selected["measured_bjd"].to_numpy(dtype=float)
    errors = selected["time_error_days"].to_numpy(dtype=float)
    if len(cycles) < 3:
        raise RuntimeError(f"{label}: need at least three accepted events")

    weight = 1.0 / errors**2
    pivot = float(np.sum(weight * cycles) / np.sum(weight))
    design = np.vstack([np.ones_like(cycles), cycles - pivot]).T
    covariance = np.linalg.inv(design.T @ (design * weight[:, None]))
    best = covariance @ (design.T @ (weight * times))

    residual = times - design @ best
    chi2 = float(np.sum(weight * residual**2))
    reduced_chi2 = chi2 / (len(cycles) - 2)
    scale = float(np.sqrt(reduced_chi2)) if reduced_chi2 > 1.0 else 1.0
    covariance = covariance * scale**2

    period = float(best[1])
    epoch_at_zero = float(best[0]) - pivot * period
    jacobian = np.array([1.0, -pivot])
    return Ephemeris(
        label=label,
        events=int(len(cycles)),
        sectors=tuple(sorted(int(s) for s in selected["sector"].unique())),
        epoch_bjd=epoch_at_zero,
        epoch_error_days=float(np.sqrt(jacobian @ covariance @ jacobian)),
        period_days=period,
        period_error_days=float(np.sqrt(covariance[1, 1])),
        covariance_days2=float(jacobian @ covariance @ np.array([0.0, 1.0])),
        reduced_chi2=reduced_chi2,
        residual_rms_minutes=float(np.sqrt(np.mean(residual**2)) * 1440.0),
        error_scale=scale,
    )


def propagation_minutes(
    epoch_error: float, period_error: float, covariance: float, cycle: float
) -> float:
    """One-sigma mid-transit uncertainty after ``cycle`` orbits."""
    variance = (
        epoch_error**2 + (cycle * period_error) ** 2 + 2.0 * cycle * covariance
    )
    return float(np.sqrt(max(0.0, variance)) * 1440.0)


def forward_table(adopted: Ephemeris) -> dict[str, dict[str, float]]:
    """Compare catalog and refined prediction windows on future dates.

    Both ephemerides are expressed at catalog cycle zero, so the same cycle
    index applies to each.  The catalog row publishes no epoch-period
    covariance, so it is taken as zero, which is the conservative reading.
    """
    table: dict[str, dict[str, float]] = {}
    for name, target_bjd in FORWARD_DATES.items():
        cycle = float(round((target_bjd - EPOCH_BJD) / PERIOD_DAYS))
        catalog = propagation_minutes(
            EPOCH_ERROR_DAYS, PERIOD_ERROR_DAYS, 0.0, cycle
        )
        refined = propagation_minutes(
            adopted.epoch_error_days,
            adopted.period_error_days,
            adopted.covariance_days2,
            cycle,
        )
        table[name] = {
            "cycle": cycle,
            "catalog_uncertainty_minutes": catalog,
            "refined_uncertainty_minutes": refined,
            "improvement_factor": catalog / refined if refined > 0 else float("nan"),
        }
    return table


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    best_curves = {s: load_light_curve(p) for s, p in curve_paths().items()}
    qlp_curves = {s: load_light_curve(p) for s, p in qlp_paths().items()}

    best_events = measure_events(best_curves, FLAT_FRACTION)
    qlp_events = measure_events(qlp_curves, FLAT_FRACTION)
    box_events = measure_events(best_curves, 1.0)
    best_events.to_csv(OUTPUT_DIR / "event_times_best_per_sector.csv", index=False)
    qlp_events.to_csv(OUTPUT_DIR / "event_times_qlp_only.csv", index=False)

    adopted = fit_ephemeris(best_events, "best pipeline per sector, trapezoid")
    qlp_only = fit_ephemeris(qlp_events, "QLP in all four sectors, trapezoid")
    box_shape = fit_ephemeris(box_events, "best pipeline per sector, box")
    spoc_subset = best_events[best_events["sector"].isin((54, 81))]
    spoc_only = fit_ephemeris(spoc_subset, "SPOC two-minute, Sectors 54 and 81")

    accepted = best_events[best_events["used_in_ephemeris"]]
    baseline_days = float(
        accepted["measured_bjd"].max() - accepted["measured_bjd"].min()
    )

    summary = {
        "target": "TOI-3505.01",
        "pipeline_choice": {
            sector: {
                "pipeline": best_curves[sector].pipeline,
                "cadence_minutes": best_curves[sector].cadence_days * 1440.0,
                "file": best_curves[sector].path.name,
                "crowdsap": best_curves[sector].crowdsap,
            }
            for sector in SECTORS
        },
        "transit_shape": {
            "model": "exposure-integrated symmetric trapezoid",
            "flat_fraction_t23_over_t14": FLAT_FRACTION,
            "from": "SPOC multi-sector radius ratio and impact parameter",
            "note": (
                "The event is close to grazing, so ingress and egress take up "
                "most of it. A box has infinitely sharp edges and reports "
                "mid-transit times that look more precise than the data allow."
            ),
        },
        "ephemeris": {
            "adopted": adopted.as_dict(),
            "qlp_only": qlp_only.as_dict(),
            "spoc_only": spoc_only.as_dict(),
            "box_shape_control": box_shape.as_dict(),
        },
        "baseline_days": baseline_days,
        "baseline_years": baseline_days / 365.25,
        "comparisons": {
            "catalog_period_days": PERIOD_DAYS,
            "catalog_period_error_days": PERIOD_ERROR_DAYS,
            "spoc_period_days": SPOC_PERIOD_DAYS,
            "spoc_period_error_days": SPOC_PERIOD_ERROR_DAYS,
            "precision_gain_over_catalog": PERIOD_ERROR_DAYS / adopted.period_error_days,
            "precision_gain_over_spoc": (
                SPOC_PERIOD_ERROR_DAYS / adopted.period_error_days
            ),
            "sigma_from_catalog": abs(adopted.period_days - PERIOD_DAYS)
            / float(np.hypot(adopted.period_error_days, PERIOD_ERROR_DAYS)),
            "sigma_from_spoc": abs(adopted.period_days - SPOC_PERIOD_DAYS)
            / float(np.hypot(adopted.period_error_days, SPOC_PERIOD_ERROR_DAYS)),
            "sigma_adopted_from_qlp_only": abs(
                adopted.period_days - qlp_only.period_days
            )
            / float(
                np.hypot(adopted.period_error_days, qlp_only.period_error_days)
            ),
        },
        "forward_propagation": forward_table(adopted),
    }
    (OUTPUT_DIR / "ephemeris_refined.json").write_text(json.dumps(summary, indent=2))

    print(json.dumps(summary["pipeline_choice"], indent=2))
    print()
    for key, fit in summary["ephemeris"].items():
        print(
            f"{key:12s} P = {fit['period_days']:.7f} +/- {fit['period_error_days']:.2e} d"
            f"   n = {fit['events']:2d}   chi2/dof = {fit['reduced_chi2']:.2f}"
            f"   RMS = {fit['residual_rms_minutes']:.1f} min"
        )
    print()
    print(f"Baseline {summary['baseline_years']:.2f} yr")
    print(
        f"Gain over catalog {summary['comparisons']['precision_gain_over_catalog']:.2f}x, "
        f"over SPOC {summary['comparisons']['precision_gain_over_spoc']:.2f}x"
    )
    print(f"\nWrote {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
