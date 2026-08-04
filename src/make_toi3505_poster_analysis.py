"""Poster-facing analysis for TOI-3505.01.

This script does not re-measure anything. It reads the frozen measurement
products already written by the pipeline in ``outputs/`` and derives the three
quantities the symposium poster reports that no earlier stage produced:

1. A four-sector linear ephemeris refit with an explicit treatment of the
   timing-residual excess, plus a forward propagation that compares the
   refined ephemeris against the ExoFOP TOI catalog ephemeris.
2. Derived planetary and stellar-system parameters, reported both as the
   uncorrected observed depth and under the documented dilution scenario.
3. Two poster figures that the pipeline never generated: the transit-timing
   observed-minus-calculated diagram and the nearby-star eclipse-depth screen
   in the standard delta-magnitude versus scatter form.

Every input path is relative to the repository root.
"""

from __future__ import annotations

import csv
import json
import math
from dataclasses import dataclass, asdict
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs" / "toi3505_poster"
OUT.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------------
# Catalog constants. Values are the ExoFOP TOI table row for TOI-3505.01 and
# the TIC v8 row for TIC 390988385, both archived under data/catalogs.
# ---------------------------------------------------------------------------

CATALOG = {
    "tic_id": 390988385,
    "toi": "3505.01",
    "tfopwg_disposition": "PC",
    "epoch_bjd_tdb": 2459793.534385,
    "epoch_error_days": 0.0020787,
    "period_days": 2.9151556,
    "period_error_days": 1.17e-05,
    "duration_hours": 2.004,
    "duration_error_hours": 0.21,
    "depth_ppt": 2.910,
    "depth_error_ppt": 0.196266,
    "planet_radius_rearth": 7.4578,
    "insolation_earth": 504.746,
    "tess_mag": 10.9372,
    "sectors_in_toi_row": (14, 41, 54),
}

STAR = {
    "teff_k": 6220.0,
    "radius_rsun": 1.33515,
    "distance_pc": 373.477,
    "tic_contamination_ratio": 0.5471626,
    "v_mag": 11.24,
    # TIC v8 lists no mass for this target. A 6220 K star with a 1.34 solar
    # radius sits on or just above the main sequence, where the empirical
    # Torres relations give about 1.25 solar masses. The value is carried as an
    # explicit assumption and the uncertainty is propagated.
    "assumed_mass_msun": 1.25,
    "assumed_mass_error_msun": 0.15,
}

# Official SPOC multi-sector geometric transit model, scope s0014-s0086, with
# sectors 54 and 81 contributing. Archived in
# outputs/toi3505_data_validation/official_multisector_tce.csv.
SPOC = {
    "period_days": 2.915145579641331,
    "period_error_days": 6.9513185e-06,
    "duration_hours": 2.7111916374936995,
    "duration_error_hours": 0.09769286,
    "depth_ppt": 3.2918889004895218,
    "depth_error_ppt": 0.118516266,
    "radius_ratio": 0.061769367661928656,
    "radius_ratio_error": 0.0015563034,
    "impact_parameter": 0.9159642085456406,
    "impact_parameter_error": 0.013228787,
    "fit_snr": 34.111057,
    "multiple_event_statistic": 27.673754,
    "odd_even_difference_sigma": 0.7784450526530438,
    "weak_secondary_mes": 2.4098403,
    "weak_secondary_depth_ppt": 0.21794022000000002,
    "bootstrap_false_alarm_probability": 1.070414881520924e-143,
    "centroid_offset_arcsec": 2.6222534,
    "centroid_offset_error_arcsec": 2.646006,
    "ghost_core_halo_ratio": 6.905730715165484,
    "observed_transits": 17,
    "all_machine_checks_pass": True,
}

# Physical constants in the unit system used below.
RSUN_REARTH = 109.076
RSUN_RJUP = 9.7311
AU_RSUN = 215.032
TEFF_SUN = 5772.0

# Known close companion from high-resolution imaging, carried through the
# pipeline as a screening quantity only.
COMPANION = {"separation_arcsec": 0.517, "delta_i_mag": 1.7}


# ---------------------------------------------------------------------------
# 1. Four-sector linear ephemeris
# ---------------------------------------------------------------------------


@dataclass
class Ephemeris:
    label: str
    events: int
    epoch_bjd_tdb: float
    epoch_error_days: float
    period_days: float
    period_error_days: float
    covariance_days2: float
    reduced_chi2: float
    residual_rms_minutes: float
    error_scale: float


def load_events() -> list[dict]:
    """Read the adopted event list from the refined ephemeris run.

    This is the trapezoid-shape, best-pipeline-per-sector list produced by
    ``refine_toi3505_ephemeris.py``.  The earlier QLP-and-box list in
    ``outputs/toi3505_tess_analysis`` is kept as the like-for-like comparison
    and is not what the posters quote.
    """
    path = (
        ROOT
        / "outputs"
        / "toi3505_ephemeris_refined"
        / "event_times_best_per_sector.csv"
    )
    with path.open() as handle:
        rows = list(csv.DictReader(handle))
    events = []
    for row in rows:
        if row["used_in_ephemeris"] != "True":
            continue
        events.append(
            {
                "sector": int(row["sector"]),
                "cycle": int(row["cycle"]),
                "pipeline": row["pipeline"],
                "measured_bjd": float(row["measured_bjd"]),
                "time_error_days": float(row["time_error_days"]),
                "depth_ppt": float(row["depth_ppt"]),
                "depth_snr": float(row["depth_snr"]),
            }
        )
    if not events:
        raise RuntimeError(f"No accepted events in {path}")
    return events


def fit_ephemeris(events: list[dict], label: str, scale_errors: bool) -> Ephemeris:
    """Weighted straight-line fit of mid-transit time against cycle number.

    The fit is referenced to the weighted-mean cycle so the reported epoch and
    period are very nearly uncorrelated. The epoch is then transformed back to
    the catalog reference cycle so the two ephemerides can be compared directly.
    """
    cycle = np.array([e["cycle"] for e in events], dtype=float)
    time = np.array([e["measured_bjd"] for e in events], dtype=float)
    sigma = np.array([e["time_error_days"] for e in events], dtype=float)

    weight = 1.0 / sigma**2
    pivot = float(np.sum(weight * cycle) / np.sum(weight))

    design = np.vstack([np.ones_like(cycle), cycle - pivot]).T
    covariance_inverse = design.T @ (design * weight[:, None])
    covariance = np.linalg.inv(covariance_inverse)
    best = covariance @ (design.T @ (weight * time))

    residual = time - design @ best
    chi2 = float(np.sum(weight * residual**2))
    dof = len(events) - 2
    reduced_chi2 = chi2 / dof
    rms_minutes = float(np.sqrt(np.mean(residual**2)) * 24.0 * 60.0)

    scale = math.sqrt(reduced_chi2) if (scale_errors and reduced_chi2 > 1.0) else 1.0
    covariance = covariance * scale**2

    # Transform the epoch from the pivot cycle to catalog cycle zero.
    period = float(best[1])
    epoch_at_pivot = float(best[0])
    epoch_at_zero = epoch_at_pivot - pivot * period
    jacobian = np.array([1.0, -pivot])
    epoch_zero_variance = float(jacobian @ covariance @ jacobian)
    covariance_epoch_period = float(jacobian @ covariance @ np.array([0.0, 1.0]))

    return Ephemeris(
        label=label,
        events=len(events),
        epoch_bjd_tdb=epoch_at_zero,
        epoch_error_days=math.sqrt(epoch_zero_variance),
        period_days=period,
        period_error_days=math.sqrt(float(covariance[1, 1])),
        covariance_days2=covariance_epoch_period,
        reduced_chi2=reduced_chi2,
        residual_rms_minutes=rms_minutes,
        error_scale=scale,
    )


def propagation_uncertainty_minutes(
    epoch_error: float, period_error: float, covariance: float, cycle: float
) -> float:
    """1-sigma mid-transit uncertainty at a given cycle, in minutes."""
    variance = epoch_error**2 + (cycle * period_error) ** 2 + 2.0 * cycle * covariance
    return math.sqrt(max(variance, 0.0)) * 24.0 * 60.0


# ---------------------------------------------------------------------------
# 2. Derived system parameters
# ---------------------------------------------------------------------------


def scaled_separation_from_duration(
    period_days: float, duration_hours: float, radius_ratio: float, impact: float
) -> float:
    """Solve the transit-duration equation for a/R*.

    Using the standard first-to-fourth-contact expression

        T = (P / pi) * arcsin( sqrt((1+k)^2 - b^2) / (a/R* * sin i) ),
        cos i = b / (a/R*),

    this inverts numerically for a/R*. Deriving the scale from the measured
    duration, period, depth, and impact parameter avoids assuming a stellar
    mass, which TIC v8 does not supply for this target.
    """
    duration_days = duration_hours / 24.0
    numerator = math.sqrt(max((1.0 + radius_ratio) ** 2 - impact**2, 1e-12))

    def duration_for(scaled: float) -> float:
        cos_i = impact / scaled
        sin_i = math.sqrt(max(1.0 - cos_i**2, 1e-12))
        argument = numerator / (scaled * sin_i)
        if argument >= 1.0:
            return float("inf")
        return (period_days / math.pi) * math.asin(argument)

    low, high = max(impact + 1e-6, 1.05), 500.0
    for _ in range(200):
        middle = 0.5 * (low + high)
        if duration_for(middle) > duration_days:
            low = middle
        else:
            high = middle
    return 0.5 * (low + high)


def transit_duration_hours(
    period_days: float, scaled: float, radius_ratio: float, impact: float
) -> float:
    """First-to-fourth-contact duration for a circular orbit, in hours."""
    cos_i = impact / scaled
    sin_i = math.sqrt(max(1.0 - cos_i**2, 1e-12))
    numerator = math.sqrt(max((1.0 + radius_ratio) ** 2 - impact**2, 0.0))
    argument = numerator / (scaled * sin_i)
    if argument >= 1.0:
        return float("inf")
    return (period_days / math.pi) * math.asin(argument) * 24.0


def scaled_separation_from_density(
    period_days: float, mass_msun: float, radius_rsun: float
) -> float:
    """a/R* implied by Kepler's third law and the adopted stellar parameters."""
    period_years = period_days / 365.25
    semimajor_au = (mass_msun * period_years**2) ** (1.0 / 3.0)
    return semimajor_au * AU_RSUN / radius_rsun


def impact_from_duration(
    period_days: float, scaled: float, radius_ratio: float, duration_hours: float
) -> float:
    """Impact parameter that reproduces a measured duration at fixed a/R*."""
    low, high = 0.0, 1.0 + radius_ratio
    if transit_duration_hours(period_days, scaled, radius_ratio, low) < duration_hours:
        return float("nan")
    for _ in range(200):
        middle = 0.5 * (low + high)
        if transit_duration_hours(period_days, scaled, radius_ratio, middle) > duration_hours:
            low = middle
        else:
            high = middle
    return 0.5 * (low + high)


def derive_system(
    radius_ratio: float,
    radius_ratio_error: float,
    impact: float,
    duration_hours: float,
    period_days: float,
    label: str,
    note: str,
) -> dict:
    """Planet and orbit parameters from the transit observables."""
    radius_rearth = radius_ratio * STAR["radius_rsun"] * RSUN_REARTH
    radius_rearth_error = radius_ratio_error * STAR["radius_rsun"] * RSUN_REARTH
    radius_rjup = radius_ratio * STAR["radius_rsun"] * RSUN_RJUP
    radius_rjup_error = radius_ratio_error * STAR["radius_rsun"] * RSUN_RJUP

    a_over_rstar = scaled_separation_from_duration(
        period_days, duration_hours, radius_ratio, impact
    )
    inclination_deg = math.degrees(math.acos(impact / a_over_rstar))
    semimajor_au = a_over_rstar * STAR["radius_rsun"] / AU_RSUN

    # Implied stellar mass from Kepler's third law, reported as a consistency
    # check on the derived geometry rather than as an input.
    period_years = period_days / 365.25
    implied_mass_msun = semimajor_au**3 / period_years**2

    equilibrium_k = STAR["teff_k"] * math.sqrt(1.0 / (2.0 * a_over_rstar))
    luminosity_lsun = (STAR["radius_rsun"] ** 2) * (STAR["teff_k"] / TEFF_SUN) ** 4
    insolation = luminosity_lsun / semimajor_au**2

    return {
        "case": label,
        "note": note,
        "radius_ratio": radius_ratio,
        "radius_ratio_error": radius_ratio_error,
        "planet_radius_rearth": radius_rearth,
        "planet_radius_rearth_error": radius_rearth_error,
        "planet_radius_rjup": radius_rjup,
        "planet_radius_rjup_error": radius_rjup_error,
        "impact_parameter": impact,
        "a_over_rstar": a_over_rstar,
        "inclination_degrees": inclination_deg,
        "semimajor_axis_au": semimajor_au,
        "implied_stellar_mass_msun": implied_mass_msun,
        "stellar_luminosity_lsun": luminosity_lsun,
        "equilibrium_temperature_k": equilibrium_k,
        "insolation_earth": insolation,
    }


# ---------------------------------------------------------------------------
# 3. Figures
# ---------------------------------------------------------------------------

PALETTE = {
    "gmu_green": "#1D6F42",
    "gmu_gold": "#FFCC33",
    "ink": "#1B2A22",
    "muted": "#6B7F73",
    "accent": "#B3391E",
    "sector": {14: "#1D6F42", 41: "#2E7FA8", 54: "#B3391E", 81: "#7A5195"},
}


def style_axis(axis) -> None:
    axis.grid(True, color="#DDE5DF", linewidth=0.8, zorder=0)
    axis.set_axisbelow(True)
    for side in ("top", "right"):
        axis.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        axis.spines[side].set_color("#8FA396")
    axis.tick_params(colors=PALETTE["ink"], labelsize=14)
    axis.xaxis.set_minor_locator(AutoMinorLocator())
    axis.yaxis.set_minor_locator(AutoMinorLocator())


def figure_timing(events: list[dict], refined: Ephemeris, path: Path) -> None:
    """Observed-minus-calculated diagram against the catalog ephemeris."""
    cycle = np.array([e["cycle"] for e in events], dtype=float)
    measured = np.array([e["measured_bjd"] for e in events], dtype=float)
    sigma_min = np.array([e["time_error_days"] for e in events]) * 24.0 * 60.0
    sectors = np.array([e["sector"] for e in events])

    catalog_predicted = (
        CATALOG["epoch_bjd_tdb"] + cycle * CATALOG["period_days"]
    )
    o_minus_c = (measured - catalog_predicted) * 24.0 * 60.0

    # Extend the horizon to the 2026 symposium epoch so the forward
    # propagation of the two ephemerides can be compared directly.
    symposium_cycle = round((2461254.0 - CATALOG["epoch_bjd_tdb"]) / CATALOG["period_days"])
    grid = np.linspace(cycle.min() - 40, symposium_cycle + 40, 700)

    figure, axis = plt.subplots(figsize=(13.6, 4.25), dpi=220)
    style_axis(axis)

    # Catalog 1-sigma propagation band.
    catalog_band = np.array(
        [
            propagation_uncertainty_minutes(
                CATALOG["epoch_error_days"], CATALOG["period_error_days"], 0.0, c
            )
            for c in grid
        ]
    )
    axis.fill_between(
        grid,
        -catalog_band,
        catalog_band,
        color=PALETTE["muted"],
        alpha=0.18,
        zorder=1,
    )

    # Refined ephemeris line and band, expressed relative to the catalog.
    refined_offset = (
        refined.epoch_bjd_tdb
        - CATALOG["epoch_bjd_tdb"]
        + grid * (refined.period_days - CATALOG["period_days"])
    ) * 24.0 * 60.0
    refined_band = np.array(
        [
            propagation_uncertainty_minutes(
                refined.epoch_error_days,
                refined.period_error_days,
                refined.covariance_days2,
                c,
            )
            for c in grid
        ]
    )
    axis.fill_between(
        grid,
        refined_offset - refined_band,
        refined_offset + refined_band,
        color=PALETTE["gmu_green"],
        alpha=0.22,
        zorder=2,
    )
    axis.plot(grid, refined_offset, color=PALETTE["gmu_green"], linewidth=2.0, zorder=3)
    axis.axhline(0.0, color=PALETTE["ink"], linewidth=1.0, zorder=3)

    for sector in (14, 41, 54, 81):
        mask = sectors == sector
        if not mask.any():
            continue
        axis.errorbar(
            cycle[mask],
            o_minus_c[mask],
            yerr=sigma_min[mask],
            fmt="o",
            markersize=6,
            color=PALETTE["sector"][sector],
            ecolor=PALETTE["sector"][sector],
            elinewidth=1.4,
            capsize=2.5,
            markeredgecolor="white",
            markeredgewidth=0.7,
            label=f"Sector {sector}",
            zorder=4,
        )

    catalog_at_symposium = propagation_uncertainty_minutes(
        CATALOG["epoch_error_days"], CATALOG["period_error_days"], 0.0, symposium_cycle
    )
    refined_at_symposium = propagation_uncertainty_minutes(
        refined.epoch_error_days,
        refined.period_error_days,
        refined.covariance_days2,
        symposium_cycle,
    )
    axis.axvline(
        symposium_cycle, color=PALETTE["ink"], linewidth=1.0, linestyle=":", zorder=3
    )
    axis.annotate(
        f"August 2026\n{catalog_at_symposium:.0f} min $\\rightarrow$ {refined_at_symposium:.0f} min",
        xy=(symposium_cycle, -34),
        ha="right",
        va="top",
        fontsize=13.5,
        color=PALETTE["ink"],
        xytext=(-12, 0),
        textcoords="offset points",
    )

    axis.set_xlabel("Transit epoch number from the catalog reference time", fontsize=16)
    axis.set_ylabel("Observed $-$ calculated (minutes)", fontsize=16)
    axis.set_title(
        "Transit timing over five years narrows the predicted transit window",
        fontsize=18,
        color=PALETTE["ink"],
        pad=10,
    )
    axis.set_ylim(-46, 46)
    axis.annotate(
        "TESS catalog, 1$\\sigma$",
        xy=(grid[-1], catalog_band[-1]),
        xytext=(-8, 5),
        textcoords="offset points",
        ha="right",
        fontsize=14,
        color=PALETTE["muted"],
    )
    axis.annotate(
        "our fit, 1$\\sigma$",
        xy=(grid[-1], refined_offset[-1] - refined_band[-1]),
        xytext=(-8, -14),
        textcoords="offset points",
        ha="right",
        fontsize=14,
        color=PALETTE["gmu_green"],
        fontweight="bold",
    )
    axis.legend(
        loc="upper left", fontsize=13.5, frameon=True, framealpha=0.96, ncol=4,
        handletextpad=0.3, columnspacing=1.0, borderpad=0.4,
    )

    sectors_used = len({e["sector"] for e in events})
    baseline_years = (
        max(e["measured_bjd"] for e in events)
        - min(e["measured_bjd"] for e in events)
    ) / 365.25
    note = (
        f"{refined.events} events, {sectors_used} sectors, "
        f"{baseline_years:.2f}-year baseline.   Refined period "
        f"{refined.period_days:.7f} $\\pm$ {refined.period_error_days*1e6:.1f}$\\times$10$^{{-6}}$ d "
        f"vs catalog {CATALOG['period_days']:.7f} $\\pm$ {CATALOG['period_error_days']*1e6:.1f}$\\times$10$^{{-6}}$ d"
    )
    figure.text(0.5, 0.012, note, ha="center", fontsize=13, color=PALETTE["muted"])
    figure.tight_layout(rect=(0, 0.05, 1, 1))
    figure.savefig(path, facecolor="white")
    figure.savefig(path.with_suffix(".svg"), facecolor="white")
    plt.close(figure)


def figure_neb(path: Path) -> dict:
    """Nearby-star screen in the delta-magnitude versus scatter form."""
    source = ROOT / "outputs" / "toi3505_ground_checks" / "nearby_star_image_measurements.csv"
    with source.open() as handle:
        rows = list(csv.DictReader(handle))

    delta = np.array([float(r["delta_tmag"]) for r in rows])
    # Depth sensitivity actually reached for each neighbour in the scheduled
    # window. A star is cleared when this is below the eclipse depth it would
    # need in order to produce the observed signal by blending into the target.
    # Two neighbours have too few valid points in the window for a depth fit;
    # they carry no limit and are shown separately as unmeasured.
    sensitivity = 3.0 * np.array(
        [
            float(r["historical_window_depth_error_ppt"])
            if r["historical_window_depth_error_ppt"]
            else np.nan
            for r in rows
        ]
    )
    unmeasured = np.isnan(sensitivity)
    required = np.array([float(r["required_eclipse_depth_ppt_simple"]) for r in rows])
    overlap = np.array([r["target_aperture_overlap"] == "True" for r in rows])
    cleared = np.array([r["transit_relevant_clearance"] == "True" for r in rows])
    separation = np.array([float(r["separation_arcsec"]) for r in rows])

    figure, axis = plt.subplots(figsize=(9.2, 5.4), dpi=220)
    style_axis(axis)
    axis.set_yscale("log")

    order = np.argsort(delta)
    axis.plot(
        delta[order],
        required[order],
        color=PALETTE["accent"],
        linewidth=2.4,
        zorder=5,
    )

    axis.scatter(
        delta[cleared],
        sensitivity[cleared],
        s=62,
        color=PALETTE["gmu_green"],
        edgecolor="white",
        linewidth=0.8,
        zorder=6,
        label=f"ruled out ({int(cleared.sum())})",
    )
    noisy = ~cleared & ~overlap & ~unmeasured
    axis.scatter(
        delta[noisy],
        sensitivity[noisy],
        s=62,
        marker="s",
        color=PALETTE["muted"],
        edgecolor="white",
        linewidth=0.8,
        zorder=6,
        label=f"too noisy ({int(noisy.sum())})",
    )
    if unmeasured.any():
        axis.scatter(
            delta[unmeasured],
            np.full(int(unmeasured.sum()), 4.2),
            s=62,
            marker="v",
            color=PALETTE["ink"],
            edgecolor="white",
            linewidth=0.8,
            zorder=6,
            label=f"no limit ({int(unmeasured.sum())})",
        )
    axis.scatter(
        delta[overlap],
        sensitivity[overlap],
        s=150,
        marker="X",
        color=PALETTE["gmu_gold"],
        edgecolor=PALETTE["ink"],
        linewidth=1.1,
        zorder=7,
        label=f"blended ({int(overlap.sum())})",
    )

    axis.set_xlabel("$\\Delta T$ relative to TOI-3505 (magnitudes)", fontsize=16)
    axis.set_ylabel("3$\\sigma$ eclipse-depth limit reached (ppt)", fontsize=16)
    axis.set_title(
        "Stars near TOI-3505 that could fake the signal",
        fontsize=18,
        color=PALETTE["ink"],
        pad=10,
    )
    axis.set_ylim(3, 6e3)
    axis.legend(
        loc="upper left", fontsize=14, frameon=True, framealpha=0.96, ncol=4,
        handletextpad=0.3, columnspacing=0.9, borderpad=0.4,
    )
    # Sit the curve label on the curve itself, angled to match how the line
    # actually renders on the log axis.
    label_x = 2.4
    label_y = float(np.interp(label_x, delta[order], required[order])) * 1.9
    axis.annotate(
        "eclipse depth needed to fake the signal",
        xy=(label_x, label_y),
        ha="center",
        fontsize=14,
        color=PALETTE["accent"],
        rotation=23,
        rotation_mode="anchor",
    )
    axis.annotate(
        "stars below the line are ruled out",
        xy=(0.985, 0.12),
        xycoords="axes fraction",
        ha="right",
        fontsize=14,
        color=PALETTE["muted"],
        style="italic",
    )
    figure.tight_layout()
    figure.savefig(path, facecolor="white")
    figure.savefig(path.with_suffix(".svg"), facecolor="white")
    plt.close(figure)

    return {
        "stars_measured": int(len(rows)),
        "cleared": int(cleared.sum()),
        "not_cleared_limit_too_shallow": int((~cleared & ~overlap & ~unmeasured).sum()),
        "not_cleared_no_limit": int(unmeasured.sum()),
        "blended_with_target": int(overlap.sum()),
        "closest_unblended_arcsec": float(separation[~overlap].min()),
        "faintest_delta_tmag_cleared": float(delta[cleared].max())
        if cleared.any()
        else float("nan"),
        "clearance_criterion": (
            "3 sigma depth limit in the scheduled window is below the eclipse "
            "depth the neighbour would need, with full window coverage and no "
            "aperture blending with the target."
        ),
    }


def figure_phase_folded(path: Path, refined: "Ephemeris") -> dict:
    """Four TESS sectors phase-folded on the refined ephemeris.

    Replaces the analysis-stage version with a poster aspect ratio and folds on
    the ephemeris measured in this work rather than the catalog one.
    """
    figure, axes = plt.subplots(2, 2, figsize=(12.0, 6.0), dpi=220, sharey=True)
    summary = []

    with (ROOT / "outputs" / "toi3505_tess_analysis" / "sector_measurements.csv").open() as handle:
        measured = {int(r["sector"]): r for r in csv.DictReader(handle)}

    for axis, sector in zip(axes.ravel(), (14, 41, 54, 81)):
        style_axis(axis)
        source = (
            ROOT
            / "outputs"
            / "toi3505_tess_analysis"
            / "clean_light_curves"
            / f"sector_{sector}_qlp.csv"
        )
        with source.open() as handle:
            rows = [r for r in csv.DictReader(handle) if r["quality_zero"] == "True"]

        time = np.array([float(r["time_bjd_tdb"]) for r in rows])
        flux = np.array([float(r["normalized_flux"]) for r in rows])

        cycle = np.round((time - refined.epoch_bjd_tdb) / refined.period_days)
        hours = (time - (refined.epoch_bjd_tdb + cycle * refined.period_days)) * 24.0
        keep = np.abs(hours) < 4.0
        hours, flux = hours[keep], flux[keep]

        axis.scatter(hours, flux, s=5, color="#B9C6BD", edgecolor="none", zorder=2)

        edges = np.arange(-4.0, 4.001, 0.25)
        index = np.digitize(hours, edges) - 1
        centre, mean, error = [], [], []
        for slot in range(len(edges) - 1):
            selected = flux[index == slot]
            if selected.size < 3:
                continue
            centre.append(0.5 * (edges[slot] + edges[slot + 1]))
            mean.append(float(np.mean(selected)))
            error.append(float(np.std(selected) / math.sqrt(selected.size)))
        axis.errorbar(
            centre,
            mean,
            yerr=error,
            fmt="o",
            markersize=4.5,
            color=PALETTE["sector"][sector],
            ecolor=PALETTE["sector"][sector],
            elinewidth=1.2,
            capsize=0,
            markeredgecolor="white",
            markeredgewidth=0.5,
            zorder=4,
        )

        record = measured[sector]
        depth = float(record["depth_ppt"]) * 1e-3
        duration = float(record["duration_hours"])
        axis.plot(
            [-4, -duration / 2, -duration / 2, duration / 2, duration / 2, 4],
            [1, 1, 1 - depth, 1 - depth, 1, 1],
            color=PALETTE["accent"],
            linewidth=1.8,
            zorder=5,
        )
        axis.axhline(1.0, color=PALETTE["ink"], linewidth=0.8, zorder=3)
        axis.set_xlim(-4, 4)
        axis.set_ylim(0.9915, 1.0055)
        axis.set_title(
            f"Sector {sector}   {float(record['depth_ppt']):.2f} $\\pm$ "
            f"{float(record['depth_error_ppt']):.2f} ppt",
            fontsize=16,
            color=PALETTE["ink"],
            pad=5,
        )
        summary.append({"sector": sector, "points": int(keep.sum())})

    for axis in axes[1]:
        axis.set_xlabel("Hours from mid-transit", fontsize=15)
    for axis in axes[:, 0]:
        axis.set_ylabel("Relative flux", fontsize=15)

    figure.suptitle(
        "The same dip shows up in four TESS sectors over five years",
        fontsize=18,
        color=PALETTE["ink"],
        y=0.985,
    )
    figure.tight_layout(rect=(0, 0, 1, 0.96))
    figure.savefig(path, facecolor="white")
    figure.savefig(path.with_suffix(".svg"), facecolor="white")
    plt.close(figure)
    return {"sectors": summary, "folded_on": "refined four-sector ephemeris"}


def figure_ground_light_curve(path: Path) -> dict:
    """Poster version of the final ground light curve.

    The AstroImageJ deliverable in ``outputs/toi3505_review_package`` carries the
    full diagnostic stack and a y-range set by those traces. For the poster the
    same adopted measurements are replotted on a flux range matched to the data,
    with the scheduled window and the injected-signal comparison shown directly.
    """
    source = (
        ROOT
        / "outputs"
        / "toi3505_final_candidate"
        / "TOI_3505.01_2022-07-22_R_final_candidate_light_curve.csv"
    )
    with source.open() as handle:
        rows = list(csv.DictReader(handle))

    used = [r for r in rows if r["used_in_primary_curve"] == "True"]
    dropped = [r for r in rows if r["used_in_primary_curve"] != "True"]

    def column(records, name):
        return np.array([float(r[name]) for r in records])

    time = column(used, "bjd_tdb") - 2459782.0
    flux = column(used, "adopted_relative_brightness")
    drop_time = column(dropped, "bjd_tdb") - 2459782.0
    drop_flux = column(dropped, "adopted_relative_brightness")

    bins = (
        ROOT
        / "outputs"
        / "toi3505_final_candidate"
        / "TOI_3505.01_2022-07-22_R_final_candidate_10min_bins.csv"
    )
    with bins.open() as handle:
        binned = list(csv.DictReader(handle))
    first = float(used[0]["bjd_tdb"]) - 2459782.0
    bin_time = column(binned, "hours") / 24.0 + first
    bin_flux = column(binned, "relative_brightness")
    bin_error = column(binned, "uncertainty")

    check = json.loads(
        (ROOT / "outputs" / "toi3505_final_candidate" / "historical_schedule_check.json").read_text()
    )
    window = check["working_interpretation"]["times"]
    ingress = window["ingress"]["bjd_tdb"] - 2459782.0
    egress = window["egress"]["bjd_tdb"] - 2459782.0
    fixed = check["fixed_window_check"]

    figure, axis = plt.subplots(figsize=(11.0, 5.5), dpi=220)
    style_axis(axis)

    axis.axvspan(
        ingress,
        egress,
        color=PALETTE["accent"],
        alpha=0.10,
        zorder=1,
    )
    for edge in (ingress, egress):
        axis.axvline(edge, color=PALETTE["accent"], linewidth=1.4, linestyle="--", zorder=2)
    axis.axhline(1.0, color=PALETTE["ink"], linewidth=1.0, zorder=2)

    # The transit that would have been seen had it occurred in this window.
    depth = CATALOG["depth_ppt"] * 1e-3
    axis.plot(
        [time.min(), ingress, ingress, egress, egress, time.max()],
        [1.0, 1.0, 1.0 - depth, 1.0 - depth, 1.0, 1.0],
        color=PALETTE["muted"],
        linewidth=2.0,
        linestyle=(0, (6, 3)),
        zorder=4,
        label="expected transit",
    )

    axis.scatter(
        time,
        flux,
        s=20,
        color="#7FA8BC",
        edgecolor="none",
        alpha=0.9,
        zorder=3,
        label=f"measurements ({len(used)})",
    )
    axis.scatter(
        drop_time,
        drop_flux,
        s=22,
        marker="x",
        color=PALETTE["muted"],
        linewidth=1.0,
        alpha=0.75,
        zorder=3,
        label=f"cut ({len(dropped)})",
    )
    axis.errorbar(
        bin_time,
        bin_flux,
        yerr=bin_error,
        fmt="o",
        markersize=8.5,
        color=PALETTE["gmu_green"],
        ecolor=PALETTE["gmu_green"],
        elinewidth=1.8,
        capsize=3,
        markeredgecolor="white",
        markeredgewidth=0.8,
        zorder=6,
        label="10-min bins",
    )

    axis.set_xlim(time.min() - 0.004, time.max() + 0.004)
    # Keep every binned point and its error bar inside the frame; the sequence
    # genuinely degrades at high airmass near the end of the night.
    lower = float(np.min(bin_flux - bin_error)) - 0.0009
    upper = max(float(np.max(bin_flux + bin_error)) + 0.0012, 1.0055)
    axis.set_ylim(lower, upper)
    axis.set_xlabel("Barycentric Julian Date (TDB) $-$ 2459782", fontsize=16)
    axis.set_ylabel("Relative brightness", fontsize=16)
    axis.set_title(
        "TOI-3505.01 from the Mason 0.8 m, UT 2022 July 22, R band, 50 s",
        fontsize=18,
        color=PALETTE["ink"],
        pad=34,
    )
    axis.annotate(
        "when the transit was scheduled",
        xy=(0.5 * (ingress + egress), upper - 0.0011),
        ha="center",
        fontsize=15,
        color=PALETTE["accent"],
    )
    # The key sits above the frame so it never covers data.
    axis.legend(
        loc="lower center",
        bbox_to_anchor=(0.5, 1.005),
        fontsize=14.5,
        frameon=False,
        ncol=4,
        handletextpad=0.3,
        columnspacing=1.4,
    )

    figure.tight_layout()
    figure.savefig(path, facecolor="white")
    figure.savefig(path.with_suffix(".svg"), facecolor="white")
    plt.close(figure)

    return {
        "points_plotted": len(used),
        "points_excluded": len(dropped),
        "window_depth_ppt": fixed["observed_depth_ppt"],
        "window_depth_error_ppt": fixed["observed_depth_error_ppt"],
        "injection_snr": fixed["injected_total_depth_snr"],
    }


def figure_field(path: Path) -> dict:
    """Plate-solved GMU field with the adopted apertures drawn on it.

    The seeing-stage figure in ``outputs/toi3505_seeing`` still shows the
    earlier 35-pixel trial aperture, so the poster figure is regenerated here
    at the adopted 25-pixel radius and 70-139-pixel sky annulus.
    """
    from astropy.io import fits
    from astropy.visualization import ZScaleInterval
    from matplotlib.patches import Circle, Rectangle

    frame = (
        ROOT
        / "data"
        / "ground"
        / "toi3505"
        / "plate_solved"
        / "TOI_3505.01_50.000s_R-0001_wcs.fits"
    )
    with fits.open(frame) as handle:
        image = handle[0].data.astype(float)

    settings = json.loads(
        (ROOT / "outputs" / "toi3505_final_candidate" / "analysis_settings.json").read_text()
    )
    aperture_file = (
        ROOT / "outputs" / "toi3505_final_candidate" / "TOI_3505.01_2022-07-22_R.apertures"
    )
    values = {}
    for line in aperture_file.read_text().splitlines():
        if "=" in line and not line.startswith("#"):
            key, _, value = line.partition("=")
            values[key.strip()] = value.strip()

    x_positions = [float(v) for v in values[".multiaperture.xapertures"].split(",")]
    y_aij = [float(v) for v in values[".multiaperture.yapertures"].split(",")]
    # AstroImageJ counts rows from the top; FITS counts from the bottom.
    height = image.shape[0]
    y_positions = [height - y for y in y_aij]

    radius = float(values[".aperture.radius"])
    back_inner = float(values[".aperture.rback1"])
    back_outer = float(values[".aperture.rback2"])
    pixel_scale = 0.3621236516728507

    interval = ZScaleInterval(contrast=0.16)
    low, high = interval.get_limits(image)

    figure, axes = plt.subplots(1, 2, figsize=(11.4, 5.7), dpi=200, width_ratios=[1.32, 1])

    axis = axes[0]
    axis.imshow(image, origin="lower", cmap="bone", vmin=low, vmax=high)
    for index, (x, y) in enumerate(zip(x_positions, y_positions)):
        is_target = index == 0
        color = PALETTE["gmu_gold"] if is_target else "#4FD1FF"
        axis.add_patch(
            Circle(
                (x, y),
                90,
                fill=False,
                edgecolor=color,
                linewidth=1.6 if is_target else 1.2,
            )
        )
        axis.text(
            x,
            y + 140,
            "T1" if is_target else f"C{index + 1}",
            ha="center",
            va="bottom",
            fontsize=14,
            color=color,
            fontweight="bold" if is_target else "normal",
        )

    zoom = 210
    axis.add_patch(
        Rectangle(
            (x_positions[0] - zoom, y_positions[0] - zoom),
            2 * zoom,
            2 * zoom,
            fill=False,
            edgecolor=PALETTE["accent"],
            linewidth=1.6,
            linestyle="--",
        )
    )
    axis.set_xlim(0, image.shape[1])
    axis.set_ylim(0, image.shape[0])
    axis.set_xticks([])
    axis.set_yticks([])
    axis.set_title(
        "Our star and the ten comparison stars",
        fontsize=16,
        color=PALETTE["ink"],
    )
    bar_pixels = 120.0 / pixel_scale
    axis.plot(
        [260, 260 + bar_pixels],
        [230, 230],
        color="white",
        linewidth=3.0,
        solid_capstyle="butt",
    )
    axis.annotate(
        "2 arcmin",
        (260 + bar_pixels / 2, 300),
        ha="center",
        fontsize=14,
        color="white",
    )

    axis = axes[1]
    axis.imshow(image, origin="lower", cmap="bone", vmin=low, vmax=high)
    for value, color, style in (
        (radius, PALETTE["gmu_gold"], "-"),
        (back_inner, "#4FD1FF", "--"),
        (back_outer, "#4FD1FF", "--"),
    ):
        axis.add_patch(
            Circle(
                (x_positions[0], y_positions[0]),
                value,
                fill=False,
                edgecolor=color,
                linewidth=1.8,
                linestyle=style,
            )
        )
    axis.set_xlim(x_positions[0] - zoom, x_positions[0] + zoom)
    axis.set_ylim(y_positions[0] - zoom, y_positions[0] + zoom)
    axis.set_xticks([])
    axis.set_yticks([])
    axis.set_title(
        f"Aperture: {radius:.0f} px star, {back_inner:.0f}-{back_outer:.0f} px sky ring",
        fontsize=16,
        color=PALETTE["ink"],
    )
    axis.annotate(
        f"{radius:.0f} px = {radius * pixel_scale:.1f} arcsec",
        (0.5, 0.035),
        xycoords="axes fraction",
        ha="center",
        fontsize=14,
        color="white",
    )

    figure.tight_layout()
    figure.savefig(path, facecolor="white")
    plt.close(figure)

    return {
        "source_radius_pixels": radius,
        "source_radius_arcsec": radius * pixel_scale,
        "sky_annulus_pixels": [back_inner, back_outer],
        "pixel_scale_arcsec": pixel_scale,
        "comparison_stars": len(x_positions) - 1,
        "settings_file": str(settings.get("source", "analysis_settings.json")),
    }


def figure_depth_consistency(path: Path) -> dict:
    """Per-sector depth against the catalog value, with the pipeline spread."""
    source = ROOT / "outputs" / "toi3505_tess_analysis" / "sector_measurements.csv"
    with source.open() as handle:
        rows = list(csv.DictReader(handle))

    sectors = [int(r["sector"]) for r in rows]
    depth = np.array([float(r["depth_ppt"]) for r in rows])
    error = np.array([float(r["depth_error_ppt"]) for r in rows])
    odd = np.array([float(r["odd_depth_ppt"]) for r in rows])
    odd_error = np.array([float(r["odd_depth_error_ppt"]) for r in rows])
    even = np.array([float(r["even_depth_ppt"]) for r in rows])
    even_error = np.array([float(r["even_depth_error_ppt"]) for r in rows])
    odd_even_sigma = np.array([float(r["odd_even_difference_sigma"]) for r in rows])
    secondary_snr = np.array([float(r["phase_0_5_snr"]) for r in rows])

    weight = 1.0 / error**2
    combined = float(np.sum(weight * depth) / np.sum(weight))
    combined_error = float(1.0 / math.sqrt(np.sum(weight)))
    chi2 = float(np.sum((depth - combined) ** 2 / error**2))

    figure, axes = plt.subplots(1, 2, figsize=(11.0, 4.6), dpi=220, width_ratios=[1.25, 1])
    for axis in axes:
        style_axis(axis)

    positions = np.arange(len(sectors), dtype=float)
    axis = axes[0]
    axis.axhspan(
        CATALOG["depth_ppt"] - CATALOG["depth_error_ppt"],
        CATALOG["depth_ppt"] + CATALOG["depth_error_ppt"],
        color=PALETTE["muted"],
        alpha=0.20,
        zorder=1,
        label="TOI catalog depth",
    )
    axis.axhline(CATALOG["depth_ppt"], color=PALETTE["muted"], linewidth=1.4, zorder=2)
    axis.axhspan(
        combined - combined_error,
        combined + combined_error,
        color=PALETTE["gmu_green"],
        alpha=0.20,
        zorder=1,
        label="This work, 4-sector weighted mean",
    )
    axis.axhline(combined, color=PALETTE["gmu_green"], linewidth=1.8, zorder=2)
    for index, sector in enumerate(sectors):
        axis.errorbar(
            positions[index],
            depth[index],
            yerr=error[index],
            fmt="o",
            markersize=9,
            color=PALETTE["sector"][sector],
            ecolor=PALETTE["sector"][sector],
            elinewidth=2.0,
            capsize=4,
            markeredgecolor="white",
            markeredgewidth=0.9,
            zorder=5,
        )
    axis.set_xticks(positions)
    axis.set_xticklabels([f"S{s}" for s in sectors], fontsize=12)
    axis.set_xlim(-0.6, len(sectors) - 0.4)
    axis.set_xlabel("TESS sector", fontsize=13)
    axis.set_ylabel("Transit depth (ppt)", fontsize=13)
    axis.set_title("Depth is stable across sectors", fontsize=14, color=PALETTE["ink"])
    axis.legend(loc="lower left", fontsize=9)

    axis = axes[1]
    axis.axhline(0.0, color=PALETTE["ink"], linewidth=1.0, zorder=3)
    axis.axhspan(-3, 3, color=PALETTE["gmu_green"], alpha=0.14, zorder=1)
    axis.errorbar(
        positions,
        odd - even,
        yerr=np.sqrt(odd_error**2 + even_error**2),
        fmt="D",
        markersize=8,
        color=PALETTE["gmu_green"],
        ecolor=PALETTE["gmu_green"],
        elinewidth=2.0,
        capsize=4,
        markeredgecolor="white",
        zorder=5,
        label="Odd $-$ even depth",
    )
    for index, sector in enumerate(sectors):
        axis.annotate(
            f"{odd_even_sigma[index]:+.2f}$\\sigma$",
            (positions[index], (odd - even)[index]),
            textcoords="offset points",
            xytext=(0, 14),
            ha="center",
            fontsize=10,
            color=PALETTE["ink"],
        )
    axis.set_xticks(positions)
    axis.set_xticklabels([f"S{s}" for s in sectors], fontsize=12)
    axis.set_xlim(-0.6, len(sectors) - 0.4)
    axis.set_ylim(-2.2, 2.6)
    axis.set_xlabel("TESS sector", fontsize=13)
    axis.set_ylabel("Depth difference (ppt)", fontsize=13)
    axis.set_title(
        "No odd-even depth difference above 1.7$\\sigma$", fontsize=14, color=PALETTE["ink"]
    )
    axis.legend(loc="lower right", fontsize=9)

    figure.tight_layout()
    figure.savefig(path, facecolor="white")
    figure.savefig(path.with_suffix(".svg"), facecolor="white")
    plt.close(figure)

    return {
        "combined_depth_ppt": combined,
        "combined_depth_error_ppt": combined_error,
        "chi2_about_combined": chi2,
        "dof": len(sectors) - 1,
        "max_absolute_odd_even_sigma": float(np.max(np.abs(odd_even_sigma))),
        "max_secondary_snr": float(np.max(secondary_snr)),
    }


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------


def main() -> None:
    events = load_events()

    raw = fit_ephemeris(events, "four-sector, formal errors", scale_errors=False)
    refined = fit_ephemeris(events, "four-sector, errors scaled", scale_errors=True)
    without_s14 = fit_ephemeris(
        [e for e in events if e["sector"] != 14],
        "sectors 41, 54, 81, errors scaled",
        scale_errors=True,
    )

    # Forward propagation to the 2026 symposium and one year beyond.
    horizons = {
        "2026-08-01": 2461254.0,
        "2027-08-01": 2461619.0,
        "2030-08-01": 2462715.0,
    }
    propagation = {}
    for name, bjd in horizons.items():
        cycle = round((bjd - CATALOG["epoch_bjd_tdb"]) / CATALOG["period_days"])
        propagation[name] = {
            "cycle": cycle,
            "catalog_uncertainty_minutes": propagation_uncertainty_minutes(
                CATALOG["epoch_error_days"], CATALOG["period_error_days"], 0.0, cycle
            ),
            "refined_uncertainty_minutes": propagation_uncertainty_minutes(
                refined.epoch_error_days,
                refined.period_error_days,
                refined.covariance_days2,
                cycle,
            ),
        }
        propagation[name]["improvement_factor"] = (
            propagation[name]["catalog_uncertainty_minutes"]
            / propagation[name]["refined_uncertainty_minutes"]
        )

    period_difference_seconds = (
        CATALOG["period_days"] - refined.period_days
    ) * 86400.0
    period_difference_sigma = abs(CATALOG["period_days"] - refined.period_days) / math.sqrt(
        CATALOG["period_error_days"] ** 2 + refined.period_error_days**2
    )
    # The catalog row is the softer comparison. SPOC's multi-sector fit is the
    # tightest published period, so the margin over it is the one that matters
    # and both belong on the board.
    spoc_difference_sigma = abs(SPOC["period_days"] - refined.period_days) / math.sqrt(
        SPOC["period_error_days"] ** 2 + refined.period_error_days**2
    )

    # Independent four-sector duration, weighted across the QLP box fits.
    with (ROOT / "outputs" / "toi3505_tess_analysis" / "sector_measurements.csv").open() as handle:
        sector_rows = list(csv.DictReader(handle))
    sector_durations = np.array([float(r["duration_hours"]) for r in sector_rows])
    sector_half_width = np.array(
        [
            0.5 * (float(r["duration_high_hours"]) - float(r["duration_low_hours"]))
            for r in sector_rows
        ]
    )
    duration_weight = 1.0 / np.maximum(sector_half_width, 1e-6) ** 2
    box_duration = float(np.sum(duration_weight * sector_durations) / np.sum(duration_weight))
    box_duration_error = float(1.0 / math.sqrt(np.sum(duration_weight)))

    # Stellar-density cross-check on the transit geometry.
    density_scaled = scaled_separation_from_density(
        SPOC["period_days"], STAR["assumed_mass_msun"], STAR["radius_rsun"]
    )
    density_scaled_low = scaled_separation_from_density(
        SPOC["period_days"],
        STAR["assumed_mass_msun"] - STAR["assumed_mass_error_msun"],
        STAR["radius_rsun"],
    )
    density_scaled_high = scaled_separation_from_density(
        SPOC["period_days"],
        STAR["assumed_mass_msun"] + STAR["assumed_mass_error_msun"],
        STAR["radius_rsun"],
    )
    duration_at_density = transit_duration_hours(
        SPOC["period_days"], density_scaled, SPOC["radius_ratio"], SPOC["impact_parameter"]
    )
    density_impact = impact_from_duration(
        SPOC["period_days"], density_scaled, SPOC["radius_ratio"], box_duration
    )
    spoc_scaled = scaled_separation_from_duration(
        SPOC["period_days"],
        SPOC["duration_hours"],
        SPOC["radius_ratio"],
        SPOC["impact_parameter"],
    )
    implied_density_ratio = (density_scaled / spoc_scaled) ** 3

    geometry = {
        "assumed_stellar_mass_msun": STAR["assumed_mass_msun"],
        "assumed_stellar_mass_error_msun": STAR["assumed_mass_error_msun"],
        "a_over_rstar_from_stellar_density": density_scaled,
        "a_over_rstar_from_stellar_density_range": [density_scaled_low, density_scaled_high],
        "a_over_rstar_from_spoc_duration": spoc_scaled,
        "spoc_duration_hours": SPOC["duration_hours"],
        "duration_predicted_at_stellar_density_hours": duration_at_density,
        "four_sector_box_duration_hours": box_duration,
        "four_sector_box_duration_error_hours": box_duration_error,
        "toi_catalog_duration_hours": CATALOG["duration_hours"],
        "spoc_minus_box_duration_hours": SPOC["duration_hours"] - box_duration,
        "spoc_density_deficit_factor": implied_density_ratio,
        "impact_parameter_at_stellar_density": density_impact,
        "spoc_impact_parameter": SPOC["impact_parameter"],
        "box_duration_caveat": (
            "A box fit to a V-shaped grazing transit measures something closer "
            "to the flat-bottom width than to first-to-fourth contact, so the "
            "box and Mandel-Agol durations are not measuring the same quantity. "
            "The duration numbers are reported side by side, not differenced as "
            "a significance test."
        ),
        "interpretation": (
            "The official solution is grazing at b = "
            f"{SPOC['impact_parameter']:.2f}. Its first-to-fourth-contact duration "
            "requires a mean stellar density about "
            f"{implied_density_ratio:.1f} times lower than the TIC v8 radius and an "
            "assumed 1.25 solar-mass star imply, which is the familiar signature of "
            "the grazing degeneracy between depth, duration, impact parameter, and "
            "a/R*. Holding the stellar density fixed and matching the measured "
            f"duration instead gives b = {density_impact:.2f} and a/R* = "
            f"{density_scaled:.2f}. Because the transit is grazing under either "
            "solution, the radius ratio is weakly constrained and the formal error "
            "understates the true uncertainty on the planet radius."
        ),
    }

    # Derived parameters from the SPOC geometric transit model.
    observed = derive_system(
        SPOC["radius_ratio"],
        SPOC["radius_ratio_error"],
        SPOC["impact_parameter"],
        SPOC["duration_hours"],
        SPOC["period_days"],
        "SPOC multi-sector fit, no dilution correction",
        "Radius ratio, impact parameter, and duration are the official "
        "Mandel-Agol fit; a/R* is inverted from the duration equation.",
    )

    # Same geometry, but with the depth rescaled to the host star under the
    # documented screening dilution scenario for the 0.517-arcsec companion.
    dilution = json.loads(
        (ROOT / "outputs" / "toi3505_tess_pixels" / "dilution_screen.json").read_text()
    )
    host_depth = dilution["if_2p91_ppt_is_an_uncorrected_observed_depth"][
        "target_host_depth_ppt"
    ]
    depth_scale = math.sqrt(host_depth / CATALOG["depth_ppt"])
    corrected = derive_system(
        SPOC["radius_ratio"] * depth_scale,
        SPOC["radius_ratio_error"] * depth_scale,
        SPOC["impact_parameter"],
        SPOC["duration_hours"],
        SPOC["period_days"],
        "screening dilution scenario, not adopted",
        "Upper bound only. Scales the radius ratio by the screening flux "
        "budget for the 0.517-arcsec companion plus catalog contamination.",
    )

    # Geometry held to the stellar density, with b solved from the independent
    # four-sector duration. This is the solution the poster adopts.
    density_constrained = derive_system(
        SPOC["radius_ratio"],
        SPOC["radius_ratio_error"],
        density_impact,
        box_duration,
        SPOC["period_days"],
        "density-constrained, four-sector duration",
        "Adopts the TIC v8 stellar radius and a 1.25 solar-mass star, then "
        "solves for the impact parameter that reproduces the measured duration.",
    )

    # Figures.
    figure_timing(events, refined, OUT / "01_transit_timing.png")
    neb = figure_neb(OUT / "02_nearby_star_screen.png")
    depth = figure_depth_consistency(OUT / "03_depth_consistency.png")
    field = figure_field(OUT / "04_field_and_aperture.png")
    ground = figure_ground_light_curve(OUT / "05_ground_light_curve.png")
    folded = figure_phase_folded(OUT / "06_phase_folded.png", refined)

    baseline_days = max(e["measured_bjd"] for e in events) - min(
        e["measured_bjd"] for e in events
    )

    summary = {
        "target": "TOI-3505.01",
        "catalog": CATALOG,
        "star": STAR,
        "spoc_multisector": SPOC,
        "companion": COMPANION,
        "ephemeris": {
            "formal": asdict(raw),
            "adopted_error_scaled": asdict(refined),
            "without_sector_14": asdict(without_s14),
            "baseline_days": baseline_days,
            "baseline_years": baseline_days / 365.25,
            "period_difference_from_catalog_seconds": period_difference_seconds,
            "period_difference_sigma": period_difference_sigma,
            "period_difference_from_spoc_sigma": spoc_difference_sigma,
            "period_precision_gain": CATALOG["period_error_days"]
            / refined.period_error_days,
            "period_precision_gain_over_spoc": SPOC["period_error_days"]
            / refined.period_error_days,
            "pipeline_per_sector": {
                sector: pipeline
                for sector, pipeline in sorted(
                    {e["sector"]: e["pipeline"] for e in events}.items()
                )
            },
            "forward_propagation": propagation,
        },
        "transit_geometry": geometry,
        "derived_parameters": {
            "observed": observed,
            "dilution_scenario": corrected,
            "density_constrained": density_constrained,
        },
        "nearby_star_screen": neb,
        "depth_consistency": depth,
        "ground_field": field,
        "ground_light_curve": ground,
        "phase_folded": folded,
    }

    (OUT / "poster_analysis.json").write_text(json.dumps(summary, indent=2))

    # Console report.
    print("=" * 78)
    print("TOI-3505.01 poster analysis")
    print("=" * 78)
    print(f"Timing baseline: {baseline_days:.1f} d ({baseline_days/365.25:.2f} yr), {len(events)} events")
    print()
    for fit in (raw, refined, without_s14):
        print(f"  {fit.label}")
        print(f"    P  = {fit.period_days:.8f} +/- {fit.period_error_days:.3e} d")
        print(f"    T0 = {fit.epoch_bjd_tdb:.6f} +/- {fit.epoch_error_days:.6f} BJD_TDB")
        print(f"    reduced chi2 = {fit.reduced_chi2:.2f}, RMS = {fit.residual_rms_minutes:.2f} min, scale = {fit.error_scale:.2f}")
    print()
    print(f"Catalog P = {CATALOG['period_days']:.8f} +/- {CATALOG['period_error_days']:.3e} d")
    print(f"Difference = {period_difference_seconds:+.3f} s ({period_difference_sigma:.2f} sigma)")
    print(f"Period precision gain = {CATALOG['period_error_days']/refined.period_error_days:.2f}x")
    print()
    for name, entry in propagation.items():
        print(
            f"  {name}: catalog {entry['catalog_uncertainty_minutes']:6.1f} min -> "
            f"refined {entry['refined_uncertainty_minutes']:6.1f} min "
            f"({entry['improvement_factor']:.1f}x tighter)"
        )
    print()
    print("Depth consistency:")
    print(f"  combined = {depth['combined_depth_ppt']:.3f} +/- {depth['combined_depth_error_ppt']:.3f} ppt")
    print(f"  chi2/dof about combined = {depth['chi2_about_combined']:.2f}/{depth['dof']}")
    print(f"  max |odd-even| = {depth['max_absolute_odd_even_sigma']:.2f} sigma")
    print(f"  max secondary SNR = {depth['max_secondary_snr']:.2f}")
    print()
    print("Transit geometry:")
    print(f"  a/R* from stellar density   = {geometry['a_over_rstar_from_stellar_density']:.2f}")
    print(f"  a/R* from SPOC duration     = {geometry['a_over_rstar_from_spoc_duration']:.2f}")
    print(f"  SPOC duration               = {SPOC['duration_hours']:.3f} +/- {SPOC['duration_error_hours']:.3f} h")
    print(f"  four-sector box duration    = {box_duration:.3f} +/- {box_duration_error:.3f} h")
    print(f"  TOI catalog duration        = {CATALOG['duration_hours']:.3f} +/- {CATALOG['duration_error_hours']:.3f} h")
    print(f"  SPOC - box                  = {geometry['spoc_minus_box_duration_hours']:+.3f} h (different definitions, not a significance test)")
    print(f"  SPOC density deficit        = {geometry['spoc_density_deficit_factor']:.1f}x below the expected stellar density")
    print(f"  b at stellar density        = {density_impact:.3f} (SPOC: {SPOC['impact_parameter']:.3f})")
    print()
    print("Derived parameters:")
    for entry in (observed, density_constrained, corrected):
        print(f"  {entry['case']}")
        print(f"    Rp/R* = {entry['radius_ratio']:.5f} +/- {entry['radius_ratio_error']:.5f}")
        print(f"    Rp    = {entry['planet_radius_rearth']:.2f} +/- {entry['planet_radius_rearth_error']:.2f} R_earth"
              f"  ({entry['planet_radius_rjup']:.3f} +/- {entry['planet_radius_rjup_error']:.3f} R_jup)")
        print(f"    b     = {entry['impact_parameter']:.3f}, i = {entry['inclination_degrees']:.2f} deg")
        print(f"    a/R*  = {entry['a_over_rstar']:.2f}, a = {entry['semimajor_axis_au']:.4f} au")
        print(f"    implied M* = {entry['implied_stellar_mass_msun']:.2f} M_sun (consistency check)")
        print(f"    Teq   = {entry['equilibrium_temperature_k']:.0f} K, S = {entry['insolation_earth']:.0f} S_earth")
    print()
    print("Nearby-star screen:", json.dumps(neb, indent=2))
    print()
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
