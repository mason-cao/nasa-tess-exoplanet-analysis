"""Search the whole 2022 TOI-3505.01 ground sequence for a transit at any time.

The fixed-window test in ``make_toi3505_ground_checks.py`` answers a narrow
question: was there a dimming inside the exact window the stale 2022 schedule
told the observers to cover?  Because that window is now known to be stale, the
more useful question is the one a reviewer asks next -- is there anything like
the expected transit anywhere in the sequence?

This script answers it directly.  The event duration is held at a published
TESS value and the midpoint is scanned across the observed span.  At each trial
midpoint a straight local baseline and one exposure-integrated box depth are
fitted, so the depth is free while the shape and timing are fixed.  The same
scan is then repeated after injecting the corresponding published depth at
each admissible midpoint.  That injection check retains the actual temporal
systematics and prevents formal white-noise covariance from being mistaken for
uniform sensitivity.

The result is a phase-sampling diagnostic for this one sequence.  It does not
provide a calibrated false-alarm probability, resolve the close companion,
correct for dilution, or convert the sequence into a detection.
"""

from __future__ import annotations

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

from make_toi3505_ground_checks import fit_linear_box, integrated_box_fraction


ROOT = Path(__file__).resolve().parents[1]
LIGHT_CURVE_PATH = (
    ROOT
    / "outputs"
    / "toi3505_final_candidate"
    / "TOI_3505.01_2022-07-22_R_final_candidate_light_curve.csv"
)
SCHEDULE_CHECK_PATH = (
    ROOT / "outputs" / "toi3505_final_candidate" / "historical_schedule_check.json"
)
TESS_SUMMARY_PATH = ROOT / "outputs" / "toi3505_tess_analysis" / "analysis_summary.json"
VALIDATION_SUMMARY_PATH = (
    ROOT / "outputs" / "toi3505_data_validation" / "analysis_summary.json"
)
OUTPUT_DIR = ROOT / "outputs" / "toi3505_ground_search"

EXPOSURE_SECONDS = 50.0
MIDPOINT_STEP_MINUTES = 2.0
MINIMUM_IN_EVENT_POINTS = 30
# A local straight baseline can only be trusted when it is anchored on both
# sides of the candidate event.  Without this, a trial box that hangs off the
# end of the sequence lets the slope trade against depth, and the fit reports a
# large spurious dimming driven entirely by the degraded end of the night.
MINIMUM_BASELINE_POINTS_EACH_SIDE = 20
# Fraction of the trial event window that must fall inside the observed span.
MINIMUM_EVENT_COVERAGE = 0.85
FORMAL_DETECTION_SIGMA = 3.0


def load_object(path: Path) -> dict[str, object]:
    """Load one JSON object and reject other top-level types."""
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError(f"Expected a JSON object in {path}")
    return payload


def duration_hypotheses() -> list[dict[str, object]]:
    """Return the published TESS durations and depths the scan holds fixed."""
    tess = load_object(TESS_SUMMARY_PATH)
    validation = load_object(VALIDATION_SUMMARY_PATH)
    multisector = validation["official_multisector_tce"]
    if not isinstance(multisector, list) or len(multisector) != 1:
        raise ValueError("Expected exactly one SPOC multi-sector TCE record")
    spoc = multisector[0]
    return [
        {
            "label": "TOI catalog",
            "duration_hours": float(tess["catalog_duration_hours"]),
            "expected_depth_ppt": float(tess["catalog_depth_ppt"]),
            "expected_depth_error_ppt": None,
            "source": "ExoFOP TOI catalog row frozen for this analysis",
        },
        {
            "label": "SPOC multi-sector",
            "duration_hours": float(spoc["fit_duration_hours"]),
            "expected_depth_ppt": float(spoc["fit_depth_ppt"]),
            "expected_depth_error_ppt": float(spoc["fit_depth_error_ppt"]),
            "source": "SPOC sectors 54 and 81 combined Mandel-Agol fit",
        },
    ]


def scan_one_duration(
    time_hours: np.ndarray,
    flux: np.ndarray,
    flux_error: np.ndarray,
    use: np.ndarray,
    duration_hours: float,
) -> pd.DataFrame:
    """Fit a free depth at every trial midpoint for one fixed duration."""
    observed = time_hours[use]
    step = MIDPOINT_STEP_MINUTES / 60.0
    span_start = float(np.min(observed))
    span_end = float(np.max(observed))
    midpoints = np.arange(span_start, span_end + 0.5 * step, step)
    rows: list[dict[str, float | int]] = []
    for midpoint in midpoints:
        event_start = float(midpoint) - duration_hours / 2.0
        event_end = float(midpoint) + duration_hours / 2.0
        inside = max(0.0, min(event_end, span_end) - max(event_start, span_start))
        coverage = inside / duration_hours
        if coverage < MINIMUM_EVENT_COVERAGE:
            continue
        box = integrated_box_fraction(
            time_hours, float(midpoint), duration_hours, EXPOSURE_SECONDS / 3600.0
        )
        in_event = int(np.sum((box > 0.0) & use))
        before = int(np.sum((box == 0.0) & use & (time_hours < event_start)))
        after = int(np.sum((box == 0.0) & use & (time_hours > event_end)))
        if (
            in_event < MINIMUM_IN_EVENT_POINTS
            or before < MINIMUM_BASELINE_POINTS_EACH_SIDE
            or after < MINIMUM_BASELINE_POINTS_EACH_SIDE
        ):
            continue
        try:
            fit = fit_linear_box(time_hours, flux, flux_error, box, use)
        except (ValueError, np.linalg.LinAlgError):
            continue
        depth_ppt = fit["depth"] * 1000.0
        error_ppt = fit["depth_error"] * 1000.0
        rows.append(
            {
                "midpoint_hours_since_first_image": float(midpoint),
                "event_coverage_fraction": coverage,
                "depth_ppt": depth_ppt,
                "depth_error_ppt": error_ppt,
                "depth_snr": depth_ppt / error_ppt if error_ppt > 0 else float("nan"),
                "formal_upper_limit_ppt": (
                    depth_ppt + FORMAL_DETECTION_SIGMA * error_ppt
                ),
                "in_event_points": in_event,
                "baseline_points_before": before,
                "baseline_points_after": after,
                "reduced_chi_square": fit["reduced_chi_square"],
                "residual_scatter_ppt": fit["residual_scatter"] * 1000.0,
            }
        )
    if not rows:
        raise RuntimeError(
            f"No trial midpoint met the coverage requirements for a "
            f"{duration_hours:.3f} h duration"
        )
    return pd.DataFrame(rows)


def assess_expected_depth_recovery(
    time_hours: np.ndarray,
    flux: np.ndarray,
    flux_error: np.ndarray,
    use: np.ndarray,
    duration_hours: float,
    expected_depth_ppt: float,
    scan: pd.DataFrame,
) -> pd.DataFrame:
    """Inject the published depth at every admissible midpoint and refit it.

    The injected curve preserves the observed residual structure.  The
    resulting fraction above the formal threshold is a phase-sampling
    diagnostic for this particular realization, not a frequentist detection
    probability or completeness estimate for an ensemble of observing nights.
    """
    assessed = scan.copy()
    total_depths: list[float] = []
    depth_errors: list[float] = []
    depth_snrs: list[float] = []
    recovered_increments: list[float] = []
    injected_depth = expected_depth_ppt / 1000.0
    exposure_hours = EXPOSURE_SECONDS / 3600.0

    for row in assessed.itertuples(index=False):
        midpoint = float(row.midpoint_hours_since_first_image)
        box = integrated_box_fraction(
            time_hours, midpoint, duration_hours, exposure_hours
        )
        injected_flux = flux * (1.0 - injected_depth * box)
        fit = fit_linear_box(time_hours, injected_flux, flux_error, box, use)
        total_depth_ppt = float(fit["depth"]) * 1000.0
        depth_error_ppt = float(fit["depth_error"]) * 1000.0
        depth_snr = (
            total_depth_ppt / depth_error_ppt if depth_error_ppt > 0.0 else float("nan")
        )
        total_depths.append(total_depth_ppt)
        depth_errors.append(depth_error_ppt)
        depth_snrs.append(depth_snr)
        recovered_increments.append(total_depth_ppt - float(row.depth_ppt))

    assessed["injected_expected_depth_ppt"] = expected_depth_ppt
    assessed["injected_total_depth_ppt"] = total_depths
    assessed["injected_depth_error_ppt"] = depth_errors
    assessed["injected_total_depth_snr"] = depth_snrs
    assessed["injected_recovered_increment_ppt"] = recovered_increments
    assessed["injected_above_formal_3sigma"] = (
        assessed["injected_total_depth_snr"] >= FORMAL_DETECTION_SIGMA
    )
    return assessed


def summarize(scan: pd.DataFrame, hypothesis: dict[str, object]) -> dict[str, object]:
    """Reduce one scan to formal fit metrics and injection recovery."""
    duration_hours = float(hypothesis["duration_hours"])
    expected_depth = float(hypothesis["expected_depth_ppt"])
    best = scan.loc[scan["depth_snr"].idxmax()]
    worst_limit = float(scan["formal_upper_limit_ppt"].max())
    recovered = scan["injected_above_formal_3sigma"].astype(bool)
    recovered_count = int(recovered.sum())
    trial_count = int(len(scan))
    return {
        "label": hypothesis["label"],
        "source": hypothesis["source"],
        "duration_hours": duration_hours,
        "searched_midpoint_range_hours": [
            float(scan["midpoint_hours_since_first_image"].min()),
            float(scan["midpoint_hours_since_first_image"].max()),
        ],
        "expected_depth_ppt": expected_depth,
        "expected_depth_error_ppt": hypothesis["expected_depth_error_ppt"],
        "trial_midpoints": trial_count,
        "midpoint_step_minutes": MIDPOINT_STEP_MINUTES,
        "best_midpoint_hours_since_first_image": float(
            best["midpoint_hours_since_first_image"]
        ),
        "best_depth_ppt": float(best["depth_ppt"]),
        "best_depth_error_ppt": float(best["depth_error_ppt"]),
        "best_depth_snr": float(best["depth_snr"]),
        "median_formal_upper_limit_ppt": float(scan["formal_upper_limit_ppt"].median()),
        "maximum_formal_upper_limit_ppt": worst_limit,
        "formal_upper_limit_below_expected_depth_everywhere": bool(
            worst_limit < expected_depth
        ),
        "formal_upper_limit_sigma": FORMAL_DETECTION_SIGMA,
        "expected_depth_formal_recovery_count": recovered_count,
        "expected_depth_formal_recovery_total": trial_count,
        "expected_depth_formal_recovery_fraction": recovered_count / trial_count,
        "expected_depth_recovered_above_formal_3sigma_everywhere": bool(
            recovered.all()
        ),
        "minimum_injected_total_depth_snr": float(
            scan["injected_total_depth_snr"].min()
        ),
        "median_injected_total_depth_snr": float(
            scan["injected_total_depth_snr"].median()
        ),
        "maximum_injected_total_depth_snr": float(
            scan["injected_total_depth_snr"].max()
        ),
        "median_injected_recovered_increment_ppt": float(
            scan["injected_recovered_increment_ppt"].median()
        ),
        "interpretation": (
            "Formal weighted-least-squares upper limits do not model temporal "
            "correlation. Injection recovery is therefore used to describe "
            "phase sensitivity, and no global exclusion is claimed."
        ),
    }


def plot_search(
    light_curve: pd.DataFrame,
    scans: dict[str, pd.DataFrame],
    summaries: list[dict[str, object]],
    schedule: dict[str, object],
    path: Path,
) -> None:
    """Draw the searched light curve and the depth-versus-midpoint scan."""
    use = light_curve["used_in_primary_curve"].to_numpy(dtype=bool)
    time = light_curve["hours_since_first_image"].to_numpy(dtype=float)
    flux = light_curve["adopted_relative_brightness"].to_numpy(dtype=float)

    figure, axes = plt.subplots(
        2, 1, figsize=(9.0, 7.2), sharex=True, height_ratios=(1.0, 1.15)
    )

    window = schedule_window_hours(schedule, light_curve)
    for axis in axes:
        if window is not None:
            axis.axvspan(
                window[0],
                window[1],
                color="#C0392B",
                alpha=0.12,
                zorder=0,
                label="2022 schedule window",
            )

    axes[0].plot(
        time[use],
        flux[use],
        ".",
        color="#2C6B5F",
        markersize=3.4,
        alpha=0.55,
        label=f"{int(use.sum())} accepted measurements",
    )
    axes[0].set_ylabel("Relative brightness")
    axes[0].set_title(
        "TOI-3505.01 ground sequence searched for a transit at any time",
        fontsize=11,
    )

    colors = {"TOI catalog": "#8A6110", "SPOC multi-sector": "#1F4E45"}
    for summary in summaries:
        label = str(summary["label"])
        scan = scans[label]
        color = colors.get(label, "#444C48")
        axes[1].plot(
            scan["midpoint_hours_since_first_image"],
            scan["depth_ppt"],
            "-",
            color=color,
            linewidth=1.4,
            label=(
                f"{label} ({summary['duration_hours']:.2f} h; injection "
                f"recovery {float(summary['expected_depth_formal_recovery_fraction']):.0%})"
            ),
        )
        axes[1].fill_between(
            scan["midpoint_hours_since_first_image"],
            scan["depth_ppt"] - scan["depth_error_ppt"],
            scan["depth_ppt"] + scan["depth_error_ppt"],
            color=color,
            alpha=0.16,
            linewidth=0,
        )
        axes[1].axhline(
            float(summary["expected_depth_ppt"]),
            color=color,
            linestyle=":",
            linewidth=1.2,
        )

    axes[1].axhline(0.0, color="#6B736E", linewidth=0.9)
    axes[1].set_xlabel("Hours since first exposure")
    axes[1].set_ylabel("Fitted depth (ppt)")
    axes[1].legend(loc="upper left", fontsize=8, framealpha=0.9)
    axes[0].legend(loc="lower left", fontsize=8, framealpha=0.9)
    axes[1].annotate(
        "Dotted lines mark the published TESS depths.\n"
        "Shading is the formal WLS 1-sigma interval.\n"
        "Positive depth means dimming.",
        xy=(0.99, 0.04),
        xycoords="axes fraction",
        ha="right",
        va="bottom",
        fontsize=8,
        color="#444C48",
    )

    figure.tight_layout()
    figure.savefig(path, dpi=200)
    svg_path = path.with_suffix(".svg")
    figure.savefig(svg_path)
    # Matplotlib leaves trailing spaces in SVG path data.  Removing them keeps
    # generated artifacts clean under ``git diff --check`` without changing
    # the rendered figure.
    svg_path.write_text(
        "\n".join(line.rstrip() for line in svg_path.read_text().splitlines()) + "\n",
        encoding="utf-8",
    )
    plt.close(figure)


def schedule_window_hours(
    schedule: dict[str, object], light_curve: pd.DataFrame
) -> tuple[float, float] | None:
    """Convert the stale schedule ingress and egress to plot hours."""
    interpretation = schedule.get("working_interpretation")
    if not isinstance(interpretation, dict):
        return None
    times = interpretation.get("times")
    if not isinstance(times, dict):
        return None
    try:
        ingress = float(times["ingress"]["bjd_tdb"])  # type: ignore[index]
        egress = float(times["egress"]["bjd_tdb"])  # type: ignore[index]
    except (KeyError, TypeError, ValueError):
        return None
    first = float(light_curve["bjd_tdb"].min())
    return ((ingress - first) * 24.0, (egress - first) * 24.0)


def write_readme(summaries: list[dict[str, object]], path: Path) -> None:
    """Record what the scan tested and what it does not establish."""
    lines = [
        "# TOI-3505.01 whole-sequence transit search",
        "",
        "The fixed-window test asks whether the 2022 schedule window contained a",
        "dimming. Because that window is now known to be stale, this scan asks the",
        "broader question: is there anything like the expected transit anywhere in",
        "the sequence? The duration is held at a published TESS value and the",
        "midpoint is scanned across the observed span with a free depth.",
        "",
        "## Result",
        "",
    ]
    for summary in summaries:
        limit = float(summary["maximum_formal_upper_limit_ppt"])
        expected = float(summary["expected_depth_ppt"])
        searched = summary["searched_midpoint_range_hours"]
        if not isinstance(searched, list) or len(searched) != 2:
            raise RuntimeError("A search summary has an invalid midpoint range")
        recovered = int(summary["expected_depth_formal_recovery_count"])
        total = int(summary["expected_depth_formal_recovery_total"])
        recovery_fraction = float(summary["expected_depth_formal_recovery_fraction"])
        lines.extend(
            [
                f"### {summary['label']} duration ({summary['duration_hours']:.3f} h)",
                "",
                f"- Searched midpoints {searched[0]:.2f} to {searched[1]:.2f} h "
                f"after the first exposure: {summary['trial_midpoints']} trials at "
                f"{summary['midpoint_step_minutes']:.0f}-minute spacing. The range "
                "is set by requiring the event to be at least "
                f"{MINIMUM_EVENT_COVERAGE:.0%} sampled with at least "
                f"{MINIMUM_BASELINE_POINTS_EACH_SIDE} baseline points on each side.",
                f"- Deepest dimming: {float(summary['best_depth_ppt']):.3f} +/- "
                f"{float(summary['best_depth_error_ppt']):.3f} ppt at "
                f"{float(summary['best_midpoint_hours_since_first_image']):.2f} h "
                f"({float(summary['best_depth_snr']):.2f} sigma).",
                f"- Median formal {summary['formal_upper_limit_sigma']}-sigma "
                "WLS upper bound "
                f"{float(summary['median_formal_upper_limit_ppt']):.3f} ppt; "
                f"weakest {limit:.3f} ppt, against a published {expected:.3f} "
                "ppt depth. These formal bounds do not include temporal "
                "covariance and are not treated as global exclusions.",
                f"- Injecting the published depth into the observed curve at "
                f"each admissible midpoint produces a formal >=3-sigma fitted "
                f"depth at {recovered}/{total} midpoints "
                f"({recovery_fraction:.1%}); the minimum injected-event formal "
                "S/N is "
                f"{float(summary['minimum_injected_total_depth_snr']):.2f}.",
                "- Because exact-depth recovery is not complete across phase, "
                "this analysis does not claim that the published depth is "
                "excluded at every midpoint.",
                "",
            ]
        )
    lines.extend(
        [
            "## Limits",
            "",
            "- The scan holds the duration fixed and assumes a symmetric box; it is",
            "  not a limb-darkened physical transit fit.",
            "- Depths are observed aperture depths. No dilution correction is",
            "  applied, and the 0.517-arcsecond companion is not resolved.",
            "- Formal WLS intervals are locally scaled by reduced chi-square but",
            "  do not model time-correlated residuals. No false-alarm probability",
            "  is quoted.",
            "- The injection-recovery percentage is a phase-sampling diagnostic",
            "  for this observed residual realization, not a calibrated detection",
            "  probability or an ensemble completeness estimate.",
            "- The sequence is only about twice the transit duration, so a fully",
            "  sampled event with baseline on both sides fits only near the middle",
            "  of the night. The search says nothing about events outside that",
            "  range, including any that fall entirely outside the sequence.",
            "- The degraded final half hour of the night, where the scatter is",
            "  about four times the mid-run value, produces a large spurious",
            "  dimming if trial events are allowed to hang off the end of the",
            "  sequence. The two-sided baseline requirement excludes those trials.",
            "",
            "## Products",
            "",
            "- `ground_search.json` - inputs, per-duration summaries, and limits.",
            "- `floating_time_scan.csv` - every trial midpoint and fitted depth.",
            "- `01_ground_transit_search.png` and `.svg` - publication figure.",
            "",
            "Regenerate with:",
            "",
            "```bash",
            ".venv/bin/python src/search_toi3505_ground_transit.py",
            "```",
            "",
        ]
    )
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    light_curve = pd.read_csv(LIGHT_CURVE_PATH)
    schedule = load_object(SCHEDULE_CHECK_PATH)

    time_hours = light_curve["hours_since_first_image"].to_numpy(dtype=float)
    flux = light_curve["adopted_relative_brightness"].to_numpy(dtype=float)
    flux_error = light_curve["raw_relative_brightness_error"].to_numpy(dtype=float)
    use = light_curve["used_in_primary_curve"].to_numpy(dtype=bool)

    scans: dict[str, pd.DataFrame] = {}
    summaries: list[dict[str, object]] = []
    frames: list[pd.DataFrame] = []
    for hypothesis in duration_hypotheses():
        label = str(hypothesis["label"])
        scan = scan_one_duration(
            time_hours, flux, flux_error, use, float(hypothesis["duration_hours"])
        )
        scan = assess_expected_depth_recovery(
            time_hours,
            flux,
            flux_error,
            use,
            float(hypothesis["duration_hours"]),
            float(hypothesis["expected_depth_ppt"]),
            scan,
        )
        scans[label] = scan
        summaries.append(summarize(scan, hypothesis))
        labelled = scan.copy()
        labelled.insert(0, "duration_label", label)
        labelled.insert(1, "duration_hours", float(hypothesis["duration_hours"]))
        frames.append(labelled)

    combined = pd.concat(frames, ignore_index=True)
    combined.to_csv(OUTPUT_DIR / "floating_time_scan.csv", index=False)

    result = {
        "target": "TOI-3505.01",
        "analysis_type": "whole-sequence floating-midpoint transit search",
        "source_files": {
            "light_curve": str(LIGHT_CURVE_PATH.relative_to(ROOT)),
            "tess_summary": str(TESS_SUMMARY_PATH.relative_to(ROOT)),
            "validation_summary": str(VALIDATION_SUMMARY_PATH.relative_to(ROOT)),
        },
        "method": (
            "At each trial midpoint a straight local baseline and one "
            "exposure-integrated box depth are fitted by weighted least squares. "
            "The duration is held at a published TESS value and the depth is free. "
            "The corresponding published depth is then injected into the actual "
            "observed curve at every admissible midpoint and refitted."
        ),
        "accepted_measurements": int(use.sum()),
        "observed_span_hours": float(time_hours[use].max() - time_hours[use].min()),
        "minimum_in_event_points": MINIMUM_IN_EVENT_POINTS,
        "minimum_baseline_points_each_side": MINIMUM_BASELINE_POINTS_EACH_SIDE,
        "minimum_event_coverage_fraction": MINIMUM_EVENT_COVERAGE,
        "durations": summaries,
        "limits": [
            "Observed aperture depths only; no dilution correction is applied.",
            "A fixed symmetric box, not a limb-darkened physical transit model.",
            "Formal WLS covariance does not model time-correlated residuals.",
            "Injection recovery is a phase diagnostic, not a detection probability.",
            "A null over the searched span says nothing about events outside it.",
        ],
    }
    (OUTPUT_DIR / "ground_search.json").write_text(
        json.dumps(result, indent=2) + "\n", encoding="utf-8"
    )
    plot_search(
        light_curve,
        scans,
        summaries,
        schedule,
        OUTPUT_DIR / "01_ground_transit_search.png",
    )
    write_readme(summaries, OUTPUT_DIR / "README.md")

    print("TOI-3505.01 whole-sequence transit search")
    for summary in summaries:
        print(
            f"  {summary['label']} ({float(summary['duration_hours']):.3f} h): "
            f"best {float(summary['best_depth_ppt']):.3f} +/- "
            f"{float(summary['best_depth_error_ppt']):.3f} ppt "
            f"({float(summary['best_depth_snr']):.2f} formal sigma)"
        )
        print(
            f"    weakest formal {FORMAL_DETECTION_SIGMA:.0f}-sigma bound "
            f"{float(summary['maximum_formal_upper_limit_ppt']):.3f} ppt vs published "
            f"{float(summary['expected_depth_ppt']):.3f} ppt"
        )
        print(
            "    exact-depth injection formal recovery "
            f"{int(summary['expected_depth_formal_recovery_count'])}/"
            f"{int(summary['expected_depth_formal_recovery_total'])} "
            f"({float(summary['expected_depth_formal_recovery_fraction']):.1%})"
        )
    print(f"Wrote {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
