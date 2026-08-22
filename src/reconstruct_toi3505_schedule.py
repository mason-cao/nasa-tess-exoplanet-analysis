"""Reconstruct the ephemeris behind the 2022 TOI-3505.01 schedule window.

The preserved internship spreadsheet row contains ingress, egress, and a
period, but not the epoch or the original formulas. Public ExoFOP reports
preserve an earlier predicted transit and both an old and revised period. This
script tests whether propagating that earlier event with the old period
reproduces the spreadsheet markers. The result is an archival consistency
test, not proof of the missing workbook's cell history.
"""

from __future__ import annotations

import csv
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
from matplotlib.ticker import MultipleLocator


ROOT = Path(__file__).resolve().parents[1]
SCHEDULE_CHECK_PATH = (
    ROOT / "outputs" / "toi3505_final_candidate" / "historical_schedule_check.json"
)
GROUND_SUMMARY_PATH = ROOT / "outputs" / "toi3505_final_candidate" / "summary.json"
EPHEMERIS_PATH = (
    ROOT / "outputs" / "toi3505_ephemeris_refined" / "ephemeris_refined.json"
)
EXOFOP_PATH = ROOT / "data" / "catalogs" / "toi3505" / "exofop_ground_followup.json"
OUTPUT_DIR = ROOT / "outputs" / "toi3505_schedule_reconstruction"

MARKERS = ("ingress", "midpoint", "egress")


def load_object(path: Path) -> dict[str, object]:
    """Load one JSON object and reject other top-level types."""
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError(f"Expected a JSON object in {path}")
    return payload


def nearest_cycle_count(reference_bjd: float, target_bjd: float, period_days: float) -> int:
    """Return the integer cycle count that brings ``reference_bjd`` nearest target."""
    if period_days <= 0.0:
        raise ValueError("period_days must be positive")
    return int(round((target_bjd - reference_bjd) / period_days))


def propagated_markers(
    source: dict[str, float], period_days: float, cycles: int
) -> dict[str, float]:
    """Propagate ingress, midpoint, and egress by an integer number of cycles."""
    return {marker: float(source[marker] + cycles * period_days) for marker in MARKERS}


def reconstruct(
    schedule: dict[str, object],
    ground: dict[str, object],
    ephemeris: dict[str, object],
    exofop: dict[str, object],
) -> dict[str, object]:
    """Build the numerical reconstruction and its explicit inference limits."""
    working = schedule["working_interpretation"]
    if not isinstance(working, dict):
        raise ValueError("working_interpretation must be an object")
    times = working["times"]
    if not isinstance(times, dict):
        raise ValueError("working_interpretation.times must be an object")

    schedule_markers = {
        "ingress": float(times["ingress"]["bjd_tdb"]),  # type: ignore[index]
        "midpoint": float(working["event_midpoint_bjd_tdb"]),
        "egress": float(times["egress"]["bjd_tdb"]),  # type: ignore[index]
    }

    reports = exofop["report_values"]
    if not isinstance(reports, dict):
        raise ValueError("report_values must be an object")
    ulmt = reports["ulmt_2021_10_15"]
    cmo = reports["cmo_2023_05_05"]
    if not isinstance(ulmt, dict) or not isinstance(cmo, dict):
        raise ValueError("ULMT and CMO report values must be objects")
    source_markers = {
        "ingress": float(ulmt["predicted_ingress_bjd_tdb"]),
        "midpoint": float(ulmt["predicted_midpoint_bjd_tdb"]),
        "egress": float(ulmt["predicted_egress_bjd_tdb"]),
    }
    old_period = float(cmo["superseded_period_days"])
    displayed_period = float(cmo["revised_period_days"])
    cycles = nearest_cycle_count(
        source_markers["midpoint"], schedule_markers["midpoint"], old_period
    )
    reconstructed = propagated_markers(source_markers, old_period, cycles)

    comparison_rows = []
    for marker in MARKERS:
        offset_seconds = (schedule_markers[marker] - reconstructed[marker]) * 86400.0
        comparison_rows.append(
            {
                "marker": marker,
                "source_bjd_tdb": source_markers[marker],
                "cycles": cycles,
                "propagated_bjd_tdb": reconstructed[marker],
                "schedule_bjd_tdb": schedule_markers[marker],
                "schedule_minus_propagated_seconds": offset_seconds,
                "absolute_offset_seconds": abs(offset_seconds),
            }
        )

    displayed_midpoint = source_markers["midpoint"] + cycles * displayed_period
    displayed_offset_hours = (
        schedule_markers["midpoint"] - displayed_midpoint
    ) * 24.0

    adopted = ephemeris["ephemeris"]["adopted"]  # type: ignore[index]
    adopted_epoch = float(adopted["epoch_bjd_tdb"])  # type: ignore[index]
    adopted_period = float(adopted["period_days"])  # type: ignore[index]
    adopted_cycle = nearest_cycle_count(
        adopted_epoch, schedule_markers["midpoint"], adopted_period
    )
    adopted_midpoint = adopted_epoch + adopted_cycle * adopted_period

    observation = ground["observation"]
    if not isinstance(observation, dict):
        raise ValueError("ground observation must be an object")
    max_offset = max(float(row["absolute_offset_seconds"]) for row in comparison_rows)

    return {
        "target": "TOI-3505.01",
        "analysis_type": "archival ephemeris consistency reconstruction",
        "source_files": {
            "schedule_check": str(SCHEDULE_CHECK_PATH.relative_to(ROOT)),
            "ground_summary": str(GROUND_SUMMARY_PATH.relative_to(ROOT)),
            "refined_ephemeris": str(EPHEMERIS_PATH.relative_to(ROOT)),
            "exofop_context": str(EXOFOP_PATH.relative_to(ROOT)),
        },
        "source_event": {
            "facility": "ULMT",
            "date_utc": "2021-10-15",
            "markers_bjd_tdb": source_markers,
            "source_url": ulmt["url"],
        },
        "periods": {
            "superseded_period_days": old_period,
            "displayed_revised_period_days": displayed_period,
            "schedule_row_period_days": float(schedule["row"]["Orbital Period"]),  # type: ignore[index]
            "adopted_tess_period_days": adopted_period,
        },
        "reconstruction": {
            "cycles": cycles,
            "marker_comparison": comparison_rows,
            "maximum_absolute_marker_offset_seconds": max_offset,
            "all_markers_agree_within_one_minute": max_offset < 60.0,
            "displayed_period_propagated_midpoint_bjd_tdb": displayed_midpoint,
            "schedule_minus_displayed_period_midpoint_hours": displayed_offset_hours,
        },
        "timeline": {
            "schedule_markers_bjd_tdb": schedule_markers,
            "reconstructed_old_ephemeris_markers_bjd_tdb": reconstructed,
            "observation_start_bjd_tdb": float(observation["start_bjd_tdb"]),
            "observation_end_bjd_tdb": float(observation["end_bjd_tdb"]),
            "adopted_ephemeris_cycle": adopted_cycle,
            "adopted_ephemeris_midpoint_bjd_tdb": adopted_midpoint,
            "schedule_minus_adopted_midpoint_hours": (
                schedule_markers["midpoint"] - adopted_midpoint
            )
            * 24.0,
        },
        "assessment": {
            "strength": "strongly supported archival reconstruction, not proof",
            "supported_statement": (
                "The 2022 schedule ingress, midpoint, and egress are each within "
                "one minute of the public 2021 prediction propagated by 96 cycles "
                "with the superseded 2.9174250-day period."
            ),
            "likely_history": (
                "The timing cells are consistent with having retained a prediction "
                "from the superseded ephemeris while the visible period cell was "
                "later changed to 2.9151488 days."
            ),
            "limits": [
                "The original workbook, formulas, revision history, and Transit Info file are unavailable.",
                "The agreement cannot distinguish a stale formula from manually copied timing cells.",
                "This reconstruction explains the schedule window but does not turn the 2022 null into a transit detection.",
            ],
        },
    }


def write_marker_csv(path: Path, rows: list[dict[str, object]]) -> None:
    """Write the marker comparison with a stable column order."""
    columns = (
        "marker",
        "source_bjd_tdb",
        "cycles",
        "propagated_bjd_tdb",
        "schedule_bjd_tdb",
        "schedule_minus_propagated_seconds",
        "absolute_offset_seconds",
    )
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)


def plot_timeline(result: dict[str, object], path: Path) -> None:
    """Plot broad timing context and a marker-level reconstruction zoom."""
    timeline = result["timeline"]
    reconstruction = result["reconstruction"]
    if not isinstance(timeline, dict) or not isinstance(reconstruction, dict):
        raise ValueError("timeline and reconstruction must be objects")
    schedule = timeline["schedule_markers_bjd_tdb"]
    old = timeline["reconstructed_old_ephemeris_markers_bjd_tdb"]
    if not isinstance(schedule, dict) or not isinstance(old, dict):
        raise ValueError("timeline markers must be objects")
    zero = float(schedule["midpoint"])

    def hours(value: float) -> float:
        return (float(value) - zero) * 24.0

    fig, (context, zoom) = plt.subplots(
        1, 2, figsize=(12.4, 5.0), dpi=240, gridspec_kw={"width_ratios": [1.35, 1.0]}
    )
    colors = {
        "ink": "#1f2933",
        "green": "#006633",
        "gold": "#ffcc33",
        "red": "#b43c35",
        "blue": "#3977a8",
        "muted": "#718096",
    }

    schedule_lo = hours(float(schedule["ingress"]))
    schedule_hi = hours(float(schedule["egress"]))
    obs_lo = hours(float(timeline["observation_start_bjd_tdb"]))
    obs_hi = hours(float(timeline["observation_end_bjd_tdb"]))
    adopted = hours(float(timeline["adopted_ephemeris_midpoint_bjd_tdb"]))
    displayed = -float(reconstruction["schedule_minus_displayed_period_midpoint_hours"])

    context.barh(3, obs_hi - obs_lo, left=obs_lo, height=0.46, color=colors["blue"])
    context.barh(
        2,
        schedule_hi - schedule_lo,
        left=schedule_lo,
        height=0.46,
        color=colors["gold"],
        edgecolor=colors["ink"],
        linewidth=0.8,
    )
    context.plot(adopted, 1, "o", color=colors["green"], markersize=8)
    context.plot(displayed, 0, "D", color=colors["red"], markersize=7)
    context.axvline(0.0, color=colors["ink"], linewidth=0.8, alpha=0.5)
    context.set_yticks(
        [0, 1, 2, 3],
        [
            "2021 epoch + revised P",
            "adopted TESS midpoint",
            "2022 schedule window",
            "GMU observing sequence",
        ],
    )
    context.set_xlabel("Hours from the 2022 schedule midpoint")
    context.set_title("Why the GMU sequence missed the current transit")
    context.grid(axis="x", color="#d7dde3", linewidth=0.7)
    context.set_axisbelow(True)
    context.annotate(
        f"{abs(adopted):.1f} h earlier",
        xy=(adopted, 1),
        xytext=(adopted + 1.2, 1.35),
        arrowprops={"arrowstyle": "-", "color": colors["green"]},
        color=colors["green"],
        fontsize=9,
    )
    context.annotate(
        f"{abs(displayed):.2f} h early",
        xy=(displayed, 0),
        xytext=(displayed - 4.6, 0.38),
        arrowprops={"arrowstyle": "-", "color": colors["red"]},
        color=colors["red"],
        fontsize=9,
    )

    marker_y = {"ingress": 2, "midpoint": 1, "egress": 0}
    rows = reconstruction["marker_comparison"]
    if not isinstance(rows, list):
        raise ValueError("marker_comparison must be a list")
    for row in rows:
        marker = str(row["marker"])
        y = marker_y[marker]
        propagated_seconds = (
            float(row["propagated_bjd_tdb"]) - float(schedule[marker])
        ) * 86400.0
        zoom.plot(0.0, y, "o", color=colors["gold"], markersize=8, label=None)
        zoom.plot(
            propagated_seconds,
            y,
            "|",
            color=colors["red"],
            markersize=16,
            markeredgewidth=2.2,
        )
        zoom.plot([0.0, propagated_seconds], [y, y], color=colors["muted"], linewidth=1.1)
        zoom.text(
            propagated_seconds,
            y - 0.18 if marker == "ingress" else y + 0.18,
            f"{propagated_seconds:+.0f} s",
            ha="center",
            va="top" if marker == "ingress" else "bottom",
            fontsize=9,
            color=colors["red"],
        )
    zoom.axvline(0.0, color=colors["ink"], linewidth=0.8, alpha=0.55)
    zoom.set_yticks([0, 1, 2], ["Egress", "Midpoint", "Ingress"])
    zoom.set_xlabel("Propagated old marker minus schedule marker (seconds)")
    zoom.set_title("Old ephemeris reproduces all three markers", pad=16)
    zoom.xaxis.set_major_locator(MultipleLocator(20))
    zoom.grid(axis="x", color="#d7dde3", linewidth=0.7)
    zoom.set_axisbelow(True)
    zoom.set_xlim(-45, 75)
    zoom.set_ylim(-0.55, 2.45)
    zoom.text(
        0.02,
        0.035,
        "Gold circles: schedule\nRed ticks: 2021 event + 96 x 2.9174250 d",
        transform=zoom.transAxes,
        fontsize=8.5,
        color=colors["muted"],
        va="bottom",
    )

    for axis in (context, zoom):
        axis.spines[["top", "right"]].set_visible(False)
        axis.tick_params(colors=colors["ink"])
    fig.suptitle(
        "TOI-3505.01: archival reconstruction of the 2022 schedule",
        fontsize=15,
        color=colors["ink"],
        y=0.995,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(path, facecolor="white", bbox_inches="tight")
    fig.savefig(path.with_suffix(".svg"), facecolor="white", bbox_inches="tight")
    plt.close(fig)


def write_readme(result: dict[str, object], path: Path) -> None:
    """Write a concise human-facing interpretation beside the machine record."""
    reconstruction = result["reconstruction"]
    timeline = result["timeline"]
    assert isinstance(reconstruction, dict) and isinstance(timeline, dict)
    path.write_text(
        f"""# TOI-3505.01 schedule reconstruction

The recovered 2022 schedule window is consistent with a superseded public
ephemeris. Propagating the ExoFOP ULMT prediction from 2021-10-15 by
`{int(reconstruction['cycles'])}` cycles at `2.9174250` days reproduces the
schedule ingress, midpoint, and egress to within
`{float(reconstruction['maximum_absolute_marker_offset_seconds']):.1f}` seconds.

The schedule row visibly contains the later `2.9151488`-day period. Applying
that period to the same 2021 midpoint instead predicts an event
`{float(reconstruction['schedule_minus_displayed_period_midpoint_hours']):.2f}`
hours before the schedule midpoint. The adopted four-sector TESS ephemeris
places the nearest event
`{float(timeline['schedule_minus_adopted_midpoint_hours']):.2f}` hours before
the schedule midpoint.

## Interpretation

The timing cells are therefore strongly consistent with a prediction retained
from the older ephemeris while the displayed period was updated. This is an
archival reconstruction, not proof: the original workbook, formulas, revision
history, and Transit Info file remain unavailable. The result explains why the
2022 observing window was stale; it does not convert the ground-based null into
a transit detection.

## Products

- `schedule_reconstruction.json` - complete inputs, calculations, and limits.
- `marker_comparison.csv` - marker-by-marker propagation residuals.
- `01_schedule_reconstruction.png` and `.svg` - publication figure.

Regenerate with:

```bash
.venv/bin/python src/reconstruct_toi3505_schedule.py
```
""",
        encoding="utf-8",
    )


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    result = reconstruct(
        load_object(SCHEDULE_CHECK_PATH),
        load_object(GROUND_SUMMARY_PATH),
        load_object(EPHEMERIS_PATH),
        load_object(EXOFOP_PATH),
    )
    (OUTPUT_DIR / "schedule_reconstruction.json").write_text(
        json.dumps(result, indent=2) + "\n", encoding="utf-8"
    )
    reconstruction = result["reconstruction"]
    assert isinstance(reconstruction, dict)
    rows = reconstruction["marker_comparison"]
    assert isinstance(rows, list)
    write_marker_csv(OUTPUT_DIR / "marker_comparison.csv", rows)
    plot_timeline(result, OUTPUT_DIR / "01_schedule_reconstruction.png")
    write_readme(result, OUTPUT_DIR / "README.md")

    print("TOI-3505.01 schedule reconstruction")
    print(f"  cycles: {reconstruction['cycles']}")
    print(
        "  maximum marker offset: "
        f"{float(reconstruction['maximum_absolute_marker_offset_seconds']):.1f} s"
    )
    print(
        "  displayed-period mismatch: "
        f"{float(reconstruction['schedule_minus_displayed_period_midpoint_hours']):.2f} h"
    )
    print(f"Wrote {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
