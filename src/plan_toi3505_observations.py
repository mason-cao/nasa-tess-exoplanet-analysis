"""Plan observable TOI-3505.01 transit windows from the adopted ephemeris.

The planner propagates the full epoch-period covariance, converts BJD_TDB
midpoints to UTC and local civil time, and evaluates the complete observing
sequence at the GMU 0.8 m site.  Its outputs are planning aids, not approved
Transit Info files; a mentor or observer must confirm a window before booking
the telescope.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from dataclasses import dataclass
from datetime import UTC, date, datetime, time, timedelta
from pathlib import Path
from zoneinfo import ZoneInfo

import numpy as np
from astropy import units as u
from astropy.coordinates import AltAz, get_body, get_sun
from astropy.time import Time
from astropy.utils import iers

from toi3505_schedule import OBSERVATORY, TARGET_COORD, WORKING_TIMEZONE


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_EPHEMERIS_PATH = (
    ROOT / "outputs" / "toi3505_ephemeris_refined" / "ephemeris_refined.json"
)
DEFAULT_VALIDATION_PATH = (
    ROOT / "outputs" / "toi3505_data_validation" / "analysis_summary.json"
)
DEFAULT_OUTPUT_DIR = ROOT / "outputs" / "toi3505_observation_plan"
DEFAULT_PLAN_DAYS = 90
DEFAULT_BASELINE_HOURS = 1.0
DEFAULT_MIN_TARGET_ALTITUDE_DEG = 25.0
DEFAULT_MAX_SUN_ALTITUDE_DEG = -18.0
GRID_SAMPLES = 49


@dataclass(frozen=True)
class PlanningLimits:
    baseline_hours: float = DEFAULT_BASELINE_HOURS
    min_target_altitude_deg: float = DEFAULT_MIN_TARGET_ALTITUDE_DEG
    max_sun_altitude_deg: float = DEFAULT_MAX_SUN_ALTITUDE_DEG

    def validate(self) -> None:
        if self.baseline_hours < 0:
            raise ValueError("Baseline hours must be non-negative")
        if not -90.0 <= self.min_target_altitude_deg <= 90.0:
            raise ValueError("Minimum target altitude must be between -90 and 90")
        if not -90.0 <= self.max_sun_altitude_deg <= 90.0:
            raise ValueError("Maximum Sun altitude must be between -90 and 90")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ephemeris", type=Path, default=DEFAULT_EPHEMERIS_PATH)
    parser.add_argument("--validation", type=Path, default=DEFAULT_VALIDATION_PATH)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument(
        "--start",
        type=date.fromisoformat,
        help="First local calendar date to consider (YYYY-MM-DD; default: today)",
    )
    parser.add_argument(
        "--end",
        type=date.fromisoformat,
        help="Exclusive local end date (YYYY-MM-DD; default: start + 90 days)",
    )
    parser.add_argument(
        "--days",
        type=int,
        default=DEFAULT_PLAN_DAYS,
        help="Planning horizon when --end is omitted (default: 90)",
    )
    parser.add_argument("--timezone", default=WORKING_TIMEZONE)
    parser.add_argument(
        "--duration-hours",
        type=float,
        help="Override the official SPOC multi-sector duration",
    )
    parser.add_argument("--baseline-hours", type=float, default=DEFAULT_BASELINE_HOURS)
    parser.add_argument(
        "--min-target-altitude-deg",
        type=float,
        default=DEFAULT_MIN_TARGET_ALTITUDE_DEG,
    )
    parser.add_argument(
        "--max-sun-altitude-deg",
        type=float,
        default=DEFAULT_MAX_SUN_ALTITUDE_DEG,
    )
    return parser.parse_args()


def configure_astropy() -> None:
    """Keep planning deterministic and usable without a network connection."""
    iers.conf.auto_download = False
    iers.conf.auto_max_age = None


def read_json(path: Path) -> dict[str, object]:
    if not path.exists():
        raise FileNotFoundError(path)
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError(f"Expected a JSON object in {path}")
    return payload


def load_adopted_ephemeris(path: Path) -> dict[str, object]:
    payload = read_json(path)
    try:
        adopted = payload["ephemeris"]["adopted"]  # type: ignore[index]
    except (KeyError, TypeError) as error:
        raise ValueError(f"No adopted ephemeris in {path}") from error
    if not isinstance(adopted, dict):
        raise ValueError(f"Invalid adopted ephemeris in {path}")

    required = (
        "events",
        "epoch_bjd_tdb",
        "epoch_error_days",
        "period_days",
        "period_error_days",
        "covariance_days2",
    )
    missing = [name for name in required if name not in adopted]
    if missing:
        raise ValueError(f"Adopted ephemeris is missing: {', '.join(missing)}")
    if float(adopted["period_days"]) <= 0:
        raise ValueError("Adopted period must be positive")
    return adopted


def load_spoc_duration_hours(path: Path) -> float:
    payload = read_json(path)
    records = payload.get("official_multisector_tce")
    if not isinstance(records, list) or len(records) != 1:
        raise ValueError(f"Expected one official multi-sector TCE in {path}")
    record = records[0]
    if not isinstance(record, dict) or "fit_duration_hours" not in record:
        raise ValueError(f"Official duration is missing from {path}")
    duration = float(record["fit_duration_hours"])
    if not 0 < duration < 24:
        raise ValueError(f"Implausible transit duration: {duration}")
    return duration


def prediction_uncertainty_minutes(ephemeris: dict[str, object], cycle: int) -> float:
    epoch_error = float(ephemeris["epoch_error_days"])
    period_error = float(ephemeris["period_error_days"])
    covariance = float(ephemeris["covariance_days2"])
    variance = epoch_error**2 + (cycle * period_error) ** 2 + 2.0 * cycle * covariance
    if variance < -1e-15:
        raise ValueError(
            f"Negative prediction variance at cycle {cycle}: {variance:.6g}"
        )
    return math.sqrt(max(0.0, variance)) * 1440.0


def bjd_tdb_to_utc(midpoint_bjd_tdb: float) -> Time:
    """Invert the target-specific barycentric correction by fixed iteration."""
    instant = Time(midpoint_bjd_tdb, format="jd", scale="tdb", location=OBSERVATORY)
    for _ in range(4):
        correction = instant.light_travel_time(TARGET_COORD, kind="barycentric")
        instant = Time(
            midpoint_bjd_tdb - correction.to_value(u.day),
            format="jd",
            scale="tdb",
            location=OBSERVATORY,
        )
    return instant.utc


def local_interval(
    start_date: date, end_date: date, timezone: ZoneInfo
) -> tuple[datetime, datetime]:
    if end_date <= start_date:
        raise ValueError("The planning end date must be after the start date")
    start = datetime.combine(start_date, time.min, tzinfo=timezone)
    end = datetime.combine(end_date, time.min, tzinfo=timezone)
    return start, end


def iso_utc(instant: Time) -> str:
    value = instant.to_datetime(timezone=UTC)
    return value.isoformat(timespec="seconds").replace("+00:00", "Z")


def iso_local(instant: Time, timezone: ZoneInfo) -> str:
    return instant.to_datetime(timezone=timezone).isoformat(timespec="seconds")


def display_local(value: str, timezone: ZoneInfo) -> str:
    instant = datetime.fromisoformat(value).astimezone(timezone)
    if instant.second >= 30:
        instant += timedelta(minutes=1)
    return instant.replace(second=0, microsecond=0).strftime("%Y-%m-%d %H:%M %Z")


def portable_path(path: Path) -> str:
    resolved = path.resolve()
    try:
        return str(resolved.relative_to(ROOT))
    except ValueError:
        return str(resolved)


def event_geometry(
    midpoint: Time,
    duration_hours: float,
    limits: PlanningLimits,
    timezone: ZoneInfo,
) -> dict[str, object]:
    half_duration = duration_hours / 2.0 * u.hour
    half_sequence = (duration_hours / 2.0 + limits.baseline_hours) * u.hour
    sequence_start = midpoint - half_sequence
    ingress = midpoint - half_duration
    egress = midpoint + half_duration
    sequence_end = midpoint + half_sequence

    grid_jd = np.linspace(sequence_start.jd, sequence_end.jd, GRID_SAMPLES)
    grid = Time(
        grid_jd,
        format="jd",
        scale="utc",
        location=OBSERVATORY,
    )
    altaz = AltAz(obstime=grid, location=OBSERVATORY)
    target_altitudes = TARGET_COORD.transform_to(altaz).alt.deg
    sun_altitudes = get_sun(grid).transform_to(altaz).alt.deg
    moon_altitudes = get_body("moon", grid, OBSERVATORY).transform_to(altaz).alt.deg

    endpoint_grid = Time(
        [sequence_start.jd, midpoint.jd, sequence_end.jd],
        format="jd",
        scale="utc",
        location=OBSERVATORY,
    )
    endpoint_altaz = AltAz(obstime=endpoint_grid, location=OBSERVATORY)
    endpoint_target = TARGET_COORD.transform_to(endpoint_altaz).alt.deg

    minimum_target = float(np.min(target_altitudes))
    maximum_target = float(np.max(target_altitudes))
    maximum_sun = float(np.max(sun_altitudes))
    maximum_moon = float(np.max(moon_altitudes))
    reasons: list[str] = []
    if minimum_target < limits.min_target_altitude_deg:
        reasons.append(f"target below {limits.min_target_altitude_deg:g} deg")
    if maximum_sun > limits.max_sun_altitude_deg:
        reasons.append(f"Sun above {limits.max_sun_altitude_deg:g} deg")

    return {
        "sequence_start_utc": iso_utc(sequence_start),
        "ingress_utc": iso_utc(ingress),
        "midpoint_utc": iso_utc(midpoint),
        "egress_utc": iso_utc(egress),
        "sequence_end_utc": iso_utc(sequence_end),
        "sequence_start_local": iso_local(sequence_start, timezone),
        "ingress_local": iso_local(ingress, timezone),
        "midpoint_local": iso_local(midpoint, timezone),
        "egress_local": iso_local(egress, timezone),
        "sequence_end_local": iso_local(sequence_end, timezone),
        "target_altitude_start_deg": float(endpoint_target[0]),
        "target_altitude_midpoint_deg": float(endpoint_target[1]),
        "target_altitude_end_deg": float(endpoint_target[2]),
        "minimum_target_altitude_deg": minimum_target,
        "maximum_target_altitude_deg": maximum_target,
        "maximum_sun_altitude_deg": maximum_sun,
        "maximum_moon_altitude_deg": maximum_moon,
        "moon_below_horizon_throughout": bool(maximum_moon < 0.0),
        "observable": not reasons,
        "rejection_reasons": "; ".join(reasons),
    }


def plan_events(
    ephemeris: dict[str, object],
    start_date: date,
    end_date: date,
    duration_hours: float,
    limits: PlanningLimits,
    timezone_name: str = WORKING_TIMEZONE,
) -> list[dict[str, object]]:
    configure_astropy()
    limits.validate()
    if not 0 < duration_hours < 24:
        raise ValueError("Transit duration must be between zero and 24 hours")
    timezone = ZoneInfo(timezone_name)
    start_local, end_local = local_interval(start_date, end_date, timezone)
    start_utc = start_local.astimezone(UTC)
    end_utc = end_local.astimezone(UTC)

    epoch = float(ephemeris["epoch_bjd_tdb"])
    period = float(ephemeris["period_days"])
    start_tdb = float(Time(start_utc).tdb.jd)
    end_tdb = float(Time(end_utc).tdb.jd)
    first_cycle = math.floor((start_tdb - epoch) / period) - 2
    last_cycle = math.ceil((end_tdb - epoch) / period) + 2

    events: list[dict[str, object]] = []
    for cycle in range(first_cycle, last_cycle + 1):
        midpoint_bjd = epoch + cycle * period
        midpoint = bjd_tdb_to_utc(midpoint_bjd)
        midpoint_datetime = midpoint.to_datetime(timezone=UTC)
        if not start_utc <= midpoint_datetime < end_utc:
            continue
        event = {
            "cycle": cycle,
            "midpoint_bjd_tdb": midpoint_bjd,
            "midpoint_uncertainty_minutes": prediction_uncertainty_minutes(
                ephemeris, cycle
            ),
            **event_geometry(midpoint, duration_hours, limits, timezone),
        }
        events.append(event)
    return events


CSV_FIELDS = (
    "cycle",
    "midpoint_bjd_tdb",
    "midpoint_uncertainty_minutes",
    "sequence_start_utc",
    "ingress_utc",
    "midpoint_utc",
    "egress_utc",
    "sequence_end_utc",
    "sequence_start_local",
    "ingress_local",
    "midpoint_local",
    "egress_local",
    "sequence_end_local",
    "target_altitude_start_deg",
    "target_altitude_midpoint_deg",
    "target_altitude_end_deg",
    "minimum_target_altitude_deg",
    "maximum_target_altitude_deg",
    "maximum_sun_altitude_deg",
    "maximum_moon_altitude_deg",
    "moon_below_horizon_throughout",
    "observable",
    "rejection_reasons",
)


def ics_escape(value: str) -> str:
    return (
        value.replace("\\", "\\\\")
        .replace(";", "\\;")
        .replace(",", "\\,")
        .replace("\n", "\\n")
    )


def ics_datetime(value: str) -> str:
    instant = datetime.fromisoformat(value.replace("Z", "+00:00"))
    return instant.astimezone(UTC).strftime("%Y%m%dT%H%M%SZ")


def fold_ics_line(line: str) -> list[str]:
    """Fold the ASCII calendar lines at RFC 5545's 75-octet boundary."""
    if not line.isascii():
        raise ValueError("Calendar output must remain ASCII")
    if len(line) <= 75:
        return [line]
    folded = [line[:75]]
    remainder = line[75:]
    while remainder:
        folded.append(" " + remainder[:74])
        remainder = remainder[74:]
    return folded


def build_calendar(
    events: list[dict[str, object]], generated_at: datetime | None = None
) -> str:
    generated = (generated_at or datetime.now(UTC)).astimezone(UTC)
    stamp = generated.strftime("%Y%m%dT%H%M%SZ")
    lines = [
        "BEGIN:VCALENDAR",
        "VERSION:2.0",
        "PRODID:-//TOI-3505.01 analysis//Observation planner//EN",
        "CALSCALE:GREGORIAN",
        "METHOD:PUBLISH",
    ]
    for event in events:
        if not bool(event["observable"]):
            continue
        start = str(event["sequence_start_utc"])
        end = str(event["sequence_end_utc"])
        description = (
            f"Predicted midpoint {event['midpoint_local']}; "
            f"1-sigma uncertainty {float(event['midpoint_uncertainty_minutes']):.2f} "
            "minutes. Planning output only; independently confirm before observing."
        )
        lines.extend(
            [
                "BEGIN:VEVENT",
                f"UID:toi3505-cycle-{event['cycle']}@nasa-tess-exoplanet-analysis",
                f"DTSTAMP:{stamp}",
                f"DTSTART:{ics_datetime(start)}",
                f"DTEND:{ics_datetime(end)}",
                "SUMMARY:TOI-3505.01 transit observation window",
                f"DESCRIPTION:{ics_escape(description)}",
                "END:VEVENT",
            ]
        )
    lines.append("END:VCALENDAR")
    folded = [part for line in lines for part in fold_ics_line(line)]
    return "\r\n".join(folded) + "\r\n"


def build_readme(
    summary: dict[str, object], timezone: ZoneInfo, ephemeris_path: Path
) -> str:
    events = summary["events"]
    assert isinstance(events, list)
    observable = [event for event in events if bool(event["observable"])]
    limits = summary["planning_limits"]
    assert isinstance(limits, dict)
    rows = []
    for event in observable:
        rows.append(
            "| {cycle} | {start} | {midpoint} | {end} | {uncertainty:.2f} | "
            "{minimum:.1f} | {moon} |".format(
                cycle=event["cycle"],
                start=display_local(str(event["sequence_start_local"]), timezone),
                midpoint=display_local(str(event["midpoint_local"]), timezone),
                end=display_local(str(event["sequence_end_local"]), timezone),
                uncertainty=float(event["midpoint_uncertainty_minutes"]),
                minimum=float(event["minimum_target_altitude_deg"]),
                moon="yes" if event["moon_below_horizon_throughout"] else "no",
            )
        )
    if not rows:
        rows.append(
            "| — | No complete windows pass the planning limits | — | — | — | — | — |"
        )

    source = portable_path(ephemeris_path)
    return f"""# TOI-3505.01 observation plan

Planning windows generated from
`{source}`. These times are not an approved Transit Info file; independently
confirm the ephemeris, observatory clock, and telescope availability before
observing.

## Observable windows

| Cycle | Sequence start | Midpoint | Sequence end | 1σ (min) | Min altitude | Moon below throughout |
| ---: | --- | --- | --- | ---: | ---: | :---: |
{chr(10).join(rows)}

The sequence includes {float(limits['baseline_hours']):g} hour of baseline on each side of a
{float(summary['duration_hours']):.3f}-hour transit. A complete window must keep the target
above {float(limits['min_target_altitude_deg']):g} degrees and the Sun below
{float(limits['max_sun_altitude_deg']):g} degrees.

## Files

- `transit_windows.csv` — every predicted midpoint in the requested range,
  including rejected windows and reasons.
- `observation_plan.json` — inputs, thresholds, adopted ephemeris, and events.
- `observable_transits.ics` — calendar events for windows that pass the limits.

Regenerate with:

```bash
.venv/bin/python src/plan_toi3505_observations.py
```
"""


def write_outputs(
    output_dir: Path,
    events: list[dict[str, object]],
    ephemeris: dict[str, object],
    ephemeris_path: Path,
    validation_path: Path,
    start_date: date,
    end_date: date,
    duration_hours: float,
    limits: PlanningLimits,
    timezone_name: str,
) -> dict[str, object]:
    output_dir.mkdir(parents=True, exist_ok=True)
    generated_at = datetime.now(UTC)
    with (output_dir / "transit_windows.csv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=CSV_FIELDS)
        writer.writeheader()
        writer.writerows(events)

    summary: dict[str, object] = {
        "target": "TOI-3505.01",
        "generated_utc": generated_at.isoformat(),
        "observatory": "GMU 0.8 m",
        "timezone": timezone_name,
        "date_range": {
            "start_local": str(start_date),
            "end_local_exclusive": str(end_date),
        },
        "duration_hours": duration_hours,
        "planning_limits": {
            "baseline_hours": limits.baseline_hours,
            "min_target_altitude_deg": limits.min_target_altitude_deg,
            "max_sun_altitude_deg": limits.max_sun_altitude_deg,
        },
        "source_files": {
            "ephemeris": portable_path(ephemeris_path),
            "duration": portable_path(validation_path),
        },
        "adopted_ephemeris": dict(ephemeris),
        "event_count": len(events),
        "observable_event_count": sum(bool(event["observable"]) for event in events),
        "events": events,
        "approval_status": "planning output; mentor confirmation required",
    }
    (output_dir / "observation_plan.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    timezone = ZoneInfo(timezone_name)
    (output_dir / "README.md").write_text(
        build_readme(summary, timezone, ephemeris_path), encoding="utf-8"
    )
    (output_dir / "observable_transits.ics").write_bytes(
        build_calendar(events, generated_at).encode("utf-8")
    )
    return summary


def main() -> None:
    args = parse_args()
    if args.days <= 0:
        raise ValueError("--days must be positive")
    timezone = ZoneInfo(args.timezone)
    start_date = args.start or datetime.now(timezone).date()
    end_date = args.end or start_date + timedelta(days=args.days)
    ephemeris_path = args.ephemeris.resolve()
    validation_path = args.validation.resolve()
    ephemeris = load_adopted_ephemeris(ephemeris_path)
    duration_hours = (
        args.duration_hours
        if args.duration_hours is not None
        else load_spoc_duration_hours(validation_path)
    )
    limits = PlanningLimits(
        baseline_hours=args.baseline_hours,
        min_target_altitude_deg=args.min_target_altitude_deg,
        max_sun_altitude_deg=args.max_sun_altitude_deg,
    )
    events = plan_events(
        ephemeris,
        start_date,
        end_date,
        duration_hours,
        limits,
        args.timezone,
    )
    summary = write_outputs(
        args.output_dir.resolve(),
        events,
        ephemeris,
        ephemeris_path,
        validation_path,
        start_date,
        end_date,
        duration_hours,
        limits,
        args.timezone,
    )
    print(
        f"Planned {summary['event_count']} transits; "
        f"{summary['observable_event_count']} pass the complete-window limits."
    )
    print(f"Wrote {args.output_dir.resolve()}")


if __name__ == "__main__":
    main()
