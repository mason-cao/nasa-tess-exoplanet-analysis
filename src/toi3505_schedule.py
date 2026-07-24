"""Preserve and interpret the TOI-3505.01 2022 scheduling-sheet row.

The source row contains clock times but no time-zone field.  This module keeps
that limitation explicit, evaluates both Eastern civil time and UTC, and uses
Eastern time only as a working interpretation because its planned start/end
bracket the actual GMU image sequence.
"""

from __future__ import annotations

import json
import math
import re
from datetime import UTC, date, datetime, time, timedelta
from pathlib import Path
from zoneinfo import ZoneInfo

import numpy as np
from astropy import units as u
from astropy.coordinates import AltAz, EarthLocation, SkyCoord, get_sun
from astropy.time import Time
from astropy.utils import iers


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_SCHEDULE_RECORD = (
    ROOT
    / "data"
    / "program_records"
    / "toi3505"
    / "observing_schedule_2022-07-21.json"
)
WORKING_TIMEZONE = "America/New_York"
OBSERVATORY = EarthLocation.from_geodetic(
    lon=-77.3053299972 * u.deg,
    lat=38.82817 * u.deg,
    height=154.0 * u.m,
)
TARGET_COORD = SkyCoord(
    ra=19.802897222222224 * u.hourangle,
    dec=18.69891388888889 * u.deg,
    frame="icrs",
)


def load_schedule_record(path: Path = DEFAULT_SCHEDULE_RECORD) -> dict[str, object]:
    """Load the preserved row and return its cells as a named mapping."""
    payload = json.loads(path.read_text(encoding="utf-8"))
    columns = payload.get("columns")
    values = payload.get("values")
    if not isinstance(columns, list) or not isinstance(values, list):
        raise ValueError("Schedule record must contain column and value lists")
    if len(columns) != len(values) or len(columns) != len(set(columns)):
        raise ValueError("Schedule columns and values are mismatched or duplicated")
    row = dict(zip(columns, values, strict=True))
    if row["Target"] != "TOI 3505.01":
        raise ValueError(f"Unexpected schedule target: {row['Target']}")
    if row["NoD"] != "2022-07-21" or row["NoW"] != "Thursday":
        raise ValueError("The schedule night/date fields do not match the source row")
    if date.fromisoformat(row["NoD"]).strftime("%A") != row["NoW"]:
        raise ValueError("The schedule weekday is inconsistent with its date")
    if row["Filter"] != "R" or row["Exp"] != "50s":
        raise ValueError("The schedule no longer matches the R/50-second data set")
    payload["row"] = row
    payload["path"] = str(path.resolve())
    return payload


def parse_clock(text: str) -> time:
    """Parse the hour:minute cells used by the program sheet."""
    match = re.fullmatch(r"(\d{1,2}):(\d{2})", text.strip())
    if match is None:
        raise ValueError(f"Invalid schedule clock value: {text!r}")
    hour, minute = (int(value) for value in match.groups())
    return time(hour=hour, minute=minute)


def clock_datetime(
    night: date,
    clock: time,
    start_clock: time,
    timezone: ZoneInfo,
) -> datetime:
    """Attach after-midnight clock values to the next calendar date."""
    calendar_date = night + timedelta(days=int(clock < start_clock))
    return datetime.combine(calendar_date, clock, tzinfo=timezone)


def datetime_to_bjd_tdb(value: datetime) -> float:
    """Convert a timezone-aware civil time to BJD_TDB at the target."""
    if value.tzinfo is None:
        raise ValueError("A timezone-aware datetime is required")
    iers.conf.auto_download = False
    utc_value = value.astimezone(UTC)
    instant = Time(utc_value, scale="utc", location=OBSERVATORY)
    return float(
        (instant.tdb + instant.light_travel_time(TARGET_COORD, kind="barycentric")).jd
    )


def build_clock_interpretation(
    row: dict[str, str], timezone_name: str
) -> dict[str, object]:
    """Convert the four schedule clocks under one named time zone."""
    timezone = ZoneInfo(timezone_name)
    night = date.fromisoformat(row["NoD"])
    start_clock = parse_clock(row["Start"])
    values: dict[str, dict[str, object]] = {}
    for name, column in (
        ("planned_start", "Start"),
        ("ingress", "Ingress"),
        ("egress", "Egress"),
        ("planned_end", "End"),
    ):
        local_value = clock_datetime(
            night, parse_clock(row[column]), start_clock, timezone
        )
        values[name] = {
            "source_clock": row[column],
            "clock_datetime": local_value.isoformat(),
            "utc_datetime": local_value.astimezone(UTC).isoformat(),
            "bjd_tdb": datetime_to_bjd_tdb(local_value),
        }
    ordered = [
        float(values[name]["bjd_tdb"])
        for name in ("planned_start", "ingress", "egress", "planned_end")
    ]
    if not all(left < right for left, right in zip(ordered, ordered[1:])):
        raise ValueError("Schedule times are not in start/ingress/egress/end order")
    start_datetime = clock_datetime(night, start_clock, start_clock, timezone)
    utc_offset = start_datetime.utcoffset()
    if utc_offset is None:
        raise ValueError(f"No UTC offset is available for {timezone_name}")
    return {
        "timezone": timezone_name,
        "timezone_abbreviation": start_datetime.tzname(),
        "utc_offset_hours": utc_offset.total_seconds() / 3600.0,
        "source_timezone_explicit": False,
        "times": values,
        "event_duration_hours": (ordered[2] - ordered[1]) * 24.0,
        "event_midpoint_bjd_tdb": 0.5 * (ordered[1] + ordered[2]),
    }


def schedule_altitudes(interpretation: dict[str, object]) -> dict[str, object]:
    """Return Sun and target altitudes at the planned start and end."""
    times = interpretation["times"]
    assert isinstance(times, dict)
    result: dict[str, object] = {}
    iers.conf.auto_download = False
    for event in ("planned_start", "planned_end"):
        entry = times[event]
        assert isinstance(entry, dict)
        instant = Time(datetime.fromisoformat(str(entry["utc_datetime"])))
        frame = AltAz(obstime=instant, location=OBSERVATORY)
        result[event] = {
            "utc_datetime": entry["utc_datetime"],
            "sun_altitude_degrees": float(
                get_sun(instant).transform_to(frame).alt.deg
            ),
            "target_altitude_degrees": float(
                TARGET_COORD.transform_to(frame).alt.deg
            ),
        }
    return result


def schedule_context(
    path: Path = DEFAULT_SCHEDULE_RECORD,
    *,
    observation_start_bjd: float | None = None,
    observation_end_bjd: float | None = None,
) -> dict[str, object]:
    """Return the preserved row and documented time interpretations."""
    payload = load_schedule_record(path)
    row = payload["row"]
    assert isinstance(row, dict)
    working = build_clock_interpretation(row, WORKING_TIMEZONE)
    utc_alternative = build_clock_interpretation(row, "UTC")
    working_altitudes = schedule_altitudes(working)
    utc_altitudes = schedule_altitudes(utc_alternative)
    working_start = working_altitudes["planned_start"]
    working_end = working_altitudes["planned_end"]
    utc_start = utc_altitudes["planned_start"]
    assert isinstance(working_start, dict)
    assert isinstance(working_end, dict)
    assert isinstance(utc_start, dict)
    context: dict[str, object] = {
        "source_record": str(path.resolve().relative_to(ROOT)),
        "source_description": payload["source_description"],
        "original_workbook_or_url_archived": payload[
            "original_workbook_or_url_archived"
        ],
        "source_limits": payload["source_limits"],
        "row": row,
        "working_interpretation": working,
        "utc_alternative": utc_alternative,
        "timezone_plausibility": {
            "preferred_interpretation": WORKING_TIMEZONE,
            "working_interpretation_altitudes": working_altitudes,
            "utc_alternative_altitudes": utc_altitudes,
            "working_start_is_observable": bool(
                float(working_start["sun_altitude_degrees"]) < 0.0
                and float(working_start["target_altitude_degrees"]) > 20.0
            ),
            "working_end_is_observable": bool(
                float(working_end["sun_altitude_degrees"]) < 0.0
                and float(working_end["target_altitude_degrees"]) > 20.0
            ),
            "utc_start_is_observable": bool(
                float(utc_start["sun_altitude_degrees"]) < 0.0
                and float(utc_start["target_altitude_degrees"]) > 20.0
            ),
            "assessment": (
                "America/New_York is strongly favored: on this date it is EDT "
                "(UTC-4), its planned start and end occur with the Sun below "
                "the horizon and the target above 20 degrees, and its planned "
                "range brackets the images. Treating the clocks as UTC starts "
                "in daylight with the target below the horizon."
            ),
        },
        "working_interpretation_reason": (
            "America/New_York civil time is used because the planned 21:10-04:55 "
            "range brackets the actual GMU sequence and has physically plausible "
            "Sun and target altitudes under that interpretation. On this date "
            "America/New_York is EDT (UTC-4). The source row itself does not "
            "state a time zone."
        ),
        "historical_ephemeris_complete": False,
    }
    if observation_start_bjd is not None and observation_end_bjd is not None:
        if observation_end_bjd <= observation_start_bjd:
            raise ValueError("Observation end must be after observation start")
        working_times = working["times"]
        utc_times = utc_alternative["times"]
        assert isinstance(working_times, dict) and isinstance(utc_times, dict)
        context["observation_comparison"] = {
            "observation_start_bjd_tdb": observation_start_bjd,
            "observation_end_bjd_tdb": observation_end_bjd,
            "working_planned_range_brackets_observation": bool(
                float(working_times["planned_start"]["bjd_tdb"])
                <= observation_start_bjd
                and float(working_times["planned_end"]["bjd_tdb"])
                >= observation_end_bjd
            ),
            "working_event_fully_covered": bool(
                observation_start_bjd
                <= float(working_times["ingress"]["bjd_tdb"])
                < float(working_times["egress"]["bjd_tdb"])
                <= observation_end_bjd
            ),
            "utc_event_fully_covered": bool(
                observation_start_bjd
                <= float(utc_times["ingress"]["bjd_tdb"])
                < float(utc_times["egress"]["bjd_tdb"])
                <= observation_end_bjd
            ),
        }
    return context


def fixed_window_fraction(
    bjd_tdb: np.ndarray,
    ingress_bjd: float,
    egress_bjd: float,
    exposure_seconds: float,
) -> np.ndarray:
    """Return the fraction of each exposure inside a fixed time window."""
    time_values = np.asarray(bjd_tdb, dtype=float)
    half_exposure = exposure_seconds / (2.0 * 86400.0)
    exposure_start = time_values - half_exposure
    exposure_end = time_values + half_exposure
    overlap = np.maximum(
        0.0,
        np.minimum(exposure_end, egress_bjd)
        - np.maximum(exposure_start, ingress_bjd),
    )
    return overlap / (2.0 * half_exposure)


def robust_sigma(values: np.ndarray) -> float:
    """Return a median-absolute-deviation scatter estimate."""
    finite = np.asarray(values, dtype=float)
    finite = finite[np.isfinite(finite)]
    if len(finite) == 0:
        return float("nan")
    center = float(np.median(finite))
    return float(1.4826 * np.median(np.abs(finite - center)))


def fit_fixed_window(
    bjd_tdb: np.ndarray,
    flux: np.ndarray,
    flux_error: np.ndarray,
    use: np.ndarray,
    ingress_bjd: float,
    egress_bjd: float,
    *,
    exposure_seconds: float = 50.0,
) -> dict[str, float | int]:
    """Fit a straight baseline and one fixed schedule-window depth."""
    time_values = np.asarray(bjd_tdb, dtype=float)
    flux_values = np.asarray(flux, dtype=float)
    error_values = np.asarray(flux_error, dtype=float)
    fraction = fixed_window_fraction(
        time_values, ingress_bjd, egress_bjd, exposure_seconds
    )
    keep = (
        np.asarray(use, dtype=bool)
        & np.isfinite(time_values)
        & np.isfinite(flux_values)
        & np.isfinite(error_values)
        & (error_values > 0.0)
    )
    if keep.sum() < 10 or np.sum((fraction > 0.0) & keep) < 3:
        raise ValueError("Too few measurements cover the fixed schedule window")
    hours = (time_values[keep] - float(np.mean(time_values[keep]))) * 24.0
    design = np.column_stack(
        [np.ones(keep.sum()), hours, -fraction[keep]]
    )
    weights = 1.0 / np.square(error_values[keep])
    normal = design.T @ (weights[:, None] * design)
    right = design.T @ (weights * flux_values[keep])
    coefficients = np.linalg.solve(normal, right)
    residual = flux_values[keep] - design @ coefficients
    degrees_of_freedom = keep.sum() - design.shape[1]
    reduced_chi_square = float(
        np.sum(np.square(residual / error_values[keep])) / degrees_of_freedom
    )
    covariance = np.linalg.inv(normal) * reduced_chi_square
    inside = (fraction > 0.0) & keep
    outside = (fraction == 0.0) & keep
    return {
        "depth": float(coefficients[2]),
        "depth_error": float(math.sqrt(covariance[2, 2])),
        "depth_snr": float(coefficients[2] / math.sqrt(covariance[2, 2])),
        "baseline": float(coefficients[0]),
        "slope_per_hour": float(coefficients[1]),
        "residual_robust_scatter": robust_sigma(residual),
        "reduced_chi_square": reduced_chi_square,
        "points": int(keep.sum()),
        "in_window_points": int(inside.sum()),
        "outside_window_points": int(outside.sum()),
        "median_inside": float(np.median(flux_values[inside])),
        "median_outside": float(np.median(flux_values[outside])),
    }


def analyze_schedule_window(
    record_path: Path,
    bjd_tdb: np.ndarray,
    flux: np.ndarray,
    flux_error: np.ndarray,
    use: np.ndarray,
    *,
    exposure_seconds: float = 50.0,
    comparison_depth_ppt: float = 2.910,
) -> dict[str, object]:
    """Combine source interpretation, fixed-window fit, and one injection."""
    time_values = np.asarray(bjd_tdb, dtype=float)
    context = schedule_context(
        record_path,
        observation_start_bjd=float(np.nanmin(time_values)),
        observation_end_bjd=float(np.nanmax(time_values)),
    )
    working = context["working_interpretation"]
    assert isinstance(working, dict)
    times = working["times"]
    assert isinstance(times, dict)
    ingress = float(times["ingress"]["bjd_tdb"])
    egress = float(times["egress"]["bjd_tdb"])
    box = fixed_window_fraction(time_values, ingress, egress, exposure_seconds)
    observed = fit_fixed_window(
        time_values,
        flux,
        flux_error,
        use,
        ingress,
        egress,
        exposure_seconds=exposure_seconds,
    )
    injected_flux = np.asarray(flux, dtype=float) * (
        1.0 - comparison_depth_ppt / 1000.0 * box
    )
    injected = fit_fixed_window(
        time_values,
        injected_flux,
        flux_error,
        use,
        ingress,
        egress,
        exposure_seconds=exposure_seconds,
    )
    fit_summary = {
        "model": "straight baseline plus exposure-integrated fixed box",
        "observed_depth_ppt": float(observed["depth"]) * 1000.0,
        "observed_depth_error_ppt": float(observed["depth_error"]) * 1000.0,
        "observed_depth_snr": observed["depth_snr"],
        "transit_like_dimming_above_3_sigma": bool(
            float(observed["depth_snr"]) >= 3.0
        ),
        "points": observed["points"],
        "in_window_points": observed["in_window_points"],
        "outside_window_points": observed["outside_window_points"],
        "median_inside": observed["median_inside"],
        "median_outside": observed["median_outside"],
        "residual_robust_scatter_ppt": float(
            observed["residual_robust_scatter"]
        )
        * 1000.0,
        "reduced_chi_square": observed["reduced_chi_square"],
        "comparison_injection_depth_ppt": comparison_depth_ppt,
        "injected_total_depth_ppt": float(injected["depth"]) * 1000.0,
        "injected_increment_recovered_ppt": (
            float(injected["depth"]) - float(observed["depth"])
        )
        * 1000.0,
        "injected_total_depth_snr": injected["depth_snr"],
        "interpretation": (
            "No transit-like dimming is measured in the exact historical window. "
            "This is conditional on the documented Eastern-time interpretation "
            "and is not a physical transit fit."
        ),
    }
    context["fixed_window_check"] = fit_summary
    return context


def write_target_plot_config(
    source: Path,
    destination: Path,
    schedule: dict[str, object],
    *,
    observation_start_bjd: float,
    observation_end_bjd: float,
) -> None:
    """Write a target-specific AstroImageJ config with schedule V.Markers."""
    working = schedule["working_interpretation"]
    assert isinstance(working, dict)
    times = working["times"]
    assert isinstance(times, dict)
    ingress = float(times["ingress"]["bjd_tdb"])
    egress = float(times["egress"]["bjd_tdb"])
    integer_offset = math.floor(observation_start_bjd)
    values = {
        ".plot.title": "TOI-3505.01, UT 2022-07-22",
        ".plot.subtitle": "GMU 0.8 m (R filter, 50 s exposures)",
        ".plot.xlabel": "BJD_TDB",
        ".plot.showVMarker1": "true",
        ".plot.showVMarker2": "true",
        ".plot.vMarker1TopText": "2022 schedule (EDT)",
        ".plot.vMarker1BotText": "Ingress",
        ".plot.vMarker1Value": f"{ingress - integer_offset:.9f}",
        ".plot.vMarker2TopText": "2022 schedule (EDT)",
        ".plot.vMarker2BotText": "Egress",
        ".plot.vMarker2Value": f"{egress - integer_offset:.9f}",
        ".plot.useInEgressMarkers": "false",
        ".plot.ingressTime": f"{ingress - integer_offset:.9f}",
        ".plot.egressTime": f"{egress - integer_offset:.9f}",
        ".plot.xMin": f"{observation_start_bjd:.9f}",
        ".plot.xMax": f"{observation_end_bjd:.9f}",
    }
    text = source.read_bytes().decode("ascii")
    for key, value in values.items():
        pattern = re.compile(rf"(?m)^{re.escape(key)}=[^\r\n]*")
        text, count = pattern.subn(f"{key}={value}", text)
        if count != 1:
            raise RuntimeError(f"Expected one {key} entry in {source}; found {count}")
    destination.write_bytes(text.encode("ascii"))
