"""Check that public TOI-3505.01 products match the canonical JSON outputs."""

from __future__ import annotations

import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
EPHEMERIS_PATH = (
    Path("outputs") / "toi3505_ephemeris_refined" / "ephemeris_refined.json"
)
RESEARCH_CONFIG_PATH = (
    Path("outputs") / "toi3505_research_record" / "frozen_analysis_config.json"
)
RESEARCH_CLAIMS_PATH = (
    Path("outputs") / "toi3505_research_record" / "claim_evidence.csv"
)
OBSERVATION_PLAN_PATH = (
    Path("outputs") / "toi3505_observation_plan" / "observation_plan.json"
)
GROUND_SUMMARY_PATH = Path("outputs") / "toi3505_final_candidate" / "summary.json"
GROUND_CHECKS_PATH = Path("outputs") / "toi3505_ground_checks" / "summary.json"
REVIEW_PACKAGE_PATH = Path("outputs") / "toi3505_review_package"
PRIVATE_PLATFORM_NAME = "dis" + "cord"
LEGACY_PACKAGE_PATH = Path("outputs") / f"toi3505_{PRIVATE_PLATFORM_NAME}_post"
REVIEW_FILES = (
    "01_TOI_3505.01_final_light_curve.png",
    "02_TOI_3505.01_data_set_fit_settings.png",
    "03_TOI_3505.01_NEB_screen.png",
)


def load_json(path: Path) -> dict[str, object]:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except FileNotFoundError:
        raise
    if not isinstance(payload, dict):
        raise ValueError(f"Expected a JSON object in {path}")
    return payload


def require_tokens(path: Path, tokens: tuple[str, ...], errors: list[str]) -> None:
    if not path.is_file():
        errors.append(f"missing file: {path}")
        return
    text = path.read_text(encoding="utf-8")
    for token in tokens:
        if token not in text:
            errors.append(f"{path}: missing canonical text {token!r}")


def compare_adopted_ephemeris(
    label: str,
    actual: object,
    expected: dict[str, object],
    errors: list[str],
) -> None:
    if not isinstance(actual, dict):
        errors.append(f"{label}: adopted ephemeris is missing or not an object")
        return
    keys = (
        "events",
        "epoch_bjd_tdb",
        "epoch_error_days",
        "period_days",
        "period_error_days",
        "covariance_days2",
    )
    for key in keys:
        if key not in actual:
            errors.append(f"{label}: missing adopted ephemeris field {key}")
            continue
        expected_value = expected[key]
        actual_value = actual[key]
        if isinstance(expected_value, (int, float)) and isinstance(
            actual_value, (int, float)
        ):
            if not math.isclose(
                float(actual_value), float(expected_value), rel_tol=1e-12, abs_tol=1e-15
            ):
                errors.append(
                    f"{label}: {key}={actual_value!r}, expected {expected_value!r}"
                )
        elif actual_value != expected_value:
            errors.append(
                f"{label}: {key}={actual_value!r}, expected {expected_value!r}"
            )


def public_text_files(root: Path) -> list[Path]:
    paths = [root / "README.md", root / ".gitignore"]
    paths.extend((root / "src").glob("*.py"))
    paths.extend((root / "src").glob("*.html"))
    paths.extend((root / "outputs").glob("*/README.md"))
    return sorted(path for path in paths if path.is_file())


def check_repository(root: Path = ROOT) -> list[str]:
    errors: list[str] = []
    ephemeris_file = root / EPHEMERIS_PATH
    if not ephemeris_file.is_file():
        return [f"missing canonical ephemeris: {ephemeris_file}"]
    payload = load_json(ephemeris_file)
    try:
        adopted = payload["ephemeris"]["adopted"]  # type: ignore[index]
        comparisons = payload["comparisons"]
    except (KeyError, TypeError) as error:
        return [f"invalid canonical ephemeris structure: {error}"]
    if not isinstance(adopted, dict) or not isinstance(comparisons, dict):
        return ["invalid adopted ephemeris or comparison object"]

    events = int(adopted["events"])
    period = float(adopted["period_days"])
    period_error = float(adopted["period_error_days"])
    gain_catalog = float(comparisons["precision_gain_over_catalog"])
    gain_spoc = float(comparisons["precision_gain_over_spoc"])
    unicode_period = f"{period:.7f} ± {period_error:.7f}"
    ascii_period = f"{period:.7f} +/- {period_error:.7f}"
    html_period = f"{period:.7f} &plusmn; {period_error:.7f}"
    ground_summary = load_json(root / GROUND_SUMMARY_PATH)
    try:
        observation_start = float(
            ground_summary["observation"]["start_bjd_tdb"]  # type: ignore[index]
        )
    except (KeyError, TypeError) as error:
        return [f"invalid ground observation summary: {error}"]
    nearest_cycle = round(
        (observation_start - float(adopted["epoch_bjd_tdb"])) / period
    )
    nearest_midpoint = float(adopted["epoch_bjd_tdb"]) + nearest_cycle * period
    ground_offset_hours = abs(observation_start - nearest_midpoint) * 24.0
    poster_offset_text = f"{ground_offset_hours:.1f} hours before"

    text_expectations = {
        root / "README.md": (unicode_period, f"{events} transit times"),
        root
        / "outputs"
        / "toi3505_ephemeris_refined"
        / "README.md": (
            ascii_period,
            f"| Events | {events} across",
            f"{gain_catalog:.2f} times tighter",
            f"{gain_spoc:.2f} times tighter",
        ),
        root
        / "outputs"
        / "toi3505_poster"
        / "README.md": (
            ascii_period,
            f"{events} mid-transit times",
            f"{gain_catalog:.2f} times tighter",
            f"{gain_spoc:.2f} times tighter",
        ),
        root
        / "src"
        / "poster_template.html": (
            html_period,
            f"{events} transits",
            poster_offset_text,
        ),
        root
        / "src"
        / "poster_template_v2.html": (
            html_period,
            f"{events} transits",
            poster_offset_text,
        ),
        root
        / "outputs"
        / "toi3505_poster"
        / "TOI-3505.01_Mason_Cao_poster.html": (
            html_period,
            f"{events} transits",
            poster_offset_text,
        ),
        root
        / "outputs"
        / "toi3505_poster"
        / "TOI-3505.01_Mason_Cao_poster_v2.html": (
            html_period,
            f"{events} transits",
            poster_offset_text,
        ),
    }
    private_record = root / "docs" / "project-record.md"
    if private_record.exists():
        text_expectations[private_record] = (
            unicode_period,
            f"{events} mid-transit times",
            f"{gain_catalog:.2f} times more precise",
            f"{gain_spoc:.2f} times more precise",
            f"{ground_offset_hours:.2f} hours before",
        )
    for path, tokens in text_expectations.items():
        require_tokens(path, tokens, errors)

    research_path = root / RESEARCH_CONFIG_PATH
    if not research_path.is_file():
        errors.append(f"missing research configuration: {research_path}")
    else:
        research = load_json(research_path)
        tess = research.get("tess")
        research_adopted = (
            tess.get("adopted_ephemeris") if isinstance(tess, dict) else None
        )
        compare_adopted_ephemeris(str(research_path), research_adopted, adopted, errors)

    ground_checks = load_json(root / GROUND_CHECKS_PATH)
    nearby = ground_checks.get("nearby_star_screen")
    nearby_tokens: tuple[str, ...] = ()
    if not isinstance(nearby, dict):
        errors.append(
            f"{root / GROUND_CHECKS_PATH}: nearby-star screen is missing or invalid"
        )
    else:
        try:
            nearby_tokens = (
                f"{int(nearby['sources_cleared_by_conditional_screen'])} of "
                f"{int(nearby['bright_enough_catalog_candidates'])}",
            )
        except (KeyError, TypeError, ValueError) as error:
            errors.append(
                f"{root / GROUND_CHECKS_PATH}: invalid nearby-star counts: {error}"
            )

    claims_path = root / RESEARCH_CLAIMS_PATH
    require_tokens(
        claims_path,
        (ascii_period, f"{events} accepted mid-transit times", *nearby_tokens),
        errors,
    )

    plan_path = root / OBSERVATION_PLAN_PATH
    if not plan_path.is_file():
        errors.append(f"missing observation plan: {plan_path}")
    else:
        plan = load_json(plan_path)
        compare_adopted_ephemeris(
            str(plan_path), plan.get("adopted_ephemeris"), adopted, errors
        )

    review_package = root / REVIEW_PACKAGE_PATH
    for name in REVIEW_FILES:
        if not (review_package / name).is_file():
            errors.append(f"missing mentor-review asset: {review_package / name}")
    legacy_package = root / LEGACY_PACKAGE_PATH
    if legacy_package.exists():
        errors.append(
            f"legacy communication-specific package still exists: {legacy_package}"
        )

    for path in public_text_files(root):
        text = path.read_text(encoding="utf-8").lower()
        if PRIVATE_PLATFORM_NAME in text:
            errors.append(
                f"public text still names a private communication platform: {path}"
            )
    return errors


def main() -> None:
    errors = check_repository()
    if errors:
        print("TOI-3505.01 consistency check failed:")
        for error in errors:
            print(f"- {error}")
        raise SystemExit(1)
    print(
        "TOI-3505.01 public products match the canonical ephemeris and package names."
    )


if __name__ == "__main__":
    main()
