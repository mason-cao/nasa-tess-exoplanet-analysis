"""Check that public TOI-3505.01 products match the canonical JSON outputs."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
EPHEMERIS_PATH = (
    Path("outputs") / "toi3505_ephemeris_refined" / "ephemeris_refined.json"
)
ROBUSTNESS_PATH = (
    Path("outputs")
    / "toi3505_ephemeris_robustness"
    / "ephemeris_robustness.json"
)
RECONSTRUCTION_PATH = (
    Path("outputs")
    / "toi3505_schedule_reconstruction"
    / "schedule_reconstruction.json"
)
EXOFOP_CONTEXT_PATH = (
    Path("data") / "catalogs" / "toi3505" / "exofop_ground_followup.json"
)
RESEARCH_CONFIG_PATH = (
    Path("outputs") / "toi3505_research_record" / "frozen_analysis_config.json"
)
RESEARCH_CLAIMS_PATH = (
    Path("outputs") / "toi3505_research_record" / "claim_evidence.csv"
)
PAPER_VALUES_PATH = Path("outputs") / "toi3505_paper" / "manuscript_values.json"
PAPER_SOURCE_PATH = Path("paper") / "TOI-3505.01_manuscript.md"
PAPER_HTML_PATH = Path("paper") / "TOI-3505.01_manuscript.html"
PAPER_PDF_PATH = Path("output") / "pdf" / "TOI-3505.01_research_paper.pdf"
FILE_MANIFEST_PATH = (
    Path("outputs") / "toi3505_research_record" / "file_manifest.csv"
)
MANIFEST_SUMMARY_PATH = (
    Path("outputs") / "toi3505_research_record" / "manifest_summary.json"
)
LITERATURE_CONTEXT_PATH = (
    Path("data") / "catalogs" / "toi3505" / "literature_context.json"
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
EXPECTED_AUTHORS = (
    "Mason Cao",
    "Annalyse Dickinson",
    "Owen Alfaro",
    "Rianne Eccleston",
    "Kasey Davidson",
    "Kevin I. Collins",
    "Peter Plavchan",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--verify-manifest",
        action="store_true",
        help="Rehash every file in the research manifest and verify its summary.",
    )
    return parser.parse_args()


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


def sha256_file(path: Path, block_size: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(block_size):
            digest.update(block)
    return digest.hexdigest()


def verify_manifest(root: Path = ROOT) -> list[str]:
    """Verify every frozen path, size, digest, and manifest-summary total."""
    errors: list[str] = []
    manifest_path = root / FILE_MANIFEST_PATH
    summary_path = root / MANIFEST_SUMMARY_PATH
    if not manifest_path.is_file():
        return [f"missing file manifest: {manifest_path}"]
    if not summary_path.is_file():
        return [f"missing manifest summary: {summary_path}"]

    with manifest_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        required_fields = {"category", "path", "size_bytes", "sha256"}
        if reader.fieldnames is None or not required_fields.issubset(reader.fieldnames):
            return [f"{manifest_path}: required columns are missing"]
        rows = list(reader)

    seen: set[str] = set()
    total_bytes = 0
    category_counts: Counter[str] = Counter()
    root_resolved = root.resolve()
    for row_number, row in enumerate(rows, start=2):
        relative_path = row["path"]
        if relative_path in seen:
            errors.append(f"{manifest_path}:{row_number}: duplicate path {relative_path}")
            continue
        seen.add(relative_path)
        candidate = (root / relative_path).resolve()
        try:
            candidate.relative_to(root_resolved)
        except ValueError:
            errors.append(
                f"{manifest_path}:{row_number}: path escapes repository: {relative_path}"
            )
            continue
        if not candidate.is_file():
            errors.append(f"{manifest_path}:{row_number}: missing file {relative_path}")
            continue
        try:
            expected_size = int(row["size_bytes"])
        except ValueError:
            errors.append(
                f"{manifest_path}:{row_number}: invalid size for {relative_path}"
            )
            continue
        actual_size = candidate.stat().st_size
        total_bytes += actual_size
        category_counts[row["category"]] += 1
        if actual_size != expected_size:
            errors.append(
                f"{relative_path}: size {actual_size}, manifest says {expected_size}"
            )
        actual_digest = sha256_file(candidate)
        if actual_digest != row["sha256"]:
            errors.append(f"{relative_path}: SHA-256 differs from manifest")

    summary = load_json(summary_path)
    if int(summary.get("files", -1)) != len(rows):
        errors.append(f"{summary_path}: file count differs from manifest")
    if int(summary.get("bytes", -1)) != total_bytes:
        errors.append(f"{summary_path}: byte total differs from current files")
    expected_counts = {
        str(key): int(value)
        for key, value in dict(summary.get("category_counts", {})).items()
    }
    if expected_counts != dict(category_counts):
        errors.append(f"{summary_path}: category counts differ from manifest")
    original_count = sum(row["category"] == "original_archive" for row in rows)
    if int(summary.get("original_archive_count", -1)) != original_count:
        errors.append(f"{summary_path}: original-archive count differs from manifest")
    if not bool(summary.get("sha256_complete")):
        errors.append(f"{summary_path}: sha256_complete must be true")
    return errors


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
    paths.extend((root / "paper").glob("*.md"))
    paths.append(root / "output" / "pdf" / "README.md")
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

    presentation = ground_summary.get("figure_presentation")
    if not isinstance(presentation, dict) or not bool(
        presentation.get("reviewer_overlap_note_resolved")
    ):
        errors.append(
            f"{root / GROUND_SUMMARY_PATH}: ground-curve legend overlap is not resolved"
        )

    robustness_path = root / ROBUSTNESS_PATH
    robustness = load_json(robustness_path)
    robustness_primary = robustness.get("primary_tess_linear")
    if not isinstance(robustness_primary, dict):
        errors.append(f"{robustness_path}: primary TESS fit is missing")
    else:
        if int(robustness_primary.get("events", -1)) != events:
            errors.append(f"{robustness_path}: primary event count differs")
        numerical_tolerances = {
            "period_days": 1e-12,
            "period_error_days": 5e-14,
            "epoch_bjd_tdb": 2e-9,
        }
        for key, tolerance in numerical_tolerances.items():
            if not math.isclose(
                float(robustness_primary.get(key, math.nan)),
                float(adopted[key]),
                rel_tol=0.0,
                abs_tol=tolerance,
            ):
                errors.append(
                    f"{robustness_path}: primary {key} differs from canonical fit"
                )
    quadratic = robustness.get("quadratic_model_control")
    if not isinstance(quadratic, dict) or float(
        quadratic.get("delta_bic_quadratic_minus_linear", -math.inf)
    ) <= 0.0:
        errors.append(f"{robustness_path}: quadratic control must favor linear BIC")

    reconstruction_path = root / RECONSTRUCTION_PATH
    reconstruction = load_json(reconstruction_path)
    recovered = reconstruction.get("reconstruction")
    reconstruction_timeline = reconstruction.get("timeline")
    if not isinstance(recovered, dict):
        errors.append(f"{reconstruction_path}: reconstruction block is missing")
    else:
        if int(recovered.get("cycles", -1)) != 96:
            errors.append(f"{reconstruction_path}: expected 96 propagated cycles")
        if not bool(recovered.get("all_markers_agree_within_one_minute")):
            errors.append(
                f"{reconstruction_path}: schedule markers do not agree within one minute"
            )
    if not isinstance(reconstruction_timeline, dict):
        errors.append(f"{reconstruction_path}: reconstruction timeline is missing")
        current_schedule_offset = math.nan
    else:
        current_schedule_offset = float(
            reconstruction_timeline.get(
                "schedule_minus_adopted_midpoint_hours", math.nan
            )
        )
    revised_schedule_offset = (
        float(recovered.get("schedule_minus_displayed_period_midpoint_hours", math.nan))
        if isinstance(recovered, dict)
        else math.nan
    )

    exofop_path = root / EXOFOP_CONTEXT_PATH
    exofop = load_json(exofop_path)
    current_status = exofop.get("current_status")
    if not isinstance(current_status, dict) or current_status.get(
        "tfopwg_disposition"
    ) != "PC":
        errors.append(f"{exofop_path}: frozen TFOPWG disposition must be PC")

    text_expectations = {
        root / "README.md": (
            unicode_period,
            f"{events} transit times",
            f"{current_schedule_offset:.2f} hours",
            f"{revised_schedule_offset:.2f}-hour offset",
        ),
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
        / "toi3505_ephemeris_robustness"
        / "README.md": (
            f"{period:.10f} +/- {period_error:.10f} days",
            "four clusters",
            "external control",
        ),
        root
        / "outputs"
        / "toi3505_schedule_reconstruction"
        / "README.md": (
            "96` cycles",
            "56.1` seconds",
            "archival reconstruction, not proof",
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

    paper_values_path = root / PAPER_VALUES_PATH
    if not paper_values_path.is_file():
        errors.append(f"missing paper values record: {paper_values_path}")
    else:
        paper_values = load_json(paper_values_path)
        compare_adopted_ephemeris(
            str(paper_values_path), paper_values.get("primary_result"), adopted, errors
        )
        if tuple(paper_values.get("authors", ())) != EXPECTED_AUTHORS:
            errors.append(
                f"{paper_values_path}: author list must contain exactly the established six authors"
            )
        status_statement = str(paper_values.get("status_statement", "")).lower()
        if "planet candidate" not in status_statement:
            errors.append(f"{paper_values_path}: candidate status boundary is missing")
        if "does not claim validation or confirmation" not in status_statement:
            errors.append(f"{paper_values_path}: validation boundary is missing")

        literature = load_json(root / LITERATURE_CONTEXT_PATH)
        expected_novelty = literature.get("novelty_assessment", {}).get(
            "supported_wording"
        )
        if paper_values.get("novelty_wording") != expected_novelty:
            errors.append(
                f"{paper_values_path}: novelty wording differs from the literature snapshot"
            )
        source_files = paper_values.get("source_files")
        if not isinstance(source_files, list) or not all(
            isinstance(path, str) for path in source_files
        ):
            errors.append(f"{paper_values_path}: source-file inventory is invalid")
        else:
            required_sources = {
                "src/build_toi3505_paper.py",
                "paper/TOI-3505.01_manuscript.md",
                "outputs/toi3505_ephemeris_refined/event_times_best_per_sector.csv",
            }
            missing_sources = required_sources - set(source_files)
            if missing_sources:
                errors.append(
                    f"{paper_values_path}: missing build sources {sorted(missing_sources)}"
                )
            if len(source_files) != len(set(source_files)):
                errors.append(f"{paper_values_path}: duplicate source-file entries")
            for relative_path in source_files:
                if not (root / relative_path).is_file():
                    errors.append(
                        f"{paper_values_path}: inventoried source is missing: {relative_path}"
                    )

    paper_source_path = root / PAPER_SOURCE_PATH
    require_tokens(
        paper_source_path,
        (
            "{{PERIOD_SHORT}}",
            "{{NOVELTY_SENTENCE}}",
            "does not validate or confirm TOI-3505.01",
        ),
        errors,
    )
    paper_html_path = root / PAPER_HTML_PATH
    require_tokens(
        paper_html_path,
        (
            unicode_period,
            "To our knowledge,",
            "does not validate or confirm TOI-3505.01",
            *EXPECTED_AUTHORS,
        ),
        errors,
    )
    if paper_html_path.is_file():
        paper_html = paper_html_path.read_text(encoding="utf-8")
        if "{{" in paper_html or "}}" in paper_html:
            errors.append(f"{paper_html_path}: unreplaced template token remains")
        if paper_html.count("data:image/svg+xml;base64,") != 6:
            errors.append(
                f"{paper_html_path}: expected exactly six embedded lossless SVG figures"
            )

    paper_pdf_path = root / PAPER_PDF_PATH
    if not paper_pdf_path.is_file():
        errors.append(f"missing paper PDF: {paper_pdf_path}")
    elif paper_pdf_path.stat().st_size < 100_000:
        errors.append(f"{paper_pdf_path}: PDF is unexpectedly small")
    else:
        with paper_pdf_path.open("rb") as handle:
            if handle.read(5) != b"%PDF-":
                errors.append(f"{paper_pdf_path}: invalid PDF signature")

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
        (
            ascii_period,
            f"{events} accepted mid-transit times",
            "56.1 seconds",
            "Delta BIC is +1.71",
            "dispositions as PC",
            *nearby_tokens,
        ),
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
    args = parse_args()
    errors = check_repository()
    if args.verify_manifest:
        errors.extend(verify_manifest())
    if errors:
        print("TOI-3505.01 consistency check failed:")
        for error in errors:
            print(f"- {error}")
        raise SystemExit(1)
    suffix = " and every manifest hash" if args.verify_manifest else ""
    print(
        "TOI-3505.01 public products match the canonical ephemeris and package names"
        f"{suffix}."
    )


if __name__ == "__main__":
    main()
