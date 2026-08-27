"""Build the self-contained TOI-3505.01 research-paper HTML and PDF.

The editable manuscript is a Markdown template. Every quantitative token in
that template is populated from a canonical JSON or CSV analysis product, and
all four figures are embedded from lossless SVG files. This keeps the prose
editable while preventing headline values from drifting away from the frozen
analysis.
"""

from __future__ import annotations

import argparse
import base64
import html
import json
import re
import shutil
import subprocess
import tempfile
import time
from datetime import date
from pathlib import Path

import markdown
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
MANUSCRIPT_PATH = ROOT / "paper" / "TOI-3505.01_manuscript.md"
HTML_PATH = ROOT / "paper" / "TOI-3505.01_manuscript.html"
VALUES_PATH = ROOT / "outputs" / "toi3505_paper" / "manuscript_values.json"
PDF_PATH = ROOT / "output" / "pdf" / "TOI-3505.01_research_paper.pdf"

EPHEMERIS_PATH = (
    ROOT / "outputs" / "toi3505_ephemeris_refined" / "ephemeris_refined.json"
)
EVENT_TIMES_PATH = (
    ROOT
    / "outputs"
    / "toi3505_ephemeris_refined"
    / "event_times_best_per_sector.csv"
)
ROBUSTNESS_PATH = (
    ROOT
    / "outputs"
    / "toi3505_ephemeris_robustness"
    / "ephemeris_robustness.json"
)
RECONSTRUCTION_PATH = (
    ROOT
    / "outputs"
    / "toi3505_schedule_reconstruction"
    / "schedule_reconstruction.json"
)
GROUND_PATH = ROOT / "outputs" / "toi3505_final_candidate" / "summary.json"
GROUND_CHECKS_PATH = ROOT / "outputs" / "toi3505_ground_checks" / "summary.json"
GROUND_SEARCH_PATH = ROOT / "outputs" / "toi3505_ground_search" / "ground_search.json"
FALSE_POSITIVE_PATH = (
    ROOT / "outputs" / "toi3505_false_positive" / "false_positive_assessment.json"
)
DILUTION_PATH = ROOT / "outputs" / "toi3505_tess_pixels" / "dilution_screen.json"
VALIDATION_PATH = (
    ROOT / "outputs" / "toi3505_data_validation" / "analysis_summary.json"
)
TESS_PATH = ROOT / "outputs" / "toi3505_tess_analysis" / "analysis_summary.json"
PIXEL_PATH = ROOT / "outputs" / "toi3505_tess_pixels" / "analysis_summary.json"
CATALOG_PATH = ROOT / "data" / "catalogs" / "toi3505" / "catalog_summary.json"
TIC_PATH = ROOT / "data" / "catalogs" / "toi3505" / "tic_v8_2p5arcmin.csv"
EXOFOP_PATH = (
    ROOT / "data" / "catalogs" / "toi3505" / "exofop_ground_followup.json"
)
LITERATURE_PATH = (
    ROOT / "data" / "catalogs" / "toi3505" / "literature_context.json"
)

FIGURES = {
    "PHASE_FIGURE": (
        ROOT
        / "outputs"
        / "toi3505_tess_analysis"
        / "02_phase_folded_sectors.svg"
    ),
    "GROUND_FIGURE": (
        ROOT
        / "outputs"
        / "toi3505_final_candidate"
        / "TOI_3505.01_2022-07-22_R_light_curve.svg"
    ),
    "ROBUSTNESS_FIGURE": (
        ROOT
        / "outputs"
        / "toi3505_ephemeris_robustness"
        / "01_ephemeris_robustness.svg"
    ),
    "SCHEDULE_FIGURE": (
        ROOT
        / "outputs"
        / "toi3505_schedule_reconstruction"
        / "01_schedule_reconstruction.svg"
    ),
    "SEARCH_FIGURE": (
        ROOT
        / "outputs"
        / "toi3505_ground_search"
        / "01_ground_transit_search.svg"
    ),
    "FP_FIGURE": (
        ROOT
        / "outputs"
        / "toi3505_false_positive"
        / "01_false_positive_tests.svg"
    ),
}

TOKEN_PATTERN = re.compile(r"{{([A-Z0-9_]+)}}")
AUTHOR_NAMES = (
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
    parser.add_argument("--manuscript", type=Path, default=MANUSCRIPT_PATH)
    parser.add_argument("--html-output", type=Path, default=HTML_PATH)
    parser.add_argument("--values-output", type=Path, default=VALUES_PATH)
    parser.add_argument(
        "--pdf-output",
        type=Path,
        default=None,
        help="Also print the rendered HTML to this PDF (use 'default' for the standard path).",
    )
    parser.add_argument(
        "--chrome",
        type=Path,
        default=None,
        help="Optional Chrome/Chromium executable used with --pdf-output.",
    )
    return parser.parse_args()


def load_json(path: Path) -> dict[str, object]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError(f"Expected a JSON object in {path}")
    return payload


def svg_data_uri(path: Path) -> str:
    """Return a lossless, self-contained SVG data URI."""
    payload = base64.b64encode(path.read_bytes()).decode("ascii")
    return f"data:image/svg+xml;base64,{payload}"


def format_ra(ra_degrees: float) -> str:
    hours_total = (ra_degrees / 15.0) % 24.0
    hours = int(hours_total)
    minutes_total = (hours_total - hours) * 60.0
    minutes = int(minutes_total)
    seconds = (minutes_total - minutes) * 60.0
    return f"{hours:02d}:{minutes:02d}:{seconds:05.2f}"


def format_dec(dec_degrees: float) -> str:
    sign = "+" if dec_degrees >= 0 else "−"
    absolute = abs(dec_degrees)
    degrees = int(absolute)
    minutes_total = (absolute - degrees) * 60.0
    minutes = int(minutes_total)
    seconds = (minutes_total - minutes) * 60.0
    return f"{sign}{degrees:02d}:{minutes:02d}:{seconds:04.1f}"


def html_table(headers: list[str], rows: list[list[str]], classes: str = "") -> str:
    class_attribute = f' class="{html.escape(classes)}"' if classes else ""
    header_html = "".join(f"<th>{cell}</th>" for cell in headers)
    row_html = "".join(
        "<tr>" + "".join(f"<td>{cell}</td>" for cell in row) + "</tr>"
        for row in rows
    )
    return (
        f"<table{class_attribute}><thead><tr>{header_html}</tr></thead>"
        f"<tbody>{row_html}</tbody></table>"
    )


def collect_values() -> tuple[dict[str, str], dict[str, object]]:
    """Collect formatted template values and a machine-readable audit record."""
    ephemeris = load_json(EPHEMERIS_PATH)
    robustness = load_json(ROBUSTNESS_PATH)
    reconstruction = load_json(RECONSTRUCTION_PATH)
    ground = load_json(GROUND_PATH)
    ground_checks = load_json(GROUND_CHECKS_PATH)
    ground_search = load_json(GROUND_SEARCH_PATH)
    false_positive = load_json(FALSE_POSITIVE_PATH)
    fp_scenarios = false_positive["scenarios"]
    chromatic = fp_scenarios["blended_eclipsing_binary"]["chromatic_depth_test"]
    centroid_sigma = fp_scenarios["nearby_eclipsing_binary"][
        "spoc_difference_image_offset_sigma"
    ]
    velocity = fp_scenarios["eclipsing_binary_on_target"][
        "eclipsing_companion_velocity_bound"
    ]
    if not bool(chromatic["consistent_with_achromatic"]):
        raise RuntimeError(
            "The MuSCAT2 depths are no longer consistent with an achromatic "
            "event; Section 4.5 and Section 5.4 must be rewritten"
        )
    if not bool(velocity["stellar_companion_disfavoured"]):
        raise RuntimeError(
            "The velocity bound no longer disfavours a stellar companion; "
            "Section 4.5 and Section 5.4 must be rewritten"
        )
    dilution = load_json(DILUTION_PATH)
    validation = load_json(VALIDATION_PATH)
    spoc_records = validation["official_multisector_tce"]
    assert isinstance(spoc_records, list) and len(spoc_records) == 1
    spoc_tce = spoc_records[0]
    tess = load_json(TESS_PATH)
    pixels = load_json(PIXEL_PATH)
    catalog = load_json(CATALOG_PATH)
    exofop = load_json(EXOFOP_PATH)
    literature = load_json(LITERATURE_PATH)
    tic = pd.read_csv(TIC_PATH)
    target_rows = tic.loc[tic["ID"] == int(catalog["tic_id"])]
    if len(target_rows) != 1:
        raise RuntimeError("The TIC snapshot must contain exactly one target row")
    target = target_rows.iloc[0]

    adopted = ephemeris["ephemeris"]["adopted"]  # type: ignore[index]
    comparisons = ephemeris["comparisons"]
    primary = robustness["primary_tess_linear"]
    anchors = robustness["four_sector_anchor_fit"]
    sector_jackknife = robustness["delete_one_sector_jackknife"]
    event_jackknife = robustness["delete_one_event_jackknife"]
    quadratic = robustness["quadratic_model_control"]
    muscat2 = robustness["external_muscat2_control"]
    schedule = reconstruction["reconstruction"]
    timeline = reconstruction["timeline"]
    fixed = ground["historical_schedule"]["fixed_window_check"]
    coverage = ground["historical_schedule"]["observation_comparison"]
    nearby = ground_checks["nearby_star_screen"]
    search_durations = ground_search["durations"]
    assert isinstance(search_durations, list)
    search_by_label = {str(entry["label"]): entry for entry in search_durations}
    search_catalog = search_by_label["TOI catalog"]
    search_spoc = search_by_label["SPOC multi-sector"]
    if int(nearby["uncleared_with_eclipse_consistent_shape"]) != 0:
        raise RuntimeError(
            "A nearby star now shows an eclipse-consistent dimming; Section 3.5 "
            "states that none does and must be rewritten"
        )
    if not all(
        bool(entry["expected_depth_excluded_everywhere"]) for entry in search_durations
    ):
        raise RuntimeError(
            "The whole-sequence search no longer excludes the published depth "
            "at every searched midpoint"
        )
    current_status = exofop["current_status"]

    if current_status["tess_disposition"] != "PC":
        raise RuntimeError("The frozen TESS disposition is no longer PC")
    if current_status["tfopwg_disposition"] != "PC":
        raise RuntimeError("The frozen TFOPWG disposition is no longer PC")
    if int(adopted["events"]) != 27:
        raise RuntimeError("The manuscript is designed for the canonical 27-event fit")
    if abs(float(primary["period_days"]) - float(adopted["period_days"])) > 1e-11:
        raise RuntimeError("Robustness and canonical periods disagree")
    if not bool(ground["figure_presentation"]["reviewer_overlap_note_resolved"]):
        raise RuntimeError("The ground figure still has an unresolved presentation note")
    if float(quadratic["delta_bic_quadratic_minus_linear"]) <= 0:
        raise RuntimeError("The frozen quadratic control no longer favors the linear model")

    event_table = pd.read_csv(EVENT_TIMES_PATH)
    accepted_events = event_table.loc[event_table["used_in_ephemeris"].astype(bool)]
    event_counts = accepted_events.groupby(["sector", "pipeline"]).size().to_dict()

    sector_dates = {
        14: "2019 Jul 18–Aug 14",
        41: "2021 Jul 24–Aug 20",
        54: "2022 Jul 09–Aug 04",
        81: "2024 Jul 15–Aug 10",
    }
    sector_rows: list[list[str]] = []
    for measurement in tess["sector_measurements"]:
        sector = int(measurement["sector"])
        selected_pipeline = ephemeris["pipeline_choice"][str(sector)]["pipeline"]
        selected_cadence = float(
            ephemeris["pipeline_choice"][str(sector)]["cadence_minutes"]
        )
        count = int(event_counts[(sector, selected_pipeline)])
        sector_rows.append(
            [
                str(sector),
                sector_dates[sector],
                html.escape(str(selected_pipeline)),
                f"{selected_cadence:g}",
                str(count),
                f"{float(measurement['depth_ppt']):.2f} ± "
                f"{float(measurement['depth_error_ppt']):.2f}",
            ]
        )

    model_rows = [
        [
            "Primary TESS linear",
            str(int(adopted["events"])),
            f"{float(adopted['period_days']):.10f}",
            f"{float(adopted['period_error_days']):.10f}",
            f"{float(adopted['reduced_chi2']):.2f}",
            "Adopted",
        ],
        [
            "Four sector anchors",
            str(int(anchors["events"])),
            f"{float(anchors['period_days']):.10f}",
            f"{float(anchors['period_error_days']):.10f}",
            f"{float(anchors['reduced_chi2']):.2f}",
            "Conservative aggregation",
        ],
        [
            "QLP only",
            str(int(ephemeris["ephemeris"]["qlp_only"]["events"])),
            f"{float(ephemeris['ephemeris']['qlp_only']['period_days']):.10f}",
            f"{float(ephemeris['ephemeris']['qlp_only']['period_error_days']):.10f}",
            f"{float(ephemeris['ephemeris']['qlp_only']['reduced_chi2']):.2f}",
            "Pipeline control",
        ],
        [
            "TESS + MuSCAT2",
            str(int(muscat2["fit"]["events"])),
            f"{float(muscat2['fit']['period_days']):.10f}",
            f"{float(muscat2['fit']['period_error_days']):.10f}",
            f"{float(muscat2['fit']['reduced_chi2']):.2f}",
            "External control only",
        ],
    ]

    marker_rows = [
        [
            str(row["marker"]).capitalize(),
            f"{float(row['source_bjd_tdb']):.7f}",
            f"{float(row['propagated_bjd_tdb']):.7f}",
            f"{float(row['schedule_bjd_tdb']):.7f}",
            f"{float(row['schedule_minus_propagated_seconds']):+.1f}",
        ]
        for row in schedule["marker_comparison"]
    ]

    public_rows: list[list[str]] = []
    for row in exofop["time_series_inventory"]["rows"]:
        scope = str(row["result_scope"])
        if "not used" in scope.lower() or "no quantitative" in scope.lower():
            role = "Context only"
        elif "external ephemeris control" in scope.lower():
            role = "Timing control only"
        else:
            role = "Report-level context"
        public_rows.append(
            [
                html.escape(str(row["date_utc"])),
                html.escape(str(row["facility"])),
                html.escape(str(row["filter"])),
                html.escape(str(row["coverage"])),
                role,
            ]
        )

    inventory = exofop["imaging_and_spectroscopy_inventory"]
    recon = exofop["reconnaissance_spectroscopy_notes"]
    reverification = exofop["status_reverification"]
    if reverification["tfopwg_disposition"] != "PC":
        raise RuntimeError(
            "The re-verified TFOPWG disposition is no longer PC; the manuscript's "
            "status boundary must be rewritten before it can be built"
        )
    followup_rows: list[list[str]] = []
    for row in inventory["high_resolution_imaging"]:
        if "companion_delta_mag" in row:
            result = (
                f"{float(row['companion_separation_arcsec']):.3f} arcsec, "
                f"delta {row['companion_delta_band']} = "
                f"{float(row['companion_delta_mag']):.2f}"
            )
        else:
            result = (
                f"{float(row['companion_separation_arcsec']):.2f} arcsec, "
                f"delta J = {float(row['companion_delta_mag_j']):.2f}, "
                f"delta Ks = {float(row['companion_delta_mag_ks']):.2f}; "
                f"{row['additional_sources']}"
            )
        followup_rows.append(
            [
                html.escape(str(row["date_utc"])),
                html.escape(f"{row['facility']} / {row['instrument']}"),
                html.escape(str(row["technique"]).title()),
                html.escape(result),
                "Companion parameters" if row["used_numerically"] else "Inventory only",
            ]
        )
    for row in inventory["spectroscopy"]:
        facility = str(row["facility"])
        if "instrument" in row:
            facility = f"{facility} / {row['instrument']}"
        result = f"{int(row['epochs'])} epochs on file"
        if str(row["facility"]) == "TRES":
            result += (
                "; velocities in phase with the photometric ephemeris, reported as "
                "consistent with about 10 Jupiter masses, on hold"
            )
        followup_rows.append(
            [
                html.escape(str(row["year"])),
                html.escape(facility),
                "Spectroscopy",
                html.escape(result),
                "Inventory only",
            ]
        )

    values = {
        "PAPER_DATE": "22 August 2026",
        "TARGET_RA_HMS": format_ra(float(catalog["center_icrs_degrees"][0])),
        "TARGET_DEC_DMS": format_dec(float(catalog["center_icrs_degrees"][1])),
        "TARGET_RA_DEG": f"{float(catalog['center_icrs_degrees'][0]):.6f}",
        "TARGET_DEC_DEG": f"{float(catalog['center_icrs_degrees'][1]):.6f}",
        "TARGET_TMAG": f"{float(catalog['target_tmag']):.4f}",
        "TARGET_TEFF": f"{float(target['Teff']):.0f}",
        "TARGET_RADIUS": f"{float(target['rad']):.3f}",
        "CATALOG_PERIOD": f"{float(comparisons['catalog_period_days']):.7f}",
        "CATALOG_PERIOD_ERROR": f"{float(comparisons['catalog_period_error_days']):.7f}",
        "CATALOG_EPOCH": f"{float(tess['catalog_ephemeris']['epoch_bjd_tdb']):.6f}",
        "CATALOG_DEPTH": f"{float(tess['catalog_depth_ppt']):.2f}",
        "CATALOG_DURATION": f"{float(tess['catalog_duration_hours']):.3f}",
        "EVENTS": str(int(adopted["events"])),
        "BASELINE_YEARS": f"{float(ephemeris['baseline_years']):.2f}",
        "PERIOD": f"{float(adopted['period_days']):.10f}",
        "PERIOD_SHORT": f"{float(adopted['period_days']):.7f}",
        "PERIOD_ERROR": f"{float(adopted['period_error_days']):.10f}",
        "PERIOD_ERROR_SHORT": f"{float(adopted['period_error_days']):.7f}",
        "EPOCH": f"{float(adopted['epoch_bjd_tdb']):.9f}",
        "EPOCH_ERROR": f"{float(adopted['epoch_error_days']):.7f}",
        "REDUCED_CHI2": f"{float(adopted['reduced_chi2']):.3f}",
        "TIMING_RMS": f"{float(adopted['residual_rms_minutes']):.2f}",
        "ERROR_SCALE": f"{float(adopted['error_scale']):.3f}",
        "GAIN_CATALOG": f"{float(comparisons['precision_gain_over_catalog']):.2f}",
        "GAIN_SPOC": f"{float(comparisons['precision_gain_over_spoc']):.2f}",
        "ANCHOR_PERIOD": f"{float(anchors['period_days']):.10f}",
        "ANCHOR_ERROR": f"{float(anchors['period_error_days']):.10f}",
        "SECTOR_JACKKNIFE_ERROR": f"{float(sector_jackknife['jackknife_standard_error_days']):.10f}",
        "SECTOR_JACKKNIFE_RATIO": f"{float(sector_jackknife['ratio_to_primary_formal_error']):.2f}",
        "EVENT_JACKKNIFE_ERROR": f"{float(event_jackknife['jackknife_standard_error_days']):.10f}",
        "EVENT_JACKKNIFE_MAX_SHIFT": f"{float(event_jackknife['maximum_absolute_period_shift_days']):.10f}",
        "QUADRATIC_TERM": f"{float(quadratic['quadratic_term_days_per_cycle2']):+.3e}",
        "QUADRATIC_ERROR": f"{float(quadratic['quadratic_term_error_days_per_cycle2']):.3e}",
        "QUADRATIC_SIGMA": f"{float(quadratic['quadratic_term_significance_sigma']):+.2f}",
        "DELTA_CHI2": f"{float(quadratic['delta_chi_square_linear_minus_quadratic']):.2f}",
        "DELTA_BIC": f"{float(quadratic['delta_bic_quadratic_minus_linear']):+.2f}",
        "MUSCAT_PERIOD": f"{float(muscat2['fit']['period_days']):.10f}",
        "MUSCAT_ERROR": f"{float(muscat2['fit']['period_error_days']):.10f}",
        "MUSCAT_PRECISION_CHANGE": f"{float(muscat2['period_precision_change_percent']):.2f}",
        "GROUND_ARCHIVED": str(int(ground["measurements_archived"])),
        "GROUND_USED": str(int(ground["measurements_used_in_primary_curve"])),
        "GROUND_EXCLUDED": str(int(ground["measurements_excluded_from_primary_curve"])),
        "GROUND_DURATION": f"{float(ground['observation']['duration_hours']):.3f}",
        "GROUND_START": f"{float(ground['observation']['start_bjd_tdb']):.9f}",
        "GROUND_END": f"{float(ground['observation']['end_bjd_tdb']):.9f}",
        "GROUND_SCATTER": f"{float(ground['adopted_primary_robust_scatter_ppt']):.3f}",
        "GROUND_DEPTH": f"{float(fixed['observed_depth_ppt']):+.3f}",
        "GROUND_DEPTH_ERROR": f"{float(fixed['observed_depth_error_ppt']):.3f}",
        "GROUND_DEPTH_SNR": f"{float(fixed['observed_depth_snr']):+.2f}",
        "INJECTION_DEPTH": f"{float(fixed['comparison_injection_depth_ppt']):.2f}",
        "INJECTION_RECOVERED": f"{float(fixed['injected_increment_recovered_ppt']):.3f}",
        "INJECTION_TOTAL_SNR": f"{float(fixed['injected_total_depth_snr']):.2f}",
        "OLD_PERIOD": f"{float(reconstruction['periods']['superseded_period_days']):.7f}",
        "REVISED_PERIOD": f"{float(reconstruction['periods']['displayed_revised_period_days']):.7f}",
        "RECONSTRUCTION_CYCLES": str(int(schedule["cycles"])),
        "MAX_MARKER_OFFSET": f"{float(schedule['maximum_absolute_marker_offset_seconds']):.1f}",
        "REVISED_OFFSET_HOURS": f"{float(schedule['schedule_minus_displayed_period_midpoint_hours']):.2f}",
        "CURRENT_OFFSET_HOURS": f"{float(timeline['schedule_minus_adopted_midpoint_hours']):.2f}",
        "CURRENT_GROUND_OFFSET_HOURS": f"{abs(float(ground['closest_predicted_transit_midpoint']['hours_from_observation_start'])):.2f}",
        "NEARBY_BRIGHT": str(int(nearby["bright_enough_catalog_candidates"])),
        "NEARBY_CLEARED": str(int(nearby["sources_cleared_by_conditional_screen"])),
        "NEARBY_OVERLAP": str(int(nearby["source_apertures_overlapping_target"])),
        "NEARBY_FAINT": str(
            int(nearby["disposition_counts"]["not cleared - flux too low"])
        ),
        "NEARBY_WEAK": str(
            int(nearby["disposition_counts"]["not cleared - depth limit too weak"])
        ),
        "TFOP_BAND_CORRECTION": f"{float(nearby['tfop_band_correction_mag']):.1f}",
        "RV_SPAN": f"{float(recon['velocity_span_km_s']):.2f}",
        "RV_MASS": "10",
        "HOST_MASS": f"{float(velocity['assumed_host_mass_solar']):.2f}",
        "CHROMATIC_MEAN": f"{float(chromatic['mean_depth_ppt']):.2f}",
        "CHROMATIC_SCATTER": f"{float(chromatic['depth_scatter_ppt']):.2f}",
        "CHROMATIC_SLOPE": f"{float(chromatic['slope_ppt_per_100nm']):+.3f}",
        "CHROMATIC_SLOPE_ERROR": f"{float(chromatic['slope_error_ppt_per_100nm']):.3f}",
        "CHROMATIC_SIGMA": f"{float(chromatic['slope_sigma']):.2f}",
        "STELLAR_MIN_K": f"{float(velocity['smallest_stellar_scenario_km_s']):.0f}",
        "STELLAR_RATIO": f"{float(velocity['ratio_smallest_stellar_to_observed']):.0f}",
        "CENTROID_SIGMA": (
            "under 1" if max(float(v) for v in centroid_sigma) < 1.0
            else f"{max(float(v) for v in centroid_sigma):.1f}"
        ),
        "SEARCH_BEST_DEPTH": f"{float(search_catalog['best_depth_ppt']):.2f}",
        "SEARCH_BEST_ERROR": f"{float(search_catalog['best_depth_error_ppt']):.2f}",
        "SEARCH_BEST_SIGMA": f"{float(search_catalog['best_depth_snr']):.2f}",
        "SEARCH_BEST_P": f"{float(search_catalog['best_trials_corrected_probability']):.3f}",
        "SEARCH_LIMIT_CATALOG": f"{float(search_catalog['median_upper_limit_ppt']):.2f}",
        "SEARCH_LIMIT_SPOC": f"{float(search_spoc['median_upper_limit_ppt']):.2f}",
        "SEARCH_RANGE_LOW": f"{float(search_catalog['searched_midpoint_range_hours'][0]):.2f}",
        "SEARCH_RANGE_HIGH": f"{float(search_catalog['searched_midpoint_range_hours'][1]):.2f}",
        "SEARCH_DURATION_CATALOG": f"{float(search_catalog['duration_hours']):.2f}",
        "SEARCH_DURATION_SPOC": f"{float(search_spoc['duration_hours']):.2f}",
        "SEARCH_COVERAGE": f"{float(ground_search['minimum_event_coverage_fraction']):.0%}",
        "CATALOG_RADIUS": f"{float(tess['catalog_planet_radius_rearth']):.2f}",
        "CATALOG_INSOLATION": f"{float(tess['catalog_insolation_earth']):.0f}",
        "SPOC_DEPTH": f"{float(spoc_tce['fit_depth_ppt']):.2f}",
        "SPOC_IMPACT": f"{float(spoc_tce['fit_impact_parameter']):.2f}",
        "COMPANION_FLUX_RATIO": f"{float(dilution['unresolved_companion']['flux_ratio_using_delta_i_as_tess_band_proxy']):.3f}",
        "TARGET_FLUX_FRACTION": f"{float(dilution['screening_target_fraction_of_total_flux']):.3f}",
        "HOST_DEPTH_SCENARIO": f"{float(dilution['if_2p91_ppt_is_an_uncorrected_observed_depth']['target_host_depth_ppt']):.2f}",
        "PLANNED_SPAN": f"{float(coverage['planned_span_hours']):.2f}",
        "LATE_START": f"{float(coverage['late_start_hours']):.2f}",
        "EARLY_END": f"{float(coverage['early_end_hours']):.2f}",
        "PLANNED_COVERAGE": f"{float(coverage['planned_window_coverage_fraction']):.0%}",
        "PRE_INGRESS_BASELINE": f"{float(coverage['pre_ingress_baseline_hours']):.2f}",
        "POST_EGRESS_BASELINE": f"{float(coverage['post_egress_baseline_hours']):.2f}",
        "SCHEDULE_ROW_NOTE": str(coverage["schedule_row_note"]),
        "COMPANION_SEPARATION": f"{float(pixels['dilution_screen']['unresolved_companion']['separation_arcsec']):.3f}",
        "COMPANION_DELTA_I": f"{float(pixels['dilution_screen']['unresolved_companion']['delta_i_mag']):.1f}",
        "TIC_CONTAMINATION": f"{float(pixels['dilution_screen']['tic_catalog_contamination_ratio']):.3f}",
        "PUBLIC_NIGHTS": str(int(exofop["time_series_inventory"]["unique_observing_nights"])),
        "TESS_STATUS": html.escape(str(current_status["tess_disposition"])),
        "TFOP_STATUS": html.escape(str(current_status["tfopwg_disposition"])),
        "NOVELTY_SENTENCE": html.escape(
            str(literature["novelty_assessment"]["supported_wording"])
        ),
        "SECTOR_TABLE": html_table(
            [
                "Sector",
                "Public light-curve span (UTC)",
                "Timing pipeline",
                "Cadence (min)",
                "Accepted timings",
                "QLP depth (ppt)",
            ],
            sector_rows,
            "compact-table sector-table",
        ),
        "MODEL_TABLE": html_table(
            ["Fit", "N", "Period (d)", "1σ error (d)", "χ²ν", "Role"],
            model_rows,
            "compact-table model-table",
        ),
        "MARKER_TABLE": html_table(
            [
                "Marker",
                "2021 BJD<sub>TDB</sub>",
                "+96 cycles",
                "2022 schedule",
                "Schedule − propagation (s)",
            ],
            marker_rows,
            "compact-table marker-table",
        ),
        "PUBLIC_TABLE": html_table(
            ["UT date", "Facility", "Filter", "Coverage", "Use here"],
            public_rows,
            "compact-table public-table",
        ),
        "FOLLOWUP_TABLE": html_table(
            ["UT date", "Facility", "Kind", "Result on file", "Use here"],
            followup_rows,
            "compact-table public-table",
        ),
    }
    values.update({name: svg_data_uri(path) for name, path in FIGURES.items()})

    record = {
        "schema_version": 1,
        "paper_date": date(2026, 8, 22).isoformat(),
        "title": "A Four-Sector TESS Ephemeris for TOI-3505.01 and the Origin of a 2022 Ground-Based Null Observation",
        "authors": list(AUTHOR_NAMES),
        "status_statement": (
            "TOI-3505.01 remains a planet candidate; the paper does not claim "
            "validation or confirmation."
        ),
        "primary_result": {
            "events": int(adopted["events"]),
            "sectors": list(adopted["sectors"]),
            "baseline_years": float(ephemeris["baseline_years"]),
            "epoch_bjd_tdb": float(adopted["epoch_bjd_tdb"]),
            "epoch_error_days": float(adopted["epoch_error_days"]),
            "period_days": float(adopted["period_days"]),
            "period_error_days": float(adopted["period_error_days"]),
            "covariance_days2": float(adopted["covariance_days2"]),
            "reduced_chi2": float(adopted["reduced_chi2"]),
            "residual_rms_minutes": float(adopted["residual_rms_minutes"]),
        },
        "ground_result": {
            "used_measurements": int(ground["measurements_used_in_primary_curve"]),
            "robust_scatter_ppt": float(ground["adopted_primary_robust_scatter_ppt"]),
            "historical_window_depth_ppt": float(fixed["observed_depth_ppt"]),
            "historical_window_depth_error_ppt": float(
                fixed["observed_depth_error_ppt"]
            ),
        },
        "schedule_reconstruction": {
            "cycles": int(schedule["cycles"]),
            "superseded_period_days": float(
                reconstruction["periods"]["superseded_period_days"]
            ),
            "maximum_marker_offset_seconds": float(
                schedule["maximum_absolute_marker_offset_seconds"]
            ),
            "strength": reconstruction["assessment"]["strength"],
        },
        "novelty_wording": literature["novelty_assessment"]["supported_wording"],
        "source_files": [
            str(path.relative_to(ROOT))
            for path in (
                Path(__file__).resolve(),
                MANUSCRIPT_PATH,
                EPHEMERIS_PATH,
                EVENT_TIMES_PATH,
                ROBUSTNESS_PATH,
                RECONSTRUCTION_PATH,
                GROUND_PATH,
                GROUND_CHECKS_PATH,
                GROUND_SEARCH_PATH,
                FALSE_POSITIVE_PATH,
                DILUTION_PATH,
                VALIDATION_PATH,
                TESS_PATH,
                PIXEL_PATH,
                CATALOG_PATH,
                TIC_PATH,
                EXOFOP_PATH,
                LITERATURE_PATH,
                *FIGURES.values(),
            )
        ],
    }
    return values, record


STYLE = r"""
:root {
  --ink: #111827;
  --muted: #4b5563;
  --rule: #9ca3af;
  --light-rule: #d1d5db;
  --accent: #164e63;
  --soft: #f3f4f6;
}
@page {
  size: Letter;
  margin: 0.62in 0.68in 0.66in;
  @bottom-center {
    content: counter(page);
    color: #6b7280;
    font-family: Arial, Helvetica, sans-serif;
    font-size: 8pt;
  }
}
* { box-sizing: border-box; }
html { background: white; }
body {
  color: var(--ink);
  background: white;
  font-family: Georgia, "Times New Roman", serif;
  font-size: 10.25pt;
  line-height: 1.36;
  margin: 0 auto;
  max-width: 7.15in;
  text-rendering: optimizeLegibility;
}
a { color: var(--accent); text-decoration: none; }
p { margin: 0 0 0.56em; text-align: justify; hyphens: auto; }
h1, h2, h3, .paper-meta, .authors, .affiliations, .keywords, figcaption,
table, .abstract-label, .status-note {
  font-family: Arial, Helvetica, sans-serif;
}
h1 {
  font-size: 21pt;
  line-height: 1.12;
  letter-spacing: -0.02em;
  margin: 0.08in auto 0.15in;
  max-width: 6.75in;
  text-align: center;
}
h2 {
  border-bottom: 0.8pt solid var(--rule);
  break-after: avoid;
  font-size: 14pt;
  line-height: 1.2;
  margin: 0.22in 0 0.09in;
  padding-bottom: 0.035in;
}
h3 {
  break-after: avoid;
  font-size: 11.25pt;
  line-height: 1.2;
  margin: 0.14in 0 0.055in;
}
.paper-meta { color: var(--muted); font-size: 8.7pt; text-align: center; }
.authors { font-size: 10.5pt; font-weight: 600; line-height: 1.35; text-align: center; }
.affiliations { color: var(--muted); font-size: 8.7pt; line-height: 1.35; text-align: center; }
.abstract {
  border-bottom: 0.8pt solid var(--rule);
  border-top: 0.8pt solid var(--rule);
  margin: 0.17in 0 0.1in;
  padding: 0.11in 0.13in 0.08in;
}
.abstract p { margin-bottom: 0.35em; text-align: left; }
.abstract-label { font-size: 9.4pt; font-weight: 700; }
.keywords { color: var(--muted); font-size: 8.8pt; margin-bottom: 0.12in; }
.status-note {
  background: var(--soft);
  border-left: 3pt solid var(--accent);
  break-inside: avoid;
  font-size: 9.1pt;
  margin: 0.1in 0;
  padding: 0.075in 0.11in;
}
.equation {
  break-inside: avoid;
  font-size: 10.6pt;
  margin: 0.11in 0;
  text-align: center;
}
figure {
  break-inside: avoid;
  margin: 0.14in auto 0.16in;
  text-align: center;
}
figure img { display: block; height: auto; margin: 0 auto; max-width: 100%; }
figure.wide img { width: 7.05in; }
figure.medium img { width: 6.45in; }
figure.ground img { width: 5.75in; }
figcaption {
  color: #1f2937;
  font-size: 8.25pt;
  line-height: 1.28;
  margin: 0.07in auto 0;
  max-width: 6.9in;
  text-align: left;
}
table {
  border-collapse: collapse;
  break-inside: avoid;
  font-size: 7.8pt;
  line-height: 1.22;
  margin: 0.11in auto 0.15in;
  width: 100%;
}
thead { display: table-header-group; }
th {
  background: var(--soft);
  border-bottom: 1pt solid var(--rule);
  border-top: 1pt solid var(--rule);
  font-weight: 700;
  padding: 0.045in 0.04in;
  text-align: left;
  vertical-align: bottom;
}
td {
  border-bottom: 0.35pt solid var(--light-rule);
  padding: 0.04in;
  text-align: left;
  vertical-align: top;
}
.model-table { font-size: 7.5pt; }
.public-table { font-size: 7.35pt; }
.table-caption {
  break-after: avoid;
  font-family: Arial, Helvetica, sans-serif;
  font-size: 8.35pt;
  line-height: 1.26;
  margin: 0.1in 0 -0.06in;
  text-align: left;
}
.references { font-size: 8.35pt; line-height: 1.25; }
.references p { margin-bottom: 0.38em; padding-left: 0.2in; text-align: left; text-indent: -0.2in; }
.no-break { break-inside: avoid; }
.page-break { break-before: page; }
ul { margin: 0.04in 0 0.08in 0.22in; padding-left: 0.12in; }
li { margin-bottom: 0.25em; }
code { font-size: 0.92em; }
"""


def render_manuscript(template_path: Path) -> tuple[str, dict[str, object]]:
    values, record = collect_values()
    source = template_path.read_text(encoding="utf-8")
    requested = set(TOKEN_PATTERN.findall(source))
    unknown = requested - values.keys()
    if unknown:
        raise RuntimeError(f"Unknown manuscript tokens: {sorted(unknown)}")
    unused = values.keys() - requested
    if unused:
        raise RuntimeError(f"Collected values are not used by manuscript: {sorted(unused)}")
    rendered = TOKEN_PATTERN.sub(lambda match: values[match.group(1)], source)
    if TOKEN_PATTERN.search(rendered):
        raise RuntimeError("Unreplaced manuscript tokens remain")

    body = markdown.markdown(
        rendered,
        extensions=("extra", "smarty", "sane_lists"),
        output_format="html5",
    )
    document = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{html.escape(str(record['title']))}</title>
<style>{STYLE}</style>
</head>
<body>
{body}
</body>
</html>
"""
    return document, record


def find_chrome(explicit: Path | None) -> Path:
    if explicit is not None:
        if explicit.is_file():
            return explicit.resolve()
        raise FileNotFoundError(f"Chrome executable not found: {explicit}")
    candidates = (
        Path("/Applications/Google Chrome.app/Contents/MacOS/Google Chrome"),
        Path("/Applications/Chromium.app/Contents/MacOS/Chromium"),
    )
    for candidate in candidates:
        if candidate.is_file():
            return candidate
    for command in ("google-chrome", "chromium", "chromium-browser"):
        resolved = shutil.which(command)
        if resolved:
            return Path(resolved)
    raise FileNotFoundError("Chrome or Chromium is required to print the PDF")


def is_valid_pdf(path: Path) -> bool:
    """Return whether *path* looks like a non-trivial PDF artifact."""
    if not path.is_file() or path.stat().st_size < 1024:
        return False
    with path.open("rb") as handle:
        return handle.read(5) == b"%PDF-"


def stop_process(process: subprocess.Popen[bytes]) -> None:
    """Stop a headless browser that remains alive after writing its PDF."""
    if process.poll() is not None:
        return
    process.terminate()
    try:
        process.wait(timeout=5)
    except subprocess.TimeoutExpired:
        process.kill()
        process.wait(timeout=5)


def print_pdf(html_path: Path, pdf_path: Path, chrome: Path | None) -> None:
    browser = find_chrome(chrome)
    pdf_path.parent.mkdir(parents=True, exist_ok=True)
    pdf_path.unlink(missing_ok=True)
    with tempfile.TemporaryDirectory(
        prefix="toi3505-paper-chrome-", ignore_cleanup_errors=True
    ) as profile:
        command = [
            str(browser),
            "--headless=new",
            "--disable-gpu",
            "--no-pdf-header-footer",
            "--allow-file-access-from-files",
            f"--user-data-dir={profile}",
            f"--print-to-pdf={pdf_path.resolve()}",
            html_path.resolve().as_uri(),
        ]
        process = subprocess.Popen(
            command,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
        deadline = time.monotonic() + 45.0
        stable_since: float | None = None
        previous_size = -1
        while process.poll() is None and time.monotonic() < deadline:
            if is_valid_pdf(pdf_path):
                current_size = pdf_path.stat().st_size
                if current_size != previous_size:
                    previous_size = current_size
                    stable_since = time.monotonic()
                elif stable_since is not None and time.monotonic() - stable_since >= 1.0:
                    break
            time.sleep(0.1)
        return_code = process.poll()
        stop_process(process)
    if not is_valid_pdf(pdf_path):
        suffix = "" if return_code in (None, 0) else f" (exit status {return_code})"
        raise RuntimeError(f"Chrome PDF generation failed{suffix}")


def main() -> None:
    args = parse_args()
    manuscript = args.manuscript.resolve()
    html_output = args.html_output.resolve()
    values_output = args.values_output.resolve()
    document, record = render_manuscript(manuscript)
    html_output.parent.mkdir(parents=True, exist_ok=True)
    values_output.parent.mkdir(parents=True, exist_ok=True)
    html_output.write_text(document, encoding="utf-8")
    values_output.write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote {html_output}")
    print(f"Wrote {values_output}")

    if args.pdf_output is not None:
        pdf_output = (
            PDF_PATH
            if str(args.pdf_output).lower() == "default"
            else args.pdf_output.resolve()
        )
        print_pdf(html_output, pdf_output, args.chrome)
        print(f"Wrote {pdf_output}")


if __name__ == "__main__":
    main()
