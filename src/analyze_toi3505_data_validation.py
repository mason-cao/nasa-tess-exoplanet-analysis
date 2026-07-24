"""Compare official SPOC reports with the project TESS measurements.

The SPOC Data Validation products are available for TOI-3505.01 in Sectors 54
and 81. A later SPOC search also combined those two sectors at 10-minute
cadence. This script reads the machine-readable report XML, checks the matching
DVT FITS headers and Exo.MAST table, verifies which sectors actually
contributed, and places a small set of official diagnostics beside the
project's QLP, SPOC-light-curve, and difference-image measurements.

This is a comparison of methods applied to the same TESS observations. It is
not an independent observation, a statistical validation, or a replacement
for resolving the known close companion.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import re
import tempfile
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Any

os.environ.setdefault(
    "MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "nasa-tess-matplotlib-cache")
)

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.io import fits


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DATA_DIR = ROOT / "data" / "tess" / "toi3505" / "data_validation"
DEFAULT_MULTISECTOR_DATA_DIR = (
    ROOT / "data" / "tess" / "toi3505" / "data_validation_multi_sector"
)
DEFAULT_EXOMAST_DIR = ROOT / "data" / "tess" / "toi3505" / "exomast"
DEFAULT_RELEASE_NOTES_DIR = (
    ROOT / "data" / "tess" / "toi3505" / "release_notes"
)
DEFAULT_TESS_ANALYSIS_DIR = ROOT / "outputs" / "toi3505_tess_analysis"
DEFAULT_PIXEL_DIR = ROOT / "outputs" / "toi3505_tess_pixels"
DEFAULT_OUTPUT_DIR = ROOT / "outputs" / "toi3505_data_validation"
DV_NAMESPACE = {"dv": "http://www.nasa.gov/2018/TESS/DV"}
EXPECTED_SECTORS = (54, 81)
MULTISECTOR_SCOPE = "s0014-s0086"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=DEFAULT_DATA_DIR)
    parser.add_argument(
        "--multisector-data-dir",
        type=Path,
        default=DEFAULT_MULTISECTOR_DATA_DIR,
    )
    parser.add_argument("--exomast-dir", type=Path, default=DEFAULT_EXOMAST_DIR)
    parser.add_argument(
        "--release-notes-dir",
        type=Path,
        default=DEFAULT_RELEASE_NOTES_DIR,
    )
    parser.add_argument(
        "--tess-analysis-dir", type=Path, default=DEFAULT_TESS_ANALYSIS_DIR
    )
    parser.add_argument("--pixel-dir", type=Path, default=DEFAULT_PIXEL_DIR)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    return parser.parse_args()


def required(parent: ET.Element, path: str) -> ET.Element:
    element = parent.find(path, DV_NAMESPACE)
    if element is None:
        raise ValueError(f"Missing required Data Validation XML element: {path}")
    return element


def float_attribute(element: ET.Element, name: str) -> float:
    if name not in element.attrib:
        tag = element.tag.split("}")[-1]
        raise ValueError(f"Missing {name!r} on Data Validation XML element {tag!r}")
    return float(element.attrib[name])


def bool_attribute(element: ET.Element, name: str) -> bool:
    value = element.attrib.get(name, "").strip().lower()
    if value not in {"true", "false"}:
        tag = element.tag.split("}")[-1]
        raise ValueError(f"Invalid {name!r} on Data Validation XML element {tag!r}")
    return value == "true"


def parameter_map(fit: ET.Element) -> dict[str, tuple[float, float]]:
    parameters: dict[str, tuple[float, float]] = {}
    container = required(fit, "dv:modelParameters")
    for parameter in container.findall("dv:modelParameter", DV_NAMESPACE):
        name = parameter.attrib.get("name")
        if not name:
            continue
        parameters[name] = (
            float_attribute(parameter, "value"),
            float_attribute(parameter, "uncertainty"),
        )
    return parameters


def required_parameter(
    parameters: dict[str, tuple[float, float]], name: str
) -> tuple[float, float]:
    if name not in parameters:
        raise ValueError(f"Missing required Data Validation fit parameter: {name}")
    return parameters[name]


def parse_dv_xml(path: Path) -> dict[str, Any]:
    """Extract the report values used in the project comparison."""
    root = ET.parse(path).getroot()
    planet = required(root, "dv:planetResults")
    all_fit = required(planet, "dv:allTransitsFit")
    parameters = parameter_map(all_fit)

    period, period_error = required_parameter(parameters, "orbitalPeriodDays")
    epoch, epoch_error = required_parameter(parameters, "transitEpochBtjd")
    duration, duration_error = required_parameter(parameters, "transitDurationHours")
    depth_ppm, depth_error_ppm = required_parameter(parameters, "transitDepthPpm")
    radius_ratio, radius_ratio_error = required_parameter(
        parameters, "ratioPlanetRadiusToStarRadius"
    )
    impact, impact_error = required_parameter(parameters, "minImpactParameter")

    binary = required(
        planet,
        "dv:binaryDiscriminationResults/"
        "dv:oddEvenTransitDepthComparisonStatistic",
    )
    odd_even_statistic = float_attribute(binary, "value")
    odd_even_significance = float_attribute(binary, "significance")

    bootstrap = required(planet, "dv:bootstrapResults")
    candidate = required(planet, "dv:planetCandidate")
    weak_secondary = required(candidate, "dv:weakSecondary")
    secondary_depth = required(weak_secondary, "dv:depthPpm")

    motion = required(
        planet, "dv:centroidResults/dv:differenceImageMotionResults"
    )
    tic_sky = required(
        motion, "dv:msTicCentroidOffsets/dv:meanSkyOffset"
    )
    control_sky = required(
        motion, "dv:msControlCentroidOffsets/dv:meanSkyOffset"
    )
    quality_summary = required(motion, "dv:summaryQualityMetric")
    difference = required(planet, "dv:differenceImageResults")
    quality = required(difference, "dv:qualityMetric")

    ghost = required(planet, "dv:ghostDiagnosticResults")
    ghost_core = required(ghost, "dv:coreApertureCorrelationStatistic")
    ghost_halo = required(ghost, "dv:haloApertureCorrelationStatistic")
    core_value = float_attribute(ghost_core, "value")
    halo_value = float_attribute(ghost_halo, "value")

    tic_offset = float_attribute(tic_sky, "value")
    tic_offset_error = float_attribute(tic_sky, "uncertainty")
    control_offset = float_attribute(control_sky, "value")
    control_offset_error = float_attribute(control_sky, "uncertainty")

    return {
        "sector": int(difference.attrib["sector"]),
        "tic_id": int(root.attrib["ticId"]),
        "toi_id": root.attrib.get("toiId", ""),
        "spoc_transit_model": all_fit.attrib.get("transitModelName", ""),
        "fit_full_convergence": bool_attribute(all_fit, "fullConvergence"),
        "fit_period_days": period,
        "fit_period_error_days": period_error,
        "fit_epoch_btjd": epoch,
        "fit_epoch_error_days": epoch_error,
        "fit_duration_hours": duration,
        "fit_duration_error_hours": duration_error,
        "fit_depth_ppt": depth_ppm / 1000.0,
        "fit_depth_error_ppt": depth_error_ppm / 1000.0,
        "fit_radius_ratio": radius_ratio,
        "fit_radius_ratio_error": radius_ratio_error,
        "fit_impact_parameter": impact,
        "fit_impact_parameter_error": impact_error,
        "fit_snr": float_attribute(all_fit, "modelFitSnr"),
        "maximum_multiple_event_statistic": float_attribute(
            candidate, "maxMultipleEventSigma"
        ),
        "observed_transits": int(candidate.attrib["observedTransitCount"]),
        "expected_transits": int(candidate.attrib["expectedTransitCount"]),
        "suspected_eclipsing_binary": bool_attribute(
            candidate, "suspectedEclipsingBinary"
        ),
        # The report's "value in sigmas" is sqrt(test statistic).
        "odd_even_difference_sigma": math.sqrt(max(odd_even_statistic, 0.0)),
        "odd_even_planet_favoring_percent": odd_even_significance * 100.0,
        "weak_secondary_mes": float_attribute(weak_secondary, "maxMes"),
        "weak_secondary_phase_days": float_attribute(
            weak_secondary, "maxMesPhaseInDays"
        ),
        "weak_secondary_depth_ppt": float_attribute(
            secondary_depth, "value"
        )
        / 1000.0,
        "weak_secondary_depth_error_ppt": float_attribute(
            secondary_depth, "uncertainty"
        )
        / 1000.0,
        "bootstrap_false_alarm_probability": float_attribute(
            bootstrap, "significance"
        ),
        "centroid_to_tic_arcsec": tic_offset,
        "centroid_to_tic_error_arcsec": tic_offset_error,
        "centroid_to_tic_sigma": (
            abs(tic_offset) / tic_offset_error if tic_offset_error > 0 else np.nan
        ),
        "centroid_to_out_of_transit_arcsec": control_offset,
        "centroid_to_out_of_transit_error_arcsec": control_offset_error,
        "centroid_to_out_of_transit_sigma": (
            abs(control_offset) / control_offset_error
            if control_offset_error > 0
            else np.nan
        ),
        "difference_image_quality": float_attribute(quality, "value"),
        "difference_images_good_fraction": float_attribute(
            quality_summary, "fractionOfGoodMetrics"
        ),
        "difference_image_transits": int(difference.attrib["numberOfTransits"]),
        "ghost_core_statistic": core_value,
        "ghost_halo_statistic": halo_value,
        "ghost_core_halo_ratio": (
            core_value / halo_value if halo_value != 0 else np.nan
        ),
    }


def parse_dvt_header(path: Path) -> dict[str, Any]:
    """Read the compact fit summary stored in the official DVT FITS file."""
    with fits.open(path, memmap=True) as hdus:
        primary = hdus[0].header
        tce_hdu = hdus["TCE_1"]
        tce = tce_hdu.header
        sector_value = primary.get("SECTOR")
        sector_vector = str(primary.get("SECTORS", ""))
        table_names = set(tce_hdu.columns.names or [])
        if "LC_DETREND" in table_names:
            finite_detrended_rows = int(
                np.count_nonzero(np.isfinite(tce_hdu.data["LC_DETREND"]))
            )
        else:
            finite_detrended_rows = None
        return {
            "sector": int(sector_value) if sector_value is not None else None,
            "tic_id": int(primary["TICID"]),
            "dvt_data_release": int(primary.get("DATA_REL", -1)),
            "dvt_sector_vector": sector_vector,
            "dvt_sectors_used": ";".join(
                str(value) for value in decode_spoc_sector_vector(sector_vector)
            ),
            "dvt_time_rows": int(len(tce_hdu.data)),
            "dvt_finite_detrended_rows": finite_detrended_rows,
            "dvt_cadence_minutes": float(tce.get("TIMEDEL", np.nan)) * 1440.0,
            "dvt_period_days": float(tce["TPERIOD"]),
            "dvt_epoch_btjd": float(tce["TEPOCH"]),
            "dvt_duration_hours": float(tce["TDUR"]),
            "dvt_depth_ppt": float(tce["TDEPTH"]) / 1000.0,
            "dvt_maximum_multiple_event_statistic": float(tce["MAXMES"]),
        }


def decode_spoc_sector_vector(vector: str) -> tuple[int, ...]:
    """Decode the SPOC bit vector whose character index is the sector number."""
    if vector and set(vector) - {"0", "1"}:
        raise ValueError("SPOC sector vector contains values other than 0 and 1")
    return tuple(index for index, value in enumerate(vector) if value == "1")


def fit_values_match(xml_row: dict[str, Any], dvt_row: dict[str, Any]) -> bool:
    """Compare fit values shared by an XML report and its DVT header."""
    return bool(
        xml_row["tic_id"] == dvt_row["tic_id"]
        and math.isclose(
            xml_row["fit_period_days"],
            dvt_row["dvt_period_days"],
            rel_tol=0.0,
            abs_tol=5e-5,
        )
        and math.isclose(
            xml_row["fit_epoch_btjd"],
            dvt_row["dvt_epoch_btjd"],
            rel_tol=0.0,
            abs_tol=5e-4,
        )
        and math.isclose(
            xml_row["fit_duration_hours"],
            dvt_row["dvt_duration_hours"],
            rel_tol=0.0,
            abs_tol=0.02,
        )
        and math.isclose(
            xml_row["fit_depth_ppt"],
            dvt_row["dvt_depth_ppt"],
            rel_tol=0.0,
            abs_tol=0.01,
        )
        and math.isclose(
            xml_row["maximum_multiple_event_statistic"],
            dvt_row["dvt_maximum_multiple_event_statistic"],
            rel_tol=0.0,
            abs_tol=0.01,
        )
    )


def dvt_matches_xml(xml_row: dict[str, Any], dvt_row: dict[str, Any]) -> bool:
    """Check that the XML and FITS belong to the same report solution."""
    return bool(
        xml_row["sector"] == dvt_row["sector"]
        and fit_values_match(xml_row, dvt_row)
    )


def relative(path: Path) -> str:
    try:
        return path.resolve().relative_to(ROOT).as_posix()
    except ValueError:
        return str(path.resolve())


def only_file(directory: Path, pattern: str, description: str) -> Path:
    matches = sorted(directory.rglob(pattern))
    if len(matches) != 1:
        raise FileNotFoundError(
            f"Expected one {description} matching {pattern}; found {len(matches)}"
        )
    return matches[0]


def parse_dr122_target_info(path: Path, tic_id: int) -> dict[str, Any]:
    """Read one target row from the official DR122 supplemental table."""
    minimum: int | None = None
    maximum: int | None = None
    matching_rows: list[list[str]] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith("#"):
            minimum_match = re.search(r"minSector:\s*(\d+)", line)
            maximum_match = re.search(r"maxSector:\s*(\d+)", line)
            if minimum_match:
                minimum = int(minimum_match.group(1))
            if maximum_match:
                maximum = int(maximum_match.group(1))
            continue
        fields = line.split()
        if fields and int(fields[0]) == tic_id:
            matching_rows.append(fields)
    if minimum is None or maximum is None:
        raise ValueError("DR122 target table does not define its sector range")
    if len(matching_rows) != 1:
        raise ValueError(
            f"Expected one DR122 target row for TIC {tic_id}; "
            f"found {len(matching_rows)}"
        )
    fields = matching_rows[0]
    if len(fields) != 4:
        raise ValueError(f"Unexpected DR122 target row: {' '.join(fields)}")
    vector = fields[1]
    expected_length = maximum - minimum + 1
    if len(vector) != expected_length or set(vector) - {"0", "1"}:
        raise ValueError("Invalid sector vector in DR122 target table")
    sectors = tuple(
        minimum + index for index, value in enumerate(vector) if value == "1"
    )
    return {
        "target_info_search_first_sector": minimum,
        "target_info_search_last_sector": maximum,
        "target_info_sector_vector": vector,
        "target_info_contributing_sectors": ";".join(
            str(value) for value in sectors
        ),
        "target_info_made_tce": bool(int(fields[2])),
        "target_info_completed_dv": bool(int(fields[3])),
    }


def parse_exomast_info(path: Path) -> dict[str, Any]:
    """Read the combined TCE header returned by the official Exo.MAST API."""
    payload = json.loads(path.read_text(encoding="utf-8"))
    header = payload.get("DV Data Header")
    if not isinstance(header, dict):
        raise ValueError("Exo.MAST response has no DV Data Header")
    return {
        "exomast_tic_id": int(header["TICID"]),
        "exomast_scope": str(header["SECTORS"]),
        "exomast_period_days": float(header["TPERIOD"]),
        "exomast_epoch_btjd": float(header["TEPOCH"]),
        "exomast_depth_ppt": float(header["TDEPTH"]) / 1000.0,
        "exomast_duration_hours": float(header["TDUR"]),
        "exomast_fit_snr": float(header["TSNR"]),
        "exomast_maximum_multiple_event_statistic": float(header["MAXMES"]),
        "exomast_observed_transits": int(header["NTRANS"]),
        "exomast_cadence_minutes": float(header["TIMEDEL"]) * 1440.0,
        "exomast_time_start_btjd": float(header["TSTART"]),
        "exomast_time_stop_btjd": float(header["TSTOP"]),
    }


def nullable_float_array(rows: list[dict[str, Any]], key: str) -> np.ndarray:
    return np.asarray(
        [np.nan if row.get(key) is None else float(row[key]) for row in rows],
        dtype=float,
    )


def parse_exomast_table(path: Path) -> dict[str, Any]:
    """Summarize the target's archived combined DV time-series table."""
    payload = json.loads(path.read_text(encoding="utf-8"))
    rows = payload.get("data")
    if not isinstance(rows, list) or not rows:
        raise ValueError("Exo.MAST response has no DV time-series rows")
    time = nullable_float_array(rows, "TIME")
    detrended = nullable_float_array(rows, "LC_DETREND")
    model = nullable_float_array(rows, "MODEL_INIT")
    finite_time = time[np.isfinite(time)]
    finite_model = model[np.isfinite(model)]
    scopes = sorted({str(row.get("SECTORS")) for row in rows})
    return {
        "exomast_table_rows": int(len(rows)),
        "exomast_table_finite_detrended_rows": int(
            np.count_nonzero(np.isfinite(detrended))
        ),
        "exomast_table_time_start_btjd": (
            float(np.min(finite_time)) if finite_time.size else None
        ),
        "exomast_table_time_stop_btjd": (
            float(np.max(finite_time)) if finite_time.size else None
        ),
        "exomast_table_model_depth_ppt": float(-np.min(finite_model) * 1000.0),
        "exomast_table_scopes": ";".join(
            value for value in scopes if value != "None"
        ),
    }


def parse_multisector_difference_images(path: Path) -> pd.DataFrame:
    """Read the per-sector difference-image entries from a combined XML."""
    root = ET.parse(path).getroot()
    planet = required(root, "dv:planetResults")
    rows: list[dict[str, Any]] = []
    for result in planet.findall("dv:differenceImageResults", DV_NAMESPACE):
        quality = required(result, "dv:qualityMetric")
        rows.append(
            {
                "sector": int(result.attrib["sector"]),
                "difference_image_transits": int(
                    result.attrib["numberOfTransits"]
                ),
                "difference_image_quality": float_attribute(quality, "value"),
                "cadences_in_transit": int(
                    result.attrib["numberOfCadencesInTransit"]
                ),
                "cadences_out_of_transit": int(
                    result.attrib["numberOfCadencesOutOfTransit"]
                ),
            }
        )
    if not rows:
        raise ValueError("Combined report has no difference-image entries")
    return pd.DataFrame(rows).sort_values("sector").reset_index(drop=True)


def exomast_info_matches_dvt(
    info: dict[str, Any], dvt: dict[str, Any]
) -> bool:
    return bool(
        info["exomast_tic_id"] == dvt["tic_id"]
        and math.isclose(
            info["exomast_period_days"],
            dvt["dvt_period_days"],
            rel_tol=0.0,
            abs_tol=1e-10,
        )
        and math.isclose(
            info["exomast_epoch_btjd"],
            dvt["dvt_epoch_btjd"],
            rel_tol=0.0,
            abs_tol=1e-9,
        )
        and math.isclose(
            info["exomast_depth_ppt"],
            dvt["dvt_depth_ppt"],
            rel_tol=0.0,
            abs_tol=1e-9,
        )
        and math.isclose(
            info["exomast_duration_hours"],
            dvt["dvt_duration_hours"],
            rel_tol=0.0,
            abs_tol=1e-9,
        )
    )


def find_multisector_products(
    data_dir: Path,
    exomast_dir: Path,
    release_notes_dir: Path,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Verify and summarize the official DR122 combined report set."""
    xml_path = only_file(data_dir, "*_dvr.xml", "combined DVR XML")
    dvt_path = only_file(data_dir, "*_dvt.fits", "combined DVT FITS")
    dvr_pdf = only_file(data_dir, "*_dvr.pdf", "combined full DVR PDF")
    dvm_pdf = only_file(data_dir, "*_dvm.pdf", "combined mini-report PDF")
    dvs_pdf = only_file(data_dir, "*_dvs.pdf", "combined summary PDF")
    info_path = exomast_dir / f"{MULTISECTOR_SCOPE}_tce1_info.json"
    table_path = exomast_dir / f"{MULTISECTOR_SCOPE}_tce1_table.json"
    tce_list_path = exomast_dir / "tce_list.json"
    products_path = exomast_dir / f"{MULTISECTOR_SCOPE}_vetting_products.json"
    target_info_path = (
        release_notes_dir
        / "tess_multisector_14_86_drn122_targetinfo_v01.txt"
    )
    release_note_path = (
        release_notes_dir / "tess_multisector_14_86_drn122_v01.pdf"
    )
    for path in (
        info_path,
        table_path,
        tce_list_path,
        products_path,
        target_info_path,
        release_note_path,
    ):
        if not path.exists():
            raise FileNotFoundError(path)

    xml_row = parse_dv_xml(xml_path)
    dvt_row = parse_dvt_header(dvt_path)
    exomast_info = parse_exomast_info(info_path)
    exomast_table = parse_exomast_table(table_path)
    target_info = parse_dr122_target_info(target_info_path, int(xml_row["tic_id"]))
    difference_images = parse_multisector_difference_images(xml_path)

    expected_sector_text = ";".join(str(value) for value in EXPECTED_SECTORS)
    checks = {
        "xml_dvt_fit": fit_values_match(xml_row, dvt_row),
        "xml_difference_sectors": tuple(difference_images["sector"])
        == EXPECTED_SECTORS,
        "dvt_sector_vector": dvt_row["dvt_sectors_used"]
        == expected_sector_text,
        "target_info_sectors": target_info[
            "target_info_contributing_sectors"
        ]
        == expected_sector_text,
        "target_info_completed": target_info["target_info_made_tce"]
        and target_info["target_info_completed_dv"],
        "exomast_info_dvt": exomast_info_matches_dvt(exomast_info, dvt_row),
        "exomast_scope": exomast_info["exomast_scope"] == MULTISECTOR_SCOPE,
        "table_rows": exomast_table["exomast_table_rows"]
        == dvt_row["dvt_time_rows"],
        "table_finite_rows": exomast_table[
            "exomast_table_finite_detrended_rows"
        ]
        == dvt_row["dvt_finite_detrended_rows"],
        "table_model_depth": math.isclose(
            exomast_table["exomast_table_model_depth_ppt"],
            dvt_row["dvt_depth_ppt"],
            rel_tol=0.0,
            abs_tol=1e-5,
        ),
    }
    failed = [name for name, passed in checks.items() if not passed]
    if failed:
        raise ValueError(
            "Combined SPOC product consistency check failed: "
            + ", ".join(failed)
        )

    tce_list = json.loads(tce_list_path.read_text(encoding="utf-8"))
    expected_tce = f"{MULTISECTOR_SCOPE}:TCE_1"
    if expected_tce not in tce_list.get("TCE", []):
        raise ValueError(f"Exo.MAST TCE list does not contain {expected_tce}")
    products = json.loads(products_path.read_text(encoding="utf-8"))
    product_names = {str(item["productSubGroupDescription"]) for item in products}
    if product_names != {"DVR", "DVM", "DVS", "DVT"}:
        raise ValueError(
            "Combined report product set is incomplete: "
            + ", ".join(sorted(product_names))
        )

    row = dict(xml_row)
    row.pop("sector", None)
    row.pop("difference_image_quality", None)
    row.pop("difference_image_transits", None)
    dvt_without_sector = dict(dvt_row)
    dvt_without_sector.pop("sector", None)
    row.update(dvt_without_sector)
    row.update(exomast_info)
    row.update(exomast_table)
    row.update(target_info)
    row.update(
        {
            "report_scope": MULTISECTOR_SCOPE,
            "contributing_sectors": expected_sector_text,
            "difference_image_sectors": expected_sector_text,
            "difference_image_transits_total": int(
                difference_images["difference_image_transits"].sum()
            ),
            "all_machine_checks_pass": True,
            "dvr_xml_file": relative(xml_path),
            "dvt_fits_file": relative(dvt_path),
            "dvr_pdf_file": relative(dvr_pdf),
            "dvm_pdf_file": relative(dvm_pdf),
            "dvs_pdf_file": relative(dvs_pdf),
            "exomast_info_file": relative(info_path),
            "exomast_table_file": relative(table_path),
            "exomast_tce_list_file": relative(tce_list_path),
            "exomast_product_list_file": relative(products_path),
            "dr122_target_info_file": relative(target_info_path),
            "dr122_release_note_file": relative(release_note_path),
        }
    )
    difference_images.insert(0, "report_scope", MULTISECTOR_SCOPE)
    difference_images["dvr_xml_file"] = relative(xml_path)
    return pd.DataFrame([row]), difference_images


def find_official_products(data_dir: Path) -> pd.DataFrame:
    xml_paths = sorted(data_dir.rglob("*_dvr.xml"))
    dvt_paths = sorted(data_dir.rglob("*_dvt.fits"))
    if len(xml_paths) != len(EXPECTED_SECTORS):
        raise FileNotFoundError(
            f"Expected {len(EXPECTED_SECTORS)} DVR XML files; found {len(xml_paths)}"
        )
    if len(dvt_paths) != len(EXPECTED_SECTORS):
        raise FileNotFoundError(
            f"Expected {len(EXPECTED_SECTORS)} DVT FITS files; found {len(dvt_paths)}"
        )

    dvt_by_sector: dict[int, tuple[Path, dict[str, Any]]] = {}
    for path in dvt_paths:
        row = parse_dvt_header(path)
        dvt_by_sector[int(row["sector"])] = (path, row)

    rows: list[dict[str, Any]] = []
    for xml_path in xml_paths:
        row = parse_dv_xml(xml_path)
        sector = int(row["sector"])
        if sector not in dvt_by_sector:
            raise FileNotFoundError(f"No DVT FITS product found for Sector {sector}")
        dvt_path, dvt_row = dvt_by_sector[sector]
        if not dvt_matches_xml(row, dvt_row):
            raise ValueError(
                f"Sector {sector} DVR XML and DVT FITS summaries do not agree"
            )

        report_pdf = xml_path.with_suffix(".pdf")
        summary_pdfs = sorted(xml_path.parent.glob("*_dvs.pdf"))
        if not report_pdf.exists() or len(summary_pdfs) != 1:
            raise FileNotFoundError(
                f"Expected one DVR and one DVS PDF beside Sector {sector} XML"
            )

        row.update(dvt_row)
        row.update(
            {
                "dvt_matches_xml_fit": True,
                "dvr_xml_file": relative(xml_path),
                "dvt_fits_file": relative(dvt_path),
                "dvr_pdf_file": relative(report_pdf),
                "dvs_pdf_file": relative(summary_pdfs[0]),
            }
        )
        rows.append(row)

    table = pd.DataFrame(rows).sort_values("sector").reset_index(drop=True)
    found_sectors = tuple(int(value) for value in table["sector"])
    if found_sectors != EXPECTED_SECTORS:
        raise ValueError(
            f"Expected official report sectors {EXPECTED_SECTORS}; found {found_sectors}"
        )
    return table


def one_row(table: pd.DataFrame, description: str, **filters: Any) -> pd.Series:
    selected = table
    for column, value in filters.items():
        selected = selected[selected[column] == value]
    if len(selected) != 1:
        raise ValueError(
            f"Expected one {description} row for {filters}; found {len(selected)}"
        )
    return selected.iloc[0]


def build_project_comparison(
    official: pd.DataFrame, tess_analysis_dir: Path, pixel_dir: Path
) -> pd.DataFrame:
    sectors = pd.read_csv(tess_analysis_dir / "sector_measurements.csv")
    periods = pd.read_csv(tess_analysis_dir / "period_search.csv")
    pipelines = pd.read_csv(tess_analysis_dir / "pipeline_comparison.csv")
    centroids = pd.read_csv(pixel_dir / "difference_image_localization.csv")

    rows: list[dict[str, Any]] = []
    for _, dv in official.iterrows():
        sector = int(dv["sector"])
        qlp = one_row(sectors, "QLP sector measurement", sector=sector, pipeline="QLP")
        period = one_row(periods, "QLP period search", sector=sector)
        spoc = one_row(
            pipelines,
            "corrected SPOC box fit",
            sector=sector,
            pipeline="SPOC",
            series="corrected",
        )
        centroid = one_row(centroids, "project centroid", sector=sector)

        rows.append(
            {
                "sector": sector,
                "project_qlp_period_search_days": float(period["best_period_days"]),
                "official_spoc_fit_period_days": float(dv["fit_period_days"]),
                "period_difference_minutes": (
                    float(period["best_period_days"]) - float(dv["fit_period_days"])
                )
                * 1440.0,
                "project_qlp_depth_ppt": float(qlp["depth_ppt"]),
                "project_qlp_depth_error_ppt": float(qlp["depth_error_ppt"]),
                "project_spoc_box_depth_ppt": float(
                    spoc["depth_ppt_fixed_2_004h"]
                ),
                "project_spoc_box_depth_error_ppt": float(
                    spoc["depth_error_ppt"]
                ),
                "official_spoc_fit_depth_ppt": float(dv["fit_depth_ppt"]),
                "official_spoc_fit_depth_error_ppt": float(
                    dv["fit_depth_error_ppt"]
                ),
                "qlp_minus_official_depth_ppt": float(qlp["depth_ppt"])
                - float(dv["fit_depth_ppt"]),
                "spoc_box_minus_official_depth_ppt": float(
                    spoc["depth_ppt_fixed_2_004h"]
                )
                - float(dv["fit_depth_ppt"]),
                "project_qlp_duration_hours": float(qlp["duration_hours"]),
                "official_spoc_fit_duration_hours": float(
                    dv["fit_duration_hours"]
                ),
                "project_qlp_odd_even_difference_sigma": abs(
                    float(qlp["odd_even_difference_sigma"])
                ),
                "official_odd_even_difference_sigma": float(
                    dv["odd_even_difference_sigma"]
                ),
                "project_qlp_phase_0_5_snr": abs(float(qlp["phase_0_5_snr"])),
                "official_weak_secondary_mes": float(dv["weak_secondary_mes"]),
                "project_centroid_offset_arcsec": float(
                    centroid["centroid_offset_arcsec"]
                ),
                "project_centroid_error_arcsec": float(
                    centroid["bootstrap_68pct_centroid_error_arcsec"]
                ),
                "official_centroid_offset_arcsec": float(
                    dv["centroid_to_tic_arcsec"]
                ),
                "official_centroid_error_arcsec": float(
                    dv["centroid_to_tic_error_arcsec"]
                ),
                "inside_one_tess_pixel": bool(
                    centroid["consistent_with_target_within_one_pixel"]
                )
                and float(dv["centroid_to_tic_arcsec"]) < 21.0,
                "comparison_limit": (
                    "Same TESS observations; this checks methods and is not "
                    "independent confirmation. Centroid error methods differ."
                ),
            }
        )
    return pd.DataFrame(rows)


def save_comparison_figure(comparison: pd.DataFrame, output_dir: Path) -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 10,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "svg.fonttype": "none",
            "svg.hashsalt": "toi3505-spoc-report-comparison",
        }
    )
    sectors = comparison["sector"].to_numpy(dtype=int)
    x = np.arange(len(sectors), dtype=float)

    fig, axes = plt.subplots(2, 1, figsize=(8.0, 8.2), sharex=True)
    depth_series = (
        (-0.22, "project_qlp_depth_ppt", "project_qlp_depth_error_ppt",
         "QLP box fit", "#2166ac", "o"),
        (0.0, "project_spoc_box_depth_ppt", "project_spoc_box_depth_error_ppt",
         "SPOC light curve, box fit", "#4d9221", "s"),
        (0.22, "official_spoc_fit_depth_ppt", "official_spoc_fit_depth_error_ppt",
         "SPOC report transit fit", "#b2182b", "D"),
    )
    for offset, value, error, label, color, marker in depth_series:
        axes[0].errorbar(
            x + offset,
            comparison[value],
            yerr=comparison[error],
            fmt=marker,
            color=color,
            markersize=7,
            capsize=4,
            linewidth=1.4,
            label=label,
        )
    axes[0].set_ylabel("Transit depth (ppt)")
    axes[0].set_title("Depth estimates")
    axes[0].legend(frameon=False, loc="best")
    axes[0].grid(axis="y", color="#dddddd", linewidth=0.8)

    centroid_series = (
        (-0.11, "project_centroid_offset_arcsec", "project_centroid_error_arcsec",
         "Project difference image", "#2166ac", "o"),
        (0.11, "official_centroid_offset_arcsec", "official_centroid_error_arcsec",
         "SPOC report", "#b2182b", "D"),
    )
    for offset, value, error, label, color, marker in centroid_series:
        axes[1].errorbar(
            x + offset,
            comparison[value],
            yerr=comparison[error],
            fmt=marker,
            color=color,
            markersize=7,
            capsize=4,
            linewidth=1.4,
            label=label,
        )
    axes[1].set_ylabel("Centroid offset from target (arcsec)")
    axes[1].set_xlabel("TESS sector")
    axes[1].set_xticks(x, [str(value) for value in sectors])
    axes[1].set_ylim(bottom=0)
    axes[1].set_title("Difference-image location")
    axes[1].legend(frameon=False, loc="best")
    axes[1].grid(axis="y", color="#dddddd", linewidth=0.8)
    axes[1].text(
        0.5,
        0.04,
        "Both methods remain within one TESS pixel (about 21 arcsec).\n"
        "The two error estimates use different methods.",
        transform=axes[1].transAxes,
        ha="center",
        va="bottom",
        fontsize=9,
    )

    fig.suptitle("TOI-3505.01: comparison with SPOC reports", fontsize=14)
    fig.text(
        0.5,
        0.01,
        "These measurements use the same TESS observations; agreement is not "
        "independent confirmation.",
        ha="center",
        va="bottom",
        fontsize=9,
    )
    fig.tight_layout(rect=(0, 0.04, 1, 0.97))
    for suffix in ("png", "svg"):
        fig.savefig(
            output_dir / f"01_spoc_report_comparison.{suffix}",
            dpi=220,
            bbox_inches="tight",
            metadata={"Creator": "TOI-3505.01 research code"},
        )
    plt.close(fig)


def json_records(table: pd.DataFrame) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    for row in table.to_dict(orient="records"):
        clean: dict[str, Any] = {}
        for key, value in row.items():
            if isinstance(value, np.generic):
                value = value.item()
            if value is None or (isinstance(value, float) and not np.isfinite(value)):
                clean[key] = None
            else:
                clean[key] = value
        records.append(clean)
    return records


def write_readme(
    output_dir: Path,
    official: pd.DataFrame,
    comparison: pd.DataFrame,
    multisector: pd.DataFrame,
    multisector_images: pd.DataFrame,
) -> None:
    s54 = one_row(official, "official Sector 54 report", sector=54)
    s81 = one_row(official, "official Sector 81 report", sector=81)
    combined = multisector.iloc[0]
    lines = [
        "# TOI-3505.01 SPOC report comparison",
        "",
        "This folder compares the official SPOC Data Validation products with the "
        "project's simpler TESS measurements. Official reports were found for "
        "Sectors 54 and 81. The MAST search did not return a SPOC target "
        "observation, and therefore no matching per-sector report, for Sectors "
        "14 or 41. SPOC later produced a combined report from the Sector 54 and "
        "81 observations.",
        "",
        "## Main results",
        "",
        f"- Sector 54 SPOC fit: {s54['fit_period_days']:.6f} days, "
        f"{s54['fit_depth_ppt']:.3f} ± {s54['fit_depth_error_ppt']:.3f} ppt, "
        f"duration {s54['fit_duration_hours']:.2f} hours.",
        f"- Sector 81 SPOC fit: {s81['fit_period_days']:.6f} days, "
        f"{s81['fit_depth_ppt']:.3f} ± {s81['fit_depth_error_ppt']:.3f} ppt, "
        f"duration {s81['fit_duration_hours']:.2f} hours.",
        f"- Combined Sector 54 and 81 SPOC fit: "
        f"{combined['fit_period_days']:.8f} days, "
        f"{combined['fit_depth_ppt']:.4f} ± "
        f"{combined['fit_depth_error_ppt']:.4f} ppt, duration "
        f"{combined['fit_duration_hours']:.3f} ± "
        f"{combined['fit_duration_error_hours']:.3f} hours, using "
        f"{int(combined['observed_transits'])} observed transits.",
        "- The combined run is labeled s0014-s0086 because that is the search "
        "range. The DR122 target table, DVT sector vector, and XML difference "
        "images all agree that only Sectors 54 and 81 contributed to this "
        "target. Its time series uses 10-minute bins.",
        f"- Odd/even differences are {s54['odd_even_difference_sigma']:.2f} "
        f"and {s81['odd_even_difference_sigma']:.2f} sigma.",
        f"- The strongest secondary-event statistics are "
        f"{s54['weak_secondary_mes']:.2f} and "
        f"{s81['weak_secondary_mes']:.2f}.",
        f"- SPOC's centroid offsets from the TIC position are "
        f"{s54['centroid_to_tic_arcsec']:.2f} ± "
        f"{s54['centroid_to_tic_error_arcsec']:.2f} arcsec and "
        f"{s81['centroid_to_tic_arcsec']:.2f} ± "
        f"{s81['centroid_to_tic_error_arcsec']:.2f} arcsec.",
        f"- The combined report gives an odd/even difference of "
        f"{combined['odd_even_difference_sigma']:.2f} sigma, a strongest "
        f"secondary statistic of {combined['weak_secondary_mes']:.2f}, and a "
        f"mean TIC-position offset of {combined['centroid_to_tic_arcsec']:.2f} "
        f"± {combined['centroid_to_tic_error_arcsec']:.2f} arcsec.",
        "- Manual review of Appendix B in all three full reports found the statement "
        "\"This target did not trigger any alerts.\"",
        "",
        "The reports support the limited conclusion that SPOC recovered the same "
        "signal and did not flag a strong odd/even difference, secondary event, "
        "or large centroid displacement in these observations. They do not prove "
        "that the signal is planetary. The project and SPOC calculations reuse "
        "the same TESS observations, and TESS cannot resolve the known "
        "0.517-arcsec companion.",
        "",
        "## Files",
        "",
        "- `official_dv_metrics.csv`: selected values read from the SPOC XML, "
        "with DVT FITS header cross-checks.",
        "- `comparison_with_project.csv`: QLP and project SPOC box fits beside "
        "the official SPOC transit fit and centroid result.",
        "- `official_multisector_tce.csv`: the combined fit with checks against "
        "the DVT FITS, Exo.MAST table, DR122 target table, and product list.",
        "- `official_multisector_difference_images.csv`: the Sector 54 and 81 "
        "difference-image entries inside the combined report.",
        "- `analysis_summary.json`: the same results in one machine-readable "
        "record.",
        "- `01_spoc_report_comparison.png` and `.svg`: compact depth and "
        "centroid comparison.",
        "",
        "The odd/even value in sigma is the square root of the test statistic, "
        "matching the full report table. SPOC centroid uncertainties include "
        "the report's systematic error floor; the project errors come from "
        "event resampling, so their error bars are not interchangeable.",
        "",
        "## Official references",
        "",
        "- [TESS data products](https://heasarc.gsfc.nasa.gov/docs/tess/data-products.html)",
        "- [MAST TESS archive](https://archive.stsci.edu/missions-and-data/tess)",
        "- [Exo.MAST light-curve services](https://exo.mast.stsci.edu/docs/dvdata_ws.html)",
        "- [TESS Data Release 122](https://archive.stsci.edu/missions/tess/doc/tess_drn/tess_multisector_14_86_drn122_v01.pdf)",
        "",
        f"Rows written: {len(official)} official report summaries and "
        f"{len(comparison)} project comparisons, plus {len(multisector)} "
        f"combined TCE row and {len(multisector_images)} combined-report "
        "difference-image rows.",
        "",
    ]
    (output_dir / "README.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    args = parse_args()
    data_dir = args.data_dir.resolve()
    multisector_data_dir = args.multisector_data_dir.resolve()
    exomast_dir = args.exomast_dir.resolve()
    release_notes_dir = args.release_notes_dir.resolve()
    tess_analysis_dir = args.tess_analysis_dir.resolve()
    pixel_dir = args.pixel_dir.resolve()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    official = find_official_products(data_dir)
    multisector, multisector_images = find_multisector_products(
        multisector_data_dir,
        exomast_dir=exomast_dir,
        release_notes_dir=release_notes_dir,
    )
    comparison = build_project_comparison(
        official, tess_analysis_dir=tess_analysis_dir, pixel_dir=pixel_dir
    )
    official.to_csv(output_dir / "official_dv_metrics.csv", index=False)
    comparison.to_csv(output_dir / "comparison_with_project.csv", index=False)
    multisector.to_csv(output_dir / "official_multisector_tce.csv", index=False)
    multisector_images.to_csv(
        output_dir / "official_multisector_difference_images.csv", index=False
    )
    save_comparison_figure(comparison, output_dir)
    write_readme(
        output_dir,
        official,
        comparison,
        multisector,
        multisector_images,
    )

    summary = {
        "target": "TOI-3505.01",
        "tic_id": 390988385,
        "report_sectors": list(EXPECTED_SECTORS),
        "official_metrics": json_records(official),
        "official_multisector_tce": json_records(multisector),
        "official_multisector_difference_images": json_records(
            multisector_images
        ),
        "project_comparison": json_records(comparison),
        "interpretation": [
            "SPOC recovered the same approximately 2.915-day transit signal.",
            "No strong odd/even or secondary-event diagnostic appears in any of the three reports.",
            "The combined report uses only Sectors 54 and 81, despite the s0014-s0086 search label.",
            "Both SPOC centroid offsets are below one TESS pixel.",
            "These are same-data method checks, not independent confirmation.",
            "The 0.517-arcsec companion remains unresolved by TESS.",
        ],
        "manual_pdf_review": {
            "per_sector_summary_pages": [1],
            "per_sector_full_report_pdf_pages": [8, 24, 38],
            "combined_summary_pdf_page": 1,
            "combined_full_report_pdf_pages": [7, 8, 27, 41],
            "dr122_release_note_pdf_pages": [4, 5, 7],
            "result": "All three report appendices state that the target did not trigger any alerts.",
        },
    }
    (output_dir / "analysis_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(f"Saved official report comparison to {output_dir}")
    for _, row in official.iterrows():
        print(
            f"Sector {int(row['sector'])}: "
            f"P={row['fit_period_days']:.6f} d, "
            f"depth={row['fit_depth_ppt']:.3f}±"
            f"{row['fit_depth_error_ppt']:.3f} ppt, "
            f"TIC centroid offset={row['centroid_to_tic_arcsec']:.2f} arcsec"
        )
    combined = multisector.iloc[0]
    print(
        "Combined Sectors 54 and 81: "
        f"P={combined['fit_period_days']:.8f} d, "
        f"depth={combined['fit_depth_ppt']:.4f}±"
        f"{combined['fit_depth_error_ppt']:.4f} ppt, "
        f"TIC centroid offset={combined['centroid_to_tic_arcsec']:.2f} arcsec"
    )


if __name__ == "__main__":
    main()
