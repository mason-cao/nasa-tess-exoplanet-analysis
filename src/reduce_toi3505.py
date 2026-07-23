"""Calibrate the GMU TOI-3505.01 image sequence for AstroImageJ.

The program-provided files are split across six ZIP archives.  This pipeline
reads them directly, builds median master darks at the two matching exposure
times, builds a gradient-corrected and normalized median R-band master flat,
and writes 32-bit calibrated science FITS images using

    calibrated = (raw science - master 50 s dark) / master R flat

No bias is subtracted because the matching darks include the camera bias, as
specified in the Schar AstroImageJ tutorial.  No exposure scaling, nonlinearity
correction, cosmic-ray filtering, plate solving, or photometry is performed.
Those are deliberately separate and auditable stages.

The flat workflow mirrors the requested AstroImageJ Data Processor settings:
each 3.5 s R flat is dark-corrected, a fitted illumination plane is divided
out ("Remove Gradient"), each flat is normalized, and the stack is combined by
the median.  The output FITS files preserve the source timing metadata and add
explicit calibration HISTORY records and checksums.
"""

from __future__ import annotations

import argparse
import io
import json
import logging
import math
import os
import re
import tempfile
import time
import zipfile
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path
from typing import Iterable

os.environ.setdefault(
    "MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "nasa-tess-matplotlib-cache")
)

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.visualization import AsinhStretch, ImageNormalize, ZScaleInterval


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_ARCHIVE_DIR = ROOT / "data_and_lectures"
DEFAULT_OUTPUT_ROOT = ROOT / "data" / "ground" / "toi3505"
DEFAULT_DIAGNOSTICS_DIR = ROOT / "outputs" / "toi3505_reduction"

ARCHIVE_GLOB = "TOI_3505.01-20260714T190730Z-1-00*.zip"
SCIENCE_RE = re.compile(
    r"(?:^|/)TOI_3505\.01_50\.000s_R-(\d{4})(?:\(1\))?\.fits$"
)
DARK_50_RE = re.compile(r"(?:^|/)Dark_50\.000s-(\d{4})\.fits$")
DARK_3_5_RE = re.compile(r"(?:^|/)Dark_3\.500s-(\d{4})\.fits$")
FLAT_RE = re.compile(r"(?:^|/)Flat_3\.500s_R-(\d{4})-final\.fits$")

EXPECTED_SHAPE = (4096, 4096)
EXPECTED_SCIENCE_COUNT = 283
EXPECTED_CALIBRATION_COUNT = 10
SAMPLE_STRIDE = 8
SPOTCHECK_THUMBNAIL_STRIDE = 16
REVIEW_THUMBNAIL_STRIDE = 32
CONTACT_SHEET_FRAMES_PER_PAGE = 42
PIPELINE_VERSION = "1.1.0"


@dataclass(frozen=True)
class MemberRef:
    """A FITS member inside one source ZIP archive."""

    archive: Path
    member: str
    name: str
    sequence: int
    crc: int
    file_size: int


@dataclass
class Inventory:
    archives: list[Path]
    science: list[MemberRef]
    dark_50: list[MemberRef]
    dark_3_5: list[MemberRef]
    flats: list[MemberRef]
    duplicate_science_members: list[dict[str, object]]

    def summary(self) -> dict[str, object]:
        return {
            "archive_count": len(self.archives),
            "archives": [path.name for path in self.archives],
            "science_frames": len(self.science),
            "science_sequence_start": self.science[0].sequence if self.science else None,
            "science_sequence_end": self.science[-1].sequence if self.science else None,
            "dark_50s": len(self.dark_50),
            "dark_3_5s": len(self.dark_3_5),
            "flat_3_5s_r": len(self.flats),
            "duplicate_science_members": self.duplicate_science_members,
        }


def configure_logging(log_path: Path, *, overwrite: bool = False) -> logging.Logger:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    logger = logging.getLogger("toi3505_reduction")
    logger.setLevel(logging.INFO)
    logger.handlers.clear()
    formatter = logging.Formatter("%(asctime)s | %(levelname)s | %(message)s")
    stream = logging.StreamHandler()
    stream.setFormatter(formatter)
    file_handler = logging.FileHandler(
        log_path,
        mode="w" if overwrite else "a",
        encoding="utf-8",
    )
    file_handler.setFormatter(formatter)
    logger.addHandler(stream)
    logger.addHandler(file_handler)
    return logger


def _choose_preferred(candidates: list[MemberRef]) -> MemberRef:
    """Prefer a normally named member over a Finder-style ``(1)`` duplicate."""
    return min(candidates, key=lambda ref: ("(1)" in ref.name, ref.member))


def _collect_category(
    matches: dict[int, list[MemberRef]],
    expected_count: int,
    category: str,
) -> tuple[list[MemberRef], list[dict[str, object]]]:
    selected: list[MemberRef] = []
    duplicates: list[dict[str, object]] = []
    for sequence, candidates in sorted(matches.items()):
        signatures = {(item.crc, item.file_size) for item in candidates}
        if len(signatures) > 1:
            raise RuntimeError(
                f"Conflicting duplicate {category} frame {sequence}: "
                + ", ".join(item.member for item in candidates)
            )
        preferred = _choose_preferred(candidates)
        selected.append(preferred)
        if len(candidates) > 1:
            duplicates.append(
                {
                    "sequence": sequence,
                    "selected": preferred.member,
                    "ignored": [item.member for item in candidates if item != preferred],
                    "crc": preferred.crc,
                    "file_size": preferred.file_size,
                }
            )
    if len(selected) != expected_count:
        raise RuntimeError(
            f"Expected {expected_count} unique {category} frames, found {len(selected)}"
        )
    return selected, duplicates


def discover_inventory(archive_dir: Path = DEFAULT_ARCHIVE_DIR) -> Inventory:
    archives = sorted(archive_dir.glob(ARCHIVE_GLOB))
    if len(archives) != 6:
        raise FileNotFoundError(
            f"Expected six TOI-3505 archives in {archive_dir}, found {len(archives)}"
        )

    science_matches: dict[int, list[MemberRef]] = {}
    dark_50_matches: dict[int, list[MemberRef]] = {}
    dark_3_5_matches: dict[int, list[MemberRef]] = {}
    flat_matches: dict[int, list[MemberRef]] = {}

    patterns = (
        (SCIENCE_RE, science_matches),
        (DARK_50_RE, dark_50_matches),
        (DARK_3_5_RE, dark_3_5_matches),
        (FLAT_RE, flat_matches),
    )
    for archive_path in archives:
        with zipfile.ZipFile(archive_path) as archive:
            for info in archive.infolist():
                for pattern, destination in patterns:
                    match = pattern.search(info.filename)
                    if match:
                        sequence = int(match.group(1))
                        destination.setdefault(sequence, []).append(
                            MemberRef(
                                archive=archive_path,
                                member=info.filename,
                                name=Path(info.filename).name,
                                sequence=sequence,
                                crc=info.CRC,
                                file_size=info.file_size,
                            )
                        )
                        break

    science, duplicates = _collect_category(
        science_matches, EXPECTED_SCIENCE_COUNT, "science"
    )
    dark_50, _ = _collect_category(
        dark_50_matches, EXPECTED_CALIBRATION_COUNT, "50 s dark"
    )
    dark_3_5, _ = _collect_category(
        dark_3_5_matches, EXPECTED_CALIBRATION_COUNT, "3.5 s dark"
    )
    flats, _ = _collect_category(
        flat_matches, EXPECTED_CALIBRATION_COUNT, "3.5 s R flat"
    )
    if [ref.sequence for ref in science] != list(range(1, EXPECTED_SCIENCE_COUNT + 1)):
        raise RuntimeError("Science sequence is not contiguous from 0001 through 0283")

    return Inventory(
        archives=archives,
        science=science,
        dark_50=dark_50,
        dark_3_5=dark_3_5,
        flats=flats,
        duplicate_science_members=duplicates,
    )


def read_member(ref: MemberRef) -> tuple[np.ndarray, fits.Header]:
    """Read a source FITS member fully, which also forces its ZIP CRC check."""
    with zipfile.ZipFile(ref.archive) as archive:
        payload = archive.read(ref.member)
    with fits.open(io.BytesIO(payload), memmap=False, uint=True) as hdul:
        if len(hdul) != 1 or hdul[0].data is None:
            raise RuntimeError(f"Expected one primary image HDU in {ref.member}")
        data = np.asarray(hdul[0].data, dtype=np.float32)
        header = hdul[0].header.copy()
    if data.shape != EXPECTED_SHAPE:
        raise RuntimeError(f"Unexpected shape {data.shape} for {ref.member}")
    return data, header


def robust_mad(values: np.ndarray) -> float:
    values = np.asarray(values, dtype=np.float64)
    finite = values[np.isfinite(values)]
    if not len(finite):
        return math.nan
    center = np.median(finite)
    return float(1.4826 * np.median(np.abs(finite - center)))


def sampled_statistics(data: np.ndarray, stride: int = SAMPLE_STRIDE) -> dict[str, float]:
    sample = np.asarray(data[::stride, ::stride], dtype=np.float64).ravel()
    finite = sample[np.isfinite(sample)]
    if not len(finite):
        return {
            "mean": math.nan,
            "median": math.nan,
            "mad": math.nan,
            "p01": math.nan,
            "p99": math.nan,
            "min": math.nan,
            "max": math.nan,
            "finite_fraction": 0.0,
        }
    p01, p99 = np.percentile(finite, [1.0, 99.0])
    return {
        "mean": float(np.mean(finite)),
        "median": float(np.median(finite)),
        "mad": robust_mad(finite),
        "p01": float(p01),
        "p99": float(p99),
        "min": float(np.min(finite)),
        "max": float(np.max(finite)),
        "finite_fraction": float(len(finite) / len(sample)),
    }


def validate_header(
    header: fits.Header,
    *,
    exposure: float,
    image_type: str,
    filter_name: str | None = None,
) -> None:
    actual_exposure = float(header.get("EXPTIME", math.nan))
    if not np.isclose(actual_exposure, exposure, atol=1e-6):
        raise RuntimeError(f"Expected EXPTIME={exposure}, found {actual_exposure}")
    actual_type = str(header.get("IMAGETYP", "")).strip().lower()
    if image_type.lower() not in actual_type:
        raise RuntimeError(f"Expected {image_type!r} image, found IMAGETYP={actual_type!r}")
    if filter_name is not None:
        actual_filter = str(header.get("FILTER", "")).strip().lower()
        if filter_name.lower() not in actual_filter:
            raise RuntimeError(
                f"Expected filter containing {filter_name!r}, found {actual_filter!r}"
            )


def calibration_header(
    source_header: fits.Header,
    *,
    image_type: str,
    combined: int,
    exposure: float,
    filter_name: str | None,
    history: Iterable[str],
) -> fits.Header:
    header = source_header.copy()
    header["IMAGETYP"] = (image_type, "Calibration product type")
    header["EXPTIME"] = (exposure, "Exposure time represented by master")
    if filter_name is not None:
        header["FILTER"] = (filter_name, "Photometric filter")
    header["NCOMBINE"] = (combined, "Number of median-combined inputs")
    header["COMBMETH"] = ("MEDIAN", "Calibration combination method")
    header["CALVERS"] = (PIPELINE_VERSION, "TOI-3505 reduction pipeline")
    header["CALDATE"] = (
        datetime.now(UTC).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "UTC calibration creation time",
    )
    for entry in history:
        header.add_history(entry)
    return header


def write_fits_atomic(
    path: Path,
    data: np.ndarray,
    header: fits.Header,
    *,
    overwrite: bool,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists() and not overwrite:
        raise FileExistsError(path)
    temporary = path.with_name(f".{path.name}.partial")
    if temporary.exists():
        temporary.unlink()
    fits.PrimaryHDU(data=np.asarray(data, dtype=np.float32), header=header).writeto(
        temporary,
        overwrite=True,
        checksum=True,
        output_verify="fix+warn",
    )
    os.replace(temporary, path)


def median_combine(
    refs: list[MemberRef],
    *,
    expected_exposure: float,
    logger: logging.Logger,
) -> tuple[np.ndarray, fits.Header, list[dict[str, object]]]:
    stack = np.empty((len(refs), *EXPECTED_SHAPE), dtype=np.float32)
    first_header: fits.Header | None = None
    rows: list[dict[str, object]] = []
    for index, ref in enumerate(refs):
        data, header = read_member(ref)
        validate_header(header, exposure=expected_exposure, image_type="dark")
        if first_header is None:
            first_header = header
        stack[index] = data
        stats = sampled_statistics(data)
        rows.append(
            {
                "kind": f"dark_{expected_exposure:g}s",
                "sequence": ref.sequence,
                "archive": ref.archive.name,
                "member": ref.member,
                **stats,
            }
        )
        logger.info("Loaded dark %d/%d: %s", index + 1, len(refs), ref.name)
    master = np.median(stack, axis=0, overwrite_input=True).astype(np.float32)
    del stack
    assert first_header is not None
    return master, first_header, rows


def fit_illumination_plane(
    image: np.ndarray,
    *,
    stride: int = SAMPLE_STRIDE,
    iterations: int = 3,
    sigma: float = 4.0,
) -> tuple[np.ndarray, dict[str, float]]:
    """Fit a robust plane to a flat and divide out only that plane gradient.

    Coordinates are normalized to [-1, 1].  Subsampling makes the fit practical
    for 4096-square frames while retaining more than 250,000 fit pixels.
    Iterative clipping prevents hot pixels and dust shadows from determining the
    illumination slope.  The returned image retains its original mean level.
    """
    sample = np.asarray(image[::stride, ::stride], dtype=np.float64)
    y_indices = np.arange(0, image.shape[0], stride, dtype=np.float64)
    x_indices = np.arange(0, image.shape[1], stride, dtype=np.float64)
    y = -1.0 + 2.0 * y_indices / (image.shape[0] - 1)
    x = -1.0 + 2.0 * x_indices / (image.shape[1] - 1)
    xx, yy = np.meshgrid(x, y)
    z = sample.ravel()
    design = np.column_stack((np.ones(z.size), xx.ravel(), yy.ravel()))
    finite = np.isfinite(z)
    if finite.sum() < 10_000:
        raise RuntimeError("Too few finite pixels to fit flat illumination plane")
    low, high = np.percentile(z[finite], [0.2, 99.8])
    mask = finite & (z >= low) & (z <= high)
    coefficients = np.zeros(3, dtype=np.float64)
    for _ in range(iterations):
        coefficients, *_ = np.linalg.lstsq(design[mask], z[mask], rcond=None)
        residual = z - design @ coefficients
        scatter = robust_mad(residual[mask])
        if not np.isfinite(scatter) or scatter <= 0:
            break
        mask = finite & (np.abs(residual) <= sigma * scatter)

    full_y = np.linspace(-1.0, 1.0, image.shape[0], dtype=np.float32)
    full_x = np.linspace(-1.0, 1.0, image.shape[1], dtype=np.float32)
    plane = (
        coefficients[0]
        + coefficients[1] * full_x[None, :]
        + coefficients[2] * full_y[:, None]
    ).astype(np.float32)
    plane_level = float(np.mean(plane, dtype=np.float64))
    if not np.isfinite(plane_level) or plane_level <= 0 or np.min(plane) <= 0:
        raise RuntimeError("Non-positive illumination plane fitted to flat")
    corrected = image * np.float32(plane_level) / plane
    details = {
        "plane_intercept_adu": float(coefficients[0]),
        "plane_x_coefficient_adu": float(coefficients[1]),
        "plane_y_coefficient_adu": float(coefficients[2]),
        "plane_peak_to_peak_fraction": float(
            (np.max(plane) - np.min(plane)) / plane_level
        ),
        "fit_pixel_count": int(mask.sum()),
        "fit_residual_mad_adu": robust_mad((z - design @ coefficients)[mask]),
    }
    return corrected.astype(np.float32, copy=False), details


def build_master_flat(
    refs: list[MemberRef],
    master_flat_dark: np.ndarray,
    *,
    remove_gradient: bool,
    logger: logging.Logger,
) -> tuple[np.ndarray, fits.Header, list[dict[str, object]]]:
    stack = np.empty((len(refs), *EXPECTED_SHAPE), dtype=np.float32)
    first_header: fits.Header | None = None
    rows: list[dict[str, object]] = []
    for index, ref in enumerate(refs):
        data, header = read_member(ref)
        validate_header(header, exposure=3.5, image_type="light", filter_name="red")
        if first_header is None:
            first_header = header
        calibrated = data - master_flat_dark
        if remove_gradient:
            calibrated, plane_details = fit_illumination_plane(calibrated)
        else:
            plane_details = {
                "plane_intercept_adu": math.nan,
                "plane_x_coefficient_adu": math.nan,
                "plane_y_coefficient_adu": math.nan,
                "plane_peak_to_peak_fraction": math.nan,
                "fit_pixel_count": 0,
                "fit_residual_mad_adu": math.nan,
            }
        finite = np.isfinite(calibrated)
        level = float(np.mean(calibrated[finite], dtype=np.float64))
        if not np.isfinite(level) or level <= 0:
            raise RuntimeError(f"Invalid normalized flat level for {ref.member}: {level}")
        normalized = calibrated / np.float32(level)
        stack[index] = normalized
        rows.append(
            {
                "kind": "flat_3.5s_R",
                "sequence": ref.sequence,
                "archive": ref.archive.name,
                "member": ref.member,
                "normalization_level_adu": level,
                **plane_details,
                **{f"normalized_{key}": value for key, value in sampled_statistics(normalized).items()},
            }
        )
        logger.info("Calibrated flat %d/%d: %s", index + 1, len(refs), ref.name)

    master = np.median(stack, axis=0, overwrite_input=True).astype(np.float32)
    del stack
    finite_positive = np.isfinite(master) & (master > 0)
    if finite_positive.sum() != master.size:
        raise RuntimeError(
            f"Master flat contains {master.size - int(finite_positive.sum())} non-positive/non-finite pixels"
        )
    master /= np.float32(np.mean(master, dtype=np.float64))
    assert first_header is not None
    return master, first_header, rows


def add_science_calibration_history(
    source_header: fits.Header,
    *,
    source_member: str,
    gradient_removed_from_flats: bool,
) -> fits.Header:
    header = source_header.copy()
    header["IMAGETYP"] = ("Reduced Light Frame", "Calibrated science image")
    header["REDUCED"] = (True, "Dark-subtracted and flat-divided")
    header["DARKCOR"] = (True, "Used matching 50 s median master dark")
    header["FLATCOR"] = (True, "Used dark-corrected median R master flat")
    header["GRADREM"] = (
        gradient_removed_from_flats,
        "Illumination plane removed from input flats",
    )
    header["DARKFILE"] = ("mdark_50.000s.fits", "Master dark filename")
    header["FLATFILE"] = ("mflat_R.fits", "Master flat filename")
    header["CALVERS"] = (PIPELINE_VERSION, "TOI-3505 reduction pipeline")
    header["CALDATE"] = (
        datetime.now(UTC).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "UTC calibration creation time",
    )
    header.add_history(f"Source ZIP member: {source_member}")
    header.add_history("Subtracted median of 10 matching 50.000 s dark frames.")
    header.add_history(
        "Divided by median of 10 dark-corrected, gradient-removed, normalized 3.500 s R flats."
        if gradient_removed_from_flats
        else "Divided by median of 10 dark-corrected, normalized 3.500 s R flats."
    )
    header.add_history("No bias scaling, cosmic-ray removal, plate solving, or photometry applied.")
    return header


def output_name(ref: MemberRef) -> str:
    return re.sub(r"(?:\(1\))?\.fits$", "_out.fits", ref.name)


def display_path(path: Path) -> str:
    try:
        return str(path.relative_to(ROOT))
    except ValueError:
        return str(path)


def load_existing_reduced(path: Path, source_member: str) -> tuple[np.ndarray, fits.Header]:
    with fits.open(path, memmap=True, checksum=True) as hdul:
        header = hdul[0].header.copy()
        if not bool(header.get("REDUCED", False)):
            raise RuntimeError(f"Existing output is not marked reduced: {path}")
        history = "\n".join(str(value) for value in header.get("HISTORY", []))
        if source_member not in history:
            raise RuntimeError(f"Existing output provenance does not match {source_member}: {path}")
        data = np.array(hdul[0].data, dtype=np.float32, copy=True)
    return data, header


def reduce_science_frames(
    refs: list[MemberRef],
    master_dark: np.ndarray,
    master_flat: np.ndarray,
    *,
    reduced_dir: Path,
    diagnostics_dir: Path,
    gradient_removed_from_flats: bool,
    overwrite: bool,
    logger: logging.Logger,
) -> tuple[
    pd.DataFrame,
    dict[int, tuple[np.ndarray, np.ndarray]],
    dict[int, np.ndarray],
]:
    reduced_dir.mkdir(parents=True, exist_ok=True)
    diagnostics_dir.mkdir(parents=True, exist_ok=True)
    jsonl_path = diagnostics_dir / "frame_manifest.jsonl"
    rows: list[dict[str, object]] = []
    selected_indices = set(np.linspace(0, len(refs) - 1, 20, dtype=int).tolist())
    selected_indices.update({0, len(refs) // 2, len(refs) - 1})
    spotcheck_thumbnails: dict[int, tuple[np.ndarray, np.ndarray]] = {}
    review_thumbnails: dict[int, np.ndarray] = {}

    if overwrite and jsonl_path.exists():
        jsonl_path.unlink()

    for index, ref in enumerate(refs):
        started = time.perf_counter()
        output_path = reduced_dir / output_name(ref)
        raw: np.ndarray | None = None
        if output_path.exists() and not overwrite:
            calibrated, header = load_existing_reduced(output_path, ref.member)
            status = "resumed-existing"
            if index in selected_indices:
                raw, _ = read_member(ref)
        else:
            raw, header = read_member(ref)
            validate_header(header, exposure=50.0, image_type="light", filter_name="red")
            if str(header.get("OBJECT", "")).strip() != "TOI_3505.01":
                raise RuntimeError(f"Unexpected OBJECT in {ref.member}: {header.get('OBJECT')}")
            calibrated = (raw - master_dark) / master_flat
            calibrated = calibrated.astype(np.float32, copy=False)
            output_header = add_science_calibration_history(
                header,
                source_member=ref.member,
                gradient_removed_from_flats=gradient_removed_from_flats,
            )
            write_fits_atomic(output_path, calibrated, output_header, overwrite=overwrite)
            header = output_header
            status = "written"

        calibrated_stats = sampled_statistics(calibrated)
        if raw is None:
            raw_stats = {
                "mean": math.nan,
                "median": math.nan,
                "mad": math.nan,
                "p01": math.nan,
                "p99": math.nan,
                "min": math.nan,
                "max": math.nan,
                "finite_fraction": math.nan,
            }
            raw_saturation_fraction = math.nan
        else:
            raw_stats = sampled_statistics(raw)
            raw_sample = raw[::SAMPLE_STRIDE, ::SAMPLE_STRIDE]
            raw_saturation_fraction = float(np.mean(raw_sample >= 65500.0))

        calibrated_sample = calibrated[::SAMPLE_STRIDE, ::SAMPLE_STRIDE]
        row = {
            "sequence": ref.sequence,
            "source_archive": ref.archive.name,
            "source_member": ref.member,
            "source_crc32": f"{ref.crc:08x}",
            "output_file": output_path.name,
            "output_relative_path": display_path(output_path),
            "output_size_bytes": output_path.stat().st_size,
            "status": status,
            "date_obs": str(header.get("DATE-OBS", "")),
            "bjd_tdb": float(header.get("BJD_TDB", math.nan)),
            "airmass": float(header.get("AIRMASS", math.nan)),
            "ccd_temp_c": float(header.get("CCD-TEMP", math.nan)),
            "exposure_s": float(header.get("EXPTIME", math.nan)),
            "filter": str(header.get("FILTER", "")).strip(),
            "object": str(header.get("OBJECT", "")).strip(),
            "has_wcs": bool(header.get("CTYPE1") and header.get("CTYPE2")),
            "raw_saturation_fraction_sampled": raw_saturation_fraction,
            "calibrated_negative_fraction_sampled": float(
                np.mean(calibrated_sample < 0)
            ),
            "elapsed_seconds": float(time.perf_counter() - started),
        }
        row.update({f"raw_{key}": value for key, value in raw_stats.items()})
        row.update(
            {f"calibrated_{key}": value for key, value in calibrated_stats.items()}
        )
        rows.append(row)
        with jsonl_path.open("a", encoding="utf-8") as handle:
            handle.write(json.dumps(row, sort_keys=True) + "\n")

        # Keep a small calibrated view of every frame so the generated contact
        # sheets support sequence-wide visual review without reopening 19 GB of
        # FITS products.  Higher-resolution raw/calibrated pairs are retained
        # only for the beginning, midpoint, and end spot checks.
        review_thumbnails[ref.sequence] = calibrated[
            ::REVIEW_THUMBNAIL_STRIDE, ::REVIEW_THUMBNAIL_STRIDE
        ].copy()
        if index in selected_indices:
            if raw is None:
                raw, _ = read_member(ref)
            spotcheck_thumbnails[ref.sequence] = (
                raw[::SPOTCHECK_THUMBNAIL_STRIDE, ::SPOTCHECK_THUMBNAIL_STRIDE].copy(),
                calibrated[
                    ::SPOTCHECK_THUMBNAIL_STRIDE, ::SPOTCHECK_THUMBNAIL_STRIDE
                ].copy(),
            )
        del calibrated
        if raw is not None:
            del raw
        if (index + 1) % 10 == 0 or index == 0 or index + 1 == len(refs):
            logger.info(
                "Reduced %d/%d science frames (latest %s, %.2f s)",
                index + 1,
                len(refs),
                ref.name,
                row["elapsed_seconds"],
            )

    table = pd.DataFrame(rows).sort_values("sequence").reset_index(drop=True)
    return table, spotcheck_thumbnails, review_thumbnails


def robust_zscore(values: pd.Series) -> np.ndarray:
    array = values.to_numpy(dtype=float)
    finite = np.isfinite(array)
    result = np.full_like(array, np.nan, dtype=float)
    if not finite.any():
        return result
    # Robust outlier labels are not meaningful for a tiny --limit smoke test.
    # Non-finite and strongly negative-pixel checks still run independently.
    if finite.sum() < 20:
        result[finite] = 0.0
        return result
    center = float(np.nanmedian(array))
    spread = robust_mad(array)
    if not np.isfinite(spread) or spread <= 0:
        result[finite] = 0.0
    else:
        result[finite] = (array[finite] - center) / spread
    return result


def apply_qa_flags(table: pd.DataFrame) -> pd.DataFrame:
    flagged = table.copy()
    metrics = [
        "calibrated_median",
        "calibrated_mad",
        "calibrated_p99",
        "calibrated_negative_fraction_sampled",
    ]
    for metric in metrics:
        flagged[f"z_{metric}"] = robust_zscore(flagged[metric])

    reasons: list[str] = []
    for _, row in flagged.iterrows():
        frame_reasons: list[str] = []
        if row["calibrated_finite_fraction"] < 1.0:
            frame_reasons.append("nonfinite_pixels")
        if abs(row["z_calibrated_median"]) > 6.0:
            frame_reasons.append("background_level_outlier")
        if abs(row["z_calibrated_mad"]) > 6.0:
            frame_reasons.append("image_scatter_outlier")
        if abs(row["z_calibrated_p99"]) > 6.0:
            frame_reasons.append("bright_level_outlier")
        if row["calibrated_negative_fraction_sampled"] > 0.10:
            frame_reasons.append("many_negative_pixels")
        reasons.append(";".join(frame_reasons))
    flagged["qa_flag"] = [bool(reason) for reason in reasons]
    flagged["qa_reasons"] = reasons
    return flagged


def display_normalization(image: np.ndarray) -> ImageNormalize:
    finite = image[np.isfinite(image)]
    if not len(finite):
        return ImageNormalize(vmin=0.0, vmax=1.0)
    try:
        vmin, vmax = ZScaleInterval(contrast=0.20).get_limits(finite)
    except Exception:
        vmin, vmax = np.percentile(finite, [1, 99.5])
    if not np.isfinite(vmin) or not np.isfinite(vmax) or vmax <= vmin:
        vmin, vmax = float(np.min(finite)), float(np.max(finite) + 1.0)
    return ImageNormalize(vmin=vmin, vmax=vmax, stretch=AsinhStretch(0.08), clip=True)


def plot_master_calibrations(
    master_dark_3_5: np.ndarray,
    master_dark_50: np.ndarray,
    master_flat: np.ndarray,
    output: Path,
) -> None:
    figure, axes = plt.subplots(2, 2, figsize=(12, 10))
    products = [
        (master_dark_3_5[::8, ::8], "Median master dark - 3.5 s", "ADU"),
        (master_dark_50[::8, ::8], "Median master dark - 50 s", "ADU"),
        (master_flat[::8, ::8], "Gradient-corrected median master flat - R", "Relative response"),
    ]
    for axis, (image, title, label) in zip(axes.flat[:3], products):
        finite = image[np.isfinite(image)]
        low, high = np.percentile(finite, [0.5, 99.5])
        rendered = axis.imshow(image, origin="lower", cmap="gray", vmin=low, vmax=high)
        axis.set_title(title)
        axis.set_xlabel("Binned CCD x")
        axis.set_ylabel("Binned CCD y")
        figure.colorbar(rendered, ax=axis, fraction=0.046, pad=0.04, label=label)

    profile_axis = axes[1, 1]
    profile_axis.plot(
        np.nanmedian(master_flat, axis=0),
        label="column median",
        linewidth=1.1,
    )
    profile_axis.plot(
        np.nanmedian(master_flat, axis=1),
        label="row median",
        linewidth=1.1,
    )
    profile_axis.axhline(1.0, color="0.3", linestyle=":")
    profile_axis.set(
        title="Master-flat spatial profiles",
        xlabel="Pixel index",
        ylabel="Relative response",
    )
    profile_axis.legend()
    profile_axis.grid(alpha=0.2)
    figure.suptitle("TOI-3505.01 calibration products", fontsize=16)
    figure.tight_layout()
    figure.savefig(output, dpi=180)
    plt.close(figure)


def plot_raw_calibrated_comparison(
    thumbnails: dict[int, tuple[np.ndarray, np.ndarray]],
    output: Path,
) -> None:
    sequences = [min(thumbnails), sorted(thumbnails)[len(thumbnails) // 2], max(thumbnails)]
    figure, axes = plt.subplots(3, 2, figsize=(10, 14))
    for row_index, sequence in enumerate(sequences):
        raw, calibrated = thumbnails[sequence]
        for column, (image, label) in enumerate(
            ((raw, "raw"), (calibrated, "dark-subtracted / flat-divided"))
        ):
            axis = axes[row_index, column]
            axis.imshow(image, origin="lower", cmap="gray", norm=display_normalization(image))
            axis.set_title(f"Frame {sequence:04d} - {label}")
            axis.set_xlabel("Binned CCD x")
            axis.set_ylabel("Binned CCD y")
    figure.suptitle("TOI-3505.01 reduction spot checks", fontsize=16)
    figure.tight_layout()
    figure.savefig(output, dpi=180)
    plt.close(figure)


def plot_contact_sheets(
    thumbnails: dict[int, np.ndarray],
    table: pd.DataFrame,
    output_dir: Path,
) -> list[Path]:
    sequences = sorted(thumbnails)
    output_dir.mkdir(parents=True, exist_ok=True)
    legacy_output = output_dir / "04_calibrated_contact_sheet.png"
    if legacy_output.exists():
        legacy_output.unlink()
    for stale in output_dir.glob("04_calibrated_contact_sheet_page_*.png"):
        stale.unlink()
    table_by_sequence = table.set_index("sequence")
    outputs: list[Path] = []
    page_count = math.ceil(len(sequences) / CONTACT_SHEET_FRAMES_PER_PAGE)
    for page_index in range(page_count):
        first = page_index * CONTACT_SHEET_FRAMES_PER_PAGE
        page_sequences = sequences[first : first + CONTACT_SHEET_FRAMES_PER_PAGE]
        columns = 7
        rows = math.ceil(len(page_sequences) / columns)
        figure, axes = plt.subplots(rows, columns, figsize=(17.5, 2.6 * rows))
        axes_array = np.atleast_1d(axes).ravel()
        for axis, sequence in zip(axes_array, page_sequences):
            image = thumbnails[sequence]
            axis.imshow(
                image,
                origin="lower",
                cmap="gray",
                norm=display_normalization(image),
            )
            record = table_by_sequence.loc[sequence]
            flag = " | QA" if bool(record["qa_flag"]) else ""
            axis.set_title(
                f"{sequence:04d} | BJD-2459782 {record['bjd_tdb'] - 2459782:.4f}{flag}",
                fontsize=7.5,
                color="#b22222" if flag else "black",
            )
            axis.set_xticks([])
            axis.set_yticks([])
        for axis in axes_array[len(page_sequences) :]:
            axis.axis("off")
        figure.suptitle(
            "TOI-3505.01 calibrated sequence visual review "
            f"({page_index + 1}/{page_count}; frames {page_sequences[0]:04d}-"
            f"{page_sequences[-1]:04d})",
            fontsize=15,
        )
        figure.tight_layout()
        output = output_dir / f"04_calibrated_contact_sheet_page_{page_index + 1:02d}.png"
        figure.savefig(output, dpi=170)
        plt.close(figure)
        outputs.append(output)
    return outputs


def plot_reduction_timeseries(table: pd.DataFrame, output: Path) -> None:
    x = table["bjd_tdb"] - 2459782.0
    figure, axes = plt.subplots(4, 1, figsize=(12, 10), sharex=True)
    metrics = [
        ("calibrated_median", "Sampled image median (ADU)"),
        ("calibrated_mad", "Sampled image robust scatter (ADU)"),
        ("airmass", "Airmass"),
        ("ccd_temp_c", "CCD temperature (C)"),
    ]
    for axis, (column, label) in zip(axes, metrics):
        good = ~table["qa_flag"]
        axis.plot(x, table[column], color="#3566a8", linewidth=0.8, alpha=0.75)
        axis.scatter(x[good], table.loc[good, column], color="#3566a8", s=8)
        axis.scatter(
            x[~good],
            table.loc[~good, column],
            color="#c43d3d",
            s=20,
            label="automatic review candidate",
        )
        axis.set_ylabel(label)
        axis.grid(alpha=0.2)
    if table["qa_flag"].any():
        axes[0].legend(loc="best")
    axes[-1].set_xlabel("BJD_TDB - 2459782")
    figure.suptitle("TOI-3505.01 calibration quality timeline", fontsize=16)
    figure.tight_layout()
    figure.savefig(output, dpi=180)
    plt.close(figure)


def write_readme(
    output_root: Path,
    diagnostics_dir: Path,
    *,
    summary: dict[str, object],
) -> None:
    content = f"""# TOI-3505.01 calibrated image set

Generated by `src/reduce_toi3505.py` (version {PIPELINE_VERSION}).

## Calibration equation

`reduced = (raw science - mdark_50.000s) / mflat_R`

- Master darks: median of 10 matching frames at each exposure time.
- Master flat: each 3.5 s R flat was dark-subtracted, its fitted illumination plane was divided out, it was normalized to mean 1, and the 10 normalized flats were median-combined.
- Output science images: 32-bit floating-point FITS with original timing metadata, calibration HISTORY cards, and FITS checksums.
- Bias subtraction: not applied because the matching dark frames contain the bias signal.
- Not applied: dark exposure scaling, nonlinearity correction, outlier/cosmic-ray filtering, plate solving, alignment, or photometry.

## Contents

- `calibration/mdark_3.500s.fits`
- `calibration/mdark_50.000s.fits`
- `calibration/mflat_R.fits`
- `calibration/master_flat_review_mask.fits` (diagnostic only; not applied)
- `reduced/*_out.fits`
- `{display_path(diagnostics_dir / 'frame_manifest.csv')}`
- `{display_path(diagnostics_dir / 'calibration_inputs.csv')}`
- `{display_path(diagnostics_dir / 'summary.json')}`
- `{display_path(diagnostics_dir / 'reduction.log')}`
- `{display_path(diagnostics_dir / '04_calibrated_contact_sheet_page_*.png')}` (every reduced frame)
- `{display_path(diagnostics_dir / '05_astroimagej_frame_0001.png')}`, `{display_path(diagnostics_dir / '06_astroimagej_frame_0282.png')}`, and `{display_path(diagnostics_dir / '07_astroimagej_frame_0283.png')}` (AstroImageJ compatibility/end-frame checks)
- `{display_path(diagnostics_dir / 'verification.json')}` (independent checksum, timing, and formula validation)

## Status

- Reduced frames: {summary['reduced_science_frames']}
- Automatically flagged for visual review: {summary['qa_flagged_frames']}
- Frames with WCS headers before plate solving: {summary['frames_with_wcs']}
- Calibrated frames represented in visual-review sheets: {summary['contact_sheet_frames']}

Automatic flags are review candidates, not rejected frames. Inspect the sequence in AstroImageJ before deciding whether to exclude any image. The next stage is plate solving or alignment, followed by seeing-profile measurement and multi-aperture photometry.

## Reproduce

From the repository root:

```bash
.venv/bin/python src/reduce_toi3505.py
.venv/bin/python src/verify_toi3505_reduction.py
```
"""
    readme_path = output_root / "README.md"
    readme_path.parent.mkdir(parents=True, exist_ok=True)
    readme_path.write_text(content, encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive-dir", type=Path, default=DEFAULT_ARCHIVE_DIR)
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument("--diagnostics-dir", type=Path, default=DEFAULT_DIAGNOSTICS_DIR)
    parser.add_argument(
        "--inventory-only",
        action="store_true",
        help="Validate archive inventory without reading all image payloads",
    )
    parser.add_argument(
        "--limit",
        type=int,
        help="Reduce only the first N science frames (for a smoke test)",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Replace existing master and reduced FITS products",
    )
    parser.add_argument(
        "--no-gradient-removal",
        action="store_true",
        help="Disable the tutorial-requested illumination-plane correction",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.archive_dir = args.archive_dir.resolve()
    args.output_root = args.output_root.resolve()
    args.diagnostics_dir = args.diagnostics_dir.resolve()
    logger = configure_logging(
        args.diagnostics_dir / "reduction.log",
        overwrite=args.overwrite,
    )
    started = time.perf_counter()

    inventory = discover_inventory(args.archive_dir)
    logger.info("Inventory validated: %s", json.dumps(inventory.summary(), sort_keys=True))
    if args.inventory_only:
        print(json.dumps(inventory.summary(), indent=2))
        return

    calibration_dir = args.output_root / "calibration"
    reduced_dir = args.output_root / "reduced"
    calibration_dir.mkdir(parents=True, exist_ok=True)
    args.diagnostics_dir.mkdir(parents=True, exist_ok=True)
    remove_gradient = not args.no_gradient_removal

    master_dark_3_5_path = calibration_dir / "mdark_3.500s.fits"
    master_dark_50_path = calibration_dir / "mdark_50.000s.fits"
    master_flat_path = calibration_dir / "mflat_R.fits"

    logger.info("Building 3.5 s median master dark")
    master_dark_3_5, dark_3_5_header, dark_3_5_rows = median_combine(
        inventory.dark_3_5,
        expected_exposure=3.5,
        logger=logger,
    )
    header = calibration_header(
        dark_3_5_header,
        image_type="Master Dark",
        combined=len(inventory.dark_3_5),
        exposure=3.5,
        filter_name=None,
        history=["Median of 10 raw 3.500 s dark frames; bias retained."],
    )
    write_fits_atomic(
        master_dark_3_5_path,
        master_dark_3_5,
        header,
        overwrite=True,
    )

    logger.info("Building 50 s median master dark")
    master_dark_50, dark_50_header, dark_50_rows = median_combine(
        inventory.dark_50,
        expected_exposure=50.0,
        logger=logger,
    )
    header = calibration_header(
        dark_50_header,
        image_type="Master Dark",
        combined=len(inventory.dark_50),
        exposure=50.0,
        filter_name=None,
        history=["Median of 10 raw 50.000 s dark frames; bias retained."],
    )
    write_fits_atomic(
        master_dark_50_path,
        master_dark_50,
        header,
        overwrite=True,
    )

    logger.info("Building dark-corrected median R master flat")
    master_flat, flat_header, flat_rows = build_master_flat(
        inventory.flats,
        master_dark_3_5,
        remove_gradient=remove_gradient,
        logger=logger,
    )
    header = calibration_header(
        flat_header,
        image_type="Master Flat",
        combined=len(inventory.flats),
        exposure=3.5,
        filter_name="Red",
        history=[
            "Each input flat was corrected with mdark_3.500s.fits.",
            "A robust illumination plane was divided from each calibrated flat."
            if remove_gradient
            else "Illumination-plane removal was disabled.",
            "Each input was normalized to mean 1 before median combination.",
        ],
    )
    header["DARKFILE"] = ("mdark_3.500s.fits", "Flat-dark master filename")
    header["GRADREM"] = (remove_gradient, "Illumination plane removed")
    write_fits_atomic(master_flat_path, master_flat, header, overwrite=True)

    flat_review_mask = (~np.isfinite(master_flat)) | (master_flat < 0.5) | (master_flat > 1.5)
    mask_header = fits.Header()
    mask_header["IMAGETYP"] = ("Flat Review Mask", "Diagnostic only; was not applied")
    mask_header["BUNIT"] = ("boolean", "1 means review master-flat pixel")
    mask_header["CALVERS"] = PIPELINE_VERSION
    fits.PrimaryHDU(data=flat_review_mask.astype(np.uint8), header=mask_header).writeto(
        calibration_dir / "master_flat_review_mask.fits",
        overwrite=True,
        checksum=True,
    )

    calibration_rows = dark_3_5_rows + dark_50_rows + flat_rows
    pd.DataFrame(calibration_rows).to_csv(
        args.diagnostics_dir / "calibration_inputs.csv", index=False
    )

    science_refs = inventory.science[: args.limit] if args.limit else inventory.science
    logger.info("Reducing %d science frames", len(science_refs))
    frame_table, spotcheck_thumbnails, review_thumbnails = reduce_science_frames(
        science_refs,
        master_dark_50,
        master_flat,
        reduced_dir=reduced_dir,
        diagnostics_dir=args.diagnostics_dir,
        gradient_removed_from_flats=remove_gradient,
        overwrite=args.overwrite,
        logger=logger,
    )
    frame_table = apply_qa_flags(frame_table)
    frame_table.to_csv(args.diagnostics_dir / "frame_manifest.csv", index=False)
    frame_table.loc[frame_table["qa_flag"]].to_csv(
        args.diagnostics_dir / "visual_review_candidates.csv", index=False
    )

    plot_master_calibrations(
        master_dark_3_5,
        master_dark_50,
        master_flat,
        args.diagnostics_dir / "01_master_calibrations.png",
    )
    plot_raw_calibrated_comparison(
        spotcheck_thumbnails,
        args.diagnostics_dir / "02_raw_vs_calibrated.png",
    )
    plot_reduction_timeseries(
        frame_table,
        args.diagnostics_dir / "03_reduction_quality_timeline.png",
    )
    contact_sheet_paths = plot_contact_sheets(
        review_thumbnails,
        frame_table,
        args.diagnostics_dir,
    )

    master_flat_stats = sampled_statistics(master_flat)
    summary = {
        "pipeline_version": PIPELINE_VERSION,
        "created_utc": datetime.now(UTC).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "calibration_formula": "(raw science - median 50 s dark) / median normalized R flat",
        "bias_subtracted": False,
        "dark_exposure_scaled": False,
        "flat_gradient_removed": remove_gradient,
        "cosmic_ray_filter_applied": False,
        "output_bitpix": -32,
        "inventory": inventory.summary(),
        "reduced_science_frames": int(len(frame_table)),
        "qa_flagged_frames": int(frame_table["qa_flag"].sum()),
        "qa_flagged_sequences": frame_table.loc[
            frame_table["qa_flag"], "sequence"
        ].astype(int).tolist(),
        "frames_with_wcs": int(frame_table["has_wcs"].sum()),
        "contact_sheet_frames": int(len(review_thumbnails)),
        "contact_sheet_pages": int(len(contact_sheet_paths)),
        "bjd_tdb_start": float(frame_table["bjd_tdb"].min()),
        "bjd_tdb_end": float(frame_table["bjd_tdb"].max()),
        "master_dark_3_5_stats": sampled_statistics(master_dark_3_5),
        "master_dark_50_stats": sampled_statistics(master_dark_50),
        "master_flat_stats": master_flat_stats,
        "master_flat_review_pixels": int(flat_review_mask.sum()),
        "master_flat_review_fraction": float(flat_review_mask.mean()),
        "reduced_directory_size_bytes": int(
            sum(path.stat().st_size for path in reduced_dir.glob("*_out.fits"))
        ),
        "elapsed_seconds": float(time.perf_counter() - started),
    }
    (args.diagnostics_dir / "summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    write_readme(args.output_root, args.diagnostics_dir, summary=summary)
    logger.info("Reduction complete: %s", json.dumps(summary, sort_keys=True))


if __name__ == "__main__":
    main()
