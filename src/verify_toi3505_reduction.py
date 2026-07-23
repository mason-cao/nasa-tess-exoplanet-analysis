"""Independently verify the complete TOI-3505.01 calibration products."""

from __future__ import annotations

import argparse
import json
import math
import time
from datetime import UTC, datetime
from pathlib import Path

import numpy as np
from astropy.io import fits

import reduce_toi3505 as reduction


def verify_fits_checksum(path: Path) -> tuple[int, int]:
    """Return FITS CHECKSUM and DATASUM verification codes (1 means valid)."""
    with fits.open(path, memmap=True, checksum=False) as hdul:
        if len(hdul) != 1 or hdul[0].data is None:
            raise RuntimeError(f"Expected one primary image HDU: {path}")
        return int(hdul[0].verify_checksum()), int(hdul[0].verify_datasum())


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive-dir", type=Path, default=reduction.DEFAULT_ARCHIVE_DIR)
    parser.add_argument("--output-root", type=Path, default=reduction.DEFAULT_OUTPUT_ROOT)
    parser.add_argument(
        "--report",
        type=Path,
        default=reduction.DEFAULT_DIAGNOSTICS_DIR / "verification.json",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    started = time.perf_counter()
    inventory = reduction.discover_inventory(args.archive_dir.resolve())
    output_root = args.output_root.resolve()
    reduced_dir = output_root / "reduced"
    calibration_dir = output_root / "calibration"
    errors: list[str] = []

    expected_paths = {
        reduced_dir / reduction.output_name(ref): ref for ref in inventory.science
    }
    actual_paths = set(reduced_dir.glob("*_out.fits"))
    missing = sorted(str(path) for path in expected_paths.keys() - actual_paths)
    unexpected = sorted(str(path) for path in actual_paths - expected_paths.keys())
    partials = sorted(str(path) for path in output_root.rglob("*.partial"))
    if missing:
        errors.append(f"Missing reduced products: {missing}")
    if unexpected:
        errors.append(f"Unexpected reduced products: {unexpected}")
    if partials:
        errors.append(f"Partial products remain: {partials}")

    calibration_products = [
        calibration_dir / "mdark_3.500s.fits",
        calibration_dir / "mdark_50.000s.fits",
        calibration_dir / "mflat_R.fits",
        calibration_dir / "master_flat_review_mask.fits",
    ]
    valid_calibration_checksums = 0
    for path in calibration_products:
        if not path.exists():
            errors.append(f"Missing calibration product: {path}")
            continue
        checksum, datasum = verify_fits_checksum(path)
        if checksum == 1 and datasum == 1:
            valid_calibration_checksums += 1
        else:
            errors.append(
                f"Invalid calibration checksum for {path.name}: "
                f"CHECKSUM={checksum}, DATASUM={datasum}"
            )

    valid_science_checksums = 0
    finite_sampled_science = 0
    bjds: list[float] = []
    for index, (path, ref) in enumerate(expected_paths.items(), start=1):
        if not path.exists():
            continue
        with fits.open(path, memmap=True, checksum=False) as hdul:
            if len(hdul) != 1 or hdul[0].data is None:
                errors.append(f"Invalid HDU layout: {path.name}")
                continue
            hdu = hdul[0]
            header = hdu.header
            checksum = int(hdu.verify_checksum())
            datasum = int(hdu.verify_datasum())
            if checksum == 1 and datasum == 1:
                valid_science_checksums += 1
            else:
                errors.append(
                    f"Invalid checksum for {path.name}: "
                    f"CHECKSUM={checksum}, DATASUM={datasum}"
                )
            if int(header.get("BITPIX", 0)) != -32:
                errors.append(f"Unexpected BITPIX for {path.name}")
            if hdu.data.shape != reduction.EXPECTED_SHAPE:
                errors.append(f"Unexpected image shape for {path.name}: {hdu.data.shape}")
            if not bool(header.get("REDUCED", False)):
                errors.append(f"Missing REDUCED marker for {path.name}")
            if not bool(header.get("DARKCOR", False)) or not bool(
                header.get("FLATCOR", False)
            ):
                errors.append(f"Missing calibration marker for {path.name}")
            if not math.isclose(float(header.get("EXPTIME", math.nan)), 50.0):
                errors.append(f"Unexpected exposure for {path.name}")
            if str(header.get("FILTER", "")).strip().lower() != "red":
                errors.append(f"Unexpected filter for {path.name}")
            if str(header.get("OBJECT", "")).strip() != "TOI_3505.01":
                errors.append(f"Unexpected target for {path.name}")
            history_value = header.get("HISTORY", [])
            if isinstance(history_value, str):
                history = history_value
            else:
                history = "\n".join(str(item) for item in history_value)
            if ref.member not in history:
                errors.append(f"Source provenance mismatch for {path.name}")
            bjd = float(header.get("BJD_TDB", math.nan))
            bjds.append(bjd)
            sample = np.asarray(hdu.data[::64, ::64])
            if np.isfinite(sample).all():
                finite_sampled_science += 1
            else:
                errors.append(f"Non-finite sampled pixels in {path.name}")
        if index % 25 == 0 or index == len(expected_paths):
            print(f"Verified {index}/{len(expected_paths)} science FITS files", flush=True)

    if not np.all(np.diff(np.asarray(bjds, dtype=float)) > 0):
        errors.append("BJD_TDB values are not strictly increasing")

    master_dark = fits.getdata(calibration_dir / "mdark_50.000s.fits").astype(
        np.float32, copy=False
    )
    master_flat = fits.getdata(calibration_dir / "mflat_R.fits").astype(
        np.float32, copy=False
    )
    formula_checks: list[dict[str, float | int]] = []
    rng = np.random.default_rng(3505)
    yy = rng.integers(0, reduction.EXPECTED_SHAPE[0], 256)
    xx = rng.integers(0, reduction.EXPECTED_SHAPE[1], 256)
    for source_index in (0, len(inventory.science) // 2, len(inventory.science) - 1):
        ref = inventory.science[source_index]
        raw, raw_header = reduction.read_member(ref)
        output_path = reduced_dir / reduction.output_name(ref)
        with fits.open(output_path, memmap=True) as hdul:
            output_values = np.asarray(hdul[0].data[yy, xx], dtype=np.float32)
            output_bjd = float(hdul[0].header["BJD_TDB"])
        expected_values = (
            (raw[yy, xx] - master_dark[yy, xx]) / master_flat[yy, xx]
        ).astype(np.float32)
        maximum_difference = float(np.max(np.abs(output_values - expected_values)))
        bjd_difference = float(abs(output_bjd - float(raw_header["BJD_TDB"])))
        formula_checks.append(
            {
                "sequence": int(ref.sequence),
                "sampled_pixels": int(len(yy)),
                "maximum_absolute_difference_adu": maximum_difference,
                "bjd_tdb_absolute_difference_days": bjd_difference,
            }
        )
        if maximum_difference != 0.0:
            errors.append(
                f"Calibration formula mismatch for sequence {ref.sequence}: "
                f"max difference {maximum_difference}"
            )
        if bjd_difference != 0.0:
            errors.append(f"BJD_TDB changed for sequence {ref.sequence}")

    report = {
        "status": "pass" if not errors else "fail",
        "verified_utc": datetime.now(UTC).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "pipeline_version": reduction.PIPELINE_VERSION,
        "science_files_expected": len(expected_paths),
        "science_files_found": len(actual_paths),
        "science_checksums_valid": valid_science_checksums,
        "calibration_checksums_valid": valid_calibration_checksums,
        "sampled_science_files_all_finite": finite_sampled_science,
        "partial_files_found": partials,
        "bjd_tdb_start": float(min(bjds)) if bjds else math.nan,
        "bjd_tdb_end": float(max(bjds)) if bjds else math.nan,
        "bjd_tdb_strictly_increasing": bool(
            len(bjds) == len(expected_paths)
            and np.all(np.diff(np.asarray(bjds, dtype=float)) > 0)
        ),
        "formula_spot_checks": formula_checks,
        "errors": errors,
        "elapsed_seconds": float(time.perf_counter() - started),
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))
    if errors:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
