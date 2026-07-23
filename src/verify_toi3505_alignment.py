"""Check every TOI-3505.01 aligned image and its recorded pixel shift."""

from __future__ import annotations

import argparse
import json
import math
import time
from datetime import UTC, datetime
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.io import fits

import align_toi3505 as alignment


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, default=alignment.DEFAULT_INPUT_DIR)
    parser.add_argument("--output-dir", type=Path, default=alignment.DEFAULT_OUTPUT_DIR)
    parser.add_argument(
        "--shift-table",
        type=Path,
        default=alignment.DEFAULT_DIAGNOSTICS_DIR / "frame_shifts.csv",
    )
    parser.add_argument(
        "--report",
        type=Path,
        default=alignment.DEFAULT_DIAGNOSTICS_DIR / "verification.json",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    started = time.perf_counter()
    input_dir = args.input_dir.resolve()
    output_dir = args.output_dir.resolve()
    frames = alignment.discover_frames(input_dir)
    shift_table = pd.read_csv(args.shift_table).set_index("sequence")
    errors: list[str] = []

    expected_outputs = {
        output_dir / alignment.aligned_name(frame.path): frame for frame in frames
    }
    actual_outputs = set(output_dir.glob("*_aligned.fits"))
    missing = sorted(path.name for path in expected_outputs.keys() - actual_outputs)
    unexpected = sorted(path.name for path in actual_outputs - expected_outputs.keys())
    partials = sorted(path.name for path in output_dir.glob("*.partial"))
    if missing:
        errors.append(f"Missing aligned images: {missing}")
    if unexpected:
        errors.append(f"Unexpected aligned images: {unexpected}")
    if partials:
        errors.append(f"Unfinished images remain: {partials}")

    valid_checksums = 0
    exact_pixel_checks = 0
    preserved_times = 0
    x_shifts: list[int] = []
    y_shifts: list[int] = []
    bjd_values: list[float] = []
    check_points = np.asarray([512, 1024, 2048, 3072, 3584], dtype=int)
    source_y, source_x = np.meshgrid(check_points, check_points, indexing="ij")

    for index, (output_path, frame) in enumerate(expected_outputs.items(), start=1):
        if not output_path.exists():
            continue
        with fits.open(frame.path, memmap=True, checksum=False) as input_hdul, fits.open(
            output_path, memmap=True, checksum=False
        ) as output_hdul:
            input_hdu = input_hdul[0]
            output_hdu = output_hdul[0]
            header = output_hdu.header
            if int(output_hdu.verify_checksum()) == 1 and int(
                output_hdu.verify_datasum()
            ) == 1:
                valid_checksums += 1
            else:
                errors.append(f"Checksum failed for {output_path.name}")
            if output_hdu.data.shape != alignment.EXPECTED_SHAPE:
                errors.append(f"Image shape changed for {output_path.name}")
            if int(header.get("BITPIX", 0)) != -32:
                errors.append(f"Image is not 32-bit for {output_path.name}")
            if not bool(header.get("ALIGNED", False)):
                errors.append(f"Alignment marker is missing from {output_path.name}")
            if str(header.get("ALNMETH", "")) != "WHOLEPIX":
                errors.append(f"Whole-pixel method is missing from {output_path.name}")

            x_shift = int(header.get("X_SHIFT", 999999))
            y_shift = int(header.get("Y_SHIFT", 999999))
            x_shifts.append(x_shift)
            y_shifts.append(y_shift)
            if frame.sequence not in shift_table.index:
                errors.append(f"Frame {frame.sequence:04d} is missing from frame_shifts.csv")
            else:
                row = shift_table.loc[frame.sequence]
                if x_shift != int(row["x_shift_pixels"]) or y_shift != int(
                    row["y_shift_pixels"]
                ):
                    errors.append(f"Recorded shift does not match frame {frame.sequence:04d}")

            output_y = source_y + y_shift
            output_x = source_x + x_shift
            input_values = np.asarray(input_hdu.data[source_y, source_x])
            output_values = np.asarray(output_hdu.data[output_y, output_x])
            if np.array_equal(input_values, output_values):
                exact_pixel_checks += 1
            else:
                errors.append(f"Pixel values changed in {output_path.name}")

            input_bjd = float(input_hdu.header.get("BJD_TDB", math.nan))
            output_bjd = float(header.get("BJD_TDB", math.nan))
            bjd_values.append(output_bjd)
            if input_bjd == output_bjd:
                preserved_times += 1
            else:
                errors.append(f"Observation time changed in {output_path.name}")

            fill_value = float(header.get("ALNFILL", math.nan))
            if x_shift > 0:
                border_value = float(output_hdu.data[2048, 0])
            elif x_shift < 0:
                border_value = float(output_hdu.data[2048, -1])
            elif y_shift > 0:
                border_value = float(output_hdu.data[0, 2048])
            elif y_shift < 0:
                border_value = float(output_hdu.data[-1, 2048])
            else:
                border_value = fill_value
            if border_value != fill_value:
                errors.append(f"Shifted edge was not filled correctly in {output_path.name}")

        if index % 25 == 0 or index == len(expected_outputs):
            print(f"Checked {index}/{len(expected_outputs)} aligned images", flush=True)

    times_in_order = bool(
        len(bjd_values) == len(expected_outputs)
        and np.all(np.diff(np.asarray(bjd_values)) > 0)
    )
    if not times_in_order:
        errors.append("The observation times are not in increasing order")

    report = {
        "status": "pass" if not errors else "fail",
        "verified_utc": datetime.now(UTC).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "alignment_version": alignment.ALIGNMENT_VERSION,
        "aligned_images_expected": len(expected_outputs),
        "aligned_images_found": len(actual_outputs),
        "valid_fits_checksums": valid_checksums,
        "exact_whole_pixel_checks": exact_pixel_checks,
        "observation_times_preserved": preserved_times,
        "observation_times_in_order": times_in_order,
        "x_shift_range_pixels": [min(x_shifts), max(x_shifts)] if x_shifts else [],
        "y_shift_range_pixels": [min(y_shifts), max(y_shifts)] if y_shifts else [],
        "unfinished_files": partials,
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
