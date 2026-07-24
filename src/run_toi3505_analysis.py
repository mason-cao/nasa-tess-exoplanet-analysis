"""Run the current TOI-3505.01 analysis in its documented order.

The default run starts from the saved 25-pixel ground table and downloaded
public TESS/catalog files.  Optional flags can refresh network data or repeat
the slower image-level ground aperture extraction.  Each stage is a normal
standalone script, so a failed step stops the run with its original error.
"""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--download",
        action="store_true",
        help="Refresh MAST, Gaia, and TIC files before analysis.",
    )
    parser.add_argument(
        "--remeasure-ground-apertures",
        action="store_true",
        help="Repeat the 281-frame multi-radius ground extraction.",
    )
    parser.add_argument(
        "--skip-nearby-images",
        action="store_true",
        help="Skip the 281-frame nearby-star measurement.",
    )
    parser.add_argument(
        "--skip-large-manifest",
        action="store_true",
        help="Do not rehash the reduced and aligned FITS collections.",
    )
    parser.add_argument("--skip-tests", action="store_true")
    parser.add_argument(
        "--plate-solve-representative",
        action="store_true",
        help="Run the five planned Astrometry.net solutions; requires the API key.",
    )
    return parser.parse_args()


def run_stage(label: str, arguments: list[str]) -> None:
    print(f"\n[{label}]", flush=True)
    subprocess.run(arguments, cwd=ROOT, check=True)


def script(name: str, *arguments: str) -> list[str]:
    return [sys.executable, str(ROOT / "src" / name), *arguments]


def require_local_inputs() -> None:
    required = (
        ROOT
        / "outputs"
        / "toi3505_aperture_check"
        / "TOI_3505.01_2022-07-22_R_measurements_25px_70-139_AIJ.tbl",
        ROOT / "data" / "catalogs" / "toi3505" / "tic_v8_2p5arcmin.csv",
        ROOT / "data" / "catalogs" / "toi3505" / "gaia_dr3_2p5arcmin.csv",
        ROOT
        / "data"
        / "program_records"
        / "toi3505"
        / "observing_schedule_2022-07-21.json",
        ROOT / "data" / "tess" / "toi3505" / "tesscut",
        ROOT / "data" / "tess" / "toi3505" / "data_validation",
        ROOT
        / "data"
        / "tess"
        / "toi3505"
        / "data_validation_multi_sector",
        ROOT
        / "data"
        / "tess"
        / "toi3505"
        / "exomast"
        / "s0014-s0086_tce1_table.json",
        ROOT
        / "data"
        / "tess"
        / "toi3505"
        / "release_notes"
        / "tess_multisector_14_86_drn122_targetinfo_v01.txt",
    )
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        raise FileNotFoundError(
            "Required local inputs are missing. Run with --download or restore: "
            + "; ".join(missing)
        )


def main() -> None:
    args = parse_args()
    if args.download:
        run_stage("Download TESS products", script("download_toi3505_tess.py"))
        run_stage("Download target catalogs", script("download_toi3505_catalogs.py"))
        run_stage(
            "Download comparison-star Gaia neighborhoods",
            script("download_toi3505_comparison_gaia.py"),
        )
    require_local_inputs()

    if args.remeasure_ground_apertures:
        run_stage(
            "Repeat ground aperture and timing checks",
            script("check_toi3505_aperture_radii.py"),
        )

    run_stage("Build final ground light curve", script("make_toi3505_final_candidate.py"))
    ground_arguments = script("make_toi3505_ground_checks.py")
    if args.skip_nearby_images:
        ground_arguments.append("--skip-nearby-images")
    run_stage("Run remaining ground checks", ground_arguments)
    run_stage("Measure all four TESS sectors", script("analyze_toi3505_tess.py"))
    run_stage("Run TESS pixel checks", script("analyze_toi3505_tess_pixels.py"))
    run_stage(
        "Compare official SPOC reports",
        script("analyze_toi3505_data_validation.py"),
    )

    plate_arguments = script("plate_solve_toi3505_representative.py")
    if args.plate_solve_representative:
        if not os.environ.get("ASTROMETRY_NET_API_KEY", "").strip():
            raise RuntimeError(
                "ASTROMETRY_NET_API_KEY is required for --plate-solve-representative"
            )
        plate_arguments.append("--run")
    run_stage("Update representative plate-solve plan", plate_arguments)

    if not args.skip_tests:
        run_stage(
            "Run tests",
            [sys.executable, "-m", "unittest", "discover", "-s", "tests", "-v"],
        )

    manifest_arguments = script("make_toi3505_research_record.py")
    if args.skip_large_manifest:
        manifest_arguments.append("--skip-large-derived")
    run_stage("Refresh research record", manifest_arguments)
    print("\nTOI-3505.01 analysis completed successfully.", flush=True)


if __name__ == "__main__":
    main()
