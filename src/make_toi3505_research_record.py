"""Build the reproducibility and claim-evidence record for TOI-3505.01.

The manifest covers the six original archives, calibration files, every
reduced and aligned ground image, downloaded TESS products, catalog snapshots,
load-bearing tables, code, and compact scientific outputs.  SHA-256 hashes are
calculated in streaming chunks so even the large FITS collection is handled
without loading whole files into memory.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import os
import platform
import plistlib
import subprocess
import sys
from datetime import UTC, datetime
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT_DIR = ROOT / "outputs" / "toi3505_research_record"
ORIGINAL_ARCHIVE_PATTERN = "TOI_3505.01-20260714T190730Z-1-*.zip"
PACKAGE_NAMES = (
    "numpy",
    "pandas",
    "astropy",
    "astroquery",
    "matplotlib",
    "scipy",
    "photutils",
    "lightkurve",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument(
        "--skip-large-derived",
        action="store_true",
        help="Skip reduced and aligned FITS hashes for a faster partial refresh.",
    )
    return parser.parse_args()


def sha256_file(path: Path, block_size: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(block_size):
            digest.update(block)
    return digest.hexdigest()


def add_group(
    groups: list[tuple[str, Path]], category: str, paths: list[Path]
) -> None:
    for path in sorted(set(paths)):
        if path.is_file():
            groups.append((category, path))


def manifest_paths(include_large_derived: bool) -> list[tuple[str, Path]]:
    groups: list[tuple[str, Path]] = []
    archives = list((ROOT / "data_and_lectures").glob(ORIGINAL_ARCHIVE_PATTERN))
    if len(archives) != 6:
        raise RuntimeError(f"Expected six original TOI-3505 archives; found {len(archives)}")
    add_group(groups, "original_archive", archives)
    add_group(
        groups,
        "ground_calibration",
        list((ROOT / "data" / "ground" / "toi3505" / "calibration").glob("*.fits")),
    )
    if include_large_derived:
        add_group(
            groups,
            "ground_reduced_image",
            list((ROOT / "data" / "ground" / "toi3505" / "reduced").glob("*.fits")),
        )
        add_group(
            groups,
            "ground_aligned_image",
            list((ROOT / "data" / "ground" / "toi3505" / "aligned").glob("*.fits")),
        )
    add_group(
        groups,
        "plate_solution",
        list((ROOT / "data" / "ground" / "toi3505" / "plate_solved").glob("*.fits")),
    )
    add_group(
        groups,
        "tess_download",
        [
            path
            for path in (ROOT / "data" / "tess" / "toi3505").rglob("*")
            if path.is_file()
            and path.suffix.lower()
            in {".fits", ".pdf", ".xml", ".json", ".txt"}
        ],
    )
    catalog_dir = ROOT / "data" / "catalogs" / "toi3505"
    add_group(
        groups,
        "catalog_snapshot",
        [path for path in catalog_dir.iterdir() if path.is_file()],
    )
    program_record_dir = ROOT / "data" / "program_records" / "toi3505"
    add_group(
        groups,
        "program_record",
        [path for path in program_record_dir.iterdir() if path.is_file()],
    )

    load_bearing = [
        ROOT
        / "outputs"
        / "toi3505_aperture_check"
        / "TOI_3505.01_2022-07-22_R_measurements_25px_70-139_AIJ.tbl",
        ROOT
        / "outputs"
        / "toi3505_final_candidate"
        / "TOI_3505.01_2022-07-22_R_final_candidate_light_curve.csv",
        ROOT / "outputs" / "toi3505_final_candidate" / "frozen_protocol.json",
        ROOT / "outputs" / "toi3505_final_candidate" / "analysis_settings.json",
        ROOT / "outputs" / "toi3505_reduction" / "frame_manifest.csv",
        ROOT / "outputs" / "toi3505_alignment" / "frame_shifts.csv",
        ROOT / "outputs" / "toi3505_aperture_check" / "aperture_radius_metrics.csv",
    ]
    add_group(groups, "load_bearing_table", load_bearing)

    evidence_folders = (
        "toi3505_final_candidate",
        "toi3505_ground_checks",
        "toi3505_tess_download",
        "toi3505_tess_analysis",
        "toi3505_tess_pixels",
        "toi3505_data_validation",
        "toi3505_representative_plate_solve",
    )
    evidence: list[Path] = []
    for folder in evidence_folders:
        directory = ROOT / "outputs" / folder
        if not directory.exists():
            continue
        evidence.extend(
            path
            for path in directory.iterdir()
            if path.is_file() and path.suffix.lower() in {".csv", ".json", ".md"}
        )
    add_group(groups, "scientific_output", evidence)
    add_group(groups, "analysis_code", list((ROOT / "src").glob("*.py")))
    add_group(groups, "analysis_test", list((ROOT / "tests").glob("test_*.py")))

    # Deduplicate paths while preserving the first, most fundamental category.
    unique: list[tuple[str, Path]] = []
    seen: set[Path] = set()
    for category, path in groups:
        resolved = path.resolve()
        if resolved in seen:
            continue
        unique.append((category, resolved))
        seen.add(resolved)
    return unique


def build_manifest(paths: list[tuple[str, Path]]) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    total = len(paths)
    for index, (category, path) in enumerate(paths, start=1):
        stat = path.stat()
        rows.append(
            {
                "category": category,
                "path": str(path.relative_to(ROOT)),
                "size_bytes": stat.st_size,
                "modified_utc": datetime.fromtimestamp(stat.st_mtime, UTC).isoformat(),
                "sha256": sha256_file(path),
            }
        )
        if index % 25 == 0 or index == total:
            print(f"Hashed {index}/{total} files", flush=True)
    return pd.DataFrame(rows)


def command_version(command: list[str]) -> str | None:
    try:
        result = subprocess.run(
            command, check=False, capture_output=True, text=True, timeout=20
        )
    except (FileNotFoundError, subprocess.TimeoutExpired):
        return None
    text = (result.stdout + "\n" + result.stderr).strip()
    return text.splitlines()[0] if text else None


def astroimagej_version() -> str | None:
    info = Path("/Applications/AstroImageJ.app/Contents/Info.plist")
    if not info.exists():
        return None
    try:
        with info.open("rb") as handle:
            data = plistlib.load(handle)
    except (OSError, plistlib.InvalidFileException):
        return None
    return str(
        data.get("CFBundleShortVersionString")
        or data.get("CFBundleVersion")
        or "installed; bundle version not listed"
    )


def astroimagej_java_version() -> str | None:
    release = Path(
        "/Applications/AstroImageJ.app/Contents/runtime/Contents/Home/release"
    )
    if not release.exists():
        return None
    for line in release.read_text(encoding="utf-8").splitlines():
        if line.startswith("JAVA_VERSION="):
            return line.split("=", 1)[1].strip().strip('"')
    return None


def git_state() -> dict[str, object]:
    commit = command_version(["git", "rev-parse", "HEAD"])
    status = subprocess.run(
        ["git", "status", "--porcelain"],
        cwd=ROOT,
        check=False,
        capture_output=True,
        text=True,
    ).stdout
    return {"commit": commit, "worktree_clean": not bool(status.strip())}


def software_record() -> dict[str, object]:
    packages = {}
    for name in PACKAGE_NAMES:
        try:
            packages[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            packages[name] = None
    return {
        "created_utc": datetime.now(UTC).isoformat(),
        "python": sys.version,
        "python_executable": sys.executable,
        "operating_system": platform.platform(),
        "machine": platform.machine(),
        "packages": packages,
        "system_java": command_version(["java", "-version"]),
        "astroimagej_bundled_java": astroimagej_java_version(),
        "astroimagej": astroimagej_version(),
        "git": git_state(),
    }


def frozen_config() -> dict[str, object]:
    protocol_path = ROOT / "outputs" / "toi3505_final_candidate" / "frozen_protocol.json"
    ground_protocol = json.loads(protocol_path.read_text(encoding="utf-8"))
    schedule_path = (
        ROOT
        / "outputs"
        / "toi3505_final_candidate"
        / "historical_schedule_check.json"
    )
    schedule = json.loads(schedule_path.read_text(encoding="utf-8"))
    return {
        "target": {"name": "TOI-3505.01", "tic_id": 390988385},
        "ground": ground_protocol,
        "historical_schedule": {
            "source_record": schedule["source_record"],
            "working_interpretation": schedule["working_interpretation"],
            "utc_alternative": schedule["utc_alternative"],
            "fixed_window_check": schedule["fixed_window_check"],
            "historical_ephemeris_complete": False,
        },
        "tess": {
            "sectors": [14, 41, 54, 81],
            "primary_pipeline": "QLP",
            "quality_rule": "finite flux and uncertainty with QUALITY equal to zero",
            "period_days": 2.9151556,
            "epoch_bjd_tdb": 2459793.534385,
            "depth_model": "exposure-integrated box with a local straight baseline per event",
            "physical_transit_model": False,
            "reference_tesscut_aperture_radius_pixels": 3.0,
            "tesscut_background_annulus_pixels": [4.0, 7.0],
            "official_spoc_report_sectors": [54, 81],
            "official_spoc_combined_report": {
                "search_scope": "s0014-s0086",
                "contributing_sectors": [54, 81],
                "cadence_minutes": 10,
                "data_release": 122,
            },
            "official_spoc_report_role": "same-observation method comparison",
        },
        "claim_limits": [
            "The GMU night is off transit under the current ephemeris, but the recovered 2022 schedule placed an event inside it.",
            "The schedule row has no explicit time zone, epoch, uncertainty, depth, or prediction source.",
            "The GMU window falls in a Sector 54 data gap.",
            "TESS-scale localization cannot separate the 0.517-arcsecond companion.",
            "The dilution screen is not a final correction.",
            "The box model is not a physical planet fit.",
        ],
    }


def claim_evidence_rows() -> list[dict[str, str]]:
    return [
        {
            "claim": "The measured ground target is TOI-3505.01.",
            "evidence": "outputs/toi3505_plate_solve/solution.json; outputs/toi3505_ground_checks/comparison_star_catalog_matches.csv",
            "status": "supported",
            "limit": "Only the first frame has an independent plate solution so far.",
        },
        {
            "claim": "The adopted ground curve uses a 25-pixel aperture and all ten comparison stars.",
            "evidence": "outputs/toi3505_aperture_check/aperture_radius_metrics.csv; outputs/toi3505_final_candidate/frozen_protocol.json",
            "status": "supported",
            "limit": "Mentor review remains requested.",
        },
        {
            "claim": "No detrending correction was justified for the posted ground curve.",
            "evidence": "outputs/toi3505_final_candidate/model_selection.csv",
            "status": "supported",
            "limit": "Applies to the tested variables and blocked checks only.",
        },
        {
            "claim": "The 2022 program schedule placed ingress and egress inside the GMU sequence under the documented Eastern-time interpretation.",
            "evidence": "data/program_records/toi3505/observing_schedule_2022-07-21.json; outputs/toi3505_final_candidate/historical_schedule_check.json",
            "status": "supported with a stated time-zone assumption",
            "limit": "The pasted row has no explicit time zone, epoch, uncertainty, depth, or original workbook/URL.",
        },
        {
            "claim": "The fixed historical-window check measures -0.658 +/- 0.395 ppt, so it finds no transit-like dimming; a 2.91-ppt same-window injection is recovered with a 2.910-ppt increment.",
            "evidence": "outputs/toi3505_final_candidate/historical_schedule_check.json; outputs/toi3505_final_candidate/historical_schedule_times.csv",
            "status": "supported under the Eastern-time schedule interpretation",
            "limit": "This is a straight-baseline exposure-integrated box check, not a physical transit fit, and the schedule timezone is unconfirmed.",
        },
        {
            "claim": "The GMU sequence does not cover a predicted transit under the current catalog ephemeris.",
            "evidence": "outputs/toi3505_final_candidate/summary.json; outputs/toi3505_tess_analysis/ephemeris_checks.json",
            "status": "supported",
            "limit": "This does not erase the conflicting historical schedule window; the schedule's full ephemeris is still unavailable.",
        },
        {
            "claim": "The TESS signal is recovered in Sectors 14, 41, 54, and 81.",
            "evidence": "outputs/toi3505_tess_analysis/sector_measurements.csv; outputs/toi3505_tess_analysis/period_search.csv",
            "status": "supported",
            "limit": "All products are measurements of the same star system, not four independent instruments.",
        },
        {
            "claim": "The two per-sector and one combined SPOC reports recover the same signal without a strong odd/even, secondary-event, or centroid warning in Sectors 54 and 81.",
            "evidence": "outputs/toi3505_data_validation/official_dv_metrics.csv; outputs/toi3505_data_validation/official_multisector_tce.csv; official DVR and DVS files listed in those tables",
            "status": "supported",
            "limit": "The s0014-s0086 label is a search range, but only Sectors 54 and 81 contributed for this target. SPOC and project measurements reuse the same TESS observations, and neither resolves the 0.517-arcsecond companion.",
        },
        {
            "claim": "No simultaneous TESS/GMU flux comparison is possible.",
            "evidence": "outputs/toi3505_tess_analysis/sector54_ground_overlap.json",
            "status": "supported",
            "limit": "The GMU sequence falls inside a Sector 54 data gap.",
        },
        {
            "claim": "The TESS difference signal is consistent with the target system at one-pixel resolution.",
            "evidence": "outputs/toi3505_tess_pixels/difference_image_localization.csv",
            "status": "supported",
            "limit": "One TESS pixel cannot resolve the 0.517-arcsecond companion.",
        },
        {
            "claim": "A final dilution-corrected planet radius has not been measured.",
            "evidence": "outputs/toi3505_tess_pixels/dilution_screen.json; outputs/toi3505_tess_analysis/pipeline_comparison.csv",
            "status": "open",
            "limit": "QLP crowding treatment and band-dependent companion contrast remain unresolved.",
        },
        {
            "claim": "The preliminary ground screen conditionally clears 32 of 44 bright-enough catalog candidates in the historical window; two target-aperture overlaps and ten incomplete or noisy cases remain.",
            "evidence": "outputs/toi3505_ground_checks/nearby_star_image_measurements.csv",
            "status": "conditional screening result",
            "limit": "It assumes Eastern local time, is not the formal program NEB procedure, and cannot resolve the 0.517-arcsecond companion.",
        },
    ]


def decision_rows() -> list[dict[str, str]]:
    return [
        {
            "date": "2026-07-23",
            "decision": "Keep the ground light curve undetrended.",
            "reason": "No tested correction passed the predeclared blocked out-of-sample rule.",
            "alternative": "Airmass, sky, width, position, comparison counts, and combined models were retained as checks.",
        },
        {
            "date": "2026-07-23",
            "decision": "Use all ten comparison stars in the primary flux sum.",
            "reason": "The ensemble was not selected by minimizing target scatter; pseudo-target checks were acceptable.",
            "alternative": "Equal-star and inverse-error ensembles were saved as robustness comparisons.",
        },
        {
            "date": "2026-07-23",
            "decision": "Do not fit a transit to the GMU sequence.",
            "reason": "The current and historical timing sources disagree, and a fixed historical-window box gives no significant dimming.",
            "alternative": "Display the historical window and report its fixed-box check without a physical transit fit.",
        },
        {
            "date": "2026-07-23",
            "decision": "Use a finite-exposure box for first-pass TESS measurements.",
            "reason": "It is transparent across cadences and does not imply unsupported physical precision.",
            "alternative": "A physical limb-darkened fit remains gated by dilution and companion information.",
        },
        {
            "date": "2026-07-23",
            "decision": "Treat nearby-star and dilution calculations as screens.",
            "reason": "The ground clearance is conditional on an unlabeled schedule timezone, the current run is not the formal program NEB procedure, and neither the ground apertures nor TESS resolve the close companion.",
            "alternative": "Confirm the schedule source, complete the formal program NEB procedure, and use resolving observations or permitted ExoFOP products.",
        },
        {
            "date": "2026-07-23",
            "decision": "Use official SPOC reports as same-observation checks, not independent validation.",
            "reason": "The two per-sector reports, the combined report, and the project measurements all use the same Sector 54 and 81 TESS pixels.",
            "alternative": "Keep the official diagnostics beside the project results while preserving the close-companion limitation.",
        },
        {
            "date": "2026-07-24",
            "decision": "Use America/New_York as the working interpretation of the schedule clocks while preserving a UTC alternative.",
            "reason": "The sheet's planned 21:10-04:55 range brackets the actual image sequence only under the Eastern civil-time interpretation.",
            "alternative": "Obtain the original workbook or Transit Info file and replace the assumption if its metadata says otherwise.",
        },
    ]


def dependence_rows() -> list[dict[str, str]]:
    return [
        {
            "product": "ground final light curve",
            "depends_on": "original archives; calibration masters; reduced and aligned images; 25-pixel AstroImageJ table; 2022 observing-schedule row",
            "shared_data_warning": "Photutils and AstroImageJ extractions use the same ground images; the schedule time zone is a documented assumption.",
        },
        {
            "product": "four-sector TESS measurements",
            "depends_on": "MAST QLP light curves; catalog ephemeris",
            "shared_data_warning": "SPOC and TESS-SPOC checks reuse the same TESS observations.",
        },
        {
            "product": "TESS custom-aperture depths and difference images",
            "depends_on": "MAST TESScut cubes; TESS box timing results",
            "shared_data_warning": "Aperture variants are not independent observations.",
        },
        {
            "product": "official SPOC report comparison",
            "depends_on": "MAST DVR, DVM, and DVS PDFs; DVR XML; DVT FITS; DR122 target table; Exo.MAST combined TCE table; project QLP and SPOC measurements",
            "shared_data_warning": "The official and project results reuse the same Sector 54 and 81 TESS observations.",
        },
        {
            "product": "dilution screen",
            "depends_on": "TIC v8; Delta-I companion measurement stated in project context; catalog depth",
            "shared_data_warning": "Delta I is only a proxy for TESS-band contrast.",
        },
        {
            "product": "ground nearby-star screen",
            "depends_on": "TIC v8; first-frame WCS; all aligned ground images; comparison counts",
            "shared_data_warning": "Clearance is conditional on the EDT schedule interpretation, and the target aperture cannot separate the 0.517-arcsecond companion.",
        },
    ]


def main() -> None:
    args = parse_args()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    include_large = not args.skip_large_derived
    paths = manifest_paths(include_large)
    manifest = build_manifest(paths)
    manifest.to_csv(output_dir / "file_manifest.csv", index=False)

    software = software_record()
    (output_dir / "software_versions.json").write_text(
        json.dumps(software, indent=2) + "\n", encoding="utf-8"
    )
    (output_dir / "frozen_analysis_config.json").write_text(
        json.dumps(frozen_config(), indent=2) + "\n", encoding="utf-8"
    )
    pd.DataFrame(claim_evidence_rows()).to_csv(
        output_dir / "claim_evidence.csv", index=False
    )
    pd.DataFrame(decision_rows()).to_csv(output_dir / "decision_log.csv", index=False)
    pd.DataFrame(dependence_rows()).to_csv(
        output_dir / "measurement_dependence.csv", index=False
    )

    category_counts = manifest.groupby("category").size().to_dict()
    category_bytes = manifest.groupby("category")["size_bytes"].sum().to_dict()
    summary = {
        "created_utc": datetime.now(UTC).isoformat(),
        "manifest_scope": "full" if include_large else "large derived images skipped",
        "files": len(manifest),
        "bytes": int(manifest["size_bytes"].sum()),
        "category_counts": {key: int(value) for key, value in category_counts.items()},
        "category_bytes": {key: int(value) for key, value in category_bytes.items()},
        "original_archive_count": int(
            (manifest["category"] == "original_archive").sum()
        ),
        "sha256_complete": True,
    }
    (output_dir / "manifest_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    readme = f"""# TOI-3505.01 research record

This folder freezes the reproducibility record for the current analysis.

- `{len(manifest)}` files are listed with size, modification time, and SHA-256.
- The six original archives are included.
- Reduced and aligned images are {'included' if include_large else 'not included in this refresh'}.
- Software and system versions are in `software_versions.json`.
- Ground and TESS choices are in `frozen_analysis_config.json`.
- Claims, decisions, and shared-data dependencies have separate CSV ledgers.

Regenerate the full record with:

```bash
.venv/bin/python src/make_toi3505_research_record.py
```

The manifest records files; it does not make unresolved provenance or
scientific limitations disappear. Those limits remain in the claim ledger.
"""
    (output_dir / "README.md").write_text(readme, encoding="utf-8")
    print(f"Saved research record to {output_dir}")


if __name__ == "__main__":
    main()
