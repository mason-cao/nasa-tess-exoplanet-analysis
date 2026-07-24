"""Plan or run plate solutions for representative TOI-3505.01 frames.

The existing plate solution uses the first calibrated exposure.  This wrapper
selects early, middle, late, and meridian-bracketing frames from the measured
night and can pass each one to ``plate_solve_toi3505.py``.  By default it only
writes the frozen frame plan.  Use ``--detect-only`` for local source-list
checks or ``--run`` when an Astrometry.net API key is available.
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
from datetime import UTC, datetime
from pathlib import Path

import numpy as np
import pandas as pd

from analyze_toi3505_photometry import load_table


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_TABLE = (
    ROOT
    / "outputs"
    / "toi3505_aperture_check"
    / "TOI_3505.01_2022-07-22_R_measurements_25px_70-139_AIJ.tbl"
)
DEFAULT_REDUCED_DIR = ROOT / "data" / "ground" / "toi3505" / "reduced"
DEFAULT_SOLVED_DIR = (
    ROOT / "data" / "ground" / "toi3505" / "plate_solved" / "representative"
)
DEFAULT_REPORT_DIR = ROOT / "outputs" / "toi3505_representative_plate_solve"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--table", type=Path, default=DEFAULT_TABLE)
    parser.add_argument("--reduced-dir", type=Path, default=DEFAULT_REDUCED_DIR)
    parser.add_argument("--solved-dir", type=Path, default=DEFAULT_SOLVED_DIR)
    parser.add_argument("--report-dir", type=Path, default=DEFAULT_REPORT_DIR)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument(
        "--detect-only",
        action="store_true",
        help="Measure local source lists for every planned frame.",
    )
    mode.add_argument(
        "--run",
        action="store_true",
        help="Submit source lists to Astrometry.net and save the solutions.",
    )
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def representative_indices(table: pd.DataFrame, bracket_frames: int = 12) -> dict[str, int]:
    """Return stable row indices for the five representative roles."""
    if len(table) < 5:
        raise ValueError("At least five frames are required")
    meridian = int(np.nanargmin(table["AIRMASS"].to_numpy(dtype=float)))
    choices = {
        "early": 0,
        "middle": len(table) // 2,
        "before_meridian": max(0, meridian - bracket_frames),
        "after_meridian": min(len(table) - 1, meridian + bracket_frames),
        "late": len(table) - 1,
    }
    # With a short sequence, roles could collide. Move a duplicate to the
    # nearest unused row so every requested check remains a separate image.
    used: set[int] = set()
    for role, index in list(choices.items()):
        if index not in used:
            used.add(index)
            continue
        alternatives = sorted(range(len(table)), key=lambda value: abs(value - index))
        replacement = next(value for value in alternatives if value not in used)
        choices[role] = replacement
        used.add(replacement)
    return choices


def reduced_name(aligned_name: str) -> str:
    suffix = "_aligned.fits"
    if not aligned_name.endswith(suffix):
        raise ValueError(f"Unexpected aligned-image name: {aligned_name}")
    return aligned_name[: -len(suffix)] + "_out.fits"


def build_plan(table: pd.DataFrame, reduced_dir: Path, solved_dir: Path) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for role, index in representative_indices(table).items():
        row = table.iloc[index]
        input_path = reduced_dir / reduced_name(str(row["Label"]))
        output_name = input_path.stem.replace("_out", "_wcs") + ".fits"
        rows.append(
            {
                "role": role,
                "row_index": index,
                "slice": int(row["slice"]),
                "bjd_tdb": float(row["BJD_TDB"]),
                "airmass": float(row["AIRMASS"]),
                "input_path": str(input_path.resolve()),
                "output_path": str((solved_dir / output_name).resolve()),
                "input_exists": input_path.exists(),
            }
        )
    return pd.DataFrame(rows).sort_values("slice").reset_index(drop=True)


def run_plan(plan: pd.DataFrame, report_dir: Path, *, detect_only: bool, overwrite: bool) -> None:
    script = ROOT / "src" / "plate_solve_toi3505.py"
    for row in plan.itertuples(index=False):
        frame_report = report_dir / f"slice_{row.slice:04d}_{row.role}"
        command = [
            sys.executable,
            str(script),
            "--input",
            row.input_path,
            "--output",
            row.output_path,
            "--report-dir",
            str(frame_report),
        ]
        if detect_only:
            command.append("--detect-only")
        if overwrite:
            command.append("--overwrite")
        print(f"Processing slice {row.slice} ({row.role})", flush=True)
        subprocess.run(command, check=True)


def main() -> None:
    args = parse_args()
    report_dir = args.report_dir.resolve()
    solved_dir = args.solved_dir.resolve()
    report_dir.mkdir(parents=True, exist_ok=True)
    solved_dir.mkdir(parents=True, exist_ok=True)
    table = load_table(
        args.table.resolve(), expected_outer_radius=139, expected_source_radius=25
    )
    plan = build_plan(table, args.reduced_dir.resolve(), solved_dir)
    if not bool(plan["input_exists"].all()):
        missing = plan.loc[~plan["input_exists"], "input_path"].tolist()
        raise FileNotFoundError(f"Missing representative inputs: {missing}")
    plan.to_csv(report_dir / "frame_plan.csv", index=False, float_format="%.10f")

    mode = "plan_only"
    if args.detect_only:
        mode = "source_detection_complete"
        run_plan(plan, report_dir, detect_only=True, overwrite=args.overwrite)
    elif args.run:
        if not os.environ.get("ASTROMETRY_NET_API_KEY", "").strip():
            raise RuntimeError(
                "ASTROMETRY_NET_API_KEY is required for --run; the key is never saved"
            )
        mode = "plate_solutions_requested"
        run_plan(plan, report_dir, detect_only=False, overwrite=args.overwrite)

    summary = {
        "target": "TOI-3505.01",
        "created_utc": datetime.now(UTC).isoformat(),
        "mode": mode,
        "representative_frames": plan.to_dict(orient="records"),
        "api_key_saved": False,
        "independent_plate_solutions_complete": bool(args.run),
    }
    (report_dir / "summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    readme = f"""# Representative TOI-3505.01 plate-solve plan

Five frames are frozen in `frame_plan.csv`: early, middle, late, and one on
each side of the minimum-airmass point. Current mode: `{mode}`.

The first frame already has a valid plate solution. The additional independent
solutions require `ASTROMETRY_NET_API_KEY` and can be run with:

```bash
.venv/bin/python src/plate_solve_toi3505_representative.py --run
```

Only measured source positions are submitted. The full science images and API
key are not uploaded or written by the scripts.
"""
    (report_dir / "README.md").write_text(readme, encoding="utf-8")
    print(f"Saved representative frame plan to {report_dir}")


if __name__ == "__main__":
    main()
