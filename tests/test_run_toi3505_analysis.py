from __future__ import annotations

import sys
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import run_toi3505_analysis as runner  # noqa: E402


class AnalysisRunnerTests(unittest.TestCase):
    def test_refinement_planning_and_consistency_run_in_order(self) -> None:
        args = SimpleNamespace(
            download=False,
            remeasure_ground_apertures=False,
            skip_nearby_images=True,
            skip_large_manifest=True,
            skip_tests=False,
            plate_solve_representative=False,
        )
        calls: list[tuple[str, list[str]]] = []

        def record(label: str, arguments: list[str]) -> None:
            calls.append((label, arguments))

        with (
            patch.object(runner, "parse_args", return_value=args),
            patch.object(runner, "require_local_inputs"),
            patch.object(runner, "run_stage", side_effect=record),
        ):
            runner.main()

        labels = [label for label, _ in calls]
        self.assertLess(
            labels.index("Measure all four TESS sectors"),
            labels.index("Refine TESS ephemeris"),
        )
        self.assertLess(
            labels.index("Compare official SPOC reports"),
            labels.index("Refine TESS ephemeris"),
        )
        self.assertLess(
            labels.index("Refine TESS ephemeris"),
            labels.index("Plan upcoming observations"),
        )
        self.assertLess(
            labels.index("Plan upcoming observations"),
            labels.index("Refresh research record"),
        )
        self.assertLess(
            labels.index("Refresh research record"),
            labels.index("Check public-product consistency"),
        )
        self.assertLess(
            labels.index("Check public-product consistency"),
            labels.index("Run tests"),
        )

        manifest_call = next(
            arguments
            for label, arguments in calls
            if label == "Refresh research record"
        )
        self.assertIn("--skip-large-derived", manifest_call)


if __name__ == "__main__":
    unittest.main()
