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
            build_paper_pdf=False,
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
            labels.index("Quantify ephemeris robustness"),
        )
        self.assertLess(
            labels.index("Quantify ephemeris robustness"),
            labels.index("Reconstruct historical schedule"),
        )
        self.assertLess(
            labels.index("Reconstruct historical schedule"),
            labels.index("Plan upcoming observations"),
        )
        self.assertLess(
            labels.index("Plan upcoming observations"),
            labels.index("Build research paper"),
        )
        self.assertLess(
            labels.index("Build research paper"),
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
        consistency_call = next(
            arguments
            for label, arguments in calls
            if label == "Check public-product consistency"
        )
        self.assertEqual(consistency_call[-1], "--verify-manifest")

        manifest_call = next(
            arguments
            for label, arguments in calls
            if label == "Refresh research record"
        )
        self.assertIn("--skip-large-derived", manifest_call)

        paper_call = next(
            arguments for label, arguments in calls if label == "Build research paper"
        )
        self.assertNotIn("--pdf-output", paper_call)

    def test_paper_pdf_flag_is_forwarded(self) -> None:
        args = SimpleNamespace(
            download=False,
            remeasure_ground_apertures=False,
            skip_nearby_images=True,
            skip_large_manifest=True,
            skip_tests=True,
            plate_solve_representative=False,
            build_paper_pdf=True,
        )
        calls: list[tuple[str, list[str]]] = []

        with (
            patch.object(runner, "parse_args", return_value=args),
            patch.object(runner, "require_local_inputs"),
            patch.object(
                runner,
                "run_stage",
                side_effect=lambda label, arguments: calls.append((label, arguments)),
            ),
        ):
            runner.main()

        paper_call = next(
            arguments for label, arguments in calls if label == "Build research paper"
        )
        self.assertEqual(paper_call[-2:], ["--pdf-output", "default"])


if __name__ == "__main__":
    unittest.main()
