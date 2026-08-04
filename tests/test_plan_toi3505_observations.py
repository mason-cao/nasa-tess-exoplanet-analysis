from __future__ import annotations

import json
import sys
import tempfile
import unittest
from datetime import date
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import plan_toi3505_observations as planner  # noqa: E402


class ObservationPlannerTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.ephemeris_path = (
            ROOT / "outputs" / "toi3505_ephemeris_refined" / "ephemeris_refined.json"
        )
        cls.validation_path = (
            ROOT / "outputs" / "toi3505_data_validation" / "analysis_summary.json"
        )
        cls.ephemeris = planner.load_adopted_ephemeris(cls.ephemeris_path)
        cls.duration = planner.load_spoc_duration_hours(cls.validation_path)
        cls.events = planner.plan_events(
            cls.ephemeris,
            date(2026, 8, 4),
            date(2026, 10, 1),
            cls.duration,
            planner.PlanningLimits(),
        )

    def test_covariance_is_used_in_prediction_uncertainty(self) -> None:
        cycle = 505
        actual = planner.prediction_uncertainty_minutes(self.ephemeris, cycle)
        without_covariance = dict(self.ephemeris)
        without_covariance["covariance_days2"] = 0.0
        uncorrelated = planner.prediction_uncertainty_minutes(without_covariance, cycle)
        self.assertAlmostEqual(actual, 3.3498, places=3)
        self.assertNotAlmostEqual(actual, uncorrelated, places=2)

    def test_august_13_window_is_observable(self) -> None:
        event = next(row for row in self.events if row["cycle"] == 505)
        self.assertTrue(event["observable"])
        self.assertEqual(event["midpoint_local"][:16], "2026-08-13T00:15")
        self.assertEqual(event["sequence_start_local"][:16], "2026-08-12T21:54")
        self.assertEqual(event["sequence_end_local"][:16], "2026-08-13T02:36")
        self.assertGreater(event["minimum_target_altitude_deg"], 45.0)
        self.assertLess(event["maximum_sun_altitude_deg"], -18.0)
        self.assertTrue(event["moon_below_horizon_throughout"])

    def test_outputs_include_csv_json_markdown_and_calendar(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            output_dir = Path(directory)
            summary = planner.write_outputs(
                output_dir,
                self.events,
                self.ephemeris,
                self.ephemeris_path,
                self.validation_path,
                date(2026, 8, 4),
                date(2026, 10, 1),
                self.duration,
                planner.PlanningLimits(),
                planner.WORKING_TIMEZONE,
            )
            self.assertGreater(summary["observable_event_count"], 0)
            for name in (
                "README.md",
                "transit_windows.csv",
                "observation_plan.json",
                "observable_transits.ics",
            ):
                self.assertTrue((output_dir / name).is_file())

            saved = json.loads(
                (output_dir / "observation_plan.json").read_text(encoding="utf-8")
            )
            self.assertEqual(saved["adopted_ephemeris"]["events"], 27)
            calendar = (output_dir / "observable_transits.ics").read_text(
                encoding="utf-8"
            )
            self.assertIn("BEGIN:VCALENDAR", calendar)
            self.assertIn("toi3505-cycle-505", calendar)
            self.assertIn("DTSTART:20260813T015416Z", calendar)
            raw_calendar = (output_dir / "observable_transits.ics").read_bytes()
            self.assertIn(b"\r\n", raw_calendar)
            self.assertLessEqual(
                max(len(line) for line in raw_calendar.split(b"\r\n")), 75
            )


if __name__ == "__main__":
    unittest.main()
