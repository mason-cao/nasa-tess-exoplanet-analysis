from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import reconstruct_toi3505_schedule as reconstruction  # noqa: E402


class ScheduleReconstructionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.result = reconstruction.reconstruct(
            reconstruction.load_object(reconstruction.SCHEDULE_CHECK_PATH),
            reconstruction.load_object(reconstruction.GROUND_SUMMARY_PATH),
            reconstruction.load_object(reconstruction.EPHEMERIS_PATH),
            reconstruction.load_object(reconstruction.EXOFOP_PATH),
        )

    def test_old_ephemeris_lands_on_the_schedule_after_96_cycles(self) -> None:
        recovered = self.result["reconstruction"]
        self.assertEqual(recovered["cycles"], 96)
        self.assertTrue(recovered["all_markers_agree_within_one_minute"])
        self.assertLess(recovered["maximum_absolute_marker_offset_seconds"], 60.0)

    def test_visible_revised_period_does_not_generate_the_timing_cells(self) -> None:
        recovered = self.result["reconstruction"]
        self.assertGreater(
            recovered["schedule_minus_displayed_period_midpoint_hours"], 5.0
        )

    def test_current_tess_event_is_outside_the_ground_sequence(self) -> None:
        timeline = self.result["timeline"]
        self.assertLess(
            timeline["adopted_ephemeris_midpoint_bjd_tdb"],
            timeline["observation_start_bjd_tdb"],
        )
        self.assertGreater(timeline["schedule_minus_adopted_midpoint_hours"], 20.0)

    def test_result_preserves_the_inference_limit(self) -> None:
        assessment = self.result["assessment"]
        self.assertIn("not proof", assessment["strength"])
        self.assertTrue(any("workbook" in item for item in assessment["limits"]))

    def test_frozen_exofop_context_is_internally_consistent(self) -> None:
        payload = json.loads(reconstruction.EXOFOP_PATH.read_text(encoding="utf-8"))
        inventory = payload["time_series_inventory"]
        self.assertEqual(inventory["table_rows"], 7)
        self.assertEqual(inventory["unique_observing_nights"], 6)
        self.assertEqual(payload["current_status"]["tfopwg_disposition"], "PC")


if __name__ == "__main__":
    unittest.main()
