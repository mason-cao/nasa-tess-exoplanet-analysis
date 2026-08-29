from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import search_toi3505_ground_transit as search  # noqa: E402


RESULT_PATH = ROOT / "outputs" / "toi3505_ground_search" / "ground_search.json"


class GroundTransitSearchTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.result = json.loads(RESULT_PATH.read_text(encoding="utf-8"))
        cls.durations = {
            str(entry["label"]): entry for entry in cls.result["durations"]
        }

    def test_both_published_durations_are_searched(self) -> None:
        self.assertEqual(set(self.durations), {"TOI catalog", "SPOC multi-sector"})

    def test_formal_bounds_are_explicitly_not_global_exclusions(self) -> None:
        for label, entry in self.durations.items():
            with self.subTest(label=label):
                self.assertTrue(
                    entry["formal_upper_limit_below_expected_depth_everywhere"]
                )
                self.assertLess(
                    entry["maximum_formal_upper_limit_ppt"],
                    entry["expected_depth_ppt"],
                )
                self.assertIn("no global exclusion", entry["interpretation"])

    def test_no_single_midpoint_reaches_formal_three_sigma(self) -> None:
        for label, entry in self.durations.items():
            with self.subTest(label=label):
                self.assertLess(entry["best_depth_snr"], 3.0)

    def test_expected_depth_recovery_is_reported_without_overstatement(self) -> None:
        scan = pd.read_csv(
            ROOT / "outputs" / "toi3505_ground_search" / "floating_time_scan.csv"
        )
        for label, entry in self.durations.items():
            with self.subTest(label=label):
                subset = scan.loc[scan["duration_label"] == label]
                recovered = subset["injected_above_formal_3sigma"].astype(bool)
                self.assertEqual(
                    int(recovered.sum()),
                    entry["expected_depth_formal_recovery_count"],
                )
                self.assertEqual(
                    len(subset), entry["expected_depth_formal_recovery_total"]
                )
                self.assertAlmostEqual(
                    float(recovered.mean()),
                    entry["expected_depth_formal_recovery_fraction"],
                )
                self.assertFalse(
                    entry["expected_depth_recovered_above_formal_3sigma_everywhere"]
                )
                self.assertGreater(
                    entry["expected_depth_formal_recovery_fraction"], 0.5
                )
                self.assertLess(entry["expected_depth_formal_recovery_fraction"], 1.0)
                self.assertAlmostEqual(
                    float(subset["injected_recovered_increment_ppt"].median()),
                    entry["expected_depth_ppt"],
                    delta=0.02,
                )

    def test_every_trial_keeps_baseline_on_both_sides(self) -> None:
        scan = pd.read_csv(
            ROOT / "outputs" / "toi3505_ground_search" / "floating_time_scan.csv"
        )
        minimum = search.MINIMUM_BASELINE_POINTS_EACH_SIDE
        self.assertTrue((scan["baseline_points_before"] >= minimum).all())
        self.assertTrue((scan["baseline_points_after"] >= minimum).all())
        self.assertTrue(
            (scan["event_coverage_fraction"] >= search.MINIMUM_EVENT_COVERAGE).all()
        )

    def test_an_injected_transit_is_recovered_at_the_injected_time(self) -> None:
        """The scan must find a transit that is genuinely present."""
        curve = pd.read_csv(search.LIGHT_CURVE_PATH)
        time = curve["hours_since_first_image"].to_numpy(dtype=float)
        flux = curve["adopted_relative_brightness"].to_numpy(dtype=float)
        error = curve["raw_relative_brightness_error"].to_numpy(dtype=float)
        use = curve["used_in_primary_curve"].to_numpy(dtype=bool)

        duration = 2.004
        injected_midpoint = 2.5
        injected_depth = 0.00291
        box = search.integrated_box_fraction(
            time, injected_midpoint, duration, search.EXPOSURE_SECONDS / 3600.0
        )
        scan = search.scan_one_duration(
            time, flux * (1.0 - injected_depth * box), error, use, duration
        )
        best = scan.loc[scan["depth_snr"].idxmax()]
        self.assertGreater(float(best["depth_snr"]), 5.0)
        self.assertLess(
            abs(float(best["midpoint_hours_since_first_image"]) - injected_midpoint),
            0.25,
        )


if __name__ == "__main__":
    unittest.main()
