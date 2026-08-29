from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import assess_toi3505_false_positive as assessment  # noqa: E402


RESULT_PATH = (
    ROOT / "outputs" / "toi3505_false_positive" / "false_positive_assessment.json"
)


class FalsePositiveAssessmentTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.result = json.loads(RESULT_PATH.read_text(encoding="utf-8"))
        cls.scenarios = cls.result["scenarios"]

    def test_every_standard_scenario_is_addressed(self) -> None:
        self.assertEqual(
            set(self.scenarios),
            {
                "nearby_eclipsing_binary",
                "blended_eclipsing_binary",
                "eclipsing_binary_on_target",
                "unresolved_close_companion",
            },
        )

    def test_every_scenario_records_what_it_cannot_exclude(self) -> None:
        for name, entry in self.scenarios.items():
            with self.subTest(scenario=name):
                self.assertTrue(str(entry["not_excluded"]).strip())

    def test_reported_depths_show_no_wavelength_trend(self) -> None:
        chromatic = self.scenarios["blended_eclipsing_binary"]["chromatic_depth_test"]
        self.assertTrue(chromatic["no_apparent_monotonic_trend"])
        self.assertFalse(chromatic["monotonic_with_wavelength"])
        self.assertIn("not a standard error", chromatic["scale_interpretation"])

    def test_a_monotonic_chromatic_pattern_is_flagged(self) -> None:
        """A depth that climbs steeply with wavelength must be identified."""
        chromatic = assessment.chromatic_depth_test(
            {"g": 2.0, "r": 3.0, "i": 4.0, "z_s": 5.0}
        )
        self.assertFalse(chromatic["no_apparent_monotonic_trend"])
        self.assertTrue(chromatic["monotonic_with_wavelength"])

    def test_stellar_companions_are_disfavoured_by_the_velocity_span(self) -> None:
        velocity = self.scenarios["eclipsing_binary_on_target"][
            "eclipsing_companion_velocity_bound"
        ]
        self.assertTrue(velocity["stellar_companion_disfavoured"])
        self.assertGreater(
            float(velocity["smallest_stellar_scenario_km_s"]),
            2.0 * float(velocity["observed_velocity_span_km_s"]),
        )

    def test_velocity_amplitude_grows_with_companion_mass(self) -> None:
        light = assessment.velocity_amplitude_km_s(0.1, 2.915, 1.25)
        heavy = assessment.velocity_amplitude_km_s(0.6, 2.915, 1.25)
        self.assertGreater(heavy, light)
        # A Jupiter-mass companion at this period stays well under 1 km/s.
        jupiter = assessment.velocity_amplitude_km_s(0.000954, 2.915, 1.25)
        self.assertLess(jupiter, 1.0)

    def test_the_unresolved_companion_is_not_claimed_to_be_addressed(self) -> None:
        entry = self.scenarios["unresolved_close_companion"]
        self.assertIn("Not addressed", str(entry["assessment"]))


if __name__ == "__main__":
    unittest.main()
