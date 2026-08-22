from __future__ import annotations

import sys
import unittest
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import analyze_toi3505_ephemeris_robustness as robustness  # noqa: E402


class EphemerisRobustnessTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.summary, cls.tables, cls.events = robustness.build_analysis()

    def test_primary_fit_exactly_reproduces_the_canonical_fit(self) -> None:
        canonical = robustness.load_json(robustness.EPHEMERIS_PATH)
        adopted = canonical["ephemeris"]["adopted"]
        primary = self.summary["primary_tess_linear"]
        self.assertEqual(primary["events"], 27)
        np.testing.assert_allclose(
            primary["period_days"], adopted["period_days"], rtol=0.0, atol=1e-12
        )
        np.testing.assert_allclose(
            primary["period_error_days"],
            adopted["period_error_days"],
            rtol=0.0,
            atol=5e-14,
        )
        np.testing.assert_allclose(
            primary["epoch_bjd_tdb"],
            adopted["epoch_bjd_tdb"],
            rtol=0.0,
            atol=2e-9,
        )

    def test_one_conservative_anchor_represents_each_sector(self) -> None:
        anchors = self.tables["sector_anchors.csv"]
        self.assertEqual([row["sector"] for row in anchors], [14, 41, 54, 81])
        for row in anchors:
            self.assertGreaterEqual(
                row["adopted_error_days"], row["formal_error_days"]
            )
            self.assertGreaterEqual(
                row["adopted_error_days"], row["empirical_error_days"]
            )

    def test_sector_anchor_fit_is_consistent_and_conservative(self) -> None:
        primary = self.summary["primary_tess_linear"]
        anchor = self.summary["four_sector_anchor_fit"]
        combined_error = np.hypot(
            primary["period_error_days"], anchor["period_error_days"]
        )
        self.assertLess(abs(primary["period_days"] - anchor["period_days"]), combined_error)
        self.assertGreater(anchor["period_error_days"], primary["period_error_days"])

    def test_four_sector_jackknife_is_labeled_and_larger_than_formal(self) -> None:
        primary = self.summary["primary_tess_linear"]
        jackknife = self.summary["delete_one_sector_jackknife"]
        self.assertEqual(jackknife["clusters"], 4)
        self.assertGreater(
            jackknife["jackknife_standard_error_days"], primary["period_error_days"]
        )
        self.assertIn("Sensitivity diagnostic only", jackknife["interpretation"])

    def test_quadratic_model_is_not_favored(self) -> None:
        quadratic = self.summary["quadratic_model_control"]
        self.assertGreater(quadratic["delta_bic_quadratic_minus_linear"], 0.0)
        self.assertLess(quadratic["delta_bic_quadratic_minus_linear"], 2.0)
        self.assertIn("not favored", quadratic["assessment"])

    def test_selection_thresholds_leave_the_period_stable(self) -> None:
        primary = self.summary["primary_tess_linear"]
        rows = self.tables["selection_sensitivity.csv"]
        self.assertEqual([row["events"] for row in rows], [27, 26, 23])
        for row in rows:
            self.assertLess(
                abs(row["period_shift_days"]), primary["period_error_days"]
            )

    def test_muscat2_is_an_external_control_not_primary(self) -> None:
        primary = self.summary["primary_tess_linear"]
        external = self.summary["external_muscat2_control"]
        self.assertEqual(external["fit"]["events"], 28)
        self.assertEqual(external["point"]["cycle"], 119)
        self.assertIn("external control only", external["point"]["role"])
        self.assertLess(
            abs(external["fit"]["period_days"] - primary["period_days"]),
            primary["period_error_days"],
        )


if __name__ == "__main__":
    unittest.main()
