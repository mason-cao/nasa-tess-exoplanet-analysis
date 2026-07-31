from __future__ import annotations

import sys
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import refine_toi3505_ephemeris as refine


def synthetic_events(
    period: float, epoch: float, cycles: np.ndarray, error_days: float
) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "sector": 54,
            "cycle": cycles,
            "measured_bjd": epoch + cycles * period,
            "time_error_days": error_days,
            "used_in_ephemeris": True,
        }
    )


class EphemerisFitTests(unittest.TestCase):
    def test_recovers_a_noiseless_line(self) -> None:
        period, epoch = 2.9151516, 2459793.5311
        cycles = np.arange(-380, 254, 37, dtype=float)
        fit = refine.fit_ephemeris(
            synthetic_events(period, epoch, cycles, 0.002), "noiseless"
        )
        self.assertAlmostEqual(fit.period_days, period, places=10)
        self.assertAlmostEqual(fit.epoch_bjd, epoch, places=6)

    def test_errors_are_inflated_exactly_once(self) -> None:
        """Regression: the tools helper already inflates by reduced chi-square.

        Scaling its output again would count the penalty twice.  This fit builds
        its own covariance, so a data set whose residuals match its error bars
        must come back with ``error_scale`` of exactly one.
        """
        period, epoch = 2.9151516, 2459793.5311
        cycles = np.arange(-30, 31, dtype=float)
        events = synthetic_events(period, epoch, cycles, 0.002)
        fit = refine.fit_ephemeris(events, "exact line")
        self.assertEqual(fit.error_scale, 1.0)
        self.assertLess(fit.reduced_chi2, 1e-12)

        # Formal error of a weighted line fit: sigma_P = sigma_t / sqrt(Sxx).
        expected = 0.002 / np.sqrt(np.sum((cycles - cycles.mean()) ** 2))
        self.assertAlmostEqual(fit.period_days, period, places=10)
        np.testing.assert_allclose(fit.period_error_days, expected, rtol=1e-9)

    def test_scatter_beyond_the_error_bars_widens_the_period(self) -> None:
        period, epoch = 2.9151516, 2459793.5311
        cycles = np.arange(-30, 31, dtype=float)
        events = synthetic_events(period, epoch, cycles, 0.002)
        rng = np.random.default_rng(3505)
        events["measured_bjd"] += rng.normal(0.0, 0.006, size=len(cycles))
        fit = refine.fit_ephemeris(events, "scattered")
        self.assertGreater(fit.reduced_chi2, 1.0)
        self.assertAlmostEqual(fit.error_scale, np.sqrt(fit.reduced_chi2), places=12)


class AdoptedResultTests(unittest.TestCase):
    """Guard the numbers the posters quote."""

    @classmethod
    def setUpClass(cls) -> None:
        path = ROOT / "outputs" / "toi3505_ephemeris_refined" / "ephemeris_refined.json"
        if not path.exists():
            raise unittest.SkipTest("run src/refine_toi3505_ephemeris.py first")
        import json

        cls.summary = json.loads(path.read_text())

    def test_partial_coverage_event_is_rejected(self) -> None:
        """Sector 81 cycle 247 straddles a gap and fits a 24 ppt event."""
        events = pd.read_csv(
            ROOT
            / "outputs"
            / "toi3505_ephemeris_refined"
            / "event_times_best_per_sector.csv"
        )
        bad = events[(events["sector"] == 81) & (events["cycle"] == 247)]
        self.assertEqual(len(bad), 1)
        self.assertLess(float(bad["window_coverage"].iloc[0]), 0.5)
        self.assertFalse(bool(bad["used_in_ephemeris"].iloc[0]))

    def test_accepted_events_are_all_well_covered(self) -> None:
        events = pd.read_csv(
            ROOT
            / "outputs"
            / "toi3505_ephemeris_refined"
            / "event_times_best_per_sector.csv"
        )
        accepted = events[events["used_in_ephemeris"]]
        self.assertGreaterEqual(float(accepted["window_coverage"].min()), 0.70)
        # No accepted event should have a depth far from the roughly 3 ppt signal.
        self.assertLess(float(accepted["depth_ppt"].max()), 8.0)

    def test_trapezoid_beats_the_box_on_goodness_of_fit(self) -> None:
        fits = self.summary["ephemeris"]
        self.assertLess(fits["adopted"]["reduced_chi2"], 2.0)
        self.assertGreater(fits["box_shape_control"]["reduced_chi2"], 10.0)

    def test_adopted_period_agrees_with_both_published_values(self) -> None:
        comparisons = self.summary["comparisons"]
        self.assertLess(comparisons["sigma_from_catalog"], 3.0)
        self.assertLess(comparisons["sigma_from_spoc"], 3.0)
        self.assertLess(comparisons["sigma_adopted_from_qlp_only"], 3.0)
        self.assertGreater(comparisons["precision_gain_over_catalog"], 1.0)
        self.assertGreater(comparisons["precision_gain_over_spoc"], 1.0)

    def test_every_sector_contributes(self) -> None:
        self.assertEqual(self.summary["ephemeris"]["adopted"]["sectors"], [14, 41, 54, 81])


if __name__ == "__main__":
    unittest.main()
