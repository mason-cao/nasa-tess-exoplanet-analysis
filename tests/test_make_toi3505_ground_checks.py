from __future__ import annotations

import sys
import unittest
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from make_toi3505_ground_checks import (
    curve_from_ensemble,
    fit_linear_box,
    integrated_box_fraction,
    noise_vs_bin_size,
)


class GroundCheckTests(unittest.TestCase):
    def test_curve_from_ensemble_is_median_normalized(self) -> None:
        target = np.array([10.0, 11.0, 9.0])
        ensemble = np.array([5.0, 5.0, 5.0])
        curve = curve_from_ensemble(target, ensemble)
        self.assertAlmostEqual(float(np.median(curve)), 1.0)
        self.assertGreater(curve[1], curve[0])

    def test_integrated_box_fraction_handles_event_edges(self) -> None:
        time = np.array([-1.0, -0.5, 0.0, 0.5, 1.0])
        fraction = integrated_box_fraction(
            time, midpoint_hours=0.0, duration_hours=1.0, exposure_hours=0.2
        )
        np.testing.assert_allclose(fraction, [0.0, 0.5, 1.0, 0.5, 0.0])

    def test_linear_box_fit_recovers_synthetic_depth_and_slope(self) -> None:
        time = np.linspace(0.0, 5.0, 301)
        box = integrated_box_fraction(
            time, midpoint_hours=2.5, duration_hours=2.0, exposure_hours=1.0 / 60.0
        )
        depth = 0.004
        flux = 1.0 + 0.001 * (time - np.mean(time)) - depth * box
        error = np.full(len(time), 0.001)
        result = fit_linear_box(time, flux, error, box, np.ones(len(time), dtype=bool))
        self.assertAlmostEqual(result["depth"], depth, places=10)
        self.assertAlmostEqual(result["slope_per_hour"], 0.001, places=10)

    def test_noise_table_has_finite_positive_values(self) -> None:
        rng = np.random.default_rng(3505)
        bjd = 2459000.0 + np.arange(300) / 1440.0
        flux = 1.0 + rng.normal(0.0, 0.002, len(bjd))
        table = noise_vs_bin_size(bjd, flux, np.ones(len(bjd), dtype=bool))
        self.assertGreater(len(table), 3)
        self.assertTrue(np.all(np.isfinite(table["beta"])))
        self.assertTrue(np.all(table["measured_robust_scatter_ppt"] > 0))


if __name__ == "__main__":
    unittest.main()
