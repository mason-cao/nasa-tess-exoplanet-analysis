from __future__ import annotations

import sys
import unittest
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from check_toi3505_aperture_radii import measure_frame, primary_radius_passes


class ApertureRadiusCheckTests(unittest.TestCase):
    def test_background_subtracted_aperture_recovers_known_signal(self) -> None:
        image = np.full((101, 101), 100.0)
        image[50, 50] += 1_000.0

        flux, background = measure_frame(
            image,
            np.array([[50.0, 50.0]]),
            radii=(3.0, 5.0),
            annulus_inner=20.0,
            annulus_outer=30.0,
        )

        np.testing.assert_allclose(background, [100.0], atol=1.0e-10)
        np.testing.assert_allclose(flux[:, 0], [1_000.0, 1_000.0], atol=1.0e-8)

    def test_primary_radius_passes_when_on_both_precision_plateaus(self) -> None:
        metrics = pd.DataFrame(
            {
                "source_radius_pixels": [30.0, 35.0, 40.0],
                "target_robust_scatter_ppt": [3.0, 3.2, 3.4],
                "median_pseudo_target_scatter_ppt": [4.2, 4.0, 3.9],
            }
        )

        self.assertTrue(primary_radius_passes(metrics))

    def test_primary_radius_fails_when_target_precision_is_too_far_from_best(self) -> None:
        metrics = pd.DataFrame(
            {
                "source_radius_pixels": [30.0, 35.0, 40.0],
                "target_robust_scatter_ppt": [3.0, 3.5, 3.4],
                "median_pseudo_target_scatter_ppt": [4.2, 4.0, 3.9],
            }
        )

        self.assertFalse(primary_radius_passes(metrics))


if __name__ == "__main__":
    unittest.main()
