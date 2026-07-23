from __future__ import annotations

import sys
import unittest
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from align_toi3505 import (
    StarMeasurement,
    aligned_name,
    measure_shift_from_stars,
    shift_whole_pixels,
)


def make_star_field() -> tuple[np.ndarray, list[StarMeasurement]]:
    yy, xx = np.indices((256, 256), dtype=float)
    image = np.full((256, 256), 500.0, dtype=np.float32)
    stars: list[StarMeasurement] = []
    for x, y, brightness in (
        (55.4, 61.7, 12_000.0),
        (94.2, 181.5, 10_000.0),
        (132.8, 115.3, 11_000.0),
        (178.4, 70.6, 9_000.0),
        (200.1, 192.2, 10_500.0),
        (71.7, 126.8, 8_500.0),
    ):
        image += brightness * np.exp(
            -((xx - x) ** 2 + (yy - y) ** 2) / (2.0 * 3.5**2)
        ).astype(np.float32)
        stars.append(
            StarMeasurement(
                x=x,
                y=y,
                flux=1.0,
                fwhm=8.24,
                peak_fraction=0.01,
            )
        )
    return image, stars


class AlignmentTests(unittest.TestCase):
    def test_whole_pixel_shift_moves_values_without_interpolation(self) -> None:
        image = np.arange(30, dtype=np.float32).reshape(5, 6)
        shifted = shift_whole_pixels(
            image,
            y_shift=-1,
            x_shift=2,
            fill_value=-1.0,
        )
        np.testing.assert_array_equal(shifted[0:4, 2:6], image[1:5, 0:4])
        self.assertTrue(np.all(shifted[:, :2] == -1.0))
        self.assertTrue(np.all(shifted[4, :] == -1.0))

    def test_star_centers_recover_integer_shift(self) -> None:
        reference, stars = make_star_field()
        moving = shift_whole_pixels(
            reference,
            y_shift=7,
            x_shift=-5,
            fill_value=500.0,
        )
        result = measure_shift_from_stars(
            moving,
            stars,
            approximate_y_shift=-8,
            approximate_x_shift=8,
        )
        self.assertIsNotNone(result)
        assert result is not None
        self.assertEqual(result["x_shift"], 5)
        self.assertEqual(result["y_shift"], -7)
        self.assertGreaterEqual(result["stars_used"], 5)
        self.assertLess(result["position_rms"], 0.6)

    def test_aligned_filename_matches_course_suffix(self) -> None:
        path = Path("TOI_3505.01_50.000s_R-0001_out.fits")
        self.assertEqual(
            aligned_name(path),
            "TOI_3505.01_50.000s_R-0001_aligned.fits",
        )


if __name__ == "__main__":
    unittest.main()
