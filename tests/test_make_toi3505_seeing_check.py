from __future__ import annotations

import csv
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from make_toi3505_seeing_check import (
    crop_with_fits_extent,
    nearest_other_source,
)


class SeeingCheckTests(unittest.TestCase):
    def test_crop_extent_uses_fits_pixel_coordinates(self) -> None:
        image = np.arange(100, dtype=float).reshape(10, 10)

        crop, extent = crop_with_fits_extent(image, 6.0, 5.0, half_width=1)

        np.testing.assert_array_equal(crop, image[3:6, 4:7])
        self.assertEqual(extent, (4.5, 7.5, 3.5, 6.5))

    def test_nearest_other_source_skips_the_target(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "sources.csv"
            with path.open("w", newline="", encoding="utf-8") as file_obj:
                writer = csv.DictWriter(
                    file_obj,
                    fieldnames=["x_pixel", "y_pixel", "width_pixels"],
                )
                writer.writeheader()
                writer.writerows(
                    [
                        {"x_pixel": 99.0, "y_pixel": 99.0, "width_pixels": 10.0},
                        {"x_pixel": 129.0, "y_pixel": 99.0, "width_pixels": 9.5},
                        {"x_pixel": 159.0, "y_pixel": 99.0, "width_pixels": 9.0},
                    ]
                )

            source = nearest_other_source(
                path,
                center_x_fits=100.0,
                center_y_fits=100.0,
            )

        self.assertAlmostEqual(source["distance_pixels"], 30.0)
        self.assertAlmostEqual(source["width_pixels"], 9.5)


if __name__ == "__main__":
    unittest.main()
