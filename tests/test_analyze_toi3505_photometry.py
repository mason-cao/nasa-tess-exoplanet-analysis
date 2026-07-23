from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from analyze_toi3505_photometry import (
    closest_predicted_midpoint,
    normalized,
    validate_native_working_table,
)


class PhotometryReviewTests(unittest.TestCase):
    def test_normalized_curve_has_median_one(self) -> None:
        curve = normalized(np.array([8.0, 10.0, 12.0]))

        self.assertAlmostEqual(float(np.median(curve)), 1.0)

    def test_stated_ephemeris_places_midpoint_before_images(self) -> None:
        result = closest_predicted_midpoint(
            start_bjd=2459782.598234811,
            end_bjd=2459782.809458706,
        )

        self.assertFalse(result["inside_observation"])
        self.assertEqual(result["epoch_number"], -4)
        self.assertAlmostEqual(result["midpoint_bjd_tdb"], 2459781.8737626)

    def test_saved_working_table_matches_selected_comparison_counts(self) -> None:
        rows = 281
        source = pd.DataFrame(
            {
                "Label": [f"image-{number:04d}.fits" for number in range(1, rows + 1)],
                "slice": np.arange(1, rows + 1),
                "Saturated": np.zeros(rows),
                "BJD_TDB": np.linspace(2459782.60, 2459782.81, rows),
                "AIRMASS": np.linspace(1.3, 1.1, rows),
                "FWHM_Mean": np.full(rows, 12.0),
                "Source_Radius": np.full(rows, 35.0),
                "Sky_Rad(min)": np.full(rows, 70.0),
                "Sky_Rad(max)": np.full(rows, 139.0),
                "Source-Sky_T1": np.linspace(500_000.0, 510_000.0, rows),
                "Source_Error_T1": np.full(rows, 2_000.0),
                "Sky/Pixel_T1": np.full(rows, 800.0),
                "Width_T1": np.full(rows, 20.0),
                "X(FITS)_T1": np.full(rows, 1_800.0),
                "Y(FITS)_T1": np.full(rows, 1_700.0),
            }
        )
        working_stars = ["C2", "C3", "C5"]
        for index, star in enumerate(working_stars, start=1):
            source[f"Source-Sky_{star}"] = np.linspace(
                300_000.0 + index * 10_000.0,
                310_000.0 + index * 10_000.0,
                rows,
            )

        comparison_counts = sum(
            source[f"Source-Sky_{star}"].to_numpy() for star in working_stars
        )
        ratio = source["Source-Sky_T1"].to_numpy() / comparison_counts
        native = source[
            [
                "Label",
                "slice",
                "Saturated",
                "BJD_TDB",
                "AIRMASS",
                "FWHM_Mean",
                "Source_Radius",
                "Sky_Rad(min)",
                "Sky_Rad(max)",
                "Source-Sky_T1",
                "Source_Error_T1",
                "Sky/Pixel_T1",
                "Width_T1",
                "X(FITS)_T1",
                "Y(FITS)_T1",
            ]
        ].copy()
        native["rel_flux_T1"] = ratio
        native["tot_C_cnts"] = comparison_counts

        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "working_measurements.xls"
            native.to_csv(path, sep="\t", index=False)
            result = validate_native_working_table(
                path,
                source,
                working_stars,
                expected_outer_radius=139,
            )

        self.assertTrue(result["present"])
        self.assertEqual(result["rows"], rows)
        self.assertLess(result["maximum_normalized_curve_difference"], 1.0e-12)


if __name__ == "__main__":
    unittest.main()
