from __future__ import annotations

import sys
import unittest
from pathlib import Path

from astropy.io import fits


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from plate_solve_toi3505 import merge_wcs_header, target_coordinates


class PlateSolveTests(unittest.TestCase):
    def test_target_coordinates_reads_ra_as_hours(self) -> None:
        header = fits.Header()
        header["RAOBJ2K"] = 19.80289722222222
        header["DECOBJ2K"] = 18.69891388888889

        target = target_coordinates(header)

        self.assertAlmostEqual(target.ra.degree, 297.0434583333333, places=10)
        self.assertAlmostEqual(target.dec.degree, 18.69891388888889, places=10)

    def test_merge_wcs_header_keeps_science_header_and_marks_solution(self) -> None:
        source = fits.Header()
        source["OBJECT"] = "TOI 3505.01"
        source["CHECKSUM"] = "old"
        source["DATASUM"] = "old"
        wcs = fits.Header()
        wcs["CTYPE1"] = "RA---TAN-SIP"
        wcs["CTYPE2"] = "DEC--TAN-SIP"
        wcs["CRVAL1"] = 297.08
        wcs["CRVAL2"] = 18.59

        merged = merge_wcs_header(source, wcs, submission_id=12345)

        self.assertEqual(merged["OBJECT"], "TOI 3505.01")
        self.assertEqual(merged["CTYPE1"], "RA---TAN-SIP")
        self.assertTrue(merged["PLTSOLVD"])
        self.assertEqual(merged["PLTSUBID"], "12345")
        self.assertNotIn("CHECKSUM", merged)
        self.assertNotIn("DATASUM", merged)


if __name__ == "__main__":
    unittest.main()
