from __future__ import annotations

import sys
import unittest
from pathlib import Path

import numpy as np
from astropy.wcs import WCS


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from analyze_toi3505_tess_pixels import (
    CutoutData,
    annulus_mask,
    circular_mask,
    cutouts_equivalent,
    local_background,
    positive_centroid,
)


class TessPixelCheckTests(unittest.TestCase):
    def test_circular_mask_has_expected_integer_pixels(self) -> None:
        mask = circular_mask((7, 7), 3.0, 3.0, 1.0)
        self.assertEqual(int(mask.sum()), 5)
        self.assertTrue(mask[3, 3])

    def test_annulus_excludes_source_circle(self) -> None:
        source = circular_mask((15, 15), 7.0, 7.0, 3.0)
        sky = annulus_mask((15, 15), 7.0, 7.0, 4.0, 6.0)
        self.assertFalse(bool(np.any(source & sky)))
        self.assertGreater(int(sky.sum()), 0)

    def test_local_background_recovers_constant_sky(self) -> None:
        cube = np.full((3, 15, 15), 100.0)
        cube[:, 7, 7] += 500.0
        errors = np.full_like(cube, 2.0)
        cutout = CutoutData(
            path=Path("synthetic.fits"),
            sector=54,
            time_bjd=np.arange(3, dtype=float),
            flux=cube,
            flux_error=errors,
            quality=np.zeros(3, dtype=np.int64),
            cadence_days=1.0 / 144.0,
            wcs=WCS(naxis=2),
            target_x=7.0,
            target_y=7.0,
            pixel_scale_arcsec=21.0,
        )
        background, background_error, mask = local_background(cutout)
        np.testing.assert_allclose(background, 100.0)
        self.assertTrue(np.all(background_error > 0))
        self.assertFalse(mask[7, 7])

    def test_positive_centroid_recovers_known_position(self) -> None:
        image = np.zeros((11, 11), dtype=float)
        image[6, 4] = 10.0
        x, y = positive_centroid(image, 4.0, 6.0, radius=2.0)
        self.assertAlmostEqual(x, 4.0)
        self.assertAlmostEqual(y, 6.0)

    def test_duplicate_cutout_requires_identical_science_arrays(self) -> None:
        first = CutoutData(
            path=Path("first.fits"),
            sector=14,
            time_bjd=np.asarray([1.0, 2.0]),
            flux=np.ones((2, 3, 3)),
            flux_error=np.full((2, 3, 3), 0.1),
            quality=np.zeros(2, dtype=np.int64),
            cadence_days=1.0 / 48.0,
            wcs=WCS(naxis=2),
            target_x=1.0,
            target_y=1.0,
            pixel_scale_arcsec=21.0,
        )
        second = CutoutData(
            path=Path("second.fits"),
            sector=first.sector,
            time_bjd=first.time_bjd.copy(),
            flux=first.flux.copy(),
            flux_error=first.flux_error.copy(),
            quality=first.quality.copy(),
            cadence_days=first.cadence_days,
            wcs=WCS(naxis=2),
            target_x=first.target_x,
            target_y=first.target_y,
            pixel_scale_arcsec=first.pixel_scale_arcsec,
        )
        self.assertTrue(cutouts_equivalent(first, second))
        second.flux[0, 0, 0] = 2.0
        self.assertFalse(cutouts_equivalent(first, second))


if __name__ == "__main__":
    unittest.main()
