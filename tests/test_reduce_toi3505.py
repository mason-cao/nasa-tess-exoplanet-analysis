from __future__ import annotations

import sys
import unittest
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from reduce_toi3505 import (
    MemberRef,
    apply_qa_flags,
    fit_illumination_plane,
    output_name,
)


class ReductionTests(unittest.TestCase):
    def test_fit_illumination_plane_removes_plane_and_preserves_level(self) -> None:
        y = np.linspace(-1.0, 1.0, 192, dtype=np.float32)
        x = np.linspace(-1.0, 1.0, 256, dtype=np.float32)
        plane = 15_000.0 + 1_200.0 * x[None, :] - 700.0 * y[:, None]
        corrected, details = fit_illumination_plane(plane, stride=2)

        self.assertTrue(np.isclose(np.mean(corrected), np.mean(plane), rtol=1e-5))
        self.assertLess(np.std(corrected) / np.mean(corrected), 1e-5)
        self.assertGreater(details["plane_peak_to_peak_fraction"], 0.20)

    def test_output_name_uses_astroimagej_style_suffix(self) -> None:
        ref = MemberRef(
            archive=Path("archive.zip"),
            member="TOI_3505.01/TOI_3505.01_50.000s_R-0034(1).fits",
            name="TOI_3505.01_50.000s_R-0034(1).fits",
            sequence=34,
            crc=0,
            file_size=1,
        )
        self.assertEqual(output_name(ref), "TOI_3505.01_50.000s_R-0034_out.fits")

    def test_tiny_smoke_test_does_not_claim_robust_outliers(self) -> None:
        table = pd.DataFrame(
            {
                "calibrated_median": [100.0, 100.0, 500.0],
                "calibrated_mad": [5.0, 5.0, 50.0],
                "calibrated_p99": [200.0, 200.0, 900.0],
                "calibrated_negative_fraction_sampled": [0.0, 0.0, 0.0],
                "calibrated_finite_fraction": [1.0, 1.0, 1.0],
            }
        )

        flagged = apply_qa_flags(table)

        self.assertFalse(flagged["qa_flag"].any())


if __name__ == "__main__":
    unittest.main()
