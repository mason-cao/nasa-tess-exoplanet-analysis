from __future__ import annotations

import hashlib
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from download_toi3505_comparison_gaia import gaia_query
from make_toi3505_research_record import claim_evidence_rows, sha256_file
from plate_solve_toi3505_representative import reduced_name, representative_indices


class ReproducibilityTests(unittest.TestCase):
    def test_streaming_sha256_matches_standard_digest(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "small.bin"
            content = b"TOI-3505 reproducibility test\n" * 17
            path.write_bytes(content)
            self.assertEqual(sha256_file(path), hashlib.sha256(content).hexdigest())

    def test_representative_frame_indices_are_distinct_and_include_ends(self) -> None:
        table = pd.DataFrame({"AIRMASS": np.r_[np.linspace(2, 1, 51), np.linspace(1.01, 2, 49)]})
        indices = representative_indices(table, bracket_frames=5)
        self.assertEqual(len(set(indices.values())), 5)
        self.assertEqual(indices["early"], 0)
        self.assertEqual(indices["late"], len(table) - 1)

    def test_reduced_name_conversion_is_exact(self) -> None:
        self.assertEqual(
            reduced_name("TOI_3505.01_50.000s_R-0001_aligned.fits"),
            "TOI_3505.01_50.000s_R-0001_out.fits",
        )
        with self.assertRaises(ValueError):
            reduced_name("wrong.fits")

    def test_targeted_gaia_query_has_one_circle_per_measured_star(self) -> None:
        positions = pd.DataFrame(
            {
                "star": ["T1", "C2"],
                "wcs_ra_deg": [297.0, 296.9],
                "wcs_dec_deg": [18.7, 18.6],
            }
        )
        query = gaia_query(positions)
        self.assertEqual(query.count("CIRCLE('ICRS'"), 2)
        self.assertIn("gaiadr3.gaia_source", query)

    def test_claim_ledger_rows_have_evidence_and_limits(self) -> None:
        rows = claim_evidence_rows()
        self.assertGreaterEqual(len(rows), 8)
        self.assertTrue(all(row["evidence"] for row in rows))
        self.assertTrue(all(row["limit"] for row in rows))


if __name__ == "__main__":
    unittest.main()
