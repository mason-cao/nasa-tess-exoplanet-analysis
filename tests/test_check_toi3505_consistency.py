from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import check_toi3505_consistency as consistency  # noqa: E402


class ConsistencyCheckTests(unittest.TestCase):
    def test_require_tokens_reports_missing_text(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "example.md"
            path.write_text("canonical value", encoding="utf-8")
            errors: list[str] = []
            consistency.require_tokens(path, ("canonical", "missing"), errors)
            self.assertEqual(len(errors), 1)
            self.assertIn("missing canonical text", errors[0])

    def test_repository_products_match_canonical_outputs(self) -> None:
        self.assertEqual(consistency.check_repository(ROOT), [])


if __name__ == "__main__":
    unittest.main()
