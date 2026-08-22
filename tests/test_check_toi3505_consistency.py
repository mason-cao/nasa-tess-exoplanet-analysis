from __future__ import annotations

import hashlib
import json
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

    def test_manifest_verifier_detects_digest_drift(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            record = root / "outputs" / "toi3505_research_record"
            record.mkdir(parents=True)
            source = root / "result.txt"
            content = b"canonical result\n"
            source.write_bytes(content)
            digest = hashlib.sha256(content).hexdigest()
            (record / "file_manifest.csv").write_text(
                "category,path,size_bytes,modified_utc,sha256\n"
                f"scientific_output,result.txt,{len(content)},2026-08-22T00:00:00Z,{digest}\n",
                encoding="utf-8",
            )
            summary = {
                "files": 1,
                "bytes": len(content),
                "category_counts": {"scientific_output": 1},
                "original_archive_count": 0,
                "sha256_complete": True,
            }
            (record / "manifest_summary.json").write_text(
                json.dumps(summary), encoding="utf-8"
            )
            self.assertEqual(consistency.verify_manifest(root), [])

            source.write_bytes(b"modified result\n")
            errors = consistency.verify_manifest(root)
            self.assertTrue(any("SHA-256 differs" in error for error in errors))


if __name__ == "__main__":
    unittest.main()
