from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from build_toi3505_paper import (  # noqa: E402
    AUTHOR_NAMES,
    PDF_PATH,
    TOKEN_PATTERN,
    collect_values,
    is_valid_pdf,
    render_manuscript,
)


class PaperBuildTests(unittest.TestCase):
    def test_canonical_values_and_status_are_frozen(self) -> None:
        _, record = collect_values()
        canonical = json.loads(
            (
                ROOT
                / "outputs"
                / "toi3505_ephemeris_refined"
                / "ephemeris_refined.json"
            ).read_text(encoding="utf-8")
        )["ephemeris"]["adopted"]
        primary = record["primary_result"]
        self.assertEqual(primary["events"], 27)
        self.assertEqual(primary["sectors"], [14, 41, 54, 81])
        self.assertAlmostEqual(primary["period_days"], canonical["period_days"], 14)
        self.assertAlmostEqual(
            primary["period_error_days"], canonical["period_error_days"], 16
        )
        self.assertIn("remains a planet candidate", record["status_statement"])

    def test_author_line_is_exactly_the_established_seven_authors(self) -> None:
        _, record = collect_values()
        self.assertEqual(record["authors"], list(AUTHOR_NAMES))
        self.assertEqual(len(record["authors"]), 7)
        # The observer named on the 2022 schedule row belongs on the author
        # list under the program's own authorship guidance.
        self.assertIn("Kevin I. Collins", record["authors"])

    def test_paper_record_inventories_build_inputs(self) -> None:
        _, record = collect_values()
        sources = record["source_files"]
        self.assertEqual(len(sources), len(set(sources)))
        for required in (
            "src/build_toi3505_paper.py",
            "paper/TOI-3505.01_manuscript.md",
            "outputs/toi3505_ephemeris_refined/event_times_best_per_sector.csv",
        ):
            self.assertIn(required, sources)
        self.assertTrue(all((ROOT / path).is_file() for path in sources))

    def test_render_is_self_contained_and_has_no_template_tokens(self) -> None:
        document, _ = render_manuscript(
            ROOT / "paper" / "TOI-3505.01_manuscript.md"
        )
        self.assertIsNone(TOKEN_PATTERN.search(document))
        self.assertNotIn("file://", document)
        self.assertEqual(document.count("data:image/svg+xml;base64,"), 6)
        self.assertIn("To our knowledge,", document)
        self.assertIn("does not validate or confirm TOI-3505.01", document)

    def test_unknown_template_token_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            template = Path(directory) / "bad.md"
            template.write_text("{{NOT_A_CANONICAL_VALUE}}", encoding="utf-8")
            with self.assertRaisesRegex(RuntimeError, "Unknown manuscript tokens"):
                render_manuscript(template)

    def test_committed_pdf_is_present_and_nontrivial(self) -> None:
        self.assertTrue(is_valid_pdf(PDF_PATH))
        self.assertGreater(PDF_PATH.stat().st_size, 100_000)

    def test_pdf_signature_rejects_small_or_non_pdf_files(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            too_small = Path(directory) / "small.pdf"
            too_small.write_bytes(b"%PDF-1.4")
            self.assertFalse(is_valid_pdf(too_small))

            not_pdf = Path(directory) / "large.pdf"
            not_pdf.write_bytes(b"plain text" * 200)
            self.assertFalse(is_valid_pdf(not_pdf))


if __name__ == "__main__":
    unittest.main()
