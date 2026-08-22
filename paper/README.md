# TOI-3505.01 manuscript

`TOI-3505.01_manuscript.md` is the editable scientific source. Quantitative
values use named `{{TOKENS}}`; `src/build_toi3505_paper.py` replaces those
tokens from canonical JSON and CSV analysis products, embeds the four SVG
figures, and writes a self-contained HTML manuscript.

Build the HTML and machine-readable paper record with:

```bash
.venv/bin/python src/build_toi3505_paper.py
```

Build the HTML and the standard PDF with:

```bash
.venv/bin/python src/build_toi3505_paper.py --pdf-output default
```

The final PDF is written to
`output/pdf/TOI-3505.01_research_paper.pdf`. Render every PDF page for visual
review before distribution:

```bash
mkdir -p tmp/pdfs
pdftoppm -png -r 150 output/pdf/TOI-3505.01_research_paper.pdf \
  tmp/pdfs/toi3505-paper
```

The paper is intentionally limited to a candidate-level timing and archival
result. It does not claim validation, confirmation, a final dilution
correction, or a planet radius. Refresh the dated literature and ExoFOP status
snapshots immediately before submission.
