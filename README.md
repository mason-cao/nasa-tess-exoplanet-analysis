# TOI-3505.01 TESS and Ground-Based Analysis

NASA data-science internship analysis combining a 2022 GMU 0.8 m R-band
observing sequence with public TESS photometry, pixel data, and SPOC validation
products.

## Results

- Ground photometry used 263 quality-reviewed exposures and achieved 2.497 ppt
  robust scatter. No clear transit was measured.
- Searching the whole sequence rather than the stale schedule window alone
  finds no event above a formal 3-sigma threshold. Exact published-depth
  injections cross that threshold at 50/62 catalog-duration midpoints (80.6%)
  and 40/42 SPOC-duration midpoints (95.2%), so the sequence does not support an
  exclusion at every phase. Details are in
  [`outputs/toi3505_ground_search`](outputs/toi3505_ground_search).
- The 2.915-day signal was recovered in TESS Sectors 14, 41, 54, and 81. The
  ground observations fall within a Sector 54 data gap.
- A fit to 27 transit times spanning 5.05 years gives
  **P = 2.9151516 ± 0.0000049 days**. Details are in
  [`outputs/toi3505_ephemeris_refined`](outputs/toi3505_ephemeris_refined).
- Sector aggregation, event and sector deletion, selection-threshold tests,
  pipeline controls, and one public MuSCAT2 external control support the
  headline period. The formal period error is 4.9e-6 days, while the
  delete-one-sector sensitivity scale is 7.3e-6 days. A quadratic timing term
  is not favored (ΔBIC = +1.71).
- Propagating a public 2021 timing prediction by 96 cycles with its superseded
  2.9174250-day period reproduces all three 2022 schedule markers within 56.1
  seconds. The adopted TESS ephemeris places the relevant midpoint 20.31 hours
  before the stale schedule; applying the schedule's visible revised period to
  the same 2021 midpoint instead gives a 5.24-hour offset.
- The standard false-positive scenarios are each constrained: no measured
  neighbour shows an eclipse-shaped event, the four reported MuSCAT2 depths
  show no apparent monotonic trend across 477-870 nm (a descriptive result
  because per-band errors are unavailable), and an eclipsing stellar companion
  would move the host by at least 12 km/s against a 2.5 km/s observed velocity
  span. The unresolved 0.517-arcsecond companion is not addressed, so this
  supports rather than validates the candidate. Details are in
  [`outputs/toi3505_false_positive`](outputs/toi3505_false_positive).
- Nearby-star screening measured all 409 catalog stars within 2.5 arcminutes
  bright enough to mimic the signal once the TFOP -0.5 mag band correction is
  applied: 288 cleared, 4 blended with the target aperture, 115 too faint for a
  decisive limit, and 2 measured without a deep enough limit. Neither of those
  last two shows an eclipse-shaped event; both are monotonic trends across the
  night.

This work refines the light curve and ephemeris but does not independently
validate the planet or resolve its known 0.517-arcsecond companion.

## Research paper

The finalized analysis manuscript is ready for mentor scientific review:

- [Final research-paper PDF](output/pdf/TOI-3505.01_research_paper.pdf)
- [Editable manuscript](paper/TOI-3505.01_manuscript.md)
- [Self-contained HTML manuscript](paper/TOI-3505.01_manuscript.html)
- [Machine-readable manuscript values](outputs/toi3505_paper/manuscript_values.json)

To our knowledge, this is the first dedicated study to combine all four
publicly available TESS sectors of TOI-3505.01 and trace the 2022 GMU null
result to a schedule window consistent with a superseded ephemeris. This is a
qualified, search-based novelty statement—not a claim of first detection,
classification, follow-up, validation, or confirmation.

The literature search, live ExoFOP status, and NASA Exoplanet Archive candidate
record were refreshed on 2026-08-29. The missing original scheduling-workbook
formulas, unresolved close companion, and absent dilution-corrected radius are
explicit scope limits; they do not block the paper's narrower timing and
archival result. Public imaging and spectroscopy for this target exist on
ExoFOP and are inventoried but not analyzed here. A final mentor review remains
appropriate before submission.

## Repository

- `src/` — analysis scripts for calibration, photometry, TESS, timing, and
  validation
- `tests/` — unit tests and reproducibility checks
- `data/` — ground observations, TESS products, catalogs, and program records
- `data_and_lectures/` — program reference material and lectures; the JWST
  lecture added on 2026-08-18 is labeled as training material with a provenance
  sidecar and is not TOI-3505 observation data
- `outputs/` — generated figures and analysis results
- `paper/` — editable manuscript and self-contained HTML build
- `LICENSE` — MIT for the code and derived products; third-party material
  redistributed here keeps its own terms
- `output/pdf/` — final rendered paper artifact
- [`outputs/toi3505_poster`](outputs/toi3505_poster) — symposium posters and
  source figures

## Reproduce

```bash
python3 -m venv .venv
.venv/bin/pip install -r requirements.txt
.venv/bin/python src/run_toi3505_analysis.py --build-paper-pdf
```

Run the tests with:

```bash
.venv/bin/python -m unittest discover -s tests
```

The analysis runner also supports `--download`, `--remeasure-ground-apertures`,
`--skip-nearby-images`, and `--skip-large-manifest`. Omit
`--build-paper-pdf` when Chrome is unavailable; the self-contained HTML and
machine-readable paper record are still rebuilt.
