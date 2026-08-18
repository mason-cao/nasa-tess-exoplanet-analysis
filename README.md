# TOI-3505.01 TESS and Ground-Based Analysis

NASA data-science internship analysis combining a 2022 GMU 0.8 m R-band
observing sequence with public TESS photometry, pixel data, and SPOC validation
products.

## Results

- Ground photometry used 263 quality-reviewed exposures and achieved 2.497 ppt
  robust scatter. No clear transit was measured.
- The 2.915-day signal was recovered in TESS Sectors 14, 41, 54, and 81. The
  ground observations fall within a Sector 54 data gap.
- A fit to 27 transit times spanning 5.05 years gives
  **P = 2.9151516 ± 0.0000049 days**. Details are in
  [`outputs/toi3505_ephemeris_refined`](outputs/toi3505_ephemeris_refined).
- Nearby-star screening measured all 280 catalog stars bright enough to mimic
  the signal within 2.5 arcminutes; 222 were cleared, while two sources inside
  the target aperture remain unresolved.

This work refines the light curve and ephemeris but does not independently
validate the planet or resolve its known 0.517-arcsecond companion.

## Current status

The project owner provided the following status update on 2026-08-18:

- The symposium posters have already been presented. Poster-version selection
  and pre-presentation approval are no longer pending tasks.
- The proposed GMU 0.8 m observation on 2026-08-12 to 2026-08-13 did not occur,
  so no new ground-based data are expected from that window.
- The internship program permits the repository to remain public and to carry
  a license.
- Final light-curve feedback has been received: the reviewer did not see a
  transit, said the curve looked good, and noted that the plot key covers some
  of the data.
- The 2022 spreadsheet row and its headers are preserved, and the project owner
  confirmed that its clocks use Eastern time. On that night this means EDT
  (UTC-4).

The remaining project work is archival and scientific review: recover or
formally document any unavailable 2022 prediction source, epoch, uncertainty,
Transit Info file, and clock-sync record; move the plot key; obtain review of
the TESS-timing choices; and run the formal AstroImageJ NEB procedure only if
the program requires it.

## Repository

- `src/` — analysis scripts for calibration, photometry, TESS, timing, and
  validation
- `tests/` — unit tests and reproducibility checks
- `data/` — ground observations, TESS products, catalogs, and program records
- `data_and_lectures/` — program reference material and lectures; the JWST
  lecture added on 2026-08-18 is labeled as training material with a provenance
  sidecar and is not TOI-3505 observation data
- `outputs/` — generated figures and analysis results
- [`outputs/toi3505_poster`](outputs/toi3505_poster) — symposium posters and
  source figures

## Reproduce

```bash
python3 -m venv .venv
.venv/bin/pip install -r requirements.txt
.venv/bin/python src/run_toi3505_analysis.py
```

Run the tests with:

```bash
.venv/bin/python -m unittest discover -s tests
```

The analysis runner also supports `--download`, `--remeasure-ground-apertures`,
`--skip-nearby-images`, and `--skip-large-manifest`.
