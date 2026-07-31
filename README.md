# TOI-3505.01 TESS and Ground-Based Analysis

This repository contains my NASA data-science internship analysis of
TOI-3505.01. It combines a 2022 GMU 0.8 m R-band observing sequence with
public TESS photometry, pixel data, and SPOC validation products.

## Current result

- 281 usable 50-second ground exposures were measured; 263 pass the documented
  image-quality review used for the final plot.
- The selected AstroImageJ settings are a 25-pixel source aperture, a
  70–139-pixel sky annulus, and comparison stars C2–C11.
- The final ground light curve has 2.497 ppt robust scatter. No detrending or
  transit model was applied because no tested correction passed the frozen
  selection rules and no clear transit was measured.
- The current ephemeris places the nearest transit midpoint 17.39 hours before
  the ground sequence. A recovered historical schedule window falls within the
  observations only under a documented local-time interpretation; its fixed
  check gives −0.658 ± 0.395 ppt, not a significant dimming.
- The 2.915-day signal is independently recovered in TESS Sectors 14, 41, 54,
  and 81. The ground sequence falls inside a Sector 54 data gap, so there is no
  simultaneous TESS comparison.
- A weighted linear fit to 27 mid-transit times across those four sectors spans
  5.05 years and gives P = 2.9151516 ± 0.0000049 days, 2.40 times more precise
  than the TOI catalog period and 1.43 times more precise than the SPOC
  multi-sector fit, consistent with both at 0.32 and 0.71 sigma. Timing errors
  were inflated by sqrt(chi2/dof) = 1.14 before quoting. See
  [`outputs/toi3505_ephemeris_refined`](outputs/toi3505_ephemeris_refined) for
  the three changes that produced it: a trapezoidal transit shape instead of a
  box, SPOC 2-minute data for Sectors 54 and 81, and the rejection of one
  partially covered event.
- All 280 catalog stars within 2.5 arcminutes that are bright enough to mimic
  the signal were measured on the ground images; 222 are cleared. Two fall
  inside the target aperture and cannot be separated.

These results support a carefully scoped light-curve and timing analysis. They
do not validate the planet or resolve the known 0.517-arcsecond companion.

## Final Discord package

The three files prepared for mentor review are in
[`outputs/toi3505_discord_post`](outputs/toi3505_discord_post):

1. `01_TOI_3505.01_final_light_curve.png`
2. `02_TOI_3505.01_data_set_fit_settings.png`
3. `03_TOI_3505.01_NEB_screen.png`

The ready-to-post message is kept locally in `docs/discord-messages.md`.

## Symposium poster

Two 48 x 36 inch boards and their figures are in
[`outputs/toi3505_poster`](outputs/toi3505_poster). `v1` is the detailed board
and `v2` the simpler one built from the AstroImageJ products the rest of the
cohort shows. See that folder's README for the regeneration steps and for the
unpublished ExoFOP follow-up this project does not analyse.

## Repository layout

- `src/` — calibration, alignment, photometry, timing, TESS, and validation
  scripts. `refine_toi3505_ephemeris.py` produces the adopted ephemeris and must
  run before the poster scripts.
- `tests/` — unit tests for the analysis and reproducibility helpers.
- `data/ground/toi3505/` — local calibrated and aligned FITS products. Large
  generated FITS files are ignored by Git.
- `data/tess/toi3505/` — public TESS products and compact archive records.
- `data/catalogs/toi3505/` — saved Gaia and TIC catalog queries.
- `data/program_records/toi3505/` — the recovered observing-schedule record.
- `outputs/` — generated results, grouped by analysis stage.
- `data_and_lectures/` — original TOI-3505 archives and program references.
- `docs/project-record.md` — concise findings, limitations, and remaining work.
  Kept locally and not tracked in Git.

## Reproduce the analysis

Create an environment and install the dependencies:

```bash
python3 -m venv .venv
.venv/bin/pip install -r requirements.txt
```

Run the current local analysis in its documented order:

```bash
.venv/bin/python src/run_toi3505_analysis.py
```

Useful options include `--download`, `--remeasure-ground-apertures`,
`--skip-nearby-images`, and `--skip-large-manifest`. Run the tests alone with:

```bash
.venv/bin/python -m unittest discover -s tests
```

## Scientific guardrails

- Keep the ground transit fit off unless a mentor confirms the historical
  timing source and the data support a fit.
- Treat the nearby-star analysis as a screen, not formal planet validation.
- Treat SPOC comparisons as checks of the same TESS observations, not
  independent detections.
- Preserve excluded frames and their written image-quality reasons.
