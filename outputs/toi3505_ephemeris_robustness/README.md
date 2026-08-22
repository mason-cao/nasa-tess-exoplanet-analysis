# TOI-3505.01 ephemeris robustness

## Primary result

The primary result remains the TESS-only, 27-event linear fit:

`P = 2.9151515671 +/- 0.0000048720 days`.

The quoted uncertainty is the event-level weighted-fit uncertainty after one
sqrt(reduced chi-square) inflation. It is not replaced by any control below.

## Robustness controls

- Four conservative sector anchors give
  `2.9151495461 +/- 0.0000051714` days.
  Each sector error is the larger of its formal weighted-mean error and the
  empirical within-sector standard error.
- Deleting one sector at a time gives a four-cluster jackknife scale of
  `0.0000072894` days. With only
  four clusters this is a sensitivity diagnostic, not a calibrated replacement
  for the primary uncertainty.
- A quadratic timing model changes chi-square by only
  `1.58` and has
  `Delta BIC = +1.71`
  relative to the linear model, so curvature is not favored.
- Adding the public MuSCAT2 timing only as an external control gives
  `P = 2.9151525056 +/- 0.0000046714` days.

The stricter depth-S/N selections, delete-one-event fits, pipeline controls,
and model table are recorded in the accompanying JSON and CSV files.

## Reporting rule

Use the primary TESS-only period in the abstract and headline result. Report
the conservative four-sector fit and delete-one-sector scale as robustness
checks. Do not describe MuSCAT2, QLP, and SPOC as mutually independent
detections: MuSCAT2 is external but tentative, while QLP and SPOC reuse TESS
observations.

Regenerate with:

```bash
.venv/bin/python src/analyze_toi3505_ephemeris_robustness.py
```
