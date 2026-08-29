# TOI-3505.01 whole-sequence transit search

The fixed-window test asks whether the 2022 schedule window contained a
dimming. Because that window is now known to be stale, this scan asks the
broader question: is there anything like the expected transit anywhere in
the sequence? The duration is held at a published TESS value and the
midpoint is scanned across the observed span with a free depth.

## Result

### TOI catalog duration (2.004 h)

- Searched midpoints 1.44 to 3.47 h after the first exposure: 62 trials at 2-minute spacing. The range is set by requiring the event to be at least 85% sampled with at least 20 baseline points on each side.
- Deepest dimming: 0.930 +/- 0.423 ppt at 1.74 h (2.20 sigma).
- Median formal 3.0-sigma WLS upper bound 0.989 ppt; weakest 2.357 ppt, against a published 2.910 ppt depth. These formal bounds do not include temporal covariance and are not treated as global exclusions.
- Injecting the published depth into the observed curve at each admissible midpoint produces a formal >=3-sigma fitted depth at 50/62 midpoints (80.6%); the minimum injected-event formal S/N is 0.68.
- Because exact-depth recovery is not complete across phase, this analysis does not claim that the published depth is excluded at every midpoint.

### SPOC multi-sector duration (2.711 h)

- Searched midpoints 1.77 to 3.14 h after the first exposure: 42 trials at 2-minute spacing. The range is set by requiring the event to be at least 85% sampled with at least 20 baseline points on each side.
- Deepest dimming: 0.422 +/- 0.397 ppt at 2.04 h (1.06 sigma).
- Median formal 3.0-sigma WLS upper bound 0.685 ppt; weakest 1.614 ppt, against a published 3.292 ppt depth. These formal bounds do not include temporal covariance and are not treated as global exclusions.
- Injecting the published depth into the observed curve at each admissible midpoint produces a formal >=3-sigma fitted depth at 40/42 midpoints (95.2%); the minimum injected-event formal S/N is 2.84.
- Because exact-depth recovery is not complete across phase, this analysis does not claim that the published depth is excluded at every midpoint.

## Limits

- The scan holds the duration fixed and assumes a symmetric box; it is
  not a limb-darkened physical transit fit.
- Depths are observed aperture depths. No dilution correction is
  applied, and the 0.517-arcsecond companion is not resolved.
- Formal WLS intervals are locally scaled by reduced chi-square but
  do not model time-correlated residuals. No false-alarm probability
  is quoted.
- The injection-recovery percentage is a phase-sampling diagnostic
  for this observed residual realization, not a calibrated detection
  probability or an ensemble completeness estimate.
- The sequence is only about twice the transit duration, so a fully
  sampled event with baseline on both sides fits only near the middle
  of the night. The search says nothing about events outside that
  range, including any that fall entirely outside the sequence.
- The degraded final half hour of the night, where the scatter is
  about four times the mid-run value, produces a large spurious
  dimming if trial events are allowed to hang off the end of the
  sequence. The two-sided baseline requirement excludes those trials.

## Products

- `ground_search.json` - inputs, per-duration summaries, and limits.
- `floating_time_scan.csv` - every trial midpoint and fitted depth.
- `01_ground_transit_search.png` and `.svg` - publication figure.

Regenerate with:

```bash
.venv/bin/python src/search_toi3505_ground_transit.py
```
