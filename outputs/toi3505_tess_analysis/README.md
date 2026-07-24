# TOI-3505.01 four-sector TESS analysis

This folder contains a reproducible first-pass measurement of the public TESS signal in Sectors 14, 41, 54, and 81.

## Method

- QLP is the common four-sector data set.
- Only finite measurements with QUALITY = 0 are used.
- The fitted shape is a box integrated across each sector's exposure time.
- Every event has its own straight local baseline.
- SPOC and TESS-SPOC products are pipeline checks of the same observations, not extra observations.
- QLP files do not contain a CROWDSAP keyword, so the QLP depths here are observed aperture depths and are not labeled dilution-corrected.

## Sector measurements

| Sector | Cadence (min) | Depth (ppt) | Duration (h) | Midpoint offset (min) | BLS period (d) |
|---:|---:|---:|---:|---:|---:|
| 14 | 30.0 | 3.728 +/- 0.419 | 1.80 | -10.1 | 2.91923 |
| 41 | 10.0 | 3.015 +/- 0.209 | 2.05 | -4.3 | 2.91394 |
| 54 | 10.0 | 2.682 +/- 0.295 | 2.05 | -2.9 | 2.91518 |
| 81 | 3.3 | 3.130 +/- 0.240 | 2.10 | -1.4 | 2.91432 |

## Timing interpretation

25 individual events pass the stated timing-quality rule. The Sector 14-only box ephemeris predicts the catalog cycle nearest the GMU night at BJD_TDB 2459783.504888 +/- 440.9 minutes. This is a modern reconstruction from the public Sector 14 light curve.

The recovered 2022 schedule row lists ingress at 00:15 and egress at 01:54. Under the documented Eastern-time interpretation, those become BJD_TDB 2459782.682401-2459782.751151, inside the GMU images. Its midpoint is 20.23 hours later than the nearest current-catalog prediction.

The schedule period is close to the later measured periods, but the source row contains no epoch or timing uncertainty. The disagreement therefore cannot be attributed to period drift alone. The historical window is a recovered scheduling prediction, not proof that a transit occurred.

The GMU sequence also falls entirely inside a Sector 54 TESS data gap. There are no quality-zero TESS points during the ground window, so a simultaneous flux comparison is not possible.

## Limits

- A box estimate does not measure limb darkening, impact parameter, or a physical planet radius.
- The reported formal errors include the local fit scatter but do not include a complete astrophysical dilution model.
- The variability periodograms are screens; QLP detrending and the simple quadratic section correction can alter long periods.
- Difference imaging and custom-aperture tests are kept in the separate pixel-analysis output.
- Official SPOC per-sector and combined-report values are kept in outputs/toi3505_data_validation. They are same-observation comparisons.

## Output tables

- `sector_measurements.csv`: depth, duration, timing, odd/even, and phase-0.5 checks.
- `period_search.csv`: independent BLS result in each sector.
- `event_times.csv`: one row per measurable transit window.
- `ephemeris_checks.json`: Sector 14-only and four-sector box timing fits.
- `pipeline_comparison.csv`: QLP, SPOC, and TESS-SPOC extraction checks.
- `variability_screen.csv`: out-of-transit residual-period search.
- `clean_light_curves/`: normalized QLP data with quality, phase, and cycle columns.
- `../toi3505_data_validation/`: official Sector 54, Sector 81, and combined-sector SPOC report comparison.

Pipeline-comparison rows: 16. Variability-screen rows: 4.
