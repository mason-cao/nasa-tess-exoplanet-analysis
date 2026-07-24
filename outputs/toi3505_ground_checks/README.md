# TOI-3505.01 ground-data checks

## Comparison ensemble

The final light curve still uses the predeclared sum of all ten comparison
stars. Its robust scatter on the final frame mask is
2.497 ppt. Equal-star and inverse-error
ensembles are saved as robustness comparisons; neither was used to retune the
published curve after looking at the target.

At a 30-minute bin width, the measured robust scatter is
1.309 ppt, compared with a
0.481-ppt white-noise expectation
(beta = 2.72). This is a descriptive time-correlation check,
not a replacement for a transit fit.

## Injection check

Two-hour box dips from 1 to 10 ppt were placed at 13 interior times. Each was
fit with a straight baseline and a fixed box. The control depth at the same
time is retained, so the saved total recovery shows the real phase-dependent
structure rather than an artificially perfect control-subtracted result.
The first tested depth with at least 90% of placements above three sigma was 5 ppt. These trials test this light curve and fitting method; they do
not simulate a target-only point-spread function or establish a corrected TESS
dilution.

## Star catalog check

11 of 11 measured positions have a Gaia source
within 3 arcseconds in the targeted query. Gaia marks 0 matched
sources as `VARIABLE`. A Gaia value of `NOT_AVAILABLE` is not proof that a star
is constant; the ground pseudo-target curves remain the direct stability
check. The table also counts Gaia sources inside the 9.05-arcsecond ground
aperture so blends are visible instead of hidden.

## Nearby-star scope

44 deduplicated TIC sources within 60 arcseconds were bright enough to mimic a 2.91-ppt event in the simple total-eclipse screen and were measured on all 281 aligned images. Under the documented Eastern-time interpretation of the 2022 schedule, 32 sources are inconsistent with the required eclipse at this screen's three-sigma level. 2 source apertures overlap the target aperture and are not cleared.

The schedule-window result uses ingress BJD_TDB
2459782.682401 and egress
2459782.751151. The source row does not state
its time zone, epoch, uncertainty, or prediction source. Clearance here is a
conservative image-level screen, not the program's formal AstroImageJ NEB
procedure or planet validation. It cannot resolve or clear the known
0.517-arcsecond companion. The current ephemeris still places its nearest
event about 17.4 hours before this sequence.

## Files

- `01_ground_light_curve_checks.png`: ensemble, binning, and injection checks.
- `comparison_ensemble_light_curves.csv` and `comparison_ensemble_metrics.csv`.
- `noise_vs_bin_size.csv`.
- `ground_light_curve_injections.csv` and `ground_injection_summary.csv`.
- `comparison_star_catalog_matches.csv`.
- `nearby_star_catalog_candidates.csv`, plus image measurements when run.
- `summary.json`: key settings, counts, and limitations.
