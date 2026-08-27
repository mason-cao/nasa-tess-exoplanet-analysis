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

409 deduplicated TIC sources within 150 arcseconds (2.5 arcminutes, the TFOP SG1 nominal radius) were bright enough to mimic a 2.91-ppt event once the TFOP -0.5 mag band correction is applied, and were measured on all 281 aligned images. Using the confirmed Eastern times from the 2022 schedule: 288 cleared, 4 blended with the target aperture, 115 too faint for a decisive limit, and 2 measured but without a limit deep enough to exclude the required eclipse. Of those 2, 0 show a dimming whose shape is consistent with an eclipse; the rest are monotonic trends across the night rather than events.

The schedule-window result uses ingress BJD_TDB
2459782.682401 and egress
2459782.751151. Mason confirmed that the
schedule clocks are Eastern time. The row does not provide the prediction
epoch, uncertainty, depth, or source. Clearance here is a conservative
image-level screen, not the program's formal AstroImageJ NEB procedure or
planet validation. It cannot resolve or clear the known
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
