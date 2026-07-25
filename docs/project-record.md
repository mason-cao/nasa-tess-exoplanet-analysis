# TOI-3505.01 Project Record

Last updated: 2026-07-24

This is the compact project record for the TOI-3505.01 analysis. Numerical
results should be traced to the saved tables and JSON files in `outputs/`, not
reconstructed from this summary.

## Ground observation

- Telescope: GMU 0.8 m.
- Filter and exposure: Johnson-Cousins R, 50 seconds.
- Delivered science frames: 283.
- Working aligned sequence: frames 0001–0281. Frames 0282–0283 were retained
  but excluded because their stellar fields are unusable after a long gap.
- Final photometry: 25-pixel source aperture, 70–139-pixel sky annulus, and
  comparison stars C2–C11.
- Measurements retained in the table: 281.
- Measurements displayed after documented image-quality review: 263.
- Robust scatter of displayed points: 2.497 ppt.
- Detrending: none; no candidate correction passed the frozen blocked tests.
- Transit model: disabled.

The final plot, AstroImageJ settings screenshot, and nearby-star field are in
`outputs/toi3505_discord_post/`.

## Timing result

The supplied archives do not contain the target's original Transit Info file.
The recovered observing-schedule row lists ingress at 00:15, egress at 01:54,
and a planned range of 21:10–04:55, but does not state a time zone, epoch,
uncertainty, or prediction source.

Interpreting those clocks as America/New_York local time places the historical
window at BJD_TDB 2459782.682400599–2459782.751150993. This interpretation is
consistent with the night-time observing range but remains an assumption.

The fixed-window check measures −0.658 ± 0.395 ppt. An injected 2.91-ppt box at
the same times is recovered at 5.71 sigma, so this simple test would have seen
a signal of the expected scale. The result is a conditional null, not a failed
physical transit fit.

Using the current period of 2.9151556 days and epoch 2459793.534385 BJD_TDB,
the nearest predicted midpoint is 2459781.8737626, 17.39 hours before the first
usable measurement. The schedule period differs by only −0.588 seconds, which
does not explain the historical timing offset without the missing epoch.

## TESS result

The 2.915-day signal is recovered separately in QLP data from Sectors 14, 41,
54, and 81. First-pass box depths are:

| Sector | Depth (ppt) | Duration (hours) |
| ---: | ---: | ---: |
| 14 | 3.728 ± 0.419 | 1.80 |
| 41 | 3.015 ± 0.209 | 2.05 |
| 54 | 2.682 ± 0.295 | 2.05 |
| 81 | 3.130 ± 0.240 | 2.10 |

Official SPOC fits for Sectors 54 and 81 give depths of 3.226 ± 0.174 and
3.409 ± 0.173 ppt. The combined SPOC fit gives a period of 2.91514558 days and
a depth of 3.2919 ± 0.1185 ppt. These are comparisons using the same TESS
observations, not independent confirmation.

The ground observation occurs inside a Sector 54 data gap. Difference imaging
places the signal within one TESS pixel of the target system but cannot resolve
the known 0.517-arcsecond companion.

## Completed pipeline

1. Raw archive and FITS-header audit.
2. Dark and flat calibration with independent verification.
3. Whole-pixel alignment of 281 frames using 25 reference stars.
4. Plate solution and target-coordinate check.
5. Seeing-profile and multi-aperture review.
6. AstroImageJ and Python differential-photometry checks.
7. Frame-quality review, comparison-ensemble checks, and injection tests.
8. Four-sector TESS light-curve analysis.
9. TESS aperture, difference-image, centroid, and dilution screens.
10. Official SPOC report comparison and reproducibility manifest.

The ordered runner is `src/run_toi3505_analysis.py`.

## Key evidence

- Final ground products: `outputs/toi3505_final_candidate/`
- Aperture selection: `outputs/toi3505_aperture_check/`
- Ground robustness and nearby stars: `outputs/toi3505_ground_checks/`
- TESS sector measurements: `outputs/toi3505_tess_analysis/`
- TESS pixel checks: `outputs/toi3505_tess_pixels/`
- SPOC comparison: `outputs/toi3505_data_validation/`
- Reproducibility record: `outputs/toi3505_research_record/`

## Remaining work

1. Post the three-image package and record the mentor response.
2. Ask the mentor to confirm the original timing source, time zone, epoch, and
   uncertainty.
3. Request review of the 25-pixel aperture and C2–C11 comparison ensemble.
4. Run the program's formal AstroImageJ NEB procedure if required.
5. Keep claims limited to the unresolved target system unless additional
   resolving evidence becomes available.

## Non-negotiable limits

- Do not report a ground transit detection from this sequence.
- Do not use the generic values in the supplied Schar plot configuration as
  TOI-3505 parameters.
- Do not describe same-observation pipeline comparisons as independent
  confirmation.
- Do not claim that seeing-limited ground data or TESS resolves the close
  companion.
- Keep all excluded measurements and their reasons available for audit.
