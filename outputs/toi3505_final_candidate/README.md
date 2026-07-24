# TOI-3505.01 R-band light curve

## Result

This light curve is ready for mentor review. A recovered observing-schedule
row lists ingress at 00:15 and egress at 01:54. Interpreting those sheet clocks
as Eastern local time gives BJD_TDB
2459782.682401-
2459782.751151, fully inside the measured
sequence. A straight-baseline, fixed-window check gives
-0.658 +/-
0.395 ppt, so it does not find a
transit-like dimming in that exact historical window. Injecting a 2.91-ppt box
at the same times gives a total recovered depth of
2.252 ppt at
5.71 sigma.

The source row does not state its time zone, epoch, timing uncertainty, or
prediction source. The Eastern-time interpretation is used because the
sheet's planned 21:10-04:55 range brackets the actual image sequence. Under
the current catalog period and epoch, the nearest predicted midpoint is
BJD_TDB
2459781.873763, or
17.39 hours before the
first exposure. No physical transit model was fitted.

## Photometry settings

- 281 R-band exposures, 50 seconds each.
- 25-pixel source aperture.
- 70–139-pixel sky annulus.
- Ten comparison stars, C2 through C11.
- All ten comparison stars are unsaturated and have stable light curves.
- The 25-pixel aperture lowered the all-frame scatter from
  3.252 to 2.689 ppt compared with
  the original 35-pixel aperture, an improvement of
  17.3%.

## Frame review

All 281 measurements remain in the saved table. The plotted light curve uses
263 measurements. The other 18 frames have
documented image problems such as clouds, poor seeing, tracking, or trailing.
No frame was excluded only because the target appeared faint. Of the
23 measurements checked by eye, 5 had no
separate image problem and remain in the light curve.

The scatter of the plotted measurements is 2.497 ppt.
The 10-minute bins are included only to make the overall shape easier to see.

## Detrending

No correction was applied. Checks against airmass, sky level, star width, image position, and comparison-star counts did not improve the light curve consistently.

## Files to attach to the Discord post

The four files are collected in `../toi3505_post/`:

- `TOI_3505.01_2022-07-22_R_light_curve.png`
- `TOI_3505.01_2022-07-22_R_measurements.xls`
- `TOI_3505.01_2022-07-22_R_light_curve.plotcfg`
- `TOI_3505.01_2022-07-22_R_seeing_profile.png`

## Supporting files

- `TOI_3505.01_2022-07-22_R_light_curve.csv`: measurements and frame notes.
- `TOI_3505.01_2022-07-22_R_10min_bins.csv`: plotted 10-minute bins.
- `comparison_star_checks.csv`: comparison-star measurements.
- `detrending_checks.csv`: results of the trend checks.
- `frame_review.csv` and the two review figures: image-by-image notes.
- `analysis_settings.json`: settings used to make the light curve.
- `historical_schedule_check.json` and `historical_schedule_times.csv`: the
  preserved schedule interpretation and fixed-window result.
- `summary.json`: short numerical summary.

## Still needed

The mentor still needs to confirm the original spreadsheet or Transit Info
source, its time zone, prediction epoch and uncertainty, and review the
25-pixel aperture and comparison stars. The stored BJD_TDB
values agree with an independent calculation to within 0.000201 seconds, but
the observatory clock-sync record has not been found. Describe this as a light
curve submitted for review, not as a transit detection.
