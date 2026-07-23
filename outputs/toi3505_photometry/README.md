# TOI-3505.01 multi-aperture photometry

## What is finished

- AstroImageJ measured the target and ten comparison stars in all 281 usable
  aligned images.
- The source radius was 35 pixels in both runs.
- Two sky rings were tested: 70-100 pixels and 70-139 pixels.
- Neither table contains a saturated measurement.
- The native AstroImageJ tables, aperture file, screenshots, review tables,
  and plots are in this folder.

## Background-area check

The 70-100 pixel ring gave 0.330% robust scatter in the target
light curve. The 70-139 pixel ring gave 0.325%. The difference is
small, but the wider ring performed slightly better while AstroImageJ's
**Remove stars from background** setting was on. I kept the 70-139 pixels ring
for the next light-curve check because it also matches the Seeing Profile.

## Comparison-star check

All ten stars gave 0.325% robust target scatter. A working check
set using C2, C3, C5, C6, C7, C10, C11 gave 0.300%. The stars set
aside for review are C4, C8, C9. This is a working selection,
not a final claim; the individual curves and the review table show why each
star should be checked with the research mentor.

There are 25 image measurements marked for review because of an
unusual target value, large star width, low comparison-star counts, or a
saturation flag. They remain in the native tables. No detrending or transit
fit has been applied.

## Timing

The table covers BJD_TDB 2459782.598235 to 2459782.809459. Using the stated
period of 2.9151556 days and epoch 2459793.534385,
the closest predicted midpoint is
2459781.873763. It is
17.39 hours before
the first image, so this data set does not contain that predicted midpoint.
This light curve can still be used to check the reduction, apertures,
comparison stars, and systematic trends, but it should not be presented as a
transit detection.

## Saved AstroImageJ review files

The working reference-star calculation is saved in
`TOI_3505.01_2022-07-22_R_working_measurements.xls`. The matching aperture
file is `TOI_3505_working.apertures`, and the matching plot settings are in
`TOI_3505.01_2022-07-22_R_working.plotcfg`. The AstroImageJ screenshots show
the measurement table, reference-star selection, Y-data settings, raw light
curve, and fit settings. The transit fit and detrending are both off.
The saved table's target curve and total comparison counts match the seven
selected comparison stars.

The `.tbl` copy has the same contents as the `.xls` table and is included so
AstroImageJ 6.0.7 can reopen the working measurements directly.

## Next Schar tutorial step

Send the working light curve, seeing plot, measurement table, aperture view,
plot configuration, and observing-condition plots to the research mentor.
Ask for confirmation of the sky ring and comparison-star selection, and ask
whether this sequence is meant to be a reduction-quality exercise or whether
there is another TOI-3505 sequence that contains the predicted transit. Add a
predicted transit model only after that timing question is resolved.
