# Next step: Multi-Aperture photometry

The next step in the Schar light-curve tutorial is to make the measurement
table in AstroImageJ.

## Images

Use the 281 aligned images ending in `_aligned.fits`. Frames 0282 and 0283
are not usable and should stay out of the stack. Only the first image has a
plate solution, so use the aligned-image settings:

- Uncheck **Use RA/Dec to locate aperture positions**.
- Check **Halt processing on WCS or centroid error**.

## Starting aperture settings

- Source radius: 35 pixels
- Background inner radius: 70 pixels
- AstroImageJ starting outer radius: 139 pixels
- Clear outer-radius test: 100 pixels
- CCD readout noise: 1.414214
- CCD dark current per second: 0.012283

The 139-pixel outer ring crosses a nearby star, so compare
it with the 70-100 pixel background area
before choosing the final setting.

## Apertures

Click TOI-3505.01 first so it is T1. Then use the 10 stars in
`comparison_star_candidates.csv` as a starting list. They have similar sizes
and brightnesses to the target and are away from the image edges. They are not
automatically the final comparison set. After the table is made, plot each
comparison star's relative flux and remove stars with large scatter or clear
trends.

Save both the aperture file and the measurement table. The measurement table
will be used with `SCHAR_Plot_config.plotcfg` for the light curve, seeing,
airmass, total comparison counts, and target-position plots.
