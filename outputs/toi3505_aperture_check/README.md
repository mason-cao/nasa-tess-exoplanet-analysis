# TOI-3505.01 aperture-radius and time check

## Aperture result

The existing 35-pixel source radius fails at least one frozen precision-plateau rule and should not be treated as finalized.

Seven source radii (20, 25, 30, 35, 40, 45, and 50 pixels) were measured from
all 281 aligned images. The sky annulus remained 70-139 pixels. Photutils used
the AstroImageJ centroid positions but independently summed the images and
estimated each local sky with a three-sigma-clipped annulus median.

The 35-pixel Photutils and AstroImageJ differential curves have correlation
0.990350; their robust point-by-point difference is
0.540 ppt. Differences are expected because the two
programs do not use identical aperture-edge or background estimators.

The frozen retention rule does not simply pick the lowest target scatter. The
35-pixel radius must be within 10% of both (1) the best target differential
scatter and (2) the best median comparison pseudo-target scatter. This guards
against tuning the aperture only to make the target look flat.

## Time result

Recomputing BJD_TDB with Astropy from each stored `JD_UTC` midpoint, the GMU
site coordinates, and the J2000 target coordinate agrees with every stored
`BJD_TDB` to within 0.000201 seconds. Reconstructing
mid-exposure from the lower-precision `DATE-OBS` strings plus 25 seconds differs
by at most 2.042 seconds. The stored barycentric
conversion is therefore internally verified; observatory clock-sync provenance
still needs confirmation from the observer or mentor.

## Files

- `01_aperture_radius_check.png`: radius metrics, AstroImageJ comparison, and curves.
- `aperture_radius_metrics.csv`: one row per tested source radius.
- `python_multi_radius_light_curves.csv`: every extracted curve.
- `time_verification.csv`: all 281 independent timing calculations.
- `summary.json`: machine-readable verdict and statistics.

This is an independent implementation of the same observation, not an
independent astrophysical data set. It does not remove the need for mentor
review of the AstroImageJ settings.
