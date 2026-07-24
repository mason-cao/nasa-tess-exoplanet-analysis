# AstroImageJ workflow practice: TOI-3718.01

This is a practice reduction of the public GMU `2023-10-12` folder. It must not
be presented as TOI-6241.01: all 270 science filenames and FITS `OBJECT` headers
identify the field as `TOI 3718.01`.

## Reduction summary

- 270 R-band science frames at 65 seconds each
- 15 dark frames and 8 R-band flat frames
- Median-combined darks; dark-corrected, normalized flats
- Seeing-profile-guided 16-pixel source aperture
- 22-32-pixel local-sky annulus
- 12 retained comparison stars
- 181 frames passed the automated quality cuts
- Approximate first-frame FWHM: 15.5 pixels
- Post-window robust scatter: 2.64 ppt

The middle of the observation is strongly clouded. A residual offset remains
between the usable pre-cloud and post-cloud segments, so the apparent depth is
not scientifically interpretable. This is an inconclusive practice light curve,
not a claimed transit detection.

## Suggested Discord attachments

1. `01_toi3718_field_identification.png`
2. `02_toi3718_seeing_profile.png`
3. `03_toi3718_practice_light_curve.png`

`04_toi3718_reduction_diagnostics.png` is useful if the mentor asks to see the
clouds, pointing changes, comparison-star stability, or rejected frames.

