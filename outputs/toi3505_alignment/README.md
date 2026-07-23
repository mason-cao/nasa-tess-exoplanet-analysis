# TOI-3505.01 image alignment

The first 281 calibrated images were matched to frame 0001. Frames
0282 and 0283 were left out because their star field was not usable in the
earlier image review.

The alignment used 25 reference stars. Every image was moved by
a whole number of pixels, so the original pixel values were not interpolated.
281 aligned FITS images were written to `data/ground/toi3505/aligned`.

- Typical star position error after alignment: 0.526 pixels
- Largest star position error after alignment: 0.819 pixels
- `frame_shifts.csv`: shift and position check for every image
- `reference_stars.csv`: stars used to measure the shifts
- `01_image_movement.png`: movement during the night
- `02_alignment_check.png`: first, middle, and last image around one star
- `summary.json`: main results and checks
- `verification.json`: file checks, observation times, and exact pixel checks

The next course step is to locate the target, save its Seeing Profile, and use
those radii for multi-aperture photometry.
