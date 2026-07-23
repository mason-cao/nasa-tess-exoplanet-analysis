# TOI-3505.01 calibrated image set

Made with `src/reduce_toi3505.py` (version 1.1.0).

## Calibration equation

`reduced = (raw science - mdark_50.000s) / mflat_R`

- Master darks: median of 10 matching frames at each exposure time.
- Master flat: each 3.5 s R flat was dark-subtracted, its fitted illumination plane was divided out, it was normalized to mean 1, and the 10 normalized flats were median-combined.
- Output science images: 32-bit floating-point FITS with original timing metadata, calibration HISTORY cards, and FITS checksums.
- Bias subtraction: not applied because the matching dark frames contain the bias signal.
- Not applied: dark exposure scaling, nonlinearity correction, outlier/cosmic-ray filtering, plate solving, alignment, or photometry.

## Contents

- `calibration/mdark_3.500s.fits`
- `calibration/mdark_50.000s.fits`
- `calibration/mflat_R.fits`
- `calibration/master_flat_review_mask.fits` (diagnostic only; not applied)
- `reduced/*_out.fits`
- `outputs/toi3505_reduction/frame_manifest.csv`
- `outputs/toi3505_reduction/calibration_inputs.csv`
- `outputs/toi3505_reduction/summary.json`
- `outputs/toi3505_reduction/reduction.log`
- `outputs/toi3505_reduction/04_calibrated_contact_sheet_page_*.png` (every reduced frame)
- `outputs/toi3505_reduction/05_astroimagej_frame_0001.png`, `outputs/toi3505_reduction/06_astroimagej_frame_0282.png`, and `outputs/toi3505_reduction/07_astroimagej_frame_0283.png` (AstroImageJ compatibility/end-frame checks)
- `outputs/toi3505_reduction/verification.json` (independent checksum, timing, and formula validation)
- `outputs/toi3505_alignment/` (movement plot, alignment check, shift table, and file checks)

## Status

- Reduced frames: 283
- Frames marked for a closer look: 3
- Frames with WCS headers before plate solving: 0
- Calibrated frames represented in visual-review sheets: 283

The marked frames are not automatically rejected. Frames 0001-0281 were aligned separately with whole-pixel shifts. Frames 0282-0283 were kept in the archive but left out of the working sequence because their star field is not usable. The next step is to locate the target, save its Seeing Profile, and use those radii for multi-aperture photometry.

## Run it again

From the repository root:

```bash
.venv/bin/python src/reduce_toi3505.py
.venv/bin/python src/verify_toi3505_reduction.py
```
