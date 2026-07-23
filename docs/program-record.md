# Schar Program Record and Progress Tracker

Last updated: 2026-07-23

This file preserves the program materials, Discord instructions, Mason's reported posts, and the evidence-based next steps for this repository. Values or requirements should not be invented. Use the supplied materials and dated announcements; ask a mentor when a required target-specific input is absent.

## Archived materials

The following three files were downloaded on 2026-07-23 and copied into `data_and_lectures/` without modifying their contents:

- `TFOP_SG1_Guidelines_Latest.pdf`
  - TFOP SG1 Observation Guidelines, revision 6.4, dated 2020-09-04.
  - 32 pages.
  - SHA-256: `24031fe16896259abc83f4f6b53c41c0fb40b329e44bfc647218a83f0a57abc9`
- `Campus Telescope TESS Follow-Up Light Curve Tutorial - Schar Program.docx.pdf`
  - Schar Program AstroImageJ reduction and light-curve tutorial.
  - 30 pages.
  - SHA-256: `7c536b9ea8a97045786dff5d31124afc254da49019416ac002c017da8a5d4a13`
- `SCHAR_Plot_config.plotcfg`
  - AstroImageJ plot configuration supplied with the 2026 instructions.
  - Internal timestamp: 2026-07-15 19:32:28 EDT.
  - SHA-256: `81d657abbdd425be2194aaa4881c84e93fc8f9faaa05d08dd32927f4fedbe546`

Important interpretation:

- The older tutorial refers to `Template.plotcfg` and an older EXOFAST limb-darkening page.
- The current Discord instructions instead say to use the supplied Schar plot configuration and the SCO LDC Calculator. Follow the current instructions for internship work.
- `SCHAR_Plot_config.plotcfg` is a reusable starting configuration, not a TOI-3505 parameter file. It contains a placeholder title (`Exoplanet Name, UTyyyy-mm-dd`), placeholder subtitle, period `1.0`, and generic marker values. Preserve the original file and replace target-specific fields only inside the working AIJ session or a target-specific saved copy.

## Discord announcements supplied by Mason

### Limb darkening calculator

Author and timestamp were not included in the copied text:

> @everyone Use this website to calculate limb darkening coefficients. The filter we used is a Johnson-Cousins (R)
>
> https://sco-ldc.com/
>
> SCO LDC Calculator
>
> Computes quadratic limb-darkening coefficients (u1, u2) for exoplanet transit light-curve fits.

Operational requirement: select Johnson-Cousins R for this dataset. Do not substitute a TESS or Sloan filter.

### Plavchan - 2026-07-17 14:20

> @everyone - as you are making light curves and trying to understand them, this is the hard part. Please post the lcs, and info like the seeing plot and xls file and plotcfg. We need to ensure you have a good aperture size and to use our plotcfg to see the systematic trends as well. Also please watch the three videos I sent out. At the end of the last one I teach a rudimentary to evaluate the statistical significance of the detection of a transit. We also need to check the timing uncertainty. When you are ready to add a transit model, please use the predicted time, depth and duration first before trying to fit the transit parameters to get a sense for what the light curve should look like. I’m seeing a lot of light curves that have a lot of systematics which tends to imply poor detrending or poor apertures or reference stars. One does need to optimize the reference star selection and detrending. I will try to record another short video on this.

### starulae - 2026-07-19 16:39

> As a reminder, please make sure to upload screenshots of all windows when asking about plotting and making your light curves. Thank you!

## Discord channel map

The three main help channels reported by Mason are:

- `data reduction questions`: calibration, darks/flats, raw science images, aperture-photometry setup, or measurement-table generation.
- `plate solving questions`: WCS, Astrometry.net, API-key, coordinate, and plate-solving problems.
- `light curve questions`: plotting, predicted markers, transit models, timing uncertainty, detrending, reference-star selection, and statistical significance.

TOI-3505 timing/model questions belong in `light curve questions`.

## Required light-curve troubleshooting package

Based on the July 17 and July 19 announcements, a plotting or light-curve question should include:

- the light-curve plot;
- the seeing-profile plot;
- the AstroImageJ measurement table (`.xls`, `.tbl`, or equivalent);
- the target-specific saved `.plotcfg`, derived from the supplied Schar configuration;
- screenshots of all relevant AstroImageJ windows, not just the plot;
- the chosen source aperture and sky-annulus radii;
- enough information to evaluate comparison/reference-star selection and detrending.

The Schar tutorial and TFOP guide add these workflow requirements:

- use `BJD_TDB` as the X data;
- initially use predicted ingress and egress as fixed V.Markers;
- begin with the target's raw, unfitted, undetrended light curve at bin size 1 and shift 0;
- plot systematics including `Width_T1`, `Sky/Pixel_T1`, inverted `AIRMASS` with `tot_C_cnts`, `X(FITS)_T1`, and `Y(FITS)_T1`;
- inspect every comparison star and remove high-scatter choices;
- only add a transit model after entering the predicted time, depth, and duration;
- calculate quadratic limb-darkening coefficients with Johnson-Cousins R when the model step is reached;
- check timing uncertainty and evaluate detection significance as demonstrated in the three videos;
- save the measurement table, aperture file, light curve, field with apertures, seeing profile, plate-solved image, `.plotcfg`, and any requested NEB/dmag products.

The three videos mentioned in the July 17 announcement are not currently stored in this repository. Their files or links still need to be obtained and reviewed.

## Mason's Discord posting record

Mason reports that the following three posts have been sent.

### Post 1 - TOI-3505.01

> i started checking the TOI-3505.01 data in AstroImageJ and found 283 50-second R-band images with the matching darks and flats. they cover BJD_TDB 2459782.598–2459782.839, but using P = 2.9151556 days and T0 = 2459793.534385 put the closest transit outside that range. did i calculate it wrong, or does this dataset not actually cover a transit?

Four AstroImageJ-related screenshots were attached. No mentor response had been received when this record was updated.

### Post 2 - AU Mic b

> i used the tess2018206045859-s0001-0000000441420236-0120-s_lc.fits dataset for AU Mic b, Sector 1. for the full light curve, i kept the whole time range and used ylim = 0.98–1.06; for the zoomed-in graph, i used t0 = 1330.402640 with xlim = t0 - 0.25 to t0 + 0.25 and ylim = 0.989–1.003.

### Post 3 - TOI-3505.01 follow-up

Posted in `light curve questions` on 2026-07-23:

> Following up on my TOI-3505.01 homework post: I don't see a Transit Info PNG with the TOI-3505.01 files. Which predicted ingress/egress times, depth, duration, and timing uncertainty should I use for the July 21, 2022 light curve?

The earlier alternative follow-up drafts were not posted.

## TOI-3505 evidence and current blocker

The TOI-3505 files are split across six `TOI_3505.01-20260714T190730Z-1-00*.zip` archives. Inspection found science FITS files, darks, flats, focus FITS files, and focus images, but no `Transit Info` PNG, `.radec` file, or other target-specific timing document.

The existing audit found:

- 283 science frames at 50 seconds in R;
- 10 matching 50-second darks;
- 10 matching 3.5-second flat darks;
- 10 matching 3.5-second R flats;
- image coverage `BJD_TDB 2459782.598234811-2459782.838668314`;
- the current NASA Exoplanet Archive ephemeris used in the original post places the adjacent predicted transit centers outside that image range.

The Schar tutorial says predicted ingress/egress and period should come from the target's `Transit Info` PNG. The July 17 announcement also requires the predicted time, depth, duration, and timing uncertainty before fitting. Because that target-specific source is absent, the supplied Schar plot configuration's generic values must not be used as TOI-3505 values.

## TOI-3505 data reduction completed - 2026-07-23

The source-independent calibration stage is complete. `src/reduce_toi3505.py` implements the Schar AstroImageJ tutorial settings as a reproducible pipeline and writes AstroImageJ-compatible 32-bit FITS products. It performed:

- a median combination of the 10 matching 3.5-second darks;
- a median combination of the 10 matching 50-second darks;
- subtraction of the 3.5-second master dark from each R flat;
- fitted illumination-plane removal from each calibrated flat, followed by normalization and median combination;
- science calibration using `(raw science - mdark_50.000s) / mflat_R`;
- preservation of the original timing metadata, including `BJD_TDB`;
- FITS provenance, `CHECKSUM`, and `DATASUM` cards for each output.

No separate bias correction was applied because the exposure-matched dark contains the bias signal. No dark scaling, cosmic-ray filtering, plate solving, alignment, photometry, detrending, or transit model was applied. Those remain separate reviewable stages.

Local products and evidence:

- `data/ground/toi3505/calibration/`: two master darks, the R master flat, and a diagnostic flat-review mask;
- `data/ground/toi3505/reduced/`: 283 calibrated science FITS files (about 18.99 GB);
- `outputs/toi3505_reduction/frame_manifest.csv`: per-frame timing, header, and sampled quality statistics;
- `outputs/toi3505_reduction/01_master_calibrations.png` through the seven `04_calibrated_contact_sheet_page_*.png` files: calibration diagnostics and a visual-review image for every frame;
- `outputs/toi3505_reduction/05_astroimagej_frame_0001.png`, `06_astroimagej_frame_0282.png`, and `07_astroimagej_frame_0283.png`: direct AstroImageJ checks confirming the 4096 x 4096, 32-bit outputs open correctly and showing both end-frame quality failures;
- `outputs/toi3505_reduction/verification.json`: independent integrity and formula verification.

Verification passed for all 283 science products and all four calibration products. All FITS checksums and data sums are valid, all sampled science images are finite, all output images are 4096 x 4096 with `BITPIX=-32`, the BJD sequence is strictly increasing, and independent pixel checks on frames 0001, 0142, and 0283 exactly reproduce the calibration equation while preserving `BJD_TDB`.

The automated review marked frames 0270, 0282, and 0283. Frame 0270 has an elevated background but its stellar field remains visible. Frames 0282-0283 occur after an approximately 41-minute gap, have a much brighter background, and direct AstroImageJ inspection confirms that neither contains a usable stellar field. Preserve them for traceability, but exclude them from the plate-solving and photometry working sequence. The master flat's low-response first detector row and one isolated low-response pixel are recorded in `master_flat_review_mask.fits`; they were not silently clipped or replaced.

## Next best course of action

1. Use frames 0001-0281 as the current working sequence and plate-solve or align them in AstroImageJ. Keep 0282-0283 in the archive but outside the working sequence. Save the plate-solved image and any target coordinate file used.

2. Create and save the seeing profile, choose an initial aperture and sky annulus from it, perform multi-aperture photometry, and save the measurement table and aperture file. Do not permanently reject other variable-background frames until their seeing, source counts, and differential-photometry behavior have been checked.

3. The request for the missing timing/model inputs is awaiting a response in the existing TOI-3505 thread in `light curve questions`. Load `SCHAR_Plot_config.plotcfg`, replace its template title/subtitle in the working session, and first produce the required raw, unfitted, undetrended light curve plus the systematics panels. Do not enter target markers or fit a transit until the target-specific predicted inputs are confirmed.

4. Test aperture sizes and comparison-star selections, preserving each meaningful working result. Once a complete package exists, post the light curve, seeing profile, measurement table, target-specific `.plotcfg`, and screenshots of all relevant windows in `light curve questions`.

5. Obtain and watch the three videos before making a significance claim or presenting a fitted transit result.
