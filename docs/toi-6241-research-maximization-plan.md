# TOI-6241.01 maximum-research blueprint

| Field | Value |
|---|---|
| Working project | Ground-based follow-up and advanced vetting of TOI-6241.01 |
| Target | TOI-6241.01 / TIC 15654898 |
| Document status | Living research design, not a statement that every analysis has already been completed |
| Last updated | 2026-07-13 |
| Program fit | NASA/TESS exoplanet-analysis internship research paper |

## Purpose

This document is the complete research blueprint for making the TOI-6241.01 project as rigorous, original, and scientifically useful as the available data allow while staying inside the internship rules.

The required paper remains a ground-based TESS Object of Interest follow-up study: determine whether the GMU observation contains a transit-like event on the intended target at the predicted time, depth, and duration. Every advanced component below strengthens that central result; none replaces it.

The project can still produce a defensible paper if the result is a nondetection, partial event, neighboring-star event, or unusable night. The outcome must be reported honestly and the title, abstract, and claims must match the evidence.

## Non-negotiable data-identity warning

The spreadsheet assigns the 2023-10-12 observation to **TOI-6241.01**, but the science FITS files currently found in the corresponding public data location identify **TOI-3718.01**. Those files must not be analyzed or described as TOI-6241.01 data.

The TOI-6241 study is blocked at the ground-based reduction stage until one of these occurs:

1. The mentor provides or identifies the correct TOI-6241.01 GMU files.
2. The mentor confirms that the headers are wrong but independently verifies that the field is TOI-6241.01 through coordinates, plate solving, logs, and finder-chart agreement.
3. The mentor approves a formal target change and the research question is rewritten for the replacement target.

This mismatch is not a minor metadata issue. Target identity is the first scientific validation test.

## Scope legend

Each proposed analysis is assigned one of five levels.

| Level | Meaning | Rule |
|---|---|---|
| **CORE** | Needed for the internship paper | Complete before pursuing advanced claims. |
| **HIGH-VALUE BONUS** | Feasible with public or existing data and likely to improve the paper | Prioritize after the first trustworthy light curve. |
| **ADVANCED** | Statistically or technically demanding | Use only if the method, assumptions, and code can be explained and validated. |
| **CONDITIONAL** | Useful only for certain outcomes or if additional data exist | Activate only when its dependency is satisfied. |
| **GATED / FUTURE** | Requires mentor approval, new observations, restricted products, or a separate paper | Describe as future work unless approval and data are obtained. |

## One-sentence research strategy

Test whether GMU R-band photometry independently recovers the predicted TOI-6241.01 transit, determine whether any signal originates on the target rather than a nearby star, quantify how robust the result is to reasonable analysis choices, and place the observation in the timing and false-positive context supplied by TESS, Gaia, and public follow-up constraints.

## Strongest defensible novelty thesis

The best potential contribution is not simply another plotted transit. It is a **source-localized, robustness-audited, ground-versus-space comparison of a still-unconfirmed TESS candidate**, with a quantitative explanation of what the GMU night can and cannot rule out.

The paper should aim to contribute four evidence layers:

1. **Recovery evidence:** whether an R-band event is present at the predicted time, depth, and duration.
2. **Localization evidence:** whether the target, a cataloged neighbor, or unresolved contamination is the most likely source.
3. **Robustness evidence:** whether the conclusion survives aperture, comparison-star, detrending, quality-cut, and pipeline choices.
4. **Timing and validation context:** whether the event timing improves the ephemeris and how the new observation changes candidate-versus-false-positive interpretations.

No document or abstract should call the work “first,” “novel,” “confirmed,” or “validated” until the corresponding checks in this blueprint are completed.

## Program-scope compliance

### What the lectures require

The lecture materials frame the paper around GMU campus-telescope follow-up of a TESS candidate and explicitly allow inconclusive results. The central test is whether a transit appears on the correct target near the expected time with a compatible depth and duration. See [Lecture 8, paper-project overview](../data_and_lectures/Lecture8_Schar2025_clean_redesign.pdf#page=13), [Lecture 8, result possibilities](../data_and_lectures/Lecture8_Schar2025_clean_redesign.pdf#page=14), and [the paper-structure guide](../data_and_lectures/Paperstructure.pdf#page=5).

### What the lectures permit as advanced work

The paper-structure slides describe statistical false-positive work such as `vespa` as beyond the normal internship scope but welcome as a bonus. That makes advanced validation permissible, not mandatory. It also means a validation result must not be attempted at the expense of a correct reduction and basic transit analysis. See [the advanced-analysis note](../data_and_lectures/Paperstructure.pdf#page=16).

### What is gated

The lecture sequence treats JWST work as a later opportunity requiring approval after a campus-telescope light curve is approved, and as a separate paper. Atmosphere/JWST work therefore appears only in the gated extension section. See [Lecture 8, JWST gate](../data_and_lectures/Lecture8_Schar2025_clean_redesign.pdf#page=4).

### AI and authorship rule

Any code, equation, model, or wording used in the final project must be understood and defensible by the student. Preserve an analysis log and disclose assistance as required by the program. See [Lecture 1, AI guidance](../data_and_lectures/Lecture1_intro.pdf#page=26).

## Target fact sheet

Values below come from downloaded NASA Exoplanet Archive and ExoFOP-style TOI snapshots. The archive export used here was released 2026-07-08; its TOI row was last updated 2024-08-22. The ExoFOP-style export was modified 2025-07-31. They must be re-snapshotted immediately before the final paper because candidate parameters and dispositions can change.

| Quantity | Snapshot value | Interpretation or caution |
|---|---:|---|
| TIC ID | 15654898 | Stable target identifier to use in every search. |
| TOI | 6241.01 | Candidate identifier. |
| TESS disposition | PC | Planet candidate, not confirmed or validated by this label. |
| TFOP disposition | PC | Follow-up disposition in the downloaded snapshot. |
| Right ascension | 01:08:19.84 / 17.082679° | Must match plate-solved GMU field. |
| Declination | +25:39:09.27 / +25.652576° | Must match plate-solved GMU field. |
| TESS magnitude | 12.0227 | More relevant to TESS photon counts than the spreadsheet V magnitude. |
| Spreadsheet V magnitude | 12.75 | Intended GMU observation metadata. |
| Period | 3.7812464 ± 0.0000165 d | Starting ephemeris, not a fixed truth. |
| Reference epoch | BJD 2459879.12449 ± 0.0025425 d | Confirm time system is BJD_TDB in the source. |
| Duration | 2.130 ± 0.332 h | Expected full-event scale. |
| Depth | 2521 ± 195 ppm | 2.521 ppt, 0.2521%, or approximately 2.741 mmag. |
| Radius ratio implied by depth | about 0.0502 | Approximation: `sqrt(0.002521)` before dilution correction. |
| Candidate radius | 4.975 ± 0.374 Earth radii | Catalog/model value, dependent on stellar radius and transit fit. |
| Host effective temperature | 5209 ± 131 K | Approximate late-G/early-K temperature range; avoid assigning a spectral subtype without a cited classification. |
| Host radius | 0.990 ± 0.060 solar radii | Input to physical-radius inference. |
| Host mass | 0.886 ± 0.109 solar masses | Catalog/model estimate. |
| Host surface gravity | log g = 4.40 ± 0.09 | Consistent with a dwarf-like host in the snapshot. |
| Distance | 278.827 ± 3.087 pc | Gaia-derived/catalog distance in the snapshot. |
| Insolation | 167.006 Earth fluxes | Model-derived candidate context. |
| Equilibrium temperature | 1001 K | Model-dependent; assumptions must be stated if used. |
| TESS sectors | 17 and 57 | Analyze separately before combining. |
| Alert date | 2023-04-03 | Candidate-history context. |
| Last TOI update in snapshot | 2024-08-22 | Recheck current archive. |
| Public follow-up counts in downloaded snapshot | 0 time series, 3 spectroscopy, 2 imaging | Counts are not scientific constraints. Inspect permissions, products, dates, and quality before use. |
| Predicted mass | 21.86 Earth masses | Prediction, not a measurement. Do not present as measured mass. |
| Predicted RV semi-amplitude | 9.7 m/s | Prediction, useful only for future-observation motivation. |
| TSM / ESM | 44 / 8.3 | Prioritization metrics from the snapshot, not atmosphere detections. |

Primary catalog starting point: [NASA Exoplanet Archive TOI table](https://exoplanetarchive.ipac.caltech.edu/cgi-bin/TblView/nph-tblView?app=ExoTbls&config=TOI).

## Intended GMU observation

| Field | Spreadsheet entry |
|---|---|
| Date | 2023-10-12, Thursday |
| Target | TOI-6241.01 |
| Predicted ingress / egress labels | 01:31 / 03:39 |
| Planned observing window | 19:35 to 06:20 |
| Filter | R |
| Exposure | 75 s |
| Priority | 3 |
| Disposition at scheduling | PC |
| Observers | Ian and Tommy |
| Listed period | 3.7812464 d |

At 75 s per exposure, a 2.13 h event contains roughly 102 exposure intervals before accounting for readout and rejected frames. The full 10 h 45 min planned window could contain roughly 516 exposure intervals before overhead. These estimates are planning checks, not the observed frame count.

## Data-identity protocol

Before calibration, create a manifest with one row per file and the following fields:

- Original path and filename.
- SHA-256 checksum.
- FITS object/target name.
- DATE-OBS, exposure time, filter, image type, binning, detector temperature, and observatory coordinates when present.
- Header right ascension and declination.
- Plate-solved field center and angular separation from TOI-6241.01.
- Science, flat, dark, bias, test, focus, or unknown classification.
- Keep/reject decision and a written reason.
- Any filename-versus-header disagreement.

Identity acceptance requires all available evidence to agree: observing log, object headers, plate solution, finder chart, target position in the field, date, filter, and exposure sequence. A filename alone is insufficient. A spreadsheet row alone is insufficient.

## Novelty verification protocol

Preliminary exact-name searches performed through 2026-07-13 for “TOI-6241,” “TOI 6241,” and “TIC 15654898” did not reveal an obvious dedicated paper in the searched title/abstract results. A downloaded follow-up snapshot also reported zero time-series observations. Those are promising leads, not proof of novelty or proof that no private/unpublished time series exists.

Before choosing any novelty language:

1. Search NASA ADS by TOI number, TIC ID, coordinates, and host aliases.
2. Search arXiv using the same identifiers and spelling variants.
3. Check the current NASA Exoplanet Archive disposition and reference fields.
4. Check current ExoFOP-TESS time-series, imaging, spectroscopy, stellar, and observation-note tabs.
5. Ask the mentor whether unpublished TFOP products, student projects, or papers in preparation already use the observation.
6. Check the GMU internship examples and repository history for a prior TOI-6241 analysis.
7. Record every query, exact search string, date, result URL, and screenshot or export in a novelty ledger.
8. Repeat the audit during final drafting and once more before submission.

Safe wording before this audit is: “We analyze a GMU follow-up observation of the TESS planet candidate TOI-6241.01.”

Potential wording after a clean audit is: “To our knowledge, this is the first publicly documented analysis of this GMU TOI-6241.01 observation.” The phrase still needs mentor approval.

## Research questions

### Primary question — CORE

Does the correctly identified GMU R-band time series show a transit-like decrease on TOI-6241.01 near the predicted window, with a depth and duration statistically compatible with the TESS candidate?

### Source-localization question — CORE / HIGH-VALUE BONUS

Is the event centered on TOI-6241.01, absent from nearby stars, and robust to aperture contamination?

### Robustness question — HIGH-VALUE BONUS

Does the scientific conclusion remain stable across defensible choices of aperture, comparison ensemble, background annulus, quality cuts, and detrending model?

### TESS consistency question — HIGH-VALUE BONUS

Do Sector 17 and Sector 57 independently recover a compatible period, depth, duration, and transit shape in SAP, PDCSAP, and a transparent custom extraction?

### Ephemeris question — CONDITIONAL

If the GMU event is well localized and precisely timed, does its mid-transit time improve or challenge the linear ephemeris when combined with individual TESS transit times?

### Detectability question — HIGH-VALUE BONUS

If the event is not recovered, was a 2.521 ppt, 2.13 h transit actually detectable in that night’s noise and coverage?

### False-positive question — ADVANCED

How do the GMU localization result, Gaia neighborhood, TESS centroid/difference-image tests, archival imaging/spectroscopy constraints, and statistical validation models change the probability of planet, eclipsing-binary, background-eclipsing-binary, and hierarchical-eclipsing-binary scenarios?

### Physical-context question — CONDITIONAL

Where would a roughly 4.98 Earth-radius, 3.78 d candidate around a roughly solar-radius, 5209 K star fall relative to short-period Neptune/sub-Saturn populations, and what follow-up measurements would be most informative?

## Hypotheses and outcome-neutral tests

The hypotheses must be written before examining the predicted GMU transit window in the final reduction.

- **H1:** The target light curve contains a negative flux event near the propagated transit window.
- **H2:** The event depth and duration are compatible with the TESS-derived values after uncertainty and bandpass differences are included.
- **H3:** Neighbor light curves do not show a deeper event capable of producing the target-aperture signal through contamination.
- **H4:** The sign and approximate size of the event survive the predeclared robustness suite.
- **H5:** A nondetection, if obtained, is informative only if injection–recovery shows high completeness for the predicted event.
- **H6:** A single timing offset is not evidence of transit-timing variations; a TTV interpretation requires multiple precise timings and a model comparison.

### Recommended primary endpoint and classification rule

The primary numerical endpoint should be the R-band transit depth from the frozen prediction-conditioned transit-plus-baseline model, with mid-time and duration allowed to vary under declared TESS-informed priors. The mentor should approve exact thresholds before the final GMU transit window is unblinded.

A strong default framework is:

- **Robust recovery:** the depth posterior excludes zero at the predeclared credibility level; the transit model improves blocked out-of-transit predictive performance or the declared model score; fitted time and duration are compatible with the propagated prediction; no control time or pseudo-target produces an equally strong event at an alarming rate; and the result remains directionally stable across the acceptable robustness set.
- **Tentative recovery:** the expected event is favored, but one major test such as red-noise significance, baseline coverage, localization, or reduction stability remains weak.
- **Informative nondetection:** injection–recovery completeness is at least the predeclared high threshold, preferably 90% or greater, for the predicted depth/duration across the allowed timing window, yet no compatible event is recovered.
- **Inconclusive:** the data fail both the robust-recovery and informative-nondetection requirements.

Do not combine many weak, correlated diagnostics as though they were independent sigma evidence. Report each component and the classification logic.

## Claim ladder

Every statement in the paper should stop at the highest rung supported by the completed evidence.

| Rung | Permitted claim | Minimum evidence |
|---|---|---|
| 0 | No interpretable constraint | Wrong target, insufficient baseline, severe systematics, or inadequate precision. |
| 1 | A dip is visible | Descriptive light curve only; not yet a transit claim. |
| 2 | A transit-like event is detected | Model-supported dip with defensible uncertainty and null-test performance. |
| 3 | Event is consistent with TOI-6241.01 prediction | Compatible time, duration, and depth. |
| 4 | Event is localized to the target at the achieved angular resolution | Neighbor-star photometry and/or difference-image/centroid evidence excludes cataloged contaminants to stated limits. |
| 5 | Candidate receives additional validation evidence | Multiple independent vetting tests and explicit remaining scenarios. |
| 6 | Statistically validated planet | Full validated-planet criteria, all required contrast/spectroscopic constraints, vetted model inputs, mentor approval, and no contradictory evidence. |
| 7 | Dynamically confirmed planet | Mass or dynamical evidence meeting accepted confirmation standards; outside the present data unless new observations exist. |

One GMU light curve by itself normally supports rungs 2–4, not rungs 6–7.

## Analysis governance and anti-bias plan

### Freeze before unblinding

Before inspecting the final predicted window, save a dated configuration containing:

- Adopted coordinates and ephemeris source.
- Expected timing uncertainty calculation.
- Calibration rules.
- Initial aperture grid and background-annulus grid.
- Comparison-star eligibility rules.
- Quality metrics and rejection thresholds.
- Primary detrending covariates.
- Primary transit model and priors.
- Detection and nondetection criteria.
- Injection–recovery criteria.
- Primary figures and tables.

If practical, optimize comparison stars and systematics on out-of-transit data while masking the predicted window. Any change after unblinding belongs in the analysis ledger with its reason and effect.

### Separate confirmatory and exploratory results

- **Confirmatory:** the frozen primary reduction and model.
- **Robustness:** predeclared alternative reasonable reductions.
- **Exploratory:** ideas added after seeing the signal.

Exploratory findings may motivate future work but must not be presented as preplanned independent evidence.

### Preserve all forks

Do not keep only the aperture, comparison ensemble, or trend model that produces the most planet-like dip. Save the full analysis grid, including null and contradictory outcomes.

## Complete data inventory

### Required ground-based inputs — CORE

- Correct TOI-6241.01 GMU science FITS sequence.
- Same-night R-band flats whenever available and valid.
- Exposure-matched dark frames or a documented scaled-dark method supported by detector behavior.
- Bias frames if required by the camera/calibration procedure.
- Observer log, weather notes, focus/pointing events, and clock information.
- Telescope, detector, gain, read-noise, pixel-scale, binning, and filter metadata.
- Finder chart and aperture file if they exist.

### TESS inputs — HIGH-VALUE BONUS

- Sector 17 and Sector 57 light-curve FITS products.
- Target Pixel Files when available.
- Full-frame-image cutouts if no suitable target product exists or to cross-check extraction.
- Data Validation report and mini-report when available.
- TESS Input Catalog and contamination-related fields such as `CROWDSAP` and `FLFRCSAP`.
- Cadence-level quality flags, centroids, background, and motion/systematics vectors.

### Catalog and follow-up inputs — HIGH-VALUE BONUS

- Current NASA Exoplanet Archive TOI snapshot.
- Current ExoFOP-TESS public metadata and products.
- Gaia DR3 cone search and astrometric/photometric columns.
- TIC neighbor table with TESS magnitudes.
- Public high-resolution imaging contrast curves.
- Public reconnaissance spectra and stellar parameters.
- Public ground-based time-series photometry, if any becomes available.

### Optional external inputs — CONDITIONAL

- Pan-STARRS, DSS, 2MASS, or similar archival images for proper-motion/background checks.
- ASAS-SN and ZTF time series for long-term variability or rotation context.
- Public SED photometry for an independent stellar-radius/temperature check.
- Public radial velocities, line-bisector measurements, or spectroscopic binary flags.
- Additional TESS sectors if the target is reobserved after the snapshot date.

External surveys must be used within their licenses and cited. Restricted TFOP products require permission.

## Ground-based reduction — CORE

### 1. Raw-file audit

1. Build the immutable manifest.
2. Inspect header consistency and exposure chronology.
3. Plate solve representative early, middle, and late frames.
4. Confirm the target remains in the field and quantify drift.
5. Plot median background, FWHM, ellipticity, airmass, centroid position, and total counts versus time.
6. Mark clouds, focus changes, meridian flip, tracking jumps, saturation, cosmic rays, and edge proximity without looking at the transit shape.

### 2. Calibration audit

1. Inspect every calibration frame and reject only with a recorded reason.
2. Create master dark, flat, and bias products with robust combination.
3. Normalize flats and inspect gradients, dust features, and temporal stability.
4. Confirm dark exposure and detector temperature compatibility.
5. Compare calibrated and uncalibrated statistics for a small frame sample.
6. Save master calibration products, uncertainty maps if supported, and processing logs.
7. Test whether flat-fielding or dark subtraction introduces structured artifacts.

### 3. Time standard

1. Preserve original DATE-OBS and document whether it marks exposure start, midpoint, or end.
2. Add half an exposure when converting a start time to mid-exposure time.
3. Confirm the observatory location and computer-clock synchronization.
4. Convert to BJD_TDB using target coordinates, observatory location, and a maintained astronomy library.
5. Retain UTC/JD columns for audit but fit transits with one declared standard.
6. Test conversion on several frames independently.

A one-minute time-handling error is scientifically important in an ephemeris study.

### 4. AstroImageJ primary reduction

AstroImageJ should remain the program-aligned primary workflow unless the mentor directs otherwise. Record:

- AstroImageJ version and operating system.
- Calibration options and master-frame paths.
- Plate-solve method.
- Target and comparison-star coordinates.
- Aperture radius and sky-annulus radii.
- Centroiding method.
- Variable-aperture settings if used.
- Comparison ensemble and weighting.
- Removed frames and reasons.
- Detrending columns and fit order.
- Time, flux, uncertainty, and quality columns exported.

Guide: [AstroImageJ documentation](https://astroimagej.com/guides/legacy/).

### 5. Independent Python cross-check

A separate, transparent Python extraction should reproduce the broad result without copying AstroImageJ outputs as its input. It should include:

- FITS ingestion and calibrated-image validation.
- Source matching or centroid tracking.
- Aperture photometry with a declared background estimator.
- Ensemble differential photometry.
- BJD_TDB calculation.
- Quality-metric table.
- Raw and detrended light curves.
- Transit fit with the same primary assumptions.

Agreement strengthens confidence. Disagreement is a result to diagnose, not a pipeline to hide.

## Differential-photometry design

### Comparison-star eligibility

Predeclare reasonable bounds based on the actual field:

- Unsaturated in every retained frame.
- Adequate signal-to-noise ratio.
- Far from bad pixels, image edges, blends, and strong gradients.
- Stable centroid and no obvious variability.
- Similar brightness when possible.
- Similar color when catalog colors are available, reducing differential extinction.
- Present throughout the full sequence.

### Ensemble construction

Evaluate these defensible approaches:

1. Equal-weight stable comparison ensemble.
2. Inverse-variance weighted ensemble.
3. Iteratively reweighted ensemble using out-of-transit scatter only.
4. Leave-one-comparison-out jackknife.
5. Each comparison star treated as a pseudo-target.

The primary ensemble must be selected using stability and field properties, not the desired transit depth.

### Differential flux

For target counts `F_target(t)` and an ensemble reference `F_ref(t)`, begin with

`f_rel(t) = F_target(t) / F_ref(t)`

and normalize using a declared out-of-transit baseline. Propagate target, sky, read, and ensemble uncertainty rather than relying only on fitted scatter.

## Full robustness matrix — HIGH-VALUE BONUS

Run and save a grid that changes one defensible decision at a time and then a limited set of combined alternatives.

### Aperture tests

- Fixed apertures spanning the empirically useful range.
- FWHM-scaled apertures, such as multiples from approximately 0.8 to 2.5 after field testing.
- Alternative background-annulus inner and outer radii.
- Median versus robust-mode background estimators.
- Apertures centered by centroiding versus fixed sky coordinates.

### Comparison-star tests

- Primary ensemble.
- Leave-one-out ensembles.
- Bright-only and color-matched subsets when enough stars exist.
- Equal versus inverse-variance weights.
- Single-comparison diagnostic light curves.
- Pseudo-target light curves for every comparison star.

### Quality-cut tests

- No discretionary cuts beyond invalid/saturated frames.
- Frozen primary thresholds.
- Slightly stricter and looser thresholds for FWHM, sky, centroid drift, ellipticity, and flux uncertainty.
- Contiguous-event checks showing whether a conclusion depends on one isolated cluster of frames.

### Detrending tests

- Constant baseline.
- Linear and quadratic time trends.
- Airmass trend.
- Centroid x/y terms.
- FWHM, sky background, and total comparison flux terms.
- Limited physically motivated combinations.
- Transit-masked systematics fit.
- Simultaneous transit-plus-systematics fit.
- Gaussian process only as an advanced alternative, with kernels and priors declared.

Do not include a large collection of unconstrained regressors merely because they deepen the event. Compare models using out-of-transit predictive performance, AIC/BIC where appropriate, residual diagnostics, and injection preservation.

### Robustness deliverable

Create a heat map or specification curve showing fitted depth, mid-time, duration, residual RMS, red-noise factor, and model score for every acceptable configuration. Report the primary result plus the range across reasonable alternatives.

## Noise and uncertainty model — HIGH-VALUE BONUS

### Theoretical noise budget

Estimate contributions from:

- Target photon noise.
- Comparison-ensemble photon noise.
- Sky-background noise.
- Read noise.
- Dark-current uncertainty.
- Flat-field uncertainty when quantifiable.
- Atmospheric scintillation using telescope diameter, exposure, airmass, altitude, and an explicitly cited formula.

Compare the theoretical precision with observed out-of-transit scatter. A large excess identifies unresolved systematics.

### Empirical white and red noise

- Unbinned residual RMS.
- RMS versus bin size compared with the white-noise expectation.
- Time-averaging `beta` factor on timescales near ingress/egress and transit duration.
- Autocorrelation and residual time-series plots.
- Allan-deviation-style or equivalent stability diagnostic.
- Residual permutation as a correlated-noise sensitivity check.
- Block bootstrap with block lengths motivated by correlation timescale.

Inflate parameter uncertainties when red noise is present. Do not quote a white-noise MCMC interval as the full uncertainty if the residuals are correlated.

### Null and negative-control tests

- Fit the transit model at many out-of-transit control times with the same duration.
- Time-scramble or block-permute residuals while preserving appropriate correlation.
- Fit comparison-star pseudo-targets at the predicted time.
- Fit positive and negative box-shaped events to test asymmetric systematics.
- Run the full search on a transit-masked series to estimate false-alarm behavior.

The predicted window is scientifically stronger only if equally convincing dips are rare elsewhere under the same analysis.

## TESS reanalysis — HIGH-VALUE BONUS

### Sector-by-sector first

Analyze Sectors 17 and 57 independently before phase-folding or jointly fitting them.

For each sector:

1. Download official light-curve and pixel products.
2. Inspect data-release notes and quality flags.
3. Compare SAP and PDCSAP flux.
4. Inspect background, centroids, spacecraft events, and momentum dumps.
5. Plot the complete sector before masking.
6. Recover the period with a transparent BLS or TLS diagnostic without assuming that the catalog period is exact.
7. Fit individual transits and the phase-folded series.
8. Record missing or partial events and data gaps.

NASA guides: [TESS data products](https://heasarc.gsfc.nasa.gov/docs/tess/data-products.html) and [Target Pixel File tutorial](https://heasarc.gsfc.nasa.gov/docs/tess/Target-Pixel-File-Tutorial.html).

### Extraction sensitivity

- Official SAP aperture.
- Smaller and larger custom pixel masks.
- Difference-image-informed aperture.
- Simple aperture photometry versus a modern public correction pipeline, if used transparently.
- Background alternatives.
- Dilution correction with and without catalog crowding values.

Plot recovered depth against aperture size. A depth that grows as contaminating pixels are included is a blend warning.

### Standard transit-vetting diagnostics

- Odd-versus-even depths and timings.
- Search for secondary eclipses at phase 0.5 and elsewhere.
- Transit-shape/V-shape assessment.
- Individual-event consistency.
- Model-shift or equivalent uniqueness test.
- Centroid motion during transit.
- In-transit minus out-of-transit difference image.
- Aperture-contamination assessment.
- Period aliases and harmonic checks.
- Search for transit-correlated background/systematic changes.

Use [DAVE](https://arxiv.org/abs/1901.07459), current TESS Data Validation products, or equivalent diagnostics as independent evidence. A tool flag is not a final classification.

### Stellar variability and flares

- Mask transits and compute a rotation/variability periodogram.
- Compare variability across sectors.
- Identify flares and test whether any overlap a transit.
- Test whether local stellar trends bias transit depth.
- Use a simple basis or Gaussian process only when justified by variability structure.
- Avoid interpreting an apparent rotation period without checking aliases and sector length.

Activity analysis supports the transit-noise model; it should not become a separate paper unless the data warrant it.

## Gaia neighborhood and blend analysis — HIGH-VALUE BONUS

The official Gaia DR3 cone search identifies the likely target as Gaia DR3 source **306365897034483712**, with approximately:

- G = 12.5721.
- BP−RP = 1.0708.
- RUWE = 0.976.
- Position = 17.0826757°, +25.6525746°.

RUWE near one is not evidence of a resolved companion, but it also does not prove the star is single.

### Preliminary neighbor screen

The following calculation assumes full contaminant throughput and treats Gaia G-band contrast as a rough proxy for the relevant photometric band. It is a screening calculation, not a validation result.

`flux ratio = 10^(-0.4 × ΔG)`

`required contaminant eclipse fraction ≈ target transit depth / flux ratio`

| Separation | Neighbor G | ΔG from target | Approximate flux ratio | Approximate eclipse needed to mimic 2.521 ppt |
|---:|---:|---:|---:|---:|
| 17.7 arcsec | 19.966 | 7.39 | 0.0011 | 229% |
| 32.8 arcsec | 16.991 | 4.42 | 0.0171 | 14.8% |
| 37.2 arcsec | 17.488 | 4.92 | 0.0108 | 23.3% |
| 44.5 arcsec | 17.927 | 5.36 | 0.0072 | 35.0% |
| 63.6 arcsec | 17.691 | 5.12 | 0.0090 | 28.1% |

Under those deliberately crude assumptions, the 17.7 arcsec cataloged neighbor is too faint to mimic the full signal even if it disappeared completely. The brighter neighbors farther away remain plausible only to the extent that the TESS or GMU apertures transmit their flux and they undergo deep eclipses.

The final analysis must replace Gaia-G approximations with TIC TESS magnitudes, measured aperture/PRF throughput, bandpass-aware contrasts, and uncertainty. TESS pixels are approximately 21 arcsec across, so several of these stars can matter to TESS even if GMU resolves them. See [NASA’s TESS crowding guide](https://heasarc.gsfc.nasa.gov/docs/tess/UnderstandingCrowdingv2.html), [TESS telescope information](https://heasarc.gsfc.nasa.gov/docs/tess/telescope_information.html), and the [Gaia DR3 archive](https://gea.esac.esa.int/archive/).

### Required localization products

- GMU finding chart with target, every Gaia/TIC neighbor, aperture, and sky annulus.
- Light curve for every neighbor bright enough to mimic the event.
- Neighbor eclipse-depth exclusion table.
- TESS pixel map with official/custom aperture and Gaia overlays.
- TESS difference image and centroid-offset result.
- Recovered depth versus TESS aperture size.
- GMU depth versus aperture size.
- Explicit unresolved-separation limit set by seeing, plate scale, and centroid precision.

### Advanced contamination extensions

- PRF-weighted contamination model for each TESS pixel.
- Image-level injection of eclipses into each neighbor to measure leakage into the target aperture.
- Bayesian model assigning event probability among target and neighbors.
- Gaia variability and non-single-star table checks.
- Public high-resolution imaging contrast-curve integration.
- Proper-motion archival-image test: inspect the target’s present position in older imaging to constrain background objects that are now hidden under it.
- Color-dependent centroid or depth checks if multi-band data become available.

## Transit modeling — CORE / ADVANCED

### Primary model

Use a physically integrated transit model such as [`batman`](https://arxiv.org/abs/1507.08285) rather than estimating every property from a box. At minimum model:

- Mid-transit time.
- Radius ratio or depth.
- Scaled semimajor axis and inclination, or duration/impact-parameter parameterization.
- Limb darkening appropriate to TESS and R band.
- Baseline/systematics terms.
- Per-dataset jitter.
- Finite exposure integration.

Use catalog values as priors with cited uncertainty, not fixed truth, unless the data cannot constrain a parameter and the fixed choice is explicitly reported.

### Ground-only fits

Produce:

1. A prediction-conditioned fit with period/shape priors from TESS.
2. A weakly constrained event fit testing whether the data prefer the expected signal.
3. A no-transit/systematics-only model.
4. A sensitivity fit with alternative limb-darkening priors.

Compare the transit and no-transit models using residual diagnostics and an appropriate information criterion or Bayesian evidence method. Avoid calling a BIC difference a universal false-alarm probability.

### TESS-only fits

- Separate sector parameters for baseline, dilution, jitter, and activity.
- Shared orbital period and geometry when justified.
- Individual mid-times for timing analysis.
- Sector-specific depths as a diagnostic before enforcing one shared depth.

### Joint TESS + GMU fit

- Shared period, epoch, geometry, and radius ratio when testing achromatic consistency.
- Separate limb-darkening coefficients by band.
- Separate dilution, baseline, noise, and detrending terms by dataset.
- Alternative fit allowing GMU and TESS depths to differ.
- Posterior predictive checks for each dataset.

R and TESS are both red optical bandpasses, so matching depths are supportive but not a strong multi-color achromaticity test.

### Inference safeguards

- Inspect parameter identifiability and prior sensitivity.
- Use multiple MCMC chains, convergence diagnostics, effective sample size, and reproducible random seeds.
- Compare with a deterministic optimizer to catch coding errors.
- Perform synthetic-data recovery with known parameters.
- Report covariance between depth, impact parameter, limb darkening, dilution, and baseline.
- Do not infer planet mass, density, composition, or atmospheric properties from transit photometry alone.

## Ephemeris and timing analysis — CONDITIONAL HIGH-VALUE BONUS

### Planning propagation

For cycle number `N`, a simple linear prediction is

`T_N = T_0 + N P`.

If epoch and period errors are treated as independent, a planning-only uncertainty is

`sigma(T_N) = sqrt[sigma(T_0)^2 + N^2 sigma(P)^2]`.

Using the snapshot ephemeris:

- The 2023-10-12/13 event is near cycle `N = 93` and predicted BJD about **2460230.780405**.
- The simple propagated 1-sigma uncertainty at that epoch is about **4.28 minutes**.
- Near 2026-07-13, the elapsed cycle is roughly `N = 359` and the same simplified uncertainty is about **9.28 minutes**.

These numbers are planning checks. The final calculation must use the actual covariance, correct BJD_TDB convention, revised archive ephemeris, and any individual timing data. It must not assume that archive errors are Gaussian or that TTVs are absent.

### Timing workflow

1. Fit every usable TESS transit time with a consistent shape model.
2. Fit the GMU time only if ingress/egress and baseline make it identifiable.
3. Assign integer epochs and verify cycle count.
4. Fit a weighted linear ephemeris with epoch-period covariance.
5. Plot observed-minus-calculated residuals.
6. Run leave-one-out fits to show each timing’s influence.
7. Compare old and refined future transit-window uncertainty.
8. Test uncertainty inflation for timing red noise.
9. Report a quadratic ephemeris only if the number, precision, baseline, and model comparison support it.

Student-led ephemeris refinement is a real contribution when the timing is trustworthy; see the [ORBYTS ephemeris-refinement study](https://arxiv.org/abs/2005.01684) and [TESS ephemeris-recovery analysis](https://arxiv.org/abs/1906.02197).

### TTV guardrail

A single nonzero O−C point can arise from noise, red systematics, partial coverage, time-standard mistakes, or an outdated linear ephemeris. Call it a timing offset, not a TTV detection. A TTV claim requires repeated coherent deviations and comparison with a no-TTV model.

## Injection–recovery and nondetection science — HIGH-VALUE BONUS

This is the strongest way to turn a null result into quantitative research.

### Injection design

- Inject physically shaped transits spanning approximately 0.5–8 ppt depth.
- Span durations around the predicted 2.13 h value and its uncertainty.
- Span mid-times across the propagated timing window and data gaps.
- Include a focused grid at the catalog depth, duration, and timing uncertainty.
- Inject into raw or minimally processed target-level data before the reduction choices being tested.
- Inject into out-of-transit segments or synthetic baselines without overwriting a possible real event.
- Repeat across realistic red-noise realizations or residual blocks.

### Recovery rule

Freeze a rule before running the final grid. A recovery might require:

- Correct event sign.
- Mid-time within a declared tolerance.
- Recovered depth significantly above zero.
- Transit model preferred to the null model by a declared metric.
- No stronger simultaneous event in control stars.

### Deliverables

- Completeness versus depth and mid-time.
- Completeness at the catalog prediction.
- False-positive recovery rate from no-injection controls.
- Bias in recovered depth and time.
- Upper limit on detectable depth for well-covered intervals.
- Comparison across primary and alternative detrending pipelines.

If completeness for a 2.521 ppt transit is low, the correct result is “the GMU data do not test the predicted event,” not “the transit did not occur.”

## Formal false-positive validation — ADVANCED

The paper can include this only after the basic light curve and source-localization work are correct. Statistical validation is conditional on the inputs and scenario model; it is not mathematical proof.

### TRICERATOPS / TRICERATOPS+

Use the current maintained version and cite its exact release. The original method explicitly models nearby-star scenarios for TESS candidates; see [Giacalone et al. 2021](https://arxiv.org/abs/2002.00691). A 2025 extension incorporates ground-based light curves in separate bandpasses; see [TRICERATOPS+](https://arxiv.org/abs/2508.02782).

Required inputs and records:

- TESS pixel data and aperture.
- Target and neighboring-star catalog properties.
- Transit fit with uncertainty.
- Contrast curves when public and permitted.
- Ground-based localization and filter response.
- Spectroscopic constraints when public.
- Software version, random seed, prior assumptions, and number of Monte Carlo draws.
- Repeated runs demonstrating numerical stability.

Report both false-positive probability and nearby-false-positive probability, plus the contribution of each scenario. Do not quote only the most favorable run.

The original TRICERATOPS study used `FPP < 0.015` and `NFPP < 10^-3` as its statistically validated region. These are method-specific published criteria, not automatic permission to use the word “validated.” Confirm whether the current version retains the same definitions and whether the available observation set meets all input assumptions.

Advanced stability audit:

- Repeat enough random seeds to estimate run-to-run scatter rather than selecting one result.
- Perturb stellar, transit, dilution, and contrast-curve inputs within uncertainty.
- Re-run with reasonable Galactic-population and binary-prior alternatives when supported.
- Test inclusion/exclusion of each observational constraint and show how much it changes FPP/NFPP.
- Verify that ground localization, contrast curves, and catalog constraints are not double counted.
- Preserve scenario samples and convergence diagnostics.
- Compare software output with the qualitative scenario evidence matrix.
- Do not apply a multiplicity boost unless an independently vetted multi-candidate architecture actually exists and the adopted method supports it.

### `vespa` historical cross-check

Use [`vespa`](https://vespa.readthedocs.io/en/latest/) only as a transparent secondary calculation with [Morton et al. 2016](https://arxiv.org/abs/1605.02825) cited. Document its scenario coverage, stellar-population assumptions, aperture radius, secondary-eclipse constraint, photometry, and contrast curves. Because nearby contamination is central for TESS, disagreement with TRICERATOPS must be investigated rather than averaged away.

### Independent vetting evidence

- TESS Data Validation diagnostics.
- DAVE or an equivalent public vetting pipeline.
- Odd/even and secondary-eclipse limits.
- Difference image and centroid constraint.
- GMU neighbor eclipse exclusion.
- Gaia catalog and archival-image constraints.
- Public high-resolution imaging.
- Public reconnaissance spectroscopy.
- Transit-shape and stellar-density consistency.
- Chromatic-depth evidence if additional filters exist.

### Scenario-level modeling

Explicitly compare:

- Planet transiting target.
- Eclipsing binary on target.
- Hierarchical eclipsing binary.
- Background or foreground eclipsing binary.
- Planet transiting a bound companion.
- Instrumental/systematic event.

State which scenarios each observation excludes and which remain. A low model FPP is not publishable validation if key nearby stars, contrast curves, or binary constraints are absent.

### Validation stop conditions

Do not claim validation if:

- The ground-based target identity is unresolved.
- A plausible neighbor remains untested.
- The TESS difference image points away from the target.
- Odd/even or secondary evidence supports an eclipsing binary.
- Required high-resolution imaging or spectroscopy is missing from the chosen validation standard.
- FPP/NFPP is unstable to reasonable priors or catalog inputs.
- The mentor has not approved a validation-level claim.

## Physical and population context — CONDITIONAL

If the transit remains plausible, compare the candidate with short-period Neptune-size and sub-Saturn populations.

### Reproducible population figure

1. Query a dated NASA Exoplanet Archive confirmed-planet table.
2. Apply declared quality cuts on period, radius, uncertainty, and host properties.
3. Plot period versus radius, coloring by equilibrium temperature, insolation, or host temperature.
4. Mark TOI-6241.01 as a candidate with a distinct symbol.
5. Show how the point moves under its radius uncertainty.
6. Avoid mixing candidate and confirmed samples without clear styling.
7. Discuss survey selection effects and the incompleteness of validation near crowded or faint targets.

Potential question: does TOI-6241.01 occupy a sparsely populated part of period–radius space sometimes associated with the hot-Neptune desert or the Neptune/sub-Saturn transition? The boundary definition must come from cited literature; the paper must not declare desert membership by visual impression.

### Stellar consistency extension

- Construct a public-photometry spectral energy distribution.
- Compare temperature, radius, extinction, and distance with catalog values.
- Use Gaia parallax carefully, including the adopted zero-point treatment if refitting.
- Compare transit-derived stellar density with independent stellar estimates.
- Investigate a large inconsistency as a possible blend, eccentricity, or model problem.

### Follow-up value

Use predicted mass, RV amplitude, TSM, and ESM only to motivate possible observations. Label all predicted quantities. A transit study cannot measure mass or atmosphere.

## Complete extension catalog

| Extension | Level | Scientific value | Dependency | Concrete output |
|---|---|---|---|---|
| Second GMU transit | GATED / FUTURE | Repeats the event and improves timing | Telescope approval and observability | Independent R-band light curve and timing |
| Multi-color ground photometry | GATED / FUTURE | Tests chromatic eclipsing-binary/blend scenarios | Approved observations in separated filters | Depth-versus-band table |
| Simultaneous multi-site observation | GATED / FUTURE | Separates local weather/systematics from astrophysical signal | Partner observers and synchronized plan | Cross-site light curves |
| Public TFOP light-curve comparison | CONDITIONAL | Independent recovery and baseline | Public/permitted product exists | Homogeneous re-fit |
| High-resolution imaging | GATED / FUTURE | Excludes close companions | Public contrast curve or new approved data | Contrast curve in validation model |
| Reconnaissance spectroscopy | GATED / FUTURE | Detects stellar binaries and refines host parameters | Public spectrum or new approved data | RV/bisector/stellar table |
| Precision radial velocities | GATED / FUTURE | Measures mass or rejects stellar companion | New instrument campaign | RV orbit and mass |
| Archival proper-motion test | HIGH-VALUE BONUS | Constrains background contaminants now hidden at target position | Sufficient target proper motion and old images | Multi-epoch image panel |
| Gaia non-single-star/variability audit | HIGH-VALUE BONUS | Adds catalog companion evidence | Current Gaia products | Astrometry/variability table |
| ASAS-SN/ZTF variability | CONDITIONAL | Checks long-period variability and rotation aliases | Adequate cadence and precision | Periodogram and phase curve |
| Pixel-level eclipse injection | ADVANCED | Quantifies neighbor leakage | TESS PRF/pixel model | Contaminant recovery map |
| Bayesian source assignment | ADVANCED | Combines centroid, flux, and neighbor evidence | Valid likelihoods and priors | Posterior source probabilities |
| Independent extraction packages | HIGH-VALUE BONUS | Tests pipeline dependence | Same raw products | Cross-pipeline comparison |
| Gaussian-process systematics | ADVANCED | Models correlated noise | Enough baseline and constrained kernel | GP versus simple-baseline comparison |
| Wavelet/celerite noise check | ADVANCED | Alternative red-noise treatment | Explainable implementation | Parameter sensitivity table |
| Transit-search algorithm comparison | HIGH-VALUE BONUS | Tests period/alias recovery | TESS light curves | BLS versus TLS recovery plot |
| Individual-transit variation | CONDITIONAL | Tests depth/time consistency | Several adequate TESS events | Event-by-event parameter plot |
| TTV search | CONDITIONAL | Could indicate additional dynamics | Multiple precise timings | O−C model comparison |
| Transit-depth variation | CONDITIONAL | Checks activity, dilution, or changing contamination | Multiple comparable events | Depth-versus-epoch plot |
| Secondary/phase-curve upper limit | CONDITIONAL | Constrains stellar false positives | Adequate TESS precision | Eclipse/phase amplitude limit |
| Flare–transit interaction audit | CONDITIONAL | Prevents activity bias | Detectable flares | Flare mask and sensitivity result |
| Catalog SED fit | ADVANCED | Independently checks host radius | Clean public photometry | SED posterior and comparison |
| Population-context analysis | HIGH-VALUE BONUS | Explains why candidate matters | Dated archive query | Period–radius figure |
| False-positive model ensemble | ADVANCED | Tests dependence on validation software | Complete inputs | TRICERATOPS/vespa evidence matrix |
| Reproducible data release | HIGH-VALUE BONUS | Makes work auditable and reusable | Permission and licenses | Code, manifests, configs, derived tables |
| Educational methods appendix | HIGH-VALUE BONUS | Demonstrates student-led reproducibility | Complete analysis log | Stepwise appendix |
| Future-transit scheduler | HIGH-VALUE BONUS | Converts ephemeris into useful planning windows | Refined covariance | Calendar/table with uncertainty windows |
| Machine-learning anomaly flagging | ADVANCED / POTENTIAL | Independent frame/systematics triage | Labeled diagnostics and interpretable model | Flag comparison with human cuts |
| Hierarchical multi-event model | ADVANCED | Separates shared transit shape from event noise | Enough individual events | Partial-pooling posterior |
| Residual search for another transiting signal | ADVANCED / SEPARATE | Could reveal a multi-candidate architecture or aliases | TOI-6241 model removed and search trials controlled | Independent BLS/TLS search and injection completeness |
| Repeatable transit-shape anomaly test | GATED / POTENTIAL | Could test spots, gravity darkening, rings, or moons | Multiple high-S/N events with repeatable structure | Competing shape-model evidence |
| Rossiter–McLaughlin observation | GATED / FUTURE | Could constrain projected orbital alignment | Confirmed planet, suitable rotation, approved spectroscopy | Spectroscopic transit model |
| Atmospheric escape spectroscopy | GATED / SEPARATE | Could test metastable helium or other escape tracers | Validated/confirmed target and suitable instrumentation | Separate spectroscopic analysis |
| JWST atmosphere analysis | GATED / SEPARATE PAPER | Atmosphere context if approved | Approved campus light curve and mentor authorization | Separate project and manuscript |

Machine learning must not classify the planet by opaque score or replace physical vetting. Its only defensible near-term role is frame-quality or anomaly triage with interpretable features and human review.

## Outcome decision tree

| Observed outcome | Best scientific interpretation | Best next analysis | Claim ceiling |
|---|---|---|---|
| Clear event on target, compatible time/depth/duration | Independent support for candidate transit | Robustness, TESS joint fit, localization, timing | Additional validation evidence |
| Event on target but depth differs | Possible bandpass, dilution, activity, or systematics issue | Aperture/dilution/chromatic/model audit | Transit-like event with discrepancy |
| Stronger event on neighbor | Likely nearby false-positive source | Quantify leakage and compare TESS centroid | Nearby false-positive evidence |
| No visible event, high injection recovery | Informative nondetection | Ephemeris, timing window, weather, prediction audit | Excludes predicted event to stated completeness |
| No visible event, low injection recovery | Inconclusive night | Noise budget and upper-limit study | No test of candidate |
| Partial transit only | Timing/depth degeneracy | Prediction-conditioned fit and injection tests | Partial event consistent or inconsistent with prediction |
| Signal depends on one reduction choice | Analysis instability | Specification curve and systematics diagnosis | No robust detection |
| Field is TOI-3718.01 | Wrong dataset for this project | Obtain correct files or approved target switch | No TOI-6241 result |
| Time metadata unreliable | Timing analysis invalid | Recover clock provenance or omit timing claims | Photometric shape only if identity is secure |

## Result-dependent title options

Use a title only after the outcome is known.

### If a robust on-target event is recovered

**Ground-Based Recovery and Source Localization of a Transit-Like Event for TOI-6241.01**

### If the main contribution is robustness and TESS comparison

**A Robustness-Audited Ground–Space Photometric Analysis of the TESS Candidate TOI-6241.01**

### If the result is an informative nondetection

**Injection–Recovery Limits from a Ground-Based Follow-Up Night of TOI-6241.01**

### If a neighbor is implicated

**Ground-Based Source Localization of a Possible False-Positive Signal near TOI-6241.01**

### If the data remain inconclusive

**Assessing the Detectability of TOI-6241.01 in GMU Ground-Based Photometry**

## Figure plan

### Main-text figures

1. Full TESS Sector 17 light curve with transit windows and quality gaps.
2. Full TESS Sector 57 light curve with the same visual conventions.
3. Sector-separated and combined phase-folded TESS transits.
4. TESS Target Pixel File image with aperture, target, and Gaia/TIC neighbors.
5. GMU field image with target, comparisons, neighbors, aperture, and north/east arrows.
6. GMU observing diagnostics: sky, airmass, FWHM, ellipticity, centroid drift, and comparison flux.
7. Primary raw and detrended GMU differential light curve with predicted and fitted contact times.
8. Best-fit transit and null models with residual panels.
9. Aperture/comparison/detrending robustness specification plot.
10. Neighbor-star light curves and eclipse-depth exclusion plot.
11. TESS difference image and centroid localization.
12. RMS-versus-bin-size red-noise diagnostic.
13. Injection–recovery completeness map or curve.
14. O−C timing diagram if timing is justified.
15. Period–radius context plot if physical context is included.

### Appendix figures

- Master calibration products and calibration diagnostics.
- Plate-solution overlays for early/middle/late frames.
- Every comparison-star pseudo-target light curve.
- Depth versus aperture radius.
- Full model-selection grid.
- Odd/even, secondary, and model-shift diagnostics.
- Sector-specific individual-transit panels.
- MCMC convergence and posterior plots.
- Prior sensitivity and posterior predictive checks.
- Null-time event-search distribution.
- Archival proper-motion image panel.
- Statistical-validation scenario probability breakdown.

Every figure must state data source, time system, flux normalization, binning, excluded points, and whether the curve is raw or detrended.

## Table plan

1. Target and host-star properties with provenance and snapshot date.
2. GMU observing log and weather/technical events.
3. File manifest summary and calibration inventory.
4. Primary photometry and detrending configuration.
5. Gaia/TIC neighbor and contamination calculations.
6. TESS and GMU transit-fit parameters.
7. Individual transit times and O−C values.
8. Robustness grid summary.
9. Injection–recovery thresholds and completeness.
10. Vetting evidence by false-positive scenario.
11. Software, package versions, data releases, and random seeds.

## Paper structure aligned with the lectures

### Abstract

Five-part structure:

1. Candidate and scientific motivation.
2. GMU and TESS data used.
3. Core reduction/model/localization method.
4. Numerical result with uncertainty and detection/completeness statement.
5. Precisely bounded implication for candidate vetting.

Write the abstract last. It must include the null or contradictory result if that is what the data show.

### Introduction

- Transit method and why TESS candidates need ground follow-up.
- False positives and TESS crowding.
- Why source-localized ground photometry is useful.
- TOI-6241.01 properties and candidate status.
- Literature/novelty audit.
- Research questions and paper organization.

### Observations and data

- GMU telescope, camera, filter, cadence, and observing window.
- Calibration frames and observing conditions.
- Data-identity proof.
- TESS sectors and products.
- Gaia/catalog/follow-up products.

### Methods

- Calibration and time conversion.
- Differential photometry.
- Comparison-star and aperture selection.
- Detrending and uncertainty.
- TESS extraction and vetting.
- Transit models.
- Robustness and injection–recovery.
- Localization and false-positive methods.

### Results

Report observations before interpretation:

- Achieved precision and coverage.
- Primary light curve and model comparison.
- Depth, duration, and timing constraints.
- Robustness range.
- Neighbor/localization results.
- Injection completeness.
- TESS consistency and timing if completed.

### Discussion

- Compatibility with predicted event.
- Possible false-positive/systematics interpretations.
- Ground-versus-TESS depth and timing.
- What is newly constrained.
- Limitations and untested scenarios.
- Candidate population and follow-up value when warranted.

The program guide specifically expects discussion of timing, depth/duration, and false-positive possibilities; see [paper discussion guidance](../data_and_lectures/Paperstructure.pdf#page=15).

### Conclusion

- One sentence answering the primary question.
- Quantitative evidence.
- Claim-ladder ceiling.
- Most valuable next observation.

### Appendices

- Full settings and manifest.
- Robustness grid.
- Comparison-star diagnostics.
- Injection protocol.
- Priors and convergence.
- Additional vetting and catalog tables.
- Reproducibility statement and AI disclosure.

## Proposed repository structure

```text
data/
  raw/
    gmu/toi6241/
    tess/toi6241/sector17/
    tess/toi6241/sector57/
  external/
    gaia/
    exoplanet_archive/
    exofop/
    archival_imaging/
  processed/
    gmu/toi6241/
    tess/toi6241/
  manifests/
configs/
  toi6241_primary.yaml
  toi6241_robustness.yaml
docs/
  toi-6241-research-maximization-plan.md
  analysis-ledger.md
  novelty-ledger.md
notebooks/
  01_data_identity.ipynb
  02_gmu_quality_control.ipynb
  03_tess_sector_analysis.ipynb
  04_neighbor_analysis.ipynb
  05_transit_model.ipynb
  06_injection_recovery.ipynb
  07_ephemeris.ipynb
  08_false_positive_validation.ipynb
src/
  toi6241/
outputs/
  toi6241/
    figures/
    tables/
    models/
    logs/
    paper/
```

Raw data should remain immutable. If raw FITS products cannot be committed because of size, licenses, or privacy, store a manifest, checksums, source/retrieval instructions, and a `.gitignore` rule. Do not silently overwrite downloaded files.

## Reproducibility standard

- Pin the Python environment and record the platform.
- Record AstroImageJ version and every manual setting.
- Save archive query text, retrieval date, and raw returned table.
- Save random seeds and sampler settings.
- Store machine-readable model priors and configuration.
- Use checksums for raw inputs and important derived tables.
- Generate paper figures from scripts, not manually edited plots.
- Keep plotting choices separate from numerical analysis.
- Write units and time standards in table column names.
- Add tests for time conversion, aperture flux on synthetic images, transit injection recovery, and ephemeris cycle counting.
- Run a clean-environment reproduction before submission.
- Preserve negative, null, and failed model runs in logs.
- Cite software and data releases.
- Obtain permission before publishing GMU or restricted follow-up data.

## Prioritized execution roadmap

### Phase 0 — Unblock and prove identity

- [ ] Obtain the correct GMU files or mentor-approved resolution.
- [ ] Build checksummed manifest.
- [ ] Plate solve early/middle/late frames.
- [ ] Confirm TOI-6241 coordinates and field.
- [ ] Inventory calibration data and logs.
- [ ] Freeze a current target-property snapshot.

**Stop rule:** do not start TOI-6241 ground-based transit interpretation until identity passes.

### Phase 1 — Complete the internship core

- [ ] Calibrate the correct night in AstroImageJ.
- [ ] Select target/comparison apertures using field quality.
- [ ] Produce raw and detrended differential light curves.
- [ ] Convert times to BJD_TDB.
- [ ] Compare observed and predicted timing, depth, and duration.
- [ ] Check nearby stars.
- [ ] Post the required approved light-curve progress update.
- [ ] Draft outcome-neutral methods and data sections.

### Phase 2 — Highest-value bonuses

- [ ] Independent Python extraction.
- [ ] Aperture/comparison/detrending specification curve.
- [ ] Red-noise and theoretical-noise budget.
- [ ] TESS Sector 17 and 57 separate reanalysis.
- [ ] TESS aperture/Gaia overlay and difference image.
- [ ] Injection–recovery completeness.
- [ ] Novelty audit and ledger.

### Phase 3 — Timing and physical context

- [ ] Individual TESS transit times.
- [ ] GMU time only if identifiable.
- [ ] Linear ephemeris and covariance.
- [ ] O−C and leave-one-out checks.
- [ ] Dated period–radius context figure.
- [ ] Future-observation window table.

### Phase 4 — Advanced false-positive analysis

- [ ] Assemble public contrast, spectra, Gaia, and centroid constraints.
- [ ] Run TRICERATOPS with repeated seeds.
- [ ] Run `vespa` only as a documented cross-check.
- [ ] Complete scenario evidence matrix.
- [ ] Review claim ceiling with mentor.

### Phase 5 — Potential new observations

- [ ] Request a second transit only if it resolves a real limitation.
- [ ] Prioritize multi-color or simultaneous coverage based on remaining scenarios.
- [ ] Seek public/new high-resolution imaging or spectroscopy only with approval.
- [ ] Treat RV and atmosphere work as later independent projects.

### Phase 6 — Paper and release

- [ ] Freeze primary analysis.
- [ ] Generate all figures and tables from reproducible scripts.
- [ ] Write result-matched title and abstract.
- [ ] Complete limitations and claim-ladder audit.
- [ ] Re-run novelty and catalog checks.
- [ ] Complete mentor review.
- [ ] Reproduce in a clean environment.
- [ ] Release only permitted data and products.

## Value-versus-effort priority

### Must do

1. Correct data identity.
2. Trustworthy calibration and differential photometry.
3. Correct BJD_TDB timing.
4. Neighbor light curves/source localization.
5. Predicted-versus-observed time, depth, and duration.
6. Honest uncertainty and result-matched claims.

### Best scientific return for added work

1. TESS sector-by-sector reanalysis.
2. Aperture/comparison/detrending robustness map.
3. Injection–recovery completeness.
4. Gaia/TESS aperture contamination analysis.
5. Independent Python reduction.
6. Linear ephemeris refinement if the GMU time is precise.

### Advanced only after the above is solid

1. TRICERATOPS+ and `vespa` model comparison.
2. Bayesian source assignment.
3. GP or hierarchical multi-event models.
4. SED refitting.
5. TTV or depth-variation analysis.

### Moonshots / separate projects

1. New multi-color network campaign.
2. High-resolution imaging.
3. Reconnaissance and precision RV spectroscopy.
4. Atmosphere/JWST analysis after explicit approval.

## Risk register

| Risk | Consequence | Early diagnostic | Mitigation |
|---|---|---|---|
| Wrong target files | Entire study invalid | Header and plate-solve mismatch | Hard identity gate |
| Missing flats/darks | Biased photometry | Calibration inventory | Ask mentor; use substitutes only with documented approval/tests |
| 2.5 ppt signal below precision | Nondetection not informative | Out-of-transit RMS and injection | Completeness analysis and honest ceiling |
| Clouds or moisture | Correlated dips | Sky/comparison flux | Quality metrics, control stars, red-noise model |
| Tracking/focus jumps | Aperture losses | Centroid/FWHM | Variable aperture, segment checks, exclusion only by frozen rules |
| Bad clock or time convention | False timing offset | Header/time audit | Mid-exposure BJD_TDB cross-check |
| Nearby eclipsing source | False target attribution | Neighbor curves/difference image | Source-localization suite |
| Detrending removes/invents transit | Biased depth | Injection and model grid | Transit-masked/simultaneous alternatives |
| Variable comparison star | Artificial target dip | Pseudo-target curves | Ensemble jackknife |
| Partial baseline | Depth/time degeneracy | Coverage plot | Strong priors or inconclusive label |
| TESS crowding/dilution error | Sector depth disagreement | Aperture sensitivity | Pixel-level contamination model |
| Catalog changes | Outdated properties | Snapshot comparison | Refresh before final draft |
| Multiple testing | Inflated significance | Analysis ledger | Frozen primary test and null controls |
| Advanced model instability | Misleading FPP/posterior | Seed/prior sensitivity | Repeat, cross-tool, simplify, disclose |
| Overclaiming validation | Scientifically incorrect conclusion | Claim-ladder audit | Mentor review and scenario matrix |
| Scope overload | Core paper unfinished | Roadmap status | Finish phases in order |
| Restricted data used improperly | Publication/ethics problem | Provenance audit | Obtain permission or omit |

## Approval gates

| Gate | Approval/evidence needed | Work unlocked |
|---|---|---|
| Correct ground data | Mentor or conclusive identity audit | GMU reduction and interpretation |
| First light curve | Program/mentor review process | Full paper analysis progression |
| Internal TFOP information | Explicit permission | Restricted constraints or notes |
| Statistical validation claim | Complete inputs plus mentor review | “Validated” language, if thresholds are met |
| New telescope observations | Mentor/telescope approval | Second night, multi-color, multi-site work |
| JWST/atmosphere work | Lecture-defined approval after campus light curve | Separate atmosphere project/paper |
| Public release | Data-owner and license permission | Repository data/products release |

## Definition of a successful project

The project succeeds if it produces a reproducible and correctly bounded answer to the primary question. Success does not require a pretty transit or planet validation.

A maximum-quality final result includes:

- Proven data identity.
- Auditable calibration and time conversion.
- A primary analysis frozen independently of the desired outcome.
- Ground and TESS light curves with source localization.
- Robust uncertainty including correlated noise.
- Quantified sensitivity for any nondetection.
- Explicit false-positive scenarios and remaining limitations.
- Reproducible code, figures, tables, and provenance.
- Claims that stop at the supported rung.
- A clear next observation selected because it resolves the largest remaining uncertainty.

## Immediate next actions

While waiting for the correct GMU data:

1. Create the current NASA Archive/ExoFOP/Gaia snapshot with query dates.
2. Download and inspect TESS Sectors 17 and 57 independently.
3. Build a TESS pixel/aperture map with the Gaia neighbors.
4. Create the analysis and novelty ledgers.
5. Write the frozen primary-analysis configuration without viewing an alleged TOI-6241 GMU transit.
6. Prepare the manifest and plate-solve checks so the correct files can be authenticated immediately.
7. Recheck whether public imaging, spectroscopy, or time-series constraints now exist.

Do not use the TOI-3718.01 files to manufacture a TOI-6241.01 result.

## Reference starting set

### Program materials

- [Lecture 1: introduction and AI guidance](../data_and_lectures/Lecture1_intro.pdf)
- [Lecture 2: NASA missions](../data_and_lectures/lecture2_NASAmissions.pdf)
- [Lecture 3: working with data](../data_and_lectures/lecture3_workingwithdatapt1%20%281%29.pdf)
- [Lecture 4](../data_and_lectures/Lecture4_Schar2025_clean_redesign.pdf)
- [Lectures 5–6](../data_and_lectures/Lecture56_clean_redesign%20%281%29.pdf)
- [Lecture 8 and research-project sequence](../data_and_lectures/Lecture8_Schar2025_clean_redesign.pdf)
- [Paper structure](../data_and_lectures/Paperstructure.pdf)
- [Example GMU TOI-5372.01 paper](https://science.gmu.edu/sites/default/files/2024-10/Ground-Based%20Follow-Up%20Observations%20of%20TESS%20Object%20of%20Interest%20%28TOI%29%205372.01.pdf)

### Mission, catalog, and follow-up documentation

- [NASA Exoplanet Archive TOI table](https://exoplanetarchive.ipac.caltech.edu/cgi-bin/TblView/nph-tblView?app=ExoTbls&config=TOI)
- [TESS data products](https://heasarc.gsfc.nasa.gov/docs/tess/data-products.html)
- [TESS Target Pixel File tutorial](https://heasarc.gsfc.nasa.gov/docs/tess/Target-Pixel-File-Tutorial.html)
- [Understanding TESS crowding](https://heasarc.gsfc.nasa.gov/docs/tess/UnderstandingCrowdingv2.html)
- [TESS telescope information](https://heasarc.gsfc.nasa.gov/docs/tess/telescope_information.html)
- [TFOP ground-based follow-up overview](https://asd.gsfc.nasa.gov/archive/tess/ground_based_followup.html)
- [Gaia DR3 archive](https://gea.esac.esa.int/archive/)
- [Gaia DR3 documentation hub](https://www.cosmos.esa.int/web/gaia/dr3)
- [AstroImageJ documentation](https://astroimagej.com/guides/legacy/)

### Methods and advanced validation

- [Kreidberg 2015: `batman` transit model](https://arxiv.org/abs/1507.08285)
- [Kostov et al. 2019: DAVE vetting pipeline](https://arxiv.org/abs/1901.07459)
- [Giacalone et al. 2021: TRICERATOPS](https://arxiv.org/abs/2002.00691)
- [Gomez Barrientos et al. 2025: TRICERATOPS+](https://arxiv.org/abs/2508.02782)
- [Morton et al. 2016: `vespa`](https://arxiv.org/abs/1605.02825)
- [Morton and Johnson 2011: false-positive probability framework](https://arxiv.org/abs/1101.5630)
- [ORBYTS student ephemeris-refinement study](https://arxiv.org/abs/2005.01684)
- [TESS ephemeris-recovery study](https://arxiv.org/abs/1906.02197)

The final bibliography should cite the exact software versions, data releases, catalog queries, and scientific definitions actually used. This list is a starting set, not permission to cite a source without reading it.
