# TOI-3505.01 false-positive evidence summary

The program's paper guidance asks that the public follow-up archive be
used to gather and present the imaging, spectroscopic, and TESS light-
curve evidence bearing on whether a candidate is a false positive. This
record collects that evidence for each standard scenario.

It is an evidence summary, not a statistical validation. No false-positive
probability is computed.

## Scenarios

### Nearby eclipsing binary (NEB)

- Evidence: Image-level screen of every catalogued neighbour within 2.5 arcmin bright enough to mimic the signal, evaluated at the 2022 schedule window with the TFOP band correction applied.
- Assessment: Disfavoured. No measured neighbour shows an eclipse-shaped event of the required depth, independent public NEB checks report the same, and the SPOC difference-image centroid stays within three sigma of the target in both sectors.
- Not excluded: Neighbours too faint for a decisive limit, and any source inside the target aperture.

### Blended eclipsing binary (BEB)

- Evidence: Wavelength dependence of the reported MuSCAT2 four-band depths.
- Assessment: Weakly disfavoured. The reported depths show no apparent monotonic trend with wavelength, but the absence of per-band errors makes this a qualitative constraint.
- Not excluded: A blend whose colour matches the target, and any blend at separations below the follow-up resolution.

### Eclipsing binary on the target (EB)

- Evidence: SPOC odd/even and weak-secondary diagnostics, plus the amplitude of the public reconnaissance velocities.
- Assessment: Disfavoured. No significant odd/even difference, no secondary eclipse detection, no SPOC eclipsing-binary flag, and a velocity span far below any stellar companion.
- Not excluded: A definitive mass requires the multi-order velocity analysis the public notes leave open.

### Unresolved close companion

- Evidence: SOAR speckle and Shane adaptive-optics imaging both resolve a companion near 0.51 arcsec that no dataset in this work can separate.
- Assessment: Not addressed. This is the limiting scenario and the reason no validation, dilution correction, or radius is claimed.
- Not excluded: Everything; the companion remains unresolved here.

## Overall

The nearby and on-target eclipsing-binary scenarios are disfavoured by the available screens, while the blended-binary chromatic check is qualitative. The unresolved 0.517-arcsecond companion is not addressed, so this summary supports continued candidate-level follow-up without validating the object.

## Limits

- An evidence summary, not a statistical validation. No false-positive probability is computed and no validation framework is run.
- External report values are used as published and are not re-reduced.
- The MuSCAT2 depths come from a tentative partial event and carry no published per-band uncertainties.
- The velocity bound assumes a host mass and treats a range of three reconnaissance velocities as a span, not a fitted amplitude.

## Products

- `false_positive_assessment.json` - every scenario and calculation.
- `01_false_positive_tests.png` and `.svg` - publication figure.

Regenerate with:

```bash
.venv/bin/python src/assess_toi3505_false_positive.py
```
