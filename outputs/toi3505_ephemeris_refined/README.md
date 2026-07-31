# Refined TOI-3505.01 linear ephemeris

Produced by `src/refine_toi3505_ephemeris.py`. This supersedes the four-sector
ephemeris in `outputs/toi3505_tess_analysis` for anything the posters quote.
The earlier fit is kept as-is; it is reproduced here as the `qlp_only` row.

## Adopted result

| | |
| --- | --- |
| Period | 2.9151516 +/- 0.0000049 d |
| Epoch, catalog cycle zero | BJD_TDB 2459793.53115 +/- 0.00088 |
| Events | 27 across Sectors 14, 41, 54, 81 |
| Baseline | 5.05 years |
| Reduced chi-square | 1.30 (errors scaled by 1.14) |
| Residual scatter | 8.5 minutes |

Agreement: 0.32 sigma from the TOI catalog period, 0.71 sigma from the SPOC
multi-sector fit, 0.68 sigma from the QLP-only version of this same fit.

Precision: 2.40 times tighter than the TOI catalog value, 1.43 times tighter
than the SPOC multi-sector value. Both comparisons belong on the poster; the
SPOC one is the harder and more honest of the two.

## Three changes from the earlier fit

**1. The transit shape is a trapezoid, not a box.** The official geometry has
b = 0.916 and Rp/R* = 0.0618, so T23/T14 = 0.378: ingress and egress take up
most of the event. A box has infinitely sharp edges, which reports mid-transit
times that look tighter than the data support and misfits the real shape. The
`box_shape_control` row is the same data through a box model and returns a
reduced chi-square of 54 against the trapezoid's 1.3.

This is why the earlier fit had to inflate its errors by 2.30. Most of that
penalty was the wrong model shape, not real timing scatter. The trapezoid needs
a factor of only 1.14.

**2. Sectors 54 and 81 use the SPOC two-minute light curves.** QLP is the only
pipeline covering Sectors 14 and 41, but SPOC published two-minute data for the
other two, and a two-minute sample constrains a mid-transit time far better than
a ten-minute Full Frame Image sample of the same transit. QLP and SPOC measure
the same photons, so only one of them may represent any given transit; the
better-sampled one is chosen per sector.

**3. One event is rejected for partial window coverage.** Sector 81 cycle 247
sits on the edge of a downlink gap. It keeps points near both ends of the fit
window, which is all the earlier rule checked, but fills less than half of it.
With the middle missing, the local baseline and the transit depth stop being
separable and the fit returns a 24 ppt event where the real signal is 3 ppt.

That event was accepted by the earlier pipeline and is in its published 25.
Accepted events now cover 75 to 100 per cent of the window and the two rejected
ones cover 34 and 46 per cent, so any threshold between those groups gives the
same answer.

## Prediction windows

| Date | Catalog | This fit |
| --- | --- | --- |
| 2026-08-01 | 9.0 min | 3.3 min |
| 2027-08-01 | 11.0 min | 4.2 min |
| 2030-08-01 | 17.1 min | 6.7 min |

## Files

| File | Contents |
| --- | --- |
| `ephemeris_refined.json` | All four fits, comparisons, and forward propagation |
| `event_times_best_per_sector.csv` | Adopted event list, with coverage and shape columns |
| `event_times_qlp_only.csv` | QLP everywhere, for the like-for-like comparison |

## Reproduce

```bash
cd src && ../.venv/bin/python refine_toi3505_ephemeris.py
```

Regression tests are in `tests/test_refine_toi3505_ephemeris.py`.

## Limits

The trapezoid is a shape approximation, not a limb-darkened model. It captures
the long ingress that a box misses, but a Mandel-Agol fit with a stellar density
prior is still the right way to pin the geometry down, and it stays on the
future-work list. No claim here depends on the ground data.
