# TOI-3505.01 TESS pixel checks

This folder contains the custom-aperture, field, difference-image, dilution-screening, and light-curve injection checks for the four TESS sectors.

## Fixed choices

- The same 15x15-pixel TESScut size is used in every sector.
- Circular radii 1.5, 2.0, 2.5, 3.0, 3.5, 4.0 pixels are compared.
- The 3.0-pixel aperture is the named reference aperture because it matches the QLP best-aperture radius in the downloaded headers.
- The depth comparison fixes the duration to 2.004 hours.

## Reference-aperture results

| Sector | Observed depth (ppt) | Difference-centroid offset (arcsec) | One-pixel check |
|---:|---:|---:|:---:|
| 14 | 1.962 +/- 0.293 | 4.5 | yes |
| 41 | 1.498 +/- 0.200 | 2.6 | yes |
| 54 | 1.510 +/- 0.194 | 1.6 | yes |
| 81 | 1.690 +/- 0.231 | 5.1 | yes |

A one-pixel match is only a TESS-scale localization. It does not distinguish the target primary from the 0.517-arcsec companion.

## Dilution screen

The TIC contamination ratio is 0.547. Using delta-I = 1.7 mag as a rough TESS-band proxy adds a companion-to-target flux ratio of 0.209. These values are kept as scenario arithmetic, not applied as a final correction, because QLP's crowding treatment and the companion's TESS-band contrast still need to be resolved.

## Injection scope

The injection tests add box-shaped dips to the extracted 3-pixel light curves at six null phases. The main recovery result keeps the control signal at each phase, so the scatter shows real phase-dependent structure. A separate increment column subtracts that control only to check the fitting arithmetic. This tests light-curve recovery, but it does not simulate target-only pixel response or repair an uncertain dilution model.

Injection summary rows: 12.

## Main files

- `depth_vs_aperture.csv`
- `difference_image_localization.csv`
- `tic_neighbor_screen_60arcsec.csv`
- `dilution_screen.json`
- `light_curve_injections.csv` and `light_curve_injection_summary.csv`
