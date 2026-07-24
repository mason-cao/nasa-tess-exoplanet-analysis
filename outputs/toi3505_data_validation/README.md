# TOI-3505.01 SPOC report comparison

This folder compares the official SPOC Data Validation products with the project's simpler TESS measurements. Official reports were found for Sectors 54 and 81. The MAST search did not return a SPOC target observation, and therefore no matching per-sector report, for Sectors 14 or 41. SPOC later produced a combined report from the Sector 54 and 81 observations.

## Main results

- Sector 54 SPOC fit: 2.914773 days, 3.226 ± 0.174 ppt, duration 2.77 hours.
- Sector 81 SPOC fit: 2.915167 days, 3.409 ± 0.173 ppt, duration 2.67 hours.
- Combined Sector 54 and 81 SPOC fit: 2.91514558 days, 3.2919 ± 0.1185 ppt, duration 2.711 ± 0.098 hours, using 17 observed transits.
- The combined run is labeled s0014-s0086 because that is the search range. The DR122 target table, DVT sector vector, and XML difference images all agree that only Sectors 54 and 81 contributed to this target. Its time series uses 10-minute bins.
- Odd/even differences are 0.22 and 0.62 sigma.
- The strongest secondary-event statistics are 2.52 and 2.11.
- SPOC's centroid offsets from the TIC position are 1.83 ± 2.55 arcsec and 1.98 ± 2.60 arcsec.
- The combined report gives an odd/even difference of 0.78 sigma, a strongest secondary statistic of 2.41, and a mean TIC-position offset of 2.62 ± 2.65 arcsec.
- Manual review of Appendix B in all three full reports found the statement "This target did not trigger any alerts."

The reports support the limited conclusion that SPOC recovered the same signal and did not flag a strong odd/even difference, secondary event, or large centroid displacement in these observations. They do not prove that the signal is planetary. The project and SPOC calculations reuse the same TESS observations, and TESS cannot resolve the known 0.517-arcsec companion.

## Files

- `official_dv_metrics.csv`: selected values read from the SPOC XML, with DVT FITS header cross-checks.
- `comparison_with_project.csv`: QLP and project SPOC box fits beside the official SPOC transit fit and centroid result.
- `official_multisector_tce.csv`: the combined fit with checks against the DVT FITS, Exo.MAST table, DR122 target table, and product list.
- `official_multisector_difference_images.csv`: the Sector 54 and 81 difference-image entries inside the combined report.
- `analysis_summary.json`: the same results in one machine-readable record.
- `01_spoc_report_comparison.png` and `.svg`: compact depth and centroid comparison.

The odd/even value in sigma is the square root of the test statistic, matching the full report table. SPOC centroid uncertainties include the report's systematic error floor; the project errors come from event resampling, so their error bars are not interchangeable.

## Official references

- [TESS data products](https://heasarc.gsfc.nasa.gov/docs/tess/data-products.html)
- [MAST TESS archive](https://archive.stsci.edu/missions-and-data/tess)
- [Exo.MAST light-curve services](https://exo.mast.stsci.edu/docs/dvdata_ws.html)
- [TESS Data Release 122](https://archive.stsci.edu/missions/tess/doc/tess_drn/tess_multisector_14_86_drn122_v01.pdf)

Rows written: 2 official report summaries and 2 project comparisons, plus 1 combined TCE row and 2 combined-report difference-image rows.
