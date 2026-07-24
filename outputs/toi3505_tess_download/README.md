# TOI-3505.01 TESS download record

This folder records the public TESS products used by the project.

- QLP light curves provide one common reduction for Sectors 14, 41, 54, and 81.
- SPOC and TESS-SPOC light curves provide same-observation pipeline checks in
  Sectors 54 and 81.
- TESScut files support the same custom-aperture test in all four sectors.
- SPOC target-pixel and per-sector Data Validation products are available for
  Sectors 54 and 81. The search found no corresponding SPOC target observation
  for Sectors 14 or 41.
- SPOC Data Release 122 also contains a combined multi-sector report. Its run
  is named `s0014-s0086`, but the target-information table and FITS sector
  vector show that only Sectors 54 and 81 contributed to TOI-3505.01. The
  combined time series uses 10-minute bins.

The checksum manifest includes FITS, PDF, XML, JSON, and text products. The
official SPOC reports are comparisons and diagnostic records, not independent
observations.

Official references:

- [QLP high-level science products](https://archive.stsci.edu/hlsp/qlp)
- [TESS data products](https://heasarc.gsfc.nasa.gov/docs/tess/data-products.html)
- [Exo.MAST light-curve services](https://exo.mast.stsci.edu/docs/dvdata_ws.html)
- [TESS Data Release 122](https://archive.stsci.edu/missions/tess/doc/tess_drn/tess_multisector_14_86_drn122_v01.pdf)
- [TESS target-pixel tutorial](https://heasarc.gsfc.nasa.gov/docs/tess/Target-Pixel-File-Tutorial.html)

The QLP archive page lists a CC BY 4.0 license and asks users to cite the QLP
papers and the TESS mission when publishing results.
