# TOI-3505.01 target and aperture check

I plate-solved the first reduced image and used the target coordinates in its
header to locate TOI-3505.01. The file I used is `TOI_3505.01_50.000s_R-0001_wcs.fits`. The plate
solution puts the target at FITS pixel
(1850.72, 1742.94). AstroImageJ centered the star at
(1850.58, 1742.82), a difference of 0.18
pixel. This confirms that the Seeing Profile was measured on the target.

The AstroImageJ Seeing Profile measured a FWHM of 13.64 pixels. Its
starting photometry settings are:

- Source radius: 35 pixels
- Background inner radius: 70 pixels
- Background outer radius: 139 pixels

The plate scale is 0.362 arcseconds per pixel. A nearby star is
137.20 pixels from the target, so it falls inside the
139-pixel outer background ring. I will also test a
70-100 pixel background area,
which stays one source radius away from that star. These are starting values
for Multi-Aperture. I still need to choose the final aperture by checking the
seeing and comparison-star trends across the full set of usable images.
