# TOI-3505.01 plate solution

I plate-solved the first reduced image using the settings in the Schar
light-curve tutorial:

- Plate scale search: 0.300-0.400 arcseconds per pixel
- Position search radius: 20 arcminutes
- Sky model order: 2
- Time keyword kept in the image: BJD_TDB

The solution found a plate scale of 0.362 arcseconds per pixel. The target
coordinates land on the measured star in the image, and the full image was not
uploaded. The local source list contained 114 stars.

The plate-solved FITS file is kept on this computer because it is too large for
the repository. `solution.json`, `source_list.csv`, and `wcs_header.txt` contain
the small results needed to check or repeat this step.
