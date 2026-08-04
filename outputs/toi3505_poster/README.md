# TOI-3505.01 symposium poster

Deliverables for the Schar Scholars symposium. Both boards are 48 x 36 inches
in landscape, matching the program template, and carry the same science.

| | v1 | v2 |
| --- | --- | --- |
| Intended as | The full, scientifically detailed board | The simpler board |
| Figures | Analysis figures made for this poster | The AstroImageJ products the rest of the cohort shows |
| Left column | Abstract, Introduction, Methods, False-alarm checks | Abstract, Methods |
| Middle | Results: 5 figures, comparison table | Results: 5 figures |
| Right column | System Parameters, Discussion, Conclusions | Discussion, Checks, Conclusions and Future Work |
| Footer | Future Work, Citations, Acknowledgements | Citations, Acknowledgements |

## Files to submit or print

| File | Use |
| --- | --- |
| `TOI-3505.01_Mason_Cao_poster.pdf` | v1 print master. One page, exactly 48 x 36 inches. |
| `TOI-3505.01_Mason_Cao_poster.png` | v1 raster backup, 6400 x 4800, about 133 dots per inch at full size. |
| `TOI-3505.01_Mason_Cao_poster.html` | v1 source. Self-contained, every image inlined. |
| `TOI-3505.01_Mason_Cao_poster_v2.pdf` | v2 print master. |
| `TOI-3505.01_Mason_Cao_poster_v2.png` | v2 raster backup. |
| `TOI-3505.01_Mason_Cao_poster_v2.html` | v2 source. |
| `v1/` | Frozen copy of the earlier v1, kept before the fork was written. |

## Figures

v1:

| File | Poster figure |
| --- | --- |
| `04_field_and_aperture.png` | Fig. 1, plate-solved field with the target, C2-C11, and the adopted aperture |
| `05_ground_light_curve.png` | Fig. 2, final 263-point light curve from the Mason 0.8 m |
| `01_transit_timing.png` | Fig. 3, observed-minus-calculated timing diagram |
| `06_phase_folded.png` | Fig. 4, four TESS sectors folded on the new period |
| `02_nearby_star_screen.png` | Fig. 5, delta magnitude against depth reached, all 280 stars |
| `03_depth_consistency.png` | Depth and odd-even comparison; on neither board, kept as backup for questions |

v2:

| File | Poster figure |
| --- | --- |
| `v2_01_seeing_profile.png` | Fig. 1, seeing profile at the adopted 25 px aperture |
| `v2_02_comparison_field.png` | Fig. 2, AstroImageJ frame with T1 and C2-C11 |
| `v2_04_neb_field.png` | Fig. 3, the 2.5 arcminute nearby-star screen on our own image |
| `../toi3505_review_package/01_TOI_3505.01_final_light_curve.png` | Fig. 4, the finished light curve |
| `v2_03_dmag_rms.png` | Fig. 5, delta magnitude against depth reached |

Three v2 figures are regenerated rather than reused. The saved AstroImageJ
seeing profile in `outputs/toi3505_seeing` was taken during the 35-pixel trial,
so it contradicts the 25-pixel aperture the analysis adopted; it is rebuilt from
the plate-solved image at the adopted settings. AstroImageJ's own nearby-star
plot was never exported, so the delta-magnitude figure and the 2.5 arcminute
field are drawn from the measured neighbour table.

Figures are inlined as PNG, not JPEG. The program's paper-writing lecture asks
for lossless figures, and for these panels PNG is also the smaller encoding in
most cases.

`poster_analysis.json` holds every derived number both boards quote.
`v2_figure_summary.json` records the v2 figure inputs.

## Regenerate

```bash
cd src
../.venv/bin/python refine_toi3505_ephemeris.py       # adopted ephemeris, must run first
../.venv/bin/python make_toi3505_ground_checks.py     # 2.5 arcmin nearby-star screen
../.venv/bin/python make_toi3505_poster_analysis.py   # v1 figures + poster_analysis.json
../.venv/bin/python make_toi3505_poster_aij_figures.py  # v2 figures
../.venv/bin/python build_toi3505_poster.py v1        # inlines images into the v1 HTML
../.venv/bin/python build_toi3505_poster.py v2        # same for v2
```

Export a print file from either HTML:

```bash
"/Applications/Google Chrome.app/Contents/MacOS/Google Chrome" --headless --disable-gpu \
  --print-to-pdf="outputs/toi3505_poster/TOI-3505.01_Mason_Cao_poster.pdf" \
  --no-pdf-header-footer --virtual-time-budget=25000 \
  "file://$PWD/outputs/toi3505_poster/TOI-3505.01_Mason_Cao_poster.html"
```

Both layouts are authored at 1600 x 1200 px, 33.3 px per inch. The print rule
scales by 2.88 so the page comes out at a real 48 x 36 inches. Body text is
about 21 pt on v1 and 22 pt on v2 at full size.

## What the posters claim

The headline is the four-sector linear ephemeris,
P = 2.9151516 +/- 0.0000049 d from 27 mid-transit times spanning 5.05 years.
It is 2.40 times tighter than the TESS Object of Interest catalog period and
1.43 times tighter than the SPOC multi-sector fit, and it agrees with both
(0.32 and 0.71 sigma).

Both boards state the two comparisons side by side and attribute the gain to
the 5.05-year baseline rather than to a better single measurement, because the
catalog row rests on three sectors and the SPOC fit on two. The earlier wording,
"twice as precise as the catalog value", was true against the catalog and
misleading against SPOC, and has been replaced.

See `../toi3505_ephemeris_refined/README.md` for how that ephemeris is built and
how it differs from the frozen four-sector fit.

## Prior work on this target

TOI-3505.01 has no publication, but ExoFOP holds unpublished follow-up that
neither board analyses. Both say so in the Discussion, and the future-work
bullets point at it rather than proposing it as new:

| Kind | What is on file |
| --- | --- |
| Ground photometry | 6 light curves: GMU 0.8 m R 2021-06-28, ULMT rp 2021-10-15, KeplerCam ip 2022-06-14, CMO-SAI g' 2023-05-05, OAA Ic 2023-07-11, TCS-MuSCAT2 g/r/i/z_s 2023-07-14 |
| Imaging | SOAR HRCam speckle 2021-10-01 (0.517", dI = 1.7); Shane ShARCS AO 2021-07-19 (0.51", dJ = 1.53, dKs = 1.59, plus five sources at 6-8") |
| Spectroscopy | TRES x3 (2021), Keck HIRES x2 (2022), including 6 radial velocities |

TFOPWG disposition is still PC.

TESS coverage is exhausted: a MAST query returns sectors 14, 41, 54, and 81 and
nothing later, so the 5.05-year baseline is the longest currently available.

Both boards keep the project's standing limits: no ground transit detection,
SPOC comparisons described as method checks rather than independent
confirmation, and no claim that ground data or TESS resolves the
0.517-arcsecond companion.
