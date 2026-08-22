# TOI-3505.01 schedule reconstruction

The recovered 2022 schedule window is consistent with a superseded public
ephemeris. Propagating the ExoFOP ULMT prediction from 2021-10-15 by
`96` cycles at `2.9174250` days reproduces the
schedule ingress, midpoint, and egress to within
`56.1` seconds.

The schedule row visibly contains the later `2.9151488`-day period. Applying
that period to the same 2021 midpoint instead predicts an event
`5.24`
hours before the schedule midpoint. The adopted four-sector TESS ephemeris
places the nearest event
`20.31` hours before
the schedule midpoint.

## Interpretation

The timing cells are therefore strongly consistent with a prediction retained
from the older ephemeris while the displayed period was updated. This is an
archival reconstruction, not proof: the original workbook, formulas, revision
history, and Transit Info file remain unavailable. The result explains why the
2022 observing window was stale; it does not convert the ground-based null into
a transit detection.

## Products

- `schedule_reconstruction.json` - complete inputs, calculations, and limits.
- `marker_comparison.csv` - marker-by-marker propagation residuals.
- `01_schedule_reconstruction.png` and `.svg` - publication figure.

Regenerate with:

```bash
.venv/bin/python src/reconstruct_toi3505_schedule.py
```
