# TOI-3505.01 observation plan

Planning windows generated from
`outputs/toi3505_ephemeris_refined/ephemeris_refined.json`. These times are not an approved Transit Info file; independently
confirm the ephemeris, observatory clock, and telescope availability before
observing.

## Observable windows

| Cycle | Sequence start | Midpoint | Sequence end | 1σ (min) | Min altitude | Moon below throughout |
| ---: | --- | --- | --- | ---: | ---: | :---: |
| 505 | 2026-08-12 21:54 EDT | 2026-08-13 00:16 EDT | 2026-08-13 02:37 EDT | 3.35 | 45.3 | yes |

The sequence includes 1 hour of baseline on each side of a
2.711-hour transit. A complete window must keep the target
above 25 degrees and the Sun below
-18 degrees.

## Files

- `transit_windows.csv` — every predicted midpoint in the requested range,
  including rejected windows and reasons.
- `observation_plan.json` — inputs, thresholds, adopted ephemeris, and events.
- `observable_transits.ics` — calendar events for windows that pass the limits.

Regenerate with:

```bash
.venv/bin/python src/plan_toi3505_observations.py
```
