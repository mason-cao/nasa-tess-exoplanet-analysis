# Representative TOI-3505.01 plate-solve plan

Five frames are frozen in `frame_plan.csv`: early, middle, late, and one on
each side of the minimum-airmass point. Current mode: `plan_only`.

The first frame already has a valid plate solution. The additional independent
solutions require `ASTROMETRY_NET_API_KEY` and can be run with:

```bash
.venv/bin/python src/plate_solve_toi3505_representative.py --run
```

Only measured source positions are submitted. The full science images and API
key are not uploaded or written by the scripts.
