# TOI-3505.01 program schedule record

`observing_schedule_2022-07-21.json` preserves the single TOI-3505.01 row
that Mason copied from the internship observing-schedule spreadsheet on
2026-07-24. The full sheet is not reproduced because the other rows are not
needed for this target and contain unrelated participant names.

The source row records:

- night: 2022-07-21;
- planned ingress and egress: 00:15 and 01:54;
- planned observing range: 21:10 to 04:55;
- R filter with 50-second exposures;
- orbital period: 2.9151488 days;
- Mason Cao as the student and Kevin / Schar Students as observers;
- the note `bad RA jumps at the beginning of the night`.

The row does not identify its time zone. The analysis uses
`America/New_York` as a documented working interpretation because the planned
21:10–04:55 interval then brackets the actual GMU sequence, which began at
22:13 EDT and ended at 03:18 EDT. The derived BJD_TDB values are never written
back into this source record. A UTC alternative is retained in the generated
timing output so the assumption remains visible.

The sky geometry strongly supports that interpretation. On 2022-07-21,
America/New_York was EDT (UTC-4). At the planned 21:10 EDT start, the Sun was
about 7.6 degrees below the horizon and TOI-3505 was about 36.9 degrees high.
At the planned 04:55 EDT end, the Sun was about 11.8 degrees below the horizon
and the target was about 35.4 degrees high. Reading the same clock cells as UTC
would place the start at 17:10 EDT, with the Sun about 36.8 degrees high and
the target about 8.0 degrees below the horizon. This makes UTC implausible for
an observing schedule, although the absent source timezone remains recorded.

This row recovers the historical scheduled window, but it is not a complete
historical ephemeris. The original workbook or URL, prediction epoch, timing
uncertainty, depth, and Transit Info file still need to be obtained from the
program.

Checks against the delivered data:

- the listed R filter and 50-second exposure match the FITS headers;
- the listed RA, 19:48:10, agrees at the row's one-second precision with the
  target coordinate used for the plate solution, 19:48:10.43;
- the alignment record contains several large early image-position steps,
  including a 53-pixel step by frame 13, consistent with the note about early
  pointing jumps; image-column motion alone is not relabeled as celestial RA;
- the listed declination is `#N/A`, so it contributes no coordinate check;
- the disposition code `VPC-+` is preserved without expanding it because the
  pasted sheet did not supply its definition.
