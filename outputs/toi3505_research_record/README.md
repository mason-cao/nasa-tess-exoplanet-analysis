# TOI-3505.01 research record

This folder freezes the reproducibility record for the current analysis.

- `812` files are listed with size, modification time, and SHA-256.
- The six original archives are included.
- Reduced and aligned images are included.
- Software and system versions are in `software_versions.json`.
- Ground and TESS choices are in `frozen_analysis_config.json`.
- The adopted refined ephemeris and observation-planning limits are frozen in
  that configuration from their canonical JSON outputs.
- The manuscript source, self-contained HTML, PDF, and machine-readable paper
  values are included in the manifest and frozen configuration.
- Claims, decisions, and shared-data dependencies have separate CSV ledgers.

Regenerate the full record with:

```bash
.venv/bin/python src/make_toi3505_research_record.py
```

The manifest records files; it does not make unresolved provenance or
scientific limitations disappear. Those limits remain in the claim ledger.
