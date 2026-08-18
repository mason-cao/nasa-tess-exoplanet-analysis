from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from toi3505_schedule import (
    DEFAULT_SCHEDULE_RECORD,
    analyze_schedule_window,
    fit_fixed_window,
    fixed_window_fraction,
    load_schedule_record,
    schedule_context,
    write_target_plot_config,
)


class ScheduleRecordTests(unittest.TestCase):
    def test_source_row_preserves_key_cells(self) -> None:
        record = load_schedule_record()
        row = record["row"]

        self.assertEqual(row["Target"], "TOI 3505.01")
        self.assertEqual(row["Ingress"], "0:15")
        self.assertEqual(row["Egress"], "1:54")
        self.assertEqual(row["Exp"], "50s")
        self.assertEqual(row["Orbital Period"], "2.9151488")
        self.assertEqual(
            row["Notes"], "bad RA jumps at the beginning of the night"
        )
        self.assertEqual(row["Apertures File"], "")
        self.assertNotIn("Finder Chart", row)
        self.assertEqual(
            record["timezone_confirmation"]["iana_timezone"],
            "America/New_York",
        )

    def test_confirmed_eastern_times_bracket_data_but_utc_does_not(self) -> None:
        context = schedule_context(
            observation_start_bjd=2459782.598234811,
            observation_end_bjd=2459782.809458706,
        )
        comparison = context["observation_comparison"]
        working = context["working_interpretation"]

        self.assertTrue(working["timezone_confirmed"])
        self.assertTrue(comparison["working_planned_range_brackets_observation"])
        self.assertTrue(comparison["working_event_fully_covered"])
        self.assertFalse(comparison["utc_event_fully_covered"])
        self.assertAlmostEqual(
            working["times"]["ingress"]["bjd_tdb"], 2459782.682400599
        )
        self.assertAlmostEqual(
            working["times"]["egress"]["bjd_tdb"], 2459782.751150993
        )

    def test_confirmed_eastern_times_match_nighttime_geometry(self) -> None:
        context = schedule_context()
        working = context["working_interpretation"]
        plausibility = context["timezone_plausibility"]

        self.assertEqual(working["timezone_abbreviation"], "EDT")
        self.assertEqual(working["utc_offset_hours"], -4.0)
        self.assertEqual(
            plausibility["confirmed_interpretation"], "America/New_York"
        )
        self.assertTrue(plausibility["working_start_is_observable"])
        self.assertTrue(plausibility["working_end_is_observable"])
        self.assertFalse(plausibility["utc_start_is_observable"])

        working_start = plausibility["working_interpretation_altitudes"][
            "planned_start"
        ]
        utc_start = plausibility["utc_alternative_altitudes"]["planned_start"]
        self.assertLess(working_start["sun_altitude_degrees"], 0.0)
        self.assertGreater(working_start["target_altitude_degrees"], 20.0)
        self.assertGreater(utc_start["sun_altitude_degrees"], 20.0)
        self.assertLess(utc_start["target_altitude_degrees"], 0.0)

    def test_fixed_window_fit_recovers_synthetic_depth(self) -> None:
        bjd = np.linspace(2459782.60, 2459782.80, 240)
        ingress = 2459782.68
        egress = 2459782.74
        fraction = fixed_window_fraction(bjd, ingress, egress, 50.0)
        hours = (bjd - np.mean(bjd)) * 24.0
        flux = 1.0 + 0.0001 * hours - 0.004 * fraction
        flux += 0.00015 * np.sin(np.arange(len(bjd)) * 0.37)
        result = fit_fixed_window(
            bjd,
            flux,
            np.full(len(bjd), 0.001),
            np.ones(len(bjd), dtype=bool),
            ingress,
            egress,
        )

        self.assertAlmostEqual(result["depth"], 0.004, delta=0.00005)
        self.assertGreater(result["in_window_points"], 60)

    def test_actual_historical_window_has_no_transit_like_dimming(self) -> None:
        curve = pd.read_csv(
            ROOT
            / "outputs"
            / "toi3505_final_candidate"
            / "TOI_3505.01_2022-07-22_R_light_curve.csv"
        )
        result = analyze_schedule_window(
            DEFAULT_SCHEDULE_RECORD,
            curve["BJD_TDB"].to_numpy(dtype=float),
            curve["Relative_Brightness"].to_numpy(dtype=float),
            curve["Flux_Error"].to_numpy(dtype=float),
            curve["Used_in_Plot"].to_numpy(dtype=bool),
        )
        fixed = result["fixed_window_check"]

        self.assertEqual(fixed["in_window_points"], 90)
        self.assertFalse(fixed["transit_like_dimming_above_3_sigma"])
        self.assertLess(fixed["observed_depth_ppt"], 0.0)
        self.assertGreater(fixed["injected_total_depth_snr"], 5.0)

    def test_plot_config_uses_historical_markers(self) -> None:
        keys = (
            ".plot.title",
            ".plot.subtitle",
            ".plot.xlabel",
            ".plot.showVMarker1",
            ".plot.showVMarker2",
            ".plot.vMarker1TopText",
            ".plot.vMarker1BotText",
            ".plot.vMarker1Value",
            ".plot.vMarker2TopText",
            ".plot.vMarker2BotText",
            ".plot.vMarker2Value",
            ".plot.useInEgressMarkers",
            ".plot.ingressTime",
            ".plot.egressTime",
            ".plot.xMin",
            ".plot.xMax",
            ".plot.orbitalPeriod0",
            ".plot.teff0",
            ".plot.priorCenter[0][1]",
            ".plot.priorCenter[0][3]",
            ".plot.useTransitFit0",
            ".plot.detrendFit0",
            ".plot.orbitalPeriod1",
            ".plot.teff1",
            ".plot.priorCenter[1][1]",
            ".plot.priorCenter[1][3]",
            ".plot.useTransitFit1",
            ".plot.detrendFit1",
        )
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory) / "source.plotcfg"
            destination = Path(directory) / "target.plotcfg"
            source.write_bytes(
                ("\r\n".join(f"{key}=old" for key in keys) + "\r\n").encode(
                    "ascii"
                )
            )
            write_target_plot_config(
                source,
                destination,
                schedule_context(),
                observation_start_bjd=2459782.598234811,
                observation_end_bjd=2459782.809458706,
                orbital_period_days=2.9151556,
                transit_center_bjd=2459781.8737626,
                catalog_depth_ppt=2.910,
                target_teff_k=6220.0,
            )
            result = destination.read_bytes()

        self.assertIn(b".plot.showVMarker1=true", result)
        self.assertIn(b".plot.title=TOI-3505.01, UT 2022-07-22", result)
        self.assertIn(
            b".plot.subtitle=GMU 0.8 m (R filter, 50 s exp, fap 25-70-139)",
            result,
        )
        self.assertIn(b".plot.vMarker1TopText=2022 schedule (EDT)", result)
        self.assertIn(b".plot.vMarker1Value=0.682400599", result)
        self.assertIn(b".plot.vMarker2Value=0.751150993", result)
        self.assertIn(b".plot.useInEgressMarkers=false", result)
        self.assertIn(b".plot.orbitalPeriod0=2.9151556", result)
        self.assertIn(b".plot.teff0=6220.0", result)
        self.assertIn(b".plot.priorCenter[0][1]=0.002910000", result)
        self.assertIn(b".plot.priorCenter[0][3]=2459781.873762600", result)
        self.assertIn(b".plot.useTransitFit0=false", result)
        self.assertIn(b".plot.detrendFit0=false", result)
        self.assertIn(b".plot.orbitalPeriod1=2.9151556", result)
        self.assertIn(b".plot.teff1=6220.0", result)
        self.assertIn(b".plot.priorCenter[1][1]=0.002910000", result)
        self.assertIn(b".plot.priorCenter[1][3]=2459781.873762600", result)
        self.assertIn(b".plot.useTransitFit1=false", result)
        self.assertIn(b".plot.detrendFit1=false", result)
        self.assertEqual(result.count(b"\r\n"), len(keys))


if __name__ == "__main__":
    unittest.main()
