from __future__ import annotations

import sys
import unittest
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from toi3505_tess_tools import (
    LightCurveData,
    event_cycles,
    grid_box_fit,
    integrated_box_fraction,
    phase_offset,
    weighted_linear_ephemeris,
)


class TessToolTests(unittest.TestCase):
    def test_phase_and_cycle_counting_at_boundaries(self) -> None:
        period = 2.0
        epoch = 100.0
        times = np.array([96.0, 98.0, 99.9, 100.0, 100.9, 102.0, 104.0])

        np.testing.assert_array_equal(
            event_cycles(times, period, epoch), np.array([-2, -1, 0, 0, 0, 1, 2])
        )
        np.testing.assert_allclose(
            phase_offset(times, period, epoch),
            np.array([0.0, 0.0, -0.1, 0.0, 0.9, 0.0, 0.0]),
        )

    def test_integrated_box_accounts_for_long_exposure(self) -> None:
        fraction = integrated_box_fraction(
            np.array([0.0, 0.05, 0.10]),
            duration_days=0.10,
            cadence_days=0.04,
        )

        np.testing.assert_allclose(fraction, np.array([1.0, 0.5, 0.0]))

    def test_grid_fit_recovers_known_box(self) -> None:
        period = 2.0
        epoch = 2459000.0
        cadence = 10.0 / 1440.0
        time = np.arange(epoch - 5.0, epoch + 5.0, cadence)
        duration = 0.09
        offset = 0.012
        depth = 0.004
        cycles = event_cycles(time, period, epoch)
        phase = phase_offset(time, period, epoch + offset)
        exposure_fraction = integrated_box_fraction(phase, duration, cadence)
        local_slope = 0.0003 * phase_offset(time, period, epoch)
        flux = 1.0 + 0.0001 * cycles + local_slope - depth * exposure_fraction
        error = np.full(len(time), 0.0005)
        curve = LightCurveData(
            path=Path("synthetic.fits"),
            sector=1,
            pipeline="QLP",
            flux_name="KSPSAP_FLUX",
            time_bjd=time,
            flux=flux,
            flux_error=error,
            quality=np.zeros(len(time), dtype=int),
            cadence_days=cadence,
            crowdsap=None,
            flfrcsap=None,
        )

        fit, _, _ = grid_box_fit(
            curve,
            period_days=period,
            epoch_bjd=epoch,
            durations_days=np.linspace(0.07, 0.11, 17),
            offsets_days=np.linspace(-0.01, 0.03, 21),
        )

        self.assertAlmostEqual(fit.depth, depth, places=5)
        self.assertAlmostEqual(fit.duration_days, duration, places=5)
        self.assertAlmostEqual(fit.time_offset_days, offset, places=5)

    def test_weighted_ephemeris_recovers_period(self) -> None:
        cycles = np.array([-10, -2, 0, 7, 15])
        epoch = 2459000.123
        period = 2.915
        times = epoch + period * cycles
        errors = np.full(len(cycles), 0.001)

        result = weighted_linear_ephemeris(cycles, times, errors)

        self.assertAlmostEqual(result["epoch_bjd"], epoch, places=8)
        self.assertAlmostEqual(result["period_days"], period, places=8)
        self.assertEqual(result["events"], len(cycles))


if __name__ == "__main__":
    unittest.main()
