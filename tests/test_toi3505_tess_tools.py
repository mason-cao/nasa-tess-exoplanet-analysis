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
    flat_fraction_from_geometry,
    grid_box_fit,
    integrated_box_fraction,
    integrated_transit_fraction,
    phase_offset,
    weighted_linear_ephemeris,
)


class TransitShapeTests(unittest.TestCase):
    """The trapezoid has to reduce to the box and conserve area."""

    def test_flat_fraction_of_one_reproduces_the_box(self) -> None:
        phase = np.linspace(-0.15, 0.15, 61)
        np.testing.assert_allclose(
            integrated_transit_fraction(phase, 0.11, 0.02, flat_fraction=1.0),
            integrated_box_fraction(phase, 0.11, 0.02),
            atol=1e-12,
        )

    def test_trapezoid_area_matches_the_analytic_value(self) -> None:
        # Summing the exposure-averaged depth over contiguous exposures
        # recovers the area of the trapezoid, (T14 + T23) / 2.
        duration, cadence, flat = 0.12, 0.001, 0.4
        centers = np.arange(-0.2, 0.2, cadence) + cadence / 2.0
        area = integrated_transit_fraction(
            centers, duration, cadence, flat_fraction=flat
        ).sum() * cadence
        expected = 0.5 * (duration + duration * flat)
        self.assertAlmostEqual(area, expected, places=6)

    def test_fully_grazing_shape_is_triangular(self) -> None:
        # With no flat bottom the profile peaks at half depth for an exposure
        # centered on mid-transit only in the limit of a long exposure; for a
        # short one it reaches full depth.
        depth_short = integrated_transit_fraction(
            np.array([0.0]), 0.10, 1e-4, flat_fraction=0.0
        )[0]
        self.assertAlmostEqual(depth_short, 1.0, places=3)
        half_way = integrated_transit_fraction(
            np.array([0.025]), 0.10, 1e-4, flat_fraction=0.0
        )[0]
        self.assertAlmostEqual(half_way, 0.5, places=3)

    def test_trapezoid_is_shallower_at_the_edges_than_a_box(self) -> None:
        edge = np.array([0.04])
        box = integrated_box_fraction(edge, 0.10, 1e-4)[0]
        trapezoid = integrated_transit_fraction(edge, 0.10, 1e-4, flat_fraction=0.4)[0]
        self.assertAlmostEqual(box, 1.0, places=3)
        self.assertLess(trapezoid, box)

    def test_flat_fraction_from_geometry(self) -> None:
        # A central transit reduces to the analytic (1 - k) / (1 + k).
        self.assertAlmostEqual(
            flat_fraction_from_geometry(0.05, 0.0), 0.95 / 1.05, places=12
        )
        # Once the planet's disk crosses the stellar limb there is no flat part.
        self.assertEqual(flat_fraction_from_geometry(0.06, 0.99), 0.0)
        # The TOI-3505.01 official geometry sits between the two.
        ratio = flat_fraction_from_geometry(0.061769367661928656, 0.9159642085456406)
        self.assertAlmostEqual(ratio, 0.37838, places=5)


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
