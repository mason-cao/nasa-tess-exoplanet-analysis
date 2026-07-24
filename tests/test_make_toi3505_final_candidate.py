from __future__ import annotations

import sys
import unittest
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from make_toi3505_final_candidate import (
    apply_condition_model,
    bin_light_curve,
    blocked_model_review,
    fit_robust_linear,
    select_promoted_model,
)


class FinalCandidateTests(unittest.TestCase):
    def test_robust_fit_resists_single_large_outlier(self) -> None:
        x = np.linspace(-2.0, 2.0, 101)
        y = 1.0 + 0.02 * x
        y[50] += 1.0
        features = pd.DataFrame({"airmass": x})

        model = fit_robust_linear(features, y, ("airmass",))
        predicted = model.predict(pd.DataFrame({"airmass": [-2.0, 2.0]}))

        self.assertAlmostEqual(float(predicted[0]), 0.96, places=3)
        self.assertAlmostEqual(float(predicted[1]), 1.04, places=3)

    def test_blocked_cv_promotes_strong_repeatable_condition_term(self) -> None:
        rng = np.random.default_rng(42)
        points = 240
        condition = np.linspace(-1.0, 1.0, points)
        raw = 1.0 + 0.025 * condition + rng.normal(0.0, 0.001, points)
        features = pd.DataFrame({"airmass": condition})
        eligible = np.ones(points, dtype=bool)

        review, _ = blocked_model_review(
            features,
            raw,
            eligible,
            candidate_models={"none": (), "airmass": ("airmass",)},
        )

        self.assertEqual(select_promoted_model(review), "airmass")
        row = review.loc[review["model"] == "airmass"].iloc[0]
        self.assertTrue(bool(row["passes_absolute_rules"]))

    def test_blocked_cv_keeps_raw_curve_when_correction_fails_rules(self) -> None:
        points = 240
        raw = 1.0 + 0.001 * np.sin(np.arange(points) * 0.7)
        features = pd.DataFrame({"airmass": np.tile([-1.0, 1.0], points // 2)})
        eligible = np.ones(points, dtype=bool)

        review, _ = blocked_model_review(
            features,
            raw,
            eligible,
            candidate_models={"none": (), "airmass": ("airmass",)},
        )

        self.assertEqual(select_promoted_model(review), "none")

    def test_no_model_returns_normalized_raw_curve(self) -> None:
        raw = np.array([0.99, 1.00, 1.01])
        features = pd.DataFrame(index=np.arange(3))
        eligible = np.ones(3, dtype=bool)

        corrected, baseline, model = apply_condition_model(
            features, raw, eligible, ()
        )

        np.testing.assert_allclose(corrected, raw)
        np.testing.assert_allclose(baseline, np.ones(3))
        self.assertEqual(model.feature_names, ())

    def test_fixed_width_bins_preserve_point_count(self) -> None:
        hours = np.array([0.00, 0.05, 0.20, 0.22])
        values = np.array([0.99, 1.01, 1.02, 0.98])

        binned = bin_light_curve(hours, values, bin_minutes=10.0)

        self.assertEqual(int(binned["points"].sum()), len(hours))
        self.assertEqual(len(binned), 2)


if __name__ == "__main__":
    unittest.main()
