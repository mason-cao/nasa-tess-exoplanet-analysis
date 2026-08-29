"""Assemble the false-positive evidence for TOI-3505.01 into one record.

The program's paper guidance asks that the public follow-up archive be used to
gather and present high-contrast imaging, spectroscopy, and TESS light-curve
tests bearing on whether a candidate is a false positive. This script collects
the evidence this project already holds, adds two calculations that the evidence
supports, and states plainly what each scenario is and is not able to exclude.

The scenarios follow the TFOP SG1 vocabulary:

* NEB, a nearby eclipsing binary resolved from the target in follow-up imaging;
* BEB, an eclipsing binary blended in both the follow-up and TESS apertures,
  usually betrayed by a transit depth that varies with wavelength;
* EB, an eclipsing stellar companion on the target itself, too deep or too
  massive to be a planet.

This is an evidence summary, not a statistical validation. It computes no
false-positive probability, and it cannot address the unresolved companion.
"""

from __future__ import annotations

import json
import os
import tempfile
from pathlib import Path

os.environ.setdefault(
    "MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "nasa-tess-matplotlib-cache")
)

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from astropy import constants as const
from astropy import units as u


ROOT = Path(__file__).resolve().parents[1]
GROUND_CHECKS_PATH = ROOT / "outputs" / "toi3505_ground_checks" / "summary.json"
VALIDATION_PATH = ROOT / "outputs" / "toi3505_data_validation" / "analysis_summary.json"
EXOFOP_PATH = ROOT / "data" / "catalogs" / "toi3505" / "exofop_ground_followup.json"
DILUTION_PATH = ROOT / "outputs" / "toi3505_tess_pixels" / "dilution_screen.json"
EPHEMERIS_PATH = (
    ROOT / "outputs" / "toi3505_ephemeris_refined" / "ephemeris_refined.json"
)
OUTPUT_DIR = ROOT / "outputs" / "toi3505_false_positive"

# Effective wavelengths of the MuSCAT2 bands, in nanometres.  Used only to test
# whether the reported depths trend with wavelength.
BAND_WAVELENGTHS_NM = {"g": 477.0, "r": 623.0, "i": 762.0, "z_s": 870.0}
# Assumed host mass for the eclipsing-companion velocity bound.  The TIC gives
# no mass, so this follows the value adopted elsewhere in the project.
HOST_MASS_SOLAR = 1.25


def load_object(path: Path) -> dict[str, object]:
    """Load one JSON object and reject other top-level types."""
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError(f"Expected a JSON object in {path}")
    return payload


def chromatic_depth_test(depths_ppt: dict[str, float]) -> dict[str, object]:
    """Describe the reported multi-band depths against wavelength.

    A blended eclipsing binary generally produces a depth that varies with
    wavelength, because the blended source and the target differ in colour.
    The public report gives no per-band uncertainties, so this function records
    a slope and residual-based scale but does not treat their ratio as a
    calibrated significance or establish statistical achromaticity.
    """
    bands = [band for band in BAND_WAVELENGTHS_NM if band in depths_ppt]
    wavelength = np.array([BAND_WAVELENGTHS_NM[band] for band in bands], dtype=float)
    depth = np.array([float(depths_ppt[band]) for band in bands], dtype=float)
    mean_depth = float(np.mean(depth))
    scatter = float(np.std(depth, ddof=1))
    centred = wavelength - float(np.mean(wavelength))
    slope = float(np.sum(centred * (depth - mean_depth)) / np.sum(centred**2))
    # The public report gives no per-band uncertainty, so the residual scatter
    # about the fitted line supplies the slope error.  Using the scatter about
    # the mean instead would be circular: a real trend inflates that scatter and
    # would hide itself.
    residual = depth - (mean_depth + slope * centred)
    degrees_of_freedom = len(depth) - 2
    residual_scatter = (
        float(np.sqrt(np.sum(residual**2) / degrees_of_freedom))
        if degrees_of_freedom > 0
        else float("nan")
    )
    slope_scale = float(residual_scatter / np.sqrt(np.sum(centred**2)))
    monotonic = bool(np.all(np.diff(depth) > 0) or np.all(np.diff(depth) < 0))
    return {
        "bands": bands,
        "wavelengths_nm": wavelength.tolist(),
        "depths_ppt": depth.tolist(),
        "mean_depth_ppt": mean_depth,
        "depth_scatter_ppt": scatter,
        "residual_scatter_about_fit_ppt": residual_scatter,
        "slope_ppt_per_100nm": slope * 100.0,
        "residual_based_slope_scale_ppt_per_100nm": slope_scale * 100.0,
        "absolute_slope_to_residual_scale_ratio": (
            abs(slope) / slope_scale if slope_scale > 0 else float("nan")
        ),
        "monotonic_with_wavelength": monotonic,
        "no_apparent_monotonic_trend": not monotonic,
        "scale_interpretation": (
            "Descriptive residual-based scale only; the public report gives no "
            "per-band measurement uncertainties, so this is not a standard "
            "error or calibrated significance."
        ),
        "limits": [
            "The reporting team described this as a tentative egress, so the "
            "depths come from a partial event at modest significance.",
            "A blend whose colour closely matches the target would not produce a "
            "measurable trend and is not excluded by this test.",
            "Four reported depths without per-band errors cannot establish "
            "statistical achromaticity.",
        ],
    }


def velocity_amplitude_km_s(
    companion_mass_solar: float, period_days: float, host_mass_solar: float
) -> float:
    """Circular-orbit radial-velocity semi-amplitude for an edge-on orbit."""
    period = (period_days * u.day).to(u.s)
    companion = (companion_mass_solar * u.M_sun).to(u.kg)
    total = ((companion_mass_solar + host_mass_solar) * u.M_sun).to(u.kg)
    amplitude = (
        (2.0 * np.pi * const.G / period) ** (1.0 / 3.0)
        * companion
        / total ** (2.0 / 3.0)
    )
    return float(amplitude.to(u.km / u.s).value)


def eclipsing_companion_bound(
    period_days: float, observed_span_km_s: float
) -> dict[str, object]:
    """Compare the observed velocity span with stellar-companion expectations.

    A stellar companion producing the observed period would move the host by
    tens of kilometres per second. The reconnaissance velocities span far less
    than that, which is the cleanest argument against an eclipsing stellar
    companion on the target itself.
    """
    scenarios = {}
    for label, mass in (
        ("0.1 solar mass", 0.1),
        ("0.3 solar mass", 0.3),
        ("0.6 solar mass", 0.6),
    ):
        scenarios[label] = velocity_amplitude_km_s(mass, period_days, HOST_MASS_SOLAR)
    smallest = min(scenarios.values())
    return {
        "assumed_host_mass_solar": HOST_MASS_SOLAR,
        "period_days": period_days,
        "predicted_semi_amplitude_km_s": scenarios,
        "observed_velocity_span_km_s": observed_span_km_s,
        "smallest_stellar_scenario_km_s": smallest,
        "ratio_smallest_stellar_to_observed": smallest / observed_span_km_s,
        "stellar_companion_disfavoured": bool(smallest > 2.0 * observed_span_km_s),
        "limits": [
            "The observed span is the range of three reconnaissance velocities, "
            "not a fitted orbital semi-amplitude.",
            "The host mass is assumed, not measured.",
            "This bounds a companion to the star the spectrograph observed. It "
            "says nothing about the unresolved 0.517-arcsecond companion.",
        ],
    }


def build_assessment() -> dict[str, object]:
    """Collect every scenario into one auditable record."""
    ground = load_object(GROUND_CHECKS_PATH)
    validation = load_object(VALIDATION_PATH)
    exofop = load_object(EXOFOP_PATH)
    dilution = load_object(DILUTION_PATH)
    ephemeris = load_object(EPHEMERIS_PATH)

    nearby = ground["nearby_star_screen"]
    assert isinstance(nearby, dict)
    dispositions = nearby["disposition_counts"]
    assert isinstance(dispositions, dict)

    muscat2 = exofop["report_values"]["muscat2_2023_07_14"]  # type: ignore[index]
    recon = exofop["reconnaissance_spectroscopy_notes"]
    assert isinstance(recon, dict)
    adopted = ephemeris["ephemeris"]["adopted"]  # type: ignore[index]

    metrics = validation["official_metrics"]
    assert isinstance(metrics, list)
    odd_even = [float(row["odd_even_difference_sigma"]) for row in metrics]
    secondary = [float(row["weak_secondary_mes"]) for row in metrics]
    centroid = [float(row["centroid_to_tic_sigma"]) for row in metrics]

    chromatic = chromatic_depth_test(muscat2["depths_ppt"])  # type: ignore[index]
    velocity = eclipsing_companion_bound(
        float(adopted["period_days"]),
        float(recon["velocity_span_km_s"]),
    )

    return {
        "target": "TOI-3505.01",
        "analysis_type": "false-positive scenario evidence summary",
        "purpose": (
            "Gather the imaging, spectroscopic, and TESS light-curve evidence "
            "bearing on each standard false-positive scenario, as the program's "
            "paper guidance directs."
        ),
        "scenarios": {
            "nearby_eclipsing_binary": {
                "evidence": (
                    "Image-level screen of every catalogued neighbour within 2.5 "
                    "arcmin bright enough to mimic the signal, evaluated at the "
                    "2022 schedule window with the TFOP band correction applied."
                ),
                "measured_sources": int(
                    nearby["bright_enough_catalog_candidates"]  # type: ignore[arg-type]
                ),
                "disposition_counts": dispositions,
                "uncleared_with_eclipse_consistent_shape": int(
                    nearby["uncleared_with_eclipse_consistent_shape"]  # type: ignore[arg-type]
                ),
                "spoc_difference_image_offset_sigma": centroid,
                "spoc_difference_image_offset_note": (
                    "Offset between the difference-image centroid and the "
                    "catalogued target position, in sigma. Values below three "
                    "indicate the transit source coincides with the target at "
                    "TESS resolution."
                ),
                "public_neb_checks_on_file": [
                    "KeplerCam 2022-06-14",
                    "ULMT 2021-10-15",
                    "MuSCAT2 2023-07-14 (133 nearby stars ruled out)",
                ],
                "assessment": (
                    "Disfavoured. No measured neighbour shows an eclipse-shaped "
                    "event of the required depth, independent public NEB checks "
                    "report the same, and the SPOC difference-image centroid "
                    "stays within three sigma of the target in both sectors."
                ),
                "not_excluded": (
                    "Neighbours too faint for a decisive limit, and any source "
                    "inside the target aperture."
                ),
            },
            "blended_eclipsing_binary": {
                "evidence": (
                    "Wavelength dependence of the reported MuSCAT2 four-band depths."
                ),
                "chromatic_depth_test": chromatic,
                "assessment": (
                    "Weakly disfavoured. The reported depths show no apparent "
                    "monotonic trend with wavelength, but the absence of "
                    "per-band errors makes this a qualitative constraint."
                ),
                "not_excluded": (
                    "A blend whose colour matches the target, and any blend at "
                    "separations below the follow-up resolution."
                ),
            },
            "eclipsing_binary_on_target": {
                "evidence": (
                    "SPOC odd/even and weak-secondary diagnostics, plus the "
                    "amplitude of the public reconnaissance velocities."
                ),
                "odd_even_difference_sigma": odd_even,
                "weak_secondary_mes": secondary,
                "suspected_eclipsing_binary_flags": [
                    bool(row["suspected_eclipsing_binary"]) for row in metrics
                ],
                "eclipsing_companion_velocity_bound": velocity,
                "reported_companion_mass_note": (
                    "The public observing notes report the velocities as in phase "
                    "with the photometric ephemeris and consistent with a "
                    "companion of about 10 Jupiter masses, with scatter, on hold."
                ),
                "assessment": (
                    "Disfavoured. No significant odd/even difference, no secondary "
                    "eclipse detection, no SPOC eclipsing-binary flag, and a "
                    "velocity span far below any stellar companion."
                ),
                "not_excluded": (
                    "A definitive mass requires the multi-order velocity analysis "
                    "the public notes leave open."
                ),
            },
            "unresolved_close_companion": {
                "evidence": (
                    "SOAR speckle and Shane adaptive-optics imaging both resolve a "
                    "companion near 0.51 arcsec that no dataset in this work can "
                    "separate."
                ),
                "separation_arcsec": float(
                    dilution["unresolved_companion"]["separation_arcsec"]  # type: ignore[index]
                ),
                "flux_ratio_tess_band_proxy": float(
                    dilution["unresolved_companion"][  # type: ignore[index]
                        "flux_ratio_using_delta_i_as_tess_band_proxy"
                    ]
                ),
                "assessment": (
                    "Not addressed. This is the limiting scenario and the reason "
                    "no validation, dilution correction, or radius is claimed."
                ),
                "not_excluded": "Everything; the companion remains unresolved here.",
            },
        },
        "overall": (
            "The nearby and on-target eclipsing-binary scenarios are disfavoured "
            "by the available screens, while the blended-binary chromatic check "
            "is qualitative. The unresolved 0.517-arcsecond companion is not "
            "addressed, so this summary supports continued candidate-level "
            "follow-up without validating the object."
        ),
        "limits": [
            "An evidence summary, not a statistical validation. No false-positive "
            "probability is computed and no validation framework is run.",
            "External report values are used as published and are not re-reduced.",
            "The MuSCAT2 depths come from a tentative partial event and carry no "
            "published per-band uncertainties.",
            "The velocity bound assumes a host mass and treats a range of three "
            "reconnaissance velocities as a span, not a fitted amplitude.",
        ],
    }


def plot_assessment(result: dict[str, object], path: Path) -> None:
    """Draw the chromatic depth test and the velocity-amplitude comparison."""
    scenarios = result["scenarios"]
    assert isinstance(scenarios, dict)
    chromatic = scenarios["blended_eclipsing_binary"]["chromatic_depth_test"]  # type: ignore[index]
    velocity = scenarios["eclipsing_binary_on_target"][  # type: ignore[index]
        "eclipsing_companion_velocity_bound"
    ]

    figure, axes = plt.subplots(1, 2, figsize=(10.0, 4.0))

    wavelengths = np.array(chromatic["wavelengths_nm"], dtype=float)
    depths = np.array(chromatic["depths_ppt"], dtype=float)
    mean_depth = float(chromatic["mean_depth_ppt"])
    scatter = float(chromatic["depth_scatter_ppt"])
    axes[0].axhspan(
        mean_depth - scatter,
        mean_depth + scatter,
        color="#2C6B5F",
        alpha=0.15,
        label="mean ± scatter",
    )
    axes[0].axhline(mean_depth, color="#2C6B5F", linewidth=1.3)
    axes[0].plot(wavelengths, depths, "o", color="#1A1F1D", markersize=7, zorder=3)
    for wavelength, depth, band in zip(
        wavelengths, depths, chromatic["bands"], strict=True
    ):
        axes[0].annotate(
            band,
            (wavelength, depth),
            textcoords="offset points",
            xytext=(0, 10),
            ha="center",
            fontsize=9,
        )
    axes[0].set_ylim(1.9, 2.85)
    axes[0].set_xlabel("Effective wavelength (nm)")
    axes[0].set_ylabel("Reported depth (ppt)")
    axes[0].set_title("No apparent monotonic trend in reported depths", fontsize=10)
    axes[0].legend(fontsize=8, loc="lower right")

    predicted = velocity["predicted_semi_amplitude_km_s"]
    labels = list(predicted)
    values = [float(predicted[label]) for label in labels]
    positions = np.arange(len(labels))
    axes[1].barh(positions, values, color="#8A6110", alpha=0.85)
    axes[1].axvline(
        float(velocity["observed_velocity_span_km_s"]),
        color="#9E2F27",
        linewidth=1.8,
        label="observed velocity span",
    )
    axes[1].set_yticks(positions)
    axes[1].set_yticklabels(labels, fontsize=9)
    axes[1].set_xscale("log")
    # Matplotlib's default log minor labels collide at this range, so the ticks
    # are set explicitly.
    ticks = [1.0, 3.0, 10.0, 30.0, 100.0]
    axes[1].set_xticks(ticks)
    axes[1].set_xticklabels([f"{tick:g}" for tick in ticks], fontsize=9)
    axes[1].set_xticks([], minor=True)
    axes[1].set_xlim(1.5, 120.0)
    axes[1].set_xlabel("Velocity semi-amplitude (km s$^{-1}$)")
    axes[1].set_title("A stellar companion would move the host far more", fontsize=10)
    axes[1].legend(fontsize=8, loc="lower right")

    figure.tight_layout()
    figure.savefig(path, dpi=200)
    svg_path = path.with_suffix(".svg")
    figure.savefig(svg_path)
    # Matplotlib leaves trailing spaces in SVG path data.  Removing them keeps
    # generated artifacts clean under ``git diff --check`` without changing
    # the rendered figure.
    svg_path.write_text(
        "\n".join(line.rstrip() for line in svg_path.read_text().splitlines()) + "\n",
        encoding="utf-8",
    )
    plt.close(figure)


def write_readme(result: dict[str, object], path: Path) -> None:
    """Record the scenarios and what each cannot reach."""
    scenarios = result["scenarios"]
    assert isinstance(scenarios, dict)
    lines = [
        "# TOI-3505.01 false-positive evidence summary",
        "",
        "The program's paper guidance asks that the public follow-up archive be",
        "used to gather and present the imaging, spectroscopic, and TESS light-",
        "curve evidence bearing on whether a candidate is a false positive. This",
        "record collects that evidence for each standard scenario.",
        "",
        "It is an evidence summary, not a statistical validation. No false-positive",
        "probability is computed.",
        "",
        "## Scenarios",
        "",
    ]
    titles = {
        "nearby_eclipsing_binary": "Nearby eclipsing binary (NEB)",
        "blended_eclipsing_binary": "Blended eclipsing binary (BEB)",
        "eclipsing_binary_on_target": "Eclipsing binary on the target (EB)",
        "unresolved_close_companion": "Unresolved close companion",
    }
    for key, title in titles.items():
        entry = scenarios[key]
        assert isinstance(entry, dict)
        lines.extend(
            [
                f"### {title}",
                "",
                f"- Evidence: {entry['evidence']}",
                f"- Assessment: {entry['assessment']}",
                f"- Not excluded: {entry['not_excluded']}",
                "",
            ]
        )
    lines.extend(
        [
            "## Overall",
            "",
            str(result["overall"]),
            "",
            "## Limits",
            "",
        ]
    )
    limits = result["limits"]
    assert isinstance(limits, list)
    lines.extend(f"- {limit}" for limit in limits)
    lines.extend(
        [
            "",
            "## Products",
            "",
            "- `false_positive_assessment.json` - every scenario and calculation.",
            "- `01_false_positive_tests.png` and `.svg` - publication figure.",
            "",
            "Regenerate with:",
            "",
            "```bash",
            ".venv/bin/python src/assess_toi3505_false_positive.py",
            "```",
            "",
        ]
    )
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    result = build_assessment()
    (OUTPUT_DIR / "false_positive_assessment.json").write_text(
        json.dumps(result, indent=2) + "\n", encoding="utf-8"
    )
    plot_assessment(result, OUTPUT_DIR / "01_false_positive_tests.png")
    write_readme(result, OUTPUT_DIR / "README.md")

    scenarios = result["scenarios"]
    assert isinstance(scenarios, dict)
    chromatic = scenarios["blended_eclipsing_binary"]["chromatic_depth_test"]  # type: ignore[index]
    velocity = scenarios["eclipsing_binary_on_target"][  # type: ignore[index]
        "eclipsing_companion_velocity_bound"
    ]
    print("TOI-3505.01 false-positive evidence summary")
    print(
        "  descriptive chromatic depth slope: "
        f"{float(chromatic['slope_ppt_per_100nm']):+.3f} ppt per 100 nm; "
        "residual-based scale "
        f"{float(chromatic['residual_based_slope_scale_ppt_per_100nm']):.3f}"
    )
    print(
        "  smallest stellar scenario: "
        f"{float(velocity['smallest_stellar_scenario_km_s']):.1f} km/s vs observed "
        f"span {float(velocity['observed_velocity_span_km_s']):.2f} km/s"
    )
    print(f"Wrote {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
