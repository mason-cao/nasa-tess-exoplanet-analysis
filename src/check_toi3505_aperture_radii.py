"""Independently check TOI-3505.01 source-aperture radii and BJD_TDB.

The AstroImageJ result currently uses a 35-pixel source radius and a 70-139
pixel sky annulus.  This script remeasures all 281 aligned images with
Photutils at seven radii, using AstroImageJ centroids but an independent
sigma-clipped annulus estimator.  It is a robustness check, not a claim that
the two implementations are independent observations.
"""

from __future__ import annotations

import argparse
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
import pandas as pd
from astropy import units as u
from astropy.coordinates import EarthLocation, SkyCoord
from astropy.io import fits
from astropy.stats import SigmaClip
from astropy.time import Time, TimeDelta
from astropy.utils import iers
from photutils.aperture import (
    ApertureStats,
    CircularAnnulus,
    CircularAperture,
    aperture_photometry,
)

from analyze_toi3505_photometry import (
    COMPARISON_STARS,
    DEFAULT_WIDE_TABLE,
    clipped_standard_deviation,
    correlation,
    load_table,
    normalized,
    robust_scatter,
)
from make_toi3505_final_candidate import robust_sigma


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_ALIGNED_DIR = ROOT / "data" / "ground" / "toi3505" / "aligned"
DEFAULT_OUTPUT_DIR = ROOT / "outputs" / "toi3505_aperture_check"
SOURCE_RADII = (20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0)
PRIMARY_RADIUS = 35.0
ANNULUS_INNER = 70.0
ANNULUS_OUTER = 139.0
PRIMARY_TOLERANCE = 0.10
STAR_NAMES = ("T1", *COMPARISON_STARS)
OBSERVATORY = EarthLocation.from_geodetic(
    lon=-77.3053299972 * u.deg,
    lat=38.82817 * u.deg,
    height=154.0 * u.m,
)
TARGET_COORD = SkyCoord(
    ra=19.802897222222224 * u.hourangle,
    dec=18.69891388888889 * u.deg,
    frame="icrs",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--aligned-dir", type=Path, default=DEFAULT_ALIGNED_DIR)
    parser.add_argument("--table", type=Path, default=DEFAULT_WIDE_TABLE)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    return parser.parse_args()


def measure_frame(
    data: np.ndarray,
    positions: np.ndarray,
    radii: tuple[float, ...] = SOURCE_RADII,
    *,
    annulus_inner: float = ANNULUS_INNER,
    annulus_outer: float = ANNULUS_OUTER,
) -> tuple[np.ndarray, np.ndarray]:
    """Return background-subtracted flux[radius, star] and local background."""
    positions = np.asarray(positions, dtype=float)
    annulus = CircularAnnulus(
        positions,
        r_in=annulus_inner,
        r_out=annulus_outer,
    )
    sigma_clip = SigmaClip(
        sigma=3.0,
        maxiters=5,
        cenfunc="median",
        stdfunc="mad_std",
    )
    annulus_stats = ApertureStats(data, annulus, sigma_clip=sigma_clip)
    background = np.asarray(annulus_stats.median, dtype=float)
    apertures = [CircularAperture(positions, r=radius) for radius in radii]
    photometry = aperture_photometry(data, apertures, method="exact")
    flux = np.empty((len(radii), len(positions)), dtype=float)
    for radius_index, aperture in enumerate(apertures):
        column = f"aperture_sum_{radius_index}"
        aperture_sum = np.asarray(photometry[column], dtype=float)
        flux[radius_index] = aperture_sum - background * float(aperture.area)
    return flux, background


def differential_curve(flux: np.ndarray) -> np.ndarray:
    """Form T1 divided by the sum of all ten comparison stars."""
    flux = np.asarray(flux, dtype=float)
    comparison_sum = np.sum(flux[:, 1:], axis=1)
    curve = np.full(len(flux), np.nan, dtype=float)
    valid = (
        np.isfinite(flux[:, 0])
        & np.isfinite(comparison_sum)
        & (flux[:, 0] > 0)
        & (comparison_sum > 0)
    )
    curve[valid] = flux[valid, 0] / comparison_sum[valid]
    return normalized(curve)


def median_pseudo_target_scatter(flux: np.ndarray) -> float:
    """Summarize comparison stability without using the target flux."""
    flux = np.asarray(flux, dtype=float)
    comparison = flux[:, 1:]
    scatters: list[float] = []
    for index in range(comparison.shape[1]):
        others = np.sum(np.delete(comparison, index, axis=1), axis=1)
        ratio = comparison[:, index] / others
        scatters.append(robust_scatter(ratio))
    return float(np.median(scatters))


def primary_radius_passes(metrics: pd.DataFrame) -> bool:
    """Require the frozen radius to sit on both precision plateaus."""
    primary = metrics.loc[np.isclose(metrics["source_radius_pixels"], PRIMARY_RADIUS)]
    if len(primary) != 1:
        raise ValueError("The metrics must contain exactly one 35-pixel row")
    target_best = float(metrics["target_robust_scatter_ppt"].min())
    pseudo_best = float(metrics["median_pseudo_target_scatter_ppt"].min())
    target_value = float(primary.iloc[0]["target_robust_scatter_ppt"])
    pseudo_value = float(primary.iloc[0]["median_pseudo_target_scatter_ppt"])
    return bool(
        target_value <= target_best * (1.0 + PRIMARY_TOLERANCE)
        and pseudo_value <= pseudo_best * (1.0 + PRIMARY_TOLERANCE)
    )


def independently_verify_times(headers: list[fits.Header]) -> pd.DataFrame:
    """Recompute barycentric times from stored mid-UTC and from DATE-OBS."""
    iers.conf.auto_download = False
    jd_utc = np.array([float(header["JD_UTC"]) for header in headers])
    stored_bjd = np.array([float(header["BJD_TDB"]) for header in headers])
    stored_midpoint = Time(jd_utc, format="jd", scale="utc", location=OBSERVATORY)
    recalculated_bjd = (
        stored_midpoint.tdb
        + stored_midpoint.light_travel_time(TARGET_COORD, kind="barycentric")
    ).jd

    start_times = Time(
        [str(header["DATE-OBS"]) for header in headers],
        format="isot",
        scale="utc",
        location=OBSERVATORY,
    )
    exposure_seconds = np.array([float(header["EXPTIME"]) for header in headers])
    date_midpoints = start_times + TimeDelta(exposure_seconds / 2.0, format="sec")
    date_bjd = (
        date_midpoints.tdb
        + date_midpoints.light_travel_time(TARGET_COORD, kind="barycentric")
    ).jd

    return pd.DataFrame(
        {
            "slice": np.arange(1, len(headers) + 1),
            "date_obs_utc_start": [str(header["DATE-OBS"]) for header in headers],
            "exposure_seconds": exposure_seconds,
            "stored_jd_utc_midpoint": jd_utc,
            "stored_bjd_tdb": stored_bjd,
            "recalculated_bjd_tdb_from_jd_utc": recalculated_bjd,
            "recalculation_minus_stored_seconds": (
                recalculated_bjd - stored_bjd
            )
            * 86400.0,
            "dateobs_midpoint_bjd_tdb": date_bjd,
            "dateobs_midpoint_minus_stored_seconds": (date_bjd - stored_bjd)
            * 86400.0,
        }
    )


def plot_aperture_check(
    output_path: Path,
    *,
    hours: np.ndarray,
    curves: dict[float, np.ndarray],
    aij_curve: np.ndarray,
    metrics: pd.DataFrame,
    primary_pass: bool,
) -> None:
    colors = plt.cm.viridis(np.linspace(0.08, 0.92, len(SOURCE_RADII)))
    fig, axes = plt.subplots(
        3,
        1,
        figsize=(11.8, 10.2),
        gridspec_kw={"height_ratios": [1.35, 1.65, 2.4]},
    )

    axes[0].plot(
        metrics["source_radius_pixels"],
        metrics["target_robust_scatter_ppt"],
        "o-",
        color="#1d5f91",
        linewidth=1.6,
        label="Target differential curve",
    )
    axes[0].plot(
        metrics["source_radius_pixels"],
        metrics["median_pseudo_target_scatter_ppt"],
        "s-",
        color="#9a5c32",
        linewidth=1.6,
        label="Median comparison pseudo-target",
    )
    axes[0].axvline(PRIMARY_RADIUS, color="#6f4b8b", linestyle="--", linewidth=1.2)
    axes[0].set_ylabel("Robust scatter (ppt)")
    axes[0].set_xlabel("Source-aperture radius (pixels)")
    axes[0].legend(loc="best", fontsize=9)
    axes[0].text(
        0.985,
        0.94,
        "35 px passes both 10% precision plateaus"
        if primary_pass
        else "35 px does not pass both precision plateaus",
        transform=axes[0].transAxes,
        ha="right",
        va="top",
        fontsize=9.4,
        color="#246b3c" if primary_pass else "#9d3927",
        bbox={"boxstyle": "round,pad=0.35", "facecolor": "white", "alpha": 0.9, "edgecolor": "#cccccc"},
    )

    primary_curve = curves[PRIMARY_RADIUS]
    axes[1].plot(
        hours,
        1000.0 * (primary_curve - aij_curve),
        ".",
        color="#6f4b8b",
        markersize=4.0,
    )
    axes[1].axhline(0.0, color="#666666", linewidth=0.8)
    axes[1].set_ylabel("Photutils − AstroImageJ\nat 35 px (ppt)")
    axes[1].set_xlabel("Hours since first exposure")
    axes[1].text(
        0.015,
        0.94,
        f"Curve correlation: {correlation(primary_curve, aij_curve):.5f}\n"
        f"Robust difference: {1000.0 * robust_sigma(primary_curve - aij_curve):.3f} ppt",
        transform=axes[1].transAxes,
        ha="left",
        va="top",
        fontsize=9.2,
        bbox={"boxstyle": "round,pad=0.35", "facecolor": "white", "alpha": 0.9, "edgecolor": "#cccccc"},
    )

    separation = 1.2
    for offset, radius, color in zip(range(len(SOURCE_RADII)), SOURCE_RADII, colors):
        percent = 100.0 * (curves[radius] - 1.0) + separation * offset
        axes[2].plot(
            hours,
            percent,
            ".",
            color=color,
            markersize=3.0,
            alpha=0.76,
        )
        axes[2].text(
            float(np.nanmax(hours)) + 0.035,
            separation * offset,
            f"{radius:.0f}px",
            color=color,
            va="center",
            fontsize=8.8,
        )
    axes[2].set_ylabel("Flux change + display offset (%)")
    axes[2].set_xlabel("Hours since first exposure (BJD TDB)")
    axes[2].set_xlim(float(np.nanmin(hours)), float(np.nanmax(hours)) + 0.23)

    for axis in axes:
        axis.grid(alpha=0.16)
    fig.suptitle(
        "TOI-3505.01 — independent source-aperture robustness check",
        fontsize=15.5,
        y=0.985,
    )
    fig.text(
        0.5,
        0.956,
        "281 aligned R-band frames · 70–139 px sigma-clipped sky annulus · all 10 comparison stars",
        ha="center",
        fontsize=10,
    )
    fig.tight_layout(rect=(0.04, 0.04, 0.97, 0.935))
    fig.savefig(output_path, dpi=220, facecolor="white")
    plt.close(fig)


def write_readme(
    output_dir: Path,
    *,
    primary_pass: bool,
    maximum_bjd_error_seconds: float,
    maximum_dateobs_offset_seconds: float,
    curve_correlation: float,
    curve_difference_ppt: float,
) -> None:
    verdict = (
        "The existing 35-pixel source radius passes the frozen target and "
        "pseudo-target precision-plateau rule."
        if primary_pass
        else "The existing 35-pixel source radius fails at least one frozen "
        "precision-plateau rule and should not be treated as finalized."
    )
    text = f"""# TOI-3505.01 aperture-radius and time check

## Aperture result

{verdict}

Seven source radii (20, 25, 30, 35, 40, 45, and 50 pixels) were measured from
all 281 aligned images. The sky annulus remained 70-139 pixels. Photutils used
the AstroImageJ centroid positions but independently summed the images and
estimated each local sky with a three-sigma-clipped annulus median.

The 35-pixel Photutils and AstroImageJ differential curves have correlation
{curve_correlation:.6f}; their robust point-by-point difference is
{curve_difference_ppt:.3f} ppt. Differences are expected because the two
programs do not use identical aperture-edge or background estimators.

The frozen retention rule does not simply pick the lowest target scatter. The
35-pixel radius must be within 10% of both (1) the best target differential
scatter and (2) the best median comparison pseudo-target scatter. This guards
against tuning the aperture only to make the target look flat.

## Time result

Recomputing BJD_TDB with Astropy from each stored `JD_UTC` midpoint, the GMU
site coordinates, and the J2000 target coordinate agrees with every stored
`BJD_TDB` to within {maximum_bjd_error_seconds:.6f} seconds. Reconstructing
mid-exposure from the lower-precision `DATE-OBS` strings plus 25 seconds differs
by at most {maximum_dateobs_offset_seconds:.3f} seconds. The stored barycentric
conversion is therefore internally verified; observatory clock-sync provenance
still needs confirmation from the observer or mentor.

## Files

- `01_aperture_radius_check.png`: radius metrics, AstroImageJ comparison, and curves.
- `aperture_radius_metrics.csv`: one row per tested source radius.
- `python_multi_radius_light_curves.csv`: every extracted curve.
- `time_verification.csv`: all 281 independent timing calculations.
- `summary.json`: machine-readable verdict and statistics.

This is an independent implementation of the same observation, not an
independent astrophysical data set. It does not remove the need for mentor
review of the AstroImageJ settings.
"""
    (output_dir / "README.md").write_text(text, encoding="utf-8")


def main() -> None:
    args = parse_args()
    aligned_dir = args.aligned_dir.resolve()
    table_path = args.table.resolve()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    table = load_table(table_path, expected_outer_radius=139)
    frame_fluxes = np.full(
        (len(table), len(SOURCE_RADII), len(STAR_NAMES)), np.nan, dtype=float
    )
    backgrounds = np.full((len(table), len(STAR_NAMES)), np.nan, dtype=float)
    headers: list[fits.Header] = []

    for row_index, row in table.iterrows():
        path = aligned_dir / str(row["Label"])
        if not path.exists():
            raise FileNotFoundError(path)
        positions = np.array(
            [
                (
                    float(row[f"X(FITS)_{star}"]) - 1.0,
                    float(row[f"Y(FITS)_{star}"]) - 1.0,
                )
                for star in STAR_NAMES
            ]
        )
        with fits.open(path, memmap=True) as hdul:
            headers.append(hdul[0].header.copy())
            flux, background = measure_frame(hdul[0].data, positions)
        frame_fluxes[row_index] = flux
        backgrounds[row_index] = background
        if (row_index + 1) % 25 == 0 or row_index + 1 == len(table):
            print(f"Measured {row_index + 1}/{len(table)} frames", flush=True)

    aij_curve = normalized(table["rel_flux_T1"])
    hours = (
        table["BJD_TDB"].to_numpy(dtype=float)
        - float(table["BJD_TDB"].iloc[0])
    ) * 24.0
    curves: dict[float, np.ndarray] = {}
    metric_rows: list[dict[str, float]] = []
    curve_columns: dict[str, np.ndarray | pd.Series] = {
        "image": table["Label"],
        "slice": table["slice"],
        "bjd_tdb": table["BJD_TDB"],
        "hours_since_first_image": hours,
        "astroimagej_35px_relative_brightness": aij_curve,
    }
    for radius_index, radius in enumerate(SOURCE_RADII):
        flux = frame_fluxes[:, radius_index, :]
        curve = differential_curve(flux)
        curves[radius] = curve
        curve_columns[f"photutils_{int(radius)}px_relative_brightness"] = curve
        metric_rows.append(
            {
                "source_radius_pixels": radius,
                "target_robust_scatter_ppt": 1000.0 * robust_scatter(curve),
                "target_clipped_standard_deviation_ppt": 1000.0
                * clipped_standard_deviation(curve),
                "target_full_standard_deviation_ppt": 1000.0
                * float(np.nanstd(curve)),
                "median_pseudo_target_scatter_ppt": 1000.0
                * median_pseudo_target_scatter(flux),
                "correlation_with_astroimagej_35px": correlation(curve, aij_curve),
                "robust_difference_from_astroimagej_35px_ppt": 1000.0
                * robust_sigma(curve - aij_curve),
            }
        )
    metrics = pd.DataFrame(metric_rows)
    primary_pass = primary_radius_passes(metrics)
    metrics["passes_target_10_percent_plateau"] = metrics[
        "target_robust_scatter_ppt"
    ] <= float(metrics["target_robust_scatter_ppt"].min()) * (
        1.0 + PRIMARY_TOLERANCE
    )
    metrics["passes_pseudo_target_10_percent_plateau"] = metrics[
        "median_pseudo_target_scatter_ppt"
    ] <= float(metrics["median_pseudo_target_scatter_ppt"].min()) * (
        1.0 + PRIMARY_TOLERANCE
    )

    time_verification = independently_verify_times(headers)
    maximum_bjd_error_seconds = float(
        np.max(np.abs(time_verification["recalculation_minus_stored_seconds"]))
    )
    maximum_dateobs_offset_seconds = float(
        np.max(np.abs(time_verification["dateobs_midpoint_minus_stored_seconds"]))
    )
    timing_pass = maximum_bjd_error_seconds < 0.1

    primary_curve = curves[PRIMARY_RADIUS]
    curve_correlation = correlation(primary_curve, aij_curve)
    curve_difference_ppt = 1000.0 * robust_sigma(primary_curve - aij_curve)
    primary_flux = frame_fluxes[:, SOURCE_RADII.index(PRIMARY_RADIUS), :]
    aij_target_counts = table["Source-Sky_T1"].to_numpy(dtype=float)
    target_count_ratio = primary_flux[:, 0] / aij_target_counts

    pd.DataFrame(curve_columns).to_csv(
        output_dir / "python_multi_radius_light_curves.csv",
        index=False,
        float_format="%.10f",
    )
    metrics.to_csv(
        output_dir / "aperture_radius_metrics.csv",
        index=False,
        float_format="%.10f",
    )
    time_verification.to_csv(
        output_dir / "time_verification.csv", index=False, float_format="%.12f"
    )
    plot_aperture_check(
        output_dir / "01_aperture_radius_check.png",
        hours=hours,
        curves=curves,
        aij_curve=aij_curve,
        metrics=metrics,
        primary_pass=primary_pass,
    )

    summary = {
        "status": "pass" if primary_pass and timing_pass else "review_required",
        "implementation": "Photutils 3 sigma-clipped median-annulus extraction",
        "images_measured": int(len(table)),
        "source_radii_pixels": list(SOURCE_RADII),
        "sky_annulus_pixels": [ANNULUS_INNER, ANNULUS_OUTER],
        "primary_radius_pixels": PRIMARY_RADIUS,
        "primary_radius_passes_frozen_precision_plateaus": primary_pass,
        "primary_radius_tolerance_fraction": PRIMARY_TOLERANCE,
        "astroimagej_photutils_35px_curve_correlation": curve_correlation,
        "astroimagej_photutils_35px_robust_curve_difference_ppt": (
            curve_difference_ppt
        ),
        "photutils_to_astroimagej_target_count_ratio": {
            "median": float(np.nanmedian(target_count_ratio)),
            "robust_scatter": robust_sigma(target_count_ratio),
        },
        "timing": {
            "bjd_tdb_recalculation_pass": timing_pass,
            "maximum_absolute_recalculation_error_seconds": (
                maximum_bjd_error_seconds
            ),
            "maximum_absolute_dateobs_midpoint_offset_seconds": (
                maximum_dateobs_offset_seconds
            ),
            "clock_sync_provenance_confirmed": False,
        },
        "metrics": metrics.to_dict(orient="records"),
    }
    (output_dir / "summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    write_readme(
        output_dir,
        primary_pass=primary_pass,
        maximum_bjd_error_seconds=maximum_bjd_error_seconds,
        maximum_dateobs_offset_seconds=maximum_dateobs_offset_seconds,
        curve_correlation=curve_correlation,
        curve_difference_ppt=curve_difference_ppt,
    )

    print(f"Wrote aperture-check products to {output_dir}")
    print(f"35-pixel primary passes: {primary_pass}")
    print(f"35-pixel AstroImageJ/Photutils curve correlation: {curve_correlation:.6f}")
    print(f"Maximum BJD_TDB recalculation error: {maximum_bjd_error_seconds:.6f} s")


if __name__ == "__main__":
    main()
