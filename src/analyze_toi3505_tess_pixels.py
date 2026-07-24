"""Run custom-aperture and difference-image checks for TOI-3505.01.

This script uses the same 15x15-pixel TESScut size in Sectors 14, 41, 54,
and 81.  It measures how the observed box depth changes with aperture radius,
checks whether the in-minus-out signal is centered within one TESS pixel of
the target, overlays the Gaia field, and runs light-curve injection tests.

The localization result cannot distinguish the target from the 0.517-arcsec
companion because one TESS pixel is about 21 arcsec across.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import tempfile
import warnings
from dataclasses import dataclass
from pathlib import Path

os.environ.setdefault(
    "MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "nasa-tess-matplotlib-cache")
)

import astropy.units as u
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.utils.exceptions import AstropyWarning
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales
from matplotlib.patches import Circle

from toi3505_tess_tools import (
    LightCurveData,
    event_cycles,
    fit_repeating_box,
    grid_box_fit,
    integrated_box_fraction,
    phase_offset,
    robust_scatter,
)


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DATA_DIR = ROOT / "data" / "tess" / "toi3505" / "tesscut"
DEFAULT_CATALOG_DIR = ROOT / "data" / "catalogs" / "toi3505"
DEFAULT_TESS_RESULTS = ROOT / "outputs" / "toi3505_tess_analysis"
DEFAULT_OUTPUT_DIR = ROOT / "outputs" / "toi3505_tess_pixels"

SECTORS = (14, 41, 54, 81)
APERTURE_RADII = (1.5, 2.0, 2.5, 3.0, 3.5, 4.0)
PRIMARY_RADIUS = 3.0
BACKGROUND_INNER_RADIUS = 4.0
BACKGROUND_OUTER_RADIUS = 7.0
PERIOD_DAYS = 2.9151556
EPOCH_BJD = 2459793.534385
DURATION_DAYS = 2.004 / 24.0
CATALOG_DEPTH = 0.002910
TARGET_RA_DEG = 297.043476
TARGET_DEC_DEG = 18.698914
SPECKLE_SEPARATION_ARCSEC = 0.517
SPECKLE_DELTA_I_MAG = 1.7
RANDOM_SEED = 3505


@dataclass
class CutoutData:
    path: Path
    sector: int
    time_bjd: np.ndarray
    flux: np.ndarray
    flux_error: np.ndarray
    quality: np.ndarray
    cadence_days: float
    wcs: WCS
    target_x: float
    target_y: float
    pixel_scale_arcsec: float


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=DEFAULT_DATA_DIR)
    parser.add_argument("--catalog-dir", type=Path, default=DEFAULT_CATALOG_DIR)
    parser.add_argument("--tess-results", type=Path, default=DEFAULT_TESS_RESULTS)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    return parser.parse_args()


def circular_mask(
    shape: tuple[int, int], center_x: float, center_y: float, radius: float
) -> np.ndarray:
    y, x = np.indices(shape, dtype=float)
    return np.hypot(x - center_x, y - center_y) <= radius


def annulus_mask(
    shape: tuple[int, int],
    center_x: float,
    center_y: float,
    inner_radius: float,
    outer_radius: float,
) -> np.ndarray:
    y, x = np.indices(shape, dtype=float)
    distance = np.hypot(x - center_x, y - center_y)
    return (distance > inner_radius) & (distance <= outer_radius)


def local_background(cutout: CutoutData) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Estimate one sky value per cadence from the QLP-style 4-7 pixel annulus."""
    mask = annulus_mask(
        cutout.flux.shape[1:],
        cutout.target_x,
        cutout.target_y,
        BACKGROUND_INNER_RADIUS,
        BACKGROUND_OUTER_RADIUS,
    )
    pixels = cutout.flux[:, mask]
    errors = cutout.flux_error[:, mask]
    background = np.nanmedian(pixels, axis=1)
    # 1.253 converts the standard error of a mean to the approximate standard
    # error of a median for normally distributed measurements.
    background_error = 1.253 * np.sqrt(np.nansum(np.square(errors), axis=1)) / mask.sum()
    return background, background_error, mask


def load_cutout(path: Path) -> CutoutData:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", AstropyWarning)
        with fits.open(path, memmap=True) as hdul:
            primary = hdul[0].header
            table_header = hdul[1].header
            table = hdul[1].data
            time_bjd = np.asarray(table["TIME"], dtype=float) + float(
                table_header.get("BJDREFI", 2457000.0)
            ) + float(table_header.get("BJDREFF", 0.0))
            flux = np.asarray(table["FLUX"], dtype=float)
            flux_error = np.asarray(table["FLUX_ERR"], dtype=float)
            quality = np.asarray(table["QUALITY"], dtype=np.int64)
            wcs = WCS(hdul[2].header)
            sector = int(primary["SECTOR"])
            cadence_days = float(table_header["TIMEDEL"])
    target = SkyCoord(TARGET_RA_DEG * u.deg, TARGET_DEC_DEG * u.deg)
    target_x, target_y = wcs.world_to_pixel(target)
    pixel_scale = float(np.mean(proj_plane_pixel_scales(wcs)) * 3600.0)
    return CutoutData(
        path=path,
        sector=sector,
        time_bjd=time_bjd,
        flux=flux,
        flux_error=flux_error,
        quality=quality,
        cadence_days=cadence_days,
        wcs=wcs,
        target_x=float(target_x),
        target_y=float(target_y),
        pixel_scale_arcsec=pixel_scale,
    )


def cutouts_equivalent(first: CutoutData, second: CutoutData) -> bool:
    """Return true only when duplicate files contain the same science data."""
    return bool(
        first.sector == second.sector
        and math.isclose(first.cadence_days, second.cadence_days)
        and math.isclose(first.target_x, second.target_x, abs_tol=1e-8)
        and math.isclose(first.target_y, second.target_y, abs_tol=1e-8)
        and math.isclose(
            first.pixel_scale_arcsec,
            second.pixel_scale_arcsec,
            abs_tol=1e-8,
        )
        and np.array_equal(first.time_bjd, second.time_bjd, equal_nan=True)
        and np.array_equal(first.flux, second.flux, equal_nan=True)
        and np.array_equal(first.flux_error, second.flux_error, equal_nan=True)
        and np.array_equal(first.quality, second.quality)
    )


def find_cutouts(data_dir: Path) -> dict[int, CutoutData]:
    cutouts: dict[int, CutoutData] = {}
    for path in sorted(data_dir.rglob("*.fits")):
        cutout = load_cutout(path)
        if cutout.sector in cutouts:
            previous = cutouts[cutout.sector]
            if not cutouts_equivalent(previous, cutout):
                raise RuntimeError(
                    f"Conflicting TESScut files for Sector {cutout.sector}: "
                    f"{previous.path} and {cutout.path}"
                )
            print(
                f"Sector {cutout.sector}: equivalent duplicate TESScut file "
                f"ignored ({cutout.path})"
            )
            continue
        cutouts[cutout.sector] = cutout
    if set(cutouts) != set(SECTORS):
        raise RuntimeError(f"TESScut sectors found: {sorted(cutouts)}; expected {SECTORS}")
    return cutouts


def extract_aperture(cutout: CutoutData, radius: float) -> tuple[LightCurveData, np.ndarray]:
    mask = circular_mask(
        cutout.flux.shape[1:], cutout.target_x, cutout.target_y, radius
    )
    background, background_error, _ = local_background(cutout)
    flux = np.nansum(cutout.flux[:, mask], axis=1) - background * mask.sum()
    error = np.sqrt(
        np.nansum(np.square(cutout.flux_error[:, mask]), axis=1)
        + np.square(background_error * mask.sum())
    )
    finite_cube = np.all(np.isfinite(cutout.flux[:, mask]), axis=1)
    quality = cutout.quality.copy()
    quality[~finite_cube] |= 1
    good = (
        (quality == 0)
        & np.isfinite(cutout.time_bjd)
        & np.isfinite(flux)
        & np.isfinite(error)
        & (error > 0)
    )
    scale = float(np.nanmedian(flux[good]))
    if not np.isfinite(scale) or scale <= 0:
        raise RuntimeError(f"Invalid Sector {cutout.sector} aperture flux")
    curve = LightCurveData(
        path=cutout.path,
        sector=cutout.sector,
        pipeline="TESScut",
        flux_name=f"circular_radius_{radius:.1f}_pixels",
        time_bjd=cutout.time_bjd,
        flux=flux / scale,
        flux_error=error / scale,
        quality=quality,
        cadence_days=cutout.cadence_days,
        crowdsap=None,
        flfrcsap=None,
    )
    return curve, mask


def measure_apertures(
    cutouts: dict[int, CutoutData], output_dir: Path
) -> tuple[pd.DataFrame, dict[int, dict[float, LightCurveData]]]:
    rows: list[dict[str, object]] = []
    curves: dict[int, dict[float, LightCurveData]] = {}
    offsets = np.linspace(-0.035, 0.035, 71)
    for sector in SECTORS:
        curves[sector] = {}
        for radius in APERTURE_RADII:
            curve, mask = extract_aperture(cutouts[sector], radius)
            curves[sector][radius] = curve
            fit, _, bounds = grid_box_fit(
                curve,
                PERIOD_DAYS,
                EPOCH_BJD,
                durations_days=np.array([DURATION_DAYS]),
                offsets_days=offsets,
            )
            rows.append(
                {
                    "sector": sector,
                    "radius_pixels": radius,
                    "aperture_pixels": int(mask.sum()),
                    "cadence_minutes": curve.cadence_days * 1440.0,
                    "quality_zero_points": int(curve.good.sum()),
                    "observed_depth_ppt": fit.depth * 1000.0,
                    "depth_error_ppt": fit.depth_error * 1000.0,
                    "midpoint_offset_minutes": fit.time_offset_days * 1440.0,
                    "midpoint_offset_low_minutes": bounds["offset_low_days"] * 1440.0,
                    "midpoint_offset_high_minutes": bounds["offset_high_days"] * 1440.0,
                    "local_residual_scatter_ppt": fit.residual_scatter * 1000.0,
                    "depth_snr": fit.depth / fit.depth_error,
                    "note": "custom extraction from the same TESS pixels; not an independent observation",
                }
            )
    table = pd.DataFrame(rows)
    table.to_csv(output_dir / "depth_vs_aperture.csv", index=False)
    return table, curves


def event_difference_images(
    cutout: CutoutData,
    midpoint_offset_days: float,
    aperture_mask: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, list[int]]:
    good = (
        (cutout.quality == 0)
        & np.isfinite(cutout.time_bjd)
        & np.all(np.isfinite(cutout.flux), axis=(1, 2))
    )
    background, _, _ = local_background(cutout)
    corrected_cube = cutout.flux - background[:, None, None]
    cycles = event_cycles(cutout.time_bjd, PERIOD_DAYS, EPOCH_BJD)
    differences: list[np.ndarray] = []
    used_cycles: list[int] = []
    for cycle in np.unique(cycles[good]):
        midpoint = EPOCH_BJD + cycle * PERIOD_DAYS + midpoint_offset_days
        dt = cutout.time_bjd - midpoint
        in_event = good & (np.abs(dt) <= DURATION_DAYS / 2.0)
        out_event = good & (np.abs(dt) >= 0.9 * DURATION_DAYS) & (
            np.abs(dt) <= 2.4 * DURATION_DAYS
        )
        if in_event.sum() < 2 or out_event.sum() < 4:
            continue
        in_image = np.nanmedian(corrected_cube[in_event], axis=0)
        out_image = np.nanmedian(corrected_cube[out_event], axis=0)
        aperture_flux = float(np.nansum(out_image[aperture_mask]))
        if not np.isfinite(aperture_flux) or aperture_flux <= 0:
            continue
        differences.append((out_image - in_image) / aperture_flux)
        used_cycles.append(int(cycle))
    if len(differences) < 3:
        raise RuntimeError(
            f"Sector {cutout.sector} has only {len(differences)} usable difference images"
        )
    stack = np.asarray(differences)
    return np.nanmedian(stack, axis=0), stack, used_cycles


def positive_centroid(
    image: np.ndarray, target_x: float, target_y: float, radius: float = 3.0
) -> tuple[float, float]:
    mask = circular_mask(image.shape, target_x, target_y, radius)
    weights = np.where(mask, np.clip(image, 0.0, None), 0.0)
    total = float(np.sum(weights))
    if total <= 0:
        return float("nan"), float("nan")
    y, x = np.indices(image.shape, dtype=float)
    return float(np.sum(x * weights) / total), float(np.sum(y * weights) / total)


def difference_localization(
    cutouts: dict[int, CutoutData],
    tess_results: Path,
    output_dir: Path,
) -> tuple[pd.DataFrame, dict[int, np.ndarray]]:
    sector_table = pd.read_csv(tess_results / "sector_measurements.csv")
    rng = np.random.default_rng(RANDOM_SEED)
    rows = []
    images: dict[int, np.ndarray] = {}
    bootstrap_count = 600
    for sector in SECTORS:
        cutout = cutouts[sector]
        offset_minutes = float(
            sector_table.loc[
                sector_table["sector"] == sector, "midpoint_offset_minutes"
            ].iloc[0]
        )
        aperture = circular_mask(
            cutout.flux.shape[1:], cutout.target_x, cutout.target_y, PRIMARY_RADIUS
        )
        image, stack, cycles = event_difference_images(
            cutout, offset_minutes / 1440.0, aperture
        )
        images[sector] = image
        centroid_x, centroid_y = positive_centroid(
            image, cutout.target_x, cutout.target_y
        )
        bootstrap_centroids = []
        for _ in range(bootstrap_count):
            chosen = rng.integers(0, len(stack), len(stack))
            sample = np.nanmedian(stack[chosen], axis=0)
            bootstrap_centroids.append(
                positive_centroid(sample, cutout.target_x, cutout.target_y)
            )
        bootstrap_centroids = np.asarray(bootstrap_centroids, dtype=float)
        finite = np.all(np.isfinite(bootstrap_centroids), axis=1)
        offsets = np.hypot(
            bootstrap_centroids[finite, 0] - centroid_x,
            bootstrap_centroids[finite, 1] - centroid_y,
        )
        centroid_error_arcsec = (
            float(np.nanpercentile(offsets, 68.27) * cutout.pixel_scale_arcsec)
            if len(offsets)
            else float("nan")
        )
        offset_pixels = float(
            np.hypot(centroid_x - cutout.target_x, centroid_y - cutout.target_y)
        )
        event_depths = np.nansum(stack[:, aperture], axis=1)
        depth_error = float(np.std(event_depths, ddof=1) / np.sqrt(len(event_depths)))
        rows.append(
            {
                "sector": sector,
                "events": len(stack),
                "cycles": ";".join(str(value) for value in cycles),
                "target_x": cutout.target_x,
                "target_y": cutout.target_y,
                "difference_centroid_x": centroid_x,
                "difference_centroid_y": centroid_y,
                "centroid_offset_pixels": offset_pixels,
                "centroid_offset_arcsec": offset_pixels * cutout.pixel_scale_arcsec,
                "bootstrap_68pct_centroid_error_arcsec": centroid_error_arcsec,
                "pixel_scale_arcsec": cutout.pixel_scale_arcsec,
                "difference_aperture_depth_ppt": float(np.median(event_depths) * 1000.0),
                "difference_aperture_depth_error_ppt": depth_error * 1000.0,
                "consistent_with_target_within_one_pixel": bool(offset_pixels <= 1.0),
                "angular_limit": "TESS cannot separate the 0.517-arcsec companion from the target",
            }
        )
    table = pd.DataFrame(rows)
    table.to_csv(output_dir / "difference_image_localization.csv", index=False)
    return table, images


def build_neighbor_screen(
    catalog_dir: Path, output_dir: Path
) -> tuple[pd.DataFrame, dict[str, object]]:
    tic = pd.read_csv(catalog_dir / "tic_v8_2p5arcmin.csv")
    target = tic.loc[tic["ID"] == 390988385].iloc[0]
    target_tmag = float(target["Tmag"])
    nearby = tic[
        np.isfinite(tic["Tmag"])
        & np.isfinite(tic["dstArcSec"])
        & (tic["dstArcSec"] <= 60.0)
    ].copy()
    nearby["delta_tmag"] = nearby["Tmag"] - target_tmag
    nearby["flux_ratio_to_target"] = np.power(10.0, -0.4 * nearby["delta_tmag"])
    nearby["eclipse_fraction_needed_simple"] = (
        CATALOG_DEPTH / nearby["flux_ratio_to_target"]
    )
    nearby["screening_verdict"] = np.where(
        nearby["eclipse_fraction_needed_simple"] <= 1.0,
        "bright enough to mimic the catalog depth in simple full-throughput arithmetic",
        "too faint to mimic the catalog depth even in a total eclipse under the simple screen",
    )
    columns = [
        "ID",
        "GAIA",
        "ra",
        "dec",
        "dstArcSec",
        "Tmag",
        "delta_tmag",
        "flux_ratio_to_target",
        "eclipse_fraction_needed_simple",
        "screening_verdict",
    ]
    nearby[columns].sort_values("dstArcSec").to_csv(
        output_dir / "tic_neighbor_screen_60arcsec.csv", index=False
    )

    companion_ratio = float(10.0 ** (-0.4 * SPECKLE_DELTA_I_MAG))
    tic_contamination = float(target["contratio"])
    total_contamination = tic_contamination + companion_ratio
    budget = {
        "target_tmag": target_tmag,
        "tic_catalog_contamination_ratio": tic_contamination,
        "unresolved_companion": {
            "separation_arcsec": SPECKLE_SEPARATION_ARCSEC,
            "delta_i_mag": SPECKLE_DELTA_I_MAG,
            "flux_ratio_using_delta_i_as_tess_band_proxy": companion_ratio,
            "warning": "Delta I is only a screening proxy for the TESS-band flux ratio.",
        },
        "screening_sum_contaminating_flux_over_target_flux": total_contamination,
        "screening_target_fraction_of_total_flux": 1.0 / (1.0 + total_contamination),
        "if_2p91_ppt_is_an_uncorrected_observed_depth": {
            "target_host_depth_ppt": CATALOG_DEPTH * (1.0 + total_contamination) * 1000.0,
            "target_host_radius_ratio": float(
                np.sqrt(CATALOG_DEPTH * (1.0 + total_contamination))
            ),
            "companion_eclipse_fraction_including_screening_contamination": CATALOG_DEPTH
            * (1.0 + total_contamination)
            / companion_ratio,
        },
        "do_not_adopt_as_final_correction": "QLP dilution treatment and band-dependent companion flux must be resolved first; this is scenario arithmetic.",
    }
    (output_dir / "dilution_screen.json").write_text(
        json.dumps(budget, indent=2) + "\n", encoding="utf-8"
    )
    return nearby[columns], budget


def run_injections(
    curves: dict[int, dict[float, LightCurveData]], output_dir: Path
) -> tuple[pd.DataFrame, pd.DataFrame]:
    depths = (0.002, 0.003, 0.004)
    phase_fractions = (0.17, 0.29, 0.41, 0.59, 0.71, 0.83)
    rows = []
    for sector in SECTORS:
        curve = curves[sector][PRIMARY_RADIUS]
        for fraction in phase_fractions:
            shifted_epoch = EPOCH_BJD + fraction * PERIOD_DAYS
            control = fit_repeating_box(
                curve,
                PERIOD_DAYS,
                shifted_epoch,
                DURATION_DAYS,
                window_days=0.22,
            )
            phase = phase_offset(curve.time_bjd, PERIOD_DAYS, shifted_epoch)
            for depth in depths:
                exposure = integrated_box_fraction(
                    phase, DURATION_DAYS, curve.cadence_days
                )
                injected_curve = LightCurveData(
                    path=curve.path,
                    sector=curve.sector,
                    pipeline=curve.pipeline,
                    flux_name=curve.flux_name,
                    time_bjd=curve.time_bjd,
                    flux=curve.flux * (1.0 - depth * exposure),
                    flux_error=curve.flux_error,
                    quality=curve.quality,
                    cadence_days=curve.cadence_days,
                    crowdsap=None,
                    flfrcsap=None,
                )
                recovered = fit_repeating_box(
                    injected_curve,
                    PERIOD_DAYS,
                    shifted_epoch,
                    DURATION_DAYS,
                    window_days=0.22,
                )
                increment = recovered.depth - control.depth
                rows.append(
                    {
                        "sector": sector,
                        "aperture_radius_pixels": PRIMARY_RADIUS,
                        "phase_fraction": fraction,
                        "injected_depth_ppt": depth * 1000.0,
                        "control_depth_ppt": control.depth * 1000.0,
                        "recovered_total_depth_ppt": recovered.depth * 1000.0,
                        "recovered_injected_increment_ppt": increment * 1000.0,
                        "total_recovery_bias_ppt": (recovered.depth - depth)
                        * 1000.0,
                        "total_recovery_fraction": recovered.depth / depth,
                        "increment_recovery_bias_ppt": (increment - depth)
                        * 1000.0,
                        "increment_recovery_fraction": increment / depth,
                        "formal_depth_error_ppt": recovered.depth_error * 1000.0,
                        "scope": "injection into extracted light curve; tests fitting bias, not target-only PRF dilution",
                    }
                )
    table = pd.DataFrame(rows)
    table.to_csv(output_dir / "light_curve_injections.csv", index=False)
    summary = (
        table.groupby(["sector", "injected_depth_ppt"], as_index=False)
        .agg(
            mean_total_recovered_depth_ppt=("recovered_total_depth_ppt", "mean"),
            total_recovered_depth_std_ppt=("recovered_total_depth_ppt", "std"),
            mean_total_recovery_fraction=("total_recovery_fraction", "mean"),
            total_recovery_fraction_std=("total_recovery_fraction", "std"),
            mean_total_bias_ppt=("total_recovery_bias_ppt", "mean"),
            total_bias_std_ppt=("total_recovery_bias_ppt", "std"),
            mean_increment_recovery_fraction=(
                "increment_recovery_fraction",
                "mean",
            ),
            increment_recovery_fraction_std=(
                "increment_recovery_fraction",
                "std",
            ),
            mean_increment_bias_ppt=("increment_recovery_bias_ppt", "mean"),
            trials=("total_recovery_fraction", "size"),
        )
    )
    summary.to_csv(output_dir / "light_curve_injection_summary.csv", index=False)
    return table, summary


def make_aperture_plot(table: pd.DataFrame, output_dir: Path) -> None:
    figure, ax = plt.subplots(figsize=(9.2, 5.7))
    markers = {14: "o", 41: "s", 54: "^", 81: "D"}
    for sector in SECTORS:
        rows = table[table["sector"] == sector]
        ax.errorbar(
            rows["radius_pixels"],
            rows["observed_depth_ppt"],
            yerr=rows["depth_error_ppt"],
            marker=markers[sector],
            markersize=6,
            capsize=3,
            linewidth=1.2,
            label=f"Sector {sector}",
        )
    ax.axvline(PRIMARY_RADIUS, color="0.35", linestyle="--", linewidth=1.0)
    ax.set_xlabel("Circular aperture radius (TESS pixels)")
    ax.set_ylabel("Observed box depth (ppt)")
    ax.set_title("TOI-3505.01 depth changes with aperture size")
    ax.legend(ncol=2)
    ax.grid(alpha=0.18)
    figure.tight_layout()
    figure.savefig(output_dir / "01_depth_vs_aperture.png", dpi=240)
    figure.savefig(output_dir / "01_depth_vs_aperture.svg")
    plt.close(figure)


def make_field_map(
    cutouts: dict[int, CutoutData], catalog_dir: Path, output_dir: Path
) -> None:
    gaia = pd.read_csv(catalog_dir / "gaia_dr3_2p5arcmin.csv")
    bright = gaia[np.isfinite(gaia["phot_g_mean_mag"]) & (gaia["phot_g_mean_mag"] < 18.0)]
    coordinates = SkyCoord(
        bright["ra"].to_numpy() * u.deg, bright["dec"].to_numpy() * u.deg
    )
    figure, axes = plt.subplots(2, 2, figsize=(10, 9.2))
    for ax, sector in zip(axes.flat, SECTORS):
        cutout = cutouts[sector]
        good = (cutout.quality == 0) & np.all(np.isfinite(cutout.flux), axis=(1, 2))
        image = np.nanmedian(cutout.flux[good], axis=0)
        low, high = np.nanpercentile(image, [10, 99.5])
        ax.imshow(image, origin="lower", cmap="gray", vmin=low, vmax=high)
        x, y = cutout.wcs.world_to_pixel(coordinates)
        inside = (x >= -0.5) & (x <= 14.5) & (y >= -0.5) & (y <= 14.5)
        sizes = np.clip(80.0 - 7.0 * (bright["phot_g_mean_mag"].to_numpy() - 11.0), 12, 70)
        ax.scatter(
            x[inside],
            y[inside],
            s=sizes[inside],
            facecolors="none",
            edgecolors="#d18435",
            linewidths=0.8,
            label="Gaia G < 18",
        )
        ax.scatter(
            cutout.target_x,
            cutout.target_y,
            marker="+",
            s=130,
            linewidths=2.0,
            color="#2d6f9e",
            label="TOI-3505.01",
        )
        ax.add_patch(
            Circle(
                (cutout.target_x, cutout.target_y),
                PRIMARY_RADIUS,
                fill=False,
                color="#2d6f9e",
                linewidth=1.3,
            )
        )
        ax.set_title(f"Sector {sector}")
        ax.set_xlabel("Cutout column (pixels)")
        ax.set_ylabel("Cutout row (pixels)")
        ax.set_xlim(-0.5, 14.5)
        ax.set_ylim(-0.5, 14.5)
    axes[0, 0].legend(loc="upper left", fontsize=8)
    figure.suptitle("TOI-3505.01 field in the four TESS sectors", fontsize=15)
    figure.tight_layout(rect=(0, 0, 1, 0.96))
    figure.savefig(output_dir / "02_tess_field_maps.png", dpi=230)
    figure.savefig(output_dir / "02_tess_field_maps.svg")
    plt.close(figure)


def make_difference_plot(
    cutouts: dict[int, CutoutData],
    localizations: pd.DataFrame,
    images: dict[int, np.ndarray],
    output_dir: Path,
) -> None:
    figure, axes = plt.subplots(2, 2, figsize=(10, 9.0))
    for ax, sector in zip(axes.flat, SECTORS):
        image = images[sector]
        limit = float(np.nanpercentile(np.abs(image), 98))
        if not np.isfinite(limit) or limit == 0:
            limit = 1.0
        ax.imshow(
            image * 1000.0,
            origin="lower",
            cmap="RdBu_r",
            vmin=-limit * 1000.0,
            vmax=limit * 1000.0,
        )
        cutout = cutouts[sector]
        row = localizations[localizations["sector"] == sector].iloc[0]
        ax.scatter(
            cutout.target_x,
            cutout.target_y,
            marker="+",
            s=130,
            linewidths=2,
            color="black",
            label="target position",
        )
        ax.scatter(
            row["difference_centroid_x"],
            row["difference_centroid_y"],
            marker="x",
            s=75,
            linewidths=1.7,
            color="#e3b13b",
            label="difference centroid",
        )
        ax.set_title(
            f"Sector {sector}: offset {row['centroid_offset_arcsec']:.1f} arcsec"
        )
        ax.set_xlabel("Cutout column (pixels)")
        ax.set_ylabel("Cutout row (pixels)")
    axes[0, 0].legend(loc="upper left", fontsize=8)
    figure.suptitle("Out-of-transit minus in-transit images", fontsize=15)
    figure.tight_layout(rect=(0, 0, 1, 0.96))
    figure.savefig(output_dir / "03_difference_images.png", dpi=230)
    figure.savefig(output_dir / "03_difference_images.svg")
    plt.close(figure)


def make_injection_plot(summary: pd.DataFrame, output_dir: Path) -> None:
    figure, ax = plt.subplots(figsize=(8.8, 5.6))
    markers = {14: "o", 41: "s", 54: "^", 81: "D"}
    for sector in SECTORS:
        rows = summary[summary["sector"] == sector]
        ax.errorbar(
            rows["injected_depth_ppt"],
            rows["mean_total_recovered_depth_ppt"],
            yerr=rows["total_recovered_depth_std_ppt"],
            marker=markers[sector],
            markersize=6,
            capsize=3,
            linewidth=1.1,
            label=f"Sector {sector}",
        )
    ax.plot([1.8, 4.2], [1.8, 4.2], color="0.35", linestyle="--", linewidth=1.0)
    ax.set_xlim(1.8, 4.2)
    ax.set_ylim(1.8, 4.2)
    ax.set_xlabel("Injected observed depth (ppt)")
    ax.set_ylabel("Recovered total depth (ppt)")
    ax.set_title("Light-curve injection recovery")
    ax.legend(ncol=2)
    ax.grid(alpha=0.18)
    figure.tight_layout()
    figure.savefig(output_dir / "04_injection_recovery.png", dpi=230)
    figure.savefig(output_dir / "04_injection_recovery.svg")
    plt.close(figure)


def write_readme(
    output_dir: Path,
    aperture_table: pd.DataFrame,
    localization: pd.DataFrame,
    injection_summary: pd.DataFrame,
    budget: dict[str, object],
) -> None:
    primary = aperture_table[aperture_table["radius_pixels"] == PRIMARY_RADIUS]
    lines = [
        "# TOI-3505.01 TESS pixel checks",
        "",
        "This folder contains the custom-aperture, field, difference-image, dilution-screening, and light-curve injection checks for the four TESS sectors.",
        "",
        "## Fixed choices",
        "",
        "- The same 15x15-pixel TESScut size is used in every sector.",
        f"- Circular radii {', '.join(str(value) for value in APERTURE_RADII)} pixels are compared.",
        f"- The 3.0-pixel aperture is the named reference aperture because it matches the QLP best-aperture radius in the downloaded headers.",
        "- The depth comparison fixes the duration to 2.004 hours.",
        "",
        "## Reference-aperture results",
        "",
        "| Sector | Observed depth (ppt) | Difference-centroid offset (arcsec) | One-pixel check |",
        "|---:|---:|---:|:---:|",
    ]
    for _, depth_row in primary.iterrows():
        loc = localization[localization["sector"] == depth_row["sector"]].iloc[0]
        lines.append(
            f"| {int(depth_row['sector'])} | {depth_row['observed_depth_ppt']:.3f} +/- {depth_row['depth_error_ppt']:.3f} | "
            f"{loc['centroid_offset_arcsec']:.1f} | {'yes' if loc['consistent_with_target_within_one_pixel'] else 'no'} |"
        )
    lines.extend(
        [
            "",
            "A one-pixel match is only a TESS-scale localization. It does not distinguish the target primary from the 0.517-arcsec companion.",
            "",
            "## Dilution screen",
            "",
            f"The TIC contamination ratio is {budget['tic_catalog_contamination_ratio']:.3f}. Using delta-I = 1.7 mag as a rough TESS-band proxy adds a companion-to-target flux ratio of {budget['unresolved_companion']['flux_ratio_using_delta_i_as_tess_band_proxy']:.3f}. These values are kept as scenario arithmetic, not applied as a final correction, because QLP's crowding treatment and the companion's TESS-band contrast still need to be resolved.",
            "",
            "## Injection scope",
            "",
            "The injection tests add box-shaped dips to the extracted 3-pixel light curves at six null phases. The main recovery result keeps the control signal at each phase, so the scatter shows real phase-dependent structure. A separate increment column subtracts that control only to check the fitting arithmetic. This tests light-curve recovery, but it does not simulate target-only pixel response or repair an uncertain dilution model.",
            "",
            f"Injection summary rows: {len(injection_summary)}.",
            "",
            "## Main files",
            "",
            "- `depth_vs_aperture.csv`",
            "- `difference_image_localization.csv`",
            "- `tic_neighbor_screen_60arcsec.csv`",
            "- `dilution_screen.json`",
            "- `light_curve_injections.csv` and `light_curve_injection_summary.csv`",
        ]
    )
    (output_dir / "README.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    args = parse_args()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    cutouts = find_cutouts(args.data_dir.resolve())
    aperture_table, curves = measure_apertures(cutouts, output_dir)
    localization, difference_images = difference_localization(
        cutouts, args.tess_results.resolve(), output_dir
    )
    _, dilution_budget = build_neighbor_screen(args.catalog_dir.resolve(), output_dir)
    _, injection_summary = run_injections(curves, output_dir)

    make_aperture_plot(aperture_table, output_dir)
    make_field_map(cutouts, args.catalog_dir.resolve(), output_dir)
    make_difference_plot(cutouts, localization, difference_images, output_dir)
    make_injection_plot(injection_summary, output_dir)
    write_readme(
        output_dir,
        aperture_table,
        localization,
        injection_summary,
        dilution_budget,
    )
    summary = {
        "target": "TOI-3505.01",
        "sectors": list(SECTORS),
        "cutout_size_pixels": [15, 15],
        "aperture_radii_pixels": list(APERTURE_RADII),
        "reference_radius_pixels": PRIMARY_RADIUS,
        "difference_image_random_seed": RANDOM_SEED,
        "localization": localization.to_dict(orient="records"),
        "dilution_screen": dilution_budget,
    }
    (output_dir / "analysis_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(f"Saved TESS pixel checks to {output_dir}")


if __name__ == "__main__":
    main()
