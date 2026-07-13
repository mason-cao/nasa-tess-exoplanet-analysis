"""Practice reduction and differential photometry for the GMU 2023-10-12 folder.

The course spreadsheet associates this date with TOI-6241.01, but every science
frame in the public folder identifies the observed field as TOI-3718.01. This
script therefore treats the observation strictly as an AstroImageJ-workflow
practice data set, never as a TOI-6241.01 research result.

The reduction follows the GMU tutorial sequence: median-combine matching darks,
dark-correct and normalize the R-band flats, dark-subtract and flat-divide the
science frames, track the field, perform circular-aperture differential
photometry, and export reproducible measurements and figures.
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.coordinates import EarthLocation, SkyCoord
from astropy.io import fits
from astropy.stats import SigmaClip
from astropy.time import Time
from astropy.utils import iers
from astropy.visualization import AsinhStretch, ImageNormalize, ZScaleInterval
import astropy.units as u
from photutils.aperture import ApertureStats, CircularAnnulus, CircularAperture, aperture_photometry
from scipy.ndimage import gaussian_filter
from skimage.registration import phase_cross_correlation


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT_DIR = ROOT / "outputs" / "astroimagej_practice"

TARGET_NAME = "TOI-3718.01"
TARGET_COORD = SkyCoord("04h22m51.28s", "+41d07m31.13s", frame="icrs")
OBSERVATORY = EarthLocation.from_geodetic(lon=-77.305345 * u.deg, lat=38.828194 * u.deg, height=148.72 * u.m)
APERTURE_RADIUS = 16.0
ANNULUS_INNER = 22.0
ANNULUS_OUTER = 32.0
EXPECTED_DEPTH = 3770e-6
EXPECTED_DURATION_HOURS = 4.263
FIRST_POINTING = SkyCoord("04h22m24s", "+41d04m25s", frame="icrs")
PLATE_FACTOR = complex(-0.9939341490381487, 0.003129725076419204)
PIXEL_SCALE_ARCSEC = 0.36
REGISTRATION_BIN = 8

# Gaia-matched positions in the first science frame. T1 is the TOI. The
# comparison stars have similar Gaia G magnitudes and avoid the frame edges.
STAR_POSITIONS = {
    "T1": (1235.59, 1895.80),
    "C2": (2274.34, 2213.09),
    "C3": (1848.34, 1673.25),
    "C4": (1438.55, 3301.73),
    "C5": (2073.10, 937.74),
    "C6": (2665.13, 2461.72),
    "C7": (1560.12, 878.11),
    "C8": (801.46, 2537.37),
    "C9": (1268.02, 1405.86),
    "C10": (3238.76, 2576.86),
    "C11": (570.55, 2559.44),
    "C12": (757.50, 2860.94),
    "C13": (1710.04, 685.38),
    "C14": (346.43, 1504.78),
    "C15": (3405.15, 3710.92),
}


def numbered_path_key(path: Path) -> int:
    match = re.search(r"-(\d{4})", path.name)
    return int(match.group(1)) if match else 0


def robust_mad(values: np.ndarray) -> float:
    values = np.asarray(values, dtype=float)
    median = np.nanmedian(values)
    return float(1.4826 * np.nanmedian(np.abs(values - median)))


def median_combine(paths: list[Path]) -> np.ndarray:
    stack = np.empty((len(paths), 4096, 4096), dtype=np.float32)
    for index, path in enumerate(paths):
        stack[index] = fits.getdata(path).astype(np.float32)
    return np.median(stack, axis=0).astype(np.float32)


def build_calibrations(dark_paths: list[Path], flat_paths: list[Path]) -> tuple[dict[float, np.ndarray], np.ndarray]:
    dark_groups: dict[float, list[Path]] = {}
    for path in dark_paths:
        exposure = float(fits.getheader(path)["EXPTIME"])
        dark_groups.setdefault(exposure, []).append(path)

    master_darks = {exposure: median_combine(paths) for exposure, paths in dark_groups.items()}
    if not master_darks:
        raise RuntimeError("No dark frames were found")

    calibrated_flats: list[np.ndarray] = []
    for path in flat_paths:
        header = fits.getheader(path)
        exposure = float(header["EXPTIME"])
        closest = min(master_darks, key=lambda value: abs(value - exposure))
        flat = fits.getdata(path).astype(np.float32) - master_darks[closest]
        central = flat[512:-512, 512:-512]
        level = float(np.nanmedian(central))
        if not np.isfinite(level) or level <= 0:
            continue
        calibrated_flats.append(flat / level)
    if len(calibrated_flats) < 3:
        raise RuntimeError("Fewer than three usable flat frames were found")
    master_flat = np.median(np.stack(calibrated_flats), axis=0).astype(np.float32)
    master_flat /= np.nanmedian(master_flat[512:-512, 512:-512])
    master_flat[~np.isfinite(master_flat) | (master_flat < 0.25)] = 1.0
    return master_darks, master_flat


def local_centroid(image: np.ndarray, x: float, y: float, half_width: int = 24) -> tuple[float, float]:
    if not np.isfinite(x) or not np.isfinite(y):
        return x, y
    x_int = int(round(x))
    y_int = int(round(y))
    x0 = max(0, x_int - half_width)
    x1 = min(image.shape[1], x_int + half_width + 1)
    y0 = max(0, y_int - half_width)
    y1 = min(image.shape[0], y_int + half_width + 1)
    if x1 - x0 < 3 or y1 - y0 < 3:
        return np.nan, np.nan
    cutout = image[y0:y1, x0:x1]
    border = np.concatenate([cutout[0], cutout[-1], cutout[:, 0], cutout[:, -1]])
    background = float(np.nanmedian(border))
    noise = robust_mad(border)
    if not np.isfinite(noise) or noise <= 0 or np.nanmax(cutout - background) < 5.0 * noise:
        return np.nan, np.nan
    weights = np.clip(cutout - background - 2.0 * noise, 0, None)
    total = float(np.nansum(weights))
    if not np.isfinite(total) or total <= 0:
        return np.nan, np.nan
    yy, xx = np.mgrid[y0:y1, x0:x1]
    return float(np.nansum(xx * weights) / total), float(np.nansum(yy * weights) / total)


def compute_bjd_tdb(jd_start: float, exposure_seconds: float) -> float:
    midpoint = Time(jd_start + exposure_seconds / 2.0 / 86400.0, format="jd", scale="utc", location=OBSERVATORY)
    return float((midpoint.tdb + midpoint.light_travel_time(TARGET_COORD)).jd)


def pointing_shift_from_header(header: fits.Header) -> np.ndarray:
    pointing = SkyCoord(
        str(header["OBJCTRA"]),
        str(header["OBJCTDEC"]),
        unit=(u.hourangle, u.deg),
        frame="icrs",
    )
    delta_ra = (pointing.ra - FIRST_POINTING.ra).wrap_at(180 * u.deg).to_value(u.deg)
    delta_dec = (pointing.dec - FIRST_POINTING.dec).to_value(u.deg)
    delta_u = -delta_ra * np.cos(FIRST_POINTING.dec.to_value(u.rad)) * 3600.0 / PIXEL_SCALE_ARCSEC
    delta_v = -delta_dec * 3600.0 / PIXEL_SCALE_ARCSEC
    shift = PLATE_FACTOR * complex(delta_u, delta_v)
    return np.array([shift.real, shift.imag], dtype=float)


def registration_image(image: np.ndarray) -> np.ndarray:
    """Build a small high-pass image for robust all-night field registration.

    The telescope pointing keywords are requested coordinates rather than a
    reliable measured WCS.  Registering every frame directly to frame 1 lets
    the reduction recover independently after the long cloudy interval.
    """
    size = image.shape[0] // REGISTRATION_BIN
    binned = image[: size * REGISTRATION_BIN, : size * REGISTRATION_BIN]
    binned = binned.reshape(size, REGISTRATION_BIN, size, REGISTRATION_BIN).mean((1, 3))
    binned = binned.astype(np.float32)
    binned -= gaussian_filter(binned, 12)
    low, high = np.nanpercentile(binned, [1.0, 99.7])
    binned = np.clip(binned, low, high)
    spread = float(np.nanstd(binned))
    if not np.isfinite(spread) or spread <= 0:
        return np.zeros_like(binned)
    return (binned - np.nanmean(binned)) / spread


def image_registration_shift(reference: np.ndarray, image: np.ndarray) -> np.ndarray:
    moving = registration_image(image)
    align_yx, _, _ = phase_cross_correlation(
        reference,
        moving,
        upsample_factor=2,
        overlap_ratio=0.30,
    )
    # phase_cross_correlation returns the shift needed to move the current
    # image back onto the reference.  Aperture positions need the inverse.
    return np.array(
        [-align_yx[1] * REGISTRATION_BIN, -align_yx[0] * REGISTRATION_BIN],
        dtype=float,
    )


def expected_window_bjd() -> tuple[float, float]:
    # The course spreadsheet gives the 2023-12-06 TOI-3718.01 window as
    # 22:35-02:49 local. Propagating that window backward by 13 periods places
    # it on this 2023-10-12 local observing night.
    period = 4.230596 * u.day
    ingress_utc = Time("2023-12-07T03:35:00", scale="utc", location=OBSERVATORY) - 13 * period
    egress_utc = Time("2023-12-07T07:49:00", scale="utc", location=OBSERVATORY) - 13 * period
    ingress_bjd = ingress_utc.tdb + ingress_utc.light_travel_time(TARGET_COORD)
    egress_bjd = egress_utc.tdb + egress_utc.light_travel_time(TARGET_COORD)
    return float(ingress_bjd.jd), float(egress_bjd.jd)


def select_comparison_stars(fluxes: np.ndarray, initial_good: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    # Columns 1 onward are comparison stars. Assess each against the median of
    # the others, then retain the most stable ensemble.
    comparison = fluxes[:, 1:]
    valid = np.isfinite(comparison) & (comparison > 0)
    coverage = np.mean(valid & initial_good[:, None], axis=0) / max(np.mean(initial_good), 1e-6)
    normalizers = np.nanmedian(np.where(valid & initial_good[:, None], comparison, np.nan), axis=0)
    normalized_comparison = comparison / normalizers
    scatters = np.full(comparison.shape[1], np.nan)
    time_index = np.linspace(-1, 1, len(comparison))
    for column in range(comparison.shape[1]):
        if coverage[column] < 0.90:
            continue
        others = np.delete(normalized_comparison, column, axis=1)
        relative = normalized_comparison[:, column] / np.nanmedian(others, axis=1)
        mask = initial_good & np.isfinite(relative) & (relative > 0)
        if mask.sum() < 20:
            continue
        coefficients = np.polyfit(time_index[mask], relative[mask], deg=2)
        residual = relative[mask] / np.polyval(coefficients, time_index[mask])
        scatters[column] = robust_mad(residual)
    finite = np.isfinite(scatters)
    if finite.sum() < 6:
        raise RuntimeError("Too few comparison stars produced usable photometry")
    center = np.nanmedian(scatters[finite])
    spread = robust_mad(scatters[finite])
    threshold = max(center + 2.5 * spread, 0.006)
    selected = finite & (scatters <= threshold)
    if selected.sum() < 8:
        best = np.argsort(np.where(finite, scatters, np.inf))[: min(10, finite.sum())]
        selected[:] = False
        selected[best] = True
    return selected, scatters


def plot_field(first_frame: np.ndarray, output: Path) -> None:
    labels = list(STAR_POSITIONS)
    positions = np.array([STAR_POSITIONS[label] for label in labels])
    interval = ZScaleInterval(contrast=0.18)
    vmin, vmax = interval.get_limits(first_frame)
    norm = ImageNormalize(vmin=vmin, vmax=vmax, stretch=AsinhStretch(0.08), clip=True)
    fig, ax = plt.subplots(figsize=(10, 10))
    ax.imshow(first_frame, origin="lower", cmap="gray", norm=norm, rasterized=True)
    for label, (x, y) in zip(labels, positions):
        color = "#ff4d4d" if label == "T1" else "#36d1dc"
        radius = 30 if label == "T1" else 22
        ax.add_patch(plt.Circle((x, y), radius, fill=False, color=color, linewidth=1.5))
        ax.text(x + radius + 8, y + radius + 8, label, color=color, fontsize=8, weight="bold")
    ax.set_title("TOI-3718.01 practice field - GMU 0.8 m telescope", fontsize=16, pad=12)
    ax.set_xlabel("CCD x pixel")
    ax.set_ylabel("CCD y pixel")
    ax.text(
        0.015,
        0.018,
        "Local night 2023-10-12 | R filter | 65 s\nT1: TOI-3718.01; cyan: candidate comparison stars",
        transform=ax.transAxes,
        color="white",
        fontsize=10,
        bbox={"facecolor": "black", "alpha": 0.70, "edgecolor": "white", "boxstyle": "round,pad=0.5"},
    )
    fig.tight_layout()
    fig.savefig(output, dpi=220)
    plt.close(fig)


def plot_seeing_profile(first_frame: np.ndarray, output: Path) -> float:
    x, y = STAR_POSITIONS["T1"]
    half = 55
    x0, y0 = int(round(x)), int(round(y))
    cut = first_frame[y0 - half : y0 + half + 1, x0 - half : x0 + half + 1]
    yy, xx = np.mgrid[-half : half + 1, -half : half + 1]
    radius = np.hypot(xx - (x - x0), yy - (y - y0))
    background = np.nanmedian(cut[(radius >= ANNULUS_INNER) & (radius <= ANNULUS_OUTER)])
    signal = cut - background
    bins = np.arange(0, 50.5, 0.5)
    centers = (bins[:-1] + bins[1:]) / 2
    profile = np.array([np.nanmedian(signal[(radius >= left) & (radius < right)]) for left, right in zip(bins[:-1], bins[1:])])
    profile /= np.nanmax(profile)
    below = np.where((centers > 1) & (profile <= 0.5))[0]
    half_max_radius = float(centers[below[0]]) if len(below) else np.nan
    fwhm = 2 * half_max_radius

    fig, (ax_image, ax_profile) = plt.subplots(1, 2, figsize=(11.5, 5.2))
    vmin, vmax = ZScaleInterval(contrast=0.25).get_limits(cut)
    ax_image.imshow(cut, origin="lower", cmap="gray", norm=ImageNormalize(vmin=vmin, vmax=vmax, stretch=AsinhStretch(0.1)))
    center = (half + (x - x0), half + (y - y0))
    for value, color, style in [
        (APERTURE_RADIUS, "#ff4d4d", "-"),
        (ANNULUS_INNER, "#36d1dc", "--"),
        (ANNULUS_OUTER, "#36d1dc", "--"),
    ]:
        ax_image.add_patch(plt.Circle(center, value, fill=False, color=color, linestyle=style, linewidth=1.5))
    ax_image.set_title("Target aperture and sky annulus")
    ax_image.set_xlabel("Cutout x pixel")
    ax_image.set_ylabel("Cutout y pixel")

    ax_profile.plot(centers, profile, color="#1f77b4", linewidth=2)
    ax_profile.axhline(0.5, color="0.45", linestyle=":", label="Half maximum")
    ax_profile.axvline(
        APERTURE_RADIUS,
        color="#ff4d4d",
        linestyle="--",
        label=f"Aperture r = {APERTURE_RADIUS:.0f} px",
    )
    ax_profile.axvspan(ANNULUS_INNER, ANNULUS_OUTER, color="#36d1dc", alpha=0.18, label="Sky annulus")
    ax_profile.set(
        xlabel="Radius (pixels)",
        ylabel="Normalized median signal",
        xlim=(0, 48),
        ylim=(-0.08, 1.08),
    )
    ax_profile.set_title(f"Radial seeing profile (FWHM about {fwhm:.1f} px)")
    ax_profile.legend(fontsize=9)
    fig.suptitle("TOI-3718.01 practice aperture selection", fontsize=15)
    fig.tight_layout()
    fig.savefig(output, dpi=220)
    plt.close(fig)
    return fwhm


def plot_light_curve(table: pd.DataFrame, ingress: float, egress: float, selected_count: int, output: Path) -> None:
    good = table["used"]
    x = table["bjd_tdb"] - 2460000.0
    fig, (ax, diagnostic) = plt.subplots(2, 1, figsize=(11.5, 8.0), sharex=True, gridspec_kw={"height_ratios": [3.3, 1]})
    ax.scatter(x[~good], table.loc[~good, "normalized_flux"], s=15, color="0.72", alpha=0.8, label="excluded quality points")
    ax.scatter(x[good], table.loc[good, "normalized_flux"], s=15, color="#1769aa", alpha=0.58, label="differential photometry")

    good_indices = np.flatnonzero(good.to_numpy())
    bin_size = 5
    bx, by = [], []
    for start in range(0, len(good_indices), bin_size):
        indices = good_indices[start : start + bin_size]
        if len(indices) >= 3:
            bx.append(float(np.nanmedian(x.iloc[indices])))
            by.append(float(np.nanmedian(table["normalized_flux"].iloc[indices])))
    ax.plot(bx, by, color="#e4572e", marker="o", markersize=3.5, linewidth=1.4, label="5-point median")

    ingress_x = ingress - 2460000.0
    egress_x = egress - 2460000.0
    ax.axvspan(ingress_x, egress_x, color="#9b7ede", alpha=0.14, label="propagated spreadsheet window")
    ax.axvline(ingress_x, color="#7559a6", linestyle="--", linewidth=1)
    ax.axvline(egress_x, color="#7559a6", linestyle="--", linewidth=1)
    ax.axhline(1.0, color="0.3", linestyle=":", linewidth=1)
    ax.axhline(1.0 - EXPECTED_DEPTH, color="#b23a48", linestyle=":", linewidth=1.3, label="TESS depth = 3.77 ppt")
    finite_flux = table.loc[good, "normalized_flux"].to_numpy()
    low, high = np.nanpercentile(finite_flux, [0.5, 99.5])
    margin = max(0.003, (high - low) * 0.25)
    ax.set_ylim(low - margin, high + margin)
    ax.set_ylabel("Normalized differential flux")
    ax.set_title("Preliminary TOI-3718.01 practice light curve", fontsize=16)
    ax.text(
        0.01,
        0.02,
        f"GMU 0.8 m | R | 65 s | r={APERTURE_RADIUS:.0f} px | {selected_count} comparison stars\n"
        "Practice data from the mismatched 2023-10-12 folder; not TOI-6241.01",
        transform=ax.transAxes,
        fontsize=9.5,
        bbox={"facecolor": "white", "alpha": 0.88, "edgecolor": "0.75", "boxstyle": "round,pad=0.45"},
    )
    ax.legend(loc="best", fontsize=8.5, ncol=2)

    diagnostic.plot(x, table["comparison_sum"] / np.nanmedian(table.loc[good, "comparison_sum"]), color="#2a9d8f", linewidth=1.0)
    diagnostic.axhline(1.0, color="0.4", linestyle=":")
    diagnostic.set_ylabel("Comparison\nensemble")
    diagnostic.set_xlabel("BJD_TDB - 2,460,000")
    diagnostic.grid(alpha=0.18)
    fig.tight_layout()
    fig.savefig(output, dpi=230)
    plt.close(fig)


def plot_diagnostics(table: pd.DataFrame, comparison_names: list[str], scatters: np.ndarray, selected: np.ndarray, output: Path) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(11.5, 7.5))
    x = table["bjd_tdb"] - 2460000.0
    axes[0, 0].plot(x, table["background"], color="#6a4c93", linewidth=1)
    axes[0, 0].set(title="Median local sky background", ylabel="ADU")
    axes[0, 1].plot(x, table["shift_x"], label="x shift", linewidth=1)
    axes[0, 1].plot(x, table["shift_y"], label="y shift", linewidth=1)
    axes[0, 1].set(title="Tracked field motion", ylabel="Pixels")
    axes[0, 1].legend()
    colors = ["#2a9d8f" if use else "0.72" for use in selected]
    axes[1, 0].bar(comparison_names, scatters * 1000, color=colors)
    axes[1, 0].tick_params(axis="x", labelrotation=55)
    axes[1, 0].set(title="Comparison-star stability", ylabel="Robust scatter (ppt)")
    axes[1, 1].scatter(x, table["centroid_rms"], c=np.where(table["used"], "#1769aa", "#d1495b"), s=12, alpha=0.7)
    axes[1, 1].set(title="Centroid consistency", ylabel="RMS among stars (pixels)")
    for ax in axes.flat:
        ax.grid(alpha=0.18)
        if ax in axes[1, :]:
            ax.set_xlabel("BJD_TDB - 2,460,000" if ax is axes[1, 1] else "Comparison star")
    fig.suptitle("TOI-3718.01 practice-reduction diagnostics", fontsize=15)
    fig.tight_layout()
    fig.savefig(output, dpi=220)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("data_dir", type=Path, help="Folder containing science/, darks/, and flats/")
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    science_paths = sorted((args.data_dir / "science").glob("*.fit"), key=numbered_path_key)
    dark_paths = sorted((args.data_dir / "darks").glob("*.fit"), key=numbered_path_key)
    flat_paths = sorted((args.data_dir / "flats").glob("*.fit"), key=numbered_path_key)
    if len(science_paths) < 20:
        raise RuntimeError(f"Only {len(science_paths)} science frames were found")

    print(f"Science={len(science_paths)}, darks={len(dark_paths)}, flats={len(flat_paths)}")
    master_darks, master_flat = build_calibrations(dark_paths, flat_paths)
    science_exposure = float(fits.getheader(science_paths[0])["EXPTIME"])
    science_dark_exposure = min(master_darks, key=lambda value: abs(value - science_exposure))
    master_science_dark = master_darks[science_dark_exposure]

    labels = list(STAR_POSITIONS)
    base_positions = np.array([STAR_POSITIONS[label] for label in labels], dtype=float)
    rows: list[dict[str, object]] = []
    all_fluxes: list[np.ndarray] = []
    first_frame: np.ndarray | None = None
    registration_reference: np.ndarray | None = None

    iers.conf.auto_download = False
    sigma_clip = SigmaClip(sigma=3.0, maxiters=5)
    for index, path in enumerate(science_paths):
        data, header = fits.getdata(path, header=True)
        calibrated = (data.astype(np.float32) - master_science_dark) / master_flat
        if first_frame is None:
            first_frame = calibrated.copy()
            registration_reference = registration_image(data.astype(np.float32))

        header_shift = pointing_shift_from_header(header)
        assert registration_reference is not None
        registered_shift = (
            np.zeros(2, dtype=float)
            if index == 0
            else image_registration_shift(registration_reference, data.astype(np.float32))
        )
        estimates = []
        for x0, y0 in base_positions:
            cx, cy = local_centroid(
                calibrated,
                x0 + registered_shift[0],
                y0 + registered_shift[1],
                half_width=14,
            )
            if np.isfinite(cx) and np.isfinite(cy):
                estimates.append((cx - x0, cy - y0))
        estimates_array = np.asarray(estimates, dtype=float)
        shift = registered_shift.copy()
        if len(estimates_array) >= 4:
            provisional = np.nanmedian(estimates_array, axis=0)
            distances = np.hypot(estimates_array[:, 0] - provisional[0], estimates_array[:, 1] - provisional[1])
            limit = max(3.0, np.nanmedian(distances) + 2.5 * robust_mad(distances))
            stable_estimates = estimates_array[distances < limit]
            if len(stable_estimates) >= 4:
                shift = np.nanmedian(stable_estimates, axis=0)

        positions = []
        centroid_valid = []
        for x0, y0 in base_positions:
            predicted_x = x0 + shift[0]
            predicted_y = y0 + shift[1]
            cx, cy = local_centroid(calibrated, predicted_x, predicted_y, half_width=12)
            valid_centroid = np.isfinite(cx) and np.isfinite(cy)
            centroid_valid.append(valid_centroid)
            positions.append((cx, cy) if valid_centroid else (predicted_x, predicted_y))
        positions_array = np.asarray(positions)
        individual_shift = positions_array - base_positions
        centroid_valid_array = np.asarray(centroid_valid, dtype=bool)
        centroid_rms = (
            float(np.sqrt(np.nanmean(np.sum((individual_shift[centroid_valid_array] - shift) ** 2, axis=1))))
            if centroid_valid_array.any()
            else np.nan
        )

        apertures = CircularAperture(positions_array, r=APERTURE_RADIUS)
        annuli = CircularAnnulus(positions_array, r_in=ANNULUS_INNER, r_out=ANNULUS_OUTER)
        sums = np.asarray(aperture_photometry(calibrated, apertures, method="exact")["aperture_sum"], dtype=float)
        sky_stats = ApertureStats(calibrated, annuli, sigma_clip=sigma_clip)
        backgrounds = np.asarray(sky_stats.median, dtype=float)
        fluxes = sums - backgrounds * apertures.area
        all_fluxes.append(fluxes)

        exposure = float(header["EXPTIME"])
        jd = float(header["JD"])
        midpoint_time = Time(jd + exposure / 2.0 / 86400.0, format="jd", scale="utc")
        rows.append(
            {
                "frame": path.name,
                "utc_midpoint": midpoint_time.isot,
                "jd_utc": midpoint_time.jd,
                "bjd_tdb": compute_bjd_tdb(jd, exposure),
                "shift_x": shift[0],
                "shift_y": shift[1],
                "centroid_rms": centroid_rms,
                "background": float(np.nanmedian(backgrounds)),
                "target_x": positions_array[0, 0],
                "target_y": positions_array[0, 1],
                "header_shift_x": header_shift[0],
                "header_shift_y": header_shift[1],
                "registration_shift_x": registered_shift[0],
                "registration_shift_y": registered_shift[1],
            }
        )
        if (index + 1) % 25 == 0 or index + 1 == len(science_paths):
            print(f"Photometry {index + 1}/{len(science_paths)}")

    assert first_frame is not None
    flux_array = np.asarray(all_fluxes)
    table = pd.DataFrame(rows)
    for column, label in enumerate(labels):
        table[f"flux_{label}"] = flux_array[:, column]

    target_finite = np.isfinite(flux_array[:, 0]) & (flux_array[:, 0] > 0)
    comparison_valid = np.isfinite(flux_array[:, 1:]) & (flux_array[:, 1:] > 0)
    enough_comparisons = comparison_valid.sum(axis=1) >= 6
    comparison_levels = np.nanmedian(np.where(comparison_valid, flux_array[:, 1:], np.nan), axis=0)
    transparency = np.nanmedian(
        np.where(comparison_valid, flux_array[:, 1:] / comparison_levels, np.nan),
        axis=1,
    )
    basic_finite = target_finite & enough_comparisons & np.isfinite(transparency)
    background_center = np.nanmedian(table.loc[basic_finite, "background"])
    background_mad = robust_mad(table.loc[basic_finite, "background"].to_numpy())
    target_inside = (
        (table["target_x"].to_numpy() > ANNULUS_OUTER + 2)
        & (table["target_x"].to_numpy() < 4096 - ANNULUS_OUTER - 2)
        & (table["target_y"].to_numpy() > ANNULUS_OUTER + 2)
        & (table["target_y"].to_numpy() < 4096 - ANNULUS_OUTER - 2)
    )
    initial_good = (
        basic_finite
        & target_inside
        & (transparency > 0.45)
        & (transparency < 1.7)
        & (table["centroid_rms"].to_numpy() < 4.0)
        & (table["background"].to_numpy() < background_center + 7.0 * max(background_mad, 1.0))
    )

    print(
        "Quality diagnostics:",
        {
            "target_finite": int(target_finite.sum()),
            "enough_comparisons": int(enough_comparisons.sum()),
            "target_inside": int(target_inside.sum()),
            "centroid_rms_percentiles": np.nanpercentile(table["centroid_rms"], [5, 50, 95]).tolist(),
            "transparency_percentiles": np.nanpercentile(transparency, [5, 50, 95]).tolist(),
            "initial_good": int(initial_good.sum()),
        },
    )
    np.savez_compressed(
        args.output_dir / "intermediate_aperture_fluxes.npz",
        fluxes=flux_array,
        initial_good=initial_good,
        labels=np.asarray(labels),
    )
    table.to_csv(args.output_dir / "intermediate_tracking.csv", index=False)

    selected, scatters = select_comparison_stars(flux_array, initial_good)
    selected_fluxes = flux_array[:, 1:][:, selected]
    comparison_sum = np.nansum(selected_fluxes, axis=1)
    relative_flux = flux_array[:, 0] / comparison_sum
    ingress, egress = expected_window_bjd()
    post_transit = initial_good & (table["bjd_tdb"].to_numpy() > egress + 0.10 / 24.0)
    normalization_mask = post_transit if post_transit.sum() >= 12 else initial_good & (table.index >= int(0.8 * len(table)))
    normalization = float(np.nanmedian(relative_flux[normalization_mask]))
    normalized_flux = relative_flux / normalization

    complete_ensemble = np.all(np.isfinite(selected_fluxes) & (selected_fluxes > 0), axis=1)
    used = initial_good & complete_ensemble
    # Remove only isolated extreme points; preserve any sustained transit-like signal.
    isolated_center = pd.Series(normalized_flux).rolling(9, center=True, min_periods=4).median().to_numpy()
    residual = normalized_flux - isolated_center
    residual_scatter = robust_mad(residual[initial_good & np.isfinite(residual)])
    used &= ~((np.abs(residual) > max(6 * residual_scatter, 0.025)) & np.isfinite(residual))

    table["comparison_sum"] = comparison_sum
    table["relative_flux"] = relative_flux
    table["normalized_flux"] = normalized_flux
    table["used"] = used
    table.to_csv(args.output_dir / "toi3718_practice_measurements.csv", index=False)

    in_transit = used & (table["bjd_tdb"].to_numpy() >= ingress) & (table["bjd_tdb"].to_numpy() <= egress)
    post = used & (table["bjd_tdb"].to_numpy() > egress + 0.10 / 24.0)
    apparent_depth = float(np.nanmedian(normalized_flux[post]) - np.nanmedian(normalized_flux[in_transit])) if post.any() and in_transit.any() else np.nan
    post_scatter = robust_mad(normalized_flux[post]) if post.any() else np.nan

    plot_field(first_frame, args.output_dir / "01_toi3718_field_identification.png")
    fwhm = plot_seeing_profile(first_frame, args.output_dir / "02_toi3718_seeing_profile.png")
    plot_light_curve(table, ingress, egress, int(selected.sum()), args.output_dir / "03_toi3718_practice_light_curve.png")
    comparison_names = labels[1:]
    plot_diagnostics(table, comparison_names, scatters, selected, args.output_dir / "04_toi3718_reduction_diagnostics.png")

    summary = {
        "classification": "practice reduction; not TOI-6241.01",
        "folder_date_local": "2023-10-12",
        "observed_target_from_headers": "TOI 3718.01",
        "science_frames": len(science_paths),
        "dark_frames": len(dark_paths),
        "flat_frames": len(flat_paths),
        "filter": str(fits.getheader(science_paths[0]).get("FILTER", "R")).strip(),
        "exposure_seconds": science_exposure,
        "aperture_radius_pixels": APERTURE_RADIUS,
        "annulus_pixels": [ANNULUS_INNER, ANNULUS_OUTER],
        "seeing_fwhm_pixels_approx": fwhm,
        "comparison_stars_selected": [name for name, use in zip(labels[1:], selected) if use],
        "frames_used": int(used.sum()),
        "expected_depth_ppt": EXPECTED_DEPTH * 1000,
        "expected_duration_hours": EXPECTED_DURATION_HOURS,
        "preliminary_apparent_depth_ppt": apparent_depth * 1000,
        "post_window_scatter_ppt": post_scatter * 1000,
        "interpretation": "Preliminary practice result only; mentor review is required before scientific interpretation.",
    }
    (args.output_dir / "practice_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
