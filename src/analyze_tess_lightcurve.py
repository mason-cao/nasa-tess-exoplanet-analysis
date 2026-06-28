from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from fits_reader import list_hdus, read_bintable, read_image

ROOT = Path(__file__).resolve().parents[1]
LC_FILE = ROOT / "data" / "raw" / "tess2018206045859-s0001-0000000441420236-0120-s_lc.fits"
DVT_FILE = ROOT / "data" / "raw" / "tess2018206190142-s0001-s0001-0000000441420236-00366_dvt.fits"
PLOT_DIR = ROOT / "outputs" / "plots"
PROCESSED_DIR = ROOT / "data" / "processed"
PLOT_DIR.mkdir(parents=True, exist_ok=True)
PROCESSED_DIR.mkdir(parents=True, exist_ok=True)


def moving_median(x, y, width_days=0.25):
    """Small, dependency-light rolling median for plotting a smooth baseline."""
    out = np.full_like(y, np.nan, dtype=float)
    half = width_days / 2
    for i, t in enumerate(x):
        m = np.isfinite(y) & (np.abs(x - t) <= half)
        if m.sum() >= 5:
            out[i] = np.nanmedian(y[m])
    return out


def main():
    # Light curve table from the FITS binary table extension.
    lc, lc_header = read_bintable(LC_FILE, hdu_index=1)
    aperture, aperture_header = read_image(LC_FILE, hdu_index=2)

    time = lc["TIME"].astype(float)
    flux = lc["PDCSAP_FLUX"].astype(float)
    flux_err = lc["PDCSAP_FLUX_ERR"].astype(float)
    quality = lc["QUALITY"].astype(int)

    # Data Validation Time Series file contains the TCE ephemeris used here.
    # The selected signal is TCE_5 because its epoch matches the transit window
    # used in the lecture demo notebook, around TBJD 1330.4.
    dvt_hdus = list_hdus(DVT_FILE)
    tce5_header = dvt_hdus[5]["header"]
    target = lc_header.get("OBJECT", "TIC 441420236").strip()
    sector = lc_header.get("SECTOR", 1)
    camera = lc_header.get("CAMERA", "")
    ccd = lc_header.get("CCD", "")
    t0 = float(tce5_header["TEPOCH"])
    period_days = float(tce5_header["TPERIOD"])
    duration_hours = float(tce5_header["TDUR"])
    duration_days = duration_hours / 24.0
    dvt_depth_ppm = float(tce5_header["TDEPTH"])

    # Clean but transparent quality selection.
    finite = np.isfinite(time) & np.isfinite(flux)
    good = finite & (quality == 0)
    flagged_finite = finite & (quality != 0)

    median_flux = np.nanmedian(flux[good])
    normalized_flux = flux / median_flux
    normalized_err = flux_err / median_flux

    # Local polynomial detrend around the requested transit.
    dt = time - t0
    zoom = good & (np.abs(dt) < 0.35)
    oot_for_fit = zoom & (np.abs(dt) > duration_days * 0.85) & (np.abs(dt) < 0.32)
    coeffs = np.polyfit(dt[oot_for_fit], normalized_flux[oot_for_fit], deg=2)
    local_baseline = np.polyval(coeffs, dt)
    detrended = normalized_flux / local_baseline
    in_transit = zoom & (np.abs(dt) <= duration_days / 2)
    out_of_transit = zoom & (np.abs(dt) > duration_days * 1.2) & (np.abs(dt) < 0.28)
    estimated_depth_ppm = (1 - np.nanmedian(detrended[in_transit])) * 1_000_000
    local_scatter_ppm = np.nanstd(detrended[out_of_transit]) * 1_000_000

    # Save clean data for transparency/reproducibility.
    processed = PROCESSED_DIR / "au_mic_tess_sector1_clean_light_curve.csv"
    with processed.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["time_tbjd", "pdcsap_flux", "pdcsap_flux_err", "normalized_flux", "quality", "used_in_clean_mask"])
        for row in zip(time, flux, flux_err, normalized_flux, quality, good):
            writer.writerow(row)

    summary = PROCESSED_DIR / "transit_summary.csv"
    with summary.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["field", "value", "unit_or_note"])
        writer.writerows([
            ["target", target, "FITS OBJECT"],
            ["sector", sector, "TESS sector"],
            ["camera", camera, "TESS camera"],
            ["ccd", ccd, "TESS CCD"],
            ["cadences_total", len(time), "rows in light-curve FITS table"],
            ["cadences_finite", int(finite.sum()), "finite PDCSAP flux"],
            ["cadences_quality_0", int(good.sum()), "finite and quality == 0"],
            ["cadences_flagged_finite", int(flagged_finite.sum()), "finite and quality != 0"],
            ["median_pdcsap_flux", f"{median_flux:.6f}", "e-/s"],
            ["selected_transit_epoch", f"{t0:.9f}", "TBJD"],
            ["dvt_period", f"{period_days:.9f}", "days"],
            ["dvt_duration", f"{duration_hours:.6f}", "hours"],
            ["dvt_depth", f"{dvt_depth_ppm:.3f}", "ppm"],
            ["estimated_local_depth", f"{estimated_depth_ppm:.3f}", "ppm, from local polynomial detrend"],
            ["local_out_of_transit_scatter", f"{local_scatter_ppm:.3f}", "ppm"],
        ])

    # 1. Full light curve with selected transit highlighted.
    fig, ax = plt.subplots(figsize=(11, 5.6))
    ax.scatter(time[good], normalized_flux[good], s=6, alpha=0.55, label="quality = 0")
    ax.scatter(time[flagged_finite], normalized_flux[flagged_finite], s=12, alpha=0.85, label="finite but flagged")
    ax.axvspan(t0 - 0.25, t0 + 0.25, alpha=0.18, label="zoomed transit window")
    ax.axvline(t0, linestyle="--", linewidth=1.2, label=f"TCE epoch = {t0:.3f} TBJD")
    ax.set_title(f"TESS Sector {sector} Light Curve for {target}")
    ax.set_xlabel("Time (TBJD = BJD - 2457000)")
    ax.set_ylabel("PDCSAP Flux / Median")
    ax.set_ylim(0.982, 1.055)
    ax.legend(loc="upper right", frameon=True)
    fig.tight_layout()
    fig.savefig(PLOT_DIR / "01_full_light_curve_quality_and_transit_window.png", dpi=220)
    plt.close(fig)

    # 2. Required zoomed plot: same time/flux-range goal as assignment.
    fig, ax = plt.subplots(figsize=(9.5, 5.8))
    ax.errorbar(time[zoom], normalized_flux[zoom], yerr=normalized_err[zoom], fmt=".", markersize=3.2, elinewidth=0.45, alpha=0.72, label="PDCSAP flux")
    dense_t = np.linspace(-0.25, 0.25, 300)
    ax.plot(t0 + dense_t, np.polyval(coeffs, dense_t), linewidth=2.0, label="local stellar baseline")
    ax.axvline(t0, linestyle="--", linewidth=1.1, label="candidate transit center")
    ax.axvspan(t0 - duration_days / 2, t0 + duration_days / 2, alpha=0.18, label=f"DVT duration ≈ {duration_hours:.2f} hr")
    ax.set_xlim(t0 - 0.25, t0 + 0.25)
    ax.set_ylim(0.989, 1.003)
    ax.set_title(f"Required Transit Zoom: {target}, TESS Sector {sector}")
    ax.set_xlabel("Time (TBJD)")
    ax.set_ylabel("Median-normalized PDCSAP Flux")
    ax.text(0.02, 0.04, f"Estimated local depth ≈ {estimated_depth_ppm:,.0f} ppm", transform=ax.transAxes, fontsize=10,
            bbox={"boxstyle": "round", "facecolor": "white", "alpha": 0.82, "edgecolor": "0.8"})
    ax.legend(loc="lower left", frameon=True)
    fig.tight_layout()
    fig.savefig(PLOT_DIR / "02_required_transit_zoom.png", dpi=240)
    plt.close(fig)

    # 3. Detrended transit plot for a cleaner scientific view.
    fig, ax = plt.subplots(figsize=(9.5, 5.5))
    ax.scatter((time[zoom] - t0) * 24, detrended[zoom], s=12, alpha=0.75, label="locally detrended flux")
    ax.axhline(1.0, linestyle="--", linewidth=1.0)
    ax.axvline(0.0, linestyle="--", linewidth=1.0, label="mid-transit")
    ax.axvspan(-duration_hours / 2, duration_hours / 2, alpha=0.18, label="DVT duration")
    ax.set_xlim(-6, 6)
    ax.set_ylim(0.9955, 1.0025)
    ax.set_title("Local Detrend of the Candidate Transit")
    ax.set_xlabel("Hours from selected transit center")
    ax.set_ylabel("Flux / Local Polynomial Baseline")
    ax.text(0.02, 0.05, f"Depth ≈ {estimated_depth_ppm:,.0f} ppm\nOOT scatter ≈ {local_scatter_ppm:,.0f} ppm", transform=ax.transAxes, fontsize=10,
            bbox={"boxstyle": "round", "facecolor": "white", "alpha": 0.82, "edgecolor": "0.8"})
    ax.legend(loc="upper right", frameon=True)
    fig.tight_layout()
    fig.savefig(PLOT_DIR / "03_local_detrended_transit.png", dpi=240)
    plt.close(fig)

    # 4. Aperture mask from the second light-curve FITS extension.
    fig, ax = plt.subplots(figsize=(6.4, 6.0))
    im = ax.imshow(aperture, origin="lower")
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("Aperture bitmask value")
    ax.set_title(f"TESS Aperture Mask for {target}")
    ax.set_xlabel("Column pixel")
    ax.set_ylabel("Row pixel")
    ax.text(0.03, 0.97, f"NPIXSAP = {aperture_header.get('NPIXSAP', 'unknown')}", transform=ax.transAxes, va="top",
            bbox={"boxstyle": "round", "facecolor": "white", "alpha": 0.82, "edgecolor": "0.8"})
    fig.tight_layout()
    fig.savefig(PLOT_DIR / "04_aperture_mask.png", dpi=220)
    plt.close(fig)

    # 5. Quality flag compact diagnostic.
    unique_flags, counts = np.unique(quality, return_counts=True)
    fig, ax = plt.subplots(figsize=(9.2, 5.2))
    ax.bar([str(v) for v in unique_flags], counts)
    ax.set_yscale("log")
    ax.set_title("TESS QUALITY Flag Distribution")
    ax.set_xlabel("QUALITY flag value")
    ax.set_ylabel("Cadence count, log scale")
    fig.tight_layout()
    fig.savefig(PLOT_DIR / "05_quality_flag_distribution.png", dpi=220)
    plt.close(fig)

    print("Analysis complete.")
    print(f"Target: {target}")
    print(f"Saved plots to: {PLOT_DIR}")
    print(f"Selected transit center: {t0:.6f} TBJD")
    print(f"Estimated local depth: {estimated_depth_ppm:.0f} ppm")


if __name__ == "__main__":
    main()
