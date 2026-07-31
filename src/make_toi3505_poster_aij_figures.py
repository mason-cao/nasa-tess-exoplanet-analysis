"""AstroImageJ-style figures for the second poster layout.

The second poster shows the same products the rest of the Schar cohort shows:
a seeing profile, the finished light curve, the comparison-star apertures, the
nearby-star field, and a delta-magnitude versus scatter plot. Two of those are
regenerated here rather than reused:

* The saved AstroImageJ seeing profile in ``outputs/toi3505_seeing`` was taken
  during the 35-pixel trial, so it disagrees with the 25-pixel aperture the
  analysis actually adopted. It is rebuilt from the plate-solved image at the
  adopted settings, in the AstroImageJ visual style.
* AstroImageJ's own nearby-star plot was never exported, so the delta-magnitude
  versus scatter figure is drawn from the measured neighbour table.

The comparison-star frame is the real AstroImageJ screenshot, cropped to the
image area so the window chrome does not eat poster space.
"""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from PIL import Image

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs" / "toi3505_poster"
OUT.mkdir(parents=True, exist_ok=True)

# Adopted photometry settings, from outputs/toi3505_final_candidate.
SOURCE_RADIUS = 25.0
BACK_INNER = 70.0
BACK_OUTER = 139.0
PIXEL_SCALE = 0.3621236516728507
OBSERVED_DEPTH_PPT = 2.910
# TFOP SG1 nominal nearby-star clearance radius.
NEB_RADIUS_ARCSEC = 150.0


def figure_seeing_profile(path: Path) -> dict:
    """Radial profile of the target drawn the way AstroImageJ draws it."""
    from astropy.io import fits

    frame = (
        ROOT
        / "data"
        / "ground"
        / "toi3505"
        / "plate_solved"
        / "TOI_3505.01_50.000s_R-0001_wcs.fits"
    )
    with fits.open(frame) as handle:
        image = handle[0].data.astype(float)

    settings = json.loads(
        (ROOT / "outputs" / "toi3505_seeing" / "aperture_settings.json").read_text()
    )
    centre_x = settings["astroimagej_center_x_fits"] - 1.0
    centre_y = settings["astroimagej_center_y_fits"] - 1.0
    fwhm = settings["fwhm_pixels"]
    hwhm = fwhm / 2.0

    half = int(BACK_OUTER) + 6
    x0, x1 = int(centre_x) - half, int(centre_x) + half
    y0, y1 = int(centre_y) - half, int(centre_y) + half
    cutout = image[y0:y1, x0:x1]

    yy, xx = np.mgrid[y0:y1, x0:x1]
    radius = np.hypot(xx - centre_x, yy - centre_y).ravel()
    value = cutout.ravel()
    keep = radius <= BACK_OUTER + 4
    radius, value = radius[keep], value[keep]

    # Mean profile, the magenta curve in the AstroImageJ window.
    edges = np.arange(0, BACK_OUTER + 4, 1.5)
    centres, means = [], []
    index = np.digitize(radius, edges) - 1
    for slot in range(len(edges) - 1):
        selected = value[index == slot]
        if selected.size:
            centres.append(0.5 * (edges[slot] + edges[slot + 1]))
            means.append(float(np.mean(selected)))

    figure, axis = plt.subplots(figsize=(7.6, 6.4), dpi=220)
    axis.set_facecolor("white")
    axis.grid(True, color="#D8D8D8", linewidth=0.8)
    axis.set_axisbelow(True)
    for side in axis.spines.values():
        side.set_color("black")

    axis.plot(
        radius,
        value,
        linestyle="none",
        marker="s",
        markersize=2.0,
        markerfacecolor="none",
        markeredgecolor="#0000FF",
        markeredgewidth=0.6,
        zorder=2,
    )
    axis.plot(centres, means, color="#FF00FF", linewidth=2.6, zorder=3)

    top = float(np.nanpercentile(value, 99.9))
    floor = float(np.nanmedian(value[radius > BACK_INNER]))
    span = top - floor
    low = floor - 0.55 * span
    high = top + 0.10 * span
    bracket = high - 0.06 * (high - low)

    def bracket_box(left: float, right: float, label: str) -> None:
        axis.plot([left, left], [low, bracket], color="#CC0000", linewidth=1.1, zorder=4)
        axis.plot([right, right], [low, bracket], color="#CC0000", linewidth=1.1, zorder=4)
        axis.plot([left, right], [bracket, bracket], color="#CC0000", linewidth=1.1, zorder=4)
        axis.annotate(
            label,
            xy=(0.5 * (left + right), bracket),
            xytext=(0, 6),
            textcoords="offset points",
            ha="center",
            fontsize=15,
            color="#CC0000",
        )

    bracket_box(0.0, SOURCE_RADIUS, "SOURCE")
    bracket_box(BACK_INNER, BACK_OUTER, "BACKGROUND")
    axis.axvline(hwhm, color="#00A000", linewidth=1.3, linestyle="--", zorder=4)

    for value_x, text, colour in (
        (hwhm, f"HWHM\n{hwhm:.2f}", "#00A000"),
        (SOURCE_RADIUS, f"Radius\n{SOURCE_RADIUS:.2f}", "#CC0000"),
        (BACK_INNER, f"Back>\n{BACK_INNER:.2f}", "#CC0000"),
        (BACK_OUTER, f"<Back\n{BACK_OUTER:.2f}", "#CC0000"),
    ):
        axis.annotate(
            text,
            xy=(value_x, low),
            xytext=(0, 6),
            textcoords="offset points",
            ha="center",
            va="bottom",
            fontsize=14,
            color=colour,
        )

    axis.set_xlim(-3, BACK_OUTER + 6)
    axis.set_ylim(low, high)
    axis.set_xlabel("Radius [pixels]", fontsize=16)
    axis.set_ylabel("ADU", fontsize=16)
    axis.tick_params(labelsize=13)
    axis.set_title(
        f"TOI-3505.01, first image\n"
        f"FITS center: ({centre_x + 1:.2f}, {centre_y + 1:.2f})\n"
        f"FWHM: {fwhm:.2f} pixels = {fwhm * PIXEL_SCALE:.2f} arcsec",
        fontsize=15,
        color="black",
    )

    figure.tight_layout()
    figure.savefig(path, facecolor="white")
    plt.close(figure)
    return {
        "fwhm_pixels": fwhm,
        "fwhm_arcsec": fwhm * PIXEL_SCALE,
        "source_radius_pixels": SOURCE_RADIUS,
        "background_pixels": [BACK_INNER, BACK_OUTER],
    }


def figure_dmag_rms(path: Path) -> dict:
    """Delta magnitude against measured scatter, the standard nearby-star plot."""
    source = (
        ROOT / "outputs" / "toi3505_ground_checks" / "nearby_star_image_measurements.csv"
    )
    with source.open() as handle:
        rows = list(csv.DictReader(handle))

    def optional(row: dict, key: str) -> float:
        """Blank means the star had no usable measurement on the night."""
        value = row[key]
        return float(value) if value else float("nan")

    dmag = np.array([float(r["delta_tmag"]) for r in rows])
    rms = np.array([optional(r, "night_robust_scatter_ppt") for r in rows])
    # The clearance test uses the depth limit reached over the whole scheduled
    # window, not the frame-to-frame scatter. Binning roughly 90 in-window
    # points buys about a factor of nine in depth, which is why stars can sit
    # above the required-depth line and still be cleared.
    limit = np.array(
        [3.0 * optional(r, "historical_window_depth_error_ppt") for r in rows]
    )
    cleared = np.array([r["transit_relevant_clearance"] == "True" for r in rows])

    grid = np.linspace(0, max(dmag.max(), 7.0) + 0.3, 400)
    required = OBSERVED_DEPTH_PPT * 10 ** (0.4 * grid)

    figure, axis = plt.subplots(figsize=(7.8, 6.2), dpi=220)
    axis.set_facecolor("white")
    axis.grid(True, color="#D8D8D8", linewidth=0.8)
    axis.set_axisbelow(True)
    for side in axis.spines.values():
        side.set_color("black")

    axis.plot(
        grid,
        required,
        color="#FF00FF",
        linewidth=2.6,
        zorder=3,
        label="depth needed to fake the signal",
    )
    axis.plot(
        dmag,
        rms,
        linestyle="none",
        marker="o",
        markersize=7,
        markerfacecolor="#CC0000",
        markeredgecolor="#7A0000",
        zorder=4,
        label="RMS of each star, single frame",
    )
    good = ~np.isnan(limit)
    axis.plot(
        dmag[good],
        limit[good],
        linestyle="none",
        marker="o",
        markersize=7,
        markerfacecolor="#00A000",
        markeredgecolor="#005500",
        zorder=5,
        label="3$\\sigma$ depth limit over the window",
    )
    for x, high, low in zip(dmag[good], rms[good], limit[good]):
        axis.plot([x, x], [low, high], color="#AAAAAA", linewidth=0.9, zorder=2)

    axis.set_xlim(-0.3, max(dmag.max(), 7.0) + 0.3)
    axis.set_ylim(1, 6000)
    axis.set_yscale("log")
    axis.set_xlabel("dmag relative to TOI-3505", fontsize=16)
    axis.set_ylabel("Depth (ppt)", fontsize=16)
    axis.tick_params(labelsize=13)
    axis.set_title(
        f"Stars within {NEB_RADIUS_ARCSEC / 60.0:.1f}' of TOI-3505.01: "
        f"{int(cleared.sum())} of {len(rows)} cleared",
        fontsize=16,
        color="black",
    )
    axis.legend(loc="upper left", fontsize=12.5, framealpha=0.96)

    figure.tight_layout()
    figure.savefig(path, facecolor="white")
    plt.close(figure)
    return {
        "stars": len(rows),
        "cleared": int(cleared.sum()),
        "not_cleared": int((~cleared).sum()),
        "observed_depth_ppt": OBSERVED_DEPTH_PPT,
    }


def figure_neb_field(path: Path) -> dict:
    """The 2.5 arcminute nearby-star screen drawn on our own image.

    TFOP SG1 asks for clearance inside a nominal 2.5 arcminute radius, which is
    wide enough to cover the TESS point response function rather than just our
    photometric aperture.  Every catalog star in that circle that is bright
    enough to fake the signal is drawn, colour-coded by whether our own light
    curve for it rules the eclipse out.
    """
    from astropy.io import fits

    frame = (
        ROOT
        / "data"
        / "ground"
        / "toi3505"
        / "plate_solved"
        / "TOI_3505.01_50.000s_R-0001_wcs.fits"
    )
    with fits.open(frame) as handle:
        image = handle[0].data.astype(float)

    candidates = {}
    with (
        ROOT / "outputs" / "toi3505_ground_checks" / "nearby_star_catalog_candidates.csv"
    ).open() as handle:
        for row in csv.DictReader(handle):
            candidates[int(row["ID"])] = (
                float(row["x_zero_indexed"]),
                float(row["y_zero_indexed"]),
            )
    with (
        ROOT / "outputs" / "toi3505_ground_checks" / "nearby_star_image_measurements.csv"
    ).open() as handle:
        measured = list(csv.DictReader(handle))

    settings = json.loads(
        (ROOT / "outputs" / "toi3505_seeing" / "aperture_settings.json").read_text()
    )
    centre_x = settings["astroimagej_center_x_fits"] - 1.0
    centre_y = settings["astroimagej_center_y_fits"] - 1.0
    radius_pixels = NEB_RADIUS_ARCSEC / PIXEL_SCALE

    half = int(radius_pixels) + 18
    x0, x1 = int(centre_x) - half, int(centre_x) + half
    y0, y1 = int(centre_y) - half, int(centre_y) + half
    cutout = image[y0:y1, x0:x1]
    low, high = np.nanpercentile(cutout, [25, 99.5])

    figure, axis = plt.subplots(figsize=(7.4, 7.4), dpi=220)
    axis.imshow(
        cutout,
        origin="lower",
        cmap="gray",
        vmin=low,
        vmax=high,
        extent=(x0, x1, y0, y1),
        interpolation="nearest",
    )

    groups = {"cleared": [], "not_cleared": [], "blended": []}
    for row in measured:
        tic = int(row["tic_id"])
        if tic not in candidates:
            continue
        if row["target_aperture_overlap"] == "True":
            groups["blended"].append(candidates[tic])
        elif row["transit_relevant_clearance"] == "True":
            groups["cleared"].append(candidates[tic])
        else:
            groups["not_cleared"].append(candidates[tic])

    style = {
        "cleared": ("#00C000", "o", 7.0, f"cleared ({len(groups['cleared'])})"),
        "not_cleared": (
            "#FF7000",
            "o",
            7.0,
            f"not cleared ({len(groups['not_cleared'])})",
        ),
        "blended": ("#FF00FF", "s", 9.0, f"blended with T1 ({len(groups['blended'])})"),
    }
    for key, points in groups.items():
        if not points:
            continue
        colour, marker, size, label = style[key]
        axis.plot(
            [p[0] for p in points],
            [p[1] for p in points],
            linestyle="none",
            marker=marker,
            markersize=size,
            markerfacecolor="none",
            markeredgecolor=colour,
            markeredgewidth=1.15,
            label=label,
        )

    circle = plt.Circle(
        (centre_x, centre_y),
        radius_pixels,
        fill=False,
        color="#CC0000",
        linewidth=1.5,
        linestyle="--",
    )
    axis.add_patch(circle)
    axis.plot(
        [centre_x],
        [centre_y],
        marker="+",
        markersize=17,
        markeredgewidth=2.4,
        color="#00CFFF",
        linestyle="none",
        label="T1",
    )

    axis.set_xlim(x0, x1)
    axis.set_ylim(y0, y1)
    axis.set_xticks([])
    axis.set_yticks([])
    axis.set_title(
        f"Stars within {NEB_RADIUS_ARCSEC / 60.0:.1f}' of TOI-3505.01: "
        f"{len(groups['cleared'])} of {len(measured)} cleared",
        fontsize=16,
        color="black",
    )
    axis.legend(
        loc="upper center",
        bbox_to_anchor=(0.5, -0.012),
        ncol=4,
        fontsize=12.5,
        frameon=False,
        handletextpad=0.35,
        columnspacing=1.1,
    )

    figure.tight_layout()
    figure.savefig(path, facecolor="white")
    plt.close(figure)
    return {
        "radius_arcsec": NEB_RADIUS_ARCSEC,
        "stars": len(measured),
        "cleared": len(groups["cleared"]),
        "not_cleared": len(groups["not_cleared"]),
        "blended": len(groups["blended"]),
    }


def figure_comparison_field(path: Path) -> dict:
    """The AstroImageJ aperture frame, cropped to the image area."""
    source = ROOT / "outputs" / "toi3505_photometry" / "06_apertures_70-139.png"
    image = Image.open(source).convert("RGB")
    # Drop the window title bar, toolbar, and the histogram strip underneath.
    box = (165, 452, 1462, 1742)
    cropped = image.crop(box)
    cropped.save(path)
    return {"source": source.name, "crop_box": box, "size": cropped.size}


def main() -> None:
    seeing = figure_seeing_profile(OUT / "v2_01_seeing_profile.png")
    field = figure_comparison_field(OUT / "v2_02_comparison_field.png")
    dmag = figure_dmag_rms(OUT / "v2_03_dmag_rms.png")
    neb = figure_neb_field(OUT / "v2_04_neb_field.png")

    summary = {
        "seeing": seeing,
        "comparison_field": field,
        "dmag_rms": dmag,
        "neb_field": neb,
    }
    (OUT / "v2_figure_summary.json").write_text(json.dumps(summary, indent=2))

    print(json.dumps(summary, indent=2))
    print(f"\nWrote four figures to {OUT}")


if __name__ == "__main__":
    main()
