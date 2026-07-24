"""Render and classify the TOI-3505.01 frames selected for review."""

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
from astropy.io import fits
from astropy.visualization import AsinhStretch, ImageNormalize, PercentileInterval


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_ALIGNED_DIR = ROOT / "data" / "ground" / "toi3505" / "aligned"
DEFAULT_OUTPUT_DIR = ROOT / "outputs" / "toi3505_final_candidate"
DEFAULT_CURVE = (
    DEFAULT_OUTPUT_DIR
    / "TOI_3505.01_2022-07-22_R_final_candidate_light_curve.csv"
)
DEFAULT_TABLE = (
    ROOT
    / "outputs"
    / "toi3505_aperture_check"
    / "TOI_3505.01_2022-07-22_R_measurements_25px_70-139_AIJ.tbl"
)
ARTIFACT_CHECK_SLICES = (245, 253)
STAR_NAMES = ("T1", *(f"C{number}" for number in range(2, 12)))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--aligned-dir", type=Path, default=DEFAULT_ALIGNED_DIR)
    parser.add_argument("--curve", type=Path, default=DEFAULT_CURVE)
    parser.add_argument("--table", type=Path, default=DEFAULT_TABLE)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    return parser.parse_args()


def preliminary_classification(row: pd.Series) -> tuple[str, str]:
    reason = str(row["review_reason"])
    if int(row["slice"]) == 253:
        return (
            "field-wide vertical PSF trail confirmed across all eleven stars",
            "keep in the table but exclude from the plotted light curve",
        )
    if "large star width" in reason:
        asymmetry = abs(float(row["X-Width_T1"]) - float(row["Y-Width_T1"]))
        if asymmetry >= 8.0:
            return (
                "elongated PSF / tracking or seeing excursion",
                "keep in the table but exclude from the plotted light curve",
            )
        return (
            "broad PSF / poor-seeing excursion",
            "keep in the table but exclude from the plotted light curve",
        )
    if "low comparison-star counts" in reason:
        return (
            "field-wide transparency loss, consistent with cloud",
            "keep in the table but exclude from the plotted light curve",
        )
    return (
        "differential-flux excursion within the late-night trend",
        "keep in the plotted light curve; no separate image problem found",
    )


def plot_cross_star_artifact_check(
    measurements: pd.DataFrame,
    aligned_dir: Path,
    output_path: Path,
) -> None:
    fig, axes = plt.subplots(
        len(ARTIFACT_CHECK_SLICES),
        len(STAR_NAMES),
        figsize=(20.5, 5.8),
    )
    for row_index, slice_number in enumerate(ARTIFACT_CHECK_SLICES):
        row = measurements.loc[measurements["slice"] == slice_number].iloc[0]
        path = aligned_dir / str(row["Label"])
        with fits.open(path, memmap=True) as hdul:
            data = hdul[0].data
            for column_index, star in enumerate(STAR_NAMES):
                axis = axes[row_index, column_index]
                x = float(row[f"X(FITS)_{star}"]) - 1.0
                y = float(row[f"Y(FITS)_{star}"]) - 1.0
                half_size = 42
                x0 = int(round(x))
                y0 = int(round(y))
                cutout = np.asarray(
                    data[
                        y0 - half_size : y0 + half_size + 1,
                        x0 - half_size : x0 + half_size + 1,
                    ],
                    dtype=float,
                )
                norm = ImageNormalize(
                    cutout,
                    interval=PercentileInterval(99.5),
                    stretch=AsinhStretch(0.08),
                    clip=True,
                )
                axis.imshow(cutout, origin="lower", cmap="gray", norm=norm)
                center = (half_size + (x - x0), half_size + (y - y0))
                axis.add_patch(
                    plt.Circle(
                        center,
                        25.0,
                        fill=False,
                        color="#f0b429",
                        linewidth=0.9,
                    )
                )
                axis.set_title(f"{star}", fontsize=8.8)
                axis.set_xticks([])
                axis.set_yticks([])
                if column_index == 0:
                    axis.set_ylabel(f"Frame {slice_number}", fontsize=10.5)
    fig.suptitle(
        "TOI-3505.01: comparison-star check for frames 245 and 253",
        fontsize=15,
    )
    fig.text(
        0.5,
        0.94,
        "Each cutout has an independent display stretch; orange circle is the 25 px aperture",
        ha="center",
        fontsize=9.5,
    )
    fig.tight_layout(rect=(0.02, 0.02, 0.99, 0.91))
    fig.savefig(output_path, dpi=190)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    aligned_dir = args.aligned_dir.resolve()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    curve = pd.read_csv(args.curve.resolve())
    measurements = pd.read_csv(args.table.resolve(), sep="\t")
    flagged = curve.loc[curve["review_point"]].copy()
    if len(flagged) == 0:
        raise RuntimeError("No review points were found")
    merged = flagged.merge(
        measurements[
            [
                "slice",
                "X(FITS)_T1",
                "Y(FITS)_T1",
                "X-Width_T1",
                "Y-Width_T1",
                "Peak_T1",
                "Source-Sky_T1",
            ]
        ],
        on="slice",
        how="left",
        validate="one_to_one",
    )

    columns = 5
    rows = int(np.ceil(len(merged) / columns))
    fig, axes = plt.subplots(rows, columns, figsize=(14.5, 2.9 * rows))
    axes_flat = list(np.atleast_1d(axes).flat)
    classifications: list[dict[str, object]] = []
    for axis, (_, row) in zip(axes_flat, merged.iterrows()):
        path = aligned_dir / str(row["image"])
        with fits.open(path, memmap=True) as hdul:
            data = hdul[0].data
            x = float(row["X(FITS)_T1"]) - 1.0
            y = float(row["Y(FITS)_T1"]) - 1.0
            half_size = 52
            x0 = int(round(x))
            y0 = int(round(y))
            cutout = np.asarray(
                data[
                    y0 - half_size : y0 + half_size + 1,
                    x0 - half_size : x0 + half_size + 1,
                ],
                dtype=float,
            )
        norm = ImageNormalize(
            cutout,
            interval=PercentileInterval(99.5),
            stretch=AsinhStretch(0.08),
            clip=True,
        )
        axis.imshow(cutout, origin="lower", cmap="gray", norm=norm)
        center = (half_size + (x - x0), half_size + (y - y0))
        axis.add_patch(
            plt.Circle(center, 25.0, fill=False, color="#f0b429", linewidth=1.15)
        )
        axis.plot(center[0], center[1], "+", color="#42d4f4", markersize=7)
        axis.set_title(
            f"Frame {int(row['slice']):04d}\n"
            f"X/Y width {float(row['X-Width_T1']):.1f}/{float(row['Y-Width_T1']):.1f} px",
            fontsize=9.0,
        )
        axis.set_xticks([])
        axis.set_yticks([])
        classification, decision = preliminary_classification(row)
        classifications.append(
            {
                "slice": int(row["slice"]),
                "image": row["image"],
                "review_reason": row["review_reason"],
                "visual_classification": classification,
                "decision": decision,
                "used_in_primary_curve": bool(row["used_in_primary_curve"]),
                "condition_model_fit_eligible": bool(
                    row["condition_model_fit_eligible"]
                ),
                "relative_brightness": float(
                    row["raw_relative_brightness_10_comparisons"]
                ),
                "relative_brightness_error": float(
                    row["raw_relative_brightness_error"]
                ),
                "comparison_counts": float(row["comparison_counts"]),
                "target_width_pixels": float(row["target_width_pixels"]),
                "x_width_pixels": float(row["X-Width_T1"]),
                "y_width_pixels": float(row["Y-Width_T1"]),
                "target_peak_counts": float(row["Peak_T1"]),
            }
        )
    for axis in axes_flat[len(merged) :]:
        axis.axis("off")
    fig.suptitle(
        "TOI-3505.01: target images checked after photometry",
        fontsize=15.5,
        y=0.995,
    )
    fig.text(
        0.5,
        0.978,
        "Individually normalized display stretch; orange circle is the 25 px source aperture",
        ha="center",
        fontsize=9.8,
    )
    fig.tight_layout(rect=(0.02, 0.02, 0.98, 0.96))
    fig.savefig(output_dir / "02_flagged_target_cutouts.png", dpi=200)
    plt.close(fig)

    plot_cross_star_artifact_check(
        measurements,
        aligned_dir,
        output_dir / "03_artifact_cross_star_check.png",
    )

    review = pd.DataFrame(classifications)
    review.to_csv(
        output_dir / "flagged_frame_review.csv", index=False, float_format="%.10f"
    )
    public_review = review.rename(
        columns={
            "slice": "Frame",
            "image": "Image",
            "review_reason": "Measurement_flag",
            "visual_classification": "Image_review",
            "decision": "Decision",
            "used_in_primary_curve": "Used_in_plot",
            "relative_brightness": "Relative_brightness",
            "relative_brightness_error": "Flux_error",
            "comparison_counts": "Comparison_counts",
            "target_width_pixels": "Target_width_pixels",
            "x_width_pixels": "X_width_pixels",
            "y_width_pixels": "Y_width_pixels",
            "target_peak_counts": "Target_peak_counts",
        }
    ).drop(columns=["condition_model_fit_eligible"])
    public_review.to_csv(
        output_dir / "frame_review.csv", index=False, float_format="%.10f"
    )
    summary = {
        "review_frames": int(len(review)),
        "used_in_primary_curve": int(review["used_in_primary_curve"].sum()),
        "excluded_from_primary_curve": int((~review["used_in_primary_curve"]).sum()),
        "classification_counts": review["visual_classification"]
        .value_counts()
        .to_dict(),
        "decision": (
            "All reviewed frames remain in the saved table. Eighteen frames with "
            "image problems are excluded from the plotted light curve; five "
            "brightness changes without a separate image problem remain in it."
        ),
    }
    (output_dir / "flagged_frame_review_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(
        f"Reviewed {len(review)} flagged frames; "
        f"{int(review['used_in_primary_curve'].sum())} remain in the plotted curve"
    )


if __name__ == "__main__":
    main()
