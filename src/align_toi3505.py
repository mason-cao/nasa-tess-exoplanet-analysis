"""Align the usable TOI-3505.01 images with whole-pixel shifts.

The Schar light-curve tutorial allows image alignment when plate solving is
skipped. It asks for at least three reference stars and recommends whole-pixel
alignment so the images are not interpolated. This script follows that route.

Frames 0282 and 0283 are left out because the earlier visual review showed
that the field is no longer usable in those images. The first 281 calibrated
images are aligned to frame 0001. A broad field match finds each approximate
shift, and the final shift is measured from the centers of reference stars.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import re
import tempfile
import time
from dataclasses import dataclass
from datetime import UTC, datetime
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
from astropy.visualization import AsinhStretch, ImageNormalize, ZScaleInterval
from scipy.fft import irfftn, rfftn
from scipy.ndimage import maximum_filter


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT_DIR = ROOT / "data" / "ground" / "toi3505" / "reduced"
DEFAULT_OUTPUT_DIR = ROOT / "data" / "ground" / "toi3505" / "aligned"
DEFAULT_DIAGNOSTICS_DIR = ROOT / "outputs" / "toi3505_alignment"
FRAME_RE = re.compile(r"-(\d{4})_out\.fits$")
EXPECTED_SHAPE = (4096, 4096)
FIRST_SEQUENCE = 1
LAST_USABLE_SEQUENCE = 281
ALIGNMENT_VERSION = "1.0.0"


@dataclass(frozen=True)
class Frame:
    sequence: int
    path: Path


@dataclass(frozen=True)
class StarMeasurement:
    x: float
    y: float
    flux: float
    fwhm: float
    peak_fraction: float


def discover_frames(
    input_dir: Path,
    *,
    first_sequence: int = FIRST_SEQUENCE,
    last_sequence: int = LAST_USABLE_SEQUENCE,
) -> list[Frame]:
    frames: list[Frame] = []
    for path in sorted(input_dir.glob("*_out.fits")):
        match = FRAME_RE.search(path.name)
        if match:
            sequence = int(match.group(1))
            if first_sequence <= sequence <= last_sequence:
                frames.append(Frame(sequence=sequence, path=path))
    expected = list(range(first_sequence, last_sequence + 1))
    found = [frame.sequence for frame in frames]
    if found != expected:
        missing = sorted(set(expected) - set(found))
        raise RuntimeError(
            f"Expected frames {first_sequence:04d}-{last_sequence:04d}; "
            f"missing {missing}"
        )
    return frames


def read_frame(path: Path) -> tuple[np.ndarray, fits.Header]:
    with fits.open(path, memmap=True, checksum=False) as hdul:
        if len(hdul) != 1 or hdul[0].data is None:
            raise RuntimeError(f"Expected one image in {path}")
        data = np.array(hdul[0].data, dtype=np.float32, copy=True)
        header = hdul[0].header.copy()
    if data.shape != EXPECTED_SHAPE:
        raise RuntimeError(f"Unexpected image shape in {path}: {data.shape}")
    if not bool(header.get("REDUCED", False)):
        raise RuntimeError(f"Input is not marked as calibrated: {path}")
    return data, header


def robust_center_and_noise(values: np.ndarray) -> tuple[float, float]:
    finite = np.asarray(values, dtype=float)
    finite = finite[np.isfinite(finite)]
    if not len(finite):
        raise RuntimeError("No finite pixels were available")
    center = float(np.median(finite))
    noise = 1.4826 * float(np.median(np.abs(finite - center)))
    if not np.isfinite(noise) or noise <= 0:
        noise = float(np.std(finite))
    if not np.isfinite(noise) or noise <= 0:
        raise RuntimeError("Could not measure the image noise")
    return center, noise


def field_signature(image: np.ndarray, *, stride: int = 8) -> np.ndarray:
    sampled = np.asarray(image[::stride, ::stride], dtype=np.float32)
    center, noise = robust_center_and_noise(sampled)
    signal = np.clip(sampled - (center + 5.0 * noise), 0.0, None)
    positive = signal[signal > 0]
    if len(positive) < 20:
        raise RuntimeError("Too few stars were found for the broad field match")
    cap = float(np.percentile(positive, 99.5))
    signal = np.minimum(signal, cap)
    signal[:2, :] = 0
    signal[-2:, :] = 0
    signal[:, :2] = 0
    signal[:, -2:] = 0
    return signal


def broad_field_shift(
    reference_fft: np.ndarray,
    moving: np.ndarray,
    *,
    signature_shape: tuple[int, int],
    stride: int = 8,
) -> tuple[int, int, float]:
    moving_signature = field_signature(moving, stride=stride)
    cross_power = reference_fft * np.conj(rfftn(moving_signature))
    cross_power /= np.maximum(np.abs(cross_power), 1e-12)
    correlation = irfftn(cross_power, s=signature_shape)
    peak_y, peak_x = np.unravel_index(np.argmax(correlation), correlation.shape)
    if peak_y > correlation.shape[0] // 2:
        peak_y -= correlation.shape[0]
    if peak_x > correlation.shape[1] // 2:
        peak_x -= correlation.shape[1]
    return (
        int(peak_y * stride),
        int(peak_x * stride),
        float(np.max(correlation)),
    )


def measure_star(
    image: np.ndarray,
    x_guess: float,
    y_guess: float,
    *,
    source_radius: int = 10,
    annulus_inner: int = 11,
    annulus_outer: int = 15,
) -> StarMeasurement | None:
    center_x = int(round(x_guess))
    center_y = int(round(y_guess))
    radius = annulus_outer
    height, width = image.shape
    if (
        center_x - radius < 0
        or center_y - radius < 0
        or center_x + radius >= width
        or center_y + radius >= height
    ):
        return None

    patch = np.asarray(
        image[
            center_y - radius : center_y + radius + 1,
            center_x - radius : center_x + radius + 1,
        ],
        dtype=float,
    )
    yy, xx = np.indices(patch.shape, dtype=float)
    local_x = xx - radius
    local_y = yy - radius
    distance = np.hypot(local_x, local_y)
    annulus = patch[(distance >= annulus_inner) & (distance <= annulus_outer)]
    finite_annulus = annulus[np.isfinite(annulus)]
    if len(finite_annulus) < 20:
        return None
    background = float(np.median(finite_annulus))
    weights = np.clip(patch - background, 0.0, None)
    weights[(distance > source_radius) | ~np.isfinite(patch)] = 0.0
    flux = float(np.sum(weights))
    if not np.isfinite(flux) or flux <= 0:
        return None

    x_offset = float(np.sum(local_x * weights) / flux)
    y_offset = float(np.sum(local_y * weights) / flux)
    x = center_x + x_offset
    y = center_y + y_offset
    x_variance = float(np.sum((local_x - x_offset) ** 2 * weights) / flux)
    y_variance = float(np.sum((local_y - y_offset) ** 2 * weights) / flux)
    fwhm = 2.355 * math.sqrt(max(0.0, 0.5 * (x_variance + y_variance)))
    peak = float(np.nanmax(patch[distance <= source_radius]) - background)
    peak_fraction = peak / flux
    return StarMeasurement(
        x=x,
        y=y,
        flux=flux,
        fwhm=fwhm,
        peak_fraction=peak_fraction,
    )


def find_reference_stars(
    image: np.ndarray,
    *,
    maximum_stars: int = 25,
    edge_margin: int = 500,
    minimum_separation: float = 220.0,
) -> list[StarMeasurement]:
    sampled = image[::8, ::8]
    center, noise = robust_center_and_noise(sampled)
    local_maximum = maximum_filter(image, size=15, mode="nearest")
    yy, xx = np.nonzero(
        (image == local_maximum)
        & (image > center + 15.0 * noise)
        & (image < 45_000.0)
    )
    peaks = image[yy, xx]
    order = np.argsort(peaks)[::-1]
    selected: list[StarMeasurement] = []
    for index in order:
        x = int(xx[index])
        y = int(yy[index])
        if (
            x < edge_margin
            or y < edge_margin
            or x >= image.shape[1] - edge_margin
            or y >= image.shape[0] - edge_margin
        ):
            continue
        star = measure_star(image, x, y)
        if star is None:
            continue
        if not (6.0 <= star.fwhm <= 13.5):
            continue
        if star.flux < 50_000.0 or star.peak_fraction >= 0.05:
            continue
        if math.hypot(star.x - x, star.y - y) > 4.5:
            continue
        if any(
            math.hypot(star.x - other.x, star.y - other.y) < minimum_separation
            for other in selected
        ):
            continue
        selected.append(star)
        if len(selected) == maximum_stars:
            break
    if len(selected) < 5:
        raise RuntimeError(
            f"Only {len(selected)} suitable reference stars were found; at least 5 are required"
        )
    return selected


def measure_shift_from_stars(
    image: np.ndarray,
    reference_stars: list[StarMeasurement],
    *,
    approximate_y_shift: float,
    approximate_x_shift: float,
) -> dict[str, float | int | list[StarMeasurement]] | None:
    moving_stars: list[StarMeasurement] = []
    reference_used: list[StarMeasurement] = []
    dx_values: list[float] = []
    dy_values: list[float] = []
    for reference in reference_stars:
        moving = measure_star(
            image,
            reference.x - approximate_x_shift,
            reference.y - approximate_y_shift,
        )
        if moving is None:
            continue
        expected_x = reference.x - approximate_x_shift
        expected_y = reference.y - approximate_y_shift
        if math.hypot(moving.x - expected_x, moving.y - expected_y) > 6.0:
            continue
        if not (5.0 <= moving.fwhm <= 15.0):
            continue
        if moving.flux < 20_000.0 or moving.peak_fraction >= 0.10:
            continue
        dx_values.append(reference.x - moving.x)
        dy_values.append(reference.y - moving.y)
        moving_stars.append(moving)
        reference_used.append(reference)

    if len(dx_values) < 5:
        return None
    dx = np.asarray(dx_values)
    dy = np.asarray(dy_values)
    median_x = float(np.median(dx))
    median_y = float(np.median(dy))
    distance = np.hypot(dx - median_x, dy - median_y)
    distance_center = float(np.median(distance))
    distance_noise = 1.4826 * float(np.median(np.abs(distance - distance_center)))
    if not np.isfinite(distance_noise):
        distance_noise = 0.0
    cutoff = max(1.5, distance_center + 4.0 * distance_noise)
    keep = distance <= cutoff
    if int(np.sum(keep)) < 5:
        return None

    final_x = float(np.median(dx[keep]))
    final_y = float(np.median(dy[keep]))
    integer_x = int(np.rint(final_x))
    integer_y = int(np.rint(final_y))
    residuals = np.hypot(dx[keep] - integer_x, dy[keep] - integer_y)
    used_moving = [star for star, use in zip(moving_stars, keep) if use]
    used_reference = [star for star, use in zip(reference_used, keep) if use]
    return {
        "x_shift": integer_x,
        "y_shift": integer_y,
        "measured_x_shift": final_x,
        "measured_y_shift": final_y,
        "stars_used": int(np.sum(keep)),
        "position_rms": float(np.sqrt(np.mean(residuals**2))),
        "median_star_width": float(np.median([star.fwhm for star in used_moving])),
        "moving_stars": used_moving,
        "reference_stars": used_reference,
    }


def choose_shift(
    image: np.ndarray,
    reference_stars: list[StarMeasurement],
    *,
    broad_y_shift: int,
    broad_x_shift: int,
    previous_y_shift: int,
    previous_x_shift: int,
) -> dict[str, float | int | list[StarMeasurement]]:
    guesses = [(broad_y_shift, broad_x_shift)]
    if (previous_y_shift, previous_x_shift) != (broad_y_shift, broad_x_shift):
        guesses.append((previous_y_shift, previous_x_shift))
    choices = [
        measure_shift_from_stars(
            image,
            reference_stars,
            approximate_y_shift=y_guess,
            approximate_x_shift=x_guess,
        )
        for y_guess, x_guess in guesses
    ]
    valid = [choice for choice in choices if choice is not None]
    if not valid:
        raise RuntimeError(
            "Reference stars could not confirm either the broad field match "
            "or the previous frame position"
        )
    return min(
        valid,
        key=lambda choice: (
            -int(choice["stars_used"]),
            float(choice["position_rms"]),
        ),
    )


def shift_whole_pixels(
    image: np.ndarray,
    *,
    y_shift: int,
    x_shift: int,
    fill_value: float,
) -> np.ndarray:
    output = np.full(image.shape, fill_value, dtype=np.float32)
    height, width = image.shape

    if y_shift >= 0:
        source_y = slice(0, height - y_shift)
        output_y = slice(y_shift, height)
    else:
        source_y = slice(-y_shift, height)
        output_y = slice(0, height + y_shift)
    if x_shift >= 0:
        source_x = slice(0, width - x_shift)
        output_x = slice(x_shift, width)
    else:
        source_x = slice(-x_shift, width)
        output_x = slice(0, width + x_shift)

    output[output_y, output_x] = image[source_y, source_x]
    return output


def aligned_name(path: Path) -> str:
    return path.name.replace("_out.fits", "_aligned.fits")


def add_alignment_header(
    header: fits.Header,
    *,
    reference_name: str,
    x_shift: int,
    y_shift: int,
    stars_used: int,
    position_rms: float,
    fill_value: float,
) -> fits.Header:
    updated = header.copy()
    for keyword in ("CHECKSUM", "DATASUM"):
        if keyword in updated:
            del updated[keyword]
    updated["ALIGNED"] = (True, "Image aligned to the first usable frame")
    updated["X_SHIFT"] = (x_shift, "Whole-pixel x shift applied")
    updated["Y_SHIFT"] = (y_shift, "Whole-pixel y shift applied")
    updated["ALNSTARS"] = (stars_used, "Reference stars used for alignment")
    updated["ALNRMS"] = (position_rms, "Star position RMS after shift, pixels")
    updated["ALNREF"] = (reference_name, "Alignment reference image")
    updated["ALNMETH"] = ("WHOLEPIX", "Whole-pixel reference-star alignment")
    updated["ALNFILL"] = (fill_value, "Background value used at shifted edges")
    updated["ALNVER"] = (ALIGNMENT_VERSION, "Alignment script version")
    updated["ALNDATE"] = (
        datetime.now(UTC).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "UTC alignment time",
    )
    updated.add_history(
        f"Aligned to {reference_name} using {stars_used} reference stars."
    )
    updated.add_history(
        f"Applied whole-pixel shift X={x_shift}, Y={y_shift}; no interpolation."
    )
    return updated


def write_fits_atomic(
    path: Path,
    image: np.ndarray,
    header: fits.Header,
    *,
    overwrite: bool,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists() and not overwrite:
        raise FileExistsError(path)
    temporary = path.with_name(f".{path.name}.partial")
    if temporary.exists():
        temporary.unlink()
    fits.PrimaryHDU(data=image, header=header).writeto(
        temporary,
        overwrite=True,
        checksum=True,
        output_verify="fix+warn",
    )
    os.replace(temporary, path)


def display_normalization(image: np.ndarray) -> ImageNormalize:
    finite = image[np.isfinite(image)]
    vmin, vmax = ZScaleInterval(contrast=0.20).get_limits(finite)
    return ImageNormalize(
        vmin=vmin,
        vmax=vmax,
        stretch=AsinhStretch(0.08),
        clip=True,
    )


def plot_movement(table: pd.DataFrame, output: Path) -> None:
    time_axis = table["bjd_tdb"] - 2459782.0
    figure, axes = plt.subplots(3, 1, figsize=(11, 8), sharex=True)
    axes[0].plot(time_axis, table["x_shift_pixels"], color="#315f9b")
    axes[0].set_ylabel("Horizontal shift (pixels)")
    axes[1].plot(time_axis, table["y_shift_pixels"], color="#9b4a31")
    axes[1].set_ylabel("Vertical shift (pixels)")
    axes[2].plot(time_axis, table["position_rms_pixels"], color="#3b7d4b")
    axes[2].axhline(1.0, color="0.4", linestyle=":", linewidth=1)
    axes[2].set_ylabel("Star position error (pixels)")
    axes[2].set_xlabel("Time (BJD TDB - 2459782)")
    for axis in axes:
        axis.grid(alpha=0.2)
    figure.suptitle("TOI-3505.01 image movement", fontsize=16)
    figure.tight_layout()
    figure.savefig(output, dpi=180)
    plt.close(figure)


def plot_alignment_check(
    snapshots: dict[int, tuple[np.ndarray, np.ndarray]],
    reference_stars: list[StarMeasurement],
    output: Path,
) -> None:
    sequences = sorted(snapshots)
    center = min(
        reference_stars,
        key=lambda star: (star.x - 2048.0) ** 2 + (star.y - 2048.0) ** 2,
    )
    crop_radius = 110
    x = int(round(center.x))
    y = int(round(center.y))
    figure, axes = plt.subplots(len(sequences), 2, figsize=(9, 4 * len(sequences)))
    axes = np.atleast_2d(axes)
    for row, sequence in enumerate(sequences):
        before, after = snapshots[sequence]
        for column, (image, label) in enumerate(
            ((before, "Before alignment"), (after, "After alignment"))
        ):
            crop = image[
                y - crop_radius : y + crop_radius + 1,
                x - crop_radius : x + crop_radius + 1,
            ]
            axis = axes[row, column]
            axis.imshow(
                crop,
                origin="lower",
                cmap="gray",
                norm=display_normalization(crop),
            )
            axis.axhline(crop_radius, color="#e6b800", linewidth=0.7, alpha=0.8)
            axis.axvline(crop_radius, color="#e6b800", linewidth=0.7, alpha=0.8)
            axis.set_title(f"Frame {sequence:04d}: {label}")
            axis.set_xlabel("Pixels across this close-up")
            axis.set_ylabel("Pixels up this close-up")
    figure.suptitle("TOI-3505.01 alignment check", fontsize=16)
    figure.tight_layout()
    figure.savefig(output, dpi=180)
    plt.close(figure)


def write_readme(
    diagnostics_dir: Path,
    *,
    frame_count: int,
    reference_count: int,
    median_rms: float,
    maximum_rms: float,
    wrote_images: bool,
) -> None:
    image_status = (
        f"{frame_count} aligned FITS images were written to "
        "`data/ground/toi3505/aligned`."
        if wrote_images
        else "This run measured the shifts without writing aligned FITS images."
    )
    content = f"""# TOI-3505.01 image alignment

The first {frame_count} calibrated images were matched to frame 0001. Frames
0282 and 0283 were left out because their star field was not usable in the
earlier image review.

The alignment used {reference_count} reference stars. Every image was moved by
a whole number of pixels, so the original pixel values were not interpolated.
{image_status}

- Typical star position error after alignment: {median_rms:.3f} pixels
- Largest star position error after alignment: {maximum_rms:.3f} pixels
- `frame_shifts.csv`: shift and position check for every image
- `reference_stars.csv`: stars used to measure the shifts
- `01_image_movement.png`: movement during the night
- `02_alignment_check.png`: first, middle, and last image around one star
- `summary.json`: main results and checks
- `verification.json`: file checks, observation times, and exact pixel checks

The next course step is to locate the target, save its Seeing Profile, and use
those radii for multi-aperture photometry.
"""
    (diagnostics_dir / "README.md").write_text(content, encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, default=DEFAULT_INPUT_DIR)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument(
        "--diagnostics-dir", type=Path, default=DEFAULT_DIAGNOSTICS_DIR
    )
    parser.add_argument("--first-sequence", type=int, default=FIRST_SEQUENCE)
    parser.add_argument("--last-sequence", type=int, default=LAST_USABLE_SEQUENCE)
    parser.add_argument("--maximum-stars", type=int, default=25)
    parser.add_argument(
        "--measure-only",
        action="store_true",
        help="Measure and check the shifts without writing aligned FITS images",
    )
    parser.add_argument(
        "--overwrite", action="store_true", help="Replace existing aligned images"
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    input_dir = args.input_dir.resolve()
    output_dir = args.output_dir.resolve()
    diagnostics_dir = args.diagnostics_dir.resolve()
    diagnostics_dir.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()

    frames = discover_frames(
        input_dir,
        first_sequence=args.first_sequence,
        last_sequence=args.last_sequence,
    )
    reference_image, reference_header = read_frame(frames[0].path)
    reference_stars = find_reference_stars(
        reference_image,
        maximum_stars=args.maximum_stars,
    )
    reference_signature = field_signature(reference_image)
    reference_fft = rfftn(reference_signature)

    pd.DataFrame(
        [
            {
                "star": index,
                "x_pixel": star.x,
                "y_pixel": star.y,
                "brightness_adu": star.flux,
                "width_pixels": star.fwhm,
            }
            for index, star in enumerate(reference_stars, start=1)
        ]
    ).to_csv(diagnostics_dir / "reference_stars.csv", index=False)

    snapshot_sequences = {
        frames[0].sequence,
        frames[len(frames) // 2].sequence,
        frames[-1].sequence,
    }
    snapshots: dict[int, tuple[np.ndarray, np.ndarray]] = {}
    rows: list[dict[str, object]] = []
    previous_y_shift = 0
    previous_x_shift = 0

    for index, frame in enumerate(frames):
        frame_started = time.perf_counter()
        if index == 0:
            image = reference_image
            header = reference_header
            broad_y_shift = 0
            broad_x_shift = 0
            correlation_peak = 1.0
        else:
            image, header = read_frame(frame.path)
            broad_y_shift, broad_x_shift, correlation_peak = broad_field_shift(
                reference_fft,
                image,
                signature_shape=reference_signature.shape,
            )

        result = choose_shift(
            image,
            reference_stars,
            broad_y_shift=broad_y_shift,
            broad_x_shift=broad_x_shift,
            previous_y_shift=previous_y_shift,
            previous_x_shift=previous_x_shift,
        )
        x_shift = int(result["x_shift"])
        y_shift = int(result["y_shift"])
        previous_x_shift = x_shift
        previous_y_shift = y_shift
        fill_value = float(np.median(image[::8, ::8]))
        output_path = output_dir / aligned_name(frame.path)

        aligned: np.ndarray | None = None
        status = "measured"
        if not args.measure_only or frame.sequence in snapshot_sequences:
            aligned = shift_whole_pixels(
                image,
                y_shift=y_shift,
                x_shift=x_shift,
                fill_value=fill_value,
            )
        if not args.measure_only:
            aligned_header = add_alignment_header(
                header,
                reference_name=frames[0].path.name,
                x_shift=x_shift,
                y_shift=y_shift,
                stars_used=int(result["stars_used"]),
                position_rms=float(result["position_rms"]),
                fill_value=fill_value,
            )
            assert aligned is not None
            write_fits_atomic(
                output_path,
                aligned,
                aligned_header,
                overwrite=args.overwrite,
            )
            status = "written"

        if frame.sequence in snapshot_sequences:
            assert aligned is not None
            snapshots[frame.sequence] = (image.copy(), aligned.copy())

        row = {
            "sequence": frame.sequence,
            "input_file": frame.path.name,
            "output_file": output_path.name,
            "status": status,
            "bjd_tdb": float(header.get("BJD_TDB", math.nan)),
            "broad_x_shift_pixels": broad_x_shift,
            "broad_y_shift_pixels": broad_y_shift,
            "x_shift_pixels": x_shift,
            "y_shift_pixels": y_shift,
            "measured_x_shift_pixels": float(result["measured_x_shift"]),
            "measured_y_shift_pixels": float(result["measured_y_shift"]),
            "reference_stars_used": int(result["stars_used"]),
            "position_rms_pixels": float(result["position_rms"]),
            "median_star_width_pixels": float(result["median_star_width"]),
            "field_match_strength": correlation_peak,
            "fill_value_adu": fill_value,
            "elapsed_seconds": float(time.perf_counter() - frame_started),
        }
        rows.append(row)
        if (index + 1) % 10 == 0 or index == 0 or index + 1 == len(frames):
            print(
                f"Aligned {index + 1}/{len(frames)} images "
                f"(frame {frame.sequence:04d}, x={x_shift}, y={y_shift}, "
                f"stars={result['stars_used']}, error={result['position_rms']:.3f})",
                flush=True,
            )
        if index != 0:
            del image
        if aligned is not None and frame.sequence not in snapshot_sequences:
            del aligned

    table = pd.DataFrame(rows)
    table.to_csv(diagnostics_dir / "frame_shifts.csv", index=False)
    plot_movement(table, diagnostics_dir / "01_image_movement.png")
    plot_alignment_check(
        snapshots,
        reference_stars,
        diagnostics_dir / "02_alignment_check.png",
    )

    summary = {
        "status": "pass",
        "alignment_version": ALIGNMENT_VERSION,
        "created_utc": datetime.now(UTC).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "reference_frame": frames[0].path.name,
        "first_sequence": frames[0].sequence,
        "last_sequence": frames[-1].sequence,
        "aligned_frames": len(frames),
        "excluded_sequences": [282, 283]
        if args.first_sequence == 1 and args.last_sequence == 281
        else [],
        "reference_stars": len(reference_stars),
        "minimum_stars_used": int(table["reference_stars_used"].min()),
        "median_stars_used": float(table["reference_stars_used"].median()),
        "whole_pixel_shifts": True,
        "interpolation_applied": False,
        "median_position_rms_pixels": float(table["position_rms_pixels"].median()),
        "maximum_position_rms_pixels": float(table["position_rms_pixels"].max()),
        "x_shift_range_pixels": [
            int(table["x_shift_pixels"].min()),
            int(table["x_shift_pixels"].max()),
        ],
        "y_shift_range_pixels": [
            int(table["y_shift_pixels"].min()),
            int(table["y_shift_pixels"].max()),
        ],
        "bjd_tdb_start": float(table["bjd_tdb"].min()),
        "bjd_tdb_end": float(table["bjd_tdb"].max()),
        "aligned_images_written": not args.measure_only,
        "elapsed_seconds": float(time.perf_counter() - started),
    }
    (diagnostics_dir / "summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    write_readme(
        diagnostics_dir,
        frame_count=len(frames),
        reference_count=len(reference_stars),
        median_rms=summary["median_position_rms_pixels"],
        maximum_rms=summary["maximum_position_rms_pixels"],
        wrote_images=not args.measure_only,
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
