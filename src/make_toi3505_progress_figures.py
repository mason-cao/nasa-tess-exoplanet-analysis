"""Make simple, reproducible progress figures for the TOI-3505.01 data set."""

from __future__ import annotations

import json
import re
import zipfile
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from fits_reader import _parse_header, _read_header_cards


ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / "data_and_lectures"
OUTPUT_DIR = ROOT / "outputs" / "toi3505_progress"

SCIENCE_PATTERN = re.compile(
    r"(?:^|/)TOI_3505\.01_50\.000s_R-(\d{4})(?:\(1\))?\.fits$"
)

# Current TESS Project Candidate values from the NASA Exoplanet Archive.
PERIOD_DAYS = 2.9151556
TRANSIT_EPOCH_BJD = 2459793.5343850


def read_primary_header(archive: zipfile.ZipFile, member: str) -> dict:
    with archive.open(member) as file_obj:
        return _parse_header(_read_header_cards(file_obj))


def science_member_priority(member: str) -> tuple[int, str]:
    """Prefer the normally named frame over a '(1)' duplicate."""
    return ("(1)" in Path(member).name, member)


def collect_headers() -> tuple[pd.DataFrame, dict]:
    archives = sorted(DATA_DIR.glob("TOI_3505.01-20260714T190730Z-1-00*.zip"))
    if not archives:
        raise FileNotFoundError("No TOI-3505.01 archives were found")

    science_candidates: dict[int, list[tuple[Path, str]]] = {}
    dark_50: set[str] = set()
    dark_3_5: set[str] = set()
    flats: set[str] = set()
    focus_fits: set[str] = set()
    focus_png: set[str] = set()

    for archive_path in archives:
        with zipfile.ZipFile(archive_path) as archive:
            for member in archive.namelist():
                name = Path(member).name
                science_match = SCIENCE_PATTERN.search(member)
                if science_match:
                    frame = int(science_match.group(1))
                    science_candidates.setdefault(frame, []).append((archive_path, member))
                elif re.fullmatch(r"Dark_50\.000s-\d{4}\.fits", name):
                    dark_50.add(name)
                elif re.fullmatch(r"Dark_3\.500s-\d{4}\.fits", name):
                    dark_3_5.add(name)
                elif re.fullmatch(r"Flat_3\.500s_R-\d{4}-final\.fits", name):
                    flats.add(name)
                elif re.fullmatch(r"FocuserImage_.*\.fits", name):
                    focus_fits.add(name)
                elif member.lower().endswith(".png") and "/focuser_images/" in member:
                    focus_png.add(name)

    rows = []
    for frame, candidates in sorted(science_candidates.items()):
        archive_path, member = min(
            candidates,
            key=lambda item: science_member_priority(item[1]),
        )
        with zipfile.ZipFile(archive_path) as archive:
            header = read_primary_header(archive, member)
        rows.append(
            {
                "frame": frame,
                "date_obs": header["DATE-OBS"],
                "bjd_tdb": float(header["BJD_TDB"]),
                "airmass": float(header["AIRMASS"]),
                "ccd_temp_c": float(header["CCD-TEMP"]),
                "exposure_s": float(header["EXPTIME"]),
                "filter": str(header["FILTER"]).strip(),
                "object": str(header["OBJECT"]).strip(),
                "archive": archive_path.name,
                "member": member,
            }
        )

    table = pd.DataFrame(rows).sort_values("frame").reset_index(drop=True)
    inventory = {
        "archive_count": len(archives),
        "science_frames": len(table),
        "duplicate_science_files": sum(
            max(0, len(candidates) - 1) for candidates in science_candidates.values()
        ),
        "dark_50s": len(dark_50),
        "dark_3_5s": len(dark_3_5),
        "flat_3_5s_r": len(flats),
        "focus_fits": len(focus_fits),
        "focus_png": len(focus_png),
        "first_bjd_tdb": float(table["bjd_tdb"].min()),
        "last_bjd_tdb": float(table["bjd_tdb"].max()),
        "span_hours": float(
            (table["bjd_tdb"].max() - table["bjd_tdb"].min()) * 24.0
        ),
    }
    return table, inventory


def clean_axes(ax: plt.Axes) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y", color="#d8d8d8", linewidth=0.7, alpha=0.7)


def plot_airmass(table: pd.DataFrame) -> None:
    x = table["bjd_tdb"].to_numpy() - 2459782.0
    y = table["airmass"].to_numpy()

    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.plot(x, y, color="#3268b2", linewidth=1.2)
    ax.scatter(x, y, color="#3268b2", s=7)
    ax.set_title("TOI-3505.01")
    ax.set_xlabel("BJD_TDB - 2459782")
    ax.set_ylabel("Airmass")
    clean_axes(ax)
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "02_airmass.png", dpi=180)
    plt.close(fig)


def plot_timing(table: pd.DataFrame) -> dict:
    observed_start = float(table["bjd_tdb"].min())
    observed_end = float(table["bjd_tdb"].max())
    observed_midpoint = 0.5 * (observed_start + observed_end)

    cycle = int(round((observed_midpoint - TRANSIT_EPOCH_BJD) / PERIOD_DAYS))
    earlier_transit = TRANSIT_EPOCH_BJD + cycle * PERIOD_DAYS
    if earlier_transit > observed_start:
        earlier_transit -= PERIOD_DAYS
    later_transit = earlier_transit + PERIOD_DAYS

    offset = 2459780.0
    fig, ax = plt.subplots(figsize=(8, 3.2))
    ax.hlines(
        0,
        observed_start - offset,
        observed_end - offset,
        color="#3268b2",
        linewidth=8,
    )
    ax.text(
        0.5 * (observed_start + observed_end) - offset,
        0.12,
        "images",
        ha="center",
        va="bottom",
    )
    for transit in (earlier_transit, later_transit):
        ax.vlines(
            transit - offset,
            -0.12,
            0.12,
            color="#c94444",
            linewidth=4,
        )
        ax.text(
            transit - offset,
            0.18,
            "transit",
            ha="center",
            va="bottom",
        )

    ax.set_title("TOI-3505.01 timing")
    ax.set_xlabel("BJD_TDB - 2459780")
    ax.set_yticks([])
    ax.set_ylim(-0.32, 0.42)
    ax.set_xlim(
        earlier_transit - offset - 0.25,
        later_transit - offset + 0.25,
    )
    ax.spines["left"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "03_timing.png", dpi=180)
    plt.close(fig)

    return {
        "period_days": PERIOD_DAYS,
        "transit_epoch_bjd": TRANSIT_EPOCH_BJD,
        "previous_transit_bjd": earlier_transit,
        "next_transit_bjd": later_transit,
    }


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 11,
            "axes.titlesize": 13,
            "axes.labelsize": 11,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
        }
    )

    table, inventory = collect_headers()
    table.to_csv(OUTPUT_DIR / "toi3505_frame_headers.csv", index=False)
    plot_airmass(table)
    timing = plot_timing(table)

    summary = {**inventory, **timing}
    (OUTPUT_DIR / "summary.json").write_text(
        json.dumps(summary, indent=2) + "\n",
        encoding="utf-8",
    )

    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
