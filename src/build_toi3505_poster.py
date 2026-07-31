"""Assemble the symposium poster for TOI-3505.01.

Reads the frozen analysis products and the poster figures, optimises each
image, and writes a self-contained print-ready HTML poster sized for a
48 x 36 inch board.
"""

from __future__ import annotations

import base64
import io
import json
import sys
from pathlib import Path

from PIL import Image

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs" / "toi3505_poster"
OUT.mkdir(parents=True, exist_ok=True)

# The program supplies the George Mason University mark on its poster
# template. It is lifted from the template rather than redrawn so the poster
# matches the other Schar Scholars boards.
LOGO_SOURCE = Path.home() / "Downloads" / "CarrieS_TOI5443.png"
LOGO_BOX = (68, 42, 428, 462)


def encode(image: Image.Image, *, jpeg: bool, quality: int = 88) -> str:
    buffer = io.BytesIO()
    if jpeg:
        image.convert("RGB").save(buffer, "JPEG", quality=quality, optimize=True, progressive=True)
        mime = "image/jpeg"
    else:
        image.save(buffer, "PNG", optimize=True)
        mime = "image/png"
    payload = base64.b64encode(buffer.getvalue()).decode("ascii")
    return f"data:{mime};base64,{payload}"


# The Schar paper-writing lecture asks for lossless figures: "Figures should be
# lossless file types - NO JPEGs. PNGs preferred." Plot panels are flat colour
# on white, so PNG is also smaller than JPEG for most of them.
def prepare(path: Path, width: int, *, jpeg: bool = False, quality: int = 92) -> str:
    image = Image.open(path)
    if image.mode in ("RGBA", "LA", "P"):
        background = Image.new("RGB", image.size, "white")
        converted = image.convert("RGBA")
        background.paste(converted, mask=converted.split()[-1])
        image = background
    if image.width > width:
        height = round(image.height * width / image.width)
        image = image.resize((width, height), Image.LANCZOS)
    return encode(image, jpeg=jpeg, quality=quality)


def prepare_logo() -> str:
    image = Image.open(LOGO_SOURCE).crop(LOGO_BOX)
    image = image.resize((360, 420), Image.LANCZOS)
    return encode(image, jpeg=False)


# Widths are the native figure resolutions. At the printed size each panel
# lands near 100-120 dots per inch on a 48 x 36 inch board.
VARIANTS = {
    "v1": {
        "template": "poster_template.html",
        "stem": "TOI-3505.01_Mason_Cao_poster",
        "figures": {
            "field": (OUT / "04_field_and_aperture.png", 2280),
            "ground": (OUT / "05_ground_light_curve.png", 2288),
            "timing": (OUT / "01_transit_timing.png", 2772),
            "phase": (OUT / "06_phase_folded.png", 2640),
            "neb": (OUT / "02_nearby_star_screen.png", 2023),
        },
    },
    # Second layout: the AstroImageJ products the rest of the cohort shows.
    "v2": {
        "template": "poster_template_v2.html",
        "stem": "TOI-3505.01_Mason_Cao_poster_v2",
        "figures": {
            "seeing": (OUT / "v2_01_seeing_profile.png", 1672),
            "compfield": (OUT / "v2_02_comparison_field.png", 1297),
            "neb": (OUT / "v2_04_neb_field.png", 1629),
            "lightcurve": (ROOT / "outputs" / "toi3505_discord_post" / "01_TOI_3505.01_final_light_curve.png", 2023),
            "dmag": (OUT / "v2_03_dmag_rms.png", 1716),
        },
    },
}


def main() -> None:
    variant = sys.argv[1] if len(sys.argv) > 1 else "v1"
    if variant not in VARIANTS:
        raise SystemExit(f"unknown variant {variant!r}; choose from {sorted(VARIANTS)}")
    config = VARIANTS[variant]

    analysis = json.loads((OUT / "poster_analysis.json").read_text())

    images = {name: prepare(path, width) for name, (path, width) in config["figures"].items()}
    images["logo"] = prepare_logo()

    template = (ROOT / "src" / config["template"]).read_text()
    for name, uri in images.items():
        template = template.replace(f"{{{{{name}}}}}", uri)

    destination = OUT / f"{config['stem']}.html"
    destination.write_text(template)

    size_mb = destination.stat().st_size / 1e6
    print(f"Wrote {destination} ({size_mb:.1f} MB)")
    for name, uri in images.items():
        print(f"  {name:<11} {len(uri) / 1e6:.2f} MB")

    remaining = [token for token in ("{{", "}}") if token in template]
    if remaining:
        print("WARNING: unreplaced template tokens remain")

    ephemeris = analysis["ephemeris"]["adopted_error_scaled"]
    print(
        f"\nHeadline: P = {ephemeris['period_days']:.7f} +/- "
        f"{ephemeris['period_error_days']:.2e} d over "
        f"{analysis['ephemeris']['baseline_years']:.2f} yr"
    )


if __name__ == "__main__":
    main()
