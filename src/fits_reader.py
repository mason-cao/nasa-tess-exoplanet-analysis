"""
Minimal FITS reader for this TESS light-curve project.

Why this exists:
The internship lecture emphasized that science data are ultimately bytes with a
header/metadata layer. Instead of hiding that with a large astronomy library,
this file reads the TESS FITS binary table using only Python + NumPy. That makes
this assignment reproducible in lightweight environments while also showing the
connection between the lecture's bits/bytes/FITS discussion and real data work.

Supported here: the FITS structures used by the included TESS light curve and
Data Validation Time Series FITS files. For general research work, astropy.io.fits
is still the standard tool.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any, Dict, List, Tuple

import numpy as np


def _read_header_cards(file_obj) -> List[str]:
    """Read FITS header blocks until the END card."""
    header = b""
    while True:
        block = file_obj.read(2880)
        if not block:
            raise EOFError("Reached end of file while reading a FITS header.")
        header += block
        for i in range(0, len(header), 80):
            card = header[i : i + 80]
            if card[:3] == b"END" and card[3:8].strip() == b"":
                return [header[j : j + 80].decode("ascii", errors="replace") for j in range(0, len(header), 80)]


def _split_value_comment(raw: str) -> str:
    """Return the FITS value field without the comment, respecting quoted strings."""
    in_quote = False
    chars = []
    for ch in raw.rstrip():
        if ch == "'":
            in_quote = not in_quote
        if ch == "/" and not in_quote:
            break
        chars.append(ch)
    return "".join(chars).strip()


def _parse_value(raw: str) -> Any:
    value = _split_value_comment(raw)
    if value.startswith("'"):
        out = []
        i = 1
        while i < len(value):
            if value[i] == "'":
                if i + 1 < len(value) and value[i + 1] == "'":
                    out.append("'")
                    i += 2
                    continue
                break
            out.append(value[i])
            i += 1
        return "".join(out)
    if value == "T":
        return True
    if value == "F":
        return False
    if value == "":
        return None
    try:
        if any(ch in value for ch in ".EeDd"):
            return float(value.replace("D", "E"))
        return int(value)
    except ValueError:
        return value


def _parse_header(cards: List[str]) -> Dict[str, Any]:
    header = {}
    for card in cards:
        key = card[:8].strip()
        if key == "END":
            break
        if card[8:10] == "= ":
            header[key] = _parse_value(card[10:80])
    return header


def _data_size(header: Dict[str, Any]) -> int:
    naxis = int(header.get("NAXIS") or 0)
    if str(header.get("XTENSION", "")).strip() == "BINTABLE":
        return int(header["NAXIS1"]) * int(header["NAXIS2"]) + int(header.get("PCOUNT") or 0)
    if naxis:
        size = abs(int(header.get("BITPIX", 8))) // 8
        for axis in range(1, naxis + 1):
            size *= int(header[f"NAXIS{axis}"])
        size *= int(header.get("GCOUNT") or 1)
        return size
    return 0


def list_hdus(path: str | Path) -> List[Dict[str, Any]]:
    """Return HDU metadata, including data offsets and parsed headers."""
    path = Path(path)
    hdus = []
    with path.open("rb") as f:
        while True:
            try:
                offset = f.tell()
                cards = _read_header_cards(f)
            except EOFError:
                break
            header = _parse_header(cards)
            size = _data_size(header)
            data_offset = f.tell()
            hdus.append({"offset": offset, "data_offset": data_offset, "size": size, "header": header, "cards": cards})
            padding = (2880 - size % 2880) % 2880
            f.seek(size + padding, 1)
    return hdus


def _parse_tform(tform: str) -> Tuple[int, str, int]:
    match = re.fullmatch(r"(\d*)([A-Z])(?:\(.*\))?", str(tform).strip())
    if not match:
        raise ValueError(f"Unsupported FITS TFORM: {tform}")
    repeat = int(match.group(1) or 1)
    code = match.group(2)
    byte_sizes = {"L": 1, "B": 1, "I": 2, "J": 4, "K": 8, "A": 1, "E": 4, "D": 8}
    if code not in byte_sizes:
        raise ValueError(f"Unsupported FITS TFORM code: {tform}")
    return repeat, code, repeat * byte_sizes[code]


def _numpy_format(code: str, repeat: int):
    formats = {"L": "?", "B": "u1", "I": ">i2", "J": ">i4", "K": ">i8", "A": "S1", "E": ">f4", "D": ">f8"}
    base = formats[code]
    return base if repeat == 1 else (base, (repeat,))


def read_bintable(path: str | Path, hdu_index: int = 1):
    """Read a FITS binary table HDU into a NumPy structured array."""
    path = Path(path)
    hdus = list_hdus(path)
    hdu = hdus[hdu_index]
    header = hdu["header"]
    nrows = int(header["NAXIS2"])
    row_size = int(header["NAXIS1"])
    nfields = int(header["TFIELDS"])

    names, formats, offsets = [], [], []
    offset = 0
    for idx in range(1, nfields + 1):
        name = str(header[f"TTYPE{idx}"]).strip()
        repeat, code, nbytes = _parse_tform(header[f"TFORM{idx}"])
        names.append(name)
        formats.append(_numpy_format(code, repeat))
        offsets.append(offset)
        offset += nbytes

    dtype = np.dtype({"names": names, "formats": formats, "offsets": offsets, "itemsize": row_size})
    with path.open("rb") as f:
        f.seek(hdu["data_offset"])
        table = np.fromfile(f, dtype=dtype, count=nrows)
    return table, header


def read_image(path: str | Path, hdu_index: int = 2):
    """Read a simple FITS image HDU into a NumPy array."""
    path = Path(path)
    hdus = list_hdus(path)
    hdu = hdus[hdu_index]
    header = hdu["header"]
    bitpix = int(header["BITPIX"])
    dims = [int(header[f"NAXIS{axis}"]) for axis in range(1, int(header["NAXIS"]) + 1)]
    dtype = {8: "u1", 16: ">i2", 32: ">i4", -32: ">f4", -64: ">f8"}[bitpix]
    with path.open("rb") as f:
        f.seek(hdu["data_offset"])
        data = np.fromfile(f, dtype=dtype, count=int(np.prod(dims)))
    return data.reshape(tuple(reversed(dims))), header
