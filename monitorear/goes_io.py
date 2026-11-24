"""List and download one ABI product prefix from the public bucket.

s3fs caches directory listings on the filesystem object. Every list
drops that cache first, so a later pass sees files that arrived after
the previous one. ``fs.ls`` has no refresh argument.
"""

from __future__ import annotations

from datetime import datetime, timedelta, timezone
from pathlib import Path

import s3fs

BUCKET = "noaa-goes19"
ROOT = Path(__file__).resolve().parents[1]


def filesystem():
    """Anonymous client for the NOAA bucket in us-east-1."""
    return s3fs.S3FileSystem(anon=True, client_kwargs={"region_name": "us-east-1"})


def parse_goes_time(token):
    """Parse the 14-digit s/e/c timestamp in an ABI filename.

    Layout is YYYYDDDHHMMSSt, with the last digit a tenth of a second.
    """
    t = token[1:] if token[:1] in "sec" else token
    t = t.replace(".nc", "")
    if len(t) != 14:
        raise ValueError(f"unexpected GOES timestamp {token!r}")
    year, doy = int(t[0:4]), int(t[4:7])
    hour, minute, second, tenth = int(t[7:9]), int(t[9:11]), int(t[11:13]), int(t[13])
    base = datetime(year, 1, 1, tzinfo=timezone.utc) + timedelta(days=doy - 1)
    return base + timedelta(
        hours=hour, minutes=minute, seconds=second, milliseconds=tenth * 100
    )


def parse_abi_name(name):
    """Split an ABI L1b object name into sector, mode, and the three times."""
    parts = name.replace(".nc", "").split("_")
    # OR, ABI-L1b-RadM1-M6C07, G19, s..., e..., c...
    product = parts[1]
    sector = "M1" if "RadM1" in product else "M2" if "RadM2" in product else "?"
    mode = "M6" if "-M6" in product else product
    return {
        "sector": sector,
        "mode": mode,
        "product": product,
        "satellite": parts[2],
        "start": parse_goes_time(parts[3]),
        "end": parse_goes_time(parts[4]),
        "created": parse_goes_time(parts[5]),
    }


def sector_label(name):
    if "RadM1" in name:
        return "M1"
    if "RadM2" in name:
        return "M2"
    if "RadC" in name:
        return "CONUS"
    if "RadF" in name:
        return "FD"
    return "?"


def list_c07(fs, prefix, t0, t1):
    fs.invalidate_cache(prefix)
    chosen = []
    for key in fs.ls(prefix):
        name = key.split("/")[-1]
        if "C07_G19_" not in name:
            continue
        meta = parse_abi_name(name)
        meta["sector"] = sector_label(name)
        meta["name"] = name
        if t0 <= meta["start"] <= t1:
            chosen.append((key, meta))
    chosen.sort(key=lambda item: item[1]["start"])
    return chosen


def download(fs, key, dest_dir):
    dest = dest_dir / key.split("/")[-1]
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    fs.get(key, str(dest))
    return dest
