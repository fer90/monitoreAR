#!/usr/bin/env python3
"""Milestone 5: serve the stored candidates as JSON and on one map.

The page reads SQLite. It does not recompute a score and it does not
download NetCDF. Each candidate's thermal crop is a PNG written when the
scan is scored, or filled in here if that file is still missing.

Run from the repository root, after the poller has created the database:

    .venv/bin/python scripts/m5_api.py

The process listens on http://127.0.0.1:8000/ until interrupted.
"""

from __future__ import annotations

import json
import math
import re
import sqlite3
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import s3fs
import uvicorn
from fastapi import FastAPI, HTTPException
from fastapi.responses import FileResponse, HTMLResponse

sys.path.insert(0, str(Path(__file__).resolve().parent))
import m1_starship_flight11 as m1
import m4_poll as m4

ROOT = m1.ROOT
CROP_DIR = ROOT / "outputs" / "crops"
REPORT_PATH = ROOT / "outputs" / "m5_api.txt"
PAGE_PATH = Path(__file__).resolve().parent / "m5_map.html"
HOST = "127.0.0.1"
PORT = 8000

# The Flight 11 figure scale. Both candidate crops use it, so a 340 K
# pixel is the same colour on each scan. See outputs/m1_inspection.txt.
VMIN_K = 283.43
VMAX_K = 350.0
CROP_NAME = re.compile(r"OR_ABI-[A-Za-z0-9._-]+\.png")

# Columns returned by GET /detections. crop_url is added beside them.
PUBLIC_COLUMNS = (
    "source_key",
    "watch",
    "scan_start",
    "scan_end",
    "sector",
    "pixel_lat",
    "pixel_lon",
    "distance_km",
    "score_k",
    "baseline_k",
    "baseline_n",
    "baseline_samples_k",
    "near_bt_k",
    "near_dqf",
    "max_bt_k",
    "p98_bt_k",
    "ring_bt_k",
    "solar_zenith_deg",
    "dqf_counts",
    "hot_pixel_dqf",
    "n_conditional_dqf",
    "n_neighbour_pixels",
    "saturated",
    "bt_ceiling_k",
    "consecutive",
    "status",
    "hold_reasons",
    "fire_file",
    "fire_nearest_mask",
    "fire_nearest_meaning",
    "fire_hits_within_20km",
)
JSON_COLUMNS = {
    "baseline_samples_k",
    "dqf_counts",
    "hold_reasons",
    "fire_hits_within_20km",
}

# The two Flight 11 emits. A later candidate may be stored beside them.
FLIGHT11_STARTS = (
    "2025-10-13 23:23:54.9",
    "2025-10-13 23:24:54.9",
)

app = FastAPI(title="monitoreAR", docs_url=None, redoc_url=None, openapi_url=None)


def crop_png_name(source_key):
    name = source_key.split("/")[-1]
    if not name.endswith(".nc"):
        raise ValueError(f"source key is not a NetCDF object: {source_key}")
    return name[:-3] + ".png"


def resolve_crop(name):
    """PNG path for a crop URL, or None when the name is not one of ours."""
    if not CROP_NAME.fullmatch(name):
        return None
    path = (CROP_DIR / name).resolve()
    if path.parent != CROP_DIR.resolve():
        return None
    return path


def normalize_since(value):
    if value is None or value.strip() == "":
        return None
    text = value.strip().replace("T", " ")
    if text.endswith("Z"):
        text = text[:-1]
    return text


def self_test_api():
    name = crop_png_name(
        "noaa-goes19/ABI-L1b-RadM/2025/286/23/"
        "OR_ABI-L1b-RadM2-M6C07_G19_s20252862323549_e20252862324019_c20252862324058.nc"
    )
    expected = "OR_ABI-L1b-RadM2-M6C07_G19_s20252862323549_e20252862324019_c20252862324058.png"
    if name != expected or resolve_crop(name) is None:
        raise SystemExit(f"crop name self-test failed: {name}")
    if resolve_crop("../" + name) is not None or resolve_crop(name + ".txt") is not None:
        raise SystemExit("crop name self-test accepted a bad name")
    if normalize_since("2025-10-13T23:24:00Z") != "2025-10-13 23:24:00":
        raise SystemExit("since self-test failed")
    if normalize_since("  ") is not None:
        raise SystemExit("blank since self-test failed")
    print("api self-test: crop name, since")


def write_crop(nc_path, source_key, score_k, status):
    """One geolocated band-7 crop. Skip the write when the PNG is already there."""
    CROP_DIR.mkdir(parents=True, exist_ok=True)
    dest = CROP_DIR / crop_png_name(source_key)
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    meta = m1.parse_abi_name(Path(nc_path).name)
    frame = m1.load_frame(Path(nc_path), meta)
    crop = frame["crop"]
    if crop is None:
        raise SystemExit(f"{Path(nc_path).name} has no Starbase crop")
    fig, ax = plt.subplots(figsize=(6.6, 5.8))
    mesh = ax.pcolormesh(
        crop["lon_edges"],
        crop["lat_edges"],
        crop["bt"],
        cmap="inferno",
        vmin=VMIN_K,
        vmax=VMAX_K,
        shading="flat",
        rasterized=True,
    )
    ax.plot(
        m1.PAD_LON,
        m1.PAD_LAT,
        marker="o",
        ms=8,
        mfc="none",
        mec="cyan",
        mew=1.5,
        label="Pad 1",
        zorder=4,
    )
    # Same mark as the Milestone 1 figure: only a hot pixel that is not
    # the pad pixel itself.
    anomalous = np.isfinite(crop["max_bt_20km"]) and crop["max_bt_20km"] > crop["median_bt"] + 10
    if anomalous and crop["max_dist_km"] > 3:
        ax.plot(
            crop["max_lon"],
            crop["max_lat"],
            marker="+",
            ms=12,
            mec="white",
            mew=1.4,
            label="hottest pixel within 20 km",
            zorder=5,
        )
    start = meta["start"]
    ax.set_title(
        f"{meta['sector']}  {start:%Y-%m-%d %H:%M:%S} UTC\n"
        f"pad pixel {crop['pad_pixel_bt']:.1f} K"
        f"   max within 20 km {crop['max_bt_20km']:.1f} K"
        f" ({crop['max_dist_km']:.1f} km)\n"
        f"score {float(score_k):.2f} K  {status}",
        fontsize=9,
    )
    ax.set_xlim(m1.CROP_LON)
    ax.set_ylim(m1.CROP_LAT)
    ax.set_aspect(1 / math.cos(math.radians(m1.PAD_LAT)))
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.grid(True, color="white", alpha=0.25, linewidth=0.4)
    ax.legend(loc="lower left", fontsize=7, framealpha=0.85)
    colorbar = fig.colorbar(mesh, ax=ax, fraction=0.046, pad=0.03)
    colorbar.set_label("Brightness temperature (K)")
    fig.subplots_adjust(left=0.14, right=0.88, top=0.80, bottom=0.16)
    fig.text(
        0.14,
        0.02,
        "GOES-19 ABI band 7. Scale 283.43–350 K, the Flight 11 figure. "
        "Pixel average, not a flame temperature.",
        fontsize=7.5,
    )
    fig.savefig(dest, dpi=120)
    plt.close(fig)
    return dest


def read_conn():
    conn = sqlite3.connect(f"file:{m4.DB_PATH}?mode=ro", uri=True)
    conn.row_factory = sqlite3.Row
    return conn


def public_row(row):
    data = {}
    for key in PUBLIC_COLUMNS:
        value = row[key]
        if key in JSON_COLUMNS and isinstance(value, str):
            value = json.loads(value)
        if key in {"score_k", "max_bt_k", "near_bt_k", "baseline_k", "p98_bt_k", "ring_bt_k", "bt_ceiling_k", "distance_km"} and value is not None:
            value = round(float(value), 2)
        if key == "saturated":
            value = bool(value)
        data[key] = value
    data["crop_url"] = "/crops/" + crop_png_name(row["source_key"])
    return data


def candidate_rows(conn, since):
    if since is None:
        rows = conn.execute(
            """
            SELECT * FROM scans
            WHERE status IN ('emit', 'held_back')
            ORDER BY scan_start DESC
            """
        ).fetchall()
    else:
        rows = conn.execute(
            """
            SELECT * FROM scans
            WHERE status IN ('emit', 'held_back') AND scan_start >= ?
            ORDER BY scan_start DESC
            """,
            (since,),
        ).fetchall()
    return [public_row(row) for row in rows]


def coverage_body(conn):
    latest = conn.execute(
        """
        SELECT scan_start, scan_end, sector, covered, gap_reason, status, score_k
        FROM scans
        WHERE watch = ?
        ORDER BY scan_start DESC
        LIMIT 1
        """,
        (m4.WATCH["name"],),
    ).fetchone()
    gaps = conn.execute(
        """
        SELECT reason, sector, gap_start, gap_end, detail
        FROM gaps
        WHERE watch = ? AND resolved_by_key IS NULL
        ORDER BY gap_start
        """,
        (m4.WATCH["name"],),
    ).fetchall()
    latest_body = None
    if latest is not None:
        latest_body = {
            "scan_start": latest["scan_start"],
            "scan_end": latest["scan_end"],
            "sector": latest["sector"],
            "covered": bool(latest["covered"]),
            "gap_reason": latest["gap_reason"],
            "status": latest["status"],
            "score_k": None if latest["score_k"] is None else round(float(latest["score_k"]), 2),
        }
    return {
        "watch": m4.WATCH["name"],
        "pad": {"lat": m1.PAD_LAT, "lon": m1.PAD_LON},
        "latest_scan": latest_body,
        "open_gaps": [dict(row) for row in gaps],
    }


def ensure_crops(conn):
    """Write a PNG for every candidate that does not have one yet."""
    rows = conn.execute(
        """
        SELECT source_key, score_k, status FROM scans
        WHERE status IN ('emit', 'held_back')
        ORDER BY scan_start
        """
    ).fetchall()
    fs = None
    written = []
    for row in rows:
        dest = CROP_DIR / crop_png_name(row["source_key"])
        if not (dest.exists() and dest.stat().st_size > 0):
            if fs is None:
                fs = s3fs.S3FileSystem(anon=True, client_kwargs={"region_name": "us-east-1"})
            path = m4.locate(fs, row["source_key"])
            dest = write_crop(path, row["source_key"], row["score_k"], row["status"])
        written.append(dest)
    return written


def accept_flight11(rows):
    found = {row["scan_start"]: row for row in rows}
    missing = [start for start in FLIGHT11_STARTS if start not in found]
    if missing:
        raise SystemExit(f"Flight 11 candidates missing from the API: {missing}")
    first, second = found[FLIGHT11_STARTS[0]], found[FLIGHT11_STARTS[1]]
    if first["status"] != "emit" or second["status"] != "emit":
        raise SystemExit("Flight 11 candidates are not both emitted")
    if first["score_k"] != 39.29 or second["score_k"] != 45.37:
        raise SystemExit(f"Flight 11 scores changed: {first['score_k']}, {second['score_k']}")
    for row in (first, second):
        path = resolve_crop(row["crop_url"].split("/")[-1])
        if path is None or not path.is_file() or path.stat().st_size == 0:
            raise SystemExit(f"missing crop for {row['scan_start']}")


def write_report(rows, coverage):
    lines = [
        "Milestone 5 — API and map, Starbase Pad 1, GOES-19 ABI band 7",
        "",
        "What this serves",
        "  GET /detections?since=YYYY-MM-DD HH:MM:SS",
        "  Candidates only: status emit or held_back, newest scan start first.",
        "  The score is the value already stored. This process does not recompute it.",
        "  GET /coverage is the latest stored scan and any open gap row.",
        "  GET / is the map. It fetches those two URLs.",
        f"  Listening on http://{HOST}:{PORT}/",
        "",
        "Crops",
        "  A candidate crop is a PNG in outputs/crops/, named from the object key.",
        "  The JSON field is crop_url. Temperatures in the JSON are rounded to 0.01 K.",
        f"  Colour scale {VMIN_K} K to {VMAX_K} K, the Flight 11 figure.",
        "",
        f"Candidates: {len(rows)}",
    ]
    for row in sorted(rows, key=lambda item: item["scan_start"]):
        path = CROP_DIR / crop_png_name(row["source_key"])
        lines.append(
            f"  {row['scan_start']}  {row['status']}  score {row['score_k']:.2f} K"
            f"  max {row['max_bt_k']:.2f} K  {row['distance_km']:.2f} km"
        )
        lines.append(f"    {path.relative_to(ROOT)}  {path.stat().st_size} bytes")
    latest = coverage["latest_scan"]
    lines.append("")
    if latest is None:
        lines.append("Latest stored scan: none")
    else:
        lines.append(
            f"Latest stored scan: {latest['scan_start']}  {latest['sector']}"
            f"  {latest['gap_reason'] or latest['status']}"
        )
    lines.append(f"Open gap rows: {len(coverage['open_gaps'])}")
    lines.extend(
        [
            "",
            "Limits visible on the page",
            "  The 20 K threshold is still the Flight 11 constant.",
            "  A coverage gap is the latest scan when the sector misses the pad.",
            "  The map draws no launch for that row.",
            "  The Cafecito link is a placeholder. The project page is not set up.",
            "  There is no account, no write API, and no alert.",
            "",
            f"Report: {REPORT_PATH.relative_to(ROOT)}",
        ]
    )
    text = "\n".join(lines) + "\n"
    REPORT_PATH.parent.mkdir(parents=True, exist_ok=True)
    REPORT_PATH.write_text(text)
    return text


@app.get("/detections")
def detections(since: str | None = None):
    conn = read_conn()
    try:
        rows = candidate_rows(conn, normalize_since(since))
    finally:
        conn.close()
    return {
        "watch": m4.WATCH["name"],
        "threshold_k": m4.m3.THRESHOLD_K,
        "neighbour_km": m4.m3.NEIGHBOUR_KM,
        "since": normalize_since(since),
        "detections": rows,
    }


@app.get("/coverage")
def coverage():
    conn = read_conn()
    try:
        return coverage_body(conn)
    finally:
        conn.close()


@app.get("/crops/{name}")
def crop_image(name: str):
    path = resolve_crop(name)
    if path is None or not path.is_file():
        raise HTTPException(status_code=404, detail="crop not found")
    return FileResponse(path, media_type="image/png")


@app.get("/")
def map_page():
    return HTMLResponse(PAGE_PATH.read_text())


def main():
    self_test_api()
    if not m4.DB_PATH.exists():
        raise SystemExit(f"no database at {m4.DB_PATH}; run scripts/m4_poll.py first")
    conn = m4.connect(m4.DB_PATH)
    try:
        ensure_crops(conn)
        rows = candidate_rows(conn, None)
        accept_flight11(rows)
        # The since filter keeps the later Flight 11 scan and drops the earlier one.
        later = {row["scan_start"] for row in candidate_rows(conn, "2025-10-13 23:24:00")}
        if FLIGHT11_STARTS[0] in later or FLIGHT11_STARTS[1] not in later:
            raise SystemExit(f"since filter failed: {sorted(later)}")
        text = write_report(rows, coverage_body(conn))
    finally:
        conn.close()
    print(text)
    print(f"listening on http://{HOST}:{PORT}/")
    uvicorn.run(app, host=HOST, port=PORT, log_level="info")


if __name__ == "__main__":
    main()
