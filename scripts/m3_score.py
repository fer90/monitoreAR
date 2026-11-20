#!/usr/bin/env python3
"""Milestone 3: score one watch site and log every candidate.

The watch list is Starbase Pad 1. For each mesoscale band-7 scan that
contains the pad, the score is the hottest pixel within 20 km minus the
median of that same maximum on the previous covering scans. A scan at or
above 20 K is a candidate. The log keeps it either way. It is held back
when the jump has already lasted 3 covering scans, or when the mesoscale
fire product marks a temporally filtered fire within 20 km of the hot pixel.

The 20 K cut and the 20 km radius are frozen from the Milestone 2
measurements. This script does not poll, and it does not serve an API.

Run from the repository root:

    .venv/bin/python scripts/m3_score.py
"""

from __future__ import annotations

import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import s3fs
import xarray as xr

sys.path.insert(0, str(Path(__file__).resolve().parent))
import m1_starship_flight11 as m1
import m2_signal as m2

ROOT = m1.ROOT
BUCKET = m1.BUCKET

# Frozen from Milestone 2. The quiet evening's peak jump was 0.21 K.
# Flight 11 inside 20 km jumped 39.29 K and then 45.41 K. The 8 K pixels
# at +130 s and +190 s were 25 km and 58 km from the pad, outside this radius.
NEIGHBOUR_KM = 20.0
THRESHOLD_K = 20.0
# Five previous scans: one or two hot scans do not move a median of five.
BASELINE_SCANS = 5
MIN_BASELINE = 3
# Flight 11 stayed hot for two scans. The measured fire stayed hot for hours.
PERSIST_SCANS = 3
# Mesoscale cadence at the pad is 60 s. A missed scan breaks the run.
GAP_BREAK_S = 90.0
# DQF 2 out of range, 3 no value, 4 focal-plane temperature. They do not enter.
DQF_EXCLUDE = (2, 3, 4)
# 10–15 are single-scan fire classes. 30–35 are the temporally filtered ones.
# Flight 11's second hot scan is itself class 10, so only 30–35 hold a candidate.
FIRE_INSTANT = set(range(10, 16))
FIRE_TEMPORAL = set(range(30, 36))
FIRE_CODES = FIRE_INSTANT | FIRE_TEMPORAL

WATCH = (
    {
        "name": "starbase-pad-1",
        "lat": m1.PAD_LAT,
        "lon": m1.PAD_LON,
        "armed": True,
        "note": "Starbase Pad 1. Only site armed in this milestone.",
    },
)

OUT_DIR = ROOT / "outputs"
REPORT_PATH = OUT_DIR / "m3_score.txt"
JSONL_PATH = OUT_DIR / "m3_candidates.jsonl"
FIRE_DIR = ROOT / "data" / "m3" / "fire"


def self_test_rule():
    """Lock the baseline arithmetic and the hold rule before any download."""
    pre = [300.14, 300.26, 300.18, 300.06, 300.01]
    if abs(float(np.median(pre)) - 300.14) > 1e-9:
        raise SystemExit("baseline self-test failed on the pre-liftoff median")
    second = [300.26, 300.18, 300.06, 300.01, 339.43]
    if abs(float(np.median(second)) - 300.18) > 1e-9:
        raise SystemExit("baseline self-test failed on the second hot scan")
    if 345.55 - float(np.median(second)) < THRESHOLD_K:
        raise SystemExit("second Flight 11 scan would fall under the frozen threshold")
    if classify(39.29, 1, False) != ("emit", []):
        raise SystemExit("a single jump must be emitted")
    if classify(45.37, 2, False) != ("emit", []):
        raise SystemExit("the second consecutive jump must be emitted")
    if classify(45.0, 2, True) != ("held_back", ["fire_product_temporal"]):
        raise SystemExit("a temporal fire mask must hold the scan back")
    if classify(40.0, 3, False) != ("held_back", ["persistence"]):
        raise SystemExit("the third consecutive jump must be held back")
    if classify(0.21, 1, False) != ("below", []):
        raise SystemExit("the quiet-evening jump must stay below the threshold")
    print("rule self-test: baseline median and hold rule ok")


def classify(score, consecutive, temporal_fire):
    if score is None:
        return "unscored", []
    if score < THRESHOLD_K:
        return "below", []
    reasons = []
    if consecutive >= PERSIST_SCANS:
        reasons.append("persistence")
    if temporal_fire:
        reasons.append("fire_product_temporal")
    return ("held_back" if reasons else "emit"), reasons


def fmt_time(when):
    return when.strftime("%Y-%m-%d %H:%M:%S.") + f"{when.microsecond // 100000:01d}"


def fmt_latlon(lat, lon):
    ns = "N" if lat >= 0 else "S"
    ew = "E" if lon >= 0 else "W"
    return f"{abs(lat):.3f} {ns}, {abs(lon):.3f} {ew}"


def rk(value, digits=2):
    if value is None or not np.isfinite(value):
        return None
    return round(float(value), digits)


def read_watch(path, meta, lat0, lon0):
    """One band-7 file at one watch point. None when the sector misses it."""
    ds = xr.open_dataset(path)
    try:
        platform = str(ds.attrs.get("platform_ID", ""))
        if platform != "G19":
            raise SystemExit(f"{path.name}: platform_ID is {platform!r}, expected G19")
        consts = m1.projection_constants(ds)
        if consts["sweep"] != "x":
            raise SystemExit(f"{path.name}: unexpected sweep axis {consts['sweep']}")
        rad = ds["Rad"]
        fk1 = float(ds["planck_fk1"].values)
        fk2 = float(ds["planck_fk2"].values)
        bc1 = float(ds["planck_bc1"].values)
        bc2 = float(ds["planck_bc2"].values)
        x = np.asarray(ds["x"].values, dtype=np.float64)
        y = np.asarray(ds["y"].values, dtype=np.float64)
        if str(ds["x"].attrs.get("units", "")).lower() not in {"rad", "radians"}:
            raise SystemExit(f"{path.name}: x units are {ds['x'].attrs.get('units')!r}")
        x0, y0 = m1.latlon_to_scan_angles(
            lat0, lon0, consts["r_eq"], consts["r_pol"], consts["h_persp"], consts["lon0"]
        )
        if not (x.min() <= x0 <= x.max() and y.min() <= y0 <= y.max()):
            return None

        ix = int(np.argmin(np.abs(x - x0)))
        iy = int(np.argmin(np.abs(y - y0)))
        r0 = max(0, iy - 30)
        r1 = min(y.size, iy + 31)
        c0 = max(0, ix - 30)
        c1 = min(x.size, ix + 31)
        bt = m1.radiance_to_bt(rad.values[r0:r1, c0:c1], fk1, fk2, bc1, bc2)
        dqf = np.rint(np.asarray(ds["DQF"].values[r0:r1, c0:c1])).astype(np.int16)
        xx, yy = np.meshgrid(x[c0:c1], y[r0:r1])
        lat, lon = m1.scan_angles_to_latlon(
            xx, yy, consts["r_eq"], consts["r_pol"], consts["h_persp"], consts["lon0"]
        )
        dist = m1.haversine_km(lat, lon, lat0, lon0)
        j = iy - r0
        i = ix - c0
        near = dist <= NEIGHBOUR_KM
        ring = (dist > m2.RING_INNER_KM) & (dist <= m2.RING_OUTER_KM)
        usable = np.isfinite(bt) & ~np.isin(dqf, DQF_EXCLUDE)
        finite_near = near & usable
        finite_ring = ring & usable
        near_dqf = dqf[near & np.isfinite(bt)]
        values, counts = np.unique(near_dqf, return_counts=True) if near_dqf.size else ([], [])
        ceiling = m2.bt_ceiling(rad, fk1, fk2, bc1, bc2)

        if finite_near.any():
            chosen = np.where(finite_near, bt, -np.inf)
            j_max, i_max = np.unravel_index(int(np.argmax(chosen)), bt.shape)
            max_bt = float(bt[j_max, i_max])
            max_dist = float(dist[j_max, i_max])
            max_lat = float(lat[j_max, i_max])
            max_lon = float(lon[j_max, i_max])
            max_dqf = int(dqf[j_max, i_max])
            p98 = float(np.nanpercentile(bt[finite_near], 98))
        else:
            max_bt = max_dist = max_lat = max_lon = p98 = float("nan")
            max_dqf = None

        ring_bt = float(np.nanmedian(bt[finite_ring])) if finite_ring.any() else float("nan")
        solar = m1.solar_zenith_deg(
            max_lat if np.isfinite(max_lat) else lat0,
            max_lon if np.isfinite(max_lon) else lon0,
            meta["start"],
        )
        return {
            "meta": meta,
            "near_bt": float(bt[j, i]) if np.isfinite(bt[j, i]) else float("nan"),
            "near_dqf": int(dqf[j, i]),
            "near_dist_km": float(dist[j, i]),
            "max_bt": max_bt,
            "max_dist_km": max_dist,
            "max_lat": max_lat,
            "max_lon": max_lon,
            "max_dqf": max_dqf,
            "p98_bt": p98,
            "ring_bt": ring_bt,
            "n_near": int(finite_near.sum()),
            "n_conditional": int(np.sum(near_dqf == 1)) if near_dqf.size else 0,
            "dqf_counts": {int(v): int(c) for v, c in zip(values, counts)},
            "bt_ceiling": None if ceiling is None else float(ceiling),
            "saturated": bool(
                ceiling is not None and np.isfinite(max_bt) and max_bt >= float(ceiling) - 0.5
            ),
            "solar_zenith_deg": float(solar),
        }
    finally:
        ds.close()


def attach_baseline(rows):
    history = []
    for row in rows:
        window = history[-BASELINE_SCANS:]
        if len(window) >= MIN_BASELINE and np.isfinite(row["max_bt"]):
            baseline = float(np.median(window))
            row["baseline_k"] = baseline
            row["baseline_n"] = len(window)
            row["baseline_samples_k"] = [round(v, 2) for v in window]
            row["score_k"] = float(row["max_bt"] - baseline)
            row["scored"] = True
        else:
            row["baseline_k"] = None
            row["baseline_n"] = len(window)
            row["baseline_samples_k"] = [round(v, 2) for v in window]
            row["score_k"] = None
            row["scored"] = False
        if np.isfinite(row["max_bt"]):
            history.append(float(row["max_bt"]))


def attach_persistence(rows):
    run = 0
    prev = None
    for row in rows:
        start = row["meta"]["start"]
        hot = row["scored"] and row["score_k"] >= THRESHOLD_K
        continuous = prev is not None and (start - prev).total_seconds() <= GAP_BREAK_S
        if hot and continuous:
            run += 1
        elif hot:
            run = 1
        else:
            run = 0
        row["consecutive"] = run
        prev = start


def fire_sector(name):
    if "FDCM1" in name:
        return "M1"
    if "FDCM2" in name:
        return "M2"
    if "FDCC" in name:
        return "CONUS"
    return "?"


def fire_prefix(row):
    sector = row["meta"]["sector"]
    if sector in {"M1", "M2"}:
        product = "ABI-L2-FDCM"
    elif sector == "CONUS":
        product = "ABI-L2-FDCC"
    else:
        return None
    start = row["meta"]["start"]
    doy = int(start.strftime("%j"))
    return f"{BUCKET}/{product}/{start.year}/{doy:03d}/{start.hour:02d}/"


def build_fire_index(fs, prefix):
    found = {}
    for key in fs.ls(prefix):
        name = key.split("/")[-1]
        if "FDCM" not in name and "FDCC" not in name:
            continue
        meta = m1.parse_abi_name(name)
        found[(fire_sector(name), meta["start"])] = key
    return found


def mask_meanings(ds):
    raw_values = ds["Mask"].attrs.get("flag_values", [])
    values = [int(v) for v in np.asarray(raw_values).tolist()]
    words = str(ds["Mask"].attrs.get("flag_meanings", "")).split()
    if len(values) != len(words):
        return {}
    return dict(zip(values, words))


def annotate_fire(fs, row, cache):
    """Record the fire mask near the hot pixel. Return True when codes 30–35 are present."""
    row["fire_file"] = "unavailable"
    row["fire_nearest_mask"] = None
    row["fire_nearest_meaning"] = None
    row["fire_hits"] = []
    prefix = fire_prefix(row)
    if prefix is None or not np.isfinite(row["max_lat"]):
        return False
    if prefix not in cache:
        print(f"  fire index {prefix}")
        cache[prefix] = build_fire_index(fs, prefix)
    key = cache[prefix].get((row["meta"]["sector"], row["meta"]["start"]))
    if key is None:
        return False
    path = m2.download(fs, key, FIRE_DIR)
    row["fire_file"] = path.name
    ds = xr.open_dataset(path)
    try:
        meanings = mask_meanings(ds)
        mask = np.rint(np.asarray(ds["Mask"].values)).astype(np.int16)
        power = np.asarray(ds["Power"].values)
        consts = m1.projection_constants(ds)
        x = np.asarray(ds["x"].values, dtype=np.float64)
        y = np.asarray(ds["y"].values, dtype=np.float64)
        x0, y0 = m1.latlon_to_scan_angles(
            row["max_lat"],
            row["max_lon"],
            consts["r_eq"],
            consts["r_pol"],
            consts["h_persp"],
            consts["lon0"],
        )
        if not (x.min() <= x0 <= x.max() and y.min() <= y0 <= y.max()):
            row["fire_file"] = "unavailable"
            return False
        ix = int(np.argmin(np.abs(x - x0)))
        iy = int(np.argmin(np.abs(y - y0)))
        row["fire_nearest_mask"] = int(mask[iy, ix])
        row["fire_nearest_meaning"] = meanings.get(int(mask[iy, ix]))
        r0 = max(0, iy - 20)
        r1 = min(y.size, iy + 21)
        c0 = max(0, ix - 20)
        c1 = min(x.size, ix + 21)
        sub = mask[r0:r1, c0:c1]
        sub_power = power[r0:r1, c0:c1]
        xx, yy = np.meshgrid(x[c0:c1], y[r0:r1])
        lat, lon = m1.scan_angles_to_latlon(
            xx, yy, consts["r_eq"], consts["r_pol"], consts["h_persp"], consts["lon0"]
        )
        dist = m1.haversine_km(lat, lon, row["max_lat"], row["max_lon"])
        hits = []
        js, is_ = np.where(np.isin(sub, list(FIRE_CODES)) & (dist <= NEIGHBOUR_KM))
        for j_hit, i_hit in zip(js, is_):
            code = int(sub[j_hit, i_hit])
            frp = float(sub_power[j_hit, i_hit])
            hits.append(
                {
                    "code": code,
                    "meaning": meanings.get(code),
                    "dist_km": float(dist[j_hit, i_hit]),
                    "lat": float(lat[j_hit, i_hit]),
                    "lon": float(lon[j_hit, i_hit]),
                    "frp_mw": None if not np.isfinite(frp) else frp,
                }
            )
        hits.sort(key=lambda item: item["dist_km"])
        row["fire_hits"] = hits
        return any(item["code"] in FIRE_TEMPORAL for item in hits)
    finally:
        ds.close()


def finalize(rows, fs, cache):
    attach_baseline(rows)
    attach_persistence(rows)
    for row in rows:
        temporal = False
        row["fire_file"] = None
        row["fire_nearest_mask"] = None
        row["fire_nearest_meaning"] = None
        row["fire_hits"] = []
        if row["scored"] and row["score_k"] >= THRESHOLD_K:
            temporal = annotate_fire(fs, row, cache)
        status, reasons = classify(
            row["score_k"] if row["scored"] else None, row["consecutive"], temporal
        )
        row["status"] = status
        row["hold_reasons"] = reasons


def gather(fs, segments, dest, watch):
    rows = []
    listed = 0
    missed = 0
    for prefix, t0, t1 in segments:
        keys = m2.list_c07(fs, prefix, t0, t1)
        listed += len(keys)
        print(f"{watch['name']} {prefix} {len(keys)} files")
        for key, meta in keys:
            path = m2.download(fs, key, dest)
            print(f"  {path.name}")
            row = read_watch(path, meta, watch["lat"], watch["lon"])
            if row is None:
                missed += 1
                continue
            row["source_key"] = key
            row["watch"] = watch["name"]
            row["armed"] = watch["armed"]
            rows.append(row)
    rows.sort(key=lambda row: row["meta"]["start"])
    return rows, listed, missed


def candidate_record(row):
    hits = []
    for item in row["fire_hits"]:
        hits.append(
            {
                "code": item["code"],
                "meaning": item["meaning"],
                "dist_km": rk(item["dist_km"]),
                "lat": rk(item["lat"], 5),
                "lon": rk(item["lon"], 5),
                "frp_mw": rk(item["frp_mw"]),
            }
        )
    return {
        "watch": row["watch"],
        "armed": row["armed"],
        "status": row["status"],
        "hold_reasons": list(row["hold_reasons"]),
        "scan_start": fmt_time(row["meta"]["start"]),
        "scan_end": fmt_time(row["meta"]["end"]),
        "sector": row["meta"]["sector"],
        "source_key": row["source_key"],
        "pixel_lat": rk(row["max_lat"], 5),
        "pixel_lon": rk(row["max_lon"], 5),
        "distance_km": rk(row["max_dist_km"]),
        "score_k": rk(row["score_k"]),
        "threshold_k": THRESHOLD_K,
        "baseline_k": rk(row["baseline_k"]),
        "baseline_n": row["baseline_n"],
        "baseline_samples_k": row["baseline_samples_k"],
        "near_bt_k": rk(row["near_bt"]),
        "near_dqf": row["near_dqf"],
        "max_bt_k": rk(row["max_bt"]),
        "p98_bt_k": rk(row["p98_bt"]),
        "ring_bt_k": rk(row["ring_bt"]),
        "solar_zenith_deg": rk(row["solar_zenith_deg"], 1),
        "dqf_counts": row["dqf_counts"],
        "hot_pixel_dqf": row["max_dqf"],
        "n_conditional_dqf": row["n_conditional"],
        "n_neighbour_pixels": row["n_near"],
        "saturated": row["saturated"],
        "bt_ceiling_k": rk(row["bt_ceiling"]),
        "consecutive": row["consecutive"],
        "fire_file": row["fire_file"],
        "fire_nearest_mask": row["fire_nearest_mask"],
        "fire_nearest_meaning": row["fire_nearest_meaning"],
        "fire_hits_within_20km": hits,
    }


def write_series(add, rows):
    add(
        "  start_utc                  sector   score_K  base_K  n   max_K"
        "   near_K  dist_km  run  status"
    )
    for row in rows:
        score = f"{row['score_k']:7.2f}" if row["scored"] else "      -"
        base = f"{row['baseline_k']:7.2f}" if row["baseline_k"] is not None else "      -"
        dist = f"{row['max_dist_km']:6.2f}" if np.isfinite(row["max_dist_km"]) else "     -"
        add(
            f"  {fmt_time(row['meta']['start'])}  {row['meta']['sector']:5s}  "
            f"{score}  {base}  {row['baseline_n']:1d}  {row['max_bt']:7.2f}  "
            f"{row['near_bt']:7.2f}  {dist}  {row['consecutive']:3d}  {row['status']}"
        )


def write_candidate(add, row):
    add(f"  {row['watch']}  {row['status']}")
    if row["hold_reasons"]:
        add(f"    hold: {', '.join(row['hold_reasons'])}")
    add(
        f"    scan {fmt_time(row['meta']['start'])} to {fmt_time(row['meta']['end'])}"
        f"  sector {row['meta']['sector']}"
    )
    add(f"    source: {row['source_key']}")
    add(
        f"    hot pixel: {fmt_latlon(row['max_lat'], row['max_lon'])},"
        f" {row['max_dist_km']:.2f} km from the watch point"
    )
    add(
        f"    score: {row['score_k']:.2f} K = {row['max_bt']:.2f} K minus"
        f" baseline {row['baseline_k']:.2f} K over {row['baseline_n']} scans"
        f" {row['baseline_samples_k']}"
    )
    add(f"    pixel nearest the watch point: {row['near_bt']:.2f} K, DQF {row['near_dqf']}")
    add(f"    98th percentile within 20 km: {row['p98_bt']:.2f} K")
    add(f"    median 10–30 km away: {row['ring_bt']:.2f} K")
    add(f"    solar zenith at the hot pixel: {row['solar_zenith_deg']:.1f} deg")
    add(f"    DQF counts within 20 km: {row['dqf_counts']}")
    add(
        f"    hot-pixel DQF: {row['max_dqf']},"
        f" conditional pixels: {row['n_conditional']},"
        f" saturated: {row['saturated']}"
    )
    add(f"    consecutive scans at or above {THRESHOLD_K:.0f} K: {row['consecutive']}")
    add(f"    fire file: {row['fire_file']}")
    if row["fire_nearest_mask"] is not None:
        add(
            f"    fire mask nearest the hot pixel: {row['fire_nearest_mask']}"
            f" {row['fire_nearest_meaning']}"
        )
    if row["fire_hits"]:
        for item in row["fire_hits"]:
            frp = "no FRP" if item["frp_mw"] is None else f"FRP {item['frp_mw']:.1f} MW"
            add(
                f"    fire code {item['code']} {item['meaning']} at"
                f" {item['dist_km']:.2f} km from the hot pixel, {frp}"
            )
    else:
        add("    fire codes 10–15 and 30–35 within 20 km of the hot pixel: none")
    add("")


def peak_row(rows):
    scored = [row for row in rows if row["scored"]]
    if not scored:
        return None
    return max(scored, key=lambda row: row["score_k"])


def write_report(launch, control, fire_rows, counts):
    lines = []
    add = lines.append
    add("Milestone 3 — score, Starbase Pad 1, GOES-19 ABI band 7")
    add("")
    add("Rule, frozen from the Milestone 2 measurements")
    add("  Watch list: starbase-pad-1 at 25.9962 N, 97.1544 W. No other site is armed.")
    add("  Product: mesoscale ABI-L1b-RadM band 7 only. A sector that misses the")
    add("  pad is a coverage gap. CONUS and full disk are not scored for the watch.")
    add("  Neighbourhood: hottest finite pixel within 20 km. DQF 2, 3, and 4 stay out.")
    add("  DQF 1 may enter and is counted.")
    add("  Baseline: median of that maximum on the previous covering scans, up to")
    add("  5, and at least 3. The score does not use the liftoff time.")
    add(f"  Threshold: {THRESHOLD_K:.0f} K. Milestone 2's quiet evening, measured against")
    add("  a fixed pre-liftoff median, peaked at 0.21 K. The two Flight 11 scans")
    add("  inside 20 km were about 39 K and 45 K. The 8 K pixels farther downrange")
    add("  were outside 20 km.")
    add(f"  Persistence: the {PERSIST_SCANS}rd consecutive covering scan at or above the")
    add(f"  threshold is held back. A gap longer than {GAP_BREAK_S:.0f} s starts a new run.")
    add("  Fire product: ABI-L2-FDCM for the same sector and the same scan start.")
    add("  Codes 30–35 within 20 km of the hot pixel hold the candidate back.")
    add("  Codes 10–15 are written in the log and do not hold it back: the second")
    add("  Flight 11 scan is itself code 10. A missing fire file stays in the log")
    add("  as fire_file=unavailable.")
    add("  The 98th percentile is logged and is not the score. On the hot scans it")
    add("  rose by about 5 K, which is under this threshold.")
    add("")
    add("Launch evening — 13 Oct 2025, 23:18–23:36 UTC")
    add(
        f"  files listed: {counts['launch_listed']},"
        f" containing the pad: {len(launch)},"
        f" elsewhere: {counts['launch_missed']}"
    )
    add("")
    write_series(add, launch)
    add("")
    add("Candidates on the launch evening")
    launch_hits = [row for row in launch if row["status"] in {"emit", "held_back"}]
    if not launch_hits:
        add("  none")
        add("")
    for row in launch_hits:
        write_candidate(add, row)
    add("Quiet evening — 25 Oct 2025, same clock window, same rule")
    add(
        f"  files listed: {counts['control_listed']},"
        f" containing the pad: {len(control)},"
        f" elsewhere: {counts['control_missed']}"
    )
    quiet_peak = peak_row(control)
    if quiet_peak is not None:
        add(
            f"  largest score: {quiet_peak['score_k']:.2f} K"
            f" at {fmt_time(quiet_peak['meta']['start'])}"
        )
    add("")
    write_series(add, control)
    add("")
    control_hits = [row for row in control if row["status"] in {"emit", "held_back"}]
    add(f"  candidates: {len(control_hits)}")
    add("")
    add("Ongoing-fire check — not a watch site")
    add("  Same score, applied to the mask-30 pixel Milestone 2 selected from")
    add("  the CONUS fire product. The pixel was already hot at the first scan,")
    add("  so its recent baseline is the fire itself.")
    add(
        f"  point: {fmt_latlon(counts['fire_lat'], counts['fire_lon'])},"
        f" product FRP {counts['fire_frp']:.1f} MW"
    )
    add(
        f"  CONUS files listed: {counts['fire_listed']},"
        f" containing the point: {len(fire_rows)},"
        f" elsewhere: {counts['fire_missed']}"
    )
    fire_peak = peak_row(fire_rows)
    if fire_peak is not None:
        add(
            f"  largest score: {fire_peak['score_k']:.2f} K"
            f" at {fmt_time(fire_peak['meta']['start'])},"
            f" {fire_peak['max_bt']:.2f} K minus baseline {fire_peak['baseline_k']:.2f} K"
        )
    fire_hits = [row for row in fire_rows if row["status"] in {"emit", "held_back"}]
    add(f"  candidates: {len(fire_hits)}")
    add("")
    write_series(add, fire_rows)
    add("")
    if fire_hits:
        for row in fire_hits:
            write_candidate(add, row)
    add("Limits visible in this run")
    add("  The threshold comes from one clear launch evening and one clear quiet")
    add("  evening twelve days later. Clouds were not in either sample.")
    add("  A fire that is already burning does not make a jump against its own")
    add("  previous scans. A fire that ignites while the site is watched can pass")
    add("  for the first two scans; the third is held back.")
    add("  The pad pixel is on the coastline. The fire algorithm can label it")
    add("  water and then miss a plume that has not moved inland.")
    add("  Code 10 on a launch scan is a real fire-product detection of the plume.")
    add("  Holding every fire code back would have dropped Flight 11.")
    add("  Mesoscale tasking still decides whether the minute is scored at all.")
    add("")
    add(f"Log: {JSONL_PATH.relative_to(ROOT)}")
    add(f"Report: {REPORT_PATH.relative_to(ROOT)}")
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    REPORT_PATH.write_text("\n".join(lines) + "\n")
    return "\n".join(lines)


def accept_flight11(launch, control):
    emit = [row for row in launch if row["status"] == "emit"]
    held = [row for row in launch if row["status"] == "held_back"]
    starts = [fmt_time(row["meta"]["start"]) for row in emit]
    expected = ["2025-10-13 23:23:54.9", "2025-10-13 23:24:54.9"]
    if starts != expected or held:
        raise SystemExit(
            "Flight 11 acceptance check failed: "
            f"emitted {starts}, held {len(held)}. See {REPORT_PATH}"
        )
    quiet = [row for row in control if row["status"] in {"emit", "held_back"}]
    if quiet:
        raise SystemExit(f"quiet evening produced {len(quiet)} candidates. See {REPORT_PATH}")
    print("acceptance: Flight 11 emits two scans, the quiet evening emits none")


def main():
    if NEIGHBOUR_KM != m2.NEIGHBOUR_KM:
        raise SystemExit("neighbourhood radius drifted from the Milestone 2 measurement")
    self_test_rule()
    m1.self_test_geolocation()
    fs = s3fs.S3FileSystem(anon=True, client_kwargs={"region_name": "us-east-1"})
    cache = {}
    watch = WATCH[0]

    launch_t0, launch_t1 = m2.clock_window(m1.LIFTOFF)
    launch, launch_listed, launch_missed = gather(
        fs,
        [(f"{BUCKET}/ABI-L1b-RadM/2025/286/23/", launch_t0, launch_t1)],
        m2.LAUNCH_DIR,
        watch,
    )
    finalize(launch, fs, cache)

    control_t0, control_t1 = m2.clock_window(m2.CONTROL_DATE)
    control, control_listed, control_missed = gather(
        fs,
        [(f"{BUCKET}/ABI-L1b-RadM/2025/{m2.CONTROL_DOY:03d}/23/", control_t0, control_t1)],
        m2.CONTROL_DIR,
        watch,
    )
    finalize(control, fs, cache)

    _path, fire, _top = m2.select_fire(fs)
    fire_watch = {
        "name": "m2-fdcc-fire-check",
        "lat": fire["lat"],
        "lon": fire["lon"],
        "armed": False,
        "note": "Ongoing fire used as a negative check. Not an armed watch site.",
    }
    segments = []
    for hour in m2.FIRE_HOURS:
        segments.append(
            (
                f"{BUCKET}/ABI-L1b-RadC/2025/286/{hour}/",
                datetime(2025, 10, 13, int(hour), 0, tzinfo=timezone.utc),
                datetime(2025, 10, 13, int(hour), 59, 59, tzinfo=timezone.utc),
            )
        )
    fire_rows, fire_listed, fire_missed = gather(fs, segments, m2.FIRE_DIR, fire_watch)
    finalize(fire_rows, fs, cache)

    records = [
        candidate_record(row)
        for row in launch + control + fire_rows
        if row["status"] in {"emit", "held_back"}
    ]
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    JSONL_PATH.write_text("".join(json.dumps(record) + "\n" for record in records))
    text = write_report(
        launch,
        control,
        fire_rows,
        {
            "launch_listed": launch_listed,
            "launch_missed": launch_missed,
            "control_listed": control_listed,
            "control_missed": control_missed,
            "fire_listed": fire_listed,
            "fire_missed": fire_missed,
            "fire_lat": fire["lat"],
            "fire_lon": fire["lon"],
            "fire_frp": fire["frp_mw"],
        },
    )
    print(text)
    accept_flight11(launch, control)


if __name__ == "__main__":
    main()
