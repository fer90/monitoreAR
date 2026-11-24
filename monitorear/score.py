"""The watch-site score.

The score is the hottest usable band-7 pixel within 20 km minus the
median of that maximum on the previous covering scans. The 20 K cut and
the 20 km radius are the constants frozen from the Flight 11 measurements.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import xarray as xr

from monitorear.calibrate import bt_ceiling, radiance_to_bt
from monitorear.geolocate import (
    PAD_LAT,
    PAD_LON,
    haversine_km,
    latlon_to_scan_angles,
    projection_constants,
    scan_angles_to_latlon,
    solar_zenith_deg,
)
from monitorear.goes_io import BUCKET, ROOT, download, parse_abi_name

# Frozen from Milestone 2. The quiet evening's peak jump was 0.21 K.
# Flight 11 inside 20 km jumped 39.29 K and then 45.41 K. The 8 K pixels
# at +130 s and +190 s were 25 km and 58 km from the pad, outside this radius.
NEIGHBOUR_KM = 20.0
# Ring used as a logged raw value. It is not the score.
RING_INNER_KM = 10.0
RING_OUTER_KM = 30.0
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
        "lat": PAD_LAT,
        "lon": PAD_LON,
        "armed": True,
        "note": "Starbase Pad 1. Only site armed in this milestone.",
    },
)

FIRE_DIR = ROOT / "data" / "m3" / "fire"
OUT_DIR = ROOT / "outputs"


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
        consts = projection_constants(ds)
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
        x0, y0 = latlon_to_scan_angles(
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
        bt = radiance_to_bt(rad.values[r0:r1, c0:c1], fk1, fk2, bc1, bc2)
        dqf = np.rint(np.asarray(ds["DQF"].values[r0:r1, c0:c1])).astype(np.int16)
        xx, yy = np.meshgrid(x[c0:c1], y[r0:r1])
        lat, lon = scan_angles_to_latlon(
            xx, yy, consts["r_eq"], consts["r_pol"], consts["h_persp"], consts["lon0"]
        )
        dist = haversine_km(lat, lon, lat0, lon0)
        j = iy - r0
        i = ix - c0
        near = dist <= NEIGHBOUR_KM
        ring = (dist > RING_INNER_KM) & (dist <= RING_OUTER_KM)
        usable = np.isfinite(bt) & ~np.isin(dqf, DQF_EXCLUDE)
        finite_near = near & usable
        finite_ring = ring & usable
        near_dqf = dqf[near & np.isfinite(bt)]
        values, counts = np.unique(near_dqf, return_counts=True) if near_dqf.size else ([], [])
        ceiling = bt_ceiling(rad, fk1, fk2, bc1, bc2)

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
        solar = solar_zenith_deg(
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
        meta = parse_abi_name(name)
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
    path = download(fs, key, FIRE_DIR)
    row["fire_file"] = path.name
    ds = xr.open_dataset(path)
    try:
        meanings = mask_meanings(ds)
        mask = np.rint(np.asarray(ds["Mask"].values)).astype(np.int16)
        power = np.asarray(ds["Power"].values)
        consts = projection_constants(ds)
        x = np.asarray(ds["x"].values, dtype=np.float64)
        y = np.asarray(ds["y"].values, dtype=np.float64)
        x0, y0 = latlon_to_scan_angles(
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
        lat, lon = scan_angles_to_latlon(
            xx, yy, consts["r_eq"], consts["r_pol"], consts["h_persp"], consts["lon0"]
        )
        dist = haversine_km(lat, lon, row["max_lat"], row["max_lon"])
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


def m3_baseline(samples):
    """Median of the recent neighbourhood maxima. None until MIN_BASELINE exist."""
    if len(samples) < MIN_BASELINE:
        return None, len(samples), []
    window = list(samples[-BASELINE_SCANS:])
    return float(np.median(window)), len(window), window
