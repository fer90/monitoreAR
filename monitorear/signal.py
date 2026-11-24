#!/usr/bin/env python3
"""Milestone 2: quantify the Flight 11 band-7 signal.

Same statistics at three places:

- Starbase Pad 1 on 13 Oct 2025 (the launch), mesoscale band 7.
- The same pad on 25 Oct 2025, the nearest evening on which a mesoscale
  sector still contained the pad at this clock time.
- The strongest good-fire pixel in the GOES-19 CONUS fire product for
  the scan that overlaps liftoff, measured in CONUS band 7.

The script reports max brightness temperature, a high percentile, and
scan-to-scan differences. It does not choose a detection threshold.

Run from the repository root:

    .venv/bin/python scripts/m2_signal.py
"""

from __future__ import annotations

from monitorear import flight11 as m1
from monitorear.calibrate import bt_ceiling
from monitorear.goes_io import download, list_c07, sector_label
from monitorear.score import NEIGHBOUR_KM, RING_INNER_KM, RING_OUTER_KM

from datetime import datetime, timedelta, timezone
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import s3fs


ROOT = m1.ROOT
BUCKET = m1.BUCKET
PAD_LAT = m1.PAD_LAT
PAD_LON = m1.PAD_LON
LIFTOFF = m1.LIFTOFF

# 18 minutes centred on the launch clock. The control day uses the same
# clock time on its own date.
WINDOW_START_CLOCK = (23, 18)
WINDOW_END_CLOCK = (23, 36)

# Mesoscale sectors on the evenings of 27 Sep–23 Oct and 26–28 Oct 2025
# were aimed elsewhere at 23:22 UTC (both sectors checked). The nearest
# later evening with a sector on the pad is 25 Oct, sector M1. Solar
# zenith that evening is about 84.6°, against 82.1° on launch evening.
CONTROL_DATE = datetime(2025, 10, 25, tzinfo=timezone.utc)
CONTROL_DOY = 298

# CONUS fire product whose scan overlaps the liftoff (23:21:17–23:23:54).
# Mask 30 is temporally_filtered_good_fire_pixel. The launch pad is
# excluded so this scan cannot select the rocket.
FIRE_PRODUCT_KEY = (
    f"{BUCKET}/ABI-L2-FDCC/2025/286/23/"
    "OR_ABI-L2-FDCC-M6_G19_s20252862321170_e20252862323543_c20252862324125.nc"
)
FIRE_MASK = 30
FIRE_PAD_EXCLUSION_KM = 100.0
# Hours of CONUS band 7 used to see whether that pixel stays hot.
FIRE_HOURS = ("21", "22", "23")

# Spatial background sits outside the hot pixel and its immediate neighbours.

LAUNCH_DIR = ROOT / "data" / "m1"
CONTROL_DIR = ROOT / "data" / "m2" / "control"
FIRE_DIR = ROOT / "data" / "m2" / "fire"
OUT_DIR = ROOT / "outputs"
REPORT_PATH = OUT_DIR / "m2_signal.txt"
FIGURE_PATH = OUT_DIR / "m2_signal.png"


def clock_window(day):
    """Return [start, end] on `day` for the 23:18–23:36 UTC window."""
    start = day.replace(
        hour=WINDOW_START_CLOCK[0],
        minute=WINDOW_START_CLOCK[1],
        second=0,
        microsecond=0,
    )
    end = day.replace(
        hour=WINDOW_END_CLOCK[0],
        minute=WINDOW_END_CLOCK[1],
        second=0,
        microsecond=0,
    )
    return start, end




def measure_point(path, meta, lat0, lon0):
    """Statistics of one band-7 file at one surface point.

    Returns None when the point lies outside this file's scan-angle box.
    Brightness temperature uses the Planck coefficients in the file.
    """
    import xarray as xr

    ds = xr.open_dataset(path)
    try:
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
            raise SystemExit(f"{path.name}: x units are {ds['x'].attrs.get('units')!r}, expected rad")
        x0, y0 = m1.latlon_to_scan_angles(
            lat0, lon0, consts["r_eq"], consts["r_pol"], consts["h_persp"], consts["lon0"]
        )
        if not (x.min() <= x0 <= x.max() and y.min() <= y0 <= y.max()):
            return None

        ix = int(np.argmin(np.abs(x - x0)))
        iy = int(np.argmin(np.abs(y - y0)))
        # ±30 pixels reaches about 70 km, past the 30 km ring and past the
        # warm pixels Milestone 1 saw 25 km and 58 km from the pad.
        r0 = max(0, iy - 30)
        r1 = min(y.size, iy + 31)
        c0 = max(0, ix - 30)
        c1 = min(x.size, ix + 31)
        bt = m1.radiance_to_bt(rad.values[r0:r1, c0:c1], fk1, fk2, bc1, bc2)
        dqf = np.asarray(ds["DQF"].values[r0:r1, c0:c1])
        xx, yy = np.meshgrid(x[c0:c1], y[r0:r1])
        lat, lon = m1.scan_angles_to_latlon(
            xx, yy, consts["r_eq"], consts["r_pol"], consts["h_persp"], consts["lon0"]
        )
        dist = m1.haversine_km(lat, lon, lat0, lon0)
        j = iy - r0
        i = ix - c0
        near = dist <= NEIGHBOUR_KM
        ring = (dist > RING_INNER_KM) & (dist <= RING_OUTER_KM)
        finite_near = near & np.isfinite(bt)
        finite_ring = ring & np.isfinite(bt)

        if finite_near.any():
            masked = np.where(finite_near, bt, -np.inf)
            j_max, i_max = np.unravel_index(int(np.argmax(masked)), bt.shape)
            max_bt = float(bt[j_max, i_max])
            max_dist = float(dist[j_max, i_max])
            max_lat = float(lat[j_max, i_max])
            max_lon = float(lon[j_max, i_max])
            p98 = float(np.nanpercentile(bt[finite_near], 98))
            med_near = float(np.nanmedian(bt[finite_near]))
        else:
            max_bt = max_dist = max_lat = max_lon = p98 = med_near = float("nan")

        ring_med = float(np.nanmedian(bt[finite_ring])) if finite_ring.any() else float("nan")
        outer = (dist > NEIGHBOUR_KM) & np.isfinite(bt)
        if outer.any():
            j_out, i_out = np.unravel_index(int(np.argmax(np.where(outer, bt, -np.inf))), bt.shape)
            outer_bt = float(bt[j_out, i_out])
            outer_dist = float(dist[j_out, i_out])
        else:
            outer_bt = outer_dist = float("nan")
        values, counts = np.unique(dqf[np.isfinite(bt)], return_counts=True)
        pix_dx = float("nan")
        if i + 1 < bt.shape[1]:
            pix_dx = float(m1.haversine_km(lat[j, i], lon[j, i], lat[j, i + 1], lon[j, i + 1]))
        pix_dy = float("nan")
        if j + 1 < bt.shape[0]:
            pix_dy = float(m1.haversine_km(lat[j, i], lon[j, i], lat[j + 1, i], lon[j + 1, i]))

        return {
            "path": path,
            "meta": meta,
            "near_bt": float(bt[j, i]),
            "near_lat": float(lat[j, i]),
            "near_lon": float(lon[j, i]),
            "near_dist_km": float(dist[j, i]),
            "max20_bt": max_bt,
            "max20_dist_km": max_dist,
            "max20_lat": max_lat,
            "max20_lon": max_lon,
            "p98_20_bt": p98,
            "median20_bt": med_near,
            "ring_bt": ring_med,
            "outer_bt": outer_bt,
            "outer_dist_km": outer_dist,
            "n_ring": int(finite_ring.sum()),
            "n_near": int(finite_near.sum()),
            "dqf_counts": {int(v): int(c) for v, c in zip(values, counts)},
            "pixel_dx_km": pix_dx,
            "pixel_dy_km": pix_dy,
            "bt_ceiling": bt_ceiling(rad, fk1, fk2, bc1, bc2),
            "planck": (fk1, fk2, bc1, bc2),
            "band_wavelength_um": float(np.asarray(ds["band_wavelength"].values).reshape(-1)[0]),
            "consts": consts,
            "subpoint_lon": float(ds["nominal_satellite_subpoint_lon"].values),
        }
    finally:
        ds.close()


def collect(fs, keys, dest_dir, lat0, lon0):
    rows = []
    missed = 0
    for key, meta in keys:
        path = download(fs, key, dest_dir)
        print(f"  {path.name}")
        row = measure_point(path, meta, lat0, lon0)
        if row is None:
            missed += 1
            continue
        rows.append(row)
    rows.sort(key=lambda row: row["meta"]["start"])
    for prev, row in zip(rows, rows[1:]):
        row["d_near"] = row["near_bt"] - prev["near_bt"]
        row["d_max20"] = row["max20_bt"] - prev["max20_bt"]
    if rows:
        rows[0]["d_near"] = float("nan")
        rows[0]["d_max20"] = float("nan")
    return rows, missed


def select_fire(fs):
    """Hottest mask-30 pixel in the pinned CONUS fire product, away from the pad."""
    import xarray as xr

    path = download(fs, FIRE_PRODUCT_KEY, FIRE_DIR)
    ds = xr.open_dataset(path)
    try:
        mask = np.asarray(ds["Mask"].values)
        power = np.asarray(ds["Power"].values)
        temp = np.asarray(ds["Temp"].values)
        area = np.asarray(ds["Area"].values)
        consts = m1.projection_constants(ds)
        x = np.asarray(ds["x"].values, dtype=np.float64)
        y = np.asarray(ds["y"].values, dtype=np.float64)
        ys, xs = np.where(mask == FIRE_MASK)
        ranked = []
        for j, i in zip(ys, xs):
            lat, lon = m1.scan_angles_to_latlon(
                x[int(i)], y[int(j)], consts["r_eq"], consts["r_pol"], consts["h_persp"], consts["lon0"]
            )
            lat, lon = float(lat), float(lon)
            if not np.isfinite(lat):
                continue
            dist = float(m1.haversine_km(PAD_LAT, PAD_LON, lat, lon))
            frp = float(power[j, i]) if np.isfinite(power[j, i]) else float("nan")
            ranked.append(
                {
                    "lat": lat,
                    "lon": lon,
                    "dist_pad_km": dist,
                    "frp_mw": frp,
                    "fire_temp_k": float(temp[j, i]) if np.isfinite(temp[j, i]) else float("nan"),
                    "area_m2": float(area[j, i]) if np.isfinite(area[j, i]) else float("nan"),
                    "mask": FIRE_MASK,
                }
            )
    finally:
        ds.close()
    ranked.sort(key=lambda item: -np.nan_to_num(item["frp_mw"], nan=-1.0))
    eligible = [item for item in ranked if item["dist_pad_km"] > FIRE_PAD_EXCLUSION_KM and np.isfinite(item["frp_mw"])]
    if not eligible:
        raise SystemExit("fire product has no mask-30 pixel far from the pad")
    return path, eligible[0], ranked[:5]


def _fmt_time(when):
    return when.strftime("%Y-%m-%d %H:%M:%S.") + f"{when.microsecond // 100000:01d}"


def _fmt_latlon(lat, lon):
    ns = "N" if lat >= 0 else "S"
    ew = "E" if lon >= 0 else "W"
    return f"{abs(lat):.3f} {ns}, {abs(lon):.3f} {ew}"


def _seconds_after(when, anchor):
    return (when - anchor).total_seconds()


def summarize_pad(rows, anchor):
    """Baseline from scans that ended at or before `anchor`, then jumps."""
    pre = [row for row in rows if row["meta"]["end"] <= anchor]
    if len(pre) < 3:
        raise SystemExit(f"need at least 3 pre-anchor scans, found {len(pre)}")
    baseline = float(np.median([row["max20_bt"] for row in pre]))
    baseline_near = float(np.median([row["near_bt"] for row in pre]))
    for row in rows:
        row["jump20"] = row["max20_bt"] - baseline
        row["jump_near"] = row["near_bt"] - baseline_near
    peak = max(rows, key=lambda row: row["jump20"])
    return {
        "n": len(rows),
        "n_pre": len(pre),
        "baseline20": baseline,
        "baseline_near": baseline_near,
        "pre_std20": float(np.std([row["max20_bt"] for row in pre], ddof=1)),
        "pre_std_near": float(np.std([row["near_bt"] for row in pre], ddof=1)),
        "peak": peak,
        "peak_jump": peak["jump20"],
    }


def quiet_noise(rows):
    values = np.array([row["max20_bt"] for row in rows], dtype=np.float64)
    steps = np.array([row["d_max20"] for row in rows[1:]], dtype=np.float64)
    near = np.array([row["near_bt"] for row in rows], dtype=np.float64)
    return {
        "n": len(rows),
        "median20": float(np.median(values)),
        "std20": float(np.std(values, ddof=1)) if len(values) > 1 else float("nan"),
        "min20": float(np.min(values)),
        "max20": float(np.max(values)),
        "median_near": float(np.median(near)),
        "std_near": float(np.std(near, ddof=1)) if len(near) > 1 else float("nan"),
        "min_near": float(np.min(near)),
        "max_near": float(np.max(near)),
        "max_step": float(np.max(steps)) if steps.size else float("nan"),
        "min_step": float(np.min(steps)) if steps.size else float("nan"),
        "std_step": float(np.std(steps, ddof=1)) if steps.size > 1 else float("nan"),
    }


def gaps_seconds(rows):
    starts = [row["meta"]["start"] for row in rows]
    return [(b - a).total_seconds() for a, b in zip(starts, starts[1:])]


def write_series(add, rows, anchor):
    add(
        "  start_utc                  sector   rel_s   near_K   p98_K  max20_K"
        "   dmax_K   ring_K  dist20"
    )
    for row in rows:
        rel = _seconds_after(row["meta"]["start"], anchor)
        dmax = row["d_max20"]
        dmax_s = f"{dmax:7.2f}" if np.isfinite(dmax) else "      -"
        add(
            f"  {_fmt_time(row['meta']['start'])}  {row['meta']['sector']:5s}  "
            f"{rel:7.1f}  {row['near_bt']:7.2f}  {row['p98_20_bt']:6.2f}  "
            f"{row['max20_bt']:7.2f}  {dmax_s}  {row['ring_bt']:7.2f}  "
            f"{row['max20_dist_km']:6.2f}"
        )


def plot_series(launch, control, fire, launch_anchor, control_anchor, snr_line):
    fig, axes = plt.subplots(1, 2, figsize=(11.8, 5.0))

    def minutes(rows, anchor_clock):
        # Minutes after 23:18 UTC on that row's own date.
        return [
            (
                row["meta"]["start"]
                - anchor_clock.replace(
                    hour=23, minute=18, second=0, microsecond=0
                )
            ).total_seconds()
            / 60.0
            for row in rows
        ]

    ax = axes[0]
    lx = minutes(launch, launch_anchor)
    cx = minutes(control, control_anchor)
    ax.plot(lx, [r["max20_bt"] for r in launch], color="#b22222", lw=1.8, marker="o", ms=3.5, label="launch, max within 20 km")
    ax.plot(lx, [r["near_bt"] for r in launch], color="#b22222", lw=1.0, ls="--", marker=".", ms=3, label="launch, pad pixel")
    ax.plot(cx, [r["max20_bt"] for r in control], color="#4d4d4d", lw=1.4, marker="o", ms=3.5, label="25 Oct, max within 20 km")
    ax.plot(cx, [r["near_bt"] for r in control], color="#4d4d4d", lw=1.0, ls="--", marker=".", ms=3, label="25 Oct, pad pixel")
    liftoff_min = (LIFTOFF - launch_anchor.replace(hour=23, minute=18, second=0, microsecond=0)).total_seconds() / 60.0
    ax.axvline(liftoff_min, color="#1f4e79", lw=0.8, ls=":", label="liftoff 23:23:45")
    ax.set_xlabel("Minutes after 23:18 UTC")
    ax.set_ylabel("Brightness temperature (K)")
    ax.set_title("Pad 1, mesoscale band 7")
    ax.grid(True, color="0.85", lw=0.4)
    ax.legend(fontsize=7, loc="upper left")

    ax = axes[1]
    times = [r["meta"]["start"] for r in fire]
    # On every scan in this window the nearest pixel is also the 20 km maximum.
    ax.plot(times, [r["near_bt"] for r in fire], color="#b22222", lw=1.6, marker="o", ms=3.5, label="fire pixel")
    ax.plot(times, [r["ring_bt"] for r in fire], color="#4d4d4d", lw=1.2, label="median 10–30 km away")
    ax.set_ylabel("Brightness temperature (K)")
    ax.set_title("Washington County fire, CONUS band 7")
    ax.grid(True, color="0.85", lw=0.4)
    ax.legend(fontsize=7, loc="best")
    ax.set_xlabel("UTC on 13 Oct 2025")
    ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter("%H:%M", tz=timezone.utc))

    fig.suptitle("GOES-19 ABI band 7 — Flight 11 signal against a quiet evening and a fire", fontsize=11)
    fig.autofmt_xdate()
    fig.subplots_adjust(left=0.07, right=0.98, bottom=0.20, top=0.86, wspace=0.28)
    fig.text(0.07, 0.02, snr_line, fontsize=6.5, ha="left", va="bottom")
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIGURE_PATH, dpi=140)
    plt.close(fig)


def main():
    m1.self_test_geolocation()
    fs = s3fs.S3FileSystem(anon=True, client_kwargs={"region_name": "us-east-1"})

    launch_t0, launch_t1 = clock_window(LIFTOFF)
    print("launch window", launch_t0, launch_t1)
    launch_keys = list_c07(fs, f"{BUCKET}/ABI-L1b-RadM/2025/286/23/", launch_t0, launch_t1)
    launch_rows, launch_miss = collect(fs, launch_keys, LAUNCH_DIR, PAD_LAT, PAD_LON)

    control_day = CONTROL_DATE
    control_t0, control_t1 = clock_window(control_day)
    # Anchor at the same clock time as liftoff, so scan-relative seconds match.
    control_anchor = control_day.replace(
        hour=LIFTOFF.hour, minute=LIFTOFF.minute, second=LIFTOFF.second, microsecond=0
    )
    print("control window", control_t0, control_t1)
    control_keys = list_c07(
        fs, f"{BUCKET}/ABI-L1b-RadM/2025/{CONTROL_DOY:03d}/23/", control_t0, control_t1
    )
    control_rows, control_miss = collect(fs, control_keys, CONTROL_DIR, PAD_LAT, PAD_LON)
    if len(control_rows) < 5:
        raise SystemExit(
            f"control evening produced {len(control_rows)} frames on the pad "
            f"({control_miss} files missed it). The sector was not on Starbase."
        )

    print("fire product", FIRE_PRODUCT_KEY.split("/")[-1])
    fire_path, fire, fire_top = select_fire(fs)
    print(
        f"  selected {fire['lat']:.3f} N, {fire['lon']:.3f} E, "
        f"FRP {fire['frp_mw']:.1f} MW, {fire['dist_pad_km']:.0f} km from the pad"
    )
    fire_rows = []
    fire_listed = 0
    fire_miss = 0
    for hour in FIRE_HOURS:
        keys = list_c07(
            fs,
            f"{BUCKET}/ABI-L1b-RadC/2025/286/{hour}/",
            datetime(2025, 10, 13, int(hour), 0, tzinfo=timezone.utc),
            datetime(2025, 10, 13, int(hour), 59, 59, tzinfo=timezone.utc),
        )
        fire_listed += len(keys)
        print(f"fire hour {hour}: {len(keys)} CONUS band-7 files")
        rows, missed = collect(fs, keys, FIRE_DIR, fire["lat"], fire["lon"])
        fire_rows.extend(rows)
        fire_miss += missed
    fire_rows.sort(key=lambda row: row["meta"]["start"])
    for prev, row in zip(fire_rows, fire_rows[1:]):
        row["d_near"] = row["near_bt"] - prev["near_bt"]
        row["d_max20"] = row["max20_bt"] - prev["max20_bt"]
    if fire_rows:
        fire_rows[0]["d_near"] = float("nan")
        fire_rows[0]["d_max20"] = float("nan")

    launch_stats = summarize_pad(launch_rows, LIFTOFF)
    # Same clock rule on the control day: scans that ended before 23:23:45.
    control_pre = summarize_pad(control_rows, control_anchor)
    noise = quiet_noise(control_rows)
    peak_jump = launch_stats["peak_jump"]
    snr_std = peak_jump / noise["std20"] if noise["std20"] and noise["std20"] > 0 else float("nan")
    snr_step = peak_jump / noise["max_step"] if noise["max_step"] and noise["max_step"] > 0 else float("nan")

    # Tail: launch scans after the peak whose jump is still positive.
    launch_jumps = sorted(((row["jump20"], row) for row in launch_rows), key=lambda item: -item[0])

    fire_anom = np.array([row["max20_bt"] - row["ring_bt"] for row in fire_rows], dtype=np.float64)
    fire_peak_i = int(np.nanargmax(fire_anom)) if fire_rows else 0
    fire_peak = fire_rows[fire_peak_i] if fire_rows else None
    hot_fire = [row for row in fire_rows if row["max20_bt"] - row["ring_bt"] > 10.0]
    if hot_fire:
        fire_span_min = (
            hot_fire[-1]["meta"]["start"] - hot_fire[0]["meta"]["start"]
        ).total_seconds() / 60.0
    else:
        fire_span_min = 0.0

    sun_launch = m1.solar_zenith_deg(PAD_LAT, PAD_LON, LIFTOFF)
    sun_control = m1.solar_zenith_deg(PAD_LAT, PAD_LON, control_anchor)
    sun_fire = m1.solar_zenith_deg(fire["lat"], fire["lon"], LIFTOFF)
    consts = launch_rows[0]["consts"]
    zenith_pad = m1.satellite_zenith_deg(
        PAD_LAT, PAD_LON, consts["r_eq"], consts["r_pol"], consts["h_persp"], consts["lon0"]
    )
    fconsts = fire_rows[0]["consts"]
    zenith_fire = m1.satellite_zenith_deg(
        fire["lat"], fire["lon"], fconsts["r_eq"], fconsts["r_pol"], fconsts["h_persp"], fconsts["lon0"]
    )

    lines = []
    add = lines.append
    add("Milestone 2 — signal comparison, Starship Flight 11, GOES-19 ABI band 7")
    add("")
    add("What is being compared")
    add("  Launch: mesoscale band 7 at Pad 1, 13 Oct 2025, 23:18–23:36 UTC.")
    add("  Quiet day: mesoscale band 7 at the same pad and clock time on 25 Oct 2025.")
    add("  Fire: CONUS band 7 at the strongest mask-30 pixel in the CONUS fire")
    add("  product for the scan starting 23:21:17 UTC on 13 Oct 2025.")
    add("  Brightness temperature uses planck_fk1, planck_fk2, planck_bc1, planck_bc2")
    add("  from each file. Reflected sunlight stays in the number.")
    add("")
    add("Statistics")
    add("  near_K     brightness temperature of the pixel nearest the point")
    add("  p98_K      98th percentile of finite pixels within 20 km")
    add("  max20_K    hottest finite pixel whose centre is within 20 km")
    add("  dmax_K     max20_K minus max20_K of the previous scan that saw the point")
    add("  ring_K     median brightness temperature of pixels 10–30 km from the point")
    add("  jump       max20_K minus the median max20_K of scans that ended by 23:23:45")
    add("  The 20 km test uses a sphere of radius 6371.0088 km, as in Milestone 1.")
    add("  No detection threshold is chosen here.")
    add("")
    add("Launch day")
    add(f"  files in the window: {len(launch_keys)}, containing the pad: {len(launch_rows)}, elsewhere: {launch_miss}")
    add(f"  solar zenith at liftoff: {sun_launch:.1f} deg")
    add(f"  satellite zenith at the pad: {zenith_pad:.1f} deg")
    add(
        f"  pixel spacing at the pad: {launch_rows[0]['pixel_dx_km']:.2f} km east-west, "
        f"{launch_rows[0]['pixel_dy_km']:.2f} km north-south"
    )
    add(
        f"  pre-liftoff scans: {launch_stats['n_pre']}, "
        f"median max-within-20 km {launch_stats['baseline20']:.2f} K, "
        f"sample std {launch_stats['pre_std20']:.2f} K"
    )
    add(
        f"  pre-liftoff pad pixel median {launch_stats['baseline_near']:.2f} K, "
        f"sample std {launch_stats['pre_std_near']:.2f} K"
    )
    peak = launch_stats["peak"]
    add(
        f"  peak jump: {peak_jump:.2f} K at {_fmt_time(peak['meta']['start'])}, "
        f"max within 20 km {peak['max20_bt']:.2f} K, "
        f"{peak['max20_dist_km']:.2f} km from the pad "
        f"({_fmt_latlon(peak['max20_lat'], peak['max20_lon'])})"
    )
    ceiling = launch_rows[0]["bt_ceiling"]
    add(f"  scaled brightness-temperature ceiling: {ceiling:.2f} K")
    add(f"  wavelength: {launch_rows[0]['band_wavelength_um']:.3f} µm")
    gaps = gaps_seconds(launch_rows)
    if gaps:
        add(
            f"  gap between pad-seeing scans: min {min(gaps):.1f} s, "
            f"median {float(np.median(gaps)):.1f} s, max {max(gaps):.1f} s"
        )
    add("")
    write_series(add, launch_rows, LIFTOFF)
    add("")
    add("  Jumps of the 20 km maximum above the pre-liftoff median, largest first")
    for jump, row in launch_jumps:
        if jump < 2.0:
            continue
        rel = _seconds_after(row["meta"]["start"], LIFTOFF)
        add(
            f"    {jump:6.2f} K at {rel:+7.1f} s  "
            f"({row['max20_bt']:.2f} K, {row['max20_dist_km']:.1f} km from the pad)"
        )
    hot_scans = [row for row in launch_rows if row["meta"]["start"] >= LIFTOFF][:2]
    if hot_scans:
        p98_peak = max(row["p98_20_bt"] for row in hot_scans)
        add(
            f"  On those two scans the 98th percentile within 20 km reaches {p98_peak:.2f} K. "
            f"On the last pre-liftoff scan it is {launch_rows[launch_stats['n_pre'] - 1]['p98_20_bt']:.2f} K."
        )
        add("  One hot pixel moves the maximum by tens of kelvin and the 98th percentile by about 5 K.")
    add("  Hottest pixel farther than 20 km, inside the ±30-pixel window:")
    outer_hits = [
        row for row in launch_rows if row["outer_bt"] > launch_stats["baseline20"] + 5.0
    ]
    if not outer_hits:
        add("    none more than 5 K above the pre-liftoff neighbourhood median")
    for row in outer_hits:
        rel = _seconds_after(row["meta"]["start"], LIFTOFF)
        add(
            f"    {row['outer_bt']:.2f} K at {rel:+7.1f} s, "
            f"{row['outer_dist_km']:.1f} km from the pad "
            f"({row['outer_bt'] - launch_stats['baseline20']:.2f} K above the 20 km baseline, "
            f"{(row['outer_bt'] - launch_stats['baseline20']) / noise['std20']:.0f} quiet-day stds)"
        )
    add("")
    add("Quiet day — 25 Oct 2025, same clock window")
    add("  Adjacent evenings were checked from 27 Sep through 28 Oct, one band-7")
    add("  file per sector nearest 23:22 UTC. Neither sector contained the pad on")
    add("  12 Oct or 14 Oct. 25 Oct sector M1 did. Solar zenith at the pad at")
    add(f"  23:23:45 UTC was {sun_control:.1f} deg, against {sun_launch:.1f} deg on launch evening.")
    add(f"  files in the window: {len(control_keys)}, containing the pad: {len(control_rows)}, elsewhere: {control_miss}")
    add(
        f"  neighbourhood max: median {noise['median20']:.2f} K, "
        f"std {noise['std20']:.2f} K, min {noise['min20']:.2f} K, max {noise['max20']:.2f} K"
    )
    add(
        f"  pad pixel: median {noise['median_near']:.2f} K, "
        f"std {noise['std_near']:.2f} K, min {noise['min_near']:.2f} K, max {noise['max_near']:.2f} K"
    )
    add(
        f"  scan-to-scan change of the neighbourhood max: "
        f"largest rise {noise['max_step']:.2f} K, largest fall {noise['min_step']:.2f} K, "
        f"std {noise['std_step']:.2f} K"
    )
    cgaps = gaps_seconds(control_rows)
    if cgaps:
        add(
            f"  gap between pad-seeing scans: min {min(cgaps):.1f} s, "
            f"median {float(np.median(cgaps)):.1f} s, max {max(cgaps):.1f} s"
        )
    add(
        f"  applying the launch-day rule (median of scans ended by 23:23:45) "
        f"gives a baseline of {control_pre['baseline20']:.2f} K and a peak jump of "
        f"{control_pre['peak_jump']:.2f} K"
    )
    add("")
    write_series(add, control_rows, control_anchor)
    add("")
    add("Signal to noise, from the numbers above")
    add(f"  launch peak jump: {peak_jump:.2f} K")
    add(f"  quiet-day std of the neighbourhood max: {noise['std20']:.2f} K over {noise['n']} scans")
    add(f"  quiet-day largest scan-to-scan rise: {noise['max_step']:.2f} K")
    add(f"  peak jump / quiet-day std: {snr_std:.1f}")
    add(f"  peak jump / quiet-day largest rise: {snr_step:.1f}")
    add(
        f"  pre-liftoff std on the launch day itself, same sector: "
        f"{launch_stats['pre_std20']:.2f} K over {launch_stats['n_pre']} scans"
    )
    # Second-largest jump, so the weak tail is in the same units.
    if len(launch_jumps) >= 2:
        second_jump, second = launch_jumps[1]
        add(
            f"  second-largest launch jump: {second_jump:.2f} K at "
            f"{_seconds_after(second['meta']['start'], LIFTOFF):+.0f} s "
            f"({second_jump / noise['std20']:.1f} quiet-day stds)"
        )
    add(
        f"  The quiet-evening median ({noise['median20']:.2f} K) sits "
        f"{launch_stats['baseline20'] - noise['median20']:.2f} K below the launch-evening baseline "
        f"({launch_stats['baseline20']:.2f} K). That day-to-day offset is larger than the "
        f"{noise['std20']:.2f} K within-evening std."
    )
    add("  These ratios describe this one launch and this one quiet evening.")
    add("  They are not a threshold.")
    add("")
    add("Fire")
    add(f"  product file: {fire_path.name}")
    add("  mask 30: temporally_filtered_good_fire_pixel")
    add(
        f"  selected pixel: {_fmt_latlon(fire['lat'], fire['lon'])}, "
        f"{fire['dist_pad_km']:.0f} km from the pad"
    )
    add(
        "  OpenStreetMap places 32.810 N, 82.759 W in Washington County, Georgia."
    )
    add(
        "  The label is a county, not a street: the band-7 pixel is about 2 km across."
    )
    add(
        f"  product fire radiative power: {fire['frp_mw']:.1f} MW, "
        f"product fire temperature: {fire['fire_temp_k']:.0f} K, "
        f"product fire area: {fire['area_m2']:.0f} m^2"
    )
    add("  The product temperature is the algorithm's sub-pixel fire temperature.")
    add("  The temperatures below are band-7 brightness temperature of the whole pixel.")
    add("  Other mask-30 pixels in that scan, by radiative power:")
    for item in fire_top:
        add(
            f"    {item['frp_mw']:7.1f} MW  {_fmt_latlon(item['lat'], item['lon'])}  "
            f"{item['dist_pad_km']:.0f} km from the pad  "
            f"T_product {item['fire_temp_k']:.0f} K"
        )
    add(f"  CONUS band-7 files listed in hours 21–23: {fire_listed}, missing the point: {fire_miss}")
    add(f"  solar zenith at the fire pixel at 23:23:45 UTC: {sun_fire:.1f} deg")
    add(f"  satellite zenith at the fire pixel: {zenith_fire:.1f} deg")
    if fire_rows:
        add(
            f"  pixel spacing: {fire_rows[0]['pixel_dx_km']:.2f} km east-west, "
            f"{fire_rows[0]['pixel_dy_km']:.2f} km north-south"
        )
        fgaps = gaps_seconds(fire_rows)
        add(
            f"  gap between CONUS scans: min {min(fgaps):.1f} s, "
            f"median {float(np.median(fgaps)):.1f} s, max {max(fgaps):.1f} s"
        )
        add(
            f"  nearest-pixel offset from the product coordinate: "
            f"{fire_rows[0]['near_dist_km']:.2f} km on the first scan"
        )
        add(
            f"  peak (neighbourhood max minus 10–30 km median): "
            f"{fire_anom[fire_peak_i]:.2f} K at {_fmt_time(fire_peak['meta']['start'])}, "
            f"nearest pixel {fire_peak['near_bt']:.2f} K, "
            f"max within 20 km {fire_peak['max20_bt']:.2f} K, "
            f"ring {fire_peak['ring_bt']:.2f} K"
        )
        add(
            f"  scans whose neighbourhood max exceeds the ring by more than 10 K: "
            f"{len(hot_fire)} of {len(fire_rows)}, "
            f"spanning {fire_span_min:.0f} min from first to last such scan"
        )
        add("  The 10 K level only counts how long the anomaly lasts. It is not a threshold.")
        add(
            f"  median anomaly over the three hours: {float(np.nanmedian(fire_anom)):.2f} K, "
            f"minimum {float(np.nanmin(fire_anom)):.2f} K"
        )
        fire_steps = [row["d_max20"] for row in fire_rows if np.isfinite(row["d_max20"])]
        if fire_steps:
            add(
                f"  largest scan-to-scan rise of the neighbourhood max: {max(fire_steps):.2f} K, "
                f"largest fall: {min(fire_steps):.2f} K"
            )
        add(f"  scaled ceiling on these CONUS files: {fire_rows[0]['bt_ceiling']:.2f} K")
    add("")
    write_series(add, fire_rows, LIFTOFF)
    add("")
    add("Limits visible in this comparison")
    add("  The quiet mesoscale sample is 12 days after the launch, because the")
    add("  sectors were tasked elsewhere on the adjacent evenings. The Sun was")
    add(f"  {sun_control - sun_launch:.1f} deg closer to the horizon.")
    add("  The fire is a 5-minute CONUS pixel. The launch is a 1-minute mesoscale")
    add("  pixel. Satellite zenith is similar. The fire pixel is in darkness")
    add(f"  (solar zenith {sun_fire:.1f} deg) while the pad is still in daylight")
    add(f"  (solar zenith {sun_launch:.1f} deg), so the backgrounds are not the same")
    add("  physical mixture. Band 7 at the pad still includes reflected sunlight.")
    add("  A 346 K brightness temperature is the average of one pixel, not the")
    add("  flame temperature. The fire product's own temperature is a different")
    add("  quantity and is printed separately.")
    add("  The 20 km neighbourhood follows a hot pixel that moves. It also")
    add("  follows any other warm pixel inside that radius.")
    add("")
    snr_line = (
        f"Peak launch jump {peak_jump:.1f} K. Quiet-evening std {noise['std20']:.2f} K "
        f"({snr_std:.0f}× that std). Quiet-evening largest rise {noise['max_step']:.2f} K. "
        "NOAA GOES-19 ABI L1b. No threshold applied."
    )
    add(f"Figure: {FIGURE_PATH.relative_to(ROOT)}")
    add(f"  {snr_line}")

    plot_series(launch_rows, control_rows, fire_rows, LIFTOFF, control_anchor, snr_line)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    REPORT_PATH.write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
