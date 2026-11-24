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

from monitorear import flight11 as m1
from monitorear import signal as m2
from monitorear.score import (
    NEIGHBOUR_KM,
    PERSIST_SCANS,
    THRESHOLD_K,
    WATCH,
    classify,
    finalize,
    fmt_latlon,
    fmt_time,
    read_watch,
    rk,
    self_test_rule,
)

import json
from datetime import datetime, timezone

import numpy as np
import s3fs


ROOT = m1.ROOT
BUCKET = m1.BUCKET

# Frozen from Milestone 2. The quiet evening's peak jump was 0.21 K.
# Flight 11 inside 20 km jumped 39.29 K and then 45.41 K. The 8 K pixels
# at +130 s and +190 s were 25 km and 58 km from the pad, outside this radius.

OUT_DIR = ROOT / "outputs"
REPORT_PATH = OUT_DIR / "m3_score.txt"
JSONL_PATH = OUT_DIR / "m3_candidates.jsonl"




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
