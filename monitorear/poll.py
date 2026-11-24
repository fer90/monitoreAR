#!/usr/bin/env python3
"""Milestone 4: poll the latest mesoscale band-7 files and store the score.

The score is the Milestone 3 rule, unchanged. One process lists the current
and previous hour of ABI-L1b-RadM on noaa-goes19, scores a file only when
Pad 1 lies inside its scan-angle grid, and writes SQLite rows. A sector that
misses the pad, a missing key, and a late file are stored as gaps. The same
object key is never scored twice. Brightness temperature is not interpolated.
A new candidate also writes its map crop under outputs/crops/.

This run replays the Flight 11 hour, then takes two live passes a minute
apart. Each pass drops the cached directory listing before it lists. The
second pass must not insert the keys the first pass already stored, and it
must see files that appeared during the sleep.

Run from the repository root:

    .venv/bin/python scripts/m4_poll.py
"""

from __future__ import annotations

from monitorear import flight11 as m1
from monitorear import score as m3
from monitorear import signal as m2
from monitorear.score import m3_baseline
from monitorear.store import (
    connect,
    find_holes,
    fmt_stamp,
    insert_measured,
    insert_outside,
    note_hole,
    note_lag,
    now_stamp,
    parse_stamp,
    previous_covering,
    previous_samples,
    resolve_holes,
    seen_keys,
    self_test_store,
)

import json
import time
from datetime import datetime, timedelta, timezone

import numpy as np
import s3fs


ROOT = m1.ROOT
BUCKET = m1.BUCKET
WATCH = m3.WATCH[0]

POLL_INTERVAL_S = 60
LIVE_CYCLES = 2

DB_PATH = ROOT / "data" / "m4" / "monitorear.sqlite"
CACHE_DIR = ROOT / "data" / "m4" / "cache"
REPORT_PATH = m3.OUT_DIR / "m4_poll.txt"
FLIGHT11_PREFIX = f"{BUCKET}/ABI-L1b-RadM/2025/286/23/"



def locate(fs, key):
    name = key.split("/")[-1]
    for folder in (CACHE_DIR, m2.LAUNCH_DIR, m2.CONTROL_DIR, m3.FIRE_DIR):
        path = folder / name
        if path.exists() and path.stat().st_size > 0:
            return path
    return m2.download(fs, key, CACHE_DIR)




def score_row(conn, fs, cache, source_key, meta, measured):
    samples = previous_samples(conn, WATCH["name"], meta["start"])
    baseline, count, window = m3_baseline(samples)
    finite = np.isfinite(measured["max_bt"])
    if baseline is None or not finite:
        score = None
        scored = False
    else:
        score = float(measured["max_bt"] - baseline)
        scored = True
    prev = previous_covering(conn, WATCH["name"], meta["start"])
    if scored and score >= m3.THRESHOLD_K:
        if (
            prev is not None
            and prev[1] is not None
            and prev[1] >= m3.THRESHOLD_K
            and (meta["start"] - prev[0]).total_seconds() <= m3.GAP_BREAK_S
        ):
            consecutive = prev[2] + 1
        else:
            consecutive = 1
    else:
        consecutive = 0
    measured["score_k"] = score
    measured["scored"] = scored
    measured["consecutive"] = consecutive
    temporal = False
    measured["fire_file"] = None
    measured["fire_nearest_mask"] = None
    measured["fire_nearest_meaning"] = None
    measured["fire_hits"] = []
    if scored and score >= m3.THRESHOLD_K:
        temporal = m3.annotate_fire(fs, measured, cache)
    status, reasons = m3.classify(score if scored else None, consecutive, temporal)
    hits = []
    for item in measured["fire_hits"]:
        hits.append(
            {
                "code": item["code"],
                "meaning": item["meaning"],
                "dist_km": m3.rk(item["dist_km"]),
                "lat": m3.rk(item["lat"], 5),
                "lon": m3.rk(item["lon"], 5),
                "frp_mw": m3.rk(item["frp_mw"]),
            }
        )
    inserted = insert_measured(
        conn,
        source_key=source_key,
        scan_start=meta["start"],
        scan_end=meta["end"],
        sector=meta["sector"],
        max_bt=float(measured["max_bt"]) if finite else None,
        score=score,
        status=status,
        consecutive=consecutive,
        extra={
            "pixel_lat": m3.rk(measured["max_lat"], 5),
            "pixel_lon": m3.rk(measured["max_lon"], 5),
            "distance_km": m3.rk(measured["max_dist_km"]),
            "baseline_k": None if baseline is None else round(baseline, 2),
            "baseline_n": count,
            "baseline_samples_k": json.dumps([round(value, 2) for value in window]),
            "near_bt_k": m3.rk(measured["near_bt"]),
            "near_dqf": measured["near_dqf"],
            "p98_bt_k": m3.rk(measured["p98_bt"]),
            "ring_bt_k": m3.rk(measured["ring_bt"]),
            "solar_zenith_deg": m3.rk(measured["solar_zenith_deg"], 1),
            "dqf_counts": json.dumps(measured["dqf_counts"]),
            "hot_pixel_dqf": measured["max_dqf"],
            "n_conditional_dqf": measured["n_conditional"],
            "n_neighbour_pixels": measured["n_near"],
            "saturated": 1 if measured["saturated"] else 0,
            "bt_ceiling_k": m3.rk(measured["bt_ceiling"]),
            "hold_reasons": json.dumps(reasons),
            "fire_file": measured["fire_file"],
            "fire_nearest_mask": measured["fire_nearest_mask"],
            "fire_nearest_meaning": measured["fire_nearest_meaning"],
            "fire_hits_within_20km": json.dumps(hits),
        },
    )
    return inserted, status


def ingest_listing(fs, conn, cache, listed, label):
    """Score listed (key, meta) pairs that are not already stored."""
    known = seen_keys(conn)
    fresh = [(key, meta) for key, meta in listed if key not in known]
    fresh.sort(key=lambda item: item[1]["start"])
    print(f"{label}: {len(listed)} listed, {len(listed) - len(fresh)} already stored, {len(fresh)} new")
    scored = 0
    outside = 0
    candidates = 0
    by_sector = {}
    for key, meta in listed:
        by_sector.setdefault(meta["sector"], []).append(meta["start"])
    for sector, starts in by_sector.items():
        starts.sort()
        for prev, nxt in find_holes(starts):
            note_hole(conn, sector, prev, nxt)
    for key, meta in fresh:
        try:
            path = locate(fs, key)
            measured = m3.read_watch(path, meta, WATCH["lat"], WATCH["lon"])
        except Exception as exc:
            print(f"  skip {key.split('/')[-1]}: {exc}")
            continue
        if measured is None:
            if insert_outside(conn, key, meta, "sector_not_covering"):
                outside += 1
            resolve_holes(conn, meta["sector"], meta["start"], key)
            continue
        measured["meta"] = meta
        measured["source_key"] = key
        inserted, status = score_row(conn, fs, cache, key, meta, measured)
        resolve_holes(conn, meta["sector"], meta["start"], key)
        if not inserted:
            continue
        scored += 1
        if status in {"emit", "held_back"}:
            candidates += 1
            # The map reads this PNG. It does not open the NetCDF.
            from monitorear import api as m5_api
            m5_api.write_crop(path, key, measured["score_k"], status)
        print(
            f"  {m3.fmt_time(meta['start'])} {meta['sector']} "
            f"{status} score={measured['score_k'] if measured['scored'] else '-'}"
        )
    conn.commit()
    return {
        "listed": len(listed),
        "already_stored": len(listed) - len(fresh),
        "newly_scored": scored,
        "newly_outside": outside,
        "candidates": candidates,
    }


def list_hour(fs, hour_start):
    hour_start = hour_start.astimezone(timezone.utc).replace(minute=0, second=0, microsecond=0)
    doy = int(hour_start.strftime("%j"))
    prefix = f"{BUCKET}/ABI-L1b-RadM/{hour_start.year}/{doy:03d}/{hour_start.hour:02d}/"
    t1 = hour_start + timedelta(hours=1) - timedelta(microseconds=1)
    try:
        # s3fs caches directory listings on this filesystem. A second pass
        # would otherwise repeat the first listing and miss files that
        # arrived during the sleep. fs.ls has no refresh argument.
        fs.invalidate_cache(prefix)
        listed = m2.list_c07(fs, prefix, hour_start, t1)
    except Exception as exc:
        print(f"list failed {prefix}: {type(exc).__name__}: {exc}")
        return prefix, None
    return prefix, listed


def hour_bounds(when):
    start = when.astimezone(timezone.utc).replace(minute=0, second=0, microsecond=0)
    return start, start + timedelta(hours=1) - timedelta(microseconds=1)


def replay_flight11(fs, conn, cache):
    t0, t1 = hour_bounds(m1.LIFTOFF)
    prefix, listed = list_hour(fs, t0)
    if listed is None:
        raise SystemExit(f"could not list {prefix}")
    # The case-study window is the whole hour. list_hour already returns it.
    del t1
    stats = ingest_listing(fs, conn, cache, listed, "flight11 replay")
    rows = conn.execute(
        """
        SELECT scan_start, score_k, status FROM scans
        WHERE watch = ? AND status IN ('emit', 'held_back') AND scan_start LIKE '2025-10-13%'
        ORDER BY scan_start
        """,
        (WATCH["name"],),
    ).fetchall()
    got = [(row["scan_start"], round(row["score_k"], 2), row["status"]) for row in rows]
    expected = [
        ("2025-10-13 23:23:54.9", 39.29, "emit"),
        ("2025-10-13 23:24:54.9", 45.37, "emit"),
    ]
    if got != expected:
        raise SystemExit(f"Flight 11 poller acceptance failed: {got}")
    outside = conn.execute(
        """
        SELECT COUNT(*) AS n FROM scans
        WHERE gap_reason = 'sector_not_covering' AND scan_start LIKE '2025-10-13 23:%'
        """
    ).fetchone()["n"]
    print(f"acceptance: Flight 11 emits the two known scans, {outside} files miss the pad")
    return stats, outside


def poll_once(fs, conn, cache, when):
    started = fmt_stamp(when)
    hours = []
    listed = []
    list_note = []
    for shift in (timedelta(hours=-1), timedelta(0)):
        prefix, hour_listed = list_hour(fs, when + shift)
        hours.append(prefix)
        if hour_listed is None:
            list_note.append(prefix)
            continue
        listed.extend(hour_listed)
        if not hour_listed:
            hour_start = (when + shift).astimezone(timezone.utc).replace(minute=0, second=0, microsecond=0)
            note_hole(conn, "", hour_start, hour_start + timedelta(hours=1))
    stats = ingest_listing(fs, conn, cache, listed, f"live {started}")
    newest = max((meta["start"] for _key, meta in listed), default=None)
    lag = None if newest is None else (when - newest).total_seconds()
    note_lag(conn, newest, -1 if lag is None else lag)
    covering = conn.execute(
        """
        SELECT MAX(scan_start) AS start_at FROM scans
        WHERE watch = ? AND covered = 1 AND scan_start >= ?
        """,
        (WATCH["name"], fmt_stamp(when - timedelta(hours=2))),
    ).fetchone()["start_at"]
    coverage_age = None
    if covering:
        coverage_age = (when - parse_stamp(covering)).total_seconds()
    conn.execute(
        """
        INSERT INTO poll_runs (
            started_at, hours, listed, already_stored, newly_scored, newly_outside,
            candidates, newest_scan_start, lag_seconds, newest_covering_start,
            coverage_age_seconds
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """,
        (
            started,
            " ".join(hours),
            stats["listed"],
            stats["already_stored"],
            stats["newly_scored"],
            stats["newly_outside"],
            stats["candidates"],
            None if newest is None else fmt_stamp(newest),
            lag,
            covering,
            coverage_age,
        ),
    )
    conn.commit()
    stats["started_at"] = started
    stats["hours"] = hours
    stats["list_failures"] = list_note
    stats["newest_scan_start"] = None if newest is None else fmt_stamp(newest)
    stats["lag_seconds"] = lag
    stats["newest_covering_start"] = covering
    stats["coverage_age_seconds"] = coverage_age
    return stats


def write_report(conn, replay_stats, replay_outside, cycles):
    lines = []
    add = lines.append
    add("Milestone 4 — poller, Starbase Pad 1, GOES-19 ABI band 7")
    add("")
    add("What this run did")
    add("  The score, the 20 K threshold, the 20 km neighbourhood, the five-scan")
    add("  median, and the fire-product hold are the Milestone 3 rule. This script")
    add("  does not recompute the threshold.")
    add(f"  Database: {DB_PATH.relative_to(ROOT)}")
    add("  A scan older than 30 min does not enter the baseline, so the Flight 11")
    add("  rows do not become the background for a later day.")
    add("  The same source key is inserted once. A later pass leaves the row as it was.")
    add("")
    add("Flight 11 hour, replayed through the poller")
    add(
        f"  listed {replay_stats['listed']}, newly scored {replay_stats['newly_scored']},"
        f" sector missed the pad {replay_outside}"
    )
    add("  candidates:")
    for row in conn.execute(
        """
        SELECT scan_start, score_k, status, distance_km, max_bt_k, source_key, hold_reasons
        FROM scans
        WHERE status IN ('emit', 'held_back') AND scan_start LIKE '2025-10-13%'
        ORDER BY scan_start
        """
    ):
        add(
            f"    {row['scan_start']}  {row['status']}  score {row['score_k']:.2f} K"
            f"  max {row['max_bt_k']:.2f} K  {row['distance_km']:.2f} km"
        )
        add(f"      {row['source_key']}")
    holes = conn.execute(
        """
        SELECT sector, gap_start, gap_end, resolved_by_key, detail FROM gaps
        WHERE reason = 'file_missing' AND gap_start LIKE '2025-10-13%'
        ORDER BY gap_start
        """
    ).fetchall()
    add(f"  missing-file gaps inside that hour: {len(holes)}")
    add("")
    add("Live passes")
    for cycle in cycles:
        lag = cycle["lag_seconds"]
        lag_s = "none listed" if lag is None else f"{lag / 60:.1f} min"
        age = cycle["coverage_age_seconds"]
        age_s = "no covering scan in the stored window" if age is None else f"{age / 60:.1f} min"
        add(f"  poll at {cycle['started_at']}")
        add(f"    hours: {' '.join(cycle['hours'])}")
        add(
            f"    listed {cycle['listed']}, already stored {cycle['already_stored']},"
            f" newly scored {cycle['newly_scored']}, sector missed the pad {cycle['newly_outside']},"
            f" new candidates {cycle['candidates']}"
        )
        add(f"    newest listed scan start: {cycle['newest_scan_start']}")
        add(f"    age of that scan at the poll clock: {lag_s}")
        add(f"    newest scan that contained the pad: {cycle['newest_covering_start']}")
        add(f"    age of pad coverage: {age_s}")
        if cycle["list_failures"]:
            add(f"    list failures: {', '.join(cycle['list_failures'])}")
    add("")
    live_candidates = conn.execute(
        """
        SELECT scan_start, status, score_k, source_key FROM scans
        WHERE status IN ('emit', 'held_back') AND scan_start NOT LIKE '2025-10-13%'
        ORDER BY scan_start
        """
    ).fetchall()
    add(f"Live candidates stored: {len(live_candidates)}")
    for row in live_candidates:
        add(f"  {row['scan_start']}  {row['status']}  {row['score_k']:.2f} K  {row['source_key']}")
    add("")
    live_holes = conn.execute(
        """
        SELECT reason, sector, gap_start, gap_end, resolved_by_key, detail FROM gaps
        WHERE gap_start NOT LIKE '2025-10-13%'
        ORDER BY noted_at
        """
    ).fetchall()
    add(f"Live gap rows: {len(live_holes)}")
    for row in live_holes[:12]:
        resolved = "open" if row["resolved_by_key"] is None else f"resolved {row['resolved_by_key']}"
        add(
            f"  {row['reason']} {row['sector'] or '-'} {row['gap_start']} .. {row['gap_end']}"
            f"  {resolved}"
        )
        if row["detail"]:
            add(f"    {row['detail']}")
    if len(live_holes) > 12:
        add(f"  ... {len(live_holes) - 12} more rows in the gaps table")
    add("")
    add("Limits visible in this run")
    add("  The threshold is still the 20 K constant chosen from Flight 11.")
    add("  A file that arrives after a later scan was already stored is scored")
    add("  from the earlier scans only. The later row is not rewritten.")
    add("  The c timestamp in a filename is product-creation time, not the time")
    add("  the object appeared in the bucket. The lag above is one poll clock")
    add("  minus the newest scan start in the listing, for this run only.")
    add("  There is no API and no map yet.")
    add("")
    add(f"Report: {REPORT_PATH.relative_to(ROOT)}")
    text = "\n".join(lines) + "\n"
    REPORT_PATH.parent.mkdir(parents=True, exist_ok=True)
    REPORT_PATH.write_text(text)
    return text


def main():
    self_test_store()
    m1.self_test_geolocation()
    fs = s3fs.S3FileSystem(anon=True, client_kwargs={"region_name": "us-east-1"})
    conn = connect(DB_PATH)
    cache = {}
    replay_stats, replay_outside = replay_flight11(fs, conn, cache)
    cycles = []
    for index in range(LIVE_CYCLES):
        if index:
            print(f"sleep {POLL_INTERVAL_S} s")
            time.sleep(POLL_INTERVAL_S)
        when = datetime.now(timezone.utc)
        cycles.append(poll_once(fs, conn, cache, when))
    text = write_report(conn, replay_stats, replay_outside, cycles)
    print(text)
    conn.close()


if __name__ == "__main__":
    main()
