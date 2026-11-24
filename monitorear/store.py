"""SQLite rows for scans, coverage gaps, and poll passes."""

from __future__ import annotations

import json
import sqlite3
from datetime import datetime, timedelta, timezone

import numpy as np

from monitorear import score as m3
from monitorear.score import m3_baseline

WATCH = m3.WATCH[0]
# A previous scan older than this is not the local background. Flight 11's
# own scans are a minute apart, so the replay is unchanged.
BASELINE_MAX_AGE_S = 30 * 60
# Same-sector starts are about 60 s apart. A wider gap is a missing file.
HOLE_GAP_S = m3.GAP_BREAK_S
# Report bucket lag once the newest listed scan is older than this.
LAG_NOTE_S = 5 * 60

SCHEMA = """
CREATE TABLE IF NOT EXISTS scans (
    source_key TEXT PRIMARY KEY,
    watch TEXT NOT NULL,
    covered INTEGER NOT NULL,
    gap_reason TEXT,
    scan_start TEXT NOT NULL,
    scan_end TEXT NOT NULL,
    sector TEXT NOT NULL,
    pixel_lat REAL,
    pixel_lon REAL,
    distance_km REAL,
    score_k REAL,
    baseline_k REAL,
    baseline_n INTEGER,
    baseline_samples_k TEXT,
    near_bt_k REAL,
    near_dqf INTEGER,
    max_bt_k REAL,
    p98_bt_k REAL,
    ring_bt_k REAL,
    solar_zenith_deg REAL,
    dqf_counts TEXT,
    hot_pixel_dqf INTEGER,
    n_conditional_dqf INTEGER,
    n_neighbour_pixels INTEGER,
    saturated INTEGER,
    bt_ceiling_k REAL,
    consecutive INTEGER,
    status TEXT NOT NULL,
    hold_reasons TEXT,
    fire_file TEXT,
    fire_nearest_mask INTEGER,
    fire_nearest_meaning TEXT,
    fire_hits_within_20km TEXT,
    stored_at TEXT NOT NULL
);
CREATE TABLE IF NOT EXISTS gaps (
    id INTEGER PRIMARY KEY,
    watch TEXT NOT NULL,
    reason TEXT NOT NULL,
    sector TEXT,
    gap_start TEXT NOT NULL,
    gap_end TEXT,
    noted_at TEXT NOT NULL,
    resolved_by_key TEXT,
    detail TEXT,
    UNIQUE (watch, reason, sector, gap_start)
);
CREATE TABLE IF NOT EXISTS poll_runs (
    id INTEGER PRIMARY KEY,
    started_at TEXT NOT NULL,
    hours TEXT NOT NULL,
    listed INTEGER NOT NULL,
    already_stored INTEGER NOT NULL,
    newly_scored INTEGER NOT NULL,
    newly_outside INTEGER NOT NULL,
    candidates INTEGER NOT NULL,
    newest_scan_start TEXT,
    lag_seconds REAL,
    newest_covering_start TEXT,
    coverage_age_seconds REAL
);
CREATE INDEX IF NOT EXISTS scans_watch_start ON scans (watch, scan_start);
"""


def connect(path):
    path.parent.mkdir(parents=True, exist_ok=True)
    conn = sqlite3.connect(path)
    conn.row_factory = sqlite3.Row
    conn.execute("PRAGMA foreign_keys = ON")
    conn.executescript(SCHEMA)
    return conn


def self_test_store():
    """Idempotent keys, a baseline that ignores stale scans, and hole detection."""
    conn = sqlite3.connect(":memory:")
    conn.row_factory = sqlite3.Row
    conn.executescript(SCHEMA)
    watch = WATCH["name"]
    now = datetime(2025, 10, 13, 23, 24, tzinfo=timezone.utc)
    stale = now - timedelta(hours=2)
    insert_measured(
        conn,
        source_key="stale",
        scan_start=stale,
        scan_end=stale,
        sector="M2",
        max_bt=250.0,
        score=None,
        status="unscored",
        consecutive=0,
    )
    samples = previous_samples(conn, watch, now)
    if samples:
        raise SystemExit("a two-hour-old scan entered the baseline")
    for index, value in enumerate((300.14, 300.26, 300.18, 300.06, 300.01)):
        start = datetime(2025, 10, 13, 23, 18, tzinfo=timezone.utc) + timedelta(minutes=index)
        insert_measured(
            conn,
            source_key=f"pre-{index}",
            scan_start=start,
            scan_end=start,
            sector="M2",
            max_bt=value,
            score=0.0,
            status="below",
            consecutive=0,
        )
    samples = previous_samples(conn, watch, datetime(2025, 10, 13, 23, 23, 54, tzinfo=timezone.utc))
    baseline, count, _window = m3_baseline(samples)
    if count != 5 or abs(baseline - 300.14) > 1e-6:
        raise SystemExit(f"poller baseline self-test failed: {baseline} n={count}")
    insert_measured(
        conn,
        source_key="pre-0",
        scan_start=datetime(2025, 10, 13, 23, 18, tzinfo=timezone.utc),
        scan_end=datetime(2025, 10, 13, 23, 18, tzinfo=timezone.utc),
        sector="M2",
        max_bt=999.0,
        score=0.0,
        status="below",
        consecutive=0,
    )
    if conn.execute("SELECT max_bt_k FROM scans WHERE source_key='pre-0'").fetchone()["max_bt_k"] != 300.14:
        raise SystemExit("a stored key was scored twice")
    holes = find_holes(
        [
            datetime(2025, 10, 13, 23, 18, tzinfo=timezone.utc),
            datetime(2025, 10, 13, 23, 19, tzinfo=timezone.utc),
            datetime(2025, 10, 13, 23, 21, tzinfo=timezone.utc),
        ]
    )
    if len(holes) != 1 or (holes[0][1] - holes[0][0]).total_seconds() != 120:
        raise SystemExit(f"hole detection self-test failed: {holes}")
    if find_holes(
        [
            datetime(2025, 10, 13, 23, 18, tzinfo=timezone.utc),
            datetime(2025, 10, 13, 23, 19, tzinfo=timezone.utc),
        ]
    ):
        raise SystemExit("a 60 s step was called a missing file")
    conn.close()
    print("store self-test: idempotent keys, fresh baseline, one missing-file hole")

def find_holes(starts):
    holes = []
    for prev, nxt in zip(starts, starts[1:]):
        gap = (nxt - prev).total_seconds()
        if gap > HOLE_GAP_S:
            holes.append((prev, nxt))
    return holes


def fmt_stamp(when):
    return m3.fmt_time(when)


def parse_stamp(text):
    # YYYY-MM-DD HH:MM:SS.t — the tenth is the last character.
    base, tenth = text[:-2], text[-1]
    when = datetime.strptime(base, "%Y-%m-%d %H:%M:%S").replace(tzinfo=timezone.utc)
    return when + timedelta(milliseconds=int(tenth) * 100)


def now_stamp():
    return fmt_stamp(datetime.now(timezone.utc))

def seen_keys(conn):
    return {row["source_key"] for row in conn.execute("SELECT source_key FROM scans")}


def previous_samples(conn, watch, start):
    earliest = fmt_stamp(start - timedelta(seconds=BASELINE_MAX_AGE_S))
    rows = conn.execute(
        """
        SELECT max_bt_k FROM scans
        WHERE watch = ? AND covered = 1 AND max_bt_k IS NOT NULL
          AND scan_start < ? AND scan_start >= ?
        ORDER BY scan_start DESC
        LIMIT ?
        """,
        (watch, fmt_stamp(start), earliest, m3.BASELINE_SCANS),
    ).fetchall()
    return [float(row["max_bt_k"]) for row in reversed(rows)]


def previous_covering(conn, watch, start):
    earliest = fmt_stamp(start - timedelta(seconds=BASELINE_MAX_AGE_S))
    row = conn.execute(
        """
        SELECT scan_start, score_k, consecutive FROM scans
        WHERE watch = ? AND covered = 1 AND scan_start < ? AND scan_start >= ?
        ORDER BY scan_start DESC
        LIMIT 1
        """,
        (watch, fmt_stamp(start), earliest),
    ).fetchone()
    if row is None:
        return None
    return parse_stamp(row["scan_start"]), row["score_k"], int(row["consecutive"] or 0)


def insert_measured(
    conn,
    source_key,
    scan_start,
    scan_end,
    sector,
    max_bt,
    score,
    status,
    consecutive,
    watch=None,
    gap_reason=None,
    extra=None,
):
    """Insert one measured scan. A duplicate key is left as it was."""
    watch = watch or WATCH["name"]
    extra = extra or {}
    cur = conn.execute(
        """
        INSERT OR IGNORE INTO scans (
            source_key, watch, covered, gap_reason, scan_start, scan_end, sector,
            pixel_lat, pixel_lon, distance_km, score_k, baseline_k, baseline_n,
            baseline_samples_k, near_bt_k, near_dqf, max_bt_k, p98_bt_k, ring_bt_k,
            solar_zenith_deg, dqf_counts, hot_pixel_dqf, n_conditional_dqf,
            n_neighbour_pixels, saturated, bt_ceiling_k, consecutive, status,
            hold_reasons, fire_file, fire_nearest_mask, fire_nearest_meaning,
            fire_hits_within_20km, stored_at
        ) VALUES (
            ?, ?, 1, ?, ?, ?, ?,
            ?, ?, ?, ?, ?, ?,
            ?, ?, ?, ?, ?, ?,
            ?, ?, ?, ?,
            ?, ?, ?, ?, ?,
            ?, ?, ?, ?,
            ?, ?
        )
        """,
        (
            source_key,
            watch,
            gap_reason,
            fmt_stamp(scan_start),
            fmt_stamp(scan_end),
            sector,
            extra.get("pixel_lat"),
            extra.get("pixel_lon"),
            extra.get("distance_km"),
            score,
            extra.get("baseline_k"),
            extra.get("baseline_n"),
            extra.get("baseline_samples_k"),
            extra.get("near_bt_k"),
            extra.get("near_dqf"),
            max_bt,
            extra.get("p98_bt_k"),
            extra.get("ring_bt_k"),
            extra.get("solar_zenith_deg"),
            extra.get("dqf_counts"),
            extra.get("hot_pixel_dqf"),
            extra.get("n_conditional_dqf"),
            extra.get("n_neighbour_pixels"),
            extra.get("saturated"),
            extra.get("bt_ceiling_k"),
            consecutive,
            status,
            extra.get("hold_reasons"),
            extra.get("fire_file"),
            extra.get("fire_nearest_mask"),
            extra.get("fire_nearest_meaning"),
            extra.get("fire_hits_within_20km"),
            now_stamp(),
        ),
    )
    return cur.rowcount == 1


def insert_outside(conn, source_key, meta, reason):
    cur = conn.execute(
        """
        INSERT OR IGNORE INTO scans (
            source_key, watch, covered, gap_reason, scan_start, scan_end, sector,
            status, stored_at
        ) VALUES (?, ?, 0, ?, ?, ?, ?, ?, ?)
        """,
        (
            source_key,
            WATCH["name"],
            reason,
            fmt_stamp(meta["start"]),
            fmt_stamp(meta["end"]),
            meta["sector"],
            "gap",
            now_stamp(),
        ),
    )
    return cur.rowcount == 1


def note_hole(conn, sector, start, end):
    conn.execute(
        """
        INSERT OR IGNORE INTO gaps (
            watch, reason, sector, gap_start, gap_end, noted_at, detail
        ) VALUES (?, 'file_missing', ?, ?, ?, ?, ?)
        """,
        (
            WATCH["name"],
            sector,
            fmt_stamp(start),
            fmt_stamp(end),
            now_stamp(),
            f"{(end - start).total_seconds():.0f} s between {sector} scan starts; no brightness temperature filled in",
        ),
    )


def resolve_holes(conn, sector, start, source_key):
    stamp = fmt_stamp(start)
    conn.execute(
        """
        UPDATE gaps
        SET resolved_by_key = ?
        WHERE watch = ? AND reason = 'file_missing' AND sector = ?
          AND resolved_by_key IS NULL AND gap_start < ? AND gap_end > ?
        """,
        (source_key, WATCH["name"], sector, stamp, stamp),
    )


def note_lag(conn, newest_start, lag_s):
    if newest_start is None or lag_s <= LAG_NOTE_S:
        conn.execute(
            """
            UPDATE gaps SET resolved_by_key = 'lag_closed', gap_end = ?
            WHERE watch = ? AND reason = 'bucket_lag' AND resolved_by_key IS NULL
            """,
            (now_stamp(), WATCH["name"]),
        )
        return
    conn.execute(
        """
        INSERT INTO gaps (watch, reason, sector, gap_start, gap_end, noted_at, detail)
        VALUES (?, 'bucket_lag', ?, ?, ?, ?, ?)
        ON CONFLICT (watch, reason, sector, gap_start) DO UPDATE SET
            gap_end = excluded.gap_end,
            noted_at = excluded.noted_at,
            detail = excluded.detail,
            resolved_by_key = NULL
        """,
        (
            WATCH["name"],
            "",
            fmt_stamp(newest_start),
            now_stamp(),
            now_stamp(),
            f"newest listed scan start is {lag_s / 60:.1f} min behind the poll clock",
        ),
    )
