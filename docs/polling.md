# How the poller stores a scan

`scripts/m4_poll.py` applies the Milestone 3 score to files as they appear, and writes the rows to SQLite. It does not recompute the 20 K threshold. It writes:

- `data/m4/monitorear.sqlite` — scans, coverage gaps, and one row per pass (gitignored with the rest of `data/`)
- `outputs/m4_poll.txt` — the Flight 11 acceptance check and the live passes from the latest run

The only armed site is Starbase Pad 1. Re-run from the repository root:

```bash
.venv/bin/python scripts/m4_poll.py
```

A store self-test runs before any download. It checks that the same object key is inserted once, that a scan two hours old stays out of the baseline, and that a 120 s hole is recorded while a 60 s step is not. The script then replays the Flight 11 hour and takes two live passes, 60 s apart, and exits. `LIVE_CYCLES` is 2. Raising it keeps the same loop going. The script does not run forever.

## One pass

Each pass lists band 7 (`C07`) under `ABI-L1b-RadM` for the current hour and the previous hour on `noaa-goes19`. The prefix is the same layout Milestone 1 uses. Before the list, the poller calls `invalidate_cache` on that hour prefix. `s3fs` keeps directory listings on the filesystem object, and `fs.ls` has no refresh argument. Without the cache drop, the second pass repeats the first listing. That happened on the first run: at 22:36:17.5 UTC the listing was still the 188 files from 22:34:16.7, newest scan start still 22:33:54.9, age 143 s. A fresh list a few seconds later already contained a scan starting 22:35:54.9. The 2.4 min age was the cache.

A key already in `scans` is left as it was, including when a file shows up after a later scan was stored. The late file is scored from the earlier rows only. The later row is not rewritten. This run did not contain a late arrival, so that path is covered by the insert rule and not by a live example.

A file whose scan-angle grid contains the pad is scored with the Milestone 3 rule. The baseline is the median of the neighbourhood maximum on the previous covering scans, up to 5, at least 3, and none of them older than 30 minutes. Flight 11's scans are a minute apart, so the replay matches Milestone 3: 39.29 K and 45.37 K. The 30-minute cutoff is what keeps those rows from becoming the background on a later day. The self-test puts a two-hour-old sample outside the window. The live hours below had no covering scan, so the cutoff was not exercised on 23 November data.

The mesoscale fire product is opened only when the score is at or above 20 K, for the same sector and the same scan start. Codes 30–35 within 20 km of the hot pixel hold the row. Codes 10–15 are stored and do not hold it. A missing fire file stays `fire_file=unavailable`.

A file whose grid misses the pad is stored with `covered=0` and `gap_reason=sector_not_covering`. The poller does not open CONUS or full-disk data for that minute. A download or a read that raises is printed and not stored, so the next pass tries the key again. A file whose `platform_ID` is not `G19` still stops the process. That check is the sanity check in `read_watch`.

## What is stored

`scans` has one row per object key. Covering scans under 20 K are stored with status `below`, and the first looks of a series with status `unscored`, because the next baseline is the median of those neighbourhood maxima. A candidate is a row with status `emit` or `held_back`. There is no separate candidates table.

`gaps` records two other holes. `file_missing` is a same-sector step longer than 90 s, or an hour prefix that listed no band-7 file. `bucket_lag` is written when the newest listed scan start is more than 5 minutes behind the poll clock, and cleared when a later pass is inside that window. A sector that simply is not aimed at the pad is a `scans` row, not a `gaps` row. `poll_runs` stores the clock, the counts, and the age for each pass.

The age in the report is the poll clock minus the newest listed scan start. It includes the scan itself, about 7 s on these files. The `c` token in the filename is product-creation time. It is not the time the object became visible on AWS.

NetCDF is reused from `data/m1/`, `data/m2/`, `data/m3/fire/`, or `data/m4/cache/` when the local file is already non-empty. New downloads go to `data/m4/cache/`.

## The Flight 11 hour

The replay lists all of hour 23 on day 286 of 2025, not only the 18-minute window Milestone 3 plotted. The database holds 120 band-7 files:

| Rows | What they are |
| --- | --- |
| 65 | `sector_not_covering` (the other mesoscale sector) |
| 3 | covering, `unscored`, fewer than 3 previous looks |
| 50 | covering, `below`, score under 20 K |
| 2 | `emit` |

The emitted scan starts are 2025-10-13 23:23:54.9 (39.29 K, 339.43 K, 1.24 km) and 23:24:54.9 (45.37 K, 345.55 K, 10.09 km). The script stops if those two rows are anything else. The confirming run listed the same 120 keys, found all of them already stored, and inserted nothing. Same-sector steps in that hour stayed inside 90 s, so the hour has no `file_missing` row.

## The live passes on 23 November 2026

The confirming run, after the listing cache was dropped, polled twice:

| Poll clock (UTC) | Listed | Already stored | New rows | Newest scan start | Age |
| --- | --- | --- | --- | --- | --- |
| 22:40:02.0 | 199 | 188 | 11 | 22:39:24.9 | 0.6 min (37 s) |
| 22:41:07.5 | 202 | 199 | 3 | 22:40:54.9 | 0.2 min (13 s) |

Both passes listed hours 21 and 22 of day 327. Every new row missed the pad. The newest scan that contained the pad was none. Live candidates: 0. `gaps`: 0 rows. Hour 21 is a full set of 120 files, all `sector_not_covering`. Hour 22 has 82 files through the scan starting 22:40:54.9, also all outside the pad. The sectors were aimed elsewhere for this whole sample. That is a coverage gap.

The second pass's newest file starts 22:40:54.9, ends 22:41:01.7, and carries creation time 22:41:04.9. The poll clock was 22:41:07.5, about 6 s after scan end and about 3 s after the creation timestamp. On the first of these two passes the newest listed file ended at 22:39:31.7. The next file, sector M2, ended at 22:40:02.0, the instant that poll started, so it was not in the listing yet. An earlier cold list, at 22:34:16.7, saw a scan that had ended at 22:34:01.8. Three looks, one evening. None of them crossed the 5-minute lag note.

## What this run does not settle

The threshold is still the constant 20 K chosen from Flight 11. The poller does not fit a new cut to live data. Catching those two scans is a consistency check. Clouds were not in the evenings that set the cut. A fire that ignites in view can still pass for two scans. The live minute can also be a coverage gap: for these two hours it was, and the score does not fall back to a 5-minute product. A file that arrives out of order was not in this sample. There is no second watch site, no API, and no map.
