# Architecture

This document records how monitoreAR is built now, and the decisions that v0.1 will keep. Decisions marked **made** are already in the Milestone 1 or Milestone 2 code, or are fixed by those runs. Decisions marked **v0.1** are the intended design and are not implemented yet. Items marked **open** still need a choice before the milestone that depends on them.

The goal is a small, readable pipeline: fetch one band, convert it in the open, score a few fixed sites, store the raw numbers, and show them. The physics and the scan pattern set the shape of that pipeline more than the web page does.

## Current system

Milestone 1 is `scripts/m1_starship_flight11.py`. It draws the Flight 11 scene. Milestone 2 is `scripts/m2_signal.py`. It compares that scene with a quiet evening and a fire. Milestone 3 is `scripts/m3_score.py`. It scores Pad 1 with the rule below and writes a candidate log. Milestone 4 is `scripts/m4_poll.py`. It polls two live passes and stores every scan in SQLite. There is still no API.

```
noaa-goes19  ABI-L1b-RadM/.../C07
        │  anonymous S3 list + download (s3fs)
        ▼
   local NetCDF under data/m1/     (gitignored)
        │  xarray, scale and offset applied by the library
        ▼
   radiance  →  brightness temperature using planck_fk1/fk2/bc1/bc2
        │
        ▼
   fixed-grid scan angles  →  latitude/longitude on the GRS80 ellipsoid
        │
        ▼
   crop around Pad 1, per-scan statistics, PNG + text report in outputs/
```

Libraries are numpy, xarray, netCDF4, s3fs, and matplotlib. SQLite comes from the standard library. There is no package layout and no web server. The equations, the Flight 11 numbers, and the poller walkthrough are in [`docs/`](docs/README.md).

The case study is Starship Flight 11, GOES-19, 13 October 2025, 23:18–23:36 UTC, ABI Mode 6, band 7 only. Pad 1 is 25.9962°N, 97.1544°W.

## Decisions already made

### Read Level 1b radiance and convert it here

**Made.** The measurement is ABI L1b band 7 (`ABI-L1b-Rad*`, channel `C07`), not the Level 2 Cloud and Moisture Imagery product. L2 CMIP already contains brightness temperature. Keeping the Planck step in our code makes the calibration visible and ties the coefficients to the same file as the radiance.

The formula is the ABI one:

```
T = (fk2 / ln(fk1 / L + 1) - bc1) / bc2
```

`L` is the scaled radiance. On the Flight 11 files the units are mW m⁻² sr⁻¹ (cm⁻¹)⁻¹, and the coefficients are scalar variables in the NetCDF, not attributes of `Rad`.

`esun` and `kappa0` are present and set to NaN on these emissive-band files, so there is no reflectance factor to apply. The conversion leaves reflected sunlight in the brightness temperature. v0.1 will not invent a solar-removal step that the file does not support.

### Use the public NOAA bucket directly

**Made.** Access is anonymous `s3fs` against `s3://noaa-goes19` in `us-east-1`. Object keys follow:

```
ABI-L1b-Rad{F|C|M}/{year}/{day-of-year}/{hour}/OR_ABI-L1b-RadM1-M6C07_G19_s..._e..._c....nc
```

`s` is the scan start, `e` the scan end, and `c` the product-creation time. Each timestamp is `YYYYDDDHHMMSSt` with the last digit in tenths of a second. The script parses those times and treats them as the instrument clock. A published liftoff time is only a label.

`goes2go` and Satpy are not dependencies. The bucket layout, the scaling, and the projection are the things this project is trying to teach.

Downloaded NetCDF stays in `data/` and is gitignored. Figures and inspection reports in `outputs/` are committed when they document a result.

### Geolocate on the fixed grid, and check it

**Made.** ABI `x` and `y` are scan angles in radians after `scale_factor` and `add_offset` are applied. They are not longitude and latitude. The projection variables (`goes_imager_projection`) give the GRS80 ellipsoid, the perspective-point height, the longitude of the projection origin, and `sweep_angle_axis = x`.

The script implements the forward and inverse mapping from PUG L1b Volume 3, section 5.1.2.8. `perspective_point_height` is the height above the ellipsoid. The satellite's distance from Earth's centre is that height plus the equatorial radius.

Two checks run every time:

- A nadir angle maps to latitude 0 and the projection longitude, and the Pad 1 coordinate round-trips to well under a pixel.
- The sector-center scan angle (`x_image`, `y_image`) matches the file's own `geospatial_lat_lon_extent` center. On the Flight 11 mesoscale file that residual was 0.2 m. That checks our math against NOAA's metadata. It does not measure the satellite's pointing error.

The fixed-grid origin on these files is 75.0°W. The nominal subsatellite longitude is 75.2°W, and the nominal height is 35786 km. Navigation uses the fixed grid, which is the grid the pixels were resampled onto. `yaw_flip_flag` was 0.

### Do not correct parallax in v0.1

**Made.** The latitude and longitude of a pixel are the intersection of the line of sight with the ellipsoid. A hot source above the surface is drawn too far from the subsatellite point.

At Pad 1 the satellite zenith angle is about 39°. The shift is about 8 km on the ground for every 10 km of altitude, toward a bearing of 322° (northwest, away from the GOES-East subpoint). The +70 s hot pixel in Flight 11 was 10 km northwest of the pad, while the pad pixel itself had already returned to the background. Motion and parallax are mixed in that offset. v0.1 stores the hot pixel's ellipsoid location and does not try to invert an altitude.

### Treat "2 km" as the nadir figure

**Made.** The file attribute says `2km at nadir`. At Starbase the measured spacing of neighbouring band-7 pixels was 2.31 km east–west and 2.47 km north–south. Watch areas farther from 75°W get larger pixels. A watch-area config should keep the measured spacing, not assume 2 km.

The brightness temperature is the Planck temperature of the average radiance in that pixel. Flight 11 peaked at 346 K. The scaled integer ceiling on the same files is 411 K, and `saturated_pixel_count` was 0. The engines are hotter and much smaller than the pixel. A detector that waits for saturation will miss a launch this size.

### Require the sector to contain the site

**Made**, and it governs the poller. GOES-East in Mode 6 produces a full disk about every 10 minutes, a CONUS sector about every 5 minutes, and two mesoscale sectors about every 60 seconds. Each mesoscale scan lasts about 7 seconds. There is no per-pixel time, so an observation is only known to fall inside that sweep.

During Flight 11 only sector M2 contained Pad 1. M1 was over the tropical Atlantic, centered near 14.4°N, 41°W. The revisit at the pad was 60 seconds, not 30. The pad pixel was hot on the scan starting 23:23:55 UTC and had cooled by the scan starting 23:24:55 UTC. CONUS band-7 files in that hour started at 23:21, 23:26, and 23:31. Full-disk files started at 23:20 and 23:30. Those products were listed and not opened. Their timestamps fall outside the two scans that carried the ~40 K anomaly.

Mesoscale sector boundaries move. M2 is not "the Texas sector" in general. Every file has to be tested by projecting the watch point into `x`/`y` and checking that it lies inside the scan-angle range. The file's lat/lon bounding box is a useful hint and is not the test, because the sector is a rectangle in scan angle.

**Made in Milestone 4.** The poller uses a mesoscale file only when that test passes. If neither sector contains the watch point, it stores `gap_reason=sector_not_covering` and does not score CONUS or full-disk data. A 5-minute product cannot see a signal that occupies one or two 1-minute scans. On 23 November 2026, hours 21 and 22 UTC, every listed mesoscale file missed Pad 1. The newest covering scan in that window was none.

### Keep every quality flag

**Made.** L1b `DQF` meanings on these files are:

| Value | Meaning |
| --- | --- |
| 0 | good pixel |
| 1 | conditionally usable |
| 2 | out of range |
| 3 | no value |
| 4 | focal-plane temperature threshold exceeded |

The Flight 11 crop was entirely flag 0. The detector will still read `flag_meanings` from the file and record the counts. Pixels with flags 2, 3, or 4 do not enter a score. Flag 1 can enter and must be counted separately. A candidate never drops its flags in order to look cleaner.

### Compare scans, because the background is already hot

**v0.1**, forced by the case study. Before liftoff the pad pixel was 299 K and the 98th percentile of the crop was about 300 K. The launch pixel reached 339 K and then 346 K. A threshold on absolute temperature, chosen anywhere near the background, either fires on ordinary warm land or sits so high that it clips the real pulse. The useful feature is the jump from the preceding scans of the same place.

Milestone 2 measured that tail. The table is [`outputs/m2_signal.txt`](outputs/m2_signal.txt) and the walkthrough is [`docs/signal-comparison.md`](docs/signal-comparison.md). Inside 20 km the only jumps above 2 K are the two launch scans: 39.29 K at +9.9 s (the pad pixel, 339.43 K) and 45.41 K at +69.9 s (345.55 K, 10.1 km from the pad), against a pre-liftoff median of 300.14 K. On 25 October, the nearest evening when a mesoscale sector contained the pad at the same clock time, that 20 km maximum has a standard deviation of 0.20 K. The two evenings differ by 1.74 K in the median, which is larger than the within-evening scatter, so an absolute temperature still does not travel from one clear evening to the next.

Farther out, the same window still sees 308.26 K at +130 s (25.4 km) and 307.54 K at +190 s (57.8 km). Those are about 8 K above the 20 km baseline, about 40 quiet-evening standard deviations, and they sit outside the neighbourhood that holds the launch. The 98th percentile within 20 km only reaches 305.02 K on the hot scans.

A CONUS fire pixel in Washington County, Georgia, the same evening stayed about 43 K above its 10–30 km ring for every scan from 21:01 to 23:56 UTC. Its peak band-7 brightness temperature was 339.91 K. The amplitude matches the launch. The duration does not: two mesoscale scans, against at least 175 minutes.

**Made in Milestone 3.** The score is the hottest usable pixel within 20 km minus the median of that maximum on the previous covering scans, up to 5 and at least 3. The threshold is 20 K. DQF values 2, 3, and 4 do not enter. The script does not use the liftoff clock. On Flight 11 it emits 39.29 K at 23:23:54.9 UTC and 45.37 K at 23:24:54.9 UTC. The quiet evening's largest score is 0.04 K. The same rule on the already-burning fire peaks at 12.36 K. The walkthrough is [`docs/scoring.md`](docs/scoring.md). The 10–30 km ring stays in the log as a raw value. It is not the threshold.

Day and night cannot share one absolute cutoff. Solar zenith at the Flight 11 liftoff was about 82°. Band 7 still sees reflected sunlight. A night launch will have a colder, cleaner background. The score stores the local solar zenith angle next to the brightness temperatures so a later threshold can use it. v0.1 does not build a separate reflectance product.

### A launch is short, a fire is not

**Made in Milestone 3**, with one change forced by the Flight 11 fire files. The log records every scan at or above 20 K. A candidate is held back when either:

- the watch site has already been at or above the threshold for 3 consecutive covering scans (a gap longer than 90 s starts the count over), or
- the mesoscale fire product `ABI-L2-FDCM`, for the same sector and the same scan start, has a temporally filtered fire (mask codes 30–35) within 20 km of the hot pixel.

Codes 10–15 are stored and do not hold the scan back. The FDCM file for 23:24:54.9 UTC marks a `good_fire_pixel` (code 10, 17.7 MW) 3.39 km from the Flight 11 hot pixel. Holding every fire code back would have dropped that candidate. The pixel nearest the pad is mask 150, the file's sea-water / bright-desert class, so a plume that stays on the pad is invisible to this product. A missing fire file is stored as `fire_file=unavailable`. The candidate is kept.

### Fixed watch list

**v0.1.** Sites live in a small table, not in a user drawing on the map. Each entry has a name, a latitude, a longitude, and a short note on why GOES-East can see it. The first entry is Starbase Pad 1, because Milestone 1 geolocated it and saw a real pulse there. Adding the Cape or Kourou means repeating the case-study check, including the local pixel size and the satellite zenith angle, before the site is armed.

Coordinates are the pad. The stored detection coordinate is the hot pixel, which can leave the pad within one scan.

### Store candidates in SQLite

**Made in Milestone 4.** The file is `data/m4/monitorear.sqlite`, written by `scripts/m4_poll.py` with the standard-library `sqlite3`. It is gitignored. The walkthrough is [`docs/polling.md`](docs/polling.md), and the latest run is [`outputs/m4_poll.txt`](outputs/m4_poll.txt).

`scans` has one row per source object key. Columns include scan start, scan end, watch-area name, pixel latitude, pixel longitude, distance from the pad, score, the pad-pixel brightness temperature, the neighbourhood maximum, the baseline and the samples that produced it, solar zenith, DQF counts, the saturated flag, the fire-product result, and the source key. Covering scans under 20 K are stored with status `below`, and scans with fewer than three previous looks with status `unscored`. The next baseline is the median of those stored maxima. A candidate is a row with status `emit` or `held_back`. The same key is inserted once. A second pass does not replace the row.

A sector that misses the pad is a `scans` row with `covered=0` and `gap_reason=sector_not_covering`. The `gaps` table is for the other two holes: `file_missing` when a same-sector step is longer than 90 s, and `bucket_lag` when the newest listed scan start is more than 5 minutes behind the poll clock. `poll_runs` stores the clock and the counts for each pass. The API, when it exists, reads these rows and does not recompute from NetCDF. The map can show a gap. It does not draw a launch for one.

### Poll the latest file, and tolerate holes

**Made in Milestone 4.** `scripts/m4_poll.py` lists the current hour and the previous hour of `ABI-L1b-RadM` band 7, scores a file when Pad 1 is inside its scan-angle grid, and exits after two passes 60 s apart. The bucket constant is `noaa-goes19`. `read_watch` stops the process when `platform_ID` is not `G19`. That check is a sanity check. `orbital_slot` is printed by the Milestone 1 inspection and is not a poller column. When NOAA moves GOES-East, the bucket constant changes in one place. The script does not discover the handoff.

`s3fs` caches directory listings, and `fs.ls` has no refresh argument. Each pass calls `invalidate_cache` on the hour prefix before listing. The first run skipped that call. Its second pass repeated 188 keys and a newest scan start of 22:33:54.9, and the report age of 2.4 min was that cache. After the fix, the confirming run at 22:40:02.0 UTC saw newest start 22:39:24.9 (age 37 s) and, one minute later, newest start 22:40:54.9 (age 13 s, about 6 s after that scan ended). Those ages are poll clock minus scan start, for this evening only.

A previous covering scan older than 30 minutes does not enter the baseline, so the Flight 11 rows cannot become the background for a later day. A file that arrives after a later scan was stored is scored from the earlier rows only. The later row is not rewritten. This sample had no late file. Missing keys are logged as `file_missing` once the same-sector step exceeds 90 s. The loop does not interpolate a brightness temperature across a gap. The `c` timestamp in the filename measured 3–5 seconds after scan end on the Flight 11 files. That is product-creation time at the ground station. It is not the time the object became visible on AWS. The 5-minute lag note was not reached on the fresh listings.

### Public surface

**v0.1.** FastAPI serves `GET /detections?since=...` from SQLite. The response is the stored columns, including the raw temperatures and the source key. One HTML page uses Leaflet to plot those points and links each point at its thermal crop.

The crop is written when the candidate is scored, as a small PNG next to the database or as bytes in a column. The page does not re-download the NetCDF. The page states the resolution, the cadence, the miss on small vehicles, and the NOAA attribution. The Cafecito control is a link with a placeholder URL until the real one exists.

There is no account system, no write API, and no alert subscription in v0.1.

### How the code will grow

**v0.1**, at packaging time. Until then, each milestone stays a script that can be run and read on its own. Milestone 6 extracts the shared pieces. The score and the poller are already settled:

| Module | Responsibility |
| --- | --- |
| `goes_io` | list and download one product prefix |
| `calibrate` | scaled radiance to brightness temperature |
| `geolocate` | fixed-grid transforms and the coverage test |
| `score` | window statistics and the conservative rule |
| `store` | SQLite writes |
| `poll` | latest-file loop |
| `api` | FastAPI and the map page |

Notebooks for Milestones 1–3 import those functions and keep the narrative. `pyproject.toml` replaces the loose `requirements.txt` at that point. Dependencies stay on the current list plus FastAPI and a SQLite driver from the standard library. Cartopy, Satpy, and goes2go stay out unless a later problem cannot be stated without them.

## Target runtime

```
cron or a one-minute loop
        │
        ▼
list ABI-L1b-RadM C07 for GOES-East
        │
        ├─ no new file ──────────────► log, sleep
        ├─ sector misses the pad ───► coverage-gap row, sleep
        ▼
download → BT → DQF mask → crop each watch point
        │
        ▼
score against the recent local background
        │
        ├─ below threshold ─────────► row stored (status below); it feeds the next baseline
        ▼
apply persistence / fire-product filter
        │
        ▼
SQLite row + thermal-crop PNG
        │
        ▼
FastAPI reads rows    Leaflet reads the JSON
```

Milestone 4 implements the list, the score, and the SQLite rows. The script takes two passes and exits. The thermal-crop file and the FastAPI page are Milestone 5.

## Non-goals

- A model trained on launch chips.
- Combining GOES-West, Himawari, or Meteosat in v0.1. The view geometry and the buckets differ, and the Flight 11 lesson is already large enough.
- Solving for plume height from the parallax shift.
- Tracking the vehicle across the Gulf as a trajectory product. The hot pixel leaves the pad pixel immediately, and later excesses of a few kelvin are not a track.
- Alerting faster than the scan. The instrument itself waits about a minute, and the public bucket adds more.

## Open points

- Whether a later launch under cloud, or a fire that flares by more than 20 K above its own recent median, still fits this rule. The frozen numbers are in [`docs/scoring.md`](docs/scoring.md). They come from one clear evening.
- The second and third watch areas, after the same geolocation and pixel-size check used for Pad 1.
- The license file, and the real Cafecito URL.
- Whether the map crop is a PNG on disk or a blob in SQLite. Either is fine at v0.1 size. Pick one in Milestone 5 and keep the API field stable.
