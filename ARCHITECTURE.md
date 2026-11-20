# Architecture

This document records how monitoreAR is built now, and the decisions that v0.1 will keep. Decisions marked **made** are already in the Milestone 1 or Milestone 2 code, or are fixed by those runs. Decisions marked **v0.1** are the intended design and are not implemented yet. Items marked **open** still need a choice before the milestone that depends on them.

The goal is a small, readable pipeline: fetch one band, convert it in the open, score a few fixed sites, store the raw numbers, and show them. The physics and the scan pattern set the shape of that pipeline more than the web page does.

## Current system

Milestone 1 is `scripts/m1_starship_flight11.py`. It draws the Flight 11 scene. Milestone 2 is `scripts/m2_signal.py`. It compares that scene with a quiet evening and a fire. Neither script detects launches.

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

Libraries are numpy, xarray, netCDF4, s3fs, and matplotlib. There is no package layout, no database, and no web server. The equations and the Flight 11 numbers are written out in [`docs/`](docs/README.md).

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

**v0.1.** The near-real-time loop uses a mesoscale file only when that test passes. If neither sector contains the watch point, the loop writes a coverage gap and does not score CONUS or full-disk data. A 5-minute product cannot see a signal that occupies one or two 1-minute scans. Falling back would create quiet false confidence.

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

A CONUS fire pixel in Washington County, Georgia, the same evening stayed about 43 K above its 10–30 km ring for every scan from 21:01 to 23:56 UTC. Its peak band-7 brightness temperature was 339.91 K. The amplitude matches the launch. The duration does not: two mesoscale scans, against at least 175 minutes. Milestone 3 still freezes the threshold, the radius, and the fire product. The 20 km window, the 10–30 km ring, and the 10 K duration count are measurements of this case.

Day and night cannot share one absolute cutoff. Solar zenith at the Flight 11 liftoff was about 82°. Band 7 still sees reflected sunlight. A night launch will have a colder, cleaner background. The score stores the local solar zenith angle next to the brightness temperatures so a later threshold can use it. v0.1 does not build a separate reflectance product.

### A launch is short, a fire is not

**v0.1.** The obvious Flight 11 anomaly lasted two mesoscale scans. A fire or an industrial source tends to stay hot. The log records the raw score either way. A candidate is held back when either:

- the same pixel stays anomalous for several consecutive scans, or
- the official GOES fire product flags that pixel.

The fire product is a filter, not the detector. If the fire file is missing, the candidate is stored and marked `fire_product=unavailable` rather than silently passed or dropped. Milestone 2 measured one real case in the CONUS product `ABI-L2-FDCC`, mask 30, on the scan starting 23:21:17 UTC. Which product the filter consults when a mesoscale sector covers the site is **open** until Milestone 3.

### Fixed watch list

**v0.1.** Sites live in a small table, not in a user drawing on the map. Each entry has a name, a latitude, a longitude, and a short note on why GOES-East can see it. The first entry is Starbase Pad 1, because Milestone 1 geolocated it and saw a real pulse there. Adding the Cape or Kourou means repeating the case-study check, including the local pixel size and the satellite zenith angle, before the site is armed.

Coordinates are the pad. The stored detection coordinate is the hot pixel, which can leave the pad within one scan.

### Store candidates in SQLite

**v0.1.** One table, one row per candidate. Columns include scan start, scan end, watch-area name, pixel latitude, pixel longitude, distance from the pad, score, pad-pixel brightness temperature, neighbourhood maximum, background value used in the score, solar zenith, DQF counts, saturated flag, fire-filter result, and the source object key.

SQLite is a single file the poller writes and the API reads. The public API does not recompute from NetCDF. JSON would be enough for a handful of rows; SQLite is the choice because the questions are by time and by site. The schema stays boring and can be dumped to JSON for the API response.

A coverage-gap row is a different table, or the same table with a null score and a reason such as `sector_not_covering` or `file_missing`. The map can show a gap. It does not draw a launch for one.

### Poll the latest file, and tolerate holes

**v0.1.** A single process wakes about once a minute, lists the current and previous hour under `ABI-L1b-RadM` for the operational GOES-East bucket, and considers band-7 objects whose scan start is newer than the last one scored.

The bucket name comes from configuration, default `noaa-goes19`, because that is the operational GOES-East satellite as of 2026. The code still reads `platform_ID` and `orbital_slot` from the file. A hardcoded "G19" check is a sanity check, not the way a future satellite handoff should be discovered. When NOAA moves GOES-East, the bucket constant changes in one place.

Missing keys are normal. The loop logs the hole, waits, and tries again. It does not interpolate a brightness temperature across a gap. The `c` timestamp in the filename measured 3–5 seconds after scan end on the Flight 11 files. That is product-creation time at the ground station. It is not the time the object became visible on AWS. The poller has to accept several extra minutes of bucket lag, and a file that appears out of order, without scoring the same key twice.

### Public surface

**v0.1.** FastAPI serves `GET /detections?since=...` from SQLite. The response is the stored columns, including the raw temperatures and the source key. One HTML page uses Leaflet to plot those points and links each point at its thermal crop.

The crop is written when the candidate is scored, as a small PNG next to the database or as bytes in a column. The page does not re-download the NetCDF. The page states the resolution, the cadence, the miss on small vehicles, and the NOAA attribution. The Cafecito control is a link with a placeholder URL until the real one exists.

There is no account system, no write API, and no alert subscription in v0.1.

### How the code will grow

**v0.1**, at packaging time. Until the detector exists, each milestone stays a script that can be run and read on its own. Milestone 6 extracts the shared pieces only after Milestones 2 and 3 have settled the math:

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
        ├─ below threshold ─────────► nothing stored (the scan time is still remembered)
        ▼
apply persistence / fire-product filter
        │
        ▼
SQLite row + thermal-crop PNG
        │
        ▼
FastAPI reads rows    Leaflet reads the JSON
```

## Non-goals

- A model trained on launch chips.
- Combining GOES-West, Himawari, or Meteosat in v0.1. The view geometry and the buckets differ, and the Flight 11 lesson is already large enough.
- Solving for plume height from the parallax shift.
- Tracking the vehicle across the Gulf as a trajectory product. The hot pixel leaves the pad pixel immediately, and later excesses of a few kelvin are not a track.
- Alerting faster than the scan. The instrument itself waits about a minute, and the public bucket adds more.

## Open points

- The numeric threshold, the background window, and the neighbourhood radius. Milestone 2 measured them and did not freeze them. On Flight 11 the 20 km maximum jumped 45.41 K and, one scan earlier, 39.29 K, against a quiet-evening standard deviation of 0.20 K and a day-to-day offset of 1.74 K. A fire the same evening held a ~43 K anomaly for at least 175 minutes. The 20 km neighbourhood, the 10–30 km ring, and the 10 K duration count are the windows used to see that. Milestone 3 freezes the rule. Details are in [`docs/signal-comparison.md`](docs/signal-comparison.md).
- Which published fire product is available at mesoscale cadence for these sites.
- The second and third watch areas, after the same geolocation and pixel-size check used for Pad 1.
- The license file, and the real Cafecito URL.
- Whether the map crop is a PNG on disk or a blob in SQLite. Either is fine at v0.1 size. Pick one in Milestone 5 and keep the API field stable.
