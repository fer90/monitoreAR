# How the page reads a candidate

`scripts/m5_api.py` serves the rows the poller already stored. It does not recompute the 20 K score, and the page does not download NetCDF. Run it from the repository root, after `scripts/m4_poll.py` has created the database:

```bash
.venv/bin/python scripts/m5_api.py
```

The process listens on [http://127.0.0.1:8000/](http://127.0.0.1:8000/) until interrupted. The startup check is `outputs/m5_api.txt`.

## The JSON

`GET /detections` returns candidates only: rows whose status is `emit` or `held_back`, newest scan start first. `since` filters `scan_start`. The stored text is `YYYY-MM-DD HH:MM:SS.t`. A `T` or a trailing `Z` is accepted and compared the same way.

```text
GET /detections?since=2025-10-13T23:24:00Z
```

That call returns the 23:24:54.9 scan and leaves out 23:23:54.9. With no `since`, both Flight 11 emits are in the list:

| Scan start | Score | Maximum | Distance | Fire product |
| --- | --- | --- | --- | --- |
| 2025-10-13 23:24:54.9 | 45.37 K | 345.55 K | 10.09 km | code 10, 3.39 km, 17.74 MW. Nearest mask 150 |
| 2025-10-13 23:23:54.9 | 39.29 K | 339.43 K | 1.24 km | no code 10–35 within 20 km. Nearest mask 150 |

Temperatures in the JSON are rounded to 0.01 K, the same precision as the reports. The response also carries the baseline samples, the pad-pixel temperature, the solar zenith, the DQF counts, the saturated flag, the source object key, and `crop_url`.

`GET /coverage` is the newest stored scan, whatever its status, plus any open row in the `gaps` table. On this database the newest scan is 2026-11-23 22:40:54.9, sector M2, `sector_not_covering`. Open gap rows: 0. The map draws the pad for that state and does not draw a launch.

## The crop

The crop is a PNG in `outputs/crops/`, named from the object key. `crop_url` points at `/crops/<that file>`. A new candidate is written by the poller from the file it just scored. This script fills in a PNG that is still missing, using a local NetCDF copy when one is already on disk.

The picture is the Milestone 1 window, geolocated on the fixed grid, with the same colour scale: 283.43 K to 350 K. The open circle is Pad 1. The white cross is the hottest pixel within 20 km when that pixel is more than 10 K above the crop median and more than 3 km from the pad, which is the mark on the Milestone 1 figure. The 23:23:54.9 hot pixel is the pad pixel, so that crop has no cross. The 23:24:54.9 hot pixel is 10.09 km away, and the cross is on it.

## The page

`GET /` is one HTML page. It fetches `/detections` and `/coverage`. The orange points are emitted hot pixels. A held candidate would be gray. The open circle is the pad, and it is not a detection. Each point's popup and each card show the crop. The limits on the page are the pixel size at Pad 1, the 60 s revisit, the coverage gap, the 20 K rule, the miss on small vehicles, and the NOAA attribution.

The basemap is OpenStreetMap, loaded in the browser along with Leaflet. The Cafecito control is a link to `https://cafecito.app/` labelled as a placeholder. The project page is not set up. There is no account, no write route, and no alert.

## What this page does not settle

The two points are the Flight 11 scans the threshold was chosen from. The newest stored scan is a later hour when neither mesoscale sector contained the pad. A fire that ignites in view can still be emitted for two scans before the persistence rule holds the third. The coordinate is the hot pixel on the ellipsoid. The 10 km offset is left as measured. The colour scale is the Flight 11 scale, so a much colder night will sit near the bottom of the same range.
