# How Milestone 3 scores a watch site

`scripts/m3_score.py` turns the Milestone 2 measurements into a rule and runs it. It writes:

- `outputs/m3_score.txt` — every scored scan, and the candidates in full
- `outputs/m3_candidates.jsonl` — one JSON object per candidate, including a candidate that was held back

The only armed site is Starbase Pad 1, 25.9962°N, 97.1544°W. The script does not know the liftoff time. Re-run from the repository root:

```bash
.venv/bin/python scripts/m3_score.py
```

A self-test of the median and the hold rule runs before any download. After the files are scored, the script checks that Flight 11 produces the two known scans and that the quiet evening produces none.

## The score

For each mesoscale band-7 scan whose fixed grid contains the pad:

1. Take the hottest finite pixel within 20 km. Pixels with DQF 2 (out of range), 3 (no value), or 4 (focal-plane temperature) stay out. DQF 1 may enter and is counted.
2. The baseline is the median of that same maximum on the previous covering scans, up to 5, and at least 3. Scans before the third prior look are logged as `unscored`.
3. The score is the current maximum minus that median.

A median of five is there so one hot scan, or two, does not become the baseline. On the second Flight 11 scan the five previous maxima are 300.26, 300.18, 300.06, 300.01, and 339.43 K. The median is 300.18 K. The 339.43 K value is in the window and does not set it.

The 98th percentile is stored and is not the score. On the two hot scans it reaches 305.02 K and 304.92 K, about 5 K above the pre-liftoff value. A 20 K rule on the percentile would have missed the launch. The maximum moves by tens of kelvin because the hot object fills one pixel.

The ring median, 10–30 km from the pad, is stored too. On both hot scans it stays 298.59 K. The score uses the earlier scans of the same neighbourhood, because that is the comparison that removes the day-to-day offset Milestone 2 measured.

## Why 20 K and 20 km

| What Milestone 2 measured | Value |
| --- | --- |
| Quiet evening, fixed pre-liftoff baseline, peak jump | 0.21 K |
| Quiet evening, standard deviation of the 20 km maximum | 0.20 K |
| Flight 11, first hot scan, inside 20 km | 39.29 K |
| Flight 11, second hot scan, inside 20 km | 45.41 K above 300.14 K |
| Pixels at +130 s and +190 s | 308.26 K at 25.4 km, 307.54 K at 57.8 km |

20 km is the radius that contains both launch pixels (1.24 km, then 10.09 km) and leaves the 8 K downrange pixels outside. 20 K sits far above the quiet evening and about halfway under the smaller launch jump. This run, which uses a trailing median instead of a baseline anchored at liftoff, scores the two scans at 39.29 K and 45.37 K. The 0.04 K difference on the second scan is the median moving from 300.14 K to 300.18 K once the earliest pre-liftoff sample drops out of the five-scan window.

On 25 October the same rule's largest score is 0.04 K, at 23:26:26 UTC. Milestone 2's 0.21 K used a baseline fixed at 23:23:45, so the warmest scan of that evening was compared with later, slightly cooler scans. The trailing baseline never builds that comparison, because those first scans are still filling the window. Nothing on the quiet evening crosses 20 K.

After the two launch scans the neighbourhood maximum returns to about 300–301 K. The largest later score is 0.75 K, at 23:34:54.9.

## What is logged

A scan at or above 20 K is a candidate. The text report and the JSON line both keep the watch name, scan start and end, sector, source object key, hot-pixel latitude and longitude, distance from the pad, score, baseline and the samples in it, the pad pixel, the 98th percentile, the 10–30 km median, solar zenith, DQF counts, the saturated flag, and the fire-product result.

Both Flight 11 candidates are emitted:

| Scan start | Score | Hot pixel | Distance | Fire product |
| --- | --- | --- | --- | --- |
| 23:23:54.9 | 39.29 K | 25.987°N, 97.161°W, 339.43 K | 1.24 km | no code 10–15 or 30–35 within 20 km |
| 23:24:54.9 | 45.37 K | 26.075°N, 97.205°W, 345.55 K | 10.09 km | code 10, 3.39 km from the hot pixel, 17.7 MW |

DQF inside 20 km is 0 on all 227 pixels of both scans. Neither pixel is within 0.5 K of the 410.92 K ceiling. Solar zenith is 82.1° and 82.3°.

## Persistence and the fire product

The third consecutive covering scan at or above 20 K is held back. A gap longer than 90 s, against a 60 s cadence at this pad, starts the count over. Flight 11 is hot for two scans, so both stay emitted. The log still records the count.

The fire file is the mesoscale product `ABI-L2-FDCM` for the same sector and the same scan start. `ABI-L2-FDCM2` exists for every M2 scan in this hour. Codes 30–35 (`temporally_filtered_*`) within 20 km of the hot pixel hold the candidate back. Codes 10–15 are written down and do not hold it back.

That split is from this launch. The scan starting 23:24:54.9 contains one `good_fire_pixel` (code 10) at 26.065°N, 97.237°W, 17.7 MW, 3.39 km from the band-7 hot pixel. The FDCM pixel nearest that hot pixel is code 150, the land-cover class the file calls sea water or bright desert, as is the pixel nearest the pad on the previous scan. The pad sits on the coastline, and the fire algorithm does not treat that class as a place a fire can be. The plume is invisible to FDCM while it is on the pad. One scan later a neighbouring cell is a good fire. Holding every fire code back would have dropped the second Flight 11 candidate. A code 30 pixel does exist in that same FDCM file, 722 km from the pad, and the 20 km test leaves it out.

A missing fire file is `fire_file=unavailable`. The candidate stays in the log.

## The ongoing fire, scored with the same rule

The mask-30 pixel from Milestone 2 (32.810°N, 82.759°W, 166.2 MW) is not a watch site. The script scores it anyway, on CONUS band 7 from 21:01 to 23:56 UTC, so the threshold can be seen against a fire that is already burning.

The pixel is 342.79 K on the first scan. By the time three previous scans exist, the baseline is already the fire. The largest score in the three hours is 12.36 K, at 22:11:17 UTC: 335.17 K minus a baseline of 322.81 K, after the pixel had cooled and warmed again. That is under 20 K, so the fire check emits nothing. The jump rule sees the start of a hot interval. It does not see a hot pixel whose recent scans are just as hot.

## What this run does not settle

The 20 K and 20 km numbers come from one clear launch and one clear evening twelve days later. Clouds were in neither sample. A fire that ignites while the site is watched can pass for two scans before the persistence count holds the third. A fire that flares by more than 20 K above its own recent median would also pass; this fire's largest flare was 12.36 K. The coast mask can hide a plume that stays on the pad. A minute with no mesoscale sector on the pad is a coverage gap, and the score does not fall back to CONUS. No second watch site is armed.
