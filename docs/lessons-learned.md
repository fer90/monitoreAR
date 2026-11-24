# Lessons from the Flight 11 case study

These are the measurements that fixed the v0.1 rule. The tables they come from are [`outputs/m1_inspection.txt`](../outputs/m1_inspection.txt), [`outputs/m2_signal.txt`](../outputs/m2_signal.txt), [`outputs/m3_score.txt`](../outputs/m3_score.txt), and [`outputs/m4_poll.txt`](../outputs/m4_poll.txt). The 20 K cut is a constant chosen from this one launch. Catching Flight 11 again is a consistency check. It is not an independent detection.

## A pixel is about 2 km, and brightness temperature is not a flame temperature

At Pad 1 the nearest band-7 pixel is about 2.3 km east–west and 2.5 km north–south. The satellite zenith angle is about 39°. Brightness temperature is the Planck temperature of the radiance in that pixel. The file's Planck coefficients do not remove reflected sunlight. A 339 K pixel is a warm average over a few square kilometres. It is not the temperature of the flame. Nothing in the launch sector reached the radiance ceiling near 411 K.

The coordinate is where the line of sight meets the ellipsoid. A hot source above the surface is left where that line meets the ground. Parallax is about 8 km of horizontal shift for each 10 km of altitude, toward the satellite, near bearing 322°. This code does not slide the plume back onto the pad. On the second hot scan the hottest pixel inside 20 km was 10.09 km from Pad 1.

## Absolute brightness temperature does not travel from one day to the next

The quiet evening of 25 October 2025, the nearest later evening with a mesoscale sector on the pad at this clock time, sat about 1.74 K away from the launch evening's background. The neighbourhood maximum on that quiet evening had a standard deviation of 0.20 K. A fixed brightness-temperature gate would move with the day. The score is a jump: the hottest usable pixel within 20 km minus the median of that maximum on the previous covering scans.

## The 20 K cut was chosen from this one launch

`THRESHOLD_K` is 20.0. It is not computed at runtime, and it is not a mean-plus-sigma cut. The quiet evening's peak jump in the Milestone 2 fixed-median sense was 0.21 K. Flight 11 jumped about 39 K and then about 45 K. Pixels about 8 K warmer sat 25 km and 58 km from the pad, outside the 20 km radius. The trailing-median scores stored for the launch are 39.29 K and 45.37 K. The same rule on the Georgia fire, which is not a watch site, peaked at 12.36 K and emitted nothing. The quiet evening's largest trailing score is 0.04 K.

The poller and the API do not refit 20 K.

## Only one mesoscale sector covered the pad

Mode 6 scans two mesoscale sectors about once a minute. On the launch evening only sector M2 contained Pad 1. Sector M1 was over the tropical Atlantic. A sector that misses the pad is a coverage gap. v0.1 does not fall back to the CONUS or full-disk product. The live passes on 23 November 2026 stored no candidate for the same reason: neither sector contained the pad.

## The fire product and the plume

ABI-L2-FDCM for the same sector and the same scan start is a hold, and only for mask codes 30–35 within 20 km of the hot pixel. Codes 10–15 are logged and do not hold. On the second Flight 11 hot scan the plume itself is code 10, 17.74 MW, 3.39 km from the hot pixel. The nearest fire-mask pixel to that hot pixel, and the pad pixel, are code 150 (sea water or bright desert). Holding every fire code would have dropped the launch. A missing fire file leaves the candidate in place.

The comparison fire was the strongest ABI-L2-FDCC mask-30 pixel at least 100 km from the pad on the CONUS file that overlaps liftoff. The product does not say wildfire versus industry. The county name is a later geocode. That pixel was already hot at the start of the download window, so the file does not give an ignition time.

## The bucket listing is cached, and the filename clock is not the arrival time

`s3fs` caches `fs.ls`. That method has no refresh argument. Each list drops the cache for the prefix first. A poll that skipped the drop repeated the previous listing; the 2.4 minute age in that run was the cache. After the drop, the same evening showed a newest scan tens of seconds old.

The `c` token in the filename is product-creation time at the ground station. It is not the time the object appeared on AWS. Scan start and scan end are the instrument clock.

## A median of five keeps the second hot scan visible

The baseline is the median of up to five previous neighbourhood maxima, and it needs at least three. History is appended after the current scan is scored. One or two hot scans do not set that median. The second Flight 11 score uses a baseline of 300.18 K, so the score is 45.37 K. Milestone 2's 45.41 K used the fixed pre-liftoff median of 300.14 K. Both numbers match their definitions. A previous sample older than 30 minutes stays out of the poller's baseline, so the launch scans do not become the background for a later day.

## Persistence of three would hold a longer burn

The third consecutive covering scan at or above 20 K is held back. A gap longer than 90 seconds breaks the run. Flight 11 was hot for exactly two scans, so both are emitted. A burn that stayed hot for a third scan would be held from that scan on. The already-burning fire stayed hot for hours and still stayed under 20 K, because the score is a jump from its own recent median, not the absolute brightness temperature.

## What this case study does not show

One clear evening does not say what a launch under cloud would do, or what a fire would do if it flared by more than 20 K above its own recent median. The colour scale on the map is the Flight 11 figure, 283.43 K to 350 K. The page is local. There is no account, no write API, and no alert.
