# How Milestone 2 compares the signal

`scripts/m2_signal.py` puts the Flight 11 band-7 pulse next to a quiet evening at the same pad and next to a fire on the same calendar evening. It writes:

- `outputs/m2_signal.png`
- `outputs/m2_signal.txt`

The script imports the Milestone 1 geolocation and Planck conversion, runs that self-test, and reuses the mesoscale files already stored under `data/m1/`. New downloads go under `data/m2/`. The numbers below are the ones that run printed. Re-run from the repository root to regenerate them:

```bash
.venv/bin/python scripts/m2_signal.py
```

No detection threshold is chosen here. Milestone 3 freezes a rule. The radii and the 10 K count in this note are the windows used to see this one case.

## The three scenes

| Scene | Product | When | Point |
| --- | --- | --- | --- |
| Launch | ABI-L1b-RadM band 7, sector M2 | 13 Oct 2025, 23:18–23:36 UTC | Starbase Pad 1, 25.9962°N, 97.1544°W |
| Quiet evening | ABI-L1b-RadM band 7, sector M1 | 25 Oct 2025, the same clock window | the same pad |
| Fire | ABI-L1b-RadC band 7 | 13 Oct 2025, 21:01–23:56 UTC, every 5 min | 32.810°N, 82.759°W |

Brightness temperature uses `planck_fk1`, `planck_fk2`, `planck_bc1`, and `planck_bc2` from each file. Reflected sunlight stays in the number. The formula is the one in [Radiance to brightness temperature](radiance-to-brightness-temperature.md).

The quiet evening is 25 October because the mesoscale sectors were tasked elsewhere on the adjacent evenings. The script's search, one band-7 file per sector nearest 23:22 UTC from 27 September through 28 October, found the pad inside a sector on 24 October (M1) and 25 October (M1). It did not find the pad on 12 October or 14 October. 25 October is the same product at the same clock time. Solar zenith at the pad at 23:23:45 UTC was 84.6° that evening and 82.1° on launch evening, so the Sun was 2.5° closer to the horizon.

The fire is the strongest mask-30 pixel at least 100 km from the pad in one pinned CONUS fire file, `OR_ABI-L2-FDCC-M6_G19_s20252862321170_e20252862323543_c20252862324125.nc`. Mask 30 is `temporally_filtered_good_fire_pixel`. The selected pixel has a product fire radiative power of 166.2 MW, a product fire temperature of 716 K, and a product fire area of 18757 m². It is 1585 km from the pad. OpenStreetMap places that coordinate in Washington County, Georgia. The label stays a county: the band-7 pixel is about 2 km across. The 716 K figure is the fire algorithm's sub-pixel temperature. Every temperature in the time series is the band-7 brightness temperature of the whole pixel.

## What each scan contributes

`measure_point` cuts a ±30-pixel window in scan-angle index space, about 70 km on a side at the pad, and then measures distance with the same 6371.0088 km sphere as Milestone 1.

| Column | Meaning |
| --- | --- |
| `near_K` | brightness temperature of the pixel nearest the point |
| `p98_K` | 98th percentile of finite pixels within 20 km |
| `max20_K` | hottest finite pixel whose centre is within 20 km |
| `dmax_K` | `max20_K` minus `max20_K` of the previous scan that saw the point |
| `ring_K` | median brightness temperature of pixels from 10 km to 30 km away |
| jump | `max20_K` minus the median `max20_K` of scans that had ended by 23:23:45 that day |

The 20 km disk at the pad is on the order of 200 pixels, so one hot pixel moves the maximum by tens of kelvin and the 98th percentile by about 5 K.

## Launch evening against the quiet evening

Five M2 scans end at or before liftoff. Their 20 km maximum has median 300.14 K and sample standard deviation 0.10 K. The pad pixel over those scans has median 299.31 K and sample standard deviation 0.04 K.

Two scans then leave that baseline:

| Scan start | Seconds after liftoff | Pad pixel | Max within 20 km | Distance |
| --- | --- | --- | --- | --- |
| 23:23:54.9 | +9.9 s | 339.43 K | 339.43 K | 1.24 km |
| 23:24:54.9 | +69.9 s | 299.69 K | 345.55 K | 10.09 km, at 26.075°N, 97.205°W |

The jump is the 20 km maximum minus 300.14 K. That is 39.29 K on the first hot scan and 45.41 K on the second. The 98th percentile within 20 km reaches 305.02 K on the first of those scans. On the last pre-liftoff scan it is 299.81 K. After the second hot scan the 20 km maximum returns to about 300–301 K. The scaled brightness-temperature ceiling on these files is 410.92 K, and the scene does not reach it.

On 25 October the same clock window has 18 pad-seeing scans, about 60 s apart. The 20 km maximum has median 298.39 K, standard deviation 0.20 K, minimum 298.20 K, and maximum 298.84 K. The largest scan-to-scan rise is 0.09 K and the largest fall is 0.22 K. Applying the launch-day rule (median of scans ended by 23:23:45, then the peak excess) gives a baseline of 298.63 K and a peak jump of 0.21 K.

Two ratios describe this pair of evenings. The 45.41 K launch jump is 231 times the quiet-evening standard deviation of 0.20 K, and 523 times the largest quiet-evening rise of 0.09 K. The quiet-evening median, 298.39 K, sits 1.74 K below the launch-evening baseline of 300.14 K. That day-to-day offset is larger than the 0.20 K within-evening scatter. An absolute brightness temperature that fits one clear evening does not fit the other, while a jump against the scans just before liftoff does.

## Pixels farther than 20 km

The ±30-pixel window also records the hottest pixel outside 20 km. Two scans exceed the 20 km baseline by more than 5 K there:

| Seconds after liftoff | Brightness temperature | Distance | Above the 20 km baseline |
| --- | --- | --- | --- |
| +129.9 s | 308.26 K | 25.4 km | 8.12 K, 41 quiet-evening standard deviations |
| +189.9 s | 307.54 K | 57.8 km | 7.40 K, 38 quiet-evening standard deviations |

Those are the same warm pixels Milestone 1 reported downrange. They are large compared with the 0.20 K quiet evening, and they are outside the 20 km neighbourhood that holds the two launch scans. This quiet evening never produced a 308 K pixel near the pad. The sample is one clear evening. It contains no cloudy frames, so the comparison does not yet show what a cold cloud or another warm spot does to that tail.

## The fire

CONUS band 7 revisited the fire pixel every 300.0 s. All 36 files in hours 21–23 contain the point. The nearest pixel is the 20 km maximum on every scan, 0.00 km from the product coordinate. Pixel spacing there is 2.10 km east–west and 2.65 km north–south. Satellite zenith is 38.9°, against 39.0° at the pad.

The anomaly used for the fire is the 20 km maximum minus the 10–30 km ring median. It peaks at 49.20 K at 23:01:17 UTC: the fire pixel is 339.91 K and the ring is 290.71 K. All 36 scans exceed the ring by more than 10 K, from 21:01 to 23:56, a span of 175 minutes. The median anomaly is 42.83 K and the minimum is 24.37 K. The largest scan-to-scan rise of the neighbourhood maximum is 7.31 K and the largest fall is 7.41 K. The ring itself cools from about 297 K to about 289 K as the night deepens. The first and last scans are still hot (342.79 K and 328.69 K), so 175 minutes is a lower bound on how long this fire was bright. The 10 K level only counts that duration.

Solar zenith at the fire pixel at 23:23:45 UTC is 95.6°. The pad at liftoff is still in daylight at 82.1°. The backgrounds are different mixtures: the pad still includes reflected sunlight, and the fire pixel is in darkness. The peak band-7 values are similar, about 346 K for the launch pixel and about 340 K for the fire pixel, while the fire product's own temperature is 716 K. What separates them in this sample is duration. The launch occupies two 60-second mesoscale scans. The fire occupies every 5-minute CONUS scan in a three-hour window.

## What stays undecided

The 20 km neighbourhood, the 10–30 km ring, the 10 K duration count, and the ratios 231 and 523 describe Flight 11, one quiet evening, and one fire. Milestone 3 still chooses the score, the radius, and which fire product to consult. This run used the CONUS fire product because that is where the bright, persistent pixel was. It did not score the mesoscale fire product.
