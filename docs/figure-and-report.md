# How the figure and the report are chosen

`main` writes two artifacts:

- `outputs/m1_starship_flight11_band7.png`, from `plot_frames`
- `outputs/m1_inspection.txt`, from `write_report`

Both use only frames whose fixed grid contains the pad. In this window that is the 18 sector-M2 scans. The colour scale, the white cross, and the two kinds of maximum are rules in the script. The numbers below are the ones that rule selected for Flight 11.

## The three panels

The left panel is the last covering frame whose scan **end** is at or before the quoted liftoff. That is the M2 file `...s20252862322549_e20252862323018_...`: scan start 23:22:54.9 UTC, scan end 23:23:01.8 UTC, pad pixel 299.27 K. Liftoff at 23:23:45 falls after that scan has finished. The next look at the pad starts at 23:23:54.9, 10 s after liftoff.

The other two panels are chosen from covering scans whose scan **start** lies from 23:23:45 to 23:27:45 inclusive. The starts in that interval are 23:23:54.9, 23:24:54.9, 23:25:54.9, and 23:26:54.9. The following scan, 23:27:54.9, starts after the window. Inside the four scans:

- the middle panel maximises `pad_pixel_bt`, the brightness temperature of the single nearest pixel
- the right panel maximises `max_bt_20km`, the hottest pixel whose centre is within 20 km of the pad

Those two maxima fall on different scans.

| Panel | Scan start | File | What it shows |
| --- | --- | --- | --- |
| before liftoff | 23:22:54.9 UTC (−50 s) | `...M2...s20252862322549_e20252862323018_c20252862323060.nc` | Pad pixel 299.3 K. Hottest pixel within 20 km is 300.0 K, 15.9 km away. |
| pad pixel heats | 23:23:54.9 UTC (+10 s) | `...M2...s20252862323549_e20252862324019_c20252862324058.nc` | Pad pixel and the 20 km maximum are the same pixel: 339.43 K, 1.24 km from the pad coordinates. |
| hot pixel has moved | 23:24:54.9 UTC (+70 s) | `...M2...s20252862324549_e20252862325018_c20252862325056.nc` | Pad pixel is back to 299.69 K. The 20 km maximum is 345.55 K at 26.075°N, 97.205°W, 10.09 km northwest of the pad. |

The titles print sector, scan start, seconds relative to liftoff, the pad-pixel temperature, and the 20 km maximum with its distance.

## Colour scale

`vmin` is the 2nd percentile of the **before** crop, 283.43 K. `vmax` is fixed at 350 K. The colormap is `inferno`.

The pre-liftoff crop already has its 98th percentile at 299.81 K, and the land sits near 299–300 K. A scale that ran from the 2nd to the 98th percentile of that crop would put the land and a 340 K plume pixel at the same top colour. The fixed top at 350 K sits above the hottest pixel in the window, 345.55 K, and it leaves the ~300 K surface and the ~340 K pixel as different colours. The cold end, 283 K, keeps the colder cloud tops visible at the bottom of the scale.

## The white cross

The cyan ring is always the pad, at `PAD_LAT`, `PAD_LON`.

The white cross is the hottest pixel within 20 km, and only when both of these hold:

```
max_bt_20km > median_bt + 10 K
max_dist_km > 3
```

A warm-background pixel is left unmarked. Applying that rule to the three panels:

| Panel | 20 km max above the median by | Distance | Cross |
| --- | --- | --- | --- |
| before | 300.01 − 298.20 = 1.8 K | 15.9 km | no |
| pad pixel heats | 339.43 − 298.20 = 41.2 K | 1.24 km | no, because the hot pixel is the pad pixel |
| hot pixel has moved | 345.55 − 298.15 = 47.4 K | 10.09 km | yes |

## Two maxima

The per-frame table reports the 20 km neighbourhood. Columns are scan start, sector, seconds relative to liftoff, pad-pixel temperature, 20 km maximum, distance of that maximum, crop median, crop 98th percentile, the file-level saturated-pixel count, and the DQF counts.

Through 23:22:54.9 the pad pixel stays at 299.2–299.3 K and the 20 km maximum stays near 300 K, about 16 km away. At +10 s both jump to 339.43 K on the pad pixel. At +70 s the pad pixel is 299.69 K and the 20 km maximum is 345.55 K at 10.09 km. From +130 s onward the 20 km maximum is 299.8–301.3 K.

The crop maximum is a second measurement. The report lists every covering frame whose hottest crop pixel exceeds 305 K:

| Scan start | Crop maximum | Where |
| --- | --- | --- |
| 23:23:54.9 (+10 s) | 339.4 K | 1.2 km, the pad pixel (25.987°N, 97.161°W) |
| 23:24:54.9 (+70 s) | 345.5 K | 10.1 km northwest (26.075°N, 97.205°W) |
| 23:25:54.9 (+130 s) | 308.3 K | 25.4 km (26.224°N, 97.126°W), outside the 20 km neighbourhood |
| 23:26:54.9 (+190 s) | 307.5 K | 57.8 km (26.409°N, 96.802°W), northeast, over the Gulf |

After those four scans the crop maximum stays near the 298–300 K background. The +130 s and +190 s pixels are a few kelvin above that background and tens of kilometres from the pad. The report prints them as crop maxima above 305 K. The 20 km column at those times is 299.81 K and 299.77 K.

A boostback reported elsewhere around 23:29 falls inside this window. The 20 km maximum on the 23:28:54.9 and 23:29:54.9 scans is 300.54 K and 300.94 K. The script records those values. It does not mark them as a second hot pixel.

## Quality lines at the bottom of the report

The closing lines compare the hottest 20 km value in the window, 345.5 K, with the scaled radiance ceiling, 410.922 K. The scene stayed below the ceiling. The last sentence of the report states what the temperature is: the pixel average, not the flame temperature. Both points are derived in [Radiance to brightness temperature](radiance-to-brightness-temperature.md).

DQF counts are printed and are all `good_pixel_qf` in this crop. They do not modify the image. CONUS and Full Disk appear in the report as start times only. The cadence, the 60 s gap, and the product-creation lag are measured from the M2 filenames, as described in [From the bucket to one frame](data-flow.md).
