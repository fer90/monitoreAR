# Pixel size, zenith, parallax, and the Sun

The fixed grid places each pixel on the ellipsoid. This note covers the four extra numbers printed under "At the pad, from the fixed grid" in [`outputs/m1_inspection.txt`](../outputs/m1_inspection.txt): how large the pixel is, how high GOES-19 sits in the sky there, how far an elevated source shifts, and where the Sun was at liftoff.

All four use the same pad, 25.9962°N, 97.1544°W, and the projection of the first frame that contains it. The functions are `haversine_km`, `satellite_zenith_deg`, `initial_bearing_deg`, and `solar_zenith_deg`.

## Distances: a sphere of radius 6371.0088 km

`haversine_km` is the great-circle distance

```
a = sin²(Δφ / 2) + cos φ1 cos φ2 sin²(Δλ / 2)
d = 2 R arcsin(sqrt(a))
```

with `R = 6371.0088` km. Latitude and longitude are in degrees on the way in. The script uses this radius for three measurements:

- pad coordinates to the nearest pixel centre: 1.238 km
- that pixel centre to its east-west neighbour, and to its north-south neighbour
- every crop pixel to the pad, including the test "is this centre within 20 km?"

The east-west neighbour of the pad pixel is 2.31 km away. The north-south neighbour is 2.47 km away. The file attribute `spatial_resolution` says "2km at nadir". At Starbase the measured spacing is larger than that nadir figure, and the two axes differ. The script reports the measured neighbour distances.

Because the nearest centre is 1.238 km from the pad and the pixel is about 2.4 km across, the pad coordinates fall inside that pixel. A 10 m uncertainty in the pad survey is far inside the same pixel.

"Max within 20 km" is the hottest finite crop pixel whose centre is at most 20 km from the pad. "Crop maximum" is the hottest pixel anywhere in the map window, which extends well past 20 km. On the +70 s scan those two are the same pixel. On the +130 s and +190 s scans the crop maximum is outside 20 km, so it shows up only in the crop-maximum list. The table in the report is the 20 km figure. Both rules are applied in [How the figure and the report are chosen](figure-and-report.md).

## Satellite zenith angle

`satellite_zenith_deg` is a spherical approximation, built from the same geocentric latitude as the fixed grid. It explains the footprint and the parallax. It is not a navigation product.

With geocentric latitude `φ_c` and geocentric radius `r_c` from [Fixed-grid geometry](fixed-grid-geometry.md), and `H` the satellite's distance from Earth's centre:

```
cos β = cos φ_c · cos(λ − λ0)
β     = arccos(cos β)
```

`β` is the angle at Earth's centre between the subsatellite point and the pad. The elevation of the satellite above the local horizontal is

```
elevation = atan2( cos β − r_c / H , sin β )
zenith    = 90° − elevation
```

At the pad this prints 39.0°. GOES-East is over the equator at about 75°W, and Starbase is 22° of longitude further west and 26° north, so the satellite is well off vertical. That is the viewing angle behind the 2.31 km by 2.47 km spacing.

## Parallax

An elevated source is seen along a slanted line of sight. The fixed grid reports the point where that line meets the ellipsoid, which is farther from the subsatellite point than the point directly below the source. For a source at altitude `h` the script estimates the ground shift as

```
shift ≈ h · tan(zenith)
```

and prints it for `h = 10` km:

```
10 km × tan(39.0°) ≈ 8.1 km
```

The direction of that shift is the initial great-circle bearing from the subsatellite point `(0°, λ0)` to the pad, from `initial_bearing_deg`:

```
y = sin(Δλ) cos φ2
x = cos φ1 sin φ2 − sin φ1 cos φ2 cos(Δλ)
bearing = atan2(y, x)    # degrees from north, wrapped to [0, 360)
```

With `λ0 = −75°` that bearing is 322°. The apparent position moves away from the subsatellite point, toward the northwest.

The +70 s hot pixel is at 26.075°N, 97.205°W, 10.09 km from the pad. That is northwest of the pad, which is the same general direction as the parallax bearing. The pixel can contain horizontal motion of the rocket and the parallax of an elevated plume at the same time. The script reports the 8.1 km scale and the measured 10.09 km offset. It does not split that offset into an altitude and a ground track.

## Solar zenith angle

`solar_zenith_deg` is the NOAA general solar-position approximation. The docstring treats it as good to about one degree, which is enough to decide whether band 7 still contains reflected sunlight.

For a UTC time and a longitude in degrees east (the pad longitude is already negative):

```
hour  = hours + minutes/60 + seconds/3600
γ     = 2π/365 · (day_of_year − 1 + (hour − 12)/24)

eqtime = 229.18 · (
    0.000075
    + 0.001868 cos γ − 0.032077 sin γ
    − 0.014615 cos 2γ − 0.040849 sin 2γ
)                                                    # minutes

decl = 0.006918
     − 0.399912 cos γ + 0.070257 sin γ
     − 0.006758 cos 2γ + 0.000907 sin 2γ
     − 0.002697 cos 3γ + 0.00148 sin 3γ             # radians

time_offset = eqtime + 4 · longitude_deg            # minutes; the zone is UTC
hour_angle  = (hour·60 + time_offset) / 4 − 180     # degrees, then converted

cos z = sin φ sin(decl) + cos φ cos(decl) cos(hour_angle)
zenith = arccos(cos z)
```

At the pad, at the quoted liftoff, this is 82.1°. The Sun is about 8° above the horizon. Band 7 brightness temperature includes the reflected solar part of the radiance for the whole window. The coefficients that would remove it are absent; that is covered in [Radiance to brightness temperature](radiance-to-brightness-temperature.md).

The background brightness temperature drifts by only a few tenths of a kelvin from 23:18 to 23:36 (the crop median goes from 298.29 K to 297.76 K). The solar zenith changes slowly on that timescale. The 40 K rise of a single pixel between the 23:22:54.9 scan and the 23:23:54.9 scan sits on top of that stable background.
