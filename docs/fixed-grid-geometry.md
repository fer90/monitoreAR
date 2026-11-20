# Fixed-grid geometry

ABI Level 1b geolocation is a fixed grid. Each pixel is a scan angle, in radians, as seen from an ideal satellite. `x` and `y` in the file are those angles. xarray has already applied their scale and offset, and `load_frame` refuses a file whose `x` units are anything other than radians. It also refuses a file whose `sweep_angle_axis` is anything other than `x`, because the equations below assume that sweep.

The projection variables live on `goes_imager_projection`. `projection_constants` reads:

| Attribute | Value in the Flight 11 files | Meaning |
| --- | --- | --- |
| `semi_major_axis` | 6378137.0 m | GRS80 equatorial radius, `r_eq` |
| `semi_minor_axis` | 6356752.31414 m | GRS80 polar radius, `r_pol` |
| `perspective_point_height` | 35786023.0 m | Height of the ideal satellite above the ellipsoid, `h_persp` |
| `longitude_of_projection_origin` | −75.0° | Longitude of the fixed grid, `λ0` |
| `sweep_angle_axis` | `x` | Which angle is measured first |
| `grid_mapping_name` | `geostationary` | |

The satellite's distance from Earth's centre in these equations is

```
H = h_persp + r_eq = 35786023 + 6378137 = 42164160 m
```

`nominal_satellite_subpoint_lon` on the same file is −75.200° (the printed value is −75.19999694824219), and `nominal_satellite_subpoint_lat` is 0. `nominal_satellite_height` is 35786.0234375 km above the ellipsoid. The fixed-grid origin is the ideal longitude −75.0°. The navigation uses `longitude_of_projection_origin` from the projection variable, which is −75.0°, together with the ellipsoid axes and `perspective_point_height`. `yaw_flip_flag` is 0.

The forward and inverse maps are `latlon_to_scan_angles` and `scan_angles_to_latlon`. Both follow the GOES-R Product User's Guide, Volume 3, L1b products, section 5.1.2.8, as implemented in the script. The intersection is with the GRS80 ellipsoid.

## Latitude and longitude to scan angles

Inputs are geodetic latitude `φ`, longitude `λ`, and the constants above.

Geocentric latitude, and the distance from Earth's centre to the ellipsoid at that latitude:

```
φ_c = arctan( (r_pol² / r_eq²) · tan φ )
e²  = 1 − (r_pol / r_eq)²
r_c = r_pol / sqrt(1 − e² cos² φ_c)
```

The vector from the ideal satellite to that surface point:

```
s_x = H − r_c cos φ_c cos(λ − λ0)
s_y = −r_c cos φ_c sin(λ − λ0)
s_z = r_c sin φ_c
```

A point past the limb fails

```
H · (H − s_x) ≥ s_y² + (r_eq² / r_pol²) · s_z²
```

and the forward function raises. For a point on the visible Earth the scan angles, with `sweep_angle_axis = x`, are

```
y = arctan(s_z / s_x)
x = arcsin( −s_y / sqrt(s_x² + s_y² + s_z²) )
```

`load_frame` uses this once for the pad. The pad is inside the sector when `x_pad` lies between the minimum and maximum of the file's `x`, and the same for `y`. That test is what drops every M1 file in this window. Their scan-angle box is over the Atlantic, so the pad's angles fall outside it.

The nearest pixel is the column with `x` closest to `x_pad` and the row with `y` closest to `y_pad`. Those two choices are independent because the coordinates are separable. On the first covering frame that pixel centre projects to 25.9870°N, 97.1614°W, which is 1.238 km from the pad coordinates. The neighbouring pixels are 2.31 km and 2.47 km away, so the pad lies inside that pixel.

## Scan angles back to latitude and longitude

`scan_angles_to_latlon` intersects the look direction with the ellipsoid. The range `r_s` along the look is the near root of `a t² + b t + c = 0`:

```
a = sin²x + cos²x · ( cos²y + (r_eq² / r_pol²) sin²y )
b = −2 H cos x cos y
c = H² − r_eq²
```

A negative discriminant means the look misses the Earth. Those pixels come back as NaN. Otherwise

```
r_s = (−b − sqrt(b² − 4ac)) / (2a)
```

The minus sign in front of the square root selects the nearer intersection. The vector from the satellite to that point is

```
s_x = r_s cos x cos y
s_y = −r_s sin x
s_z = r_s cos x sin y
```

and the geodetic coordinates are

```
lat = arctan( (r_eq² / r_pol²) · s_z / sqrt( (H − s_x)² + s_y² ) )
lon = λ0 − arctan2( s_y , H − s_x )
```

The script calls this three times on a covering frame: once for every pixel centre (to crop and to measure distances), once for the cell edges (to draw the mesh), and once for the sector-centre check below.

Nadir is the special case `(x, y) = (0, 0)`. It maps to latitude 0 and longitude `λ0`.

## Two checks, and what each one measures

`self_test_geolocation` runs before the download. It uses constants written at the top of the script, not the file:

```
_TEST_R_EQ = 6378137.0
_TEST_R_POL = 6356752.31414
_TEST_H    = 35786023.0          # passed in as h_persp
_TEST_LON0 = -75.2               # nominal subpoint, for the test only
```

It checks that `(0, 0)` returns latitude 0 and longitude −75.2°, then projects the pad to scan angles and back. The round trip on those constants is 0.00 m. The script stops if that distance exceeds 200 m.

`_center_residual_m` uses the file. It converts the file's `x_image` and `y_image` with the file's own projection constants and compares the result with `geospatial_lat_center` and `geospatial_lon_center`. On the first covering frame the difference is 0.182 m.

The 0.182 m figure says the scan-angle equations reproduce the latitude and longitude NOAA stored for the centre of this sector. Spacecraft pointing error is a separate residual, and this milestone does not measure it. Both residuals are far smaller than the 2 km pixel.

## Where an elevated source is drawn

The inverse intersection is with the ellipsoid. A plume at altitude is recorded where its line of sight meets that ellipsoid. The marker in the figure is the pad on that same ellipsoid. The script leaves the hot pixel where the fixed grid puts it. The size and direction of the shift are estimated in [Pixel size, zenith, parallax, and the Sun](viewing-geometry.md). The estimate is a scale, and it is not applied as a correction.

## Cell edges in the figure

`edges_1d` turns a row of pixel-centre scan angles into cell boundaries: the interior boundaries are midpoints, and each end is extended by half of the neighbouring step. `plot_frames` converts those edge angles to latitude and longitude and draws `pcolormesh` with `shading="flat"`. The cyan pad marker is plotted at the true pad latitude and longitude, so it falls in the cell whose centre is the nearest pixel.

The axes are equal in degrees after `set_aspect(1 / cos(pad latitude))`, which shortens a degree of longitude the way a local map does. The coastline in the picture is the boundary between warm land and cooler water. The script draws no shoreline.
