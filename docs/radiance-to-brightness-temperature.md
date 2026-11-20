# Radiance to brightness temperature

ABI Level 1b stores radiance. Brightness temperature is computed in `radiance_to_bt` from coefficients stored in the same file. The Level 2 Cloud and Moisture Imagery product already contains brightness temperature; this repository starts from radiance so that step stays visible.

## What the file stores

On the first M2 frame that contains the pad (`OR_ABI-L1b-RadM2-M6C07_G19_s20252862318549_...`):

| Field | Value |
| --- | --- |
| `band_id` | 7 |
| `band_wavelength` | 3.890000104904175 µm |
| `Rad` units | mW m⁻² sr⁻¹ (cm⁻¹)⁻¹ |
| `Rad` `standard_name` | `toa_outgoing_radiance_per_unit_wavenumber` |
| `planck_fk1` | 200248.0, in the same radiance units |
| `planck_fk2` | 3685.8701171875 K |
| `planck_bc1` | 0.4898500144481659 K |
| `planck_bc2` | 0.9993100166320801 |
| `esun` | NaN |
| `kappa0` | NaN |

`planck_fk1`, `planck_fk2`, `planck_bc1`, and `planck_bc2` are scalar variables in the NetCDF file. On these GOES-19 files they are not attributes of `Rad`. `load_frame` reads them with `float(ds["planck_fk1"].values)` and the same for the other three.

The variable attributes name the coefficients:

- `planck_fk1` is the file's factor \(2hc^{2}\nu^{3}\), already expressed in mW m⁻² sr⁻¹ (cm⁻¹)⁻¹. \(\nu\) is the band's central wavenumber.
- `planck_fk2` is the file's factor \(hc\nu/k\), in kelvin. The attribute text writes this as `h c nu/b`, with \(b\) Boltzmann's constant.
- `planck_bc1` and `planck_bc2` are the spectral bandpass correction. The attribute writes the corrected temperature as \((B(\nu) - b_{c1}) / b_{c2}\), where \(B(\nu)\) is the monochromatic Planck temperature. \(b_{c2}\) is already 0.99931, so the correction is a few tenths of a kelvin.

xarray applies `scale_factor` and `add_offset` before `.values` is returned. `Rad.values` is therefore physical radiance in the units above. The scale and offset still matter for the ceiling, below.

## The formula

```
T = (fk2 / ln(fk1 / L + 1) - bc1) / bc2
```

`L` is one pixel of `Rad`. The logarithm is natural. `radiance_to_bt` writes NaN wherever `L` is non-finite or `L <= 0`, and applies the formula everywhere else.

The monochromatic piece is the Planck function solved for temperature. A blackbody spectral radiance per unit wavenumber has the shape

```
Bν(T) = 2hc²ν³ / (exp(hcν / kT) - 1)
```

Set `Bν = L`, insert `fk1` for `2hc²ν³` and `fk2` for `hcν/k`, and solve:

```
T_mono = fk2 / ln(fk1 / L + 1)
T      = (T_mono - bc1) / bc2
```

The second line is the bandpass correction from the file attributes. Band 7 is a finite filter, so the monochromatic inversion at the central wavenumber is adjusted by `bc1` and `bc2`.

The same coefficients turn the measured brightness temperatures back into radiance. That inverse is how the formula behaves; the script stores temperature, and the inspection report does not print `L`:

| Brightness temperature | Radiance from the formula, same units as `Rad` |
| --- | --- |
| 299.27 K, pad pixel on the last pre-liftoff scan | 0.907 |
| 339.43 K, pad pixel at +10 s | 3.883 |
| 345.55 K, hottest pixel within 20 km at +70 s | 4.704 |

A 40 K rise at this wavelength is about a fourfold rise in radiance. The Planck curve is steep at 3.9 µm, which is why a hot object that fills only part of a pixel can still move the pixel a long way in temperature.

## What 339 K is

339.43 K is the brightness temperature of the whole pixel: the temperature of a blackbody that would send the satellite the radiance that pixel recorded. The Raptor engines are much hotter than that, and they are much smaller than the pixel that contains them. The pixel is 2.31 km by 2.47 km at the pad (see [Pixel size, zenith, parallax, and the Sun](viewing-geometry.md)), so most of its area is ground, smoke, and sky. The number is the average radiance of that mixture, expressed as a temperature.

The hottest near-pad value in the window is 345.55 K. The scaled encoding can represent up to 410.922 K. The scene stayed below that ceiling, and the file-level `saturated_pixel_count` is 0 on every frame that contains the pad.

## The 411 K ceiling

The integer stored on disk has `valid_range` `[0, 16382]` (14 bits; 16383 is `_FillValue`). After xarray's scaling the largest finite radiance is

```
L_max = 16382 × 0.001564351 + (−0.0376) = 25.589598
```

`scale_factor` is 0.001564351 and `add_offset` is −0.0376, both from the variable encoding. `load_frame` builds `L_max` from `valid_range`, `scale_factor`, and `add_offset`, then passes that one radiance through `radiance_to_bt`. The inspection prints the result: 410.9220468980713 K. A hotter source can only be reported as this value. Flight 11 did not reach it.

## Sunlight stays in the number

`esun` would be the bandpass-weighted solar irradiance at the mean Earth–Sun distance, in W m⁻² µm⁻¹. `kappa0` would be the factor that turns radiance into a reflectance. Both are NaN on this file. The attribute text says they are the reflectance conversion, and the inspection says the same: NaN means this emissive-band file does not carry them.

`radiance_to_bt` therefore inverts the observed radiance as it stands. Reflected sunlight is part of that radiance. At liftoff the solar zenith angle was 82.1°, so the Sun was about 8° above the horizon, and band 7 still contains a solar contribution. How that angle is computed is in [Pixel size, zenith, parallax, and the Sun](viewing-geometry.md).

The background in the crop stays near 298–300 K across the whole window, including the slow drift of a few tenths of a kelvin. The pad pixel moves by about 40 K on one scan and is back near 300 K on the next. That jump is a change in the pixel, on top of a background that already includes the sunlight.

## Data quality flags

`DQF` uses flag values 0–4:

| Value | Meaning in the file |
| --- | --- |
| 0 | `good_pixel_qf` |
| 1 | `conditionally_usable_pixel_qf` |
| 2 | `out_of_range_pixel_qf` |
| 3 | `no_value_pixel_qf` |
| 4 | `focal_plane_temperature_threshold_exceeded_qf` |

The script counts those flags on the finite brightness temperatures inside the crop. It does not blank the image with them. Every covering frame in this window counted `{0: 7848}`: 7848 crop pixels, all good. The full mesoscale array is 500×500; 7848 is the crop, not the sector.
