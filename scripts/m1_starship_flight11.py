#!/usr/bin/env python3
"""Milestone 1: GOES-19 ABI band 7 around Starship Flight 11.

Downloads public L1b mesoscale radiance from the NOAA Open Data bucket,
converts it to brightness temperature, and plots the scene with the
Starbase pad on the ABI fixed grid.

Run from the repository root:

    .venv/bin/python scripts/m1_starship_flight11.py

What this script deliberately does not do: detect launches, compare a
control day, or poll for new files. Those are later milestones.
"""

from __future__ import annotations

import math
from datetime import datetime, timedelta, timezone
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import s3fs
import xarray as xr

# Liftoff 6:23:45 pm CDT = 23:23:45 UTC on 13 Oct 2025
# (NASASpaceflight.com flight timeline). File scan times below are the
# instrument's clock, not this quoted liftoff.
LIFTOFF = datetime(2025, 10, 13, 23, 23, 45, tzinfo=timezone.utc)

# Orbital Launch Pad 1 (OLM-A). Launch Library 2 lists 25.9962 N,
# 97.1544 W, which matches Wikidata 25°59'46.2"N, 97°9'16.0"W to ~10 m.
# A band-7 pixel here is a few kilometres wide, so that uncertainty
# does not move the marker by a pixel.
PAD_LAT = 25.9962
PAD_LON = -97.1544

# Scene around the pad and the first minutes of flight over the Gulf.
CROP_LAT = (25.35, 26.85)
CROP_LON = (-97.85, -95.70)

BUCKET = "noaa-goes19"
# Day 286 is 13 Oct 2025. Mesoscale L1b, both sectors, band 7 only.
PREFIX = f"{BUCKET}/ABI-L1b-RadM/2025/286/23/"
WINDOW_START = datetime(2025, 10, 13, 23, 18, tzinfo=timezone.utc)
WINDOW_END = datetime(2025, 10, 13, 23, 36, tzinfo=timezone.utc)

ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / "data" / "m1"
OUT_DIR = ROOT / "outputs"
REPORT_PATH = OUT_DIR / "m1_inspection.txt"
FIGURE_PATH = OUT_DIR / "m1_starship_flight11_band7.png"

# GRS80 values used only by the offline geolocation self-test.
# Real files carry their own semi-axes and perspective height.
_TEST_R_EQ = 6378137.0
_TEST_R_POL = 6356752.31414
_TEST_H = 35786023.0
_TEST_LON0 = -75.2


def haversine_km(lat1, lon1, lat2, lon2):
    """Great-circle distance in kilometres. Inputs are degrees."""
    r = 6371.0088
    p1 = np.radians(lat1)
    p2 = np.radians(lat2)
    dphi = np.radians(lat2 - lat1)
    dlmb = np.radians(lon2 - lon1)
    a = np.sin(dphi / 2) ** 2 + np.cos(p1) * np.cos(p2) * np.sin(dlmb / 2) ** 2
    return 2 * r * np.arcsin(np.sqrt(np.clip(a, 0, 1)))


def latlon_to_scan_angles(lat_deg, lon_deg, r_eq, r_pol, h_persp, lon0_deg):
    """Geodetic lat/lon to ABI fixed-grid scan angles (radians).

    PUG L1b Volume 3, section 5.1.2.8. The fixed grid is the view from an
    ideal satellite. ``h_persp`` is ``perspective_point_height`` (height
    above the ellipsoid). The satellite's distance from Earth's centre is
    that height plus the equatorial radius.
    """
    lat = np.radians(lat_deg)
    lon = np.radians(lon_deg)
    lam0 = np.radians(lon0_deg)
    h_sat = h_persp + r_eq
    e2 = 1.0 - (r_pol / r_eq) ** 2

    # Geocentric latitude and the distance from Earth's centre to the
    # ellipsoid at that latitude.
    phi_c = np.arctan((r_pol**2 / r_eq**2) * np.tan(lat))
    r_c = r_pol / np.sqrt(1.0 - e2 * np.cos(phi_c) ** 2)

    s_x = h_sat - r_c * np.cos(phi_c) * np.cos(lon - lam0)
    s_y = -r_c * np.cos(phi_c) * np.sin(lon - lam0)
    s_z = r_c * np.sin(phi_c)

    # Discriminant of the line-of-sight / ellipsoid intersection.
    # Negative means the point is past the limb.
    visible = h_sat * (h_sat - s_x) >= (s_y**2 + (r_eq**2 / r_pol**2) * s_z**2)
    if not np.all(visible):
        raise ValueError("point is beyond the satellite limb")

    # sweep_angle_axis = x, which is how ABI defines the fixed grid.
    y = np.arctan(s_z / s_x)
    x = np.arcsin(-s_y / np.sqrt(s_x**2 + s_y**2 + s_z**2))
    return x, y


def scan_angles_to_latlon(x_rad, y_rad, r_eq, r_pol, h_persp, lon0_deg):
    """ABI fixed-grid scan angles (radians) to geodetic lat/lon.

    Inverse of ``latlon_to_scan_angles``. Pixels that look past the Earth
    come back as NaN. The intersection is with the GRS80 ellipsoid, so a
    hot source at altitude is placed where the line of sight hits the
    ellipsoid, not at the source (parallax).
    """
    x = np.asarray(x_rad, dtype=np.float64)
    y = np.asarray(y_rad, dtype=np.float64)
    h_sat = h_persp + r_eq
    lam0 = np.radians(lon0_deg)

    a = np.sin(x) ** 2 + np.cos(x) ** 2 * (
        np.cos(y) ** 2 + (r_eq**2 / r_pol**2) * np.sin(y) ** 2
    )
    b = -2.0 * h_sat * np.cos(x) * np.cos(y)
    c = h_sat**2 - r_eq**2
    disc = b**2 - 4.0 * a * c
    good = disc >= 0

    rs = np.full(np.broadcast(x, y).shape, np.nan, dtype=np.float64)
    rs[good] = (-b[good] - np.sqrt(disc[good])) / (2.0 * a[good])

    s_x = rs * np.cos(x) * np.cos(y)
    s_y = -rs * np.sin(x)
    s_z = rs * np.cos(x) * np.sin(y)

    lat = np.full(rs.shape, np.nan, dtype=np.float64)
    lon = np.full(rs.shape, np.nan, dtype=np.float64)
    horiz = np.sqrt((h_sat - s_x) ** 2 + s_y**2)
    ok = good & (horiz > 0)
    lat[ok] = np.degrees(
        np.arctan((r_eq**2 / r_pol**2) * (s_z[ok] / horiz[ok]))
    )
    lon[ok] = np.degrees(lam0 - np.arctan2(s_y[ok], h_sat - s_x[ok]))
    return lat, lon


def satellite_zenith_deg(lat_deg, lon_deg, r_eq, r_pol, h_persp, lon0_deg):
    """Approximate geocentric satellite zenith angle at one surface point.

    Enough to explain why the pixel footprint grows and why an elevated
    plume is displaced. Not a navigation product.
    """
    lat = math.radians(lat_deg)
    # Geocentric latitude, consistent with the fixed-grid equations.
    phi_c = math.atan((r_pol**2 / r_eq**2) * math.tan(lat))
    e2 = 1.0 - (r_pol / r_eq) ** 2
    r_c = r_pol / math.sqrt(1.0 - e2 * math.cos(phi_c) ** 2)
    h_sat = h_persp + r_eq
    cos_beta = math.cos(phi_c) * math.cos(math.radians(lon_deg - lon0_deg))
    cos_beta = min(1.0, max(-1.0, cos_beta))
    beta = math.acos(cos_beta)
    # Elevation of the satellite above the local horizontal, spherical.
    elev = math.atan2(math.cos(beta) - r_c / h_sat, math.sin(beta))
    return 90.0 - math.degrees(elev)


def solar_zenith_deg(lat_deg, lon_deg, when):
    """Solar zenith angle in degrees. NOAA general solar-position series.

    Accurate to about a degree, which is enough to know whether band 7
    still contains reflected sunlight. Longitude is degrees east.
    """
    when = when.astimezone(timezone.utc)
    doy = when.timetuple().tm_yday
    hour = when.hour + when.minute / 60 + when.second / 3600
    gamma = 2 * math.pi / 365 * (doy - 1 + (hour - 12) / 24)
    eqtime = 229.18 * (
        0.000075
        + 0.001868 * math.cos(gamma)
        - 0.032077 * math.sin(gamma)
        - 0.014615 * math.cos(2 * gamma)
        - 0.040849 * math.sin(2 * gamma)
    )
    decl = (
        0.006918
        - 0.399912 * math.cos(gamma)
        + 0.070257 * math.sin(gamma)
        - 0.006758 * math.cos(2 * gamma)
        + 0.000907 * math.sin(2 * gamma)
        - 0.002697 * math.cos(3 * gamma)
        + 0.00148 * math.sin(3 * gamma)
    )
    time_offset = eqtime + 4 * lon_deg  # minutes; timezone is UTC
    tst = hour * 60 + time_offset
    ha = math.radians(tst / 4 - 180)
    lat = math.radians(lat_deg)
    cos_z = math.sin(lat) * math.sin(decl) + math.cos(lat) * math.cos(decl) * math.cos(ha)
    cos_z = min(1.0, max(-1.0, cos_z))
    return math.degrees(math.acos(cos_z))


def initial_bearing_deg(lat1, lon1, lat2, lon2):
    """Initial great-circle bearing from point 1 to point 2, degrees from north."""
    p1, p2 = math.radians(lat1), math.radians(lat2)
    dl = math.radians(lon2 - lon1)
    y = math.sin(dl) * math.cos(p2)
    x = math.cos(p1) * math.sin(p2) - math.sin(p1) * math.cos(p2) * math.cos(dl)
    return (math.degrees(math.atan2(y, x)) + 360.0) % 360.0


def radiance_to_bt(rad, fk1, fk2, bc1, bc2):
    """ABI L1b radiance to brightness temperature (K).

    T = [fk2 / ln(fk1 / L + 1) - bc1] / bc2

    The coefficients live on the radiance variable and already match that
    file's radiance units. This is a Planck inversion of the observed
    radiance. It does not remove reflected sunlight, which band 7
    (3.9 µm) still sees while the sun is up.
    """
    rad = np.asarray(rad, dtype=np.float64)
    bt = np.full(rad.shape, np.nan, dtype=np.float64)
    ok = np.isfinite(rad) & (rad > 0)
    bt[ok] = (fk2 / np.log(fk1 / rad[ok] + 1.0) - bc1) / bc2
    return bt


def parse_goes_time(token):
    """Parse the 14-digit s/e/c timestamp in an ABI filename.

    Layout is YYYYDDDHHMMSSt, with the last digit a tenth of a second.
    """
    t = token[1:] if token[:1] in "sec" else token
    t = t.replace(".nc", "")
    if len(t) != 14:
        raise ValueError(f"unexpected GOES timestamp {token!r}")
    year, doy = int(t[0:4]), int(t[4:7])
    hour, minute, second, tenth = int(t[7:9]), int(t[9:11]), int(t[11:13]), int(t[13])
    base = datetime(year, 1, 1, tzinfo=timezone.utc) + timedelta(days=doy - 1)
    return base + timedelta(
        hours=hour, minutes=minute, seconds=second, milliseconds=tenth * 100
    )


def parse_abi_name(name):
    """Split an ABI L1b object name into sector, mode, and the three times."""
    parts = name.replace(".nc", "").split("_")
    # OR, ABI-L1b-RadM1-M6C07, G19, s..., e..., c...
    product = parts[1]
    sector = "M1" if "RadM1" in product else "M2" if "RadM2" in product else "?"
    mode = "M6" if "-M6" in product else product
    return {
        "sector": sector,
        "mode": mode,
        "product": product,
        "satellite": parts[2],
        "start": parse_goes_time(parts[3]),
        "end": parse_goes_time(parts[4]),
        "created": parse_goes_time(parts[5]),
    }


def edges_1d(centers):
    """Pixel-center coordinates to cell edges. Spacing may vary slightly."""
    v = np.asarray(centers, dtype=np.float64)
    d = np.diff(v)
    out = np.empty(v.size + 1, dtype=np.float64)
    out[1:-1] = (v[:-1] + v[1:]) / 2
    out[0] = v[0] - d[0] / 2
    out[-1] = v[-1] + d[-1] / 2
    return out


def projection_constants(ds):
    proj = ds["goes_imager_projection"]
    return {
        "r_eq": float(proj.attrs["semi_major_axis"]),
        "r_pol": float(proj.attrs["semi_minor_axis"]),
        "h_persp": float(proj.attrs["perspective_point_height"]),
        "lon0": float(proj.attrs["longitude_of_projection_origin"]),
        "sweep": str(proj.attrs.get("sweep_angle_axis", "x")),
    }


def self_test_geolocation():
    """Nadir point and a pad round-trip, before any download."""
    lat, lon = scan_angles_to_latlon(
        0.0, 0.0, _TEST_R_EQ, _TEST_R_POL, _TEST_H, _TEST_LON0
    )
    if abs(float(lat)) > 1e-6 or abs(float(lon) - _TEST_LON0) > 1e-4:
        raise SystemExit(f"nadir geolocation failed: lat={lat} lon={lon}")
    x, y = latlon_to_scan_angles(
        PAD_LAT, PAD_LON, _TEST_R_EQ, _TEST_R_POL, _TEST_H, _TEST_LON0
    )
    lat2, lon2 = scan_angles_to_latlon(x, y, _TEST_R_EQ, _TEST_R_POL, _TEST_H, _TEST_LON0)
    err_m = float(haversine_km(PAD_LAT, PAD_LON, float(lat2), float(lon2)) * 1000)
    print(f"geolocation self-test: nadir ok, pad round-trip {err_m:.2f} m")
    if err_m > 200:
        raise SystemExit("fixed-grid round trip is off by more than 200 m")
    return err_m


def list_band7_keys(fs):
    keys = fs.ls(PREFIX)
    chosen = []
    for key in keys:
        name = key.split("/")[-1]
        if "C07_G19_" not in name:
            continue
        meta = parse_abi_name(name)
        if WINDOW_START <= meta["start"] <= WINDOW_END:
            chosen.append((key, meta))
    chosen.sort(key=lambda item: item[1]["start"])
    return chosen


def download(fs, key):
    dest = DATA_DIR / key.split("/")[-1]
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    fs.get(key, str(dest))
    return dest


def list_other_sector_times(fs, product_prefix):
    """Start times of band-7 files in one hour, for the cadence comparison.

    Listing only. CONUS and Full Disk are not downloaded in this milestone.
    """
    keys = fs.ls(product_prefix)
    times = []
    for key in keys:
        name = key.split("/")[-1]
        if "C07_G19_" not in name:
            continue
        start = parse_goes_time(name.split("_")[3])
        if WINDOW_START <= start <= WINDOW_END:
            times.append(start)
    return sorted(times)


def _center_residual_m(ds, consts):
    lat, lon = scan_angles_to_latlon(
        float(ds["x_image"].values),
        float(ds["y_image"].values),
        consts["r_eq"],
        consts["r_pol"],
        consts["h_persp"],
        consts["lon0"],
    )
    extent = ds["geospatial_lat_lon_extent"].attrs
    return float(
        haversine_km(
            float(lat),
            float(lon),
            float(extent["geospatial_lat_center"]),
            float(extent["geospatial_lon_center"]),
        )
        * 1000
    )


def load_frame(path, meta):
    """Open one L1b file and reduce it to the Starbase crop."""
    ds = xr.open_dataset(path)
    try:
        consts = projection_constants(ds)
        if consts["sweep"] != "x":
            raise SystemExit(f"{path.name}: unexpected sweep axis {consts['sweep']}")

        rad = ds["Rad"]
        # On current GOES-19 L1b files these are scalar variables, not
        # attributes of Rad. They are paired with radiance in
        # mW m-2 sr-1 (cm-1)-1.
        fk1 = float(ds["planck_fk1"].values)
        fk2 = float(ds["planck_fk2"].values)
        bc1 = float(ds["planck_bc1"].values)
        bc2 = float(ds["planck_bc2"].values)
        # x and y are scaled integers. xarray applies scale_factor, so
        # these values are scan angles in radians, not pixel indices.
        x = np.asarray(ds["x"].values, dtype=np.float64)
        y = np.asarray(ds["y"].values, dtype=np.float64)
        if str(ds["x"].attrs.get("units", "")).lower() not in {"rad", "radians"}:
            raise SystemExit(f"{path.name}: x units are {ds['x'].attrs.get('units')!r}, expected rad")

        x_pad, y_pad = latlon_to_scan_angles(
            PAD_LAT, PAD_LON, consts["r_eq"], consts["r_pol"], consts["h_persp"], consts["lon0"]
        )
        contains = (x.min() <= x_pad <= x.max()) and (y.min() <= y_pad <= y.max())

        bt = radiance_to_bt(rad.values, fk1, fk2, bc1, bc2)
        dqf = np.asarray(ds["DQF"].values)
        dqf_attrs = {
            "flag_values": np.asarray(ds["DQF"].attrs.get("flag_values", [])).tolist(),
            "flag_meanings": str(ds["DQF"].attrs.get("flag_meanings", "")),
            "flag_masks": np.asarray(ds["DQF"].attrs.get("flag_masks", [])).tolist()
            if "flag_masks" in ds["DQF"].attrs
            else None,
        }

        # Saturation reference: brightness temperature of the highest
        # radiance the scaled integer can represent.
        enc = rad.encoding
        valid = enc.get("valid_range", rad.attrs.get("valid_range"))
        bt_ceiling = None
        if valid is not None:
            scale = float(enc.get("scale_factor", rad.attrs.get("scale_factor", 1.0)))
            offset = float(enc.get("add_offset", rad.attrs.get("add_offset", 0.0)))
            rad_max = float(valid[1]) * scale + offset
            if rad_max > 0:
                bt_ceiling = float(radiance_to_bt(np.array([rad_max]), fk1, fk2, bc1, bc2)[0])

        file_info = {
            "attrs": {
                k: ds.attrs.get(k)
                for k in (
                    "platform_ID",
                    "instrument_type",
                    "spatial_resolution",
                    "scene_id",
                    "timeline_id",
                    "dataset_name",
                    "title",
                    "time_coverage_start",
                    "time_coverage_end",
                    "date_created",
                    "orbital_slot",
                    "production_site",
                )
            },
            "rad_units": rad.attrs.get("units"),
            "rad_standard_name": rad.attrs.get("standard_name"),
            "band_id": int(np.asarray(ds["band_id"].values).reshape(-1)[0]),
            "band_wavelength_um": float(np.asarray(ds["band_wavelength"].values).reshape(-1)[0]),
            "planck": (fk1, fk2, bc1, bc2),
            "dqf_attrs": dqf_attrs,
            "bt_ceiling": bt_ceiling,
            "consts": consts,
            "x_units": ds["x"].attrs.get("units"),
            "shape": tuple(int(v) for v in bt.shape),
            # The fixed-grid origin is an ideal longitude. The real
            # subsatellite point is a separate variable and can differ.
            "subpoint_lat": float(ds["nominal_satellite_subpoint_lat"].values),
            "subpoint_lon": float(ds["nominal_satellite_subpoint_lon"].values),
            "subpoint_height": float(ds["nominal_satellite_height"].values),
            "subpoint_height_units": ds["nominal_satellite_height"].attrs.get("units"),
            "yaw_flip_flag": float(ds["yaw_flip_flag"].values),
            "esun": float(ds["esun"].values),
            "kappa0": float(ds["kappa0"].values),
            "saturated_pixel_count": int(ds["saturated_pixel_count"].values),
            "valid_pixel_count": int(ds["valid_pixel_count"].values),
            # Check our scan-angle conversion against the lat/lon NOAA
            # already stored for the center of this sector.
            "center_residual_m": _center_residual_m(ds, consts),
            "sector_center_lat": float(ds["geospatial_lat_lon_extent"].attrs["geospatial_lat_center"]),
            "sector_center_lon": float(ds["geospatial_lat_lon_extent"].attrs["geospatial_lon_center"]),
            "sector_bounds": (
                float(ds["geospatial_lat_lon_extent"].attrs["geospatial_southbound_latitude"]),
                float(ds["geospatial_lat_lon_extent"].attrs["geospatial_northbound_latitude"]),
                float(ds["geospatial_lat_lon_extent"].attrs["geospatial_westbound_longitude"]),
                float(ds["geospatial_lat_lon_extent"].attrs["geospatial_eastbound_longitude"]),
            ),
        }

        record = {
            "path": path,
            "meta": meta,
            "contains_pad": bool(contains),
            "file_info": file_info,
            "crop": None,
        }
        if not contains:
            return record

        xx, yy = np.meshgrid(x, y)
        lat, lon = scan_angles_to_latlon(
            xx, yy, consts["r_eq"], consts["r_pol"], consts["h_persp"], consts["lon0"]
        )
        in_window = (
            (lat >= CROP_LAT[0])
            & (lat <= CROP_LAT[1])
            & (lon >= CROP_LON[0])
            & (lon <= CROP_LON[1])
        )
        rows = np.where(in_window.any(axis=1))[0]
        cols = np.where(in_window.any(axis=0))[0]
        if rows.size == 0 or cols.size == 0:
            record["contains_pad"] = False
            return record
        r0, r1 = int(rows[0]), int(rows[-1]) + 1
        c0, c1 = int(cols[0]), int(cols[-1]) + 1

        x_e = edges_1d(x[c0:c1])
        y_e = edges_1d(y[r0:r1])
        xe, ye = np.meshgrid(x_e, y_e)
        lat_e, lon_e = scan_angles_to_latlon(
            xe, ye, consts["r_eq"], consts["r_pol"], consts["h_persp"], consts["lon0"]
        )

        iy = int(np.argmin(np.abs(y - y_pad)))
        ix = int(np.argmin(np.abs(x - x_pad)))
        pad_lat_px = float(lat[iy, ix])
        pad_lon_px = float(lon[iy, ix])

        # Footprint of the pixel that contains the pad, from its neighbours.
        def _km(i0, j0, i1, j1):
            return float(haversine_km(lat[i0, j0], lon[i0, j0], lat[i1, j1], lon[i1, j1]))

        pix_dx = _km(iy, ix, iy, min(ix + 1, lon.shape[1] - 1))
        pix_dy = _km(iy, ix, min(iy + 1, lat.shape[0] - 1), ix)

        crop_bt = bt[r0:r1, c0:c1]
        crop_lat = lat[r0:r1, c0:c1]
        crop_lon = lon[r0:r1, c0:c1]
        crop_dqf = dqf[r0:r1, c0:c1]
        dist = haversine_km(crop_lat, crop_lon, PAD_LAT, PAD_LON)
        near = dist <= 20.0
        finite_near = near & np.isfinite(crop_bt)
        if finite_near.any():
            masked = np.where(finite_near, crop_bt, -np.inf)
            j, i = np.unravel_index(int(np.argmax(masked)), crop_bt.shape)
            max_bt = float(crop_bt[j, i])
            max_dist = float(dist[j, i])
            max_lat = float(crop_lat[j, i])
            max_lon = float(crop_lon[j, i])
        else:
            max_bt = max_dist = max_lat = max_lon = float("nan")

        flat = np.nanargmax(crop_bt)
        j_all, i_all = np.unravel_index(int(flat), crop_bt.shape)

        values, counts = np.unique(crop_dqf[np.isfinite(crop_bt)], return_counts=True)
        record["crop"] = {
            "bt": crop_bt,
            "lat_edges": lat_e,
            "lon_edges": lon_e,
            "dqf_counts": {int(v): int(c) for v, c in zip(values, counts)},
            "max_bt_20km": max_bt,
            "max_dist_km": max_dist,
            "max_lat": max_lat,
            "max_lon": max_lon,
            "crop_max_bt": float(crop_bt[j_all, i_all]),
            "crop_max_dist_km": float(dist[j_all, i_all]),
            "crop_max_lat": float(crop_lat[j_all, i_all]),
            "crop_max_lon": float(crop_lon[j_all, i_all]),
            "median_bt": float(np.nanmedian(crop_bt)),
            "p02_bt": float(np.nanpercentile(crop_bt, 2)),
            "p98_bt": float(np.nanpercentile(crop_bt, 98)),
            "pad_pixel_lat": pad_lat_px,
            "pad_pixel_lon": pad_lon_px,
            "pad_offset_km": float(haversine_km(PAD_LAT, PAD_LON, pad_lat_px, pad_lon_px)),
            "pad_pixel_bt": float(bt[iy, ix]),
            "pixel_dx_km": pix_dx,
            "pixel_dy_km": pix_dy,
            "x_pad": float(x_pad),
            "y_pad": float(y_pad),
        }
        return record
    finally:
        ds.close()


def fmt_time(when):
    return when.strftime("%Y-%m-%d %H:%M:%S.") + f"{when.microsecond // 100000:01d}"


def write_report(lines):
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    REPORT_PATH.write_text("\n".join(lines) + "\n")


def plot_frames(frames, labels, vmin, vmax):
    fig, axes = plt.subplots(
        1, 3, figsize=(13.6, 5.4), sharex=True, sharey=True, constrained_layout=True
    )
    mesh = None
    for ax, frame, label in zip(axes, frames, labels):
        crop = frame["crop"]
        mesh = ax.pcolormesh(
            crop["lon_edges"],
            crop["lat_edges"],
            crop["bt"],
            cmap="inferno",
            vmin=vmin,
            vmax=vmax,
            shading="flat",
            rasterized=True,
        )
        ax.plot(PAD_LON, PAD_LAT, marker="o", ms=8, mfc="none", mec="cyan", mew=1.5, label="Pad 1", zorder=4)
        # A warm-background pixel is not a plume. Mark the hot pixel only
        # when it is far above the scene and not the pad pixel itself.
        anomalous = np.isfinite(crop["max_bt_20km"]) and crop["max_bt_20km"] > crop["median_bt"] + 10
        if anomalous and crop["max_dist_km"] > 3:
            ax.plot(
                crop["max_lon"],
                crop["max_lat"],
                marker="+",
                ms=11,
                mec="white",
                mew=1.4,
                label="hottest pixel within 20 km",
                zorder=5,
            )
        start = frame["meta"]["start"]
        rel = (start - LIFTOFF).total_seconds()
        ax.set_title(
            f"{label}\n{frame['meta']['sector']}  {start:%H:%M:%S} UTC  ({rel:+.0f} s)\n"
            f"pad pixel {crop['pad_pixel_bt']:.1f} K"
            f"   max within 20 km {crop['max_bt_20km']:.1f} K"
            f" ({crop['max_dist_km']:.1f} km)",
            fontsize=9,
        )
        ax.set_xlim(CROP_LON)
        ax.set_ylim(CROP_LAT)
        ax.set_aspect(1 / math.cos(math.radians(PAD_LAT)))
        ax.set_xlabel("Longitude")
        ax.grid(True, color="white", alpha=0.25, linewidth=0.4)
    axes[0].set_ylabel("Latitude")
    handles, legend_labels = [], []
    for ax in axes:
        found_handles, found_labels = ax.get_legend_handles_labels()
        for handle, text in zip(found_handles, found_labels):
            if text not in legend_labels:
                handles.append(handle)
                legend_labels.append(text)
    axes[0].legend(handles, legend_labels, loc="lower left", fontsize=7, framealpha=0.85)
    cbar = fig.colorbar(mesh, ax=axes, shrink=0.86, pad=0.02)
    cbar.set_label("Brightness temperature (K)")
    fig.suptitle(
        "GOES-19 ABI band 7 (3.9 µm), Starship Flight 11 — 13 Oct 2025",
        fontsize=12,
    )
    fig.text(
        0.01,
        0.005,
        "NOAA GOES-19 ABI L1b radiance, converted with the Planck coefficients stored in each file. "
        "Color scale runs from the pre-launch 2nd percentile to 350 K, so a ~300 K surface and a ~340 K plume pixel are not the same color. "
        "The pad is geolocated on the fixed grid (GRS80 ellipsoid). No parallax correction.",
        fontsize=6.5,
        ha="left",
        va="bottom",
    )
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIGURE_PATH, dpi=140)
    plt.close(fig)


def main():
    err_m = self_test_geolocation()
    DATA_DIR.mkdir(parents=True, exist_ok=True)

    fs = s3fs.S3FileSystem(anon=True, client_kwargs={"region_name": "us-east-1"})
    chosen = list_band7_keys(fs)
    if not chosen:
        raise SystemExit(f"no band-7 mesoscale files under {PREFIX}")
    print(f"downloading {len(chosen)} mesoscale band-7 files")

    frames = []
    for key, meta in chosen:
        path = download(fs, key)
        print(f"  {path.name}")
        frames.append(load_frame(path, meta))

    covering = [f for f in frames if f["contains_pad"]]
    missing = [f for f in frames if not f["contains_pad"]]
    if not covering:
        raise SystemExit(
            "neither mesoscale sector contained Starbase in this window. "
            "The sectors were aimed somewhere else; CONUS would be the next place to look."
        )

    info = covering[0]["file_info"]
    consts = info["consts"]
    zenith = satellite_zenith_deg(
        PAD_LAT, PAD_LON, consts["r_eq"], consts["r_pol"], consts["h_persp"], consts["lon0"]
    )
    sun_z = solar_zenith_deg(PAD_LAT, PAD_LON, LIFTOFF)
    bearing = initial_bearing_deg(0.0, consts["lon0"], PAD_LAT, PAD_LON)
    parallax_per_10km = 10.0 * math.tan(math.radians(zenith))

    # Same-hour CONUS and Full Disk start times, so the cadence comparison
    # is measured rather than remembered.
    conus_times = list_other_sector_times(fs, f"{BUCKET}/ABI-L1b-RadC/2025/286/23/")
    fulldisk_times = list_other_sector_times(fs, f"{BUCKET}/ABI-L1b-RadF/2025/286/23/")

    pre = [f for f in covering if f["meta"]["end"] <= LIFTOFF]
    if not pre:
        raise SystemExit("no pre-liftoff frame contained the pad")
    before = pre[-1]
    after = [
        f
        for f in covering
        if f["meta"]["start"] >= LIFTOFF and f["meta"]["start"] <= LIFTOFF + timedelta(minutes=4)
    ]
    if not after:
        raise SystemExit("no post-liftoff frame in the first 4 minutes contained the pad")
    # The pad pixel and the neighbourhood maximum are not the same scan:
    # the first post-liftoff scan heats the pad pixel, and the next scan's
    # hottest pixel has moved off it.
    on_pad = max(after, key=lambda f: f["crop"]["pad_pixel_bt"])
    displaced = max(after, key=lambda f: f["crop"]["max_bt_20km"])
    # 350 K is above the observed plume (~346 K) and well above the ~300 K
    # background, so the hot pixel is not clipped to the same color as the land.
    plot_frames(
        [before, on_pad, displaced],
        ["before liftoff", "pad pixel heats", "hot pixel has moved"],
        vmin=before["crop"]["p02_bt"],
        vmax=350.0,
    )

    lines = []
    add = lines.append
    add("Milestone 1 inspection — Starship Flight 11, GOES-19 ABI band 7")
    add("")
    add(f"Liftoff reference: {LIFTOFF.isoformat()}  (quoted; not read from the satellite)")
    add(f"Pad: {PAD_LAT:.4f} N, {PAD_LON:.4f} W")
    add(f"Geolocation self-test round-trip at the pad: {err_m:.2f} m")
    add("")
    add("Instrument file (first frame that contains the pad)")
    add(f"  file: {covering[0]['path'].name}")
    for key, value in info["attrs"].items():
        add(f"  {key}: {value}")
    add(f"  band_id: {info['band_id']}")
    add(f"  band_wavelength_um: {info['band_wavelength_um']}")
    add(f"  radiance units: {info['rad_units']}")
    add(f"  radiance standard_name: {info['rad_standard_name']}")
    add(f"  planck fk1, fk2, bc1, bc2: {info['planck']}")
    add(f"  esun (W m-2 um-1): {info['esun']}")
    add(f"  kappa0: {info['kappa0']}")
    add("  esun and kappa0 are the reflectance conversion. NaN means this emissive-band file does not carry them.")
    add(f"  BT at the top of the scaled radiance range: {info['bt_ceiling']}")
    add(f"  array shape (y, x): {info['shape']}")
    add(f"  x units: {info['x_units']}")
    add(f"  fixed-grid lon0 (deg): {consts['lon0']}")
    add(f"  actual subsatellite lat, lon: {info['subpoint_lat']}, {info['subpoint_lon']}")
    add(f"  nominal_satellite_height: {info['subpoint_height']} {info['subpoint_height_units']}")
    add(f"  yaw_flip_flag: {info['yaw_flip_flag']}")
    add(f"  perspective_point_height (m): {consts['h_persp']}")
    add(f"  sweep_angle_axis: {consts['sweep']}")
    add(f"  DQF flag_values: {info['dqf_attrs']['flag_values']}")
    add(f"  DQF flag_meanings: {info['dqf_attrs']['flag_meanings']}")
    add("")
    crop0 = covering[0]["crop"]
    add("At the pad, from the fixed grid")
    add(f"  nearest pixel center: {crop0['pad_pixel_lat']:.4f} N, {crop0['pad_pixel_lon']:.4f} W")
    add(f"  pad-to-pixel-center distance: {crop0['pad_offset_km']:.3f} km")
    add(f"  scan-angle conversion vs file geospatial center: {info['center_residual_m']:.3f} m")
    add(f"  neighbouring-pixel spacing: {crop0['pixel_dx_km']:.2f} km (east-west), {crop0['pixel_dy_km']:.2f} km (north-south)")
    add(f"  satellite zenith angle (approx): {zenith:.1f} deg")
    add(f"  parallax: about {parallax_per_10km:.1f} km of apparent shift per 10 km of source altitude")
    add(f"  that shift is away from the subsatellite point, bearing {bearing:.0f} deg from north")
    add(f"  solar zenith angle at liftoff (NOAA approximation): {sun_z:.1f} deg")
    add("  band 7 BT includes reflected solar radiance while the solar zenith is below ~90 deg")
    add("")
    add("Sectors in the 23:18–23:36 UTC window")
    add(f"  files listed: {len(frames)}")
    add(f"  files whose fixed grid contains the pad: {len(covering)}")
    if missing:
        sectors = sorted({f["meta"]["sector"] for f in missing})
        add(f"  files that do not contain the pad: {len(missing)} (sectors {', '.join(sectors)})")
        # One example per sector: mesoscale domains are tasked and can sit far from the pad.
        seen = set()
        for f in missing:
            sector = f["meta"]["sector"]
            if sector in seen:
                continue
            seen.add(sector)
            s, n, w, e = f["file_info"]["sector_bounds"]
            add(
                f"    {sector} example center {f['file_info']['sector_center_lat']:.2f} N, "
                f"{f['file_info']['sector_center_lon']:.2f} E; "
                f"bounds lat {s:.2f} to {n:.2f}, lon {w:.2f} to {e:.2f} "
                f"({f['path'].name})"
            )
    durations = [(f["meta"]["end"] - f["meta"]["start"]).total_seconds() for f in covering]
    lags = [(f["meta"]["created"] - f["meta"]["end"]).total_seconds() for f in covering]
    starts = [f["meta"]["start"] for f in covering]
    gaps = [(b - a).total_seconds() for a, b in zip(starts, starts[1:])]
    add(f"  scan duration (end-start): min {min(durations):.1f} s, max {max(durations):.1f} s")
    add(f"  gap between successive scans that see the pad: min {min(gaps):.1f} s, median {float(np.median(gaps)):.1f} s, max {max(gaps):.1f} s")
    add(f"  filename creation minus scan end: min {min(lags):.1f} s, median {float(np.median(lags)):.1f} s, max {max(lags):.1f} s")
    add("  The creation stamp is ground-system product time, not the moment the object appeared on AWS.")
    add(f"  CONUS band-7 files in this window (listed, not opened): {len(conus_times)}")
    if conus_times:
        add("    " + ", ".join(t.strftime("%H:%M:%S") for t in conus_times))
    add(f"  Full Disk band-7 files in this window (listed, not opened): {len(fulldisk_times)}")
    if fulldisk_times:
        add("    " + ", ".join(t.strftime("%H:%M:%S") for t in fulldisk_times))
    add("")
    add("Per-frame brightness temperature for frames that contain the pad")
    add("  start_utc             sector  rel_s  pad_K  max20_K  dist_km  median_K  p98_K  saturated_in_file  dqf_counts")
    for f in covering:
        c = f["crop"]
        rel = (f["meta"]["start"] - LIFTOFF).total_seconds()
        add(
            f"  {fmt_time(f['meta']['start'])}  {f['meta']['sector']}  "
            f"{rel:7.1f}  {c['pad_pixel_bt']:6.2f}  {c['max_bt_20km']:7.2f}  {c['max_dist_km']:6.2f}  "
            f"{c['median_bt']:7.2f}  {c['p98_bt']:6.2f}  "
            f"{f['file_info']['saturated_pixel_count']:5d}  {c['dqf_counts']}"
        )
    add("")
    add("Frames whose crop maximum exceeds 305 K")
    hot_frames = [f for f in covering if f["crop"]["crop_max_bt"] > 305]
    if not hot_frames:
        add("  none")
    for f in hot_frames:
        c = f["crop"]
        rel = (f["meta"]["start"] - LIFTOFF).total_seconds()
        add(
            f"  {fmt_time(f['meta']['start'])}  {rel:+.0f} s  "
            f"{c['crop_max_bt']:.1f} K at {c['crop_max_dist_km']:.1f} km from the pad "
            f"({c['crop_max_lat']:.3f} N, {c['crop_max_lon']:.3f} W)"
        )
    add("  Background in this crop sits near 298–300 K. After these scans the crop maximum stays near that background.")
    add("")
    add("Figure frames")
    for label, f in ("before", before), ("pad_pixel_heats", on_pad), ("hot_pixel_moved", displaced):
        add(f"  {label}: {f['path'].name}")
    add(f"  color scale: {before['crop']['p02_bt']:.2f} K to 350 K")
    add(f"  wrote {FIGURE_PATH.relative_to(ROOT)}")
    add("")
    ceiling = info["bt_ceiling"]
    peak = max(f["crop"]["max_bt_20km"] for f in covering)
    if ceiling is not None and peak > ceiling - 0.5:
        add(f"Saturation: a near-pad pixel reached the scaled radiance ceiling ({ceiling:.1f} K).")
    else:
        add(f"Saturation: near-pad maximum {peak:.1f} K stayed below the scaled ceiling ({ceiling}).")
    add("The reported temperature is the pixel average, not the flame temperature.")

    write_report(lines)
    print("\n".join(lines))


if __name__ == "__main__":
    main()
