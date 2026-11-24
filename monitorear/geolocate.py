"""Fixed-grid geolocation for GOES-R ABI.

Scan angles in the file are radians. The conversion places the line of
sight on the GRS80 ellipsoid. A hot source above the surface is left
where that line of sight meets the ellipsoid (parallax).
"""

from __future__ import annotations

import math
from datetime import datetime, timezone

import numpy as np

# Liftoff 6:23:45 pm CDT = 23:23:45 UTC on 13 Oct 2025
# (NASASpaceflight.com flight timeline). File scan times are the
# instrument clock, not this quoted liftoff.
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
