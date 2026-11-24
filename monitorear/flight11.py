"""Flight 11 band-7 scene.

Run from the repository root:

    .venv/bin/python -m monitorear.flight11
"""

from __future__ import annotations

import math
from datetime import datetime, timedelta, timezone

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import s3fs
import xarray as xr

from monitorear.calibrate import radiance_to_bt
from monitorear.geolocate import (
    CROP_LAT,
    CROP_LON,
    LIFTOFF,
    PAD_LAT,
    PAD_LON,
    edges_1d,
    haversine_km,
    initial_bearing_deg,
    latlon_to_scan_angles,
    projection_constants,
    satellite_zenith_deg,
    scan_angles_to_latlon,
    self_test_geolocation,
    solar_zenith_deg,
)
from monitorear.goes_io import BUCKET, ROOT, parse_abi_name, parse_goes_time

PREFIX = f"{BUCKET}/ABI-L1b-RadM/2025/286/23/"
WINDOW_START = datetime(2025, 10, 13, 23, 18, tzinfo=timezone.utc)
WINDOW_END = datetime(2025, 10, 13, 23, 36, tzinfo=timezone.utc)

DATA_DIR = ROOT / "data" / "m1"
OUT_DIR = ROOT / "outputs"
REPORT_PATH = OUT_DIR / "m1_inspection.txt"
FIGURE_PATH = OUT_DIR / "m1_starship_flight11_band7.png"


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
