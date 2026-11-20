# How Milestone 1 calculates its numbers

These notes walk through `scripts/m1_starship_flight11.py` and the measurements in [`outputs/m1_inspection.txt`](../outputs/m1_inspection.txt). They cover the Starship Flight 11 case study only. The detector, the control day, and the public API are later milestones and are described in [`ARCHITECTURE.md`](../ARCHITECTURE.md).

The launch used as the example is Starship Flight 11, Starbase Pad 1, 13 October 2025, 23:23:45 UTC. The imagery is NOAA GOES-19 ABI Level 1b band 7.

Read them in this order:

1. [From the bucket to one frame](data-flow.md) — filenames, the two mesoscale sectors, and what the script keeps.
2. [Radiance to brightness temperature](radiance-to-brightness-temperature.md) — the Planck coefficients in the file, the 411 K ceiling, and what 339 K means.
3. [Fixed-grid geometry](fixed-grid-geometry.md) — scan angles to latitude and longitude, and the checks on that conversion.
4. [Pixel size, zenith, parallax, and the Sun](viewing-geometry.md) — the 2.3 km pixel, the 39° view, the 8 km parallax scale, and the solar zenith at liftoff.
5. [How the figure and the report are chosen](figure-and-report.md) — the three panels, the colour scale, the white cross, and the two maxima.

Every temperature and distance quoted here is printed by that script. Re-run it from the repository root to regenerate the report:

```bash
.venv/bin/python scripts/m1_starship_flight11.py
```
