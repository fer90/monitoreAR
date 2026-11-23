# How the case study calculates its numbers

These notes walk through the Starship Flight 11 case study and the poller that stores its score. The first five cover `scripts/m1_starship_flight11.py` and [`outputs/m1_inspection.txt`](../outputs/m1_inspection.txt). The sixth covers `scripts/m2_signal.py` and [`outputs/m2_signal.txt`](../outputs/m2_signal.txt). The seventh covers `scripts/m3_score.py` and [`outputs/m3_score.txt`](../outputs/m3_score.txt). The eighth covers `scripts/m4_poll.py` and [`outputs/m4_poll.txt`](../outputs/m4_poll.txt). The public API is a later milestone and is described in [`ARCHITECTURE.md`](../ARCHITECTURE.md).

The launch used as the example is Starship Flight 11, Starbase Pad 1, 13 October 2025, 23:23:45 UTC. The imagery is NOAA GOES-19 ABI Level 1b band 7.

Read them in this order:

1. [From the bucket to one frame](data-flow.md) — filenames, the two mesoscale sectors, and what the script keeps.
2. [Radiance to brightness temperature](radiance-to-brightness-temperature.md) — the Planck coefficients in the file, the 411 K ceiling, and what 339 K means.
3. [Fixed-grid geometry](fixed-grid-geometry.md) — scan angles to latitude and longitude, and the checks on that conversion.
4. [Pixel size, zenith, parallax, and the Sun](viewing-geometry.md) — the 2.3 km pixel, the 39° view, the 8 km parallax scale, and the solar zenith at liftoff.
5. [How the figure and the report are chosen](figure-and-report.md) — the three panels, the colour scale, the white cross, and the two maxima.
6. [How the signal is compared](signal-comparison.md) — the quiet evening, the fire, the 45 K jump, and why none of that is a threshold yet.
7. [How a watch site is scored](scoring.md) — the 20 K rule, the two Flight 11 candidates, and the fire-product hold.
8. [How the poller stores a scan](polling.md) — the SQLite rows, the listing cache, and the live coverage gap.

Every temperature and distance quoted here is printed by those scripts. Re-run them from the repository root to regenerate the reports:

```bash
.venv/bin/python scripts/m1_starship_flight11.py
.venv/bin/python scripts/m2_signal.py
.venv/bin/python scripts/m3_score.py
.venv/bin/python scripts/m4_poll.py
```
