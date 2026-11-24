# Running monitoreAR

The measurement code is the `monitorear` package. The files in `scripts/` call the same entry points, so the milestone commands still work. Install once from the repository root:

```bash
python3 -m venv .venv
.venv/bin/pip install -e .
```

On Debian or Ubuntu, `python3 -m venv` needs the `python3-venv` package. Jupyter is optional and is not required to run the detector:

```bash
.venv/bin/pip install -e ".[notebooks]"
```

`data/` holds the NetCDF cache and `data/m4/monitorear.sqlite`. Both stay gitignored. `.venv/` stays gitignored too.

## Offline checks

These do not download anything. From the repository root:

```bash
.venv/bin/python -c "from monitorear.geolocate import self_test_geolocation; self_test_geolocation()"
.venv/bin/python -c "from monitorear.score import self_test_rule; self_test_rule()"
.venv/bin/python -c "from monitorear.store import self_test_store; self_test_store()"
.venv/bin/python -c "from monitorear.api import self_test_api; self_test_api()"
```

The geolocation check is a pad round-trip on the fixed grid. The rule check locks the trailing median and the hold rule, including the frozen 20 K cut. The store check uses an in-memory database. The API check rejects a crop name that escapes the crop directory.

## Case study

These commands read the public NOAA bucket when a file is not already under `data/`. A second run reuses files that are already there.

```bash
.venv/bin/python scripts/m1_starship_flight11.py
.venv/bin/python scripts/m2_signal.py
.venv/bin/python scripts/m3_score.py
.venv/bin/python scripts/m4_poll.py
.venv/bin/python scripts/m5_api.py
```

The same entry points are `python -m monitorear.flight11`, `signal`, `scoring`, `poll`, and `api`. After `pip install -e .`, the console scripts are `monitorear-flight11`, `monitorear-signal`, `monitorear-scoring`, `monitorear-poll`, and `monitorear-api`.

`scripts/m4_poll.py` replays Flight 11 hour 23, then takes two live passes 60 seconds apart, and exits. `LIVE_CYCLES` is 2. It is not a forever loop. `scripts/m5_api.py` listens on http://127.0.0.1:8000/ until interrupted. It reads `data/m4/monitorear.sqlite` and does not recompute the score. The page it serves is `monitorear/map.html`. `scripts/m5_map.html` is that same file.

## Notebooks

`notebooks/01_flight11_scene.ipynb`, `02_signal_comparison.ipynb`, and `03_scoring.ipynb` narrate Milestones 1–3 and run the offline self-tests. They read reports that are already in `outputs/`. They do not list the bucket and they do not choose a new threshold. Open them from the repository root.

## License

The code is under the MIT license in [`LICENSE`](../LICENSE). Copyright (c) 2026 Fernando Pereira. A copy of the program keeps the copyright notice and the license text. The NOAA data files are not covered by that license.
