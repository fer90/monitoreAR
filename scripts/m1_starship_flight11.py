#!/usr/bin/env python3
"""Run the Flight 11 band-7 scene.

The measurement code lives in monitorear.flight11. From the repository root:

    .venv/bin/python scripts/m1_starship_flight11.py
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from monitorear.flight11 import main

if __name__ == "__main__":
    main()
