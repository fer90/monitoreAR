#!/usr/bin/env python3
"""Run the Flight 11 watch-site score.

The measurement code lives in monitorear.scoring. From the repository root:

    .venv/bin/python scripts/m3_score.py
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from monitorear.scoring import main

if __name__ == "__main__":
    main()
