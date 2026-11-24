#!/usr/bin/env python3
"""Run the Flight 11 signal comparison.

The measurement code lives in monitorear.signal. From the repository root:

    .venv/bin/python scripts/m2_signal.py
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from monitorear.signal import main

if __name__ == "__main__":
    main()
