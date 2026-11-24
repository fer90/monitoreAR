#!/usr/bin/env python3
"""Serve the stored candidates.

The application lives in monitorear.api. From the repository root:

    .venv/bin/python scripts/m5_api.py
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from monitorear.api import main

if __name__ == "__main__":
    main()
