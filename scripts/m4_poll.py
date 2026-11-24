#!/usr/bin/env python3
"""Run the mesoscale poller.

The loop lives in monitorear.poll. It takes two passes and exits.
From the repository root:

    .venv/bin/python scripts/m4_poll.py
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from monitorear.poll import main

if __name__ == "__main__":
    main()
