#!/usr/bin/env python3
"""Script entry point kept for the documented command `python python/run_outer.py ...`;
the code is in the package module slevolve.run_outer (console script: see docs/src/options.md)."""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from slevolve.run_outer import main  # noqa: E402

if __name__ == "__main__":
    main()
