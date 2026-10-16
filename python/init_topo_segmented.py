#!/usr/bin/env python3
"""Script entry point kept for the documented command `python python/init_topo_segmented.py ...`;
the code is in the package module slevolve.init_topo_segmented (console script: see docs/src/options.md)."""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from slevolve.init_topo_segmented import main  # noqa: E402

if __name__ == "__main__":
    main()
