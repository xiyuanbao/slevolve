#!/usr/bin/env python3
"""Script entry point for slevolve.tools.write_lovenum (kept so `python python/tools/write_lovenum.py ...` works)."""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from slevolve.tools import write_lovenum as _m  # noqa: E402

if __name__ == "__main__":
    _m.main()
else:
    sys.modules[__name__] = _m
