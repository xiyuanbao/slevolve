#!/usr/bin/env python3
"""Script entry point for slevolve.tools.make_reference (kept so `python python/tools/make_reference.py ...` works)."""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from slevolve.tools import make_reference as _m  # noqa: E402

if __name__ == "__main__":
    _m.main()
else:
    sys.modules[__name__] = _m
