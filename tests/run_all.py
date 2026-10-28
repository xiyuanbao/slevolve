#!/usr/bin/env python3
"""Run all release tests: python tests/run_all.py  (numpy required, ducc0 optional; no pytest)."""
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
TESTS = ["test_recursion.py", "test_sht_roundtrip.py", "test_small_case.py", "test_extrap_formula.py",
         "test_small_case_vardt.py", "test_floating_ice.py", "test_ramp.py"]


def main():
    t0 = time.perf_counter()
    failed = []
    for t in TESTS:
        print(f"== {t}", flush=True)
        r = subprocess.run([sys.executable, os.path.join(HERE, t)])
        if r.returncode != 0:
            failed.append(t)
    print(f"== {len(TESTS) - len(failed)}/{len(TESTS)} passed in {time.perf_counter() - t0:.1f} s"
          + (f"; FAILED: {', '.join(failed)}" if failed else ""))
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
