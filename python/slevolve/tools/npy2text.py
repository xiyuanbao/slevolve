#!/usr/bin/env python3
"""npy2text.py -- .npy grids -> legacy grid text formats.

  npy2text.py [--format es16.9|row-e16.7|row-e.18|row-E18.10] [--ncol 1024] [--outdir DIR] FILE.npy ...

es16.9      %16.9E one per line (ES16.9E2 one-value-per-line output format)  [default]
row-e16.7   %16.7e, 1024 per line, no separator (OSF archive topo/tgrid)
row-e.18    %.18e, 1024 per line, ' ' separator (np.savetxt default, OSF ice)
row-E18.10  %18.10E, 1024 per line, ' ' separator (nadi-case ice files)
Values are written in C-flat (stream) order.  Output: <basename without .npy>.
"""
import argparse, os, sys
import numpy as np  # noqa: E402
from .. import npyio  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("files", nargs="+")
    ap.add_argument("--format", default="es16.9", choices=sorted(npyio.TEXT_FORMATS))
    ap.add_argument("--ncol", type=int, default=None)
    ap.add_argument("--outdir", default=None)
    a = ap.parse_args()
    for f in a.files:
        arr = np.load(f)
        base = os.path.basename(f)
        base = base[:-4] if base.endswith(".npy") else base
        out = os.path.join(a.outdir or os.path.dirname(f), base)
        npyio.write_text_grid(out, np.asarray(arr, dtype=np.float64), a.format, a.ncol)
    print(f"npy2text: wrote {len(a.files)} files ({a.format})")


if __name__ == "__main__":
    main()
