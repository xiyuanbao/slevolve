#!/usr/bin/env python3
"""text2npy.py -- legacy list-directed grid text -> .npy (shape (nlat, nlon), C order, stream order).

  text2npy.py [--f4] [--nlat 512] [--nlon 1024] [--outdir DIR] FILE [FILE ...]
  text2npy.py [--f4] [--outdir DIR] --dir DIR --pattern 'AISload_gd_Run85_0_0_*'

--f4 stores '<f4' when the file is exactly representable in float32, else '<f8'
(with a warning).  Output: <outdir or input dir>/<basename>.npy.
Same conventions as fortran/tools/text2npy.f90 (see fortran/IO_FORMATS.md).
"""
import argparse, fnmatch, os, sys
from .. import npyio  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("files", nargs="*")
    ap.add_argument("--f4", action="store_true")
    ap.add_argument("--nlat", type=int, default=512)
    ap.add_argument("--nlon", type=int, default=1024)
    ap.add_argument("--outdir", default=None)
    ap.add_argument("--dir", default=None)
    ap.add_argument("--pattern", default="*")
    a = ap.parse_args()
    files = list(a.files)
    if a.dir:
        files += sorted(os.path.join(a.dir, f) for f in os.listdir(a.dir)
                        if fnmatch.fnmatch(f, a.pattern) and not f.endswith((".npy", ".bin"))
                        and os.path.isfile(os.path.join(a.dir, f)))
    if not files:
        ap.error("no input files")
    nf4 = 0
    for f in files:
        g = npyio.read_text_grid(f, a.nlat, a.nlon)
        out = os.path.join(a.outdir or os.path.dirname(f), os.path.basename(f) + ".npy")
        dt = npyio.save_grid(out, g, "f4" if a.f4 else "f8")
        if a.f4 and dt != "f4":
            print(f"text2npy: {os.path.basename(f)} is not exactly representable in float32 -> stored as <f8",
                  file=sys.stderr)
        nf4 += dt == "f4"
    print(f"text2npy: converted {len(files)} files ({nf4} as <f4)")


if __name__ == "__main__":
    main()
