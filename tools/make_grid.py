#!/usr/bin/env python3
"""make_grid.py -- write the grid files GLlat_<N>.txt and GLlon_<N>.txt read by the Fortran solver.

    python tools/make_grid.py --outdir INPUT_FILES/grid            # N = 512 (512 x 1024 grid)
    python tools/make_grid.py --outdir INPUT_FILES/grid --nglv 512

GLlat_<N>.txt : the N Gauss-Legendre latitudes in DEGREES, north to south (first value +89.73 for
                N = 512), one value per line, printf format '%16.7e'.
GLlon_<N>.txt : the 2N equispaced longitudes i*360/(2N) in DEGREES east, i = 0..2N-1 (0 to 359.65
                for N = 512), same format.
These are the latitude rows and longitude columns of every grid file (ice, topography, outputs):
row j of a grid is latitude GLlat[j], column i is longitude GLlon[i].

For N = 512, compared with GLlat_512.txt and GLlon_512.txt of the Han et al. (2022) OSF archive:
the longitude file is identical; 58 of the 512 latitude lines differ in the last printed digit
(maximum difference 5.3e-7 degrees; the archived values agree with these latitudes rounded to 6
decimals except at +-84.819211). The solvers do not read these files. Needs numpy only.
"""
import argparse
import os

import numpy as np


def gl_lat_lon(nglv):
    """(lat_deg north to south, lon_deg) of the nglv x 2*nglv Gauss-Legendre grid."""
    x, _ = np.polynomial.legendre.leggauss(nglv)  # ascending sin(lat): south to north
    lat = np.degrees(np.arcsin(x[::-1]))
    lon = np.arange(2 * nglv) * 360.0 / (2 * nglv)
    return lat, lon


def main():
    ap = argparse.ArgumentParser(description="write GLlat_<N>.txt and GLlon_<N>.txt")
    ap.add_argument("--outdir", required=True, help="directory to write into (e.g. INPUT_FILES/grid)")
    ap.add_argument("--nglv", type=int, default=512, help="number of Gauss-Legendre latitudes N (default 512)")
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    lat, lon = gl_lat_lon(a.nglv)
    for name, v in ((f"GLlat_{a.nglv}.txt", lat), (f"GLlon_{a.nglv}.txt", lon)):
        p = os.path.join(a.outdir, name)
        with open(p, "w") as fh:
            fh.write("".join("%16.7e\n" % x for x in v))
        print("wrote", p, len(v), "values")


if __name__ == "__main__":
    main()
