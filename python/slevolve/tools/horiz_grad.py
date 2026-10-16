#!/usr/bin/env python3
"""
horiz_grad.py -- grid components of the horizontal displacement from the
spectral scalar Psi written by either solver (horiz<N>.npy, '<c16' (L+1, L+1),
spharmt U convention; Fortran with SLEVOLVE_HORIZ>=1, Python with --horizontal).

    u_theta = dPsi/dtheta   (theta = colatitude, positive SOUTHWARD)
    u_phi   = (1/sin theta) dPsi/dphi   (positive EASTWARD)
    u_north = -u_theta,  u_east = u_phi   (metres)

Usage:
    horiz_grad.py [--nlat 512] [--nlon 1024] [--threads 4] [--outdir D] horiz550.npy ...
        -> utheta550.npy, uphi550.npy  ('<f8' (nlat, nlon), rows north -> south)
    horiz_grad.py --selftest
Needs ducc0 (python/slevolve/horiz.py).
"""
import argparse
import os
import re
import sys

import numpy as np

from ..horiz import gradient  # noqa: E402


def selftest():
    from ..sht_ducc import DuccSHT

    nlat, nlon, L = 64, 128, 63
    sh = DuccSHT(nlat, nlon, L)
    x, _ = np.polynomial.legendre.leggauss(nlat)
    th = np.arccos(x[::-1])[:, None]
    ph = (2 * np.pi * np.arange(nlon) / nlon)[None, :]
    s, c = np.sin(th), np.cos(th)
    psi = c**3 + s**2 * c * np.sin(2 * ph) + s * np.cos(ph)
    dth = -3 * c**2 * s + (2 * s * c**2 - s**3) * np.sin(2 * ph) + c * np.cos(ph)
    dph = 2 * s * c * np.cos(2 * ph) - np.sin(ph)
    ut, up = gradient(sh.analysis(psi).copy(), nlat, nlon)
    et, ep = np.abs(ut - dth).max(), np.abs(up - dph).max()
    print(f"selftest: max|u_theta err| {et:.2e}, max|u_phi err| {ep:.2e}")
    return et < 1e-10 and ep < 1e-10


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("files", nargs="*")
    ap.add_argument("--nlat", type=int, default=512)
    ap.add_argument("--nlon", type=int, default=1024)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--outdir", default=None)
    ap.add_argument("--selftest", action="store_true",
                    help="verification only: check the gradient against analytic fields and exit")
    a = ap.parse_args()
    if a.selftest:
        sys.exit(0 if selftest() else 1)
    for f in a.files:
        U = np.load(f)
        ut, up = gradient(U, a.nlat, a.nlon, nthreads=a.threads)
        d = a.outdir or os.path.dirname(os.path.abspath(f))
        m = re.match(r"horiz(.*)\.npy$", os.path.basename(f))
        tag = m.group(1) if m else os.path.splitext(os.path.basename(f))[0]
        np.save(os.path.join(d, f"utheta{tag}.npy"), ut)
        np.save(os.path.join(d, f"uphi{tag}.npy"), up)
        print(f"{f} -> utheta{tag}.npy, uphi{tag}.npy  (max |u_h| {np.hypot(ut, up).max():.4e} m)")


if __name__ == "__main__":
    main()
