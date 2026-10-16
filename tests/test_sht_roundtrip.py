#!/usr/bin/env python3
"""Spherical-harmonic analysis -> synthesis round trip on the 512 x 1024 Gauss grid.

A band-limited field is synthesised from random coefficients (degrees <= lmax), analysed and
synthesised again; the relative maximum difference of the coefficients and of the grid is checked.
numpy backend (python/slevolve/sht.py) and, if importable, ducc0 (python/slevolve/sht_ducc.py)."""
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "python"))
from slevolve.sht import SHT  # noqa: E402

NLAT, NLON = 512, 1024


def roundtrip(sht, lmax, rng):
    nl = lmax + 1
    c = np.tril(rng.normal(size=(nl, nl)) + 1j * rng.normal(size=(nl, nl)))
    c[:, 0] = c[:, 0].real  # m = 0 coefficients of a real field are real
    g = sht.synthesis(c, out=np.zeros((NLAT, NLON)))
    c2 = sht.analysis(g, out=np.zeros((nl, nl), complex))
    g2 = sht.synthesis(c2, out=np.zeros((NLAT, NLON)))
    ec = np.abs(c2 - c).max() / np.abs(c).max()
    eg = np.abs(g2 - g).max() / np.abs(g).max()
    return ec, eg


def main():
    rng = np.random.default_rng(7)
    worst = 0.0
    backends = [("numpy", lambda L: SHT(NLAT, NLON, L, workers=4))]
    try:
        import ducc0  # noqa: F401
        from slevolve.sht_ducc import make_sht

        backends.append(("ducc0", lambda L: make_sht("ducc", NLAT, NLON, L, workers=4)))
    except ImportError:
        print("test_sht_roundtrip: ducc0 not importable, numpy backend only")
    for name, mk in backends:
        for lmax in (63, 511):
            t0 = time.perf_counter()
            ec, eg = roundtrip(mk(lmax), lmax, rng)
            worst = max(worst, ec, eg)
            print(f"test_sht_roundtrip: {name:5s} lmax {lmax:3d}: coefficients {ec:.3e}, "
                  f"grid {eg:.3e}  ({time.perf_counter() - t0:.2f} s)")
    assert worst < 1e-12, worst
    print("test_sht_roundtrip: PASS")


if __name__ == "__main__":
    main()
