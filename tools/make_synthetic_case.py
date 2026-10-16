#!/usr/bin/env python3
"""make_synthetic_case.py -- write a complete INPUT_FILES tree for a synthetic Earth.

    python tools/make_synthetic_case.py --outdir SYNTH/INPUT_FILES                 # 512 x 1024 grid
    python tools/make_synthetic_case.py --outdir QS/INPUT_FILES --nglv 128          # 128 x 256 grid
    python tools/make_synthetic_case.py --outdir X/INPUT_FILES --topography bed.npy # real bedrock grid
    python tools/make_synthetic_case.py --write-seed reference/synthetic/topography_seed.npy  # maintainers

The synthetic Earth
-------------------
* topography (others/synth_topo.npy): an Earth-like random planet, built once on the 512 x 1024 grid by
  build_topography_seed() below and shipped as reference/synthetic/topography_seed.npy (rounded to
  whole millimetres).  Both generators read that file, so the Python and the Fortran twin
  tools/make_synthetic_case.f90 write identical topographies without the Fortran program needing a
  spherical-harmonic synthesis or numpy's random-number generator.  The construction:

  1. a Gaussian random field on the sphere by spectral synthesis with the shipped SHT
     (python/slevolve/sht.py), random coefficients from numpy.random.default_rng(SEED), degrees 2..200,
     power per degree proportional to l^-2 (variance of each coefficient proportional to l^-3;
     with l^-2 per coefficient the field is dominated by small scales and the "continents" are
     speckle), normalised to unit variance; two broad Gaussian bumps are added so that there is
     land at the south pole and around 60 N 130 E (the two ice-sheet sites);
  2. a hypsometric transfer: the area-weighted quantile q of the field is mapped to elevation:
     29 % land (q > 0.71: +5 m at the coast rising to a few hundred metres, with a tail to about
     2 km), a shelf of 3 % of the area (q 0.68-0.71: -200 m rising to the coast, mostly -100 to
     -200 m), a continental slope (q 0.62-0.68: -200 m to -3400 m) and the ocean floor (q < 0.62:
     -4500 m to -3700 m with q) with mid-ocean-ridge-like bumps up to about -2500 m where a
     low-pass version (degrees <= 30) of the field is near its median over the ocean floor.  The
     map is continuous, so the coast is smooth;
  3. low-amplitude fractal detail: the elevation is multiplied by 1 + a d, where d is a second
     unit-variance random field (degrees 20..250, coefficient variance l^-2, clipped to +-3) and
     a = 0.08 on land, 0.025 below sea level (the sign of the elevation, and so the coast, is kept).

  For other grids (--nglv N != 512) the seed is sampled at the nearest grid point (latitude row
  nearest in latitude; longitude column (2 i 512 + N) div (2 N)).
* ice (icemodel/synth_ice_<k>.npy, k = 0..8): two Vialov (perfectly plastic / Glen n = 3 equilibrium)
  profiles H = Hk (1 - (r / R)^(4/3))^(3/8), r = angular distance from the centre:
  a polar ice sheet centred on the south pole (R = 17 deg, H0 = 3000 m) and a mid-latitude ice sheet
  centred at 56 N 132 E (R = 12 deg, H0 = 2500 m), both only where the initial topography is above
  -500 m (ice on land or shelf, none over the deep ocean).  Both melt over 8 steps at fixed radius
  (the profile is scaled by Hk; POLAR and MIDLAT below; printed by the script).
* earthmodel/prem_512.l96C.ump3.lm1: copied from lovenum/ (Case 3 Earth model).
* grid/GLlat_<N>.txt, GLlon_<N>.txt: as tools/make_grid.py (the solvers do not read them).

--topography FILE replaces the synthetic topography by a real bedrock grid (an (N, 2N) .npy in the
SLEvolve layout, N = --nglv); the ice geometry and schedule stay the same.  Not used by default.

Arrays are (N, 2N) little-endian float64 .npy, rows north to south (Gauss-Legendre latitudes),
columns at longitudes i*360/(2N).  Topography and ice are rounded to whole millimetres: the values
then survive the ES16.9 text round trip of the initial topography in the Fortran forward driver
(legacy behaviour), and the Python and Fortran generators, whose elementary functions may differ in
the final digit, write identical files.  Run the case with

    NSTEPS=8 ICEMODEL=synth_ice_ TOPOMODEL=synth_topo PLANETMODEL=prem_512.l96C.ump3.lm1 \\
        bash fortran/build.sh                               # (tests/test_fortran_python.sh does this)
    python python/run_ais.py --steps 8 --lmax 511 --icedir SYNTH/INPUT_FILES/icemodel \\
        --icepfx synth_ice_ --topo SYNTH/INPUT_FILES/others/synth_topo.npy \\
        --earth SYNTH/INPUT_FILES/earthmodel/prem_512.l96C.ump3.lm1 --out OUT_synth

Needs numpy only (--write-seed also needs python/slevolve/sht.py, which needs numpy).
"""
import argparse
import os
import shutil
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..")
LOVE = "prem_512.l96C.ump3.lm1"
SEED_FILE = os.path.join(ROOT, "reference", "synthetic", "topography_seed.npy")
SEED_NGLV = 512
SEED = 20261013
POLAR = [3000.0, 2625.0, 2250.0, 1875.0, 1500.0, 1125.0, 750.0, 375.0, 0.0]   # polar sheet, steps 0..8
MIDLAT = [2500.0, 2400.0, 2150.0, 1800.0, 1350.0, 900.0, 500.0, 200.0, 0.0]   # mid-latitude sheet
HMAX = POLAR                     # backwards-compatible name (polar ice thickness at the pole)
# (name, centre latitude, centre longitude, radius in degrees) of the two ice sheets
SHEETS = (("polar", -90.0, 0.0, 17.0), ("mid-latitude", 56.0, 132.0, 12.0))
CAP_RADIUS = SHEETS[0][3]
ICE_MIN_BED = -500.0           # ice only where the initial topography is above -500 m
D2R = np.pi / 180.0
# hypsometric transfer (area quantiles of the random field)
Q_LAND, Q_SHELF, Q_FLOOR = 0.71, 0.68, 0.62


def gl_lat_lon(nglv):
    """(lat_deg north to south, lon_deg) of the nglv x 2*nglv Gauss-Legendre grid (tools/make_grid.py)."""
    x, _ = np.polynomial.legendre.leggauss(nglv)
    return np.arcsin(x[::-1]) * (180.0 / np.pi), np.arange(2 * nglv) * 360.0 / (2 * nglv)


def smoothstep(x):
    x = np.clip(x, 0.0, 1.0)
    return x * x * (3.0 - 2.0 * x)


def dist(LAT, LON, lat0, lon0):
    """Angular distance (deg) of every grid point from (lat0, lon0)."""
    p, p0, dl = LAT * D2R, lat0 * D2R, (LON - lon0) * D2R
    c = np.sin(p) * np.sin(p0) + np.cos(p) * np.cos(p0) * np.cos(dl)
    return np.arccos(np.clip(c, -1.0, 1.0)) / D2R


def round_mm(f):
    """Round to whole millimetres (np.round(1000 f) / 1000; the Fortran twin uses anint)."""
    return np.round(f * 1000.0) / 1000.0


def build_topography_seed(seed=SEED, nglv=SEED_NGLV):
    """Earth-like random topography on the nglv x 2 nglv grid (see the module docstring), rounded to mm."""
    sys.path.insert(0, os.path.join(ROOT, "python"))
    from slevolve.sht import SHT
    L = 255
    lat, lon = gl_lat_lon(nglv)
    LON, LAT = np.meshgrid(lon, lat)
    W = np.polynomial.legendre.leggauss(nglv)[1][::-1][:, None] * np.ones((1, 2 * nglv))
    rng = np.random.default_rng(seed)
    sht = SHT(nglv, 2 * nglv, L)
    l = np.arange(L + 1)[:, None]
    m = np.arange(L + 1)[None, :]

    def unit(g):
        g = g - (W * g).sum() / W.sum()
        return g / np.sqrt((W * g * g).sum() / W.sum())

    def field(l0, l1, p, lcut=None):
        s = np.where((l >= l0) & (l <= l1) & (m <= l), np.maximum(l, 1.0) ** (-p / 2.0), 0.0)
        c = (rng.standard_normal(s.shape) + 1j * rng.standard_normal(s.shape)) * s / np.sqrt(2.0)
        c[:, 0] = c[:, 0].real * np.sqrt(2.0)
        g = unit(sht.synthesis(c))
        if lcut is None:
            return g
        c[lcut + 1:] = 0.0
        return g, unit(sht.synthesis(c))

    g, glow = field(2, 200, 3.0, lcut=30)
    g = g + 3.5 * np.exp(-(dist(LAT, LON, -90.0, 0.0) / 22.0) ** 2) \
          + 2.5 * np.exp(-(dist(LAT, LON, 60.0, 130.0) / 15.0) ** 2)
    detail = field(20, 250, 2.0)
    # area-weighted quantile of every grid point
    o = np.argsort(g, axis=None, kind="stable")
    q = np.empty(g.size)
    q[o] = np.cumsum(W.ravel()[o]) / W.sum()
    q = q.reshape(g.shape) - 0.5 * W / W.sum()
    s = np.clip((q - Q_LAND) / (1.0 - Q_LAND), 0.0, 1.0)
    land = 5.0 + 700.0 * s ** 1.6 + 1400.0 * s ** 12
    shelf = -200.0 + 195.0 * np.clip((q - Q_SHELF) / (Q_LAND - Q_SHELF), 0.0, 1.0) ** 3
    slope = -3400.0 + 3200.0 * smoothstep((q - Q_FLOOR) / (Q_SHELF - Q_FLOOR))
    oc = q < Q_FLOOR
    med, sd = np.median(glow[oc]), glow[oc].std()
    floor = -4500.0 + 800.0 * (q / Q_FLOOR) + 1900.0 * np.exp(-((glow - med) / (0.3 * sd)) ** 2)
    e = np.where(q >= Q_LAND, land, np.where(q >= Q_SHELF, shelf, np.where(q >= Q_FLOOR, slope, floor)))
    e = e * (1.0 + np.where(e > 0.0, 0.08, 0.025) * np.clip(detail, -3.0, 3.0))
    return round_mm(e)


def sample_seed(seed, nglv):
    """Nearest-grid-point sample of the 512 x 1024 seed on the nglv x 2 nglv grid (as the Fortran twin)."""
    ns = seed.shape[0]
    if nglv == ns:
        return seed.copy()
    lat_s, _ = gl_lat_lon(ns)
    lat, _ = gl_lat_lon(nglv)
    rows = np.array([int(np.argmin(np.abs(lat_s - la))) for la in lat])
    cols = (2 * np.arange(2 * nglv) * ns + nglv) // (2 * nglv) % (2 * ns)
    return np.ascontiguousarray(seed[rows][:, cols])


def vialov(d, R):
    """Vialov profile shape (1 - (d / R)^(4/3))^(3/8) inside d < R, 0 outside."""
    return np.clip(1.0 - (np.minimum(d, R) / R) ** (4.0 / 3.0), 0.0, None) ** 0.375


def synthetic_earth(nglv, topography=None, seed_file=SEED_FILE):
    """Return (lat, lon, topo, [ice_0..ice_8]) on the nglv x 2*nglv grid."""
    lat, lon = gl_lat_lon(nglv)
    LON, LAT = np.meshgrid(lon, lat)
    if topography is None:
        topo = sample_seed(np.load(seed_file).astype(np.float64), nglv)
    else:
        topo = np.load(topography).astype(np.float64)
        if topo.shape != LAT.shape:
            raise SystemExit(f"--topography: shape {topo.shape}, expected {LAT.shape}")
    (_, la1, lo1, R1), (_, la2, lo2, R2) = SHEETS
    polar = vialov(90.0 + LAT, R1)          # angular distance from the south pole
    mid = vialov(dist(LAT, LON, la2, lo2), R2)
    grounded = topo > ICE_MIN_BED       # ice on land or shelf only
    ice = [round_mm((hp * polar + hm * mid) * grounded) for hp, hm in zip(POLAR, MIDLAT)]
    return lat, lon, topo, ice


def describe():
    s = [f"seed {SEED}: random topography (reference/synthetic/topography_seed.npy), two Vialov ice sheets"]
    for (name, la, lo, R), H in zip(SHEETS, (POLAR, MIDLAT)):
        s.append(f"  {name} ice sheet: centre {la:.1f} N {lo:.1f} E, radius {R:.1f} deg, H0 {H[0]:.0f} m")
    s.append("  step:            " + " ".join(f"{k:6d}" for k in range(len(POLAR))))
    s.append("  polar ice (m):   " + " ".join(f"{h:6.0f}" for h in POLAR))
    s.append("  mid-lat ice (m): " + " ".join(f"{h:6.0f}" for h in MIDLAT))
    return "\n".join(s)


def write(outdir, nglv, topography=None):
    lat, lon, topo, ice = synthetic_earth(nglv, topography)
    for d in ("grid", "icemodel", "earthmodel", "others"):
        os.makedirs(os.path.join(outdir, d), exist_ok=True)
    for name, v in ((f"GLlat_{nglv}.txt", lat), (f"GLlon_{nglv}.txt", lon)):
        with open(os.path.join(outdir, "grid", name), "w") as fh:
            fh.write("".join("%16.7e\n" % x for x in v))
    np.save(os.path.join(outdir, "others", "synth_topo.npy"), np.ascontiguousarray(topo, dtype="<f8"))
    for k, h in enumerate(ice):
        np.save(os.path.join(outdir, "icemodel", f"synth_ice_{k}.npy"), np.ascontiguousarray(h, dtype="<f8"))
    shutil.copy(os.path.join(ROOT, "lovenum", LOVE), os.path.join(outdir, "earthmodel", LOVE))
    return lat, lon, topo, ice


def main():
    ap = argparse.ArgumentParser(description="write a synthetic-Earth INPUT_FILES tree (random Earth-like "
                                             "topography, two melting ice sheets, 8 steps)")
    ap.add_argument("--outdir", help="INPUT_FILES directory to create")
    ap.add_argument("--nglv", type=int, default=512,
                    help="number of Gauss-Legendre latitudes N; grid N x 2N (default 512, the Fortran grid)")
    ap.add_argument("--topography", metavar="FILE", default=None,
                    help="use this (N, 2N) .npy bedrock/topography grid instead of the synthetic one")
    ap.add_argument("--write-seed", metavar="FILE", default=None,
                    help="(maintainers) rebuild the 512 x 1024 topography seed and write it to FILE")
    a = ap.parse_args()
    if a.write_seed:
        np.save(a.write_seed, build_topography_seed().astype("<f8"))
        print("wrote", a.write_seed)
        if not a.outdir:
            return
    if not a.outdir:
        ap.error("--outdir is required")
    lat, lon, topo, ice = write(a.outdir, a.nglv, a.topography)
    w = np.polynomial.legendre.leggauss(a.nglv)[1][::-1][:, None] * np.ones((1, 2 * a.nglv))
    land = (w * (topo >= 0)).sum() / w.sum()
    print(describe())
    print(f"wrote {a.outdir}: grid {a.nglv}x{2 * a.nglv}, land fraction {land:.3f}, "
          f"{len(ice)} ice files synth_ice_0..{len(ice) - 1}, topography others/synth_topo.npy"
          f"{' (from ' + a.topography + ')' if a.topography else ''}, earthmodel/{LOVE}")


if __name__ == "__main__":
    main()
