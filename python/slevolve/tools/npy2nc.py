#!/usr/bin/env python3
"""npy2nc.py -- export a directory of per-step .npy grid snapshots (tgrid<N>.npy,
optionally ocean/beta/R/G) written by the Fortran or Python solver into one NetCDF4 file.

  npy2nc.py RUNDIR --out run.nc [--fields tgrid,ocean] [--times FILE | --dt YEARS]
            [--earth NAME] [--solver NAME] [--attr key=value ...] [--zlib]

Coordinates
  lat  Gauss-Legendre latitudes of the solver's grid, degrees north, row order of
       the .npy arrays (north -> south): asin(mu_j) from the solver's own gaulw
       (python/slevolve/sht.py), identical to INPUT_FILES/grid/GLlat_512.txt to < 1e-6 deg.
  lon  the solver's longitude of column i, i*360/nlon (INPUT_FILES/grid/GLlon_512.txt;
       the FFT puts column 0 at 0 deg).  --lon-convention sampled gives (i+1)*360/nlon,
       the one-cell-shifted mapping of the archived ice regridding
       (tools/regrid.py); the default stays with the solver grid.
  time years since the start of the run: OUTPUT_SLM/times (Fortran) if present or
       --times FILE, else step*--dt.
Global attributes: solver, Earth model, planet constants (Python solver.Planet,
overridden by run_log.json rhoi/radius when present), conversion provenance.
"""
import argparse
import datetime
import json
import os
import re
import sys

import numpy as np
import warnings

warnings.filterwarnings("ignore", category=DeprecationWarning)  # netCDF4 vs numpy 2.5 shape setter

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # python/

UNITS = {"tgrid": ("m", "topography (bedrock + grounded ice), relative to sea level"),
         "ocean": ("1", "ocean function (1 = ocean)"),
         "beta": ("1", "grounded-ice mask function beta (1 = ice free)"),
         "R": ("m", "radial displacement of the solid surface"),
         "G": ("m", "geoid (sea-surface) perturbation")}


def gl_lat_deg(nlat):
    from ..sht import gaulw
    z, _ = gaulw(nlat)
    return np.degrees(np.arcsin(z))


def steps_of(rundir, field):
    pat = re.compile(rf"^{re.escape(field)}(\d+)\.npy$")
    return sorted(int(m.group(1)) for f in os.listdir(rundir) if (m := pat.match(f)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("rundir", help="directory holding tgrid<N>.npy (e.g. run/OUTPUT_SLM)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--fields", default="tgrid")
    ap.add_argument("--times", default=None, help="text file, one time (years) per step 0..N")
    ap.add_argument("--dt", type=float, default=1.0, help="years per step if no times file")
    ap.add_argument("--earth", default=None, help="Earth model file name (attribute)")
    ap.add_argument("--solver", default=None, help="solver/version string (attribute)")
    ap.add_argument("--lon-convention", default="solver", choices=("solver", "sampled"))
    ap.add_argument("--attr", action="append", default=[], help="extra global attribute key=value")
    ap.add_argument("--zlib", action="store_true", help="deflate-compress the variables")
    a = ap.parse_args()
    import netCDF4

    fields = [f for f in a.fields.split(",") if f]
    steps = steps_of(a.rundir, fields[0])
    if not steps:
        sys.exit(f"no {fields[0]}<N>.npy in {a.rundir}")
    for f in fields[1:]:
        s = steps_of(a.rundir, f)
        if s != steps:
            sys.exit(f"field {f} has steps {s[:5]}... differing from {fields[0]}")
    a0 = np.load(os.path.join(a.rundir, f"{fields[0]}{steps[0]}.npy"), mmap_mode="r")
    nlat, nlon = a0.shape

    tfile = a.times or (os.path.join(a.rundir, "times") if os.path.exists(os.path.join(a.rundir, "times")) else None)
    if tfile:
        tall = np.loadtxt(tfile, ndmin=1)
        time = tall[steps]
        tsrc = f"file {os.path.abspath(tfile)} (row = step)"
    else:
        time = np.asarray(steps, dtype=float) * a.dt
        tsrc = f"step * {a.dt} yr"

    from ..solver import Planet
    p = Planet()
    log = {}
    lj = os.path.join(a.rundir, "run_log.json")
    if os.path.exists(lj):
        log = json.load(open(lj))
        if log.get("rhoi"):
            p.rhoi = float(log["rhoi"])
        if log.get("radius"):
            p.radius = float(log["radius"])
            p.moiA = 0.3296145 * p.mass * p.radius**2
            p.moiC = 0.3307007 * p.mass * p.radius**2
    earth = a.earth or (os.path.basename(log["earth"]) if log.get("earth") else "unknown")
    if a.solver:
        solver = a.solver
    elif os.path.exists(os.path.join(a.rundir, "snap_index.txt")):
        solver = "SLEvolve (Fortran)"
    elif log:
        solver = f"SLEvolve (Python), opt={','.join(log.get('opt', []))}"
    else:
        solver = "unknown"

    nc = netCDF4.Dataset(a.out, "w", format="NETCDF4")
    nc.createDimension("time", len(steps))
    nc.createDimension("lat", nlat)
    nc.createDimension("lon", nlon)
    vt = nc.createVariable("time", "f8", ("time",))
    vt.units = "years since start of simulation"
    vt.long_name = "model time"
    vt.source = tsrc
    vs = nc.createVariable("step", "i4", ("time",))
    vs.long_name = "solver time-step index (file suffix N)"
    vla = nc.createVariable("lat", "f8", ("lat",))
    vla.units = "degrees_north"
    vla.long_name = "Gauss-Legendre latitude (row j of the grid, north to south)"
    vlo = nc.createVariable("lon", "f8", ("lon",))
    vlo.units = "degrees_east"
    vlo.long_name = "longitude (column i of the grid)"
    vlo.convention = ("i*360/nlon (solver SHT grid, GLlon file)" if a.lon_convention == "solver"
                      else "(i+1)*360/nlon (one-cell-shifted ice-sampling convention)")
    vt[:] = time
    vs[:] = steps
    vla[:] = gl_lat_deg(nlat)
    off = 0 if a.lon_convention == "solver" else 1
    vlo[:] = (np.arange(nlon) + off) * 360.0 / nlon
    for f in fields:
        v = nc.createVariable(f, "f8", ("time", "lat", "lon"), zlib=a.zlib, complevel=4,
                              chunksizes=(1, nlat, nlon))
        u, ln = UNITS.get(f, ("1", f))
        v.units, v.long_name = u, ln
        for k, s in enumerate(steps):
            v[k, :, :] = np.load(os.path.join(a.rundir, f"{f}{s}.npy")).astype(np.float64)
    g = {
        "title": "sea-level model snapshots",
        "Conventions": "CF-1.8",
        "solver": solver,
        "earth_model": earth,
        "grid": f"Gauss-Legendre {nlat} x {nlon}",
        "planet_radius_m": p.radius, "planet_mass_kg": p.mass,
        "rho_ice_kg_m3": p.rhoi, "rho_water_kg_m3": p.rhow,
        "g_m_s2": p.gacc, "omega_rad_s": p.omega,
        "moi_A_kg_m2": p.moiA, "moi_C_kg_m2": p.moiC, "fluid_love_kf": p.kf,
        "constants_source": "slevolve/solver.py Planet (planets_mod::earth_init)" + (" + run_log.json" if log else ""),
        "source_directory": os.path.abspath(a.rundir),
        "history": f"{datetime.datetime.now().isoformat(timespec='seconds')} npy2nc.py "
                   + " ".join(sys.argv[1:]),
    }
    for kv in a.attr:
        k, _, v = kv.partition("=")
        g[k] = v
    nc.setncatts(g)
    nc.close()
    print(f"npy2nc: {len(steps)} steps x {fields} -> {a.out}")


if __name__ == "__main__":
    main()
