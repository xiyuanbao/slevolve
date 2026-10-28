#!/usr/bin/env python3
"""
run_ais.py -- drive the NumPy exact-memory sea-level solver on the Antarctic
forward benchmark of Han et al. (2022) (ice files AISload_gd_Run85_0_0_k, k = 0..550), or on any
ice history <icepfx>0..N with a uniform step --dt (years; default 1) or the non-uniform steps
of a times file (--times: one time in years per ice file).

Threads are pinned *before* NumPy is imported.

Example
    python3 run_ais.py --steps 50 --threads 4 --out /path/to/output \
        --earth /path/to/INPUT_FILES/earthmodel/<file> --icedir /path/to/INPUT_FILES/icemodel \
        --topo /path/to/INPUT_FILES/others/topo0_archive
"""

import argparse
import os
import sys


def _parse_argv():
    ap = argparse.ArgumentParser()
    ap.add_argument("--steps", type=int, default=50, help="number of time steps N; reads ice files 0..N (default 50)")
    ap.add_argument("--dt", type=float, default=1.0,
                    help="step length in YEARS (default 1); must equal the time spacing of the "
                         "ice files <icepfx>0, <icepfx>1, ... (file k is the ice at time k*dt)")
    ap.add_argument("--times", default=None,
                    help="optional times file: one time t_k in years per line for ice files 0..N, strictly "
                         "increasing; step k has length t_k - t_(k-1) and --dt is ignored")
    ap.add_argument("--threads", type=int, default=1, help="threads for the transforms and numba kernels (set before NumPy is imported; default 1)")
    ap.add_argument("--out", default=None, help="output directory (pass explicitly; development default <root>/results/slevolve_nNN_tT)")
    ap.add_argument("--root", default=None, help="root directory with INPUT_FILES/ (default: parent of python/)")
    ap.add_argument("--norder", type=int, default=512, help="spectral truncation (default 512; the benchmark cases use 511)")
    ap.add_argument("--lmax", type=int, default=None,
                    help="spectral truncation; alias for --norder (ducc backend needs <= nglv-1)")
    ap.add_argument("--sht", default="numpy", choices=("numpy", "ducc"),
                    help="SHT backend (default numpy)")
    ap.add_argument("--nglv", type=int, default=512, help="number of Gauss-Legendre latitudes of the grid (nglv x 2*nglv; default 512)")
    ap.add_argument("--tpw", default="T", help="T|F rotational feedback (default T)")
    ap.add_argument("--checkmarine", default="T", help="T|F marine check: ice below flotation is removed from the load (default T)")
    ap.add_argument("--calcRG", default="T", help="T|F compute the radial displacement R and geoid G grids (default T)")
    ap.add_argument("--epsilon1", type=float, default=1.0e-5, help="inner-iteration tolerance (default 1e-5)")
    ap.add_argument("--profile", action="store_true", help="print a per-stage timing profile (log only)")
    ap.add_argument("--cprofile", action="store_true", help="run under cProfile, write cprofile.txt to --out")
    ap.add_argument("--save-fields", default="tgrid,ocean,beta",
                    help="comma list of per-step outputs: tgrid,ocean,beta,dS,R,G,dSLdt,dRdt (default "
                         "tgrid,ocean,beta); dSLdt/dRdt are the analytic rates (m/yr) at the end of each "
                         "step, with TPWrate<N>.txt when --tpw T")
    ap.add_argument("--ramp", action="store_true",
                    help="second-order load placement: each load increment accrues linearly over its "
                         "step (default off: step placement of Kendall et al. 2005)")
    ap.add_argument("--horizontal", action="store_true",
                    help="also write the horizontal-displacement scalar Psi_lm (horiz<N>.npy, "
                         "u_h = grad_h Psi); output only, default off")
    ap.add_argument("--horizontal-grid", action="store_true",
                    help="with --horizontal: also write utheta<N>.npy / uphi<N>.npy (ducc0)")
    ap.add_argument("--cache-only", action="store_true", help="only build the ice .npy cache")
    ap.add_argument("--rebuild-cache", action="store_true",
                    help="rebuild the ice cache even if its manifest (directory, prefix, file count, "
                         "size and mtime of every ice file) matches")
    ap.add_argument("--chunk", type=int, default=16, help="block size of the transform batches (performance only; default 16)")
    ap.add_argument("--mem-threads", type=int, default=1,
                    help="threads for the hereditary advance (identical results; default 1)")
    ap.add_argument("--opt", default="",
                    help="optional optimisations: comma list of raw,diag,dsl,numba,sums "
                         "or 'all'; default '' = none")
    # ---- iteration accelerations (all off by default) ----------------------
    ap.add_argument("--inner", default="off", choices=("off", "extrap", "anderson"),
                    help="inner-loop acceleration (default off)")
    ap.add_argument("--anderson-depth", type=int, default=3, help="with --inner anderson: mixing depth m (default 3)")
    ap.add_argument("--anderson-freeze", type=int, default=2,
                    help="freeze the ocean mask after this plain iteration (1 or 2)")
    ap.add_argument("--count-sht", action="store_true",
                    help="count analyses/syntheses during the time loop (log only)")
    ap.add_argument("--anderson-keep", action="store_true",
                    help="keep the Anderson history across a mask refresh")
    ap.add_argument("--rhoi", type=float, default=None, help="override Planet.rhoi (default 920)")
    ap.add_argument("--radius", type=float, default=None,
                    help="override Planet.radius (default 6.371e6); moiA/moiC follow")
    ap.add_argument("--icedir", default=None,
                    help="ice directory: read <icedir>/<icepfx><k>{.npy,.bin,text} directly (npy preferred), no cache; pass explicitly")
    ap.add_argument("--icepfx", default="AISload_gd_Run85_0_0_", help="ice file prefix: files <icedir>/<icepfx><k> (default: the Case 1 name)")
    ap.add_argument("--topo", default=None, help="initial topography (.npy, .bin LE f8, or text); pass explicitly")
    ap.add_argument("--io-compat", default="fortran", choices=("fortran", "legacy"),
                    help="fortran: ocean/beta saved as <f8 and dS as fortran_order <c16, exactly the "
                         "Fortran solver's .npy files (fortran/IO_FORMATS.md); legacy: ocean/beta as uint8")
    ap.add_argument("--earth", default=None,
                    help="Love-number file (default <root>/INPUT_FILES/earthmodel/prem_512.l60K2C.sum18p6.dum19p2.tz19p4.lm22, "
                         "pass --earth, --icedir and --topo explicitly)")
    return ap.parse_args()


ARGS = _parse_argv()
_NT = str(ARGS.threads)
for _v in (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "MKL_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "NUMBA_NUM_THREADS",
):
    os.environ[_v] = _NT

import json  # noqa: E402
import subprocess  # noqa: E402
import time  # noqa: E402

import numpy as np  # noqa: E402

from .sht_ducc import make_sht  # noqa: E402
from .solver import LoveNumbers, Planet, SeaLevelSolver, read_grid_text  # noqa: E402
from .solver_opt import OptSeaLevelSolver, parse_opts  # noqa: E402
from . import npyio  # noqa: E402
from . import _paths  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = ARGS.root or _paths.ROOT
ICEDIR = os.path.join(ROOT, "INPUT_FILES", "icemodel")
ICEPFX = "AISload_gd_Run85_0_0_"
TOPO = os.path.join(ROOT, "INPUT_FILES", "others", "topo0_archive")
EARTH = os.path.join(
    ROOT, "INPUT_FILES", "earthmodel", "prem_512.l60K2C.sum18p6.dum19p2.tz19p4.lm22"
)
CACHE = _paths.CACHE
if ARGS.topo:
    TOPO = ARGS.topo
if ARGS.earth:
    EARTH = ARGS.earth


def make_planet(rhoi=None, radius=None):
    """Planet() with optional rhoi / radius overrides (moiA, moiC recomputed
    exactly as planets_mod does: 0.3296145*mass*radius**2 etc.)."""
    p = Planet()
    if rhoi is not None:
        p.rhoi = float(rhoi)
    if radius is not None:
        p.radius = float(radius)
        p.moiA = 0.3296145 * p.mass * p.radius**2
        p.moiC = 0.3307007 * p.mass * p.radius**2
    return p


def read_grid_any(path, nlat, nlon):
    # .npy / .bin / text (and <path>.npy / <path>.bin next to a text path)
    return npyio.read_grid(path, nlat, nlon)


class BinIceDir:
    """ice[k] -> <dir>/<pfx><k>.bin read on demand (little-endian float64)."""

    def __init__(self, d, pfx, nlat, nlon):
        self.d, self.pfx, self.shape = d, pfx, (nlat, nlon)

    def __getitem__(self, k):
        return np.fromfile(os.path.join(self.d, f"{self.pfx}{k}.bin"), dtype="<f8").reshape(self.shape)


def uptime():
    try:
        return subprocess.check_output(["uptime"], text=True).strip()
    except Exception:
        return "?"


def loadavg():
    try:
        return os.getloadavg()
    except Exception:
        return (float("nan"),) * 3


def build_ice_cache(nmax, nglv):
    """Parse the text ice files once into a single memmapped float64 .npy stack."""
    os.makedirs(CACHE, exist_ok=True)
    path = os.path.join(CACHE, f"ice_{ICEPFX}0_{nmax}.npy")
    man = npyio.ice_manifest(ICEDIR, ICEPFX, nmax + 1)
    # reuse an existing cache (this length or longer) whose manifest matches the ice files
    if not ARGS.rebuild_cache:
        cands = []
        for f in sorted(os.listdir(CACHE)):
            if f.startswith(f"ice_{ICEPFX}0_") and f.endswith(".npy"):
                try:
                    n = int(f[len(f"ice_{ICEPFX}0_") : -4])
                except ValueError:
                    continue
                if n >= nmax:
                    cands.append(os.path.join(CACHE, f))
        reason = "no cache"
        for c in cands:
            ok, why = npyio.cache_matches(c, man)
            if ok:
                return c, 0.0
            reason = why
        if cands:
            print(f"[slevolve] ice cache out of date ({reason}): rebuilding {path}")
    else:
        print(f"[slevolve] --rebuild-cache: rebuilding {path}")
    t0 = time.perf_counter()
    arr = np.lib.format.open_memmap(
        path + ".tmp", mode="w+", dtype=np.float64, shape=(nmax + 1, nglv, 2 * nglv)
    )
    # Fast path (no text parsing): reuse the flat little-endian float64 binary cache
    # of the Fortran solver (INPUT_FILES/icemodel/<name>.bin or .npy, same
    # element order as the text files).  Falls back to the text parser per file.
    bindir = ICEDIR
    nbin = 0
    for k in range(nmax + 1):
        b = npyio.grid_source(bindir, f"{ICEPFX}{k}")
        if b is not None and b.endswith((".npy", ".bin")):
            arr[k] = npyio.read_grid(b, nglv, 2 * nglv)  # .npy (f4/f8) or .bin (LE f8)
            nbin += 1
        else:
            arr[k] = read_grid_text(os.path.join(ICEDIR, f"{ICEPFX}{k}"), nglv, 2 * nglv)
    print(f"[slevolve] ice cache built from {nbin} binary + {nmax + 1 - nbin} text files")
    arr.flush()
    del arr
    os.replace(path + ".tmp", path)
    npyio.write_manifest(path, man)
    return path, time.perf_counter() - t0


def main():
    a = ARGS
    nglv = a.nglv
    norder = a.lmax if a.lmax is not None else a.norder
    tpw = a.tpw.upper().startswith("T")
    cm = a.checkmarine.upper().startswith("T")
    rg = a.calcRG.upper().startswith("T")
    out = a.out or os.path.join(ROOT, "results", f"slevolve_n{a.steps}_t{a.threads}")
    os.makedirs(out, exist_ok=True)
    fields = set(f for f in a.save_fields.split(",") if f)

    times = npyio.read_times(a.times, a.steps) if a.times else None
    dts = np.diff(times) if times is not None else None
    log = {
        "argv": sys.argv,
        "threads": a.threads,
        "steps": a.steps,
        "dt_years": a.dt if times is None else "variable",
        "times_file": os.path.abspath(a.times) if a.times else None,
        "times_years": [float(x) for x in times] if times is not None else None,
        "norder": norder,
        "sht_backend": a.sht,
        "mem_threads": a.mem_threads,
        "opt": sorted(parse_opts(a.opt)),
        "nglv": nglv,
        "tpw": tpw,
        "checkmarine": cm,
        "calcRG": rg,
        "epsilon1": a.epsilon1,
        "ramp": bool(a.ramp),
        "numpy": np.__version__,
        "uptime_start": uptime(),
        "loadavg_start": loadavg(),
    }
    print(
        f"[slevolve] backend={a.sht} lmax={norder} threads={a.threads} steps={a.steps} "
        + (f"dt={a.dt:g} yr " if times is None else
           f"times file {a.times}: t0={times[0]:g} yr, dt {dts.min():g}-{dts.max():g} yr ")
        + f"opt={sorted(parse_opts(a.opt))}  {log['uptime_start']}"
    )

    if a.icedir:
        ice = npyio.IceDir(a.icedir, a.icepfx, nglv, 2 * nglv)
        log["icedir"] = a.icedir
    else:
        icepath, tcache = build_ice_cache(a.steps, nglv)
        log["ice_cache_seconds"] = tcache
        log["ice_cache_path"] = icepath
        print(f"[slevolve] ice cache: {tcache:.1f} s -> {icepath}")
        if a.cache_only:
            return
        ice = np.load(icepath, mmap_mode="r")
    log.update({"inner": a.inner, "anderson_depth": a.anderson_depth, "rhoi": a.rhoi,
                "radius": a.radius, "topo": TOPO, "earth": EARTH})

    t0 = time.perf_counter()
    topo0 = read_grid_any(TOPO, nglv, 2 * nglv)
    t_topo = time.perf_counter() - t0
    t0 = time.perf_counter()
    love = LoveNumbers(EARTH, norder)
    t_love = time.perf_counter() - t0
    t0 = time.perf_counter()
    opts = parse_opts(a.opt)
    sht = make_sht(a.sht, nglv, 2 * nglv, norder, workers=a.threads, chunk=a.chunk,
                   raw="raw" in opts)
    if "numba" in opts:
        t_nb = time.perf_counter()
        from . import kernels

        kernels.warmup()
        print(f"[slevolve] numba warm-up: {time.perf_counter()-t_nb:.2f} s")
        log["numba_warmup_s"] = time.perf_counter() - t_nb
    t_sht = time.perf_counter() - t0
    print(
        f"[slevolve] setup: topo {t_topo:.2f}s  love {t_love:.2f}s  "
        f"SHT[{a.sht}] setup {t_sht:.1f}s (tables {sht.table_bytes/1e6:.0f} MB)"
    )
    log.update(
        {
            "setup_topo_s": t_topo,
            "setup_love_s": t_love,
            "setup_sht_s": t_sht,
            "sht_table_bytes": int(sht.table_bytes),
            "love_max_modes": love.max_modes,
        }
    )

    t0 = time.perf_counter()
    ice0 = np.ascontiguousarray(ice[0])
    acc = a.inner != "off"
    cls = OptSeaLevelSolver if (opts or acc) else SeaLevelSolver
    extra = {"opt": opts} if opts else {}
    if acc:
        extra.update({"opt": opts, "inner": a.inner, "anderson_depth": a.anderson_depth,
                      "anderson_freeze": a.anderson_freeze, "anderson_keep": a.anderson_keep})
    if a.rhoi is not None or a.radius is not None:
        extra["planet"] = make_planet(a.rhoi, a.radius)
    slv = cls(
        sht,
        love,
        topo0,
        ice0,
        norder=norder,
        dt=a.dt if times is None else float(dts[0]),
        tpw=tpw,
        checkmarine=cm,
        calcRG=rg,
        epsilon1=a.epsilon1,
        mem_threads=a.mem_threads,
        ramp=a.ramp,
        **extra,
    )
    t_init = time.perf_counter() - t0
    if a.horizontal or a.horizontal_grid:
        slv.horizontal = True
        log["horizontal"] = "grid" if a.horizontal_grid else "coef"
    log["solver_init_s"] = t_init
    print(f"[slevolve] sl_solver_init: {t_init:.2f} s")

    # step-0 outputs
    def save_mask(path, m):
        if a.io_compat == "fortran":
            np.save(path, np.ascontiguousarray(m, dtype="<f8"))
        else:
            np.save(path, np.asarray(m).astype(np.uint8))

    np.save(os.path.join(out, "tgrid0.npy"), slv.tinit_0)
    save_mask(os.path.join(out, "ocean0.npy"), slv.cxy0)
    save_mask(os.path.join(out, "beta0.npy"), slv.beta0)
    scal0 = {
        "ocean_area": slv.ocean_area0,
        "oceanBeta_area": slv.oceanBeta_area0,
        "ice_volume": slv.ice_volume0,
    }

    numiter = []
    walls = []
    scalars = [scal0]
    prof_tot = {}

    if a.cprofile:
        import cProfile
        import pstats

        pr = cProfile.Profile()
        pr.enable()

    nsht = {"analysis": 0, "synthesis": 0}
    if a.count_sht:
        for _nm in ("analysis", "synthesis"):
            _f = getattr(sht, _nm)

            def _wrap(*aa, _f=_f, _nm=_nm, **kk):
                nsht[_nm] += 1
                return _f(*aa, **kk)

            setattr(sht, _nm, _wrap)
    want_rates = bool(fields & {"dSLdt", "dRdt"})
    t_run = time.perf_counter()
    for k in range(1, a.steps + 1):
        icek = np.ascontiguousarray(ice[k])
        r = slv.step(icek, profile=a.profile, dt=None if dts is None else float(dts[k - 1]),
                     rates=want_rates)
        numiter.append(r["numiter"])
        walls.append(r["wall"])
        if a.profile:
            for kk, vv in r["prof"].items():
                prof_tot[kk] = prof_tot.get(kk, 0.0) + vv
        if "tgrid" in fields:
            np.save(os.path.join(out, f"tgrid{k}.npy"), r["tgrid"])
        if "ocean" in fields:
            save_mask(os.path.join(out, f"ocean{k}.npy"), r["ocean"])
        if "beta" in fields:
            save_mask(os.path.join(out, f"beta{k}.npy"), r["beta"])
        if "R" in fields and r["R"] is not None:
            np.save(os.path.join(out, f"R{k}.npy"), r["R"])
        if "G" in fields and r["G"] is not None:
            np.save(os.path.join(out, f"G{k}.npy"), r["G"])
        if "dS" in fields:
            if a.io_compat == "fortran":
                npyio.save_spectral(os.path.join(out, f"dS_converged{k}.npy"), r["dS_converged"])
            else:
                np.save(os.path.join(out, f"dS_converged{k}.npy"), r["dS_converged"])
        if r.get("rates") is not None:
            npyio.save_rates(out, f"{k}", r["rates"], fields)
        if r.get("horiz") is not None:
            npyio.save_spectral(os.path.join(out, f"horiz{k}.npy"), r["horiz"])
            if a.horizontal_grid:
                from .horiz import gradient

                ut, up = gradient(r["horiz"], nglv, 2 * nglv, nthreads=a.threads)
                np.save(os.path.join(out, f"utheta{k}.npy"), ut)
                np.save(os.path.join(out, f"uphi{k}.npy"), up)
        scalars.append(
            {
                kk: r[kk]
                for kk in ("ocean_area", "oceanBeta_area", "gmslc_ocn", "gmslc_ocnBeta", "ice_volume")
            }
        )
        if k % 10 == 0 or k <= 3:
            print(
                f"[slevolve] step {k:4d}  numiter={r['numiter']}  {r['wall']:.3f} s"
                f"   load={loadavg()[0]:.2f}",
                flush=True,
            )
    t_run = time.perf_counter() - t_run

    if a.cprofile:
        pr.disable()
        with open(os.path.join(out, "cprofile.txt"), "w") as fh:
            st = pstats.Stats(pr, stream=fh)
            st.sort_stats("cumulative").print_stats(40)
        print(f"[slevolve] cProfile written to {out}/cprofile.txt")

    np.save(os.path.join(out, "numiter.npy"), np.asarray(numiter, dtype=np.int32))
    np.save(os.path.join(out, "step_wall_seconds.npy"), np.asarray(walls))
    log.update(
        {
            "total_solver_seconds": t_run,
            "per_step_mean": float(np.mean(walls)) if walls else None,
            "per_step_median": float(np.median(walls)) if walls else None,
            "per_step_first10": [float(x) for x in walls[:10]],
            "per_step_last10": [float(x) for x in walls[-10:]],
            "numiter_total": int(np.sum(numiter)),
            "numiter_mean": float(np.mean(numiter)) if numiter else None,
            "inner_stats": getattr(slv, "inner_stats", None),
            "sht_counts": nsht if a.count_sht else None,
            "profile_seconds": prof_tot,
            "uptime_end": uptime(),
            "loadavg_end": loadavg(),
            "scalars": scalars,
        }
    )
    with open(os.path.join(out, "run_log.json"), "w") as fh:
        json.dump(log, fh, indent=1)
    print(
        f"[slevolve] DONE {a.steps} steps in {t_run:.1f} s "
        f"({t_run/max(a.steps,1):.3f} s/step), numiter total {sum(numiter)}"
    )
    print(f"[slevolve] {log['uptime_end']}")
    if a.profile:
        tot = prof_tot.get("total", t_run)
        for kk in sorted(prof_tot, key=lambda x: -prof_tot[x]):
            print(f"[slevolve] profile {kk:10s} {prof_tot[kk]:8.2f} s  {100*prof_tot[kk]/tot:5.1f} %")


if __name__ == "__main__":
    main()
