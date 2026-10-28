#!/usr/bin/env python3
"""
run_outer.py -- outer (initial-topography) iteration of Han et al. (2022) around the
Python solver, with an optional per-cell secant update.

The legacy Fortran (sl_solver_init, initial_topo
= .false.) does, for outer iteration itersl:
    itersl == 1 : tinit_0 = truetopo                                   (eq. 48)
    itersl  > 1 : tinit_0 = tgrid0_<itersl-1> + (truetopo - pred_pres_topo_<itersl-1>)
and at the last step (current_time == L_sim) writes pred_pres_topo_<itersl> =
topoxy and tgrid0_<itersl> = tinit_0.  --outer picard reproduces exactly that
(in memory; the files are also written as .npy).

--outer secant (from outer iteration 3 on), per cell, with r_k = pred_k - obs:
    s      = (r_k - r_{k-1}) / (T0_k - T0_{k-1})   where |T0_k - T0_{k-1}| > --secant-tol
    s      = 1 (Picard) where the denominator is tiny or s is not finite
    s      = clip(s, 0.5, 1.5)                      (Picard is s = 1)
    T0_{k+1} = T0_k - r_k / s

Example (ICE-NADI-1, L = 511):
    python3 run_outer.py --iters 5 --inputs INPUT_FILES --lmax 511 \
        --earth INPUT_FILES/earthmodel/prem_512.l96C.ump3.lm1 --out OUT
"""
import argparse
import os
import sys


def _args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outer", default="picard", choices=("picard", "secant"), help="outer update of the initial topography: picard (legacy) or per-cell secant from pass 3")
    ap.add_argument("--iters", type=int, default=5, help="number of outer passes (default 5)")
    ap.add_argument("--inputs", required=True,
                    help="INPUT_FILES directory (icemodel/, others/)")
    ap.add_argument("--icepfx", default="nadi_", help="ice file prefix (default nadi_)")
    ap.add_argument("--nfiles", type=int, default=51, help="number of ice files 0..nfiles-1 = steps + 1 (default 51)")
    ap.add_argument("--topo-pres", default=None, help="default <inputs>/others/topo_pres")
    ap.add_argument("--topo0", default=None,
                    help="initial guess for T0 of outer iteration 1 (text grid or .npy); default topo_pres")
    ap.add_argument("--earth", required=True, help="Love-number file (no default)")
    ap.add_argument("--lmax", type=int, required=True, help="spectral truncation (511 for the L=511 runs)")
    ap.add_argument("--nglv", type=int, default=512, help="grid latitudes (default 512)")
    ap.add_argument("--dt", type=float, default=500.0, help="step length in YEARS (one fixed value; time spacing of the ice files)")
    ap.add_argument("--times", default=None,
                    help="optional times file: one time t_k in years per line for the ice files (or --ice-seq "
                         "entries) 0..N, strictly increasing; step k has length t_k - t_(k-1) and --dt is ignored")
    ap.add_argument("--rhoi", type=float, default=920.0, help="ice density (default 920)")
    ap.add_argument("--radius", type=float, default=6.371e6, help="Earth radius in m (default 6.371e6)")
    ap.add_argument("--tpw", default="T", help="T|F rotational feedback (default T)")
    ap.add_argument("--checkmarine", default="T", help="T|F marine check (default T)")
    ap.add_argument("--epsilon1", type=float, default=1e-5, help="inner-iteration tolerance (default 1e-5)")
    ap.add_argument("--threads", type=int, default=4,
                    help="ducc and numba threads (NUMBA_NUM_THREADS); with --opt numba the state advance "
                         "uses the numba thread pool, so there is no separate --mem-threads")
    ap.add_argument("--opt", default="raw,diag,dsl,numba,sums", help="optimisations, as run_ais.py --opt (default all)")
    ap.add_argument("--inner", default="off", choices=("off", "extrap", "anderson"), help="inner-iteration start: off, extrap (recommended) or anderson (default off)")
    ap.add_argument("--secant-tol", type=float, default=1e-3, help="m; smaller |dT0| -> Picard")
    ap.add_argument("--smin", type=float, default=0.5, help="with --outer secant: lower bound of the secant factor (default 0.5)")
    ap.add_argument("--smax", type=float, default=1.5, help="with --outer secant: upper bound of the secant factor (default 1.5)")
    ap.add_argument("--save-steps", default="", help="comma list of steps whose tgrid to save")
    ap.add_argument("--ramp", action="store_true",
                    help="second-order load placement: each load increment accrues linearly over its "
                         "step (default off: step placement of Kendall et al. 2005)")
    ap.add_argument("--rates", action="store_true",
                    help="write the analytic rates dSLdt<N>_it<k>.npy, dRdt<N>_it<k>.npy (m/yr) and "
                         "TPWrate<N>_it<k>.txt at --save-steps and the last step; output only")
    ap.add_argument("--horizontal", action="store_true",
                    help="also write the horizontal-displacement scalar Psi_lm (horiz<N>_it<k>.npy at --save-steps and the last step, "
                         "u_h = grad_h Psi); output only, default off")
    ap.add_argument("--horizontal-grid", action="store_true",
                    help="with --horizontal: also write utheta<N>.npy / uphi<N>.npy (ducc0)")
    ap.add_argument("--reuse", default="off", choices=("off", "warm", "warm+ice", "warm+corr", "warm+corr2"),
                    help="warm start across outer passes: warm = start each step from the converged load increment and ocean mask of the previous pass (recommended for --iters > 1); warm+ice, warm+corr, warm+corr2 experimental; default off")
    ap.add_argument("--count-sht", action="store_true", help="count transforms per outer iteration")
    ap.add_argument("--ice-seq", default="",
                    help="custom ice history built from the cached stack: comma list of 'i' or 'i:j:w' "
                         "((1-w) I_i + w I_j); '@FILE' reads the list from FILE. Default: files 0..nfiles-1")
    ap.add_argument("--rebuild-cache", action="store_true",
                    help="rebuild the ice cache even if its manifest (directory, prefix, file count, "
                         "size and mtime of every ice file) matches")
    ap.add_argument("--out", required=True, help="output directory")
    return ap.parse_args()


A = _args()
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS",
           "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ[_v] = str(A.threads)

import json  # noqa: E402
import time  # noqa: E402

import numpy as np  # noqa: E402

from .sht_ducc import make_sht  # noqa: E402
from .solver import LoveNumbers, Planet, read_grid_text  # noqa: E402
from .solver_opt import OptSeaLevelSolver, parse_opts  # noqa: E402
from . import npyio  # noqa: E402
from . import _paths  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))


def planet(rhoi, radius):
    p = Planet()
    p.rhoi = rhoi
    p.radius = radius
    p.moiA = 0.3296145 * p.mass * p.radius**2
    p.moiC = 0.3307007 * p.mass * p.radius**2
    return p


def ice_stack(a):
    nl, no = a.nglv, 2 * a.nglv
    tag = os.path.basename(os.path.dirname(os.path.normpath(a.inputs))) or "case"
    path = os.path.join(_paths.CACHE, f"ice_{tag}_{a.icepfx}0_{a.nfiles-1}.npy")
    icedir = os.path.join(a.inputs, "icemodel")
    man = npyio.ice_manifest(icedir, a.icepfx, a.nfiles)
    ok, why = npyio.cache_matches(path, man)
    if a.rebuild_cache or not ok:
        if a.rebuild_cache:
            print(f"[outer] --rebuild-cache: rebuilding {path}", flush=True)
        elif os.path.exists(path):
            print(f"[outer] ice cache out of date ({why}): rebuilding {path}", flush=True)
        os.makedirs(os.path.dirname(path), exist_ok=True)
        arr = np.lib.format.open_memmap(path + ".tmp", "w+", np.float64, (a.nfiles, nl, no))
        for k in range(a.nfiles):
            # <name>.npy / <name>.bin next to the text file are preferred (npyio)
            arr[k] = npyio.find_grid(icedir, f"{a.icepfx}{k}", nl, no)
        arr.flush()
        del arr
        os.replace(path + ".tmp", path)
        npyio.write_manifest(path, man)
    return np.load(path, mmap_mode="r")


def aw_rms(f, w):
    W = np.broadcast_to(w[:, None], f.shape)
    return float(np.sqrt((W * f**2).sum() / W.sum()))


def main():
    a = A
    os.makedirs(a.out, exist_ok=True)
    print(f"[outer] step length dt = {a.dt:g} yr, ice files {a.icepfx}0..{a.nfiles - 1}", flush=True)
    nl, no = a.nglv, 2 * a.nglv
    obs = npyio.read_grid(a.topo_pres or os.path.join(a.inputs, "others", "topo_pres"), nl, no)
    earth = a.earth
    ice = ice_stack(a)
    if a.ice_seq:
        spec = open(a.ice_seq[1:]).read() if a.ice_seq.startswith("@") else a.ice_seq
        items = [x.strip() for x in spec.replace("\n", ",").split(",") if x.strip()]

        class _Seq:
            def __init__(self, base, items):
                self.b, self.it = base, items

            def __len__(self):
                return len(self.it)

            def __getitem__(self, k):
                p = self.it[k].split(":")
                if len(p) == 1:
                    return self.b[int(p[0])]
                i, j, w = int(p[0]), int(p[1]), float(p[2])
                return (1.0 - w) * self.b[i] + w * self.b[j]

        ice = _Seq(ice, items)
    nsteps = len(ice) - 1
    dts = None
    if a.times:
        times = npyio.read_times(a.times)
        if times.size != nsteps + 1:
            raise SystemExit(f"--times {a.times}: {times.size} times for {nsteps + 1} ice states")
        dts = np.diff(times)
        print(f"[outer] times file {a.times}: t0={times[0]:g} yr, dt {dts.min():g}-{dts.max():g} yr "
              "(--dt ignored)", flush=True)
    love = LoveNumbers(earth, a.lmax)
    opts = parse_opts(a.opt)
    sht = make_sht("ducc", nl, no, a.lmax, workers=a.threads, raw="raw" in opts)
    if "numba" in opts:
        from . import kernels
        kernels.warmup()
    _, gw = np.polynomial.legendre.leggauss(nl)
    gw = gw[::-1]
    save_steps = set(int(x) for x in a.save_steps.split(",") if x)
    tpw = a.tpw.upper().startswith("T")
    cm = a.checkmarine.upper().startswith("T")

    T0 = obs.copy()  # itersl = 1  (eq. 48)
    if a.topo0:
        T0 = npyio.read_grid(a.topo0, nl, no)
    nsht = {"analysis": 0, "synthesis": 0}
    if a.count_sht:
        for _nm in ("analysis", "synthesis"):
            _f = getattr(sht, _nm)

            def _wrap(*aa, _f=_f, _nm=_nm, **kk):
                nsht[_nm] += 1
                return _f(*aa, **kk)

            setattr(sht, _nm, _wrap)
    store = {}  # --reuse: per-step converged state, survives across outer iterations
    hist = []  # (T0_k, r_k)
    log = {"argv": sys.argv, "outer": a.outer, "iters": [], "ramp": bool(a.ramp),
           "times_years": [float(x) for x in times] if dts is not None else None}
    for it in range(1, a.iters + 1):
        t0 = time.perf_counter()
        slv = OptSeaLevelSolver(
            sht, love, T0, np.ascontiguousarray(ice[0]), norder=a.lmax, dt=a.dt if dts is None else float(dts[0]), tpw=tpw,
            checkmarine=cm, calcRG=True, epsilon1=a.epsilon1, planet=planet(a.rhoi, a.radius),
            opt=opts, inner=a.inner, reuse=a.reuse, reuse_store=store, ramp=a.ramp)
        slv.horizontal = bool(a.horizontal or a.horizontal_grid)
        nit = 0
        nits = []
        n0 = dict(nsht)
        for k in range(1, nsteps + 1):
            want = a.rates and (k in save_steps or k == nsteps)
            r = slv.step(np.ascontiguousarray(ice[k]), dt=None if dts is None else float(dts[k - 1]),
                         rates=want)
            if want:
                npyio.save_rates(a.out, f"{k}_it{it}", r["rates"])
            nit += r["numiter"]
            nits.append(r["numiter"])
            if k in save_steps:
                np.save(os.path.join(a.out, f"tgrid{k}_it{it}.npy"), r["tgrid"])
            if r.get("horiz") is not None and (k in save_steps or k == nsteps):
                npyio.save_spectral(os.path.join(a.out, f"horiz{k}_it{it}.npy"), r["horiz"])
                if a.horizontal_grid:
                    from .horiz import gradient

                    ut, up = gradient(r["horiz"], nl, no, nthreads=a.threads)
                    np.save(os.path.join(a.out, f"utheta{k}_it{it}.npy"), ut)
                    np.save(os.path.join(a.out, f"uphi{k}_it{it}.npy"), up)
        pred = r["tgrid"].copy()
        wall = time.perf_counter() - t0
        np.save(os.path.join(a.out, f"pred_pres_topo_{it}.npy"), pred)
        np.save(os.path.join(a.out, f"tgrid0_{it}.npy"), T0)
        res = pred - obs
        rec = {"it": it, "rms_pred_minus_obs": aw_rms(res, gw),
               "max_abs": float(np.abs(res).max()), "rms_T0_minus_obs": aw_rms(T0 - obs, gw),
               "numiter": nit, "wall_s": wall, "numiter_steps": nits}
        if a.count_sht:
            rec["sht"] = {kk: nsht[kk] - n0[kk] for kk in nsht}
        if a.reuse != "off":
            rec["reuse_stats"] = dict(slv.reuse_stats)
            rec["store_bytes"] = int(sum(v.nbytes for e in store.values() for v in e.values()))
        # ---- next initial topography ------------------------------------
        hist.append((T0, res))
        if a.outer == "secant" and len(hist) >= 2:
            (T0p, rp), (T0c, rc) = hist[-2], hist[-1]
            dT = T0c - T0p
            ok = np.abs(dT) > a.secant_tol
            s = np.ones_like(dT)
            with np.errstate(divide="ignore", invalid="ignore"):
                s[ok] = (rc[ok] - rp[ok]) / dT[ok]
            s[~np.isfinite(s)] = 1.0
            ncl = int(((s < a.smin) | (s > a.smax)).sum())
            s = np.clip(s, a.smin, a.smax)
            T0 = T0c - rc / s
            rec.update({"secant_cells_picard": int((~ok).sum()), "secant_cells_clamped": ncl,
                        "secant_s_mean": float(s[ok].mean()) if ok.any() else 1.0})
        else:
            T0 = T0 - res  # Picard: tgrid0 + (truetopo - pred)
        hist = hist[-2:]
        log["iters"].append(rec)
        print(json.dumps({kk: vv for kk, vv in rec.items() if kk != "numiter_steps"}), flush=True)
    np.save(os.path.join(a.out, "T0_next.npy"), T0)
    with open(os.path.join(a.out, "outer_log.json"), "w") as fh:
        json.dump(log, fh, indent=1)


if __name__ == "__main__":
    main()
