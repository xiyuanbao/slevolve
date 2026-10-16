#!/usr/bin/env python3
"""Reference quantities of the three benchmark cases from a Python output directory.

    python python/tools/make_reference.py --case case1 --dir OUT_case1_py --refdir reference
    python python/tools/make_reference.py --case case2 --dir OUT_case2_py --refdir reference
    python python/tools/make_reference.py --case case3 --dir OUT_case3_py --refdir reference \
        --topo-pres INPUT_FILES/others/topo_pres

Writes reference/<case>.json (scalar summaries, site values, inner iterations, Case 3 misfits) and
full-resolution fields reference/<case>/<name>.npy (512 x 1024, '<f8'):
  case1: dtopo_<k>.npy = tgrid<k> - tgrid0 and dsl_<k>.npy = -(tgrid<k> - tgrid0) on ocean<k> (0 on land)
         for k = 1, 100, 550 (run_ais.py output, --save-fields tgrid,ocean)
  case2: dtopo_<k>.npy = tgrid<k>_it1 - tgrid0_1 for k = 1, 600, 1200 (run_outer.py --iters 1
         --save-steps 1,600,1200)
  case3: dtopo_p1_s50.npy, dtopo_p<N>_s50.npy = tgrid50_it<p> - tgrid0_<p> for the first and the final
         pass, and T0_minus_obs.npy = tgrid0_<N> - T_obs (converged initial topography of the final pass)
Relative sea-level change is -(topography change); it is given only where the ocean mask is saved.
"""
import argparse
import json
import os

import numpy as np

NLAT, NLON = 512, 1024
SITES = [  # name, lat (deg N), lon (deg E, -180..180)
    ("Hudson Bay centre", 60.0, -85.0), ("James Bay", 53.0, -80.5), ("Barbados", 13.1, -59.6),
    ("Bermuda", 32.3, -64.8), ("Fennoscandia centre (Gulf of Bothnia)", 63.0, 20.0),
    ("North Sea", 56.0, 3.0), ("Ross Ice Shelf front", -77.5, 180.0),
    ("Ronne Ice Shelf front", -75.5, -55.0), ("Amundsen Sea sector", -73.0, -112.0),
    ("Thwaites grounding zone", -75.5, -106.5), ("Greenland centre", 72.5, -40.0),
    ("Greenland SE coast", 65.6, -37.6), ("Florida", 27.5, -81.0), ("Patagonia", -50.0, -72.0),
    ("Tahiti", -17.6, -149.5), ("Japan (Tokyo Bay)", 35.6, 139.8), ("Antarctic Peninsula", -65.0, -63.0),
    ("Weddell Sea", -70.0, -40.0), ("Hudson Strait", 62.5, -70.0), ("Barents Sea", 75.0, 40.0)]
STEPS = {"case1": [1, 100, 550], "case2": [1, 600, 1200]}


def grid():
    x, w = np.polynomial.legendre.leggauss(NLAT)
    lat = np.degrees(np.arcsin(x[::-1]))  # row 0 = northernmost (GLlat_512.txt order)
    return lat, np.arange(NLON) * 360.0 / NLON, w[::-1]


def site_cells():
    lat, lon, _ = grid()
    out = []
    for name, la, lo in SITES:
        i = int(np.argmin(np.abs(lat - la)))
        j = int(np.argmin(np.abs(((lon - lo % 360.0) + 180.0) % 360.0 - 180.0)))
        out.append({"name": name, "lat": la, "lon": lo, "i": i, "j": j,
                    "cell_lat": round(float(lat[i]), 4), "cell_lon": round(float(lon[j]), 4)})
    return out


def stats(f, mask=None):
    _, _, w = grid()
    W = np.broadcast_to(w[:, None], f.shape)
    if mask is not None:
        W = W * mask
    return {"rms": float(np.sqrt((W * f * f).sum() / W.sum())),
            "maxabs": float(np.abs(f[W > 0]).max())}


def load(d, name):
    return np.load(os.path.join(d, name)).astype("<f8")


def fields(case, d, topo_pres=None):
    """{name: (field, mask or None)} and metadata, from an output directory."""
    F, meta = {}, {}
    if case == "case1":
        t0 = load(d, "tgrid0.npy")
        for k in STEPS[case]:
            if not os.path.exists(os.path.join(d, f"tgrid{k}.npy")):
                continue  # e.g. the Fortran example writes only tgrid550 (SLEVOLVE_SNAP=550)
            dt = load(d, f"tgrid{k}.npy") - t0
            F[f"dtopo_{k}"] = (dt, None)
            p = os.path.join(d, f"ocean{k}.npy")
            if os.path.exists(p):
                oc = (np.load(p) > 0.5).astype("<f8")
                F[f"dsl_{k}"] = (-dt * oc, oc)
        if os.path.exists(os.path.join(d, "numiter.npy")):  # Python run_ais.py
            meta["numiter"] = [int(x) for x in np.load(os.path.join(d, "numiter.npy"))]
        elif os.path.exists(os.path.join(d, "numiter")):  # Fortran OUTPUT_SLM/numiter (text)
            meta["numiter"] = [int(x) for x in open(os.path.join(d, "numiter")).read().split()]
    elif case == "case2":
        t0 = load(d, "tgrid0_1.npy")
        for k in STEPS[case]:
            F[f"dtopo_{k}"] = (load(d, f"tgrid{k}_it1.npy") - t0, None)
        log = json.load(open(os.path.join(d, "outer_log.json")))
        meta["numiter"] = log["iters"][0]["numiter_steps"]
        meta["misfit_rms"] = [log["iters"][0]["rms_pred_minus_obs"]]
    else:
        log = json.load(open(os.path.join(d, "outer_log.json")))
        n = len(log["iters"])
        for p in sorted({1, n}):
            F[f"dtopo_p{p}_s50"] = (load(d, f"tgrid50_it{p}.npy") - load(d, f"tgrid0_{p}.npy"), None)
        if topo_pres is None:
            raise SystemExit("case3 needs --topo-pres (observed present topography)")
        import sys
        from .. import npyio
        obs = npyio.read_grid(topo_pres, NLAT, NLON)
        F["T0_minus_obs"] = (load(d, f"tgrid0_{n}.npy") - obs, None)
        meta["numiter"] = [r["numiter_steps"] for r in log["iters"]]
        meta["misfit_rms"] = [r["rms_pred_minus_obs"] for r in log["iters"]]
        meta["misfit_max"] = [r["max_abs"] for r in log["iters"]]
    return F, meta


def summary(case, F, meta):
    cells = site_cells()
    out = {"case": case, "units": "m", "grid": "512x1024 Gauss-Legendre, row 0 north, lon 0..359.65 E",
           "sites": cells, "fields": {}}
    for name, (f, m) in F.items():
        s = stats(f, m)
        s["sites"] = [float(f[c["i"], c["j"]]) for c in cells]
        out["fields"][name] = s
    out.update(meta)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--case", required=True, choices=("case1", "case2", "case3"))
    ap.add_argument("--dir", required=True, help="output directory of the run")
    ap.add_argument("--refdir", required=True, help="reference directory (writes <case>.json, <case>/)")
    ap.add_argument("--topo-pres", default=None, help="case3: observed present topography")
    ap.add_argument("--note", default="", help="free text stored in the JSON (command, run time)")
    a = ap.parse_args()
    F, meta = fields(a.case, a.dir, a.topo_pres)
    out = summary(a.case, F, meta)
    out["note"] = a.note
    os.makedirs(os.path.join(a.refdir, a.case), exist_ok=True)
    for name, (f, _) in F.items():
        np.save(os.path.join(a.refdir, a.case, name + ".npy"), np.ascontiguousarray(f, dtype="<f8"))
    with open(os.path.join(a.refdir, a.case + ".json"), "w") as fh:
        json.dump(out, fh, indent=1)
    print(f"wrote {a.refdir}/{a.case}.json and {len(F)} fields in {a.refdir}/{a.case}/")


if __name__ == "__main__":
    main()
