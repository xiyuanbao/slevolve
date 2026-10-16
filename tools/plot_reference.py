#!/usr/bin/env python3
"""plot_reference.py -- maps of the shipped reference results (reference/<case>/*.npy, <case>.json).

    python tools/plot_reference.py --outdir docs/src/img --topo1 T1 --topo2 T2 --topo3 T3   # all cases
    python tools/plot_reference.py --outdir figs --case case1 --topo1 $INPUT1/others/topo0_archive
    python tools/plot_reference.py --outdir figs --case synthetic       # needs no input data

Writes
  ref_case1_dtopo550.png   Case 1: topography change at step 550 (south polar stereographic + global)
  ref_case2_dtopo1200.png  Case 2: topography change at step 1200 of the pass from T_obs (global)
  ref_case3_dtopo50.png    Case 3: topography change at step 50 of pass 5 (global)
  ref_case3_misfit.png     Case 3: outer misfit (area-weighted RMS of T_pred - T_obs) per pass
  synthetic_test.png       synthetic Earth of the Fortran tests: inputs and tgrid8 - tgrid0
                           (reference/synthetic/tgrid8.npy)

Maps are drawn with python/tools/plot_fields.py (Mollweide; south polar stereographic).  The coast is
the 0 m contour of the topography given with --topo1 (Case 1: the initial topography topo0_archive),
--topo2 (Case 2: T_obs = etopo2_512_ice6gC) and --topo3 (Case 3: topo_pres); .npy, .bin or legacy text
grids are accepted.  Without it the map is drawn without a coast and a warning is printed.  To plot
your own run in the same way, build the same field first, e.g. for Case 1 dtopo_550 = tgrid550 -
tgrid0 (python/tools/make_reference.py --case case1 --dir OUT --refdir MYREF writes it), then pass
--refdir MYREF.  Needs numpy and matplotlib.
"""
import argparse
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..")
sys.path.insert(0, os.path.join(ROOT, "python"))
sys.path.insert(0, HERE)
from slevolve.tools.plot_fields import CHANGE_CMAP, change_colorbar, mollweide_map, nice_limit, south_polar_map  # noqa: E402


def load_topo(path, shape):
    if not path:
        print("warning: no topography given, map drawn without coastline")
        return np.full(shape, -1.0)
    from slevolve.npyio import read_grid
    base = path[:-4] if path.endswith(".npy") else path
    return read_grid(base if os.path.exists(base + ".npy") or not os.path.exists(path) else path,
                     nlat=shape[0], nlon=shape[1])


def report(name, f, lim):
    print(f"{name}: colour limits -{lim:g}..{lim:g} m (99th percentile of |field|: "
          f"{np.percentile(np.abs(f), 99.0):.1f} m; max |field| {np.abs(f).max():.1f} m)")


def global_fig(plt, f, topo, title, path, vmax, label="topography change (m)"):
    fig = plt.figure(figsize=(7.4, 3.6), constrained_layout=True)
    ax = fig.add_subplot(1, 1, 1, projection="mollweide")
    pc = mollweide_map(ax, f, topo, cmap=CHANGE_CMAP, vmin=-vmax, vmax=vmax)
    ax.set_title(title, fontsize=9)
    change_colorbar(fig, pc, vmax, label=label, ax=ax, shrink=0.8)
    fig.savefig(path, dpi=110)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description="maps of the reference results")
    ap.add_argument("--outdir", required=True, help="directory for the PNG files")
    ap.add_argument("--refdir", default=os.path.join(ROOT, "reference"), help="reference directory (default reference/)")
    ap.add_argument("--case", choices=("case1", "case2", "case3", "synthetic", "all"), default="all")
    ap.add_argument("--topo1", help="Case 1 initial topography (coastline)")
    ap.add_argument("--topo2", help="Case 2 present topography T_obs (coastline)")
    ap.add_argument("--topo3", help="Case 3 present topography topo_pres (coastline)")
    a = ap.parse_args()
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    os.makedirs(a.outdir, exist_ok=True)
    cases = ("case1", "case2", "case3", "synthetic") if a.case == "all" else (a.case,)
    R = a.refdir
    if "case1" in cases:
        f = np.load(os.path.join(R, "case1", "dtopo_550.npy"))
        t = load_topo(a.topo1, f.shape)
        fig = plt.figure(figsize=(10.0, 3.9), constrained_layout=True)
        gs = fig.add_gridspec(1, 3, width_ratios=[1.0, 1.55, 0.04])
        ax = fig.add_subplot(gs[0])
        lim = 400.0                       # as Fig. 3a of the paper
        report("case1 dtopo_550", f, lim)
        pc = south_polar_map(ax, f, t, latmax=-60.0, cmap=CHANGE_CMAP, vmin=-lim, vmax=lim)
        ax.set_title("south polar view (60 S to the pole)", fontsize=9)
        ax = fig.add_subplot(gs[1], projection="mollweide")
        mollweide_map(ax, f, t, cmap=CHANGE_CMAP, vmin=-lim, vmax=lim)
        ax.set_title("global", fontsize=9)
        change_colorbar(fig, pc, lim, label="topography change (m)", cax=fig.add_subplot(gs[2]))
        fig.suptitle("Case 1: topography change, step 550 (m); black: initial coast", fontsize=10)
        fig.savefig(os.path.join(a.outdir, "ref_case1_dtopo550.png"), dpi=110)
        plt.close(fig)
    if "case2" in cases:
        f = np.load(os.path.join(R, "case2", "dtopo_1200.npy"))
        report("case2 dtopo_1200", f, nice_limit(f))
        global_fig(plt, f, load_topo(a.topo2, f.shape),
                   "Case 2: topography change, step 1200 of the pass from T_obs (m); black: present coast",
                   os.path.join(a.outdir, "ref_case2_dtopo1200.png"), nice_limit(f))
    if "case3" in cases:
        f = np.load(os.path.join(R, "case3", "dtopo_p5_s50.npy"))
        report("case3 dtopo_p5_s50", f, nice_limit(f))
        global_fig(plt, f, load_topo(a.topo3, f.shape),
                   "Case 3: topography change, step 50 of pass 5 (m); black: present coast",
                   os.path.join(a.outdir, "ref_case3_dtopo50.png"), nice_limit(f))
        ref = json.load(open(os.path.join(R, "case3.json")))
        m = ref["misfit_rms"]
        fig, ax = plt.subplots(figsize=(4.6, 3.0), constrained_layout=True)
        ax.bar(np.arange(1, len(m) + 1), m, color="C0")
        ax.set_yscale("log")
        for k, v in enumerate(m, 1):
            ax.text(k, v * 1.15, f"{v:.3g}", ha="center", fontsize=8)
        ax.set_xlabel("outer pass")
        ax.set_ylabel("RMS of T_pred - T_obs (m)")
        ax.set_title("Case 3: outer misfit per pass", fontsize=10)
        ax.set_ylim(min(m) / 2, max(m) * 3)
        fig.savefig(os.path.join(a.outdir, "ref_case3_misfit.png"), dpi=110)
        plt.close(fig)
    if "synthetic" in cases:
        from make_synthetic_case import synthetic_earth
        tg8 = np.load(os.path.join(R, "synthetic", "tgrid8.npy"))
        _, _, topo, ice = synthetic_earth(tg8.shape[0])
        d = tg8 - topo
        nan0 = np.where(ice[0] > 0, ice[0], np.nan)
        fig = plt.figure(figsize=(11.0, 6.4), constrained_layout=True)
        ax = fig.add_subplot(2, 2, 1, projection="mollweide")
        pc = mollweide_map(ax, topo, topo, cmap="terrain", vmin=-5000, vmax=2500)
        ax.set_title("synthetic topography tgrid0 (m), 512 x 1024; black: coast", fontsize=9)
        fig.colorbar(pc, ax=ax, shrink=0.8)
        ax = fig.add_subplot(2, 2, 2, projection="mollweide")
        pc = mollweide_map(ax, nan0, topo, ice=ice[0], cmap="Blues", vmin=0, vmax=3000)
        ax.set_title("ice at step 0 (m); dashed: ice margin", fontsize=9)
        fig.colorbar(pc, ax=ax, shrink=0.8)
        ax = fig.add_subplot(2, 2, 3, projection="mollweide")
        lim = nice_limit(d)
        report("synthetic tgrid8 - tgrid0", d, lim)
        pc = mollweide_map(ax, d, topo, cmap=CHANGE_CMAP, vmin=-lim, vmax=lim)
        ax.set_title("tgrid8 - tgrid0 (m), 8 steps of 1 yr (reference/synthetic/tgrid8.npy)", fontsize=9)
        change_colorbar(fig, pc, lim, label="topography change (m)", ax=ax, shrink=0.8)
        ax = fig.add_subplot(2, 2, 4)
        pc = south_polar_map(ax, d, topo, ice=ice[0], latmax=-55.0, cmap=CHANGE_CMAP, vmin=-lim, vmax=lim)
        ax.set_title("tgrid8 - tgrid0, south polar view (m)", fontsize=9)
        change_colorbar(fig, pc, lim, label="topography change (m)", ax=ax, shrink=0.8)
        fig.savefig(os.path.join(a.outdir, "synthetic_test.png"), dpi=110)
        plt.close(fig)
    print("wrote figures to", a.outdir)


if __name__ == "__main__":
    main()
