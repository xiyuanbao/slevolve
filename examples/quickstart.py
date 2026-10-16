#!/usr/bin/env python3
"""Quick start (Python solver only): two melting ice sheets on a synthetic Earth, 128 x 256 grid.

The synthetic Earth of tools/make_synthetic_case.py (the Fortran tests use the same Earth at
512 x 1024): an Earth-like random topography (29 % land in continents of varied size with irregular
coasts, shelves at -100 to -200 m, ocean floor at -3700 to -4500 m with mid-ocean ridges), a polar
ice sheet (Vialov profile, 3000 m at the south pole, radius 17 degrees) and a mid-latitude ice sheet
(Vialov profile, 2500 m at 56 N 132 E, radius 12 degrees) that both melt to zero over 8 steps.
Here the steps are 1 kyr long so that the viscous response shows; maximum degree 63, Earth model
lovenum/prem_512.l96C.ump3.lm1 (first 63 degrees), rotational feedback on, marine (flotation) check
on.  A second run with the marine check off is made for comparison.
numpy only (scipy is used for the FFT if installed; matplotlib for the figures); a few seconds.

    python examples/quickstart.py [--outdir DIR] [--no-plot]

Prints the inner iterations per step and the relative sea-level change at three sites and, if
matplotlib is installed, writes quickstart_inputs.png, quickstart_topo.png, quickstart_map.png and
quickstart_series.png to --outdir (default: the current directory).  The maps are drawn with
python/tools/plot_fields.py (Mollweide and south polar stereographic; black line: initial coast).
"""
import argparse
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..")
sys.path.insert(0, os.path.join(ROOT, "python"))
sys.path.insert(0, os.path.join(ROOT, "tools"))
from slevolve.sht import SHT  # noqa: E402
from slevolve.solver import LoveNumbers, SeaLevelSolver  # noqa: E402
from make_synthetic_case import MIDLAT, POLAR, synthetic_earth  # noqa: E402

NLAT, NLON, LMAX, DT = 128, 256, 63, 1000.0
EARTH = os.path.join(ROOT, "lovenum", "prem_512.l96C.ump3.lm1")
# (name, latitude, longitude) of the three time-series sites
SITES = [("near field (82S 270E, under the polar ice)", -82.0, 270.0),
         ("intermediate (55S 0E, ocean)", -55.0, 0.0),
         ("far field (10N 180E, ocean)", 10.0, 180.0)]


def run(sht, love, topo, ice, checkmarine):
    slv = SeaLevelSolver(sht, love, topo, ice[0], norder=LMAX, dt=DT, tpw=True,
                         checkmarine=checkmarine, calcRG=True, epsilon1=1.0e-5)
    rsl, tg, nit = [np.zeros_like(topo)], [topo], []
    for k in range(1, len(ice)):
        r = slv.step(np.ascontiguousarray(ice[k]))
        tg.append(r["tgrid"].copy())
        rsl.append(-(r["tgrid"] - topo))  # relative sea-level change = -(topography change)
        nit.append(r["numiter"])
    return np.array(rsl), np.array(tg), nit


def main():
    ap = argparse.ArgumentParser(description="quick start: synthetic Earth, two melting ice sheets (Python solver)")
    ap.add_argument("--outdir", default=".", help="directory for the figures (default: current directory)")
    ap.add_argument("--no-plot", action="store_true", help="do not write figures")
    a = ap.parse_args()
    t0 = time.perf_counter()
    sht = SHT(NLAT, NLON, LMAX)
    love = LoveNumbers(EARTH, LMAX)
    lat, lon, topo, ice = synthetic_earth(NLAT)
    rsl, tg, nit = run(sht, love, topo, ice, True)
    rsl_nm, _, nit_nm = run(sht, love, topo, ice, False)
    t_run = time.perf_counter() - t0
    w = np.polynomial.legendre.leggauss(NLAT)[1][:, None] * np.ones((1, NLON))
    ocean = topo < 0
    idx = {name: (int(np.argmin(np.abs(lat - la))), int(np.argmin(np.abs(lon - lo)))) for name, la, lo in SITES}
    far = SITES[2][0]
    print(f"quickstart: synthetic Earth, grid {NLAT}x{NLON}, lmax {LMAX}, {len(nit)} steps of {DT:.0f} yr, "
          f"rotation on, marine check on")
    for k in range(1, len(nit) + 1):
        print(f"step {k}: ice {POLAR[k]:5.0f} m (pole) {MIDLAT[k]:5.0f} m (mid-lat)  inner iterations {nit[k-1]}  "
              f"ocean-mean RSL change {(w * rsl[k])[ocean].sum() / w[ocean].sum():+7.3f} m  "
              f"far field {rsl[k][idx[far]]:+7.3f} m")
    print(f"inner iterations per step: {nit} (marine check off: {nit_nm})")
    for name, ij in idx.items():
        print(f"final RSL change, {name}: {rsl[-1][ij]:+9.3f} m (marine check off {rsl_nm[-1][ij]:+9.3f} m)")
    print(f"quickstart: {t_run:.1f} s")
    if a.no_plot:
        return
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not available: no figures written")
        return
    from slevolve.tools.plot_fields import CHANGE_CMAP, change_colorbar, mollweide_map, nice_limit, south_polar_map
    os.makedirs(a.outdir, exist_ok=True)
    t = np.arange(len(nit) + 1) * DT / 1000.0
    nan0 = lambda f: np.where(f > 0, f, np.nan)  # noqa: E731

    def mark_sites(ax, polar):
        for name, ij in idx.items():
            la, lo = lat[ij[0]], lon[ij[1]]
            if polar:
                if la < -50:
                    from slevolve.tools.plot_fields import stereo_xy
                    x, y = stereo_xy(np.array([la]), np.array([lo]))
                    ax.plot(x, y, "o", mfc="none", mec="m", ms=6)
            else:
                ax.plot(np.radians(lo - 360.0 if lo >= 180.0 else lo), np.radians(la), "o", mfc="none", mec="m", ms=5)

    # inputs: topography, ice at the start, melt schedule
    fig = plt.figure(figsize=(11.0, 6.4), constrained_layout=True)
    ax = fig.add_subplot(2, 2, 1, projection="mollweide")
    pc = mollweide_map(ax, topo, topo, cmap="terrain", vmin=-5000, vmax=2500)
    ax.set_title("initial topography (m); black: coast", fontsize=9)
    fig.colorbar(pc, ax=ax, shrink=0.8)
    ax = fig.add_subplot(2, 2, 2, projection="mollweide")
    pc = mollweide_map(ax, nan0(ice[0]), topo, cmap="Blues", vmin=0, vmax=3000)
    ax.set_title("ice thickness at step 0 (m); the ice is gone at step 8", fontsize=9)
    fig.colorbar(pc, ax=ax, shrink=0.8)
    ax = fig.add_subplot(2, 2, 3)
    pc = south_polar_map(ax, nan0(ice[0]), topo, cmap="Blues", vmin=0, vmax=3000, latmax=-55.0)
    ax.set_title("polar ice sheet at step 0 (m), south polar view", fontsize=9)
    fig.colorbar(pc, ax=ax, shrink=0.8)
    ax = fig.add_subplot(2, 2, 4)
    ax.plot(t, POLAR, "o-", ms=3, label="polar sheet, at the pole")
    ax.plot(t, MIDLAT, "s-", ms=3, label="mid-latitude sheet, at its centre")
    ax.set_xlabel("time (kyr)")
    ax.set_ylabel("ice thickness (m)")
    ax.legend(fontsize=8)
    ax.set_title("melt schedule (Vialov profiles H (1 - (r / R)^(4/3))^(3/8))", fontsize=9)
    fig.savefig(os.path.join(a.outdir, "quickstart_inputs.png"), dpi=110)
    plt.close(fig)

    # final topography and its change
    fig = plt.figure(figsize=(11.0, 6.4), constrained_layout=True)
    ax = fig.add_subplot(2, 2, 1, projection="mollweide")
    pc = mollweide_map(ax, tg[-1], topo, cmap="terrain", vmin=-5000, vmax=2500)
    ax.set_title(f"final topography after {len(nit)} kyr (m); black: initial coast", fontsize=9)
    fig.colorbar(pc, ax=ax, shrink=0.8)
    ax = fig.add_subplot(2, 2, 2, projection="mollweide")
    dt_lim = nice_limit(tg[-1] - topo)      # 99th percentile of |change|, rounded up to a nice number
    pc = mollweide_map(ax, tg[-1] - topo, topo, cmap=CHANGE_CMAP, vmin=-dt_lim, vmax=dt_lim)
    ax.set_title("topography change (m)", fontsize=9)
    change_colorbar(fig, pc, dt_lim, label="topography change (m)", ax=ax, shrink=0.8)
    ax = fig.add_subplot(2, 2, 3)
    pc = south_polar_map(ax, tg[-1], topo, cmap="terrain", vmin=-5000, vmax=2500, latmax=-55.0)
    ax.set_title("final topography, south polar view (m)", fontsize=9)
    fig.colorbar(pc, ax=ax, shrink=0.8)
    ax = fig.add_subplot(2, 2, 4)
    pc = south_polar_map(ax, tg[-1] - topo, topo, cmap=CHANGE_CMAP, vmin=-dt_lim, vmax=dt_lim, latmax=-55.0)
    ax.set_title("topography change, south polar view (m)", fontsize=9)
    change_colorbar(fig, pc, dt_lim, label="topography change (m)", ax=ax, shrink=0.8)
    fig.savefig(os.path.join(a.outdir, "quickstart_topo.png"), dpi=110)
    plt.close(fig)

    # sea-level change maps: (a) RSL over the ocean, (b) its departure from the ocean mean
    # (c) topography change everywhere
    ocean_mean = float((w * rsl[-1])[ocean].sum() / w[ocean].sum())
    rsl_oc = np.where(ocean, rsl[-1], np.nan)
    fp = rsl_oc - ocean_mean
    rsl_lim, fp_lim = nice_limit(rsl_oc), nice_limit(fp)
    fig = plt.figure(figsize=(15.0, 3.4), constrained_layout=True)
    panels = [(rsl_oc, rsl_lim, f"(a) relative sea-level change over the ocean after {len(nit)} kyr (m)",
               "relative sea-level change (m)", True),
              (fp, fp_lim, f"(b) sea-level change minus the ocean mean ({ocean_mean:+.0f} m)",
               "departure from ocean mean (m)", True),
              (tg[-1] - topo, dt_lim, "(c) topography change (uplift > 0) (m)", "topography change (m)", False)]
    for k, (f, lim, title, label, sites) in enumerate(panels):
        ax = fig.add_subplot(1, 3, k + 1, projection="mollweide")
        pc = mollweide_map(ax, f, topo, ice=ice[0], cmap=CHANGE_CMAP, vmin=-lim, vmax=lim)
        if sites:
            mark_sites(ax, False)
        ax.set_title(title, fontsize=8)
        change_colorbar(fig, pc, lim, label=label, ax=ax, orientation="horizontal", shrink=0.7, pad=0.02)
    fig.savefig(os.path.join(a.outdir, "quickstart_map.png"), dpi=110)
    plt.close(fig)

    # time series
    fig, axs = plt.subplots(1, len(SITES), figsize=(12.0, 3.0), constrained_layout=True)
    for ax, (name, ij) in zip(axs, idx.items()):
        ax.plot(t, rsl[:, ij[0], ij[1]], "o-", color="C0", ms=3, label="marine check on")
        ax.plot(t, rsl_nm[:, ij[0], ij[1]], "--", color="C1", label="marine check off")
        ax.set_title(name, fontsize=8)
        ax.set_xlabel("time (kyr)")
        ax.axhline(0.0, color="0.7", lw=0.6)
        ax.axhline(ocean_mean, color="k", lw=0.8, ls=":", label=f"ocean mean at {t[-1]:.0f} kyr ({ocean_mean:+.0f} m)")
    axs[0].set_ylabel("relative sea-level change (m)")
    axs[0].legend(fontsize=8)
    fig.savefig(os.path.join(a.outdir, "quickstart_series.png"), dpi=110)
    plt.close(fig)
    print(f"figures written to {a.outdir} (colour limits: topography change +-{dt_lim:g} m, "
          f"ocean sea-level change +-{rsl_lim:g} m, departure from the ocean mean +-{fp_lim:g} m; ocean mean {ocean_mean:+.3f} m)")


if __name__ == "__main__":
    main()
