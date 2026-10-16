#!/usr/bin/env python3
"""plot_fields.py -- maps of SLEvolve grid fields: global Mollweide and Antarctic polar stereographic.

    python python/tools/plot_fields.py FIELD.npy --topo TOPO0.npy --out map.png
    python python/tools/plot_fields.py OUT/tgrid550.npy --minus TOPO0.npy --topo TOPO0.npy \\
        --view south --ice ICE.npy --sym 400 --title "topography change (m)" --out d.png
    python python/tools/plot_fields.py OUT/tgrid8.npy --minus TOPO0.npy --topo TOPO0.npy --sym auto --out d.png

FIELD is an (N, 2N) grid in the SLEvolve layout (rows north to south on Gauss-Legendre latitudes,
columns at longitudes i*360/(2N); 512 x 1024 or 128 x 256, any N works).  Coastlines are drawn on
every map as the 0 m contour of the topography passed with --topo (the present-day or initial
topography of the run); with --ice the grounding/ice margin is drawn as the contour ice = 0 (dashed).

Views: --view global (Mollweide, matplotlib's built-in projection='mollweide', longitude wrapped to
-180..180, north up), --view south (south polar stereographic, x, y computed here from latitude and
longitude, 0 E at the top), --view both (side by side).  Needs numpy and matplotlib only; if cartopy is
importable it is not required and not used for the projection itself, so the maps look the same with
and without it.

Colour convention for changes (topography change, sea-level change), as in Fig. 3a of the SLEvolve
paper: the diverging colormap CHANGE_CMAP = 'PuOr_r' with symmetric limits -L..L fixed per figure
and arrows on the colour bar (extend='both') for values beyond them, instead of saturating silently.
L is a fixed value (Case 1: 400 m, as in the paper) or nice_limit(field): the 99th percentile of
|field| rounded up to a nice number (1, 1.5, 2, 2.5, 3, 4, 5, 6, 8 x 10^k); --sym auto does this.

As a module:  from plot_fields import grid_latlon, mollweide_map, south_polar_map, nice_limit,
              change_colorbar, CHANGE_CMAP
"""
import argparse

import numpy as np

CHANGE_CMAP = "PuOr_r"     # diverging map of the paper's topography-change panels (Fig. 3a)


def nice_limit(f, pct=99.0):
    """pct-th percentile of |f| (finite values) rounded up to 1, 1.5, 2, 2.5, 3, 4, 5, 6 or 8 x 10^k."""
    v = float(np.percentile(np.abs(f[np.isfinite(f)]), pct))
    if v <= 0.0:
        return 1.0
    e = 10.0 ** np.floor(np.log10(v))
    for m in (1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0):
        if m * e >= v * (1.0 - 1e-12):
            return float(m * e)
    return float(10.0 * e)


def change_colorbar(fig, pc, lim, label="m", **kw):
    """Colour bar of a symmetric change map: extend arrows both ends, ticks at -lim, -lim/2, 0, lim/2, lim."""
    cb = fig.colorbar(pc, extend="both", label=label, **kw)
    cb.set_ticks(np.linspace(-lim, lim, 5))
    cb.outline.set_linewidth(0.4)
    return cb


def grid_latlon(nlat):
    """(lat_deg north to south, lon_deg) of the nlat x 2*nlat Gauss-Legendre grid."""
    x, _ = np.polynomial.legendre.leggauss(nlat)
    return np.degrees(np.arcsin(x[::-1])), np.arange(2 * nlat) * 360.0 / (2 * nlat)


def _wrap(f, lon):
    """Reorder columns to longitudes -180..180 and close the seam (adds one column)."""
    lw = np.where(lon >= 180.0, lon - 360.0, lon)
    o = np.argsort(lw, kind="stable")
    lw, g = lw[o], f[:, o]
    return np.append(lw, lw[0] + 360.0), np.concatenate([g, g[:, :1]], axis=1)


def mollweide_map(ax, f, topo, ice=None, coast_lw=0.6, **kw):
    """Draw f on a matplotlib 'mollweide' axes; coast = topo 0 m contour; optional ice = 0 contour."""
    lat, lon = grid_latlon(f.shape[0])
    lw, ff = _wrap(f, lon)
    X, Y = np.meshgrid(np.radians(np.minimum(lw, 180.0)), np.radians(lat))
    kw.setdefault("rasterized", True)
    pc = ax.pcolormesh(X, Y, ff, shading="nearest", **kw)
    _, tt = _wrap(topo, lon)
    ax.contour(X, Y, tt, levels=[0.0], colors="k", linewidths=coast_lw)
    if ice is not None:
        _, ii = _wrap(ice, lon)
        ax.contour(X, Y, ii, levels=[1.0], colors="0.25", linewidths=coast_lw, linestyles="--")
    ax.set_xticklabels([])
    ax.tick_params(axis="y", labelsize=6)
    ax.grid(True, lw=0.3, color="0.6")
    return pc


def stereo_xy(lat, lon):
    """South polar stereographic coordinates (unit sphere; 0 E up, 90 E right)."""
    rho = 2.0 * np.tan(np.radians(90.0 + lat) / 2.0)
    return rho * np.sin(np.radians(lon)), rho * np.cos(np.radians(lon))


def south_polar_map(ax, f, topo, ice=None, latmax=-50.0, coast_lw=0.6, **kw):
    """Draw f south of latmax in south polar stereographic projection; coast = topo 0 m contour."""
    lat, lon = grid_latlon(f.shape[0])
    rows = lat <= latmax
    lo = np.append(lon, 360.0)
    LON, LAT = np.meshgrid(lo, lat[rows])
    X, Y = stereo_xy(LAT, LON)
    close = lambda g: np.concatenate([g[rows], g[rows][:, :1]], axis=1)  # noqa: E731
    # cell edges (curvilinear in x, y): midpoints in latitude (the pole at the southern end), and
    # longitudes half a cell either side
    le = np.concatenate([[90.0], 0.5 * (lat[1:] + lat[:-1]), [-90.0]])
    i0 = int(np.argmax(rows))
    dl = lon[1] - lon[0]
    LONe, LATe = np.meshgrid(np.append(lon - 0.5 * dl, 360.0 - 0.5 * dl), le[i0:])
    Xe, Ye = stereo_xy(LATe, LONe)
    kw.setdefault("rasterized", True)
    pc = ax.pcolormesh(Xe, Ye, f[rows], shading="flat", **kw)
    ax.contour(X, Y, close(topo), levels=[0.0], colors="k", linewidths=coast_lw)
    if ice is not None:
        ax.contour(X, Y, close(ice), levels=[1.0], colors="0.25", linewidths=coast_lw, linestyles="--")
    # graticule: latitude circles every 10 degrees, meridians every 30 degrees
    t = np.linspace(0.0, 360.0, 361)
    for la in np.arange(-80.0, latmax + 0.1, 10.0):
        x, y = stereo_xy(np.full_like(t, la), t)
        ax.plot(x, y, lw=0.3, color="0.6")
    for lm in range(0, 360, 30):
        x, y = stereo_xy(np.array([-90.0, latmax]), np.array([lm, lm]))
        ax.plot(x, y, lw=0.3, color="0.6")
    rmax = stereo_xy(np.array([latmax]), np.array([0.0]))[1][0]
    for lm, lab in ((0, "0"), (90, "90E"), (180, "180"), (270, "90W")):
        x, y = stereo_xy(np.array([latmax]), np.array([float(lm)]))
        ax.text(1.07 * x[0], 1.07 * y[0], lab, ha="center", va="center", fontsize=6)
    ax.set_xlim(-1.15 * rmax, 1.15 * rmax)
    ax.set_ylim(-1.15 * rmax, 1.15 * rmax)
    ax.set_aspect("equal")
    ax.axis("off")
    return pc


def main():
    ap = argparse.ArgumentParser(description="global Mollweide / south polar stereographic map of a grid field")
    ap.add_argument("field", help="(N, 2N) .npy grid")
    ap.add_argument("--topo", required=True, help="topography grid whose 0 m contour is drawn as the coast")
    ap.add_argument("--minus", help="subtract this grid from FIELD (e.g. tgrid0 for a topography change)")
    ap.add_argument("--ice", help="ice-thickness grid: draw the ice margin (contour ice = 0)")
    ap.add_argument("--view", choices=("global", "south", "both"), default="global")
    ap.add_argument("--latmax", type=float, default=-50.0, help="northern edge of the south polar view (deg)")
    ap.add_argument("--cmap", default=CHANGE_CMAP, help=f"colormap (default {CHANGE_CMAP})")
    ap.add_argument("--sym", help="symmetric colour range -SYM..SYM with extend arrows on the colour bar; "
                                  "'auto': 99th percentile of |FIELD| rounded up to a nice number")
    ap.add_argument("--vmin", type=float)
    ap.add_argument("--vmax", type=float)
    ap.add_argument("--title", default="")
    ap.add_argument("--label", default="m", help="colour-bar label")
    ap.add_argument("--out", required=True, help="output image (png, pdf, ...)")
    ap.add_argument("--dpi", type=int, default=130)
    a = ap.parse_args()
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    f = np.load(a.field).astype(float)
    if a.minus:
        f = f - np.load(a.minus)
    topo = np.load(a.topo)
    ice = np.load(a.ice) if a.ice else None
    kw = {"cmap": a.cmap}
    lim = None
    if a.sym is not None:
        lim = nice_limit(f) if a.sym == "auto" else float(a.sym)
        kw.update(vmin=-lim, vmax=lim)
        print(f"colour limits -{lim:g}..{lim:g}")
    if a.vmin is not None:
        kw["vmin"] = a.vmin
    if a.vmax is not None:
        kw["vmax"] = a.vmax
    views = ("global", "south") if a.view == "both" else (a.view,)
    fig = plt.figure(figsize=(4.6 * len(views) + (1.4 if "global" in views else 0.4), 3.2), constrained_layout=True)
    for k, v in enumerate(views):
        if v == "global":
            ax = fig.add_subplot(1, len(views), k + 1, projection="mollweide")
            pc = mollweide_map(ax, f, topo, ice, **kw)
        else:
            ax = fig.add_subplot(1, len(views), k + 1)
            pc = south_polar_map(ax, f, topo, ice, latmax=a.latmax, **kw)
        if lim is not None:
            change_colorbar(fig, pc, lim, label=a.label, ax=ax, shrink=0.8)
        else:
            fig.colorbar(pc, ax=ax, shrink=0.8, label=a.label)
    if a.title:
        fig.suptitle(a.title, fontsize=10)
    fig.savefig(a.out, dpi=a.dpi)
    print("wrote", a.out)


if __name__ == "__main__":
    main()
