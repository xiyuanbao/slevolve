#!/usr/bin/env python3
"""
init_topo_segmented.py -- partially relaxed (segmented) initial-topography guess T0 for outer pass 1.

Segments.  V_n = sum(I_n * area) for the ice files n = 0..N.  V is smoothed with a 5-point running
mean (edge-padded); the extrema of the smoothed series plus both ends delimit monotonic legs.  Legs
whose smoothed volume change is below THRESHOLD (default 0.01 = 1 %) of the smoothed range are
pruned iteratively, smallest first (an end point is never removed, the interior extremum next to it
is).  Segment j = files a_j..b_j; tau_j = sum |dV_n| age_n / sum |dV_n| over the steps n in
(a_j, b_j], age_n = t_N - t_n (= dt (N - n) for uniform steps; t_n from --times when given).

Guess (variant B, self-consistent sea-level equation with an intermediate Earth response):
    dSL = c + S,  S = -sum_j GmR_j(rho_i (L_j - L_{j+1})) - GmR_w(-rho_w dS)
  L_j   = boundary ice load I(a_j) * grounded(I(a_j), T0); first = I_0 g_0, last = I_N g_N(present),
  GmR_j = (G - R) operator with Love numbers relaxed for tau_j (elastic + sum_k r_k/s_k (1-exp(-s_k tau))),
  GmR_w = same with the |dV_j|-weighted mean (h, k) of the segments (water load),
  dS    = dSL * C_p - T0 * (C_p - C_0) (water load), c from exact water-volume conservation (bisection).
  Iterated until the area-weighted RMS change is < 1e-3 m (max 40 iterations).
One segment reduces to the single-tau variant B.

CLI:
  python3 init_topo_segmented.py --inputs INPUT_FILES --icepfx nadi_ --nfiles 51 --dt 500 \\
      --earth INPUT_FILES/earthmodel/prem_512.l96C.ump3.lm1 --lmax 511 --out T0_segmented.npy
  (--start T.npy : starting T0 of the fixed point, default T_pres; --c0 : starting uniform term.)
"""
import argparse
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
from .solver import LoveNumbers, Planet  # noqa: E402
from . import npyio  # noqa: E402


def grounded(I, T, rhoi, rhow):
    return (I > 0) & ((T > 0) | (rhoi * I > rhow * np.abs(T)))


def segments(V, dt, threshold=0.01):
    """Monotonic segments of the ice-volume series V (files 0..N); dt in kyr: a scalar (uniform
    steps) or the times t_0..t_N of the files in kyr (non-uniform steps)."""
    V = np.asarray(V, float); n = len(V); N = n - 1
    if np.ndim(dt) == 0:  # uniform steps: the original expression, unchanged
        age_of = lambda i: dt * (N - i)  # noqa: E731
    else:
        tk = np.asarray(dt, float)
        assert tk.shape == (n,)
        age_of = lambda i: float(tk[N] - tk[i]) if np.ndim(i) == 0 else tk[N] - tk[i]  # noqa: E731
    vp = np.pad(V, 2, mode="edge"); sm = np.convolve(vp, np.ones(5) / 5, "valid")
    d = np.sign(np.diff(sm)); ext = [0] + [i for i in range(1, n - 1) if d[i - 1] * d[i] < 0] + [n - 1]
    thr = threshold * (sm.max() - sm.min())
    while len(ext) > 2:
        legs = [abs(sm[ext[i + 1]] - sm[ext[i]]) for i in range(len(ext) - 1)]; i = int(np.argmin(legs))
        if legs[i] >= thr: break
        if i == 0: del ext[1]
        elif i == len(ext) - 2: del ext[-2]
        else: del ext[i:i + 2]
    segs = []
    for a, b in zip(ext[:-1], ext[1:]):
        st = np.arange(a + 1, b + 1); dv = np.abs(V[st] - V[st - 1]); age = age_of(st)
        segs.append(dict(a=int(a), b=int(b), t_a_ka=age_of(a), t_b_ka=age_of(b),
                         tau=float((dv * age).sum() / dv.sum()), dV=float(V[b] - V[a])))
    return segs


class _Case:
    """Grid, present state and conservation constraint (as in run_outer.py)."""

    def __init__(self, Tp, I0, Ip, rhoi=920.0, rhow=1000.0, radius=6.371e6):
        nl, no = Tp.shape
        self.Tp, self.I0, self.Ip, self.rhoi, self.rhow, self.a = Tp, I0, Ip, rhoi, rhow, radius
        x, gw = np.polynomial.legendre.leggauss(nl); gw = gw[::-1]
        self.area = np.broadcast_to((gw * 2 * np.pi / no * radius**2)[:, None], (nl, no)).copy()
        self.gp = grounded(Ip, Tp, rhoi, rhow)
        self.Cp = (Tp < 0) & ~self.gp
        self.glob = np.ones_like(Tp, bool)

    def rms(self, f, m):
        return float(np.sqrt((self.area[m] * f[m]**2).sum() / self.area[m].sum()))

    def conserv_resid(self, T0):
        g0 = grounded(self.I0, T0, self.rhoi, self.rhow); C0 = (T0 < 0) & ~g0
        V0 = (np.clip(-T0, 0, None) * C0 * self.area).sum(); Vp = (np.clip(-self.Tp, 0, None) * self.Cp * self.area).sum()
        dI = self.I0 * g0 - self.Ip * self.gp
        return V0 - Vp + self.rhoi / self.rhow * (dI * self.area).sum(), g0, C0, dI

    def solve_c(self, S, c, half=80.0):
        lo, hi = c - half, c + half
        for _ in range(60):
            mid = 0.5 * (lo + hi)
            f = self.conserv_resid(self.Tp + mid + S)[0]
            lo, hi = (mid, hi) if f > 0 else (lo, mid)
        return 0.5 * (lo + hi)


def segmented_guess(Tp, ice, earth, lmax, dt, threshold=0.01, start=None, c0=62.0, rhoi=920.0,
                    rhow=1000.0, radius=6.371e6, workers=4, c_half=80.0, tol=1e-3, maxit=40):
    """Return (T0, info).  Tp: present topography (nlat, nlon); ice: sequence of ice grids 0..N
    (present = last); earth: Love-number file; dt: kyr per step, or the file times t_0..t_N in kyr; start: starting T0 (default Tp)."""
    from .sht_ducc import make_sht
    nl, no = Tp.shape; N = len(ice) - 1
    cs = _Case(Tp, np.asarray(ice[0]), np.asarray(ice[N]), rhoi, rhow, radius)
    V = np.array([(np.asarray(ice[n]) * cs.area).sum() for n in range(N + 1)])
    segs = segments(V, dt, threshold)
    love = LoveNumbers(earth, lmax)
    p = Planet(); p.radius = radius
    l = np.arange(lmax + 1, dtype=float); ttl = np.zeros_like(l)
    ttl[1:] = 4 * np.pi * radius**3 / (p.mass * (2 * l[1:] + 1))
    sht = make_sht("ducc", nl, no, lmax, workers=workers)

    def love_tau(tau):
        h, k = love.he.copy(), love.ke.copy()
        for ll in range(1, lmax + 1):
            nm = love.nmod[ll]; s = love.s[:nm, ll]; f = 1.0 - np.exp(-s * tau)
            h[ll] += (love.r[:nm, ll] / s * f).sum(); k[ll] += (love.rprime[:nm, ll] / s * f).sum()
        return h, k

    def GmR_factory(h, k):
        fac = ttl * (1 + k - h)
        def GmR(load):
            r = sht.synthesis(sht.analysis(np.ascontiguousarray(load)).copy() * fac[:, None]).copy()
            return r - (r * cs.area).sum() / cs.area.sum()
        return GmR

    hk = [love_tau(sg["tau"]) for sg in segs]; GmR = [GmR_factory(h, k) for h, k in hk]
    w = np.abs([sg["dV"] for sg in segs]); w = w / w.sum()
    GmRw = GmR[0] if len(segs) == 1 else GmR_factory(sum(wi * h for wi, (h, k) in zip(w, hk)),
                                                      sum(wi * k for wi, (h, k) in zip(w, hk)))
    Ib = [np.asarray(ice[sg["a"]]) for sg in segs] + [cs.Ip]
    dSL = np.zeros_like(Tp) if start is None else np.asarray(start, float) - Tp
    c = c0
    for it in range(maxit):
        T0 = Tp + dSL; _, g0, C0, dI = cs.conserv_resid(T0)
        dS = dSL * cs.Cp - T0 * (cs.Cp.astype(float) - C0)
        L = [cs.I0 * g0] + [I * grounded(I, T0, rhoi, rhow) for I in Ib[1:-1]] + [cs.Ip * cs.gp]
        Sx = (-sum(G_(rhoi * (L[j] - L[j + 1])) for j, G_ in enumerate(GmR)) - GmRw(-rhow * dS)) if len(segs) > 1 \
            else -GmR[0](rhoi * dI - rhow * dS)
        cn = cs.solve_c(Sx, c, c_half); new = cn + Sx; ch = cs.rms(new - dSL, cs.glob); dSL = new; c = cn
        if ch < tol: break
    return Tp + dSL, dict(segments=segs, c=c, iters=it + 1, last_change_m=ch, threshold=threshold,
                          dt_kyr=dt if np.ndim(dt) == 0 else [float(x) for x in dt],
                          earth=os.path.abspath(earth), lmax=lmax)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--inputs", required=True, help="INPUT_FILES directory (icemodel/, others/)")
    ap.add_argument("--icepfx", required=True, help="ice file prefix, e.g. nadi_")
    ap.add_argument("--nfiles", type=int, required=True, help="number of ice files (0..nfiles-1; last = present)")
    ap.add_argument("--dt", type=float, default=None,
                    help="step length in YEARS (time spacing of the ice files; converted to kyr internally); "
                         "required unless --times is given")
    ap.add_argument("--times", default=None,
                    help="optional times file: one time t_k in years per line for files 0..nfiles-1, strictly "
                         "increasing; the segment ages use these times (--dt ignored)")
    ap.add_argument("--earth", required=True, help="Love-number file")
    ap.add_argument("--lmax", type=int, required=True, help="spectral truncation (511 for the benchmark cases)")
    ap.add_argument("--topo-pres", default=None, help="present topography (default <inputs>/others/topo_pres)")
    ap.add_argument("--nglv", type=int, default=512, help="grid latitudes (default 512)")
    ap.add_argument("--threshold", type=float, default=0.01, help="leg pruning threshold, fraction of range (default 0.01)")
    ap.add_argument("--start", default=None, help="starting T0 of the fixed point (default T_pres)")
    ap.add_argument("--c0", type=float, default=62.0, help="starting uniform term c (m)")
    ap.add_argument("--rhoi", type=float, default=920.0, help="ice density (default 920)")
    ap.add_argument("--rhow", type=float, default=1000.0, help="water density (default 1000)")
    ap.add_argument("--radius", type=float, default=6.371e6, help="Earth radius in m (default 6.371e6)")
    ap.add_argument("--workers", type=int, default=4, help="ducc threads")
    ap.add_argument("--out", required=True, help="output .npy (a .json with the segments is written next to it)")
    a = ap.parse_args()
    if a.times:
        tk = npyio.read_times(a.times)
        if tk.size != a.nfiles:
            raise SystemExit(f"--times {a.times}: {tk.size} times for {a.nfiles} ice files")
        dtk = tk / 1000.0
        print(f"[init_topo] times file {a.times}: t0={tk[0]:g} yr, dt {np.diff(tk).min():g}-"
              f"{np.diff(tk).max():g} yr (--dt ignored)", flush=True)
    elif a.dt is None:
        raise SystemExit("--dt or --times is required")
    else:
        dtk = a.dt / 1000.0
        print(f"[init_topo] step length dt = {a.dt:g} yr = {a.dt / 1000.0:g} kyr", flush=True)
    nl, no = a.nglv, 2 * a.nglv
    Tp = npyio.read_grid(a.topo_pres or os.path.join(a.inputs, "others", "topo_pres"), nl, no)
    ice = [npyio.find_grid(os.path.join(a.inputs, "icemodel"), f"{a.icepfx}{k}", nl, no) for k in range(a.nfiles)]
    start = npyio.read_grid(a.start, nl, no) if a.start else None
    T0, info = segmented_guess(Tp, ice, a.earth, a.lmax, dtk, a.threshold, start, a.c0, a.rhoi, a.rhow, a.radius, a.workers)
    np.save(a.out, T0)
    json.dump(info, open(os.path.splitext(a.out)[0] + ".json", "w"), indent=1)
    print("segments", [(s["a"], s["b"], round(s["tau"], 3)) for s in info["segments"]], "c %.4f iters %d" % (info["c"], info["iters"]))


if __name__ == "__main__":
    main()
