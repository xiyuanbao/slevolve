#!/usr/bin/env python3
"""Present-day J2-dot by two independent routes, both pinned to a closed-form value.

WHY.  In releases up to 723d95f the axial Liouville equation used +1/C, inherited from the rotation
block of the legacy code; Mitrovica & Wahr (2011, AREPS 39, eq. 15) have m3 = -Delta I33 / C.  With
+1/C the inertia route to J2-dot had the opposite sign to the gravity route, and the round-off
reproduction of the legacy code could not catch it, because the legacy code had the same sign.  A
downstream J2-dot read-out with a matching sign slip ((1 - beta) instead of (1 + beta) in the
removal of the centrifugal term) then made J2-dot look right.  A sign in one route compensated
elsewhere is caught only by pinning each route to a value computed outside the solver; this test
does that, and also requires the two routes to agree.

Routes (analytic rates, second-order placement --ramp, as the documentation recommends):
  route 1 (gravity):  J2-dot = -sqrt(5) (G20-dot - Lambda20-dot / g) / a, G = SL + R from the rate
                      grids dSLdt + dRdt (degree (2,0) taken here by an independent Gauss-Legendre
                      quadrature), Lambda20-dot from the rotational rates;
  route 2 (inertia):  J2-dot = -3 C m3-dot / (2 M a^2) (axial Liouville equation, MW2011 eq. 15;
                      for a degree-2 mass-conserving change the inertia trace is unchanged, so
                      Delta J2 = (3/2) Delta I33 / (M a^2)).

Set-up (test-only input; solver unchanged).  32 x 64 Gauss grid, lmax 10, rotation on, no marine
check, Earth model lovenum/prem_512.l96C.ump3.lm1, 12 steps of 1 kyr.  The land (+5000 m) is covered
everywhere by ice of 3000 m except two Gauss latitude rows nearest the nodes of P2 (mu = +-1/sqrt 3),
which are ice-free ocean (-1000 m).  The ice changes by h(t) f(mu), f = P2 + c P4 with c chosen so
that f = 0 on the ocean rows; h is piecewise linear in time (ice grows, melts, keeps melting slowly
in the last step).  The change is zonal, symmetric about the equator and of zero mean, so the sea
level is the same on both ocean rows and mass conservation makes the ocean load exactly zero: the
load is the degree-2 ice load rho_i h(t) P2 (the P4 part does not enter J2).

Known value (closed form, independent code).  The Love-number file is parsed here directly (degree-2
block: elastic k^E, inverse relaxation times s_k, residues r_k of the load and tidal k).  Unit
Heaviside load response k(t) = k^E + sum_k (r_k / s_k)(1 - exp(-s_k t)), so for a piecewise-linear
load with slopes h'_j on [t_j, t_j+1]
    J2_load-dot(t) = -(4 pi a^2 rho_i / (5 M)) [ (1 + k^E) h'(t-)
                     + sum_j h'_j sum_k (r_k / s_k) (exp(-s_k (t - t_j+1)) - exp(-s_k (t - t_j))) ],
    m3 = -(2/3)(M a^2 / C) J2_load  (MW2011 eq. 15),
    J2_rot-dot  = (2 a omega^2 / (3 g)) d/dt (k_T * m3)   (rotational deformation, tidal k_T;
                  the convolution of the closed-form m3-dot with the modes of k_T is integrated with
                  40-point Gauss-Legendre rules per step).
  Route 1 measures J2_load + J2_rot; route 2 measures J2_load: the axial equation of the solver
  (as MW2011 eq. 15) has no rotational deformation in Delta I33.

Tolerances.
  route 2 vs J2_load-dot: 1e-8 relative.  m3 is a linear functional of the piecewise-linear load,
      which the ramp placement and the analytic rate integrate exactly; the ocean load is zero, so the
      inner iteration does not enter; what remains is the accuracy of the (2,0) function of the
      legacy transform, 4-pi-normalised to a relative 1.7e-9 (documented in outputs.md), which
      enters through the analysis of the ice load (measured 1.7e-9).
  route 1 vs J2_load-dot + J2_rot-dot: 5e-5 relative.  Lambda is not piecewise linear in time (m3
      carries the viscous memory), so the ramp placement of the rotational increments and the rate
      read-out of the rotational sea level have an O(dt^2) error in a term that is itself ~0.5 % of
      J2-dot: measured 3.2e-3 of J2_rot (1.6e-5 of J2-dot) at 1-kyr steps, 1.4e-3 and 3.7e-4 at 250
      and 125 yr.
  route 1 vs route 2: |difference| / |J2-dot| <= 1 %.  The difference is J2_rot, the rotational
      deformation that the axial equation omits in Delta I33 (MW2011 eq. 15): measured 0.49 % here,
      ~0.3 % in production runs; the ratio of the totals is bounded by the fluid limit
      (4/9)(a omega^2 / g)(M a^2 / C) k_T^fluid = 0.43 %, a ratio of rates can exceed it slightly
      because the rotational deformation lags.  The difference must also equal the closed-form
      J2_rot-dot to 1e-2 relative (the O(dt^2) error above).
A sign error of m3 flips route 2 (relative error 2) and changes route 1 by 2 J2_rot (~1 %), so each
check fails on its own.

Option on (axial_rotation = "deg2", "full"; same set-up).  The axial equation then contains the rotational
part of I33 (tidal k_T acting on the solver's own Lambda_20) and, for "full", the degree-0 trace term:
    m3 + Q (k_T * m3) + tau_c m3 = m3_L,  Q = (4/9)(a omega^2 / g)(M a^2 / C),  m3_L = -(2/3)(M a^2 / C) J2_load,
with tau_c = 0 for "deg2".  Known value (closed form, independent code): this linear equation together with the
modal states of k_T and of the load k and the piecewise-linear ice height is a linear time-invariant system per
1-kyr segment, integrated exactly with a matrix exponential (scaling and squaring, NumPy only); the left
derivative at the end gives m3-dot.  The quadratic terms of Lambda(m), relative size |m3| ~ 1e-9, are left out.
  route 2 ("deg2"): -3 C m3-dot / (2 M a^2), now the full J2-dot.  Route 2 ("full"): the trace term changes
      tr(I) but not J2, so route 2 is -(3 / (2 M a^2)) (C m3-dot + T Lambda00-dot), T = (3/2) tau_c C / (a^2 omega^2);
      the closed form is then J2-dot = -(3 / (2 M a^2)) C (1 + tau_c) m3-dot.
  Tolerances.
  route 1 vs route 2: 1e-8 relative.  Both routes now measure the same total; what remains is the accuracy of
      the (2,0) function of the legacy transform (relative 1.7e-9, as for route 2 above); measured 1.7e-9.  With
      the option off this difference is the 0.485 % above.
  route 1 and route 2 vs the closed form: 5e-5 relative.  The rotational part is not piecewise linear in time,
      so the ramp placement has the O(dt^2) error of route 1 above, now in both routes: measured 1.3e-5 at 1-kyr
      steps, a partial cancellation (2.2e-5, 6.7e-6 and 1.75e-6 at 500, 250 and 125 yr).
  m3-dot must also be smaller in magnitude than with the option off (the feedback is stabilising), by less
  than the fluid-limit ratio Q k_T^fluid (+ tau_c) = 0.43 % (+ 0.05 %): measured 0.48 % and 0.53 % (a ratio of
  rates can exceed the ratio of totals slightly, as above), so the bound is checked at twice that ratio.
"""
import math
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "python"))
from slevolve.sht import SHT  # noqa: E402
from slevolve.solver import LoveNumbers, Planet, SeaLevelSolver, TAU_C_PREM  # noqa: E402

NLAT, NLON, LMAX, DT = 32, 64, 10, 1000.0
EARTH = os.path.join(HERE, "..", "lovenum", "prem_512.l96C.ump3.lm1")
H = [0.0, 250.0, 500.0, 750.0, 1000.0, 1000.0, 700.0, 400.0, 100.0, 0.0, 0.0, 0.0, -20.0]  # h(t_n), m
ICE0, LAND, SEA = 3000.0, 5000.0, -1000.0
# model constants (Planet() of the solver; asserted equal below)
A_E, M_E, G_E, OMEGA, RHOI = 6.371e6, 5.976e24, 9.80665, 7.292e-5, 920.0
C_E = 0.3307007 * M_E * A_E**2
TOL_R2, TOL_R1, TOL_12, TOL_ROT = 1e-8, 5e-5, 1e-2, 1e-2
TOL_AX12, TOL_AXK = 1e-8, 5e-5          # option on: route 1 vs route 2, each route vs the closed form


def p2(x):
    return 0.5 * (3 * x * x - 1)


def p4(x):
    return (35 * x**4 - 30 * x * x + 3) / 8


def read_degree2(path):
    """Independent reader: elastic k, s_k (1/kyr) and k residues (1/kyr) of the degree-2 load and
    tidal blocks (k stored multiplied by the degree)."""
    tok = open(path).read().split("\n", 1)[1].split()
    i = 0
    while True:
        l, n = int(tok[i]), int(tok[i + 1])
        i += 2
        s = np.array(tok[i:i + n], float); i += n
        el = np.array(tok[i:i + 5], float); i += 10          # elastic, fluid
        rk = np.array(tok[i + 2 * n:i + 3 * n], float); i += 3 * n
        if l == 1:
            continue
        tel = np.array(tok[i:i + 5], float)
        tfl = np.array(tok[i + 5:i + 10], float); i += 10
        trk = np.array(tok[i + 2 * n:i + 3 * n], float)
        return {"s": s, "kE": el[4] / 2, "rk": rk / 2, "kTE": tel[4] / 2, "rkT": trk / 2,
                "kTF": tfl[4] / 2}


def j2load_rate(t, tn, slope, kE, s, rk):
    """Closed-form d/dt J2_load at t (yr), left derivative; times converted to kyr for s."""
    fac = -4 * math.pi * A_E**2 * RHOI / (5 * M_E)
    j = np.searchsorted(tn, t, side="left") - 1          # step containing t (left limit)
    tot = (1 + kE) * slope[j]
    for jj in range(j + 1):
        te = min(tn[jj + 1], t)
        tot += slope[jj] * np.sum((rk / s) * (np.exp(-s * (t - te) / 1e3) - np.exp(-s * (t - tn[jj]) / 1e3)))
    return fac * tot


def known_values(love):
    tn = DT * np.arange(len(H))
    slope = np.diff(H) / DT
    s, t_end = love["s"], tn[-1]
    j2l = j2load_rate(t_end, tn, slope, love["kE"], s, love["rk"])
    # m3-dot(t') = -(2/3)(M a^2 / C) J2_load-dot(t'); d/dt (k_T * m3) = k_T^E m3-dot(t) +
    # sum_k (r^T_k / s_k) s_k int exp(-s_k (t - t')) m3-dot(t') dt'
    cm = -(2.0 / 3.0) * M_E * A_E**2 / C_E
    xg, wg = np.polynomial.legendre.leggauss(40)
    conv = 0.0
    for j in range(len(tn) - 1):
        tq = tn[j] + (xg + 1) * DT / 2
        md = np.array([cm * j2load_rate(x, tn, slope, love["kE"], s, love["rk"]) for x in tq])
        ker = (love["rkT"][None, :] * np.exp(-s[None, :] * (t_end - tq[:, None]) / 1e3)).sum(1) / 1e3
        conv += (wg * DT / 2) @ (ker * md)
    kt_rate = love["kTE"] * cm * j2l + conv
    j2r = 2 * A_E * OMEGA**2 / (3 * G_E) * kt_rate
    bound = (4.0 / 9.0) * (A_E * OMEGA**2 / G_E) * (M_E * A_E**2 / C_E) * love["kTF"]
    return j2l, j2r, bound


def expm(a):
    """Matrix exponential by scaling and squaring of a Taylor series (NumPy only)."""
    nrm = np.abs(a).sum(1).max()
    sq = max(0, int(math.ceil(math.log2(nrm / 0.25)))) if nrm > 0 else 0
    b = a / 2.0**sq
    e = np.eye(len(a))
    term = np.eye(len(a))
    for k in range(1, 30):
        term = term @ b / k
        e = e + term
    for _ in range(sq):
        e = e @ e
    return e


def known_selfconsistent(love, tau_c):
    """J2-dot (and m3-dot, per yr) at the end of the history for m3 + Q (k_T * m3) + tau_c m3 = m3_L (linear).
    State X = [y_k (k_T modes), z_k (load-k modes), h, h'] with y_k' = -s_k y_k + m3, z_k' = -s_k z_k + h,
    m3 = (cL ((1 + k^E) h + sum r_k z_k) - Q sum r^T_k y_k) / (1 + Q k^T_E + tau_c), time in kyr."""
    s, rk, rkT, kE, kTE = love["s"], love["rk"], love["rkT"], love["kE"], love["kTE"]
    K = len(s)
    Q = (4.0 / 9.0) * (A_E * OMEGA**2 / G_E) * (M_E * A_E**2 / C_E)
    cL = -(2.0 / 3.0) * (M_E * A_E**2 / C_E) * (-4 * math.pi * A_E**2 * RHOI / (5 * M_E))  # m3_L per unit load
    D = 1.0 + Q * kTE + tau_c
    n = 2 * K + 2
    iy, iz, ih, isl = slice(0, K), slice(K, 2 * K), 2 * K, 2 * K + 1
    m3row = np.zeros(n)
    m3row[ih] = cL * (1 + kE) / D
    m3row[iz] = cL * rk / D
    m3row[iy] = -Q * rkT / D
    M = np.zeros((n, n))
    M[iy, iy] = -np.diag(s)
    M[iy, :] += m3row[None, :]
    M[iz, iz] = -np.diag(s)
    M[iz, ih] = 1.0
    M[ih, isl] = 1.0
    tn = DT * np.arange(len(H)) / 1e3
    slope = np.diff(H) / (DT / 1e3)
    X = np.zeros(n)
    X[ih] = H[0]
    for j in range(len(tn) - 1):
        X[isl] = slope[j]
        X = expm(M * (tn[j + 1] - tn[j])) @ X
    m3dot = (m3row @ (M @ X)) / 1e3          # left derivative, per yr
    return -3 * C_E * (1.0 + tau_c) * m3dot / (2 * M_E * A_E**2), m3dot


def run_solver(sht, axial_rotation="off"):
    mu = sht.gaulats
    rows = [int(np.argmin(np.abs(mu - 1 / math.sqrt(3)))), int(np.argmin(np.abs(mu + 1 / math.sqrt(3))))]
    mb = mu[rows[0]]
    c = -p2(mb) / p4(mb)
    f = (p2(mu) + c * p4(mu))[:, None] * np.ones((1, NLON))
    sea = np.zeros((NLAT, NLON), bool)
    sea[rows] = True
    topo = np.where(sea, SEA, LAND)
    ice = [np.where(sea, 0.0, ICE0 + h * f) for h in H]
    assert min(x.min() for x in ice) >= 0.0 and np.abs(f[sea]).max() < 1e-15
    love = LoveNumbers(EARTH, LMAX)
    slv = SeaLevelSolver(sht, love, topo, ice[0], norder=LMAX, dt=DT, tpw=True, checkmarine=False,
                         calcRG=True, epsilon1=1e-5, ramp=True, axial_rotation=axial_rotation)
    for k in range(1, len(ice)):
        r = slv.step(np.ascontiguousarray(ice[k]), rates=True)
    return r, sea


def routes(sht, r, tau_c=0.0):
    """route 1 (geoid) and route 2 (inertia; with the trace correction when tau_c > 0) and m3-dot."""
    rt = r["rates"]
    # independent degree-(2,0) coefficient (4-pi normalisation) of G-dot = dSLdt + dRdt
    x, w = np.polynomial.legendre.leggauss(NLAT)
    order = np.argsort(sht.gaulats)
    assert np.allclose(x, sht.gaulats[order], atol=1e-14)
    wq = np.empty(NLAT)
    wq[order] = w
    gdot = (rt["dSLdt"] + rt["dRdt"]).mean(1)
    g20 = math.sqrt(5.0) / 2 * np.sum(wq * p2(sht.gaulats) * gdot)
    lam20 = rt["lambdadot"][2, 0].real
    route1 = -math.sqrt(5.0) * (g20 - lam20 / G_E) / A_E
    tdot3 = 1.5 * tau_c * C_E / (A_E**2 * OMEGA**2) * rt["lambdadot"][0, 0].real
    route2 = -3 * (C_E * rt["mdot"][2] + tdot3) / (2 * M_E * A_E**2)
    return route1, route2, rt["mdot"][2]


def main():
    t0 = time.perf_counter()
    P = Planet()
    assert (P.radius, P.mass, P.gacc, P.omega, P.rhoi, P.moiC) == (A_E, M_E, G_E, OMEGA, RHOI, C_E)
    sht = SHT(NLAT, NLON, LMAX)
    r, sea = run_solver(sht)
    route1, route2, mdot_off = routes(sht, r)

    love = read_degree2(EARTH)
    j2l, j2r, bound = known_values(love)
    known = j2l + j2r
    ocean_load = float(np.abs(r["dS_converged"]).max())
    e2 = abs(route2 - j2l) / abs(j2l)
    e1 = abs(route1 - known) / abs(known)
    d12 = abs(route1 - route2) / abs(known)
    erot = abs((route1 - route2) - j2r) / abs(j2r)
    print(f"test_j2dot_two_routes: known (closed form, k2 of {os.path.basename(EARTH)}): "
          f"J2-dot {known:.6e} /yr = load {j2l:.6e} + rotational deformation {j2r:.6e}")
    print(f"test_j2dot_two_routes: route 1 (geoid, Lambda removed)  {route1:.6e} /yr, rel. error {e1:.1e} "
          f"(tol {TOL_R1:g})")
    print(f"test_j2dot_two_routes: route 2 (-3 C m3-dot / 2 M a^2)  {route2:.6e} /yr, rel. error vs load "
          f"part {e2:.1e} (tol {TOL_R2:g})")
    print(f"test_j2dot_two_routes: route 1 - route 2 = {d12 * 100:.3f} % of J2-dot (tol {TOL_12 * 100:g} %, "
          f"fluid-limit ratio {bound * 100:.2f} %); vs closed-form rotational term {erot:.1e} (tol {TOL_ROT:g}); "
          f"max |ocean load| {ocean_load:.1e} m")
    ok = e2 < TOL_R2 and e1 < TOL_R1 and d12 < TOL_12 and erot < TOL_ROT and ocean_load < 1e-9

    # ---- option on: the rotational part of I33 (deg2) and the trace term (full)
    for mode, tau in (("deg2", 0.0), ("full", TAU_C_PREM)):
        r, _ = run_solver(sht, mode)
        a1, a2, md = routes(sht, r, tau)
        kn, kmd = known_selfconsistent(love, tau)
        f1, f2, g12 = abs(a1 - kn) / abs(kn), abs(a2 - kn) / abs(kn), abs(a1 - a2) / abs(kn)
        dm = md / mdot_off - 1.0
        print(f"test_j2dot_two_routes: [{mode}] known (closed form with the rotational I33"
              + (", trace term" if tau else "") + f") J2-dot {kn:.6e} /yr; route 1 {a1:.6e} rel. error {f1:.1e}, "
              f"route 2 {a2:.6e} rel. error {f2:.1e} (tol {TOL_AXK:g}); route 1 - route 2 = {g12:.1e} of J2-dot "
              f"(tol {TOL_AX12:g}); m3-dot {dm * 100:+.3f} % vs option off, rel. error vs closed form "
              f"{abs(kmd / md - 1.0):.1e}")
        ok &= f1 < TOL_AXK and f2 < TOL_AXK and g12 < TOL_AX12 and abs(kmd / md - 1.0) < TOL_AXK
        ok &= -2 * (bound + tau) < dm < 0.0
    print(f"test_j2dot_two_routes: ({time.perf_counter() - t0:.1f} s)")
    print("test_j2dot_two_routes: " + ("PASS" if ok else "FAIL"))
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
