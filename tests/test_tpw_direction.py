#!/usr/bin/env python3
"""Direction of true polar wander from the rotation records, pinned to the physical sign.

WHY.  The spectral coefficients of the solver use the U basis, complex harmonics without the
Condon-Shortley phase.  The rotation block, inherited from the rotation block of the legacy code,
used the order-1 expressions of Milne & Mitrovica (1998, GJI 133, eqs. 33c and 43a,b), which are
written for harmonics with that phase (their eq. 33e).  In releases up to 2ca4705 the inertia
products il13, il23 and the centrifugal coefficient Lambda21 therefore both had the wrong sign.  The
two errors cancelled in every field, so the fields, J2-dot, the TPW speed and the round-off
reproduction of the legacy code (which shares the sign) could not reveal them; only the direction of
the polar-motion record (m1, m2 and their rates) was reversed by 180 degrees.  This test pins that
direction to quantities computed here, outside the solver, so it fails on such a sign error.

Physical sign (Mitrovica & Wahr 2011, AREPS 39, eq. 10): m1 + i m2 = (I13 + i I23) / (C - A) for the
elastic and initial response, with I_i3 = -int sigma x_i x_3 dA the perturbation of the inertia
products by the surface load sigma; (m1, m2) is the displacement of the rotation pole toward 0 E and
90 E.  Mass added at 45 N, 0 E gives I13 < 0, m1 < 0: the pole moves away from the load, toward 180 E.

Set-up (test-only input; solver unchanged).  32 x 64 Gauss grid, lmax 10, rotation on, no marine
check, second-order placement (ramp) and analytic rates, Earth model lovenum/prem_512.l96C.ump3.lm1.
Land (+5000 m) everywhere except a zonal ocean cap of 30 degrees radius around the south pole
(-4000 m; the shoreline does not move).  An ice cap of 10 degrees radius centred at 45 N and lon0
grows linearly: 1 m in a first step of 1 yr, then 100 m per kyr over 6 steps of 1 kyr.
  case A: lon0 = 0 E.   The rasterised cap is symmetric about the 0 E - 180 E meridian plane of the
          grid, so I23 = 0 and the expected direction of the pole is exactly 180 E.
  case B: lon0 = 120 E. Not a grid meridian, so the rasterised cap is slightly asymmetric; the
          expected direction is that of (I13, I23) from the quadrature here (300.009 E).

Checks (each case):
  (a) sign and closed form after the first step (1 yr; viscous part small):
      (m1, m2) = alpha (1 + k2E) (I13, I23) / (C - A), alpha = 1 / (1 - k2TE / kf), with I13, I23 by
      Gauss-Legendre quadrature of the load (ice plus ocean, ocean load = water depth change over
      the ocean grid cells) and k2E, k2TE parsed here from the Love-number file;
      relative tolerance 2e-3: the viscous response accumulated within 1 yr (measured 2.1e-4);
  (b) TPW direction from the rate record: write TPWrate with the release writer (slevolve.npyio,
      9 significant digits), read it back, atan2(m2-dot, m1-dot) in degrees east must equal the
      direction of (I13-dot, I23-dot) of the quadrature of the load rate over the last step, which
      points away from the added mass (case A: 180 E): tolerance 1e-3 degrees (the record has 9
      digits; in case B the ocean part of the load and its response can differ in phase from the
      ice: measured 1e-6 degree);
  (c) the same direction from the geoid: N21 = G21 - Lambda21 / g, with the (2,1) coefficient of the
      rate grids G-dot = dSLdt + dRdt taken here by quadrature (U basis, 4-pi normalised Y21 =
      sqrt(15/2) sin(th) cos(th) e^{i phi}); (m1-dot, m2-dot) must be parallel to
      -(Re N21-dot, -Im N21-dot) within 1e-3 degrees (geoid of a mass excess is high over it);
  (d) Lambda21 of the record against the quadrature of Lambda = (|omega|^2 r^2 - (omega.r)^2) / 2 on
      the grid for omega = Omega (m1, m2, 1 + m3) of the record (U basis): relative 1e-9 (the field is
      a polynomial of degree 2, which the quadrature integrates exactly: measured 1.5e-14).
A reversed sign of (m1, m2) fails (a), (b) and (c) by 180 degrees or a relative error of 2; a
reversed Lambda21 alone fails (d).
"""
import math
import os
import sys
import tempfile
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "python"))
from slevolve.npyio import save_rates  # noqa: E402
from slevolve.sht import SHT  # noqa: E402
from slevolve.solver import LoveNumbers, Planet, SeaLevelSolver  # noqa: E402

NLAT, NLON, LMAX = 32, 64, 10
EARTH = os.path.join(HERE, "..", "lovenum", "prem_512.l96C.ump3.lm1")
LAND, SEA = 5000.0, -4000.0
DTS = [1.0] + [1000.0] * 6                       # step lengths, yr
H = np.concatenate([[0.0, 1.0], 1.0 + 100.0 * np.arange(1, 7)])   # ice thickness on the cap, m
# model constants (Planet() of the solver; asserted equal below)
A_E, M_E, G_E, OMEGA, RHOI, RHOW = 6.371e6, 5.976e24, 9.80665, 7.292e-5, 920.0, 1000.0
MOI_A, MOI_C, KF = 0.3296145 * M_E * A_E**2, 0.3307007 * M_E * A_E**2, 0.9342 + 0.008
TOL_A, TOL_DIR, TOL_LAM = 2e-3, 1e-3, 1e-9


def read_k2_elastic(path):
    """Independent reader: elastic degree-2 load k and tidal k (stored multiplied by the degree)."""
    tok = open(path).read().split("\n", 1)[1].split()
    i = 0
    while True:
        l, n = int(tok[i]), int(tok[i + 1])
        i += 2 + n
        el = np.array(tok[i:i + 5], float); i += 10 + 3 * n
        if l == 1:
            continue
        tel = np.array(tok[i:i + 5], float)
        return el[4] / 2, tel[4] / 2


class Quad:
    """Gauss-Legendre x uniform-longitude quadrature on the solver grid (rows north -> south)."""

    def __init__(self, mu):
        x, w = np.polynomial.legendre.leggauss(NLAT)
        order = np.argsort(mu)
        assert np.allclose(x, mu[order], atol=1e-14)
        self.w = np.empty(NLAT)
        self.w[order] = w
        self.mu = mu
        self.st = np.sqrt(1 - mu**2)
        self.lon = 2 * math.pi * np.arange(NLON) / NLON
        # dA / a^2 = dmu dphi
        self.dA = (self.w[:, None] * (2 * math.pi / NLON)) * np.ones((1, NLON))

    def inertia13_23(self, sigma):
        """I_i3 = -int sigma x_i x_3 dA (kg m^2) for a surface density sigma (kg/m^2)."""
        xz = A_E**2 * (self.st * self.mu)[:, None] * np.cos(self.lon)[None, :]
        yz = A_E**2 * (self.st * self.mu)[:, None] * np.sin(self.lon)[None, :]
        return (-A_E**2 * np.sum(sigma * xz * self.dA), -A_E**2 * np.sum(sigma * yz * self.dA))

    def u21(self, f):
        """U-basis (2,1) coefficient (1/4pi) int f conj(Y21) dOmega, Y21 = sqrt(15/2) sin cos e^{i phi}."""
        y21c = math.sqrt(7.5) * (self.st * self.mu)[:, None] * np.exp(-1j * self.lon)[None, :]
        return np.sum(f * y21c * self.dA) / (4 * math.pi)

    def lambda21(self, m):
        """U-basis (2,1) coefficient of Lambda = (|w|^2 a^2 - (w.r)^2) / 2 at r = a, w = Omega (m1, m2, 1 + m3)."""
        w = OMEGA * np.array([m[0], m[1], 1 + m[2]])
        rx = A_E * self.st[:, None] * np.cos(self.lon)[None, :]
        ry = A_E * self.st[:, None] * np.sin(self.lon)[None, :]
        rz = A_E * self.mu[:, None] * np.ones((1, NLON))
        lam = 0.5 * (w @ w * A_E**2 - (w[0] * rx + w[1] * ry + w[2] * rz) ** 2)
        return self.u21(lam)


def cap(mu, lon0):
    lat = np.degrees(np.arcsin(mu))[:, None]
    lon = (360.0 * np.arange(NLON) / NLON)[None, :]
    c = (np.sin(np.radians(45.0)) * np.sin(np.radians(lat))
         + np.cos(np.radians(45.0)) * np.cos(np.radians(lat)) * np.cos(np.radians(lon - lon0)))
    return c > np.cos(np.radians(10.0))


def angle_deg(x, y):
    return math.degrees(math.atan2(y, x)) % 360.0


def dangle(a, b):
    return abs((a - b + 180.0) % 360.0 - 180.0)


def run_case(lon0, sht, love, q, kE, kTE):
    mu = sht.gaulats
    sea = np.repeat((np.degrees(np.arcsin(mu)) < -60.0)[:, None], NLON, 1)
    topo = np.where(sea, SEA, LAND)
    cp = cap(mu, lon0)
    assert not (cp & sea).any()
    ice = [np.where(cp, h, 0.0) for h in H]
    slv = SeaLevelSolver(sht, love, topo, ice[0], norder=LMAX, dt=DTS[0], tpw=True, checkmarine=False,
                         calcRG=True, epsilon1=1e-12, ramp=True)
    tgrid_prev = topo.copy()
    out = {}
    for k in range(1, len(ice)):
        r = slv.step(np.ascontiguousarray(ice[k]), dt=DTS[k - 1], rates=True)
        assert np.array_equal(r["ocean"] > 0.5, sea), "the shoreline moved"
        # load since step 0: ice + water over the (fixed) ocean; water depth change = -(tgrid - topo)
        sig = RHOI * ice[k] + RHOW * np.where(sea, -(r["tgrid"] - topo), 0.0)
        sig_prev = RHOI * ice[k - 1] + RHOW * np.where(sea, -(tgrid_prev - topo), 0.0)
        tgrid_prev = r["tgrid"].copy()
        if k == 1:
            i13, i23 = q.inertia13_23(sig)
            al = 1.0 / (1.0 - kTE / KF)
            pred = al * (1 + kE) * np.array([i13, i23]) / (MOI_C - MOI_A)
            mm = r["mm"]
            out["a"] = (mm[:2].copy(), pred, (i13, i23),
                        float(np.hypot(*(mm[:2] - pred)) / np.hypot(*pred)))
    # (b) the rate record, through the release writer
    rt = r["rates"]
    with tempfile.TemporaryDirectory() as td:
        save_rates(td, "X", rt, fields=())
        rec = open(os.path.join(td, "TPWrate" + "X.txt")).read().split("\n")
    mdot = np.array(rec[1].split(), float)
    # load rate over the last step (piecewise linear in time under the ramp placement)
    i13d, i23d = (np.array(q.inertia13_23(sig)) - np.array(q.inertia13_23(sig_prev))) / DTS[-1]
    dir_rec = angle_deg(mdot[0], mdot[1])
    dir_load = angle_deg(i13d, i23d)          # MW2011 eq. 10: m-dot parallel to (I13-dot, I23-dot)
    # (c) geoid rate, centrifugal term removed
    n21 = q.u21(rt["dSLdt"] + rt["dRdt"]) - rt["lambdadot"][2, 1] / G_E
    dir_geoid = angle_deg(-n21.real, n21.imag)
    # (d) Lambda21 of the record vs quadrature for the record's m
    lam_q = q.lambda21(r["mm"])
    lam_r = r["lambda"][2, 1]
    out.update(b=(dir_rec, dir_load), c=dir_geoid, d=(lam_r, lam_q, abs(lam_r - lam_q) / abs(lam_q)),
               mdot=mdot)
    return out


def main():
    t0 = time.perf_counter()
    P = Planet()
    assert (P.radius, P.mass, P.gacc, P.omega, P.rhoi, P.rhow, P.moiA, P.moiC, P.kf) == \
        (A_E, M_E, G_E, OMEGA, RHOI, RHOW, MOI_A, MOI_C, KF)
    sht = SHT(NLAT, NLON, LMAX)
    love = LoveNumbers(EARTH, LMAX)
    kE, kTE = read_k2_elastic(EARTH)
    q = Quad(sht.gaulats)
    ok = True
    for name, lon0, expect in (("A", 0.0, 180.0), ("B", 120.0, None)):
        o = run_case(lon0, sht, love, q, kE, kTE)
        m, pred, (i13, i23), ea = o["a"]
        dir_rec, dir_load = o["b"]
        exp = dir_load if expect is None else expect
        big = np.abs(pred) > 1e-6 * np.hypot(*pred)  # signs of the non-zero components (case A: m2 = 0)
        ok_a = ea < TOL_A and bool(np.all(np.sign(m[big]) == np.sign(pred[big])))
        if expect is not None:                       # case A: added mass at 0 E, I13 < 0, m1 < 0
            ok_a = ok_a and m[0] < 0 and i13 < 0
        ok_b = dangle(dir_rec, exp) < TOL_DIR and dangle(dir_load, exp) < TOL_DIR
        ok_c = dangle(o["c"], dir_rec) < TOL_DIR
        ok_d = o["d"][2] < TOL_LAM
        ok = ok and ok_a and ok_b and ok_c and ok_d
        print(f"test_tpw_direction: case {name} (cap at 45 N, {lon0:g} E): step 1 (1 yr) I13, I23 = "
              f"{i13:.4e}, {i23:.4e} kg m^2 (quadrature); m1, m2 = {m[0]:.6e}, {m[1]:.6e}, closed form "
              f"{pred[0]:.6e}, {pred[1]:.6e}, rel. diff {ea:.1e} (tol {TOL_A:g}) "
              + ("ok" if ok_a else "FAIL"))
        print(f"test_tpw_direction: case {name}: TPW direction atan2(m2-dot, m1-dot) of TPWrate "
              f"{dir_rec:.6f} E; (I13-dot, I23-dot) of the load {dir_load:.6f} E"
              + (f" (expected {expect:g} E)" if expect is not None else "")
              + f" (tol {TOL_DIR:g} deg) " + ("ok" if ok_b else "FAIL"))
        print(f"test_tpw_direction: case {name}: direction from the geoid rate -(Re N21-dot, -Im N21-dot) "
              f"{o['c']:.6f} E (tol {TOL_DIR:g} deg) " + ("ok" if ok_c else "FAIL"))
        print(f"test_tpw_direction: case {name}: Lambda21 record {o['d'][0]:.6e}, quadrature of Lambda(omega) "
              f"{o['d'][1]:.6e}, rel. diff {o['d'][2]:.1e} (tol {TOL_LAM:g}) " + ("ok" if ok_d else "FAIL"))
    print(f"test_tpw_direction: ({time.perf_counter() - t0:.1f} s)")
    print("test_tpw_direction: " + ("PASS" if ok else "FAIL"))
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
