#!/usr/bin/env python3
"""Second-order load placement (ramp) and the analytic rate read-out, on the idealised small case of
test_small_case.py (128 x 256 Gauss grid, lmax 63, rotation on, marine check on).

  (1) ramp factor: phi(0) = 1, phi + (1 - phi) = 1, series and closed form agree at the switch.
  (2) default unchanged: ramp=False and rates=True (output only) give exactly the default run.
  (3) recursion: the ramp solver with the recursive modal state against the same solver with an
      explicit full-history sum of the ramp kernel,
          sum_{j<n} dL_j sum_k a_k (1 - phi(s_k dt_j) exp(-s_k (t_n - t_j))),
      for uniform and non-uniform steps (|diff| < 1e-9 m, topography and radial displacement).
  (4) plain vs optimised solver with the ramp (numba kernels and the raw ducc0 backend when
      installed; extrapolated inner start): topography and rates agree.
  (5) analytic rates against a fine left difference on a piecewise-linear prescribed load
      (the step before t_n is split at t_n - eps): non-rotational read-out (ExactMemory.read_rate)
      and the rotational rates of SeaLevelSolver._rotation_rate (m-dot, sea-level term), each to
      1e-5 relative; inserting the extra point changes the ramp response by < 1e-9 relative.
"""
import math
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "python"))
from slevolve.sht import SHT  # noqa: E402
from slevolve.solver import ExactMemory, LoveNumbers, Planet, SeaLevelSolver, ramp_factor  # noqa: E402
from slevolve.solver_opt import OptSeaLevelSolver  # noqa: E402
from test_small_case import DT, EARTH, EPS1, LMAX, NLAT, NLON, setup  # noqa: E402

VSTEPS = [1000.0, 500.0, 2000.0, 250.0, 1000.0, 750.0, 1500.0, 500.0]


class HistoryMemoryRamp(ExactMemory):
    """Explicit full-history evaluation of the ramp kernel (step lengths through set_dt)."""

    def __init__(self, love, lmax, dt_years, tpw):
        super().__init__(love, lmax, dt_years, tpw)
        self.R = None
        self.t = 0.0
        self.hist = []  # (end time of the step, step length, load, il, m, lambda)

    def advance_and_read(self):
        self.t += self.dt_years
        nl = self.lmax + 1
        visc = np.zeros((nl, nl), complex)
        rad = np.zeros((nl, nl), complex)
        sil, sm = np.zeros((3, 3)), np.zeros(3)
        stsl, strr = np.zeros(3, complex), np.zeros(3, complex)
        k2 = self.nmod2
        for tj, hj, L, il, m, lam in self.hist:
            phi, _ = ramp_factor(self._s * hj / 1000.0)
            f = np.where(self._s > 0, 1.0 - phi * np.exp(-self._s * (self.t - tj) / 1000.0), 0.0)
            visc += (self.a_sl * f).sum(0)[:, None] * L
            rad += (self.a_rr * f).sum(0)[:, None] * L
            f2 = f[:k2, 2]
            sil += (self.a_il[:k2] * f2).sum() * il
            sm += (self.a_m[:k2] * f2).sum() * m
            stsl += (self.a_tsl[:k2] * f2).sum() * lam
            strr += (self.a_trr[:k2] * f2).sum() * lam
        return visc, rad, sil, sm, stsl, strr

    def _ramp_commit_load(self, current_load):
        pass  # no recursive state

    def commit(self, current_load, dil=None, dm=None, dlambda=None):
        z = np.zeros
        self.hist.append((self.t, self.dt_years, current_load.copy(),
                          z((3, 3)) if dil is None else dil.copy(), z(3) if dm is None else dm.copy(),
                          z(3, complex) if dlambda is None else dlambda.copy()))
        super().commit(current_load, dil, dm, dlambda)


def run(cls, sht, love, topo, ice, dts=None, memory=None, rates=False, **kw):
    dt0 = DT if dts is None else dts[0]
    slv = cls(sht, love, topo, ice[0], norder=LMAX, dt=dt0, tpw=True, checkmarine=True, calcRG=True,
              epsilon1=EPS1, **kw)
    if memory is not None:
        slv.state = memory(love, LMAX, dt0, True)
        slv.state.ramp = slv.ramp
    out = []
    for k in range(1, len(ice)):
        r = slv.step(np.ascontiguousarray(ice[k]), dt=None if dts is None else dts[k - 1], rates=rates)
        out.append({"tgrid": r["tgrid"].copy(), "R": r["R"].copy(), "numiter": r["numiter"],
                    "rates": r["rates"]})
    return out


def maxdiff(x, y, keys=("tgrid", "R")):
    return max(max(float(np.abs(a[k] - b[k]).max()) for k in keys) for a, b in zip(x, y))


def check_factor():
    x = np.array([0.0, 1e-9, 1e-5, 1e-4 * 0.999, 1e-4 * 1.001, 0.1, 1.0, 10.0])
    phi, omp = ramp_factor(x)
    ok = phi[0] == 1.0 and omp[0] == 0.0 and np.allclose(phi + omp, 1.0, atol=1e-15)
    ok &= np.allclose(phi[1:], -np.expm1(-x[1:]) / x[1:], rtol=1e-12) and abs(phi[3] - phi[4]) < 1e-7
    print(f"test_ramp: (1) ramp factor limits {'ok' if ok else 'FAILED'}")
    return ok


def check_rates_left_difference(love2):
    """(5) on a piecewise-linear prescribed load at lmax 2 (no sea-level equation)."""
    P = Planet()
    tc = 500.0 * np.arange(51)
    amp = np.where(tc < 5000, 1.0, np.clip((tc - 17000) / -12000.0, 0, 1)) * 3e4 + 50 * np.sin(tc / 900.0)

    def load(t):
        a = np.interp(t, tc, amp)
        L = np.zeros((len(t), 3, 3), complex)
        L[:, 1, 0], L[:, 1, 1] = 0.7 * a, 0.2 * a * (1 + 0.3j)
        L[:, 2, 0], L[:, 2, 1], L[:, 2, 2] = a, 0.3 * a * (1 - 0.5j), 0.1 * a * (0.4 + 1j)
        return L - L[0]

    def integrate(t, want_rate):
        st = ExactMemory(love2, 2, float(t[1] - t[0]), True)
        st.ramp = True
        kTE2, hTE2, kE2 = love2.kTE[2], love2.hTE[2], love2.ke[2]
        stub = object.__new__(SeaLevelSolver)
        stub.p, stub.state, stub.love = P, st, love2
        stub.il, stub.mm, stub.lam = np.zeros((3, 3)), np.zeros(3), np.zeros((3, 3), complex)
        r2w2 = P.radius**2 * P.omega**2
        L = load(t)
        e_sl = 1 + love2.ke[1:3] - love2.he[1:3]
        res, rot = [], []
        rates = None
        for j in range(1, len(t)):
            st.set_dt(float(t[j] - t[j - 1]))
            stub.dt = float(t[j] - t[j - 1])
            visc, _, sum_il, sum_m, tsl, _ = st.advance_and_read()
            cur = L[j] - L[j - 1]
            cur[0, :] = 0.0
            resp = e_sl[:, None] * L[j][1:3] + visc[1:3] + st.c_sl[1:3, None] * cur[1:3]
            oldil, oldm, lamp = st.p_il.copy(), st.p_m.copy(), st.p_lambda.copy()
            z20, z21 = L[j][2, 0], L[j][2, 1]
            il = stub.il
            il[2, 2] = (-8 * math.pi * P.radius**4 / (3 * math.sqrt(5.0))) * z20.real
            il[0, 2] = (8 * math.pi * P.radius**4 / math.sqrt(30.0)) * z21.real
            il[1, 2] = (-8 * math.pi * P.radius**4 / math.sqrt(30.0)) * (-1.0j * z21).real
            sil = sum_il + st.c_il * (il - oldil)
            al = 1 / (1 - kTE2 / P.kf)
            for i in (0, 1):
                stub.mm[i] = al * ((1 / (P.moiC - P.moiA)) * ((1 + kE2) * il[i, 2] + sil[i, 2])
                                   + (1 / P.kf) * (sum_m[i] - st.c_m * oldm[i])) / (1 - al * st.c_m / P.kf)
            stub.mm[2] = (1.0 / P.moiC) * ((1 + kE2) * il[2, 2] + sil[2, 2])
            m1, m2, m3 = stub.mm
            stub.lam[2, 0] = (r2w2 / (6.0 * math.sqrt(5.0))) * (m1**2 + m2**2 - 2.0 * m3**2 - 4.0 * m3)
            stub.lam[2, 1] = (r2w2 / math.sqrt(30.0)) * ((m1 - 1j * m2) * (1.0 + m3))
            stub.lam[2, 2] = (r2w2 / math.sqrt(5.0 * 24.0)) * (m2**2 - m1**2 + 2.0j * m1 * m2)
            tot = stub.lam[2, 0:3].copy()
            rsl = ((1 + kTE2 - hTE2) * tot + tsl + st.c_tsl * (tot - lamp)) / P.gacc
            res.append(resp)
            rot.append((stub.mm.copy(), rsl))
            if want_rate and j == len(t) - 1:
                rv, _ = st.read_rate(cur)
                lin = e_sl[:, None] * cur[1:3] / stub.dt + rv[1:3]
                rs, _, rr = SeaLevelSolver._rotation_rate(stub, oldil, oldm, lamp)
                rates = (lin, rr["mdot"], rs[2, 0:3])
            st.commit(cur, il - oldil, stub.mm - oldm, tot - lamp)
        return res, rot, rates

    res, rot, (lin, md, rsl) = integrate(tc, True)
    eps = 1e-2
    t2 = np.concatenate([tc[:-1], [tc[-1] - eps], tc[-1:]])
    res2, rot2, _ = integrate(t2, False)
    same = max(float(np.abs(res2[-1] - res[-1]).max() / np.abs(res[-1]).max()),
               float(np.abs(rot2[-1][1] - rot[-1][1]).max() / np.abs(rot[-1][1]).max()))
    fd_lin = (res2[-1] - res2[-2]) / eps
    fd_m = (rot2[-1][0] - rot2[-2][0]) / eps
    fd_sl = (rot2[-1][1] - rot2[-2][1]) / eps
    e_lin = float(np.abs(fd_lin - lin).max() / np.abs(lin).max())
    e_m = float(np.abs(fd_m[:2] - md[:2]).max() / np.abs(md[:2]).max())
    e_sl = float(np.abs(fd_sl - rsl).max() / np.abs(rsl).max())
    print(f"test_ramp: (5) analytic rate vs left difference (eps {eps} yr): non-rotational {e_lin:.2e}, "
          f"m-dot {e_m:.2e}, rotational sea level {e_sl:.2e} (relative); extra point changes the "
          f"response by {same:.1e}")
    return e_lin < 1e-5 and e_m < 1e-5 and e_sl < 1e-5 and same < 1e-9


def main():
    t_start = time.perf_counter()
    sht = SHT(NLAT, NLON, LMAX)
    love = LoveNumbers(EARTH, LMAX)
    topo, ice = setup(sht)
    ok = check_factor()

    a = run(SeaLevelSolver, sht, love, topo, ice)
    a0 = run(SeaLevelSolver, sht, love, topo, ice, ramp=False, rates=True)
    d_def = maxdiff(a, a0)
    print(f"test_ramp: (2) ramp=False with rates vs default run: max |diff| {d_def:.1e} m")
    ok &= d_def == 0.0

    r = run(SeaLevelSolver, sht, love, topo, ice, ramp=True, rates=True)
    h = run(SeaLevelSolver, sht, love, topo, ice, ramp=True, memory=HistoryMemoryRamp)
    rv = run(SeaLevelSolver, sht, love, topo, ice, dts=VSTEPS, ramp=True)
    hv = run(SeaLevelSolver, sht, love, topo, ice, dts=VSTEPS, ramp=True, memory=HistoryMemoryRamp)
    e_h, e_hv = maxdiff(r, h), maxdiff(rv, hv)
    shift = max(float(np.abs(x["tgrid"] - y["tgrid"]).max()) for x, y in zip(a, r))
    print(f"test_ramp: (3) ramp, recursive vs explicit history: max |diff| {e_h:.2e} m (uniform), "
          f"{e_hv:.2e} m (steps {VSTEPS}); ramp - step placement up to {shift:.1f} m")
    ok &= e_h < 1e-9 and e_hv < 1e-9 and shift > 1.0

    opt, sht2, tag = set(), sht, "numpy"
    try:
        import numba  # noqa: F401

        opt = {"numba"}
        tag = "numba"
        try:
            from slevolve.sht_ducc import make_sht

            sht2 = make_sht("ducc", NLAT, NLON, LMAX, workers=2, raw=True)
            opt = {"raw", "diag", "dsl", "numba", "sums"}
            tag = "numba + raw ducc0"
        except ImportError:
            pass
    except ImportError:
        pass
    b = run(OptSeaLevelSolver, sht2, love, topo, ice, ramp=True, rates=True, opt=opt)
    bx = run(OptSeaLevelSolver, sht2, love, topo, ice, ramp=True, opt=opt, inner="extrap")
    e_b = maxdiff(r, b)
    e_rate = max(float(np.abs(x["rates"][k] - y["rates"][k]).max()) for x, y in zip(r, b)
                 for k in ("dSLdt", "dRdt"))
    e_md = max(float(np.abs(x["rates"]["mdot"] - y["rates"]["mdot"]).max() / np.abs(x["rates"]["mdot"]).max())
               for x, y in zip(r, b))
    e_x = max(float(np.abs(x["tgrid"] - y["tgrid"]).max()) for x, y in zip(r, bx))
    dsl = max(float(np.abs(x["tgrid"] - topo).max()) for x in r)
    print(f"test_ramp: (4) ramp, plain vs optimised ({tag}): topography/R {e_b:.2e} m, rates "
          f"{e_rate:.2e} m/yr, m-dot {e_md:.1e} relative; extrapolated start {e_x:.2e} m "
          f"(max |sea-level change| {dsl:.0f} m, epsilon1 {EPS1:g})")
    ok &= e_b < 1e-9 and e_rate < 1e-12 and e_md < 1e-9 and e_x / dsl < 1e-3

    ok &= check_rates_left_difference(LoveNumbers(EARTH, 2))
    print(f"test_ramp: ({time.perf_counter() - t_start:.1f} s)")
    print("test_ramp: " + ("PASS" if ok else "FAIL"))
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
