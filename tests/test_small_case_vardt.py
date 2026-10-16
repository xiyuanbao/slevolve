#!/usr/bin/env python3
"""Idealised small case of test_small_case.py run with NON-UNIFORM time steps from a times file.

Same grid, Earth model, topography and ice (8 steps, the same ice per step), but the steps are
1, 2, 5, 10, 20, 50, 100, 200 yr, given as a times file (t_k = 0, 1, 3, 8, ...) read with
npyio.read_times, as run_ais.py --times does; step k is advanced with dt = t_k - t_(k-1).

  (a) recursive solver (solver.SeaLevelSolver.step(ice, dt=dt_k))
  (c) the same with the modal state replaced by an explicit full-history sum,
      sum_{j<n} dL_j sum_k a_k (1 - exp(-s_k (t_n - t_j) / 1000))       (class HistoryMemoryVar)
  (b) solver_opt.OptSeaLevelSolver with the extrapolated inner start; the step lengths passed to
      the guess routine are recorded and must be (dt_j, dt_(j-1), dt_(j-2)) from step 3 on.
  (d) regression: a uniform times file (t_k = k * 1000 yr) must reproduce the run without a times
      file (constant dt = 1000 yr) exactly, for (a) and (b).

Checks: |a - c| < 1e-9 m (topography and radial displacement, every step); (d) max diff 0.0.
"""
import os
import sys
import tempfile
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "python"))
from slevolve import npyio  # noqa: E402
from slevolve import solver_opt  # noqa: E402
from slevolve.sht import SHT  # noqa: E402
from slevolve.solver import ExactMemory, LoveNumbers, SeaLevelSolver  # noqa: E402
from slevolve.solver_opt import OptSeaLevelSolver  # noqa: E402
from test_small_case import EARTH, EPS1, LMAX, NLAT, setup  # noqa: E402

STEPS = [1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0, 200.0]


class HistoryMemoryVar(ExactMemory):
    """Explicit full-history evaluation with the step lengths set through set_dt."""

    def __init__(self, love, lmax, dt_years, tpw):
        super().__init__(love, lmax, dt_years, tpw)
        self.R = None
        self.s_kl = np.zeros((self.kmax, lmax + 1))
        for l in range(1, lmax + 1):
            nm = int(love.nmod[l])
            self.s_kl[:nm, l] = love.s[:nm, l]
        self.t = 0.0
        self.hist = []  # (time of commit, load, il, m, lambda)

    def advance_and_read(self):
        self.t += self.dt_years  # end of the current step
        nl = self.lmax + 1
        visc = np.zeros((nl, nl), complex)
        rad = np.zeros((nl, nl), complex)
        sil, sm = np.zeros((3, 3)), np.zeros(3)
        stsl, strr = np.zeros(3, complex), np.zeros(3, complex)
        k2 = self.nmod2
        for tj, L, il, m, lam in self.hist:
            f = np.where(self.s_kl > 0, 1.0 - np.exp(-self.s_kl * (self.t - tj) / 1000.0), 0.0)
            visc += (self.a_sl * f).sum(0)[:, None] * L
            rad += (self.a_rr * f).sum(0)[:, None] * L
            f2 = f[:k2, 2]
            sil += (self.a_il[:k2] * f2).sum() * il
            sm += (self.a_m[:k2] * f2).sum() * m
            stsl += (self.a_tsl[:k2] * f2).sum() * lam
            strr += (self.a_trr[:k2] * f2).sum() * lam
        return visc, rad, sil, sm, stsl, strr

    def commit(self, current_load, dil=None, dm=None, dlambda=None):
        z = np.zeros
        self.hist.append((self.t, current_load.copy(),
                          z((3, 3)) if dil is None else dil.copy(),
                          z(3) if dm is None else dm.copy(),
                          z(3, complex) if dlambda is None else dlambda.copy()))
        super().commit(current_load, dil, dm, dlambda)


def run(cls, sht, love, topo, ice, dts, memory=None, **kw):
    """dts: per-step lengths (years) or None (constructor dt = 1000 yr, no per-step dt)."""
    dt0 = 1000.0 if dts is None else float(dts[0])
    slv = cls(sht, love, topo, ice[0], norder=LMAX, dt=dt0, tpw=True, checkmarine=True,
              calcRG=True, epsilon1=EPS1, **kw)
    if memory is not None:
        slv.state = memory(love, LMAX, dt0, True)
    out = []
    for k in range(1, len(ice)):
        r = slv.step(np.ascontiguousarray(ice[k]), dt=None if dts is None else float(dts[k - 1]))
        out.append({"tgrid": r["tgrid"].copy(), "R": r["R"].copy(), "numiter": r["numiter"]})
    return out


def maxdiff(x, y):
    return max(max(np.abs(a["tgrid"] - b["tgrid"]).max(), np.abs(a["R"] - b["R"]).max()) for a, b in zip(x, y))


def main():
    t_start = time.perf_counter()
    sht = SHT(NLAT, 2 * NLAT, LMAX)
    love = LoveNumbers(EARTH, LMAX)
    topo, ice = setup(sht)
    with tempfile.TemporaryDirectory() as d:
        p = os.path.join(d, "times")
        np.savetxt(p, np.concatenate(([0.0], np.cumsum(STEPS))), fmt="%.17g")
        dts = np.diff(npyio.read_times(p, len(STEPS)))
        pu = os.path.join(d, "times_uniform")
        np.savetxt(pu, 1000.0 * np.arange(len(STEPS) + 1), fmt="%.17g")
        dtu = np.diff(npyio.read_times(pu, len(STEPS)))
    assert list(dts) == STEPS

    a = run(SeaLevelSolver, sht, love, topo, ice, dts)
    c = run(SeaLevelSolver, sht, love, topo, ice, dts, memory=HistoryMemoryVar)
    calls = []
    orig = solver_opt.extrap_guess

    def rec(h1, h2, dt0, dt1, dt2, out=None):
        calls.append((dt0, dt1, dt2))
        return orig(h1, h2, dt0, dt1, dt2, out=out)

    solver_opt.extrap_guess = rec
    try:
        b = run(OptSeaLevelSolver, sht, love, topo, ice, dts, inner="extrap")
    finally:
        solver_opt.extrap_guess = orig
    expect = [(STEPS[j], STEPS[j - 1], STEPS[j - 2]) for j in range(2, len(STEPS))]

    eac = maxdiff(a, c)
    eab = max(np.abs(x["tgrid"] - y["tgrid"]).max() for x, y in zip(a, b))
    dsl = max(np.abs(x["tgrid"] - topo).max() for x in a)
    print(f"test_small_case_vardt: grid {NLAT}x{2 * NLAT}, lmax {LMAX}, steps {STEPS} yr")
    print(f"test_small_case_vardt: inner iterations (a) {[x['numiter'] for x in a]}, (b) {[x['numiter'] for x in b]}")
    print(f"test_small_case_vardt: (a) recursive vs (c) explicit history: max |diff| {eac:.3e} m "
          f"(topography and radial displacement, all steps; max |sea-level change| {dsl:.1f} m)")
    print(f"test_small_case_vardt: (a) vs (b) extrapolated start: max |diff| {eab:.3e} m")
    print(f"test_small_case_vardt: step lengths passed to the extrapolation (steps 3..8): {calls}")

    a0 = run(SeaLevelSolver, sht, love, topo, ice, None)
    a1 = run(SeaLevelSolver, sht, love, topo, ice, dtu)
    b0 = run(OptSeaLevelSolver, sht, love, topo, ice, None, inner="extrap")
    b1 = run(OptSeaLevelSolver, sht, love, topo, ice, dtu, inner="extrap")
    ra, rb = maxdiff(a0, a1), maxdiff(b0, b1)
    print(f"test_small_case_vardt: (d) uniform times file vs no times file: max |diff| (a) {float(ra):.1f} m, (b) {float(rb):.1f} m")
    print(f"test_small_case_vardt: ({time.perf_counter() - t_start:.1f} s)")
    assert eac < 1e-9, eac
    assert eab / dsl < 1e-3, eab
    assert calls == expect, (calls, expect)
    assert ra == 0.0 and rb == 0.0, (ra, rb)
    print("test_small_case_vardt: PASS")


if __name__ == "__main__":
    main()
