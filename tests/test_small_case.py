#!/usr/bin/env python3
"""Idealised axisymmetric ice age case on a reduced grid (128 x 256 Gauss grid, lmax 63).

Earth model: lovenum/prem_512.l96C.ump3.lm1 (first 63 degrees).  Topography: deep ocean (-4000 m)
with a polar continent whose surface slopes from +300 m at the north pole to -800 m at 20 degrees
colatitude; the ice is a parabolic cap of up to 1000 m and 14 degrees radius centred at the pole,
grown and melted over 8 steps of 1 kyr.  The outer part of the cap is afloat (removed by the
marine check), and the shoreline moves with the loading.  Rotational feedback on.

  (a) plain recursive solver (solver.SeaLevelSolver)
  (b) the same with the extrapolated inner start (solver_opt.OptSeaLevelSolver, inner="extrap")
  (c) the solver of (a) with the recursive modal state replaced by an explicit full-history
      evaluation of the viscous terms (class HistoryMemory below)

Checks: |a - c| < 1e-9 m for topography and radial displacement at every step; |a - b| within the
inner tolerance.  Also reports an independent fixed-point residual: for the converged step (run with
rotation off, so that the ocean-load update is a function of dS alone), the ocean-load update is
applied once more to the converged iterate.  The convergence criterion of the solver is the ratio
of Kendall et al. (2005) eq. 83, | sum|dS_new| - sum|dS_old| | / sum|dS_old| <= epsilon1, inherited
from the legacy code; it bounds a change of the coefficient sums, not the residual itself, so the
residual is reported for information.
"""
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "python"))
from slevolve.sht import SHT  # noqa: E402
from slevolve.solver import EPS0, ExactMemory, LoveNumbers, SeaLevelSolver  # noqa: E402
from slevolve.solver_opt import OptSeaLevelSolver  # noqa: E402

NLAT, NLON, LMAX, DT, EPS1 = 128, 256, 63, 1000.0, 1.0e-5
EARTH = os.path.join(HERE, "..", "lovenum", "prem_512.l96C.ump3.lm1")
H0 = [0.0, 300.0, 700.0, 1000.0, 1000.0, 800.0, 400.0, 150.0, 0.0]  # cap thickness, steps 0..8


class HistoryMemory(ExactMemory):
    """Same interface as ExactMemory, but the viscous terms are evaluated as the explicit sum over
    all committed load increments, sum_{j<n} dL_j sum_k a_k (1 - exp(-s_k (t_n - t_j)))."""

    def __init__(self, love, lmax, dt_years, tpw):
        super().__init__(love, lmax, dt_years, tpw)
        self.R = None
        self.s_kl = np.zeros((self.kmax, lmax + 1))
        for l in range(1, lmax + 1):
            nm = int(love.nmod[l])
            self.s_kl[:nm, l] = love.s[:nm, l]
        self.dt_years = dt_years
        self.hist = []  # committed increments (load, il, m, lambda)
        self.n = 0

    def advance_and_read(self):
        self.n += 1
        nl = self.lmax + 1
        visc = np.zeros((nl, nl), complex)
        rad = np.zeros((nl, nl), complex)
        sil, sm = np.zeros((3, 3)), np.zeros(3)
        stsl, strr = np.zeros(3, complex), np.zeros(3, complex)
        k2 = self.nmod2
        for j, (L, il, m, lam) in enumerate(self.hist, start=1):
            age = (self.n - j) * self.dt_years
            f = np.where(self.s_kl > 0, 1.0 - np.exp(-self.s_kl * age / 1000.0), 0.0)
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
        self.hist.append((current_load.copy(),
                          z((3, 3)) if dil is None else dil.copy(),
                          z(3) if dm is None else dm.copy(),
                          z(3, complex) if dlambda is None else dlambda.copy()))
        super().commit(current_load, dil, dm, dlambda)


class RecordingMemory(ExactMemory):
    """ExactMemory that keeps the viscous term of the last advance (for the residual check)."""

    def advance_and_read(self):
        out = super().advance_and_read()
        self.last_visc = out[0].copy()
        return out


def setup(sht):
    colat = np.degrees(np.arccos(sht.gaulats))[:, None] * np.ones((1, NLON))
    topo = np.where(colat <= 20.0, 300.0 - 1100.0 * colat / 20.0, -800.0)
    topo = np.where(colat > 25.0, -4000.0, np.where(colat > 20.0, -800.0 - 640.0 * (colat - 20.0), topo))
    shape = np.sqrt(np.clip(1.0 - (colat / 14.0) ** 2, 0.0, None))
    ice = [h * shape for h in H0]
    return topo, ice


def run(cls, sht, love, topo, ice, tpw=True, memory=None, **kw):
    slv = cls(sht, love, topo, ice[0], norder=LMAX, dt=DT, tpw=tpw, checkmarine=True,
              calcRG=True, epsilon1=EPS1, **kw)
    if memory is not None:
        slv.state = memory(love, LMAX, DT, tpw)
    out = []
    for k in range(1, len(ice)):
        r = slv.step(np.ascontiguousarray(ice[k]))
        out.append({"tgrid": r["tgrid"].copy(), "R": r["R"].copy(), "ocean": r["ocean"].copy(),
                    "numiter": r["numiter"]})
    return slv, out


def residual(slv, ice_prev_state):
    """Apply the ocean-load update once more to the converged dS of the last step (tpw off)."""
    p, sht, st = slv.p, slv.sht, slv.state
    dSprev, dI, dS = ice_prev_state, slv.deltaicestarlm, slv.load_inc_lm
    A = p.rhoi * dI + p.rhow * dSprev + p.rhow * dS
    dsl = A * slv._ttlekhl + st.last_visc * slv._ttl
    g = sht.synthesis(dsl, out=np.zeros((NLAT, NLON)))
    cstar = (slv.tinit - g < 0.0) * slv.beta  # ocean function of the converged iterate
    nl = LMAX + 1
    rO = sht.analysis(g * cstar, out=np.zeros((nl, nl), complex))
    cl = sht.analysis(cstar, out=np.zeros((nl, nl), complex))
    tO = sht.analysis(slv.tinit * (cstar - slv.cstar0), out=np.zeros((nl, nl), complex))
    conserv = ((1 / cl[0, 0]) * (-(p.rhoi / p.rhow) * dI[0, 0] - rO[0, 0] + tO[0, 0])).real
    new = -dSprev + rO + cl * conserv - tO
    dgrid = sht.synthesis(new - dS, out=np.zeros((NLAT, NLON)))
    xi = abs(np.abs(new).sum() - np.abs(dS).sum()) / np.abs(dS).sum()
    rel = np.abs(new - dS).sum() / np.abs(dS).sum()
    return float(np.abs(dgrid).max()), float(rel), float(xi)


def main():
    t_start = time.perf_counter()
    sht = SHT(NLAT, NLON, LMAX)
    love = LoveNumbers(EARTH, LMAX)
    topo, ice = setup(sht)

    _, a = run(SeaLevelSolver, sht, love, topo, ice)
    _, b = run(OptSeaLevelSolver, sht, love, topo, ice, inner="extrap")
    _, c = run(SeaLevelSolver, sht, love, topo, ice, memory=HistoryMemory)

    ocean0 = (topo < 0).astype(float)
    nshore = [int((x["ocean"] != ocean0).sum()) for x in a]
    print(f"test_small_case: grid {NLAT}x{NLON}, lmax {LMAX}, {len(a)} steps of {DT:.0f} yr; "
          f"cells whose ocean function changed per step: {nshore}")
    floating = int(((ice[3] > 0) & (topo < 0) & (920.0 * ice[3] <= 1000.0 * np.abs(topo))).sum())
    print(f"test_small_case: cells with floating ice at the maximum (removed by the marine check): {floating}")
    print(f"test_small_case: inner iterations (a) {[x['numiter'] for x in a]}")
    print(f"test_small_case: inner iterations (b) {[x['numiter'] for x in b]}")
    eac = max(max(np.abs(x["tgrid"] - y["tgrid"]).max(), np.abs(x["R"] - y["R"]).max()) for x, y in zip(a, c))
    eab = max(np.abs(x["tgrid"] - y["tgrid"]).max() for x, y in zip(a, b))
    dsl = max(np.abs(x["tgrid"] - topo).max() for x in a)
    print(f"test_small_case: (a) recursive vs (c) explicit history: max |diff| {eac:.3e} m "
          f"(topography and radial displacement, all steps)")
    print(f"test_small_case: (a) plain vs (b) extrapolated start: max |diff| {eab:.3e} m "
          f"(max |sea-level change| {dsl:.1f} m, ratio {eab / dsl:.1e}; epsilon1 = {EPS1:g})")

    # fixed-point residual of the converged iterate (rotation off), at every step
    slv = SeaLevelSolver(sht, love, topo, ice[0], norder=LMAX, dt=DT, tpw=False, checkmarine=True,
                         calcRG=True, epsilon1=EPS1)
    slv.state = RecordingMemory(love, LMAX, DT, False)
    worst = (0.0, 0.0, 0.0)
    for k in range(1, len(ice)):
        dSprev = slv.state.p_ocean.copy()
        slv.step(np.ascontiguousarray(ice[k]))
        r = residual(slv, dSprev)
        worst = tuple(max(u, v) for u, v in zip(worst, r))
    print(f"test_small_case: fixed-point residual (one more ocean-load update, rotation off): "
          f"max on grid {worst[0]:.3e} m, sum|dS_new - dS| / sum|dS| {worst[1]:.3e}, "
          f"eq. 83 ratio {worst[2]:.3e}")
    print(f"test_small_case: ({time.perf_counter() - t_start:.1f} s)")
    assert eac < 1e-9, eac
    assert eab / dsl < 1e-3, eab
    assert all(np.isfinite(worst))
    assert sum(nshore) > 0 and floating > 0
    print("test_small_case: PASS")


if __name__ == "__main__":
    main()
