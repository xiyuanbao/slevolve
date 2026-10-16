#!/usr/bin/env python3
"""Recursive modal-state update vs the explicit history sum (synthetic spectrum, non-uniform steps).

The release class solver.ExactMemory advances R_k <- d R_k + (1 - d) P with d = exp(-s_k dt / 1000)
(s_k in 1/kyr, dt in yr).  For non-uniform steps the same update holds with the decay of each step,
so the test sets the step length with ExactMemory.set_dt (the per-step path the solvers use with a
times file; checked to be exactly equal to exp(-s dt / 1000) built here) and compares the response read after each advance,
sum_k a_k R_k, with the explicit sum over the committed load increments,
sum_{j<n} dL_j sum_k a_k (1 - exp(-s_k (t_n - t_j) / 1000)).
"""
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "python"))
from slevolve.solver import ExactMemory, NPAM  # noqa: E402


class SyntheticLove:
    """Minimal stand-in for solver.LoveNumbers: random decay rates and residues."""

    def __init__(self, lmax, kmax, rng):
        nl = lmax + 1
        self.nmod = np.zeros(nl, dtype=np.int64)
        self.s = np.zeros((NPAM, nl))
        self.r = np.zeros((NPAM, nl))
        self.rprime = np.zeros((NPAM, nl))
        self.rT = np.zeros((NPAM, nl))
        self.rprimeT = np.zeros((NPAM, nl))
        for l in range(1, nl):
            nm = int(rng.integers(1, kmax + 1))
            self.nmod[l] = nm
            self.s[:nm, l] = 10.0 ** rng.uniform(-2, 2, nm)  # 1/kyr: relaxation times 10 yr - 100 kyr
            for a in (self.r, self.rprime, self.rT, self.rprimeT):
                a[:nm, l] = rng.normal(size=nm)
        self.nmod[2] = kmax  # all modes at degree 2 (rotational terms)
        self.s[:kmax, 2] = 10.0 ** rng.uniform(-2, 2, kmax)
        for a in (self.r, self.rprime, self.rT, self.rprimeT):
            a[:kmax, 2] = rng.normal(size=kmax)
        self.max_modes = int(self.nmod.max())


def main():
    t_start = time.perf_counter()
    rng = np.random.default_rng(12345)
    lmax, kmax, nsteps = 12, 8, 60
    love = SyntheticLove(lmax, kmax, rng)
    st = ExactMemory(love, lmax, 1.0, True)
    s = np.zeros((st.kmax, lmax + 1))
    for l in range(1, lmax + 1):
        s[: love.nmod[l], l] = love.s[: love.nmod[l], l]
    exists = np.zeros_like(s, dtype=bool)
    for l in range(1, lmax + 1):
        exists[: love.nmod[l], l] = True
    dts = rng.uniform(1.0, 2000.0, nsteps)  # non-uniform step lengths, yr
    t = np.concatenate(([0.0], np.cumsum(dts)))  # t[n] = end of step n
    nl = lmax + 1
    dL, dIL, dM, dLam = [], [], [], []
    worst = 0.0
    for n in range(1, nsteps + 1):
        d = np.where(exists, np.exp(-s * dts[n - 1] / 1000.0), 0.0)
        st.set_dt(dts[n - 1])
        assert np.array_equal(st.decay, d) and np.array_equal(st.omd, 1.0 - d)
        visc, rad, sil, sm, stsl, strr = st.advance_and_read()
        # explicit history sum over the increments committed at steps j < n
        ev = np.zeros((nl, nl), complex)
        er = np.zeros((nl, nl), complex)
        eil = np.zeros((3, 3))
        em = np.zeros(3)
        etsl = np.zeros(3, complex)
        etrr = np.zeros(3, complex)
        for j in range(1, n):
            age = t[n] - t[j]  # committed at the end of step j, read after the advance of step n
            f = np.where(exists, 1.0 - np.exp(-s * age / 1000.0), 0.0)
            ev += (st.a_sl * f).sum(0)[:, None] * dL[j - 1]
            er += (st.a_rr * f).sum(0)[:, None] * dL[j - 1]
            f2 = f[: st.nmod2, 2]
            eil += (st.a_il[: st.nmod2] * f2).sum() * dIL[j - 1]
            em += (st.a_m[: st.nmod2] * f2).sum() * dM[j - 1]
            etsl += (st.a_tsl[: st.nmod2] * f2).sum() * dLam[j - 1]
            etrr += (st.a_trr[: st.nmod2] * f2).sum() * dLam[j - 1]
        for rec, ex in ((visc, ev), (rad, er), (sil, eil), (sm, em), (stsl, etsl), (strr, etrr)):
            sc = max(np.abs(ex).max(), 1e-300)
            if n > 1:
                worst = max(worst, np.abs(rec - ex).max() / sc)
        # random step-wise load increment (lower triangle m <= l, l >= 1)
        L = rng.normal(size=(nl, nl)) + 1j * rng.normal(size=(nl, nl))
        L = np.tril(L)
        L[0, :] = 0.0
        il = np.zeros((3, 3))
        il[:, 2] = rng.normal(size=3)
        m = rng.normal(size=3)
        lam = rng.normal(size=3) + 1j * rng.normal(size=3)
        st.commit(L, il, m, lam)
        dL.append(L), dIL.append(il), dM.append(m), dLam.append(lam)
    wall = time.perf_counter() - t_start
    print(f"test_recursion: {nsteps} non-uniform steps (dt {dts.min():.1f}-{dts.max():.1f} yr), "
          f"lmax {lmax}, up to {kmax} modes: max relative error {worst:.3e}  ({wall:.2f} s)")
    assert worst < 1e-12, worst
    print("test_recursion: PASS")


if __name__ == "__main__":
    main()
