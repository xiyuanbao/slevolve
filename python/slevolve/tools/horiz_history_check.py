#!/usr/bin/env python3
"""
horiz_history_check.py -- verification only: checks the implementation of the
recursive horizontal read-out against the explicit full-history sum over past load increments (Kendall et al. 2005 form):

    sum_k (l_k/s_k) R_k,lm(t_n)  ==  sum_{j<n} dL_lm(j) sum_k (l_k/s_k) (1 - exp(-s_k (t_n - t_j)))

(and the same for the degree-2 rotational term with the tidal l residues).
Runs the Python solver (Case 1 layout) for --steps steps, records the committed
increments, and compares at every step.  Short: 20 steps at L = 511 take seconds.
"""
import argparse
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))

from ..solver import LoveNumbers, Planet  # noqa: E402
from ..solver_opt import OptSeaLevelSolver, parse_opts  # noqa: E402
from ..sht_ducc import make_sht  # noqa: E402
from .. import npyio  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description="Verification only: compares the recursive horizontal read-out against the explicit full-history sum.")
    ap.add_argument("--inputs", required=True, help="INPUT_FILES tree (Case 1)")
    ap.add_argument("--earth", required=True)
    ap.add_argument("--icepfx", default="AISload_gd_Run85_0_0_")
    ap.add_argument("--topo", default="others/topo0_archive")
    ap.add_argument("--steps", type=int, default=20)
    ap.add_argument("--lmax", type=int, default=511)
    ap.add_argument("--threads", type=int, default=8)
    a = ap.parse_args()

    nlat, nlon = 512, 1024
    ice = npyio.IceDir(os.path.join(a.inputs, "icemodel"), a.icepfx, nlat, nlon)
    topo = npyio.read_grid(os.path.join(a.inputs, a.topo), nlat, nlon)
    love = LoveNumbers(a.earth, a.lmax)
    opts = parse_opts("all")
    sht = make_sht("ducc", nlat, nlon, a.lmax, workers=a.threads, raw=True)
    p = Planet()  # Han et al. (2022) constants, as run_ais.py --rhoi 920 --radius 6.371e6
    p.rhoi, p.radius = 920.0, 6.371e6
    p.moiA = 0.3296145 * p.mass * p.radius**2
    p.moiC = 0.3307007 * p.mass * p.radius**2
    slv = OptSeaLevelSolver(sht, love, topo, np.ascontiguousarray(ice[0]), norder=a.lmax, dt=1,
                            tpw=True, checkmarine=True, calcRG=True, planet=p, opt=opts)
    slv.horizontal = True
    st = slv.state
    s_kl = love.s[: st.kmax, : a.lmax + 1]  # inverse relaxation times, 1/kyr (0 where no mode)
    dt = 1.0
    P_hist = [np.zeros_like(st.p_load)]
    L_hist = [np.zeros(3, dtype=np.complex128)]
    worst = worst_rot = 0.0
    for n in range(1, a.steps + 1):
        slv.step(np.ascontiguousarray(ice[n]))
        rec, rec_thl = st.read_horizontal()  # states of step n (commit does not touch R_k)
        # explicit history sum over increments committed at steps j < n
        brute = np.zeros_like(rec)
        brute_thl = np.zeros(3, dtype=np.complex128)
        for j in range(1, n):
            dL = P_hist[j] - P_hist[j - 1]
            dLam = L_hist[j] - L_hist[j - 1]
            age = (n - j) * dt
            w = (st.a_hl * (1.0 - np.exp(-s_kl * age / 1000.0))).sum(axis=0)
            brute += w[:, None] * dL
            wt = (st.a_thl * (1.0 - np.exp(-s_kl[:, 2] * age / 1000.0))).sum()
            brute_thl += wt * dLam
        scale = max(np.abs(brute).max(), 1e-300)
        e = np.abs(rec - brute).max() / scale if n > 1 else np.abs(rec).max()
        et = (np.abs(rec_thl - brute_thl).max() / max(np.abs(brute_thl).max(), 1e-300)) if n > 1 else 0.0
        worst, worst_rot = max(worst, e), max(worst_rot, et)
        print(f"step {n:3d}  max|rec-brute|/max|brute| = {e:.3e}   rotational {et:.3e}")
        P_hist.append(st.p_load.copy())
        L_hist.append(st.p_lambda.copy())
    print(f"WORST relative difference: load {worst:.3e}, rotational {worst_rot:.3e}")


if __name__ == "__main__":
    main()
