#!/usr/bin/env python3
"""Floating-ice convention on the idealised case of test_small_case.py.

beta = 1 where the grounded ice I* is zero (Kendall et al. 2005, eqs 43-44), flotation tested
against the current topography iterate.  Checks: the plain solver (solver.SeaLevelSolver) and the
optimised one (solver_opt.OptSeaLevelSolver, with the numba kernels when numba is installed) agree
to < 1e-9 m; beta is 0 exactly on the grounded ice; the ocean under the floating part of the cap
is water-loading ocean (beta = 1 under floating ice).
"""
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "python"))
import test_small_case as tsc  # noqa: E402
from slevolve.sht import SHT  # noqa: E402
from slevolve.solver import EPS0, LoveNumbers, SeaLevelSolver  # noqa: E402
from slevolve.solver_opt import OptSeaLevelSolver  # noqa: E402


def main():
    t_start = time.perf_counter()
    sht = SHT(tsc.NLAT, tsc.NLON, tsc.LMAX)
    love = LoveNumbers(tsc.EARTH, tsc.LMAX)
    topo, ice = tsc.setup(sht)

    # ---- plain vs optimised solver ------------------------------------------------
    sa, a = tsc.run(SeaLevelSolver, sht, love, topo, ice)
    opt = ()
    try:
        import numba  # noqa: F401
        opt = ("numba",)
    except ImportError:
        pass
    sb, b = tsc.run(OptSeaLevelSolver, sht, love, topo, ice, opt=opt)
    eab = max(np.abs(x["tgrid"] - y["tgrid"]).max() for x, y in zip(a, b))
    print(f"test_floating_ice: inner iterations {[x['numiter'] for x in a]} "
          f"(optimised {[x['numiter'] for x in b]}); grounded-ice re-guesses / ice analyses redone "
          f"{sa.gr_counts['reguess']} / {sa.gr_counts['ice_analyses']}")
    print(f"test_floating_ice: plain vs optimised{'+numba' if opt else ''} solver: max |diff| {eab:.3e} m")

    # step 3 (ice maximum): floating cells carry water load, beta = 0 exactly on grounded ice
    float_cells = int(((ice[3] > 0) & (topo < 0) & (920.0 * ice[3] <= 1000.0 * np.abs(topo))).sum())
    s3 = SeaLevelSolver(sht, love, topo, ice[0], norder=tsc.LMAX, dt=tsc.DT, checkmarine=True,
                        epsilon1=tsc.EPS1)
    for k in range(1, 4):
        s3.step(np.ascontiguousarray(ice[k]))
    under_float = int(((ice[3] > 0) & (s3.icestarxy == 0) & (s3.beta == 1)).sum())
    beta_ok = bool(np.array_equal(s3.beta == 0, s3.icestarxy >= EPS0))
    print(f"test_floating_ice: step 3 (ice maximum): cells with floating ice {float_cells}, of which "
          f"water-loading ocean (beta = 1) {under_float}; beta = 0 exactly on grounded ice: {beta_ok}")
    print(f"test_floating_ice: ({time.perf_counter() - t_start:.1f} s)")
    assert eab < 1e-9, eab
    assert beta_ok and under_float > 0
    print("test_floating_ice: PASS")


if __name__ == "__main__":
    main()
