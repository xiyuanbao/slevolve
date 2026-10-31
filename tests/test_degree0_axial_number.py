#!/usr/bin/env python3
"""Elastic axial number tau_c: tools/degree0_axial_number.py against closed forms, and the solvers' constant.

(1) Homogeneous self-gravitating compressible sphere (solid, and fluid with mu = 0), g = (4/3) pi G rho r: the
    degree-0 problem M (div u)' = -(16/3) pi G rho^2 U - rho dPhi/dr, Phi = (r/a)^2, M = lambda + 2 mu, has the
    solution U = c r + d j1(k r), k^2 = (16/3) pi G rho^2 / M, c = -3 / (8 pi G rho a^2), d from
    sigma(a) = M div u - 4 mu U / r = 0.  int rho U r^3 dr from the tool (two-node deck, step 20 km)
    must match it to 1e-9.
(2) The deck of the shipped Love numbers (lovenum/prem_continental.deck): tau_c at steps 20 km and 10 km must
    agree to 1e-10 (convergence in the radial step).
(3) tau_c at 10 km must equal the constant of the Python solver (solver.TAU_C_PREM) to 1e-9, and the Fortran
    default (fortran/src/support_mod.f90 axial_tau_c, fortran/build.sh AXIAL_TAU_C) must be the same number, so
    the constants cannot drift from the tool.
"""
import math
import os
import re
import sys
import tempfile
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..")
sys.path.insert(0, os.path.join(ROOT, "python"))
sys.path.insert(0, os.path.join(ROOT, "tools"))
import degree0_axial_number as d0  # noqa: E402
from slevolve.solver import TAU_C_PREM  # noqa: E402

TOL_SPHERE, TOL_CONV, TOL_CONST = 1e-9, 1e-10, 1e-9


def sphere_closed_form(a, rho, vp, vs, G=d0.G_NEWTON):
    """int_0^a rho U r^3 dr for the homogeneous sphere under Phi = (r/a)^2."""
    mu = rho * vs * vs
    M = rho * vp * vp
    k = math.sqrt(16.0 / 3.0 * math.pi * G * rho * rho / M)
    c = -3.0 / (8.0 * math.pi * G * rho * a * a)
    j0 = lambda x: math.sin(x) / x  # noqa: E731
    j1 = lambda x: math.sin(x) / x**2 - math.cos(x) / x  # noqa: E731
    sig_p = (3 * M - 4 * mu) * c
    sig_h = M * k * j0(k * a) - 4 * mu * j1(k * a) / a
    dd = -sig_p / sig_h
    x, w = np.polynomial.legendre.leggauss(64)
    r = 0.5 * a * (x + 1)
    U = c * r + dd * np.array([j1(k * ri) for ri in r])
    return float(0.5 * a * np.sum(w * rho * U * r**3))


def main():
    t0 = time.perf_counter()
    ok = True
    a = 6.371e6
    with tempfile.TemporaryDirectory() as tmp:
        for name, rho, vp, vs in (("solid", 5500.0, 9000.0, 5000.0), ("fluid", 5500.0, 9000.0, 0.0)):
            path = os.path.join(tmp, name + ".deck")
            g = 4.0 / 3.0 * math.pi * d0.G_NEWTON * rho * a
            np.savetxt(path, [[0.0, rho, vp, vs, 0.0], [a, rho, vp, vs, g]])
            layers, rs = d0.deck_model(path)
            got = d0.solve(layers, rs, 20e3)["I0"]
            ref = sphere_closed_form(a, rho, vp, vs)
            e = abs(got / ref - 1.0)
            print(f"test_degree0_axial_number: homogeneous {name} sphere: int rho U r^3 dr {got:.12e}, closed form "
                  f"{ref:.12e}, rel. error {e:.1e} (tol {TOL_SPHERE:g})")
            ok &= e < TOL_SPHERE
    t40 = d0.tau_c_of_deck(h=20e3)
    t20 = d0.tau_c_of_deck(h=10e3)
    ec = abs(t40 / t20 - 1.0)
    print(f"test_degree0_axial_number: {os.path.basename(d0.DEFAULT_DECK)}: tau_c {t40:.15e} (step 20 km), "
          f"{t20:.15e} (10 km), rel. difference {ec:.1e} (tol {TOL_CONV:g})")
    ok &= ec < TOL_CONV
    es = abs(t20 / TAU_C_PREM - 1.0)
    print(f"test_degree0_axial_number: solver.TAU_C_PREM {TAU_C_PREM:.15e}, rel. difference {es:.1e} "
          f"(tol {TOL_CONST:g})")
    ok &= es < TOL_CONST
    for f, pat in (("fortran/src/support_mod.f90", r"real :: axial_tau_c=([0-9.Ee+-]+)"),
                   ("fortran/build.sh", r'AXIAL_TAU_C="\$\{AXIAL_TAU_C:-([0-9.Ee+-]+)\}"')):
        m = re.search(pat, open(os.path.join(ROOT, f)).read())
        same = m is not None and float(m.group(1)) == TAU_C_PREM
        print(f"test_degree0_axial_number: {f} default {m.group(1) if m else 'not found'}: "
              + ("equal to solver.TAU_C_PREM" if same else "DIFFERS from solver.TAU_C_PREM"))
        ok &= same
    print(f"test_degree0_axial_number: ({time.perf_counter() - t0:.1f} s)")
    print("test_degree0_axial_number: " + ("PASS" if ok else "FAIL"))
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
