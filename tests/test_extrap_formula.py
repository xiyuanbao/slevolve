#!/usr/bin/env python3
"""Extrapolated inner start for uniform and non-uniform steps (paper Sect. 2.5).

    dS0_j = dt_j [ Sdot_{j-1} + (Sdot_{j-1} - Sdot_{j-2}) (dt_j + dt_{j-1}) / (dt_{j-1} + dt_{j-2}) ],
    Sdot_i = dS_i / dt_i,   = 2 dS_{j-1} - dS_{j-2} for equal steps.

A known sequence of increments dS_i (random spectral arrays) and step lengths dt_i is fed through
(1) the Python guess routine solver_opt.extrap_guess, used by OptSeaLevelSolver.step, and
(2) the Fortran routine modal_state_mod::extrap_coeffs, used by sl_model_mod::sl_solver (compiled
    from fortran/src/modal_state_mod.f90 with the release flags by $FC or gfortran; skipped with a
    message when no working compiler is found),
and compared with the formula evaluated independently here, for a uniform and a non-uniform
schedule.  For the uniform schedule the result must be exactly 2 dS_{j-1} - dS_{j-2}.
"""
import os
import shutil
import subprocess
import sys
import tempfile

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "python"))
from slevolve.solver_opt import extrap_guess  # noqa: E402

SCHEDULES = {
    "uniform": [1.0] * 9,
    "non-uniform": [1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0, 200.0, 37.5],
}


def formula(h1, h2, dtj, dtm1, dtm2):
    s1, s2 = h1 / dtm1, h2 / dtm2
    return dtj * (s1 + (s1 - s2) * (dtj + dtm1) / (dtm1 + dtm2))


FDRIVER = """program t
  use modal_state_mod, only: dp, extrap_coeffs
  implicit none
  real(dp) :: a, b, c, c1, c2
  integer :: ios
  do
    read(*,*,iostat=ios) a, b, c
    if (ios /= 0) exit
    call extrap_coeffs(a, b, c, c1, c2)
    write(*,'(2ES26.17E3)') c1, c2
  enddo
end program t
"""


def fortran_coeffs(triples):
    fc = shutil.which(os.environ.get("FC", "gfortran"))
    if fc is None:
        return None
    try:
        return _fortran_coeffs(fc, triples)
    except (subprocess.CalledProcessError, OSError) as e:
        print(f"test_extrap_formula: Fortran compile/run with {fc} failed ({e.__class__.__name__}); "
              "activate the Fortran toolchain or set FC to test the Fortran routine")
        return None


def _fortran_coeffs(fc, triples):
    src = os.path.join(HERE, "..", "fortran", "src", "modal_state_mod.f90")
    with tempfile.TemporaryDirectory() as d:
        drv = os.path.join(d, "t.f90")
        open(drv, "w").write(FDRIVER)
        flags = ["-O3", "-ffree-line-length-none", "-fdefault-real-8", "-fdefault-double-8", "-cpp"]
        subprocess.run([fc, *flags, "-c", src, "-o", os.path.join(d, "m.o"), "-J", d], check=True, cwd=d, capture_output=True)
        subprocess.run([fc, *flags, drv, os.path.join(d, "m.o"), "-I", d, "-o", os.path.join(d, "t")],
                       check=True, cwd=d, capture_output=True)
        inp = "".join(f"{a!r} {b!r} {c!r}\n" for a, b, c in triples)
        out = subprocess.run([os.path.join(d, "t")], input=inp, capture_output=True, text=True, check=True).stdout
    return [tuple(float(x) for x in line.split()) for line in out.strip().splitlines()]


def main():
    rng = np.random.default_rng(7)
    worst_py, worst_f = 0.0, 0.0
    fortran_ran = True
    for name, dts in SCHEDULES.items():
        n = len(dts)
        dS = [rng.normal(size=(6, 6)) + 1j * rng.normal(size=(6, 6)) for _ in range(n)]
        # steps j = 3..n (1-based) use the guess; steps 1, 2 use the eustatic start
        triples = [(dts[j - 1], dts[j - 2], dts[j - 3]) for j in range(3, n + 1)]
        fco = fortran_coeffs(triples)
        if fco is None:
            fortran_ran = False
        for idx, j in enumerate(range(3, n + 1)):
            h1, h2 = dS[j - 2], dS[j - 3]  # dS_{j-1}, dS_{j-2}
            dtj, dtm1, dtm2 = triples[idx]
            ref = formula(h1, h2, dtj, dtm1, dtm2)
            got = extrap_guess(h1, h2, dtj, dtm1, dtm2)
            sc = np.abs(ref).max()
            e = np.abs(got - ref).max() / sc
            worst_py = max(worst_py, e)
            if name == "uniform":
                assert np.array_equal(got, 2.0 * h1 - h2), "uniform steps must give 2 dS_{j-1} - dS_{j-2} exactly"
            if fco is not None:
                c1, c2 = fco[idx]
                if name == "uniform":
                    assert (c1, c2) == (2.0, -1.0), (c1, c2)
                ef = np.abs(c1 * h1 + c2 * h2 - ref).max() / sc
                worst_f = max(worst_f, ef)
        print(f"test_extrap_formula: {name:11s} steps {dts}: steps 3..{n} checked")
    print(f"test_extrap_formula: Python extrap_guess vs formula: max relative error {worst_py:.3e}")
    if fortran_ran:
        print(f"test_extrap_formula: Fortran extrap_coeffs vs formula: max relative error {worst_f:.3e}")
    else:
        print("test_extrap_formula: Fortran routine skipped (no working gfortran; set FC or activate the toolchain)")
    assert worst_py < 1e-14 and worst_f < 1e-14
    print("test_extrap_formula: PASS")


if __name__ == "__main__":
    main()
