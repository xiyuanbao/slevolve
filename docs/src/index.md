# SLEvolve v1.0

![SLEvolve](img/logo.svg){ width="360" }

SLEvolve was developed by [Xiyuan Bao](https://xiyuanbao.github.io/). It combines an exact recursive reformulation of the viscoelastic memory, engineering optimizations of the solver, and physics-inspired iteration starts; the implementation was AI-agent-assisted, with the formulation, the algorithms and the acceptance criteria set by the author and every change verified against the legacy solver.

## What SLEvolve solves

SLEvolve solves the gravitationally self-consistent sea-level equation on a spherically symmetric
viscoelastic Earth, including shoreline migration, marine ice grounding, rotational feedback and the
outer iteration on the initial topography of ice age runs. The explicit load-history sum of the
pseudo-spectral formulation is replaced by an exact recursive modal state, so the cost per time step
does not grow with the length of the load history and the total cost grows as O(N) instead of O(N^2)
in the number of time steps N. A Fortran solver derived from the legacy code of Han et al. (2022) and
a Python solver written from scratch give results identical to round-off (below 1e-10 m in the benchmark cases).

## Install and first run

Get the code:

```bash
git clone https://github.com/xiyuanbao/slevolve && cd slevolve
```

**Python** (no install step; the quick start runs the Python solver):

```bash
conda create -n slevolve -c conda-forge python=3.12 numpy numba ducc0 scipy netcdf4 matplotlib
conda activate slevolve
python examples/quickstart.py
```

**Fortran** (build, then the Fortran test, which also needs the Python env above for the comparison):

```bash
conda create -n slevolve-fortran -c conda-forge gfortran fftw openblas llvm-openmp
conda activate slevolve-fortran
CASE=case1 bash fortran/build.sh            # last line: built: .../fortran/build/case1/exe (...)
PYTHON="$(conda run -n slevolve which python)" bash tests/test_fortran_python.sh   # prints PASS
```

**Fortran only** (no Python at all; toolchain from conda, Homebrew or apt, see
[Fortran without Python](install.md#fortran-without-python)):

```bash
bash tests/test_fortran_synthetic.sh        # builds, generates the synthetic Earth in Fortran, prints PASS
```

Details: [Installation](install.md), [Verifying the installation](verify.md).

## How to cite

Bao, X., Mousavi, S. M., Latychev, K., and Mitrovica, J. X.: An exact finite-state representation of
viscoelastic memory for accelerated sea-level modeling: SLEvolve v1.0, Geoscientific Model
Development, submitted, 2026.

Equation, table and section numbers in these pages refer to that paper. Source code:
<https://github.com/xiyuanbao/slevolve>.

## The manual

- **[Getting started](install.md)**: installation, a first run, and checking the installation.
- **[User guide](benchmarks.md)**: running the benchmark cases, setting up a new case, Love numbers,
  horizontal displacement, file conversion and reproducing the published results.
- **[Reference](options.md)**: every option and switch, input and output files, file formats, code
  structure and a glossary.
- **[Background](method.md)**: the method, benchmark provenance, degree-1 Love numbers, release notes,
  and the [lineage and references](lineage.md).

**Next:** [Installation](install.md).
