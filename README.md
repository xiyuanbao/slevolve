<p align="center"><img src="docs/src/img/logo_small.png" alt="SLEvolve" width="320"></p>

# SLEvolve v1.0

<p align="center"><b>📖 Manual: <a href="https://xiyuanbao.github.io/slevolve/">https://xiyuanbao.github.io/slevolve/</a></b> — installation, quick start, benchmark cases, options, method<br>
<a href="https://doi.org/10.5281/zenodo.23422437"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.23422437.svg" alt="DOI"></a></p>

SLEvolve was developed by [Xiyuan Bao](https://xiyuanbao.github.io/). It combines an exact recursive reformulation of the viscoelastic memory, engineering optimizations of the solver, and physics-inspired iteration starts; the implementation was AI-agent-assisted, with the formulation, the algorithms and the acceptance criteria set by the author and every change verified against the legacy solver.
SLEvolve is a 1-D pseudo-spectral, gravitationally self-consistent sea-level solver for a spherically
symmetric viscoelastic Earth, based on the normal-mode (Love-number) formulation of Kendall et al.
(2005), with shoreline migration, marine ice grounding, rotational feedback, forward runs and ice age
runs (outer iteration on the initial topography). The viscoelastic memory is carried in a recursive
modal state instead of a sum over the full load history, so the total cost grows as O(N) instead of
O(N^2) in the number of time steps N. A Fortran solver derived from the legacy code of Han et al.
(2022) and a Python solver written from scratch; the paper's version (Zenodo v1.0.0, doi 10.5281/zenodo.23422438) reproduces the legacy solver to round-off (below 1e-11 m in the three
benchmark cases, with the Fortran and Python solvers identical to round-off); the current code additionally applies the corrected floating-ice convention and the corrected sign of the axial rotation component m3 (see
[Provenance](https://xiyuanbao.github.io/slevolve/provenance.html)).

## What is new relative to the legacy code

- Recursive modal state for all six viscous history terms (sea level, radial displacement, four
  degree-2 rotational terms), advanced by one decay factor per step; the recursion is exact.
- Placement rule: *advance-and-read* before the inner iteration, *commit* after convergence.
- Extrapolated inner start: the first guess of each step's ocean-load increment is extrapolated
  linearly from the two previous converged increments.
- Ice age runs: partially relaxed first-pass initial topography, warm-started outer passes, and an
  in-process outer-iteration driver.
- Faster spherical-harmonic transforms (on-the-fly Legendre recurrence or BLAS table, FFTW, N-S
  symmetry, OpenMP); binary `.npy` I/O with text <-> `.npy` converters in Fortran and Python.
- Horizontal surface displacement output (scalar potential Psi_lm, gradient on the grid).
- Floating-ice convention corrected: the ocean under floating ice shelves carries water load
  (grounded-ice mask of Kendall et al. 2005, eqs 43-44, flotation on the current topography).

Details: [Method overview](https://xiyuanbao.github.io/slevolve/method.html).

## Install and first run

```bash
git clone https://github.com/xiyuanbao/slevolve && cd slevolve
```

**Fortran** (gfortran, FFTW, OpenMP; BLAS from Apple Accelerate on macOS, OpenBLAS elsewhere):

```bash
conda create -n slevolve-fortran -c conda-forge gfortran fftw openblas llvm-openmp
conda activate slevolve-fortran
CASE=case1 bash fortran/build.sh     # -> fortran/build/case1/exe (about 10 s)
bash tests/test_fortran_synthetic.sh # Fortran-only test (no Python): vs reference/synthetic, "PASS"
```

`CASE=case1|case2|case3` selects a benchmark case's compiled-in settings; without `CASE` the build
goes to `fortran/build/custom/`. All switches: [Installation](https://xiyuanbao.github.io/slevolve/install.html) and
[Command-line options](https://xiyuanbao.github.io/slevolve/options.html).

**Python** (no install step needed; run the scripts in `python/`):

```bash
conda create -n slevolve -c conda-forge python=3.12 numpy numba ducc0 scipy netcdf4 matplotlib
conda activate slevolve
python examples/quickstart.py        # synthetic Earth, two melting ice sheets, Python solver, < 1 s
```

Or install the Python solver as the package `slevolve` (the wheel holds only the Python package):

```bash
pip install '.[fast]'                # numpy, scipy, numba, ducc0; editable: pip install -e '.[all]'
slevolve-forward --help              # = python python/run_ais.py; also slevolve-outer, slevolve-init-topo,
                                     #   slevolve-text2npy, slevolve-npy2text, slevolve-npy2nc
```

Tests:

- `python tests/run_all.py` runs the Python tests ("== 7/7 passed").
- `bash tests/test_fortran_python.sh` (Fortran vs Python, synthetic Earth),
  `bash tests/test_fortran_python_vardt.sh` (non-uniform steps) and
  `bash tests/test_fortran_python_ramp.sh` (optional second-order placement and rates) print "PASS".

See [Verifying the installation](https://xiyuanbao.github.io/slevolve/verify.html) and the tested environment in
[Installation](https://xiyuanbao.github.io/slevolve/install.html).

## Benchmark cases

The input data are not redistributed; scripts in `tools/` convert the archived files, and the
Love-number files of all three cases are in `lovenum/`.

- **Case 1**, Antarctic forward run (550 x 1 yr): ice from Yousefi et al. (2022),
  doi:10.5061/dryad.ncjsxkswm; Earth model and topography from Han et al. (2022), doi:10.17605/OSF.IO/8PTFM.
- **Case 2**, two glacial cycles (1200 x 0.2 kyr): ice history and topography from the OSF archive of
  Han et al. (2022), doi:10.17605/OSF.IO/8PTFM.
- **Case 3**, ICE-NADI-1 deglaciation (50 x 0.5 kyr; Lucas et al. 2026): ice history from Zenodo,
  doi:10.5281/zenodo.21464139.

Commands, expected results and run times: [Running the benchmark cases](https://xiyuanbao.github.io/slevolve/benchmarks.html);
data provenance: [Benchmark provenance](https://xiyuanbao.github.io/slevolve/provenance.html). Reference outputs of the three cases
are in `reference/`; check a run with `python tests/check_reference.py --case case1 --dir OUT_case1_py`.

## Known limitations

- SLEvolve reads Love numbers but does not compute them (`python/tools/write_lovenum.py` writes your
  own in the solver's layout).
- The Fortran grid is fixed at 512 x 1024 and the input file names are compiled in.
- Steps are uniform, or non-uniform from an optional times file ([Input files](https://xiyuanbao.github.io/slevolve/inputs.html)).

Full list: [Release notes](https://xiyuanbao.github.io/slevolve/release_notes.html).

## Documentation

Manual: https://xiyuanbao.github.io/slevolve/ (sources in `docs/`, MkDocs Material; `mkdocs serve`).

- **Getting started**: [installation](https://xiyuanbao.github.io/slevolve/install.html), [quick start](https://xiyuanbao.github.io/slevolve/quickstart.html),
  [verifying the installation](https://xiyuanbao.github.io/slevolve/verify.html).
- **User guide**: [benchmark cases](https://xiyuanbao.github.io/slevolve/benchmarks.html), [new case](https://xiyuanbao.github.io/slevolve/own_case.html),
  [Love numbers](https://xiyuanbao.github.io/slevolve/lovenum.html), [horizontal displacement](https://xiyuanbao.github.io/slevolve/horizontal.html),
  [file conversion](https://xiyuanbao.github.io/slevolve/conversion.html), [reproducing the published results](https://xiyuanbao.github.io/slevolve/reproduce.html).
- **Reference**: [options](https://xiyuanbao.github.io/slevolve/options.html), [inputs](https://xiyuanbao.github.io/slevolve/inputs.html),
  [outputs](https://xiyuanbao.github.io/slevolve/outputs.html), [formats](https://xiyuanbao.github.io/slevolve/formats.html),
  [code structure](https://xiyuanbao.github.io/slevolve/structure.html), [glossary](https://xiyuanbao.github.io/slevolve/reference/glossary.html).
- **Background**: [method](https://xiyuanbao.github.io/slevolve/method.html), [provenance](https://xiyuanbao.github.io/slevolve/provenance.html),
  [degree-1 Love numbers](https://xiyuanbao.github.io/slevolve/degree1.html), [release notes](https://xiyuanbao.github.io/slevolve/release_notes.html),
  [lineage and references](https://xiyuanbao.github.io/slevolve/lineage.html).

## Citation

If you use SLEvolve, please cite:

Bao, X., Mousavi, S. M., Latychev, K., and Mitrovica, J. X.: An exact finite-state representation of
viscoelastic memory for accelerated sea-level modeling: SLEvolve v1.0, Geoscientific Model
Development, submitted, 2026.

Software archive: Zenodo, https://doi.org/10.5281/zenodo.23422437 (all versions; v1.0.0:
https://doi.org/10.5281/zenodo.23422438). Machine-readable metadata: `CITATION.cff`.

## Lineage

SLEvolve builds on the sea-level theory of Kendall et al. (2005) and is a further development of the
code of Han et al. (2022); see [Lineage and references](https://xiyuanbao.github.io/slevolve/lineage.html).

- Kendall, R. A., Mitrovica, J. X., and Milne, G. A.: On post-glacial sea level - II. Numerical
  formulation and comparative results on spherically symmetric models, Geophys. J. Int., 161, 679-706,
  https://doi.org/10.1111/j.1365-246X.2005.02553.x, 2005.
- Han, H. K., Gomez, N., and Wan, J. X. W.: Capturing the interactions between ice sheets, sea level
  and the solid Earth on a range of timescales: a new "time window" algorithm, Geosci. Model Dev., 15,
  1355-1373, https://doi.org/10.5194/gmd-15-1355-2022, 2022.

## License

GPL-3.0-or-later (see `LICENSE`). The parts of `fortran/` derived from the code of Han et al. (2022)
remain under the MIT license of that code, reproduced in `LICENSE-upstream-MIT`.

## Contact

Please report problems and questions through GitHub issues at https://github.com/xiyuanbao/slevolve.
Author: [Xiyuan Bao](https://xiyuanbao.github.io/).
