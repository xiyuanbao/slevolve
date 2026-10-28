# Installation

This page lists the tested environment and how to build the Fortran solver and set up the Python solver.

## Tested environment

Only the environment below has been tested.

### Fortran

| item | version | how recorded |
|---|---|---|
| machine | Apple M4 Max, macOS 15.1 (24B2082) | `sysctl -n machdep.cpu.brand_string`, `sw_vers` |
| compiler | GNU Fortran (conda-forge gcc 16.2.0-5) 16.2.0 | `gfortran --version` |
| libgfortran | 16.2.0 (conda-forge) | `conda list` |
| FFTW 3 | 3.3.11 (conda-forge), `libfftw3.3.dylib` | `conda list`; `otool -L exe` |
| OpenMP runtime | llvm-openmp 23.1.1 (`libomp.dylib`, conda-forge) | `conda list`; `otool -L exe` |
| BLAS | Apple Accelerate framework (system) | `otool -L exe` |

`otool -L` of a `CASE=case1` build lists `libfftw3.3.dylib`, `libgfortran.5.dylib`,
`Accelerate.framework`, `libomp.dylib`, `libquadmath.0.dylib`, `libSystem.B.dylib`.

### Python

| package | version | used for |
|---|---|---|
| Python | 3.12.14 | |
| numpy | 2.5.3 | everything |
| numba | 0.67.0 | `--opt numba` kernels (`kernels.py`), state advance threading |
| ducc0 | 0.41.0 | `--sht ducc`, `run_outer.py`, `init_topo_segmented.py` transforms (`sht_ducc.py`), `horiz.py` |
| scipy | 1.18.1 | `sht.py` (`--sht numpy` backend: `scipy.fft`; without scipy it falls back to `numpy.fft`), `tools/regrid.py` |
| matplotlib | 3.11.2 | optional: the figures of `examples/quickstart.py`, `tools/plot_reference.py` |
| netCDF4 | 1.7.4 | optional: `python/tools/npy2nc.py` and the input-preparation scripts |

## Get the code

```bash
git clone https://github.com/xiyuanbao/slevolve && cd slevolve
```

All commands in this manual are run from this directory (the release root).

## Fortran

**Toolchain.** gfortran, FFTW 3, OpenMP and (outside macOS) OpenBLAS, all from conda-forge:

```bash
conda create -n slevolve-fortran -c conda-forge gfortran fftw openblas llvm-openmp
conda activate slevolve-fortran
```

This environment was created and used to build and pass the Fortran test on the tested machine. On
macOS the build links Apple Accelerate for BLAS; elsewhere it links the OpenBLAS of the environment
(`BLASLIB=openblas`, the default on Linux; not tested). Instead of activating the environment you can
pass `CONDA_ENV=slevolve-fortran` to the build script.

**Build.**

```bash
CASE=case1 bash fortran/build.sh
```

takes about 10 s and ends with the line

```
built: /path/to/slevolve/fortran/build/case1/exe (case=case1 driver=run_from_zero_plus sht=otf ...)
```

The tools `make_cache`, `text2npy` and `npy2text` are built into the same directory. The linker may
print `ld: warning: duplicate -rpath ... ignored` (once per linked program); this is harmless.
Running `./exe` without arguments prints the arguments the driver expects, e.g.
`usage: ./exe LAST   (driver run_from_zero_plus: forward run of steps 1..LAST; the step length is the build-time DT)`.
To check the build, run the Fortran test ([Verifying the installation](verify.md#fortran-test)).

**Build directory.** `fortran/build/<CASE>/` for `CASE=case1|case2|case3`, `fortran/build/custom/`
without `CASE`. A switch that differs from its default adds a suffix (`SHT=blas` gives
`case1_blas`, `OMP=0` gives `case1_serial`, `NSTEPS=20` gives `case1_N20`); `TAG=<name>` sets the
name directly.

**Positional arguments** `fortran/build.sh [TPW] [CHECKMARINE] [NSTEPS]`. Normally none is given.

| argument | values | default | meaning |
|---|---|---|---|
| `TPW` | `T`, `F` | `T` | rotational feedback (true polar wander) |
| `CHECKMARINE` | `T`, `F` | `T` | marine check: ice that would float on the current topography iterate is removed from the load |
| `NSTEPS` | integer | from `CASE` (550 without `CASE`) | number of time steps compiled in (`L_sim = DT * NSTEPS`); also settable as `NSTEPS=` |

**Switches and when to choose them** (environment variables; full table in
[Command-line options](options.md#fortran-buildsh)):

| switch | default | when to change it |
|---|---|---|
| `CASE=case1\|case2\|case3` | unset | to build one of the benchmark cases: sets the file names, `DT`, `NSTEPS`, `INITIAL_TOPO` and the driver of that case. Unset for your own case, with the file names given explicitly ([Setting up a new case](own_case.md#6-running)) |
| `DRIVER` | from `CASE`; `run_from_zero_plus` | `run_from_zero_plus`: forward run, `./exe LAST`, step length `DT`. `run_outer_plus`: one outer pass per process, `./exe LAST [ITERSL]`. `run_outer_mem`: several outer passes in one process, `./exe LAST K0 K1` (fastest for ice age runs) |
| `DT` | from `CASE`; 1 | step length in integer years; must equal the time spacing of your ice files. The forward driver uses the build-time `DT` (any integer); non-uniform steps come from a times file ([Input files](inputs.md#times-file-optional-non-uniform-time-steps)) |
| `INITIAL_TOPO` | from `CASE`; `T` | `T`: forward run from a known initial topography (`TOPO_INITIAL`); `F`: ice age run, outer iteration on the initial topography from the present topography `TOPOMODEL` |
| `PLANETMODEL`, `ICEMODEL`, `TOPOMODEL`, `TOPO_INITIAL` | from `CASE`; Case 1 names | your file names in `INPUT_FILES/earthmodel/`, `icemodel/` (prefix), `others/` |
| `NORDER` | 511 | maximum degree. 511 is the largest degree resolved by the 512 Gauss latitudes and is used by all benchmark cases; 512 (the legacy default) remains selectable |
| `SHT` | `otf` | Legendre transform. `otf` (on-the-fly recurrence) is the fastest from about 4 threads up and needs half the memory of `blas`; `blas` (precomputed table and DGEMM) can be faster on 1-2 threads. `exact` and `ref` are for verification only (they reproduce the legacy transform) |
| `OMP` | 8 | any value other than 0 compiles with OpenMP; the number of threads is set at run time with `OMP_NUM_THREADS`. `OMP=0` builds a serial executable |
| `FFT` | `fftw` | `fft991` uses the FFT routine of the legacy code (no FFTW needed; slower) |
| `SYM` | 1 | north-south symmetry in the transform (faster; same result to round-off). `SYM=0` for checks only |
| `BLASLIB` | `accelerate` on macOS, else `openblas` | the BLAS library (`OPENBLAS_PREFIX` if not in the conda env) |
| `EPS1` | 1.0E-5 | inner-iteration tolerance epsilon1 |
| `LEGACY_CONSTANTS=1` | unset | ice density 910 kg m^-3 and radius 6.371229e6 m of the older legacy code, instead of 920 and 6.371e6 (Han et al. 2022) |
| `RAMP=1` | unset (`0`) | second-order load placement: each load increment accrues linearly over its step instead of being applied at the end of it ([Method](method.md#second-order-load-placement-optional)); directory suffix `_ramp`. Leave unset to reproduce the paper and `reference/`. The analytic rates are a run-time switch, `SLEVOLVE_RATES=1`, in either build |
| `FC`, `FFTW_PREFIX`, `OPENBLAS_PREFIX`, `CONDA_ENV` | `gfortran`, `$CONDA_PREFIX` | toolchain without an activated conda env (not tested) |

Compiler flags: `-O3 -ffree-line-length-none -fdefault-real-8 -fdefault-double-8 -cpp`.

**Build all three benchmark cases**

```bash
CASE=case1 bash fortran/build.sh    # -> fortran/build/case1/exe   forward run, ./exe 550
CASE=case2 bash fortran/build.sh    # -> fortran/build/case2/exe   one outer pass per process, ./exe 1200 1
CASE=case3 bash fortran/build.sh    # -> fortran/build/case3/exe   in-process outer iteration, ./exe 50 K0 K1
```

## Fortran without Python

The Fortran solver, its tools and the [Fortran-only test](verify.md#fortran-only-test) need no Python.
The ice age wrapper `examples/run_outer_mem.sh` needs no Python either. The segmented first-pass
initial topography (`python/init_topo_segmented.py`) and the input-preparation scripts in `tools/` do; without
Python, start the outer iteration from the present topography (`TOPO0=pres`) or supply your own `.npy`.
Instead of conda, the toolchain can come from the system package manager. **Only the macOS/conda route
above was tested here**; the commands below are given for convenience and have not been run on this
release.

**macOS (Homebrew):**

```bash
brew install gcc fftw openblas libomp
FC=gfortran-15 FFTW_PREFIX=$(brew --prefix fftw) \
  NSTEPS=8 ICEMODEL=synth_ice_ TOPOMODEL=synth_topo PLANETMODEL=prem_512.l96C.ump3.lm1 \
  bash fortran/build.sh
```

Homebrew installs gfortran with its major version in the name (`gfortran-15` or similar; check
`ls $(brew --prefix)/bin/gfortran*`). On macOS the build links Apple Accelerate for BLAS; to use the
Homebrew OpenBLAS instead add `BLASLIB=openblas OPENBLAS_PREFIX=$(brew --prefix openblas)`. gfortran
brings its own OpenMP runtime (`libgomp`); `libomp` is only needed if you build with clang-based
tools.

**Debian/Ubuntu (apt):**

```bash
sudo apt install gfortran libfftw3-dev libopenblas-dev
FC=gfortran FFTW_PREFIX=/usr BLASLIB=openblas OPENBLAS_PREFIX=/usr \
  NSTEPS=8 ICEMODEL=synth_ice_ TOPOMODEL=synth_topo PLANETMODEL=prem_512.l96C.ump3.lm1 \
  bash fortran/build.sh
```

`build.sh` links `-L$FFTW_PREFIX/lib -lfftw3` and `-L$OPENBLAS_PREFIX/lib -lopenblas`. On Debian/Ubuntu
the libraries are in `/usr/lib/x86_64-linux-gnu` (or `/usr/lib/aarch64-linux-gnu`), which the linker
searches by default, so `/usr` works even though `/usr/lib` itself does not hold them.

| variable | meaning | Homebrew | apt |
|---|---|---|---|
| `FC` | Fortran compiler | `gfortran-15` (versioned name) | `gfortran` |
| `FFTW_PREFIX` | prefix of `lib/libfftw3` | `$(brew --prefix fftw)` | `/usr` |
| `BLASLIB` | `accelerate` or `openblas` | `accelerate` (default on macOS) | `openblas` (default elsewhere) |
| `OPENBLAS_PREFIX` | prefix of `lib/libopenblas` (with `BLASLIB=openblas`) | `$(brew --prefix openblas)` | `/usr` |

Then run the Fortran-only test with the same variables in the environment:

```bash
FC=... FFTW_PREFIX=... bash tests/test_fortran_synthetic.sh
```

It builds the solver, generates the synthetic Earth with the Fortran generator
`tools/make_synthetic_case.f90` (same Earth and files as the Python generator), runs 8 steps and
compares the result with the shipped reference `reference/synthetic/tgrid8.npy`.

## Python

* Python 3 with `numpy`, `numba`, `ducc0` and `scipy` (for the `--sht numpy` backend). Optional:
  `netCDF4` for `python/tools/npy2nc.py` and the input-preparation scripts. Tested versions under
  [Tested environment](#tested-environment).
* Two routes, same code and options: run the scripts in `python/` directly (no install step), or
  install the package `slevolve` with pip ([Python: pip install](#python-pip-install)) and use the
  console scripts.
* `run_ais.py` and `run_outer.py` pin thread counts before NumPy is imported.
* The default SHT backend of `run_ais.py` is `numpy` (table, `sht.py::SHT`), and its optimisations
  and the extrapolated inner start are off by default. For production runs use
  `--sht ducc --opt all --inner extrap` (the configuration of all benchmark runs).
* `run_outer.py` and `init_topo_segmented.py` always use ducc0.

Environment:

```bash
conda create -n slevolve -c conda-forge python=3.12 numpy numba ducc0 scipy netcdf4 matplotlib
conda activate slevolve
```

Check it with `python tests/run_all.py` ([Verifying the installation](verify.md#python-tests)).

## Python: pip install

The Python solver is also the package `slevolve` (sources in `python/slevolve/`, metadata in
`pyproject.toml` at the release root). From a clone:

```bash
pip install .                 # numpy, scipy
pip install '.[fast]'         # plus numba and ducc0 (needed for --opt numba, --sht ducc, run_outer, init-topo)
pip install -e '.[all]'       # editable install for development; all = fast + plot (matplotlib) + netcdf (netCDF4)
```

Extras: `fast` (numba, ducc0), `plot` (matplotlib), `netcdf` (netCDF4), `all`. Python >= 3.10.

Console scripts (same options as the scripts; `--help` lists them):

| console script | script in the repository |
|---|---|
| `slevolve-forward` | `python/run_ais.py` |
| `slevolve-outer` | `python/run_outer.py` |
| `slevolve-init-topo` | `python/init_topo_segmented.py` |
| `slevolve-text2npy` | `python/tools/text2npy.py` |
| `slevolve-npy2text` | `python/tools/npy2text.py` |
| `slevolve-npy2nc` | `python/tools/npy2nc.py` |

For example the Case 1 forward run of [Benchmarks](benchmarks.md) becomes
`slevolve-forward --steps 550 --sht ducc --opt all --inner extrap --lmax 511 --icedir ... --topo ... --earth lovenum/... --out ...`.
In Python: `import slevolve` (`slevolve.__version__`, `slevolve.SeaLevelSolver`, `slevolve.OptSeaLevelSolver`,
`slevolve.LoveNumbers`, `slevolve.Planet`, `slevolve.SHT`, `slevolve.make_sht`; loaded on first use, so
`import slevolve` does not import NumPy).

The wheel contains only the Python package. The Fortran solver, the Love-number files (`lovenum/`), the
tests, the examples, the input-preparation tools (`tools/`) and the reference data (`reference/`) are in
the repository, not in the wheel; pass the Love-number file with `--earth`. Ice caches built by
`slevolve-forward` without `--icedir` and by `slevolve-outer` go to `python/cache/` in a checkout
(also for an editable install), else to `~/.cache/slevolve`; `SLEVOLVE_CACHE_DIR` overrides both.

Manual: built with mkdocs 1.6.1 and mkdocs-material 9.7.7 (`mkdocs build --strict` in `docs/`).

**Next:** [Quick start](quickstart.md): a first run in under a second.
