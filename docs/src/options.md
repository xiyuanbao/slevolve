# Command-line options

This page lists in one place every build switch, environment variable and command-line option of the two solvers, with its default as set in the code.

## Fortran: `build.sh`

Usage: `fortran/build.sh [TPW] [CHECKMARINE] [NSTEPS]` (how to choose the switches:
[Installation](install.md#fortran)). The executable is `fortran/build/<CASE or custom><suffixes>/exe`;
a suffix is added only for a switch that differs from its default.

* Positional arguments (normally none is given):
  * `TPW` = `T|F`: rotational feedback on or off (default `T`).
  * `CHECKMARINE` = `T|F`: marine check on or off; with `T` ice that would float on the current
    topography is removed from the load (default `T`).
  * `NSTEPS`: number of steps compiled in (default from `CASE`, 550 without `CASE`); `L_sim = Ldt1 = DT * NSTEPS`.

### Case presets

`CASE=`; any variable set explicitly wins.

| `CASE` | sets | driver |
|---|---|---|
| `case1` | `TOPOMODEL=topo0_archive`, `NSTEPS=550`, `DT=1` | `run_from_zero_plus` |
| `case2` | `DT=200 NSTEPS=1200 INITIAL_TOPO=F PLANETMODEL=prem_512.l120C.ump5.lm5 ICEMODEL=iceload TOPOMODEL=etopo2_512_ice6gC` | `run_outer_plus` |
| `case3` | `DT=500 NSTEPS=50 INITIAL_TOPO=F PLANETMODEL=prem_512.l96C.ump3.lm1 ICEMODEL=nadi_ TOPOMODEL=topo_pres` | `run_outer_mem` |
| unset | the defaults of the table below (Case 1 file names, `NSTEPS=550`, `DT=1`); build directory `custom` | `run_from_zero_plus` |

### Build switches

| var | values | default | effect |
|---|---|---|---|
| `DRIVER` | `run_from_zero_plus`, `run_outer_plus`, `run_outer_mem` | from `CASE`, else `run_from_zero_plus` | the main program: forward run `./exe LAST` (uniform step = build-time `DT`, any length; non-uniform steps from a times file); one outer pass per process `./exe LAST [ITERSL]`; outer passes K0..K1 in one process `./exe LAST K0 K1`. `./exe` with wrong arguments prints this usage line |
| `NSTEPS` | integer | from `CASE`, else 550 | number of time steps of the run |
| `DT` | integer years | from `CASE`, else 1 | uniform step length `dt1..dt4`; must equal the time spacing of the ice files. With a times file (`INPUT_FILES/others/times` or `SLEVOLVE_TIMES`) the steps come from the file and `DT` only counts steps: build with `DT` unset |
| `INITIAL_TOPO` | `T`, `F` | from `CASE`, else `T` | `T`: forward run from the known topography `TOPO_INITIAL`; `F`: outer iteration on the initial topography, starting from `TOPOMODEL` |
| `PLANETMODEL`, `ICEMODEL`, `TOPOMODEL`, `TOPO_INITIAL` | file names | from `CASE`, else the Case 1 names; `TOPO_INITIAL = TOPOMODEL` | the input file names compiled in ([Input files](inputs.md)) |
| `NORDER` | integer | 511 | maximum spherical-harmonic degree; 512 (the legacy default) remains selectable |
| `EPS1` | float | `1.0E-5` | inner-iteration tolerance `epsilon1` |
| `LEGACY_CONSTANTS` | `1`, unset | unset | `1`: ice density 910 kg m^-3 and radius 6.371229e6 m of the older legacy code; unset: 920 and 6.371e6 (Han et al. 2022), used by all benchmark cases |
| `SHT` | `otf`, `blas`, `exact`, `ref` | `otf` | Legendre transform: `otf` computes the associated Legendre functions on the fly (fastest from about 4 threads, half the memory); `blas` uses a precomputed table and DGEMM; `exact` and `ref`: verification only (they reproduce the legacy transform) |
| `FFT` | `fftw`, `fft991` | `fftw` | longitude FFT: FFTW, or the routine of the legacy code (no FFTW needed) |
| `SYM` | `1`, `0` | `1` | `1` uses the north-south symmetry of the Legendre functions (fewer operations); `0` for checks only |
| `OMP` | `0`, `N` | `8` | any value other than 0 compiles with OpenMP (threads set at run time by `OMP_NUM_THREADS`; `N` itself is not used); `0` builds a serial executable |
| `BLASLIB` | `accelerate`, `openblas` | `accelerate` on macOS, else `openblas` | BLAS library linked |
| `TAG` | name | derived from `CASE` and the switches | build directory `fortran/build/<TAG>/` |
| `FC` | compiler | `gfortran` | Fortran compiler |
| `CONDA_ENV` | env name | unset | activate this conda env before building (`CONDA_BASE` defaults to `conda info --base`) |
| `FFTW_PREFIX` | dir | `$CONDA_PREFIX`, else `/usr/local` | FFTW 3 prefix (`lib/libfftw3`) |
| `OPENBLAS_PREFIX` | dir | `$CONDA_PREFIX`, else `/usr/local` | OpenBLAS prefix with `BLASLIB=openblas` |

## Fortran: run-time environment variables (read at run time by the executable)

The last column says whether the switch changes the numerical result.

| variable | values | default | effect | changes result? |
|---|---|---|---|---|
| `SLEVOLVE_TIMES` | path, or `off` | unset: `INPUT_FILES/others/times` if it exists | times file for non-uniform steps ([format](inputs.md#times-file-optional-non-uniform-time-steps)): one time in years per ice file `0..N`, strictly increasing; step `j` has length `t_j - t_(j-1)`, used for the modal-state decay factors, the extrapolated inner start and the time stamps in `OUTPUT_SLM/times`. `off` (or `none`, `0`) ignores the file, giving the uniform build-time `DT`. Read once at start-up | yes (different steps); a uniform times file `t_k = k * DT` gives the identical result |
| `SLEVOLVE_INNER` | `off`/`extrap`/`anderson` | `off` | start of the inner iteration at each step: `extrap` extrapolates the ocean-load increment linearly from the two previous steps (recommended; fewer iterations); `anderson` adds Anderson mixing (experimental) | within `epsilon1` |
| `SLEVOLVE_AA_M` | 1..5 | 3 | with `anderson`: number of previous iterates mixed | within `epsilon1` |
| `SLEVOLVE_AA_MASK` | `frozen`/`refresh` | `frozen` | with `anderson`: `refresh` continues iterating when the frozen ocean mask has become stale instead of stopping | within `epsilon1` |
| `SLEVOLVE_INNER_LOG` | 0/1 | 0 | writes the residual and the number of ocean-mask changes of every inner iteration to the log | no |
| `SLEVOLVE_REUSE` | `off`/`warm` | `off` | outer iteration: each pass starts every step from the converged load increment and ocean mask of the same step of the previous pass (recommended for ice age runs) | within `epsilon1` |
| `SLEVOLVE_OUTER` | `picard`/`secant` | unset | outer iteration: `picard` is the legacy update of the initial topography; `secant` uses a per-cell secant update from pass 3 (experimental). Setting either also writes the outer-iteration files `pred_pres_topo_<k>`, `tgrid0_<k>` in snapshot mode | `secant` changes the iterates |
| `SLEVOLVE_ICECACHE` | 0/1 | 1 | `run_outer_mem` only: keeps the ice grids in memory across passes instead of re-reading them | no |
| `SLEVOLVE_SNAP` | `all`/`none`/list (`1,2,5`, `1-20,50`) | unset = every step | writes the per-step grids only at the listed steps (step 0 is always added); `none` writes none | no |
| `SLEVOLVE_SNAP_FIELDS` | `tgrid[,ocean,beta,dS,R,G,horiz,Psi]`/`all` | `tgrid` | fields written at the steps of `SLEVOLVE_SNAP` | no |
| `SLEVOLVE_HORIZ` | 0/1/2 | 0 | horizontal displacement ([Horizontal displacement output](horizontal.md)): 1 writes `horiz<N>.npy` (Psi_lm), 2 also `Psi<N>.npy` | no (output only) |
| `SLEVOLVE_FMT` | `npy`/`bin` | `npy` | binary output container: `.npy` or headerless little-endian `.bin` | no |
| `SLEVOLVE_TEXT` | 0/1 | 0 | also writes the grids as legacy text files (slow, large) | no |
| `SLEVOLVE_BIN` | 0/1 | 1 | `0` suppresses the binary grid output | no |
| `SLEVOLVE_BININ` | 0/1 | 1 | `1` reads `<name>.npy` or `<name>.bin` before the text file `<name>`; `0` reads text only | no |
| `SLEVOLVE_CKPT` | K | 0 | writes the modal state to `exact_memory_state_<i>.bin` every K steps (0 = never) | no |
| `SLEVOLVE_OPT` | `all`/`none`/list | `all` | transform reductions, each replacing a spherical-harmonic transform by an equivalent cheaper operation: `diag` (degree-0 diagnostics by quadrature), `dsl` (total sea-level grid without a synthesis), `lazyrg` (skip the R synthesis when R and G are not written), `icebatch` (ice analysis in the first-guess sweep), `rt` (one analysis of rO - tO), `merge` (rO analysis and the next ocean-function guess in one sweep) | round-off only |
| `SLEVOLVE_BATCH` | 0/1 | 1 | batches independent transforms into one sweep; `0` runs them separately | no |
| `SLEVOLVE_FFT` | `fftw`/`fft991` | `fftw` in FFTW builds, else `fft991` | FFT used at run time (`fftw` in a build without FFTW stops with an error) | round-off only |
| `SLEVOLVE_ZERO` | `full`/`min` | `min` | `full` clears all working arrays at every step, as the legacy code; `min` clears only those that are read before being set | no |
| `SLEVOLVE_COMPAT` | bit flags | 0 | verification only: re-enables legacy per-step work (1 re-read step-0 fields, 2 re-read the Earth model, 4 save and reload the modal state, 8 recompute `t0lm`) | no |
| `SLEVOLVE_SHTPROF` | 0/1 | 0 | prints the time spent in each part of the transforms at the end | no |

Note: `SLEVOLVE_INNER=anderson` requires `SLEVOLVE_OPT` to include `merge` and `rt`.

**Threading and stack**

* `OMP_NUM_THREADS`: OpenMP threads.
* `VECLIB_MAXIMUM_THREADS`: Accelerate BLAS threads.
* `OMP_STACKSIZE=32M`: per-thread stack, set by the scripts in `examples/`; set it for your own runs too.

## Fortran: compiled-in parameters inherited from the legacy code (`user_specs_mod`)

| parameter | meaning | how to set |
|---|---|---|
| `tpw` | rotational feedback (true polar wander) | build arg 1 |
| `checkmarine` | flotation check of grounded ice | build arg 2 |
| `epsilon1` / `epsilon2` | inner convergence tolerance / outer tolerance (`epsilon2` is not used by the code) | `EPS1` (epsilon2 fixed at 1e-5) |
| `norder`, `nglv` | max degree, grid size (nglv=512 → 512×1024) | `NORDER` |
| `L_sim`, `dt1..dt4`, `Ldt1..Ldt4` | run length and step sizes | `DT` sets `dt1..dt4`; `NSTEPS` sets `L_sim = Ldt1 = DT*NSTEPS`. One build-time step size (`Ldt2..4 = 0`); these values index the steps (`L_sim/dt1` = number of steps). Non-uniform steps come from a times file at run time (`SLEVOLVE_TIMES`), not from `dt2..dt4` |
| `initial_topo` | forward run (`T`) or outer topography iteration (`F`) | `INITIAL_TOPO` |
| `planetmodel`, `icemodel`, `topomodel`, `topo_initial` | input file names | `PLANETMODEL`, `ICEMODEL`, `TOPOMODEL`, `TOPO_INITIAL` |
| `rhoi`, `radius` | ice density, Earth radius | default 920, 6.371E6; `LEGACY_CONSTANTS=1`: 910, 6.371229E6 |
| `iceVolume`, `calcRG` | legacy switches | fixed by the sed (`.true.`) |

`nglv` is fixed at 512 in the support module (`support_mod.f90`); `build.sh` has no switch for it.

## Python command lines

Each script has a console script with the same options after `pip install .`
([Installation](install.md#python-pip-install)): `python/run_ais.py` = `slevolve-forward`,
`python/run_outer.py` = `slevolve-outer`, `python/init_topo_segmented.py` = `slevolve-init-topo`,
`python/tools/text2npy.py` = `slevolve-text2npy`, `python/tools/npy2text.py` = `slevolve-npy2text`,
`python/tools/npy2nc.py` = `slevolve-npy2nc`.

Threads: `run_ais.py` and `run_outer.py` set `OMP_NUM_THREADS`, `OPENBLAS_NUM_THREADS`,
`VECLIB_MAXIMUM_THREADS`, `MKL_NUM_THREADS`, `NUMEXPR_NUM_THREADS` and `NUMBA_NUM_THREADS` to `--threads`
before NumPy is imported. With `--opt numba` the modal-state advance runs in the numba thread pool;
`run_ais.py --mem-threads` is then ignored (`solver_opt.py` sets `mem_threads = 1`), and `run_outer.py`
has no such option.

**`python/run_ais.py`** / `slevolve-forward` (forward run; any ice history `<icepfx>0..N` with uniform step `--dt` or the steps of a times file `--times`; the
defaults of `--icepfx`, `--earth`, `--topo` are the Case 1 names):

| option | default | meaning |
|---|---|---|
| `--steps` | 50 | number of steps N; ice files 0..N are read |
| `--dt` | 1 | step length in **years**; must equal the time spacing of the ice files. Printed at start-up and stored in `run_log.json` (`dt_years`) |
| `--times` | none | times file ([format](inputs.md#times-file-optional-non-uniform-time-steps)): time in years of ice files `0..N`, one per line, strictly increasing (the first `--steps + 1` values are used). Step `k` has length `t_k - t_(k-1)`; `--dt` is ignored. Stored in `run_log.json` (`times_file`, `times_years`; `dt_years` = `"variable"`) |
| `--threads` | 1 | threads of the transforms, BLAS and numba kernels (set before NumPy is imported) |
| `--mem-threads` | 1 | threads for the modal-state advance in the plain NumPy path (result unchanged; ignored with `--opt numba`) |
| `--sht` | `numpy` | transform backend: `numpy` (Legendre table, `sht.py`) or `ducc` (ducc0, much faster; recommended) |
| `--lmax` / `--norder` | 512 | maximum degree (ducc needs ≤ nglv−1; the benchmark cases use 511) |
| `--nglv` | 512 | number of Gauss-Legendre latitudes; the grid is nglv x 2 nglv |
| `--tpw`, `--checkmarine` | `T` | rotational feedback; marine check (as in Fortran) |
| `--calcRG` | `T` | computes the radial displacement R and geoid G grids (output only; `F` skips them) |
| `--epsilon1` | 1e-5 | inner-iteration tolerance |
| `--opt` | `''` | optimisations (round-off only), comma list or `all`: `raw` (ducc0 transforms without a separate rescaling pass), `diag` (degree-0 diagnostics by quadrature instead of four analyses), `dsl` (total sea-level grid without a synthesis), `numba` (compiled fused kernels; needs numba), `sums` (convergence sums carried between iterations and computed in parallel). Recommended: `all` |
| `--inner` | `off` | start of the inner iteration: `extrap` extrapolates the ocean-load increment from the two previous steps (recommended; result within `epsilon1`); `anderson` adds Anderson mixing (experimental) |
| `--anderson-depth`, `--anderson-freeze`, `--anderson-keep` | 3, 2, off | with `--inner anderson`: number of iterates mixed; plain iteration after which the ocean mask is frozen; keep the mixing history when the mask is refreshed |
| `--rhoi`, `--radius` | 920, 6.371e6 | ice density and Earth radius (Han et al. 2022) |
| `--icedir`, `--icepfx` | development cache / `AISload_gd_Run85_0_0_` | ice files `<icedir>/<icepfx><k>` (`.npy`, `.bin` or text), read directly. Pass explicitly |
| `--topo` | development tree | initial topography (`.npy`, `.bin` or text). Pass explicitly |
| `--earth` | development tree | Love-number file. Pass explicitly |
| `--save-fields` | `tgrid,ocean,beta` | fields written at every step: any of `tgrid,ocean,beta,dS,R,G` |
| `--rebuild-cache` | off | rebuild the binary ice cache (`python/cache/`). The cache is keyed by a JSON manifest next to it (ice directory, prefix, file count, size and mtime of every ice file) and is rebuilt automatically, with a message, when any of these change |
| `--horizontal` | off | write `horiz<N>.npy` (Psi_lm, [Horizontal displacement output](horizontal.md)) every step; output only |
| `--horizontal-grid` | off | also `utheta<N>.npy`, `uphi<N>.npy` (ducc0) |
| `--io-compat` | `fortran` | `fortran`: ocean/beta `<f8`, dS Fortran-order `<c16`; `legacy`: uint8 masks |
| `--out` | development path | output directory. Pass explicitly |
| `--root` | parent of `python/` | base directory of the development defaults; not needed when `--earth`, `--icedir`, `--topo` and `--out` are given |
| `--chunk` | 16 | block size of the batched transforms (performance only) |
| `--cache-only` | off | only build the ice cache in `python/cache/` (without `--icedir`), then stop |
| `--profile`, `--cprofile`, `--count-sht` | off | timing per stage; a cProfile report `cprofile.txt`; transform counts in the log. Log only |

**`python/run_outer.py`** / `slevolve-outer` (ice age run with outer iteration):

| option | default | meaning |
|---|---|---|
| `--iters` | 5 | number of outer passes |
| `--outer` | `picard` | update of the initial topography between passes: `picard` (legacy) or `secant`, a per-cell secant step from pass 3 (experimental; `--secant-tol` below which a cell uses Picard, `--smin`, `--smax` bounds of the step factor; defaults 1e-3 m, 0.5, 1.5) |
| `--reuse` | `off` | `warm`: each pass starts every step from the converged load increment and ocean mask of the previous pass (recommended; result within `epsilon1`). Experimental, not used in the paper: `warm+ice` also reuses the stored ice analysis; `warm+corr` adds the pass-to-pass change of the previous step; `warm+corr2` extrapolates that change linearly |
| `--inner` | `off` | as for `run_ais.py` |
| `--topo0` | `topo_pres` | T0 of pass 1 (text or `.npy`, e.g. a partially relaxed guess) |
| `--topo-pres` | `<inputs>/others/topo_pres` | observed present topography |
| `--inputs` | **required** | `INPUT_FILES` directory (`icemodel/`, `others/`) |
| `--icepfx`, `--nfiles` | `nadi_`, 51 | ice files `<icepfx>0 .. <icepfx><nfiles-1>` in `<inputs>/icemodel/`; `nfiles` = steps + 1 |
| `--ice-seq` | files `0..nfiles-1` | a custom ice history built from the cached files: comma list of `i` or `i:j:w` (the blend (1-w) I_i + w I_j), or `@FILE` |
| `--earth` | **required** | Love numbers |
| `--rebuild-cache` | off | rebuild the binary ice cache (`python/cache/`). The cache is keyed by a JSON manifest next to it (ice directory, prefix, file count, size and mtime of every ice file) and is rebuilt automatically, with a message, when any of these change |
| `--lmax` | **required** | truncation (511 for the benchmark cases) |
| `--dt` | 500 | step length in **years** (one fixed value; the spacing of the ice files). Printed at start-up |
| `--times` | none | times file: one time in years per ice state `0..N` (the ice files, or the `--ice-seq` entries), exactly N + 1 values, strictly increasing; step `k` has length `t_k - t_(k-1)`, `--dt` is ignored. Stored in `outer_log.json` (`times_years`) |
| `--rhoi`, `--radius` | 920, 6.371e6 | ice density and Earth radius |
| `--tpw`, `--checkmarine`, `--epsilon1` | T, T, 1e-5 | as for `run_ais.py` |
| `--threads` | 4 | ducc and numba threads (also `NUMBA_NUM_THREADS`) |
| `--opt` | `raw,diag,dsl,numba,sums` | as for `run_ais.py` |
| `--save-steps` | none | comma list of steps whose topography `tgrid<N>_it<k>.npy` is written in every pass |
| `--count-sht` | off | transform counts per pass in the log |
| `--out` | **required** | output directory (`outer_log.json`, `pred_pres_topo_<k>.npy`, `tgrid0_<k>.npy`, `T0_next.npy`, saved steps) |
| `--horizontal`, `--horizontal-grid` | off | `horiz<N>_it<k>.npy` (and `utheta/uphi<N>_it<k>.npy`) at `--save-steps` and the last step |

`run_outer.py` has no `--sht` and no `--mem-threads` option. It always uses ducc. With `--opt numba`
(default) the modal-state advance runs in the numba thread pool, sized by `--threads`; `solver_opt.py`
then sets `mem_threads = 1`, so `run_ais.py --mem-threads` has no effect with `--opt numba` either.

**`python/init_topo_segmented.py`** / `slevolve-init-topo` (partially relaxed, segmented first-pass initial topography for outer pass 1; defaults from the code):

| option | default | meaning |
|---|---|---|
| `--inputs` | **required** | `INPUT_FILES` directory (`icemodel/`, `others/`) |
| `--icepfx` | **required** | ice file prefix, e.g. `nadi_` |
| `--nfiles` | **required** | number of ice files (`0..nfiles-1`; last = present) |
| `--dt` | required unless `--times` | step length in **years** (converted to kyr internally) |
| `--times` | none | times file, one time in years per ice file (`nfiles` values, strictly increasing); the segment ages are `t_N - t_n` (instead of `dt (N - n)`); `--dt` is ignored |
| `--earth` | **required** | Love-number file |
| `--lmax` | **required** | truncation |
| `--topo-pres` | `<inputs>/others/topo_pres` | present topography |
| `--nglv` | 512 | number of Gauss-Legendre latitudes; the grid is nglv x 2 nglv |
| `--threshold` | 0.01 | monotonic legs of the (smoothed) ice-volume history whose volume change is below this fraction of the total range are merged into their neighbours |
| `--start` | present topography | starting T0 of the fixed point |
| `--c0` | 62.0 | starting uniform term c (m) |
| `--rhoi`, `--rhow`, `--radius` | 920, 1000, 6.371e6 | constants |
| `--workers` | 4 | ducc threads |
| `--out` | **required** | output `.npy`; a `.json` with the segments is written next to it |

**Next:** [Input files](inputs.md).
