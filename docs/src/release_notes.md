# Release notes and known limitations

This page lists the contents of v1.0, its differences from the legacy code and its known limitations.

## v1.0

**Note of 2026-10-30: sign of the axial rotation component m3 (correction).** The axial component
of the linearised Liouville equation used +1/C; the correct factor is -1/C,
\(m_3=-[(1+k^E_2)\,I_{33}+\text{viscous terms}]/C\) (Mitrovica & Wahr 2011, Annu. Rev. Earth Planet.
Sci. 39, 577-616, eq. 15; the same sign is in J. X. Mitrovica's sea-level code). A load that moves mass
towards the rotation axis now speeds up the spin (\(m_3>0\)). The error was inherited from the rotation
block of the legacy code. Both solvers are corrected, including the analytic rate of m3 in
`TPWrate<N>`; the equatorial (polar-motion) equations and the rotational potential were already
correct. Size of the change: \(\dot J_2\) about 1 % (Case 3 at the present, from the geoid rate:
-7.139e-12 to -7.062e-12 /yr), the sign of the length-of-day rate reversed (Case 3: +0.124 to
-0.124 ms per century), the polar-motion rate changed by a relative 2e-6, and relative sea level at the
millimetre level (max |change| of the reference fields: Case 1 6.4 mm, Case 2 4.1 mm, Case 3 0.9 mm;
0.15 m on the 220 m of the 8-yr synthetic test case). With the correction, \(-3C\dot m_3/(2Ma^2)\) and
\(\dot J_2\) read from the geoid agree in sign and to 0.4 %. `reference/` and the expected outputs of the
tests and of the quick start were recomputed; the documentation figures do not change visibly and
were kept. The legacy code has the same sign, so the round-off reproduction of the legacy code in the
v1.0.0 paper is unaffected as a statement about that code. See
[Benchmark provenance](provenance.md#differences-from-the-runs-of-han-et-al-2022).

**Note of 2026-10-28: optional second-order load placement and analytic rates.** Two opt-in
additions, both off by default, so every default output, `reference/` and the results of the paper
are unchanged:

- *Second-order load placement* (Fortran build switch `RAMP=1`, Python `--ramp` for `run_ais.py` and
  `run_outer.py`): each load increment accrues linearly over its step instead of being applied at its
  end. Each mode picks up the factor phi(s dt) = (1 - exp(-s dt)) / (s dt) in its memory and the
  same-step weight 1 - phi; the elastic and direct sea-level terms stay at the end of the step, and the
  polar-motion equation is solved in closed form with the implicit current increment. The time
  stepping becomes second order where the load is smooth in time; works with non-uniform steps, the
  extrapolated inner start and the warm start. See
  [Method](method.md#second-order-load-placement-optional).
- *Analytic rates* (Fortran `SLEVOLVE_RATES=1`, Python `--save-fields ...,dSLdt,dRdt` or
  `run_outer.py --rates`): the rates of sea level and radial displacement at the end of a step
  (`dSLdt<N>`, `dRdt<N>`, m/yr) and of the rotational quantities (`TPWrate<N>`), read from the modal
  states instead of a two-point difference; output only. Available in either build, but more
  accurate than a two-point difference only together with the ramp. See
  [Output fields](outputs.md#analytic-rates-optional).

The placement was derived and verified in a separate study. New tests: `tests/test_ramp.py` (in
`run_all.py`) and `tests/test_fortran_python_ramp.sh` ([Verifying the installation](verify.md#second-order-load-placement-and-rates)).

**Note of 2026-10-26: floating-ice convention.** The ice function \(\beta\) is now the grounded-ice
mask of Kendall et al. (2005, eqs 43-44): the ocean under floating ice shelves carries water load and
belongs to the ocean area of the conservation term, and the flotation check uses the current
topography iterate of the inner iteration (mask updated at most three times per step). The legacy
code, and the code of the paper (Zenodo v1.0.0, doi 10.5281/zenodo.23422438), excluded every
ice-covered cell, floating shelves included, from the water-loading ocean and tested flotation
against the initial topography of the pass. The change applies to both solvers and has no switch;
`reference/`, the expected outputs of the tests and the documentation figures were recomputed with
it. New test: `tests/test_floating_ice.py`. See
[Benchmark provenance](provenance.md#differences-from-the-runs-of-han-et-al-2022).

**Contents of the release**

- `fortran/`: the Fortran solver (`src/`), the build script `build.sh` with presets for the three
  benchmark cases, and the tools `make_cache`, `text2npy` and `npy2text` (`tools/`).
- `python/`: the Python solver (`run_ais.py` forward run, `run_outer.py` ice age run with outer
  iteration, `init_topo_segmented.py` first-pass initial topography) and the converters and
  post-processing tools in `python/tools/`.
- `pyproject.toml`: the Python solver as the pip-installable package `slevolve` 1.0.0 (sources in
  `python/slevolve/`, subpackage `slevolve.tools`; extras `fast`, `plot`, `netcdf`, `all`) with the
  console scripts `slevolve-forward`, `slevolve-outer`, `slevolve-init-topo`, `slevolve-text2npy`,
  `slevolve-npy2text`, `slevolve-npy2nc`. The scripts in `python/` remain the no-install route.
- `tools/`: input preparation for the benchmark cases (`build_case1_total_ice.py`, `regrid.py`,
  `prepare_case2_inputs.sh`, `prepare_case3_nadi.py`). No input data are redistributed.
- `lovenum/`: the Love-number files of the three cases, in the degree-1 frame CE (k_1 = 0), and
  `DEGREE1_NOTES.md` on degree-1 Love numbers. Case 1 (`prem_512.l60K2C.sum18p6.dum19p2.tz19p4.lm22`)
  is the file distributed with Han et al. (2022), included unchanged; it is CE in its elastic and fluid
  limits but carries a degree-1 initialisation defect in its residues, which changes the degree-1
  R, G and horizontal outputs but not sea level or topography. Cases 2 and 3
  (`prem_512.l120C.ump5.lm5`, `prem_512.l96C.ump3.lm1`) were computed with the collocation method of
  Mitrovica and Peltier (1992) (after Peltier 1985).
- `tools/make_grid.py` (grid files `GLlat_512.txt`, `GLlon_512.txt`) and
  `python/tools/write_lovenum.py` (Love numbers of your own Earth model in the solver's layout).
- `reference/`: reference outputs of the three cases, checked with `tests/check_reference.py`.
- `tests/`: the Python tests `python tests/run_all.py` (recursion vs explicit history sum, transform
  round trip, idealised small case, the extrapolated inner start for uniform and non-uniform steps,
  the small case with non-uniform steps), the Fortran test `tests/test_fortran_python.sh` (Fortran vs
  Python on a synthetic Earth), its non-uniform-step twin `tests/test_fortran_python_vardt.sh`
  (times file; also checks that a uniform times file reproduces the run without one exactly), the Fortran-only test `tests/test_fortran_synthetic.sh` (no Python;
  against `reference/synthetic/tgrid8.npy`) and `tests/check_reference.py`.
- `tools/make_synthetic_case.py` and its Fortran twin `tools/make_synthetic_case.f90` (synthetic-Earth
  input tree: an Earth-like random topography from a fixed seed, shipped as
  `reference/synthetic/topography_seed.npy`, and two melting Vialov ice sheets),
  `python/tools/plot_fields.py` (Mollweide and south polar stereographic maps with coastlines) and
  `tools/plot_reference.py` (maps of the reference results).
- `examples/`: the quick start `quickstart.py`, run scripts for the three benchmark cases
  (`case1_forward.sh`, `case2_one_pass.sh`, `case3_ice_age.sh`, each for Fortran or Python) and the
  in-process outer-iteration wrapper `run_outer_mem.sh`.
- `docs/`: this manual.

**Differences from the legacy code of Han et al. (2022)**

- The viscous history terms (sea level, radial displacement and the four degree-2 rotational terms)
  are carried in a recursive normal-mode state, advanced by one decay factor per step, instead of an
  explicit sum over the load history. The recursion is exact.
- One process runs all time steps with the Earth model, modal state and transform tables resident;
  binary `.npy` input and output, with text <-> `.npy` converters in Fortran and Python.
- Fewer spherical-harmonic transforms per step; Legendre transform by an on-the-fly recurrence or a
  BLAS table; FFTW; north-south symmetry; OpenMP.
- Optional extrapolated start of the inner iteration and warm start of outer passes (these change only
  the starting point of an iteration; the converged result agrees within the inner tolerance).
- Ice age runs: partially relaxed (segmented) first-pass initial topography, and an in-process
  outer-iteration driver.
- Optional output of the horizontal surface displacement.
- The longitude FFT routine `fft99b` of the transform module is used in a corrected form.

The code lineage is described in [Lineage and references](lineage.md).

**Known limitations**

* **Ice-sheet coupling.** The ice-sheet coupling interface inherited from the legacy code (the `coupling` flag and the coupled
  ice-input path in `fortran/src/sl_model_mod.f90`) is retained in the Fortran source but is
  not used by any benchmark case, is absent from the Python version, and has not been tested with the
  recursive state in v1.0.
* **Time steps.** Uniform steps (the benchmark cases use 1 yr, 200 yr and 500 yr): Python `--dt` in
  years (`run_ais.py`, default 1; `run_outer.py`; `init_topo_segmented.py`), printed at start-up;
  Fortran compiled in by `build.sh DT=<years>` (the forward driver `run_from_zero_plus` uses the build-time `DT`,
  verified against Python with `DT=1` and `DT=5`). Non-uniform steps: an optional times file, one time in years per ice file
  ([format](inputs.md#times-file-optional-non-uniform-time-steps)), read by Fortran at run time from
  `INPUT_FILES/others/times` or `SLEVOLVE_TIMES=<path>` (build with `DT` unset) and by Python with
  `--times <file>`. The benchmark cases use uniform steps; the non-uniform mode is verified on the
  idealised small case and the synthetic Earth (`tests/test_small_case_vardt.py`,
  `tests/test_fortran_python_vardt.sh`), not on a benchmark case.

  What the times file changes (step `j` has length `dt_j = t_j - t_(j-1)`):

    1. the decay factors of the modal-state update, `d_k = exp(-s_k dt_j)`, recomputed whenever the
       step length changes (cached otherwise), in both solvers; the recursion remains exact for the
       step-wise load history (`tests/test_recursion.py`, `tests/test_small_case_vardt.py`);
    2. the extrapolated inner start: from step 3 on,
       `dS0_j = dt_j [Sdot_(j-1) + (Sdot_(j-1) - Sdot_(j-2)) (dt_j + dt_(j-1)) / (dt_(j-1) + dt_(j-2))]`,
       `Sdot_i = dS_i / dt_i`, which is exactly `2 dS_(j-1) - dS_(j-2)` for equal steps
       (`tests/test_extrap_formula.py`);
    3. the time stamps: Fortran `OUTPUT_SLM/times` holds `t_0 .. t_N` (written as `ES24.16E3` instead of
       the uniform-step `ES14.4E2`, so that non-integer times are kept); Python `run_log.json` and
       `outer_log.json` store `times_years` (`run_log.json` `dt_years` = `"variable"`);
    4. Python `init_topo_segmented.py --times`: the segment ages are `t_N - t_n` instead of `dt (N - n)`.

  Unchanged: the step indexing (ice file k is read at step k; the Fortran `NSTEPS`, `L_sim/dt1`, the
  time window `TIMEWINDOW`, the checkpoint cadence and the `--save-steps`/`SLEVOLVE_SNAP` step numbers
  count steps, not years). Without a times file the code paths are those of uniform steps, and a
  uniform times file reproduces them exactly. The legacy `input_times`/`timearray` flag of
  `user_specs_mod` is not used (it still stops the exact-memory path when set).
* **Fortran grid and names.** `nglv` is fixed at 512, and the input file names are compiled in: every
  new case needs a rebuild.
* **Love numbers** are read, not computed ([Love numbers](lovenum.md)).
* **Python defaults.** `run_ais.py` has development defaults for `--earth`, `--icedir`, `--topo` and
  `--out`, and `run_outer.py` the Case 3 defaults `--icepfx nadi_ --nfiles 51 --dt 500`; always pass
  them explicitly.
* **Converters.** Self-describing metadata (grid, units, time axis) exist only via `npy2nc.py`. The
  checkpoint and `reuse_<k>.bin` record files are not converted to `.npy`.
* **Build portability.** Only the macOS / Apple-silicon configuration of [Installation](install.md#tested-environment)
  has been tested; the non-conda (`FFTW_PREFIX`) and OpenBLAS paths have not.
* Anderson mixing, the secant outer update and the `--reuse` variants other than `warm` are
  experimental and were not used in the paper.
* The example scripts in `examples/` are templates for the benchmark cases; they expect the input
  trees produced by the scripts in `tools/`.

## Names changed from the development versions

Files, environment variables and identifiers were renamed for the release; scripts written for the
development versions need these changes.

| development name | release name |
|---|---|
| `fortran/build6.sh`, `fortran/src6/` | `fortran/build.sh`, `fortran/src/` |
| build directory `fortran/build/pp6_case1_c2021_otf_fftw_sym1_omp8_n511_tpwT_cmT_N550/` | `fortran/build/case1/` (suffixes only for non-default switches) |
| build defaults `SHT=blas`, `OMP=0`, `NORDER=512` | `SHT=otf`, `OMP=8`, `NORDER=511` |
| `C2021=1` (default constants), `C2021=0` | removed (the default); `LEGACY_CONSTANTS=1` |
| `sl_model_mod_pureplus.f90` | `sl_model_mod.f90` |
| `exact_response_memory_mod_stage2.f90`, module `exact_response_memory_mod` | `modal_state_mod.f90`, module `modal_state_mod` |
| `pureplus_io_mod.f90`, module `pureplus_io_mod` | `io_mod.f90`, module `io_mod` |
| `support_runtime_real512_orig.f90` | `support_mod.f90` |
| `PUREPLUS3_OPT` | `SLEVOLVE_OPT` |
| `PUREPLUS_INNER`, `_INNER_LOG`, `_AA_M`, `_AA_MASK`, `_REUSE`, `_OUTER`, `_ICECACHE`, `_SNAP`, `_SNAP_FIELDS`, `_HORIZ`, `_FMT`, `_TEXT`, `_BIN`, `_BININ`, `_CKPT`, `_BATCH`, `_FFT`, `_ZERO`, `_COMPAT`, `_SHTPROF` | `SLEVOLVE_INNER`, ... (same suffix) |
| Fortran `dSlm`, `dsllm`, `deltasllm` (Python `dSlm`, ...) | `load_inc_lm`, `sl_inc_lm`, `sl_tot_lm` |
| Fortran prefixes `pp3_` (transform reductions), `pp4_` (inner iteration), `pp4_dtm*/dsm*/c1/c2` (extrapolation), `pp4_ru_`, `pp4_reuse_` (warm start), `pp4_outer`, `pp4_secant` | `opt_`, `inner_`, `extrap_`, `warm_`, `outer_mode`, `outer_secant` |
| `examples/quickstart_disc.py` | `examples/quickstart.py` (synthetic Earth) |

## Checks made during development

These checks were made during development with tools that are not part of this release; they are
recorded here and are not meant to be repeated.

**Exactness against the legacy solver** (code of the paper, Zenodo v1.0.0, before the floating-ice
correction of 2026-10-26)

* Run both solvers on the same inputs and configuration with `SLEVOLVE_INNER=off` and no warm start.
  The legacy reference build uses the corrected `fft99b`. The legacy build scripts and the comparison
  scripts used during development are not part of this release.
* **Pass criteria:**
  * final topography differs by **< 1e-11 m** between exact SLEvolve variants and the legacy solver
    (binary output);
  * identical `numiter` at every step;
  * identical `ocean` and `beta` masks.
* Comparisons against legacy *text* output are limited by the ES16.9E2 print precision (~4e-6 m).

**Acceptance criteria used during development**

* Every stage was accepted only after it reproduced the preceding stage numerically on the benchmark
  cases (paper Sect. 2.8): same inner-iteration counts and masks, and topography within round-off.
* Switches that only change the iteration start (extrap, warm, T0 guess) were accepted on agreement
  within `epsilon1` and on the converged outer misfit.
* Fortran and Python are cross-checked against each other.
* Repeat runs give identical output.

**File converters**


* Converters, on one file of each text layout (`row-e.18`, `row-e16.7`, `row-E18.10`):
  * Fortran- and Python-written `.npy` files are identical.
  * `np.load` of the Fortran file equals the text parse.
  * text -> npy -> text reproduces the original file exactly in its own format, in both languages.
  * `es16.9` output is identical between Fortran and Python and round-trips exactly.
  * The Fortran reader also handles fortran_order=True and `|i1` files.
* Solver, 20 Case 1 steps with all snapshot fields: reading `<f4` `.npy` ice and a `<f8` `.npy`
  topography gives output identical to reading the headerless `.bin` caches, and all scalar text
  outputs are identical (except CPU time).
* Python, 20 steps, lmax 511: reading the Fortran-written `.npy` inputs gives output identical to
  reading the text files; the `.npy` headers of `tgrid`, `ocean` and `dS_converged` are identical to
  the Fortran files'.

**Next:** [Lineage and references](lineage.md).
