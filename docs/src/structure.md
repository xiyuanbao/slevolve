# Code structure

This page describes how the release sources are organised and how one time step is carried out in
the code. Routine and class names are given with the file that contains them. The code excerpts are
copied from the release sources. `...` marks lines left out between complete statements.

## 1. Repository layout

```
README.md                          overview, build, quick start, references
LICENSE                            GPL-3.0 text
LICENSE-upstream-MIT               MIT license of the upstream Fortran code
CITATION.cff                       citation metadata
fortran/
  build.sh                         release build: generates user_specs from the support module (sed),
                                   compiles src/ with the chosen DRIVER/SHT/FFT/OMP, builds the tools
  IO_FORMATS.md                    binary file formats (.npy, .bin)
  src/
    run_from_zero_plus.f90         driver: forward run, ./exe LAST (step = DT or times file)
    run_outer_plus.f90             driver: one outer pass per process, ./exe LAST [ITERSL]
    run_outer_mem.f90              driver: outer passes K0..K1 in one process, ./exe LAST K0 K1
    sl_model_mod.f90               module sl_model_mod: the solver (step 0, time step, inner
                                   iteration, outer update, warm-start store, outputs)
    modal_state_mod.f90            module modal_state_mod: recursive modal state
                                   (advance_and_read, read_horizontal, commit, checkpoint)
    io_mod.f90                     module io_mod: run-time options (SLEVOLVE_* environment
                                   variables), grid/spectral read and write, snapshot selection
    npy_io_mod.f90                 module npy_io_mod: dependency-free .npy reader/writer
    spharmt_otf.f90                module spharmt: spherical-harmonic transforms (BLAS table,
                                   on-the-fly recurrence, verification backends), FFTW/FFT991
    otf_kernels.inc                on-the-fly Legendre recurrence kernels, included by spharmt_otf.f90
                                   for SHT=otf
    spharmt_plus.f90               alternative module spharmt without the on-the-fly backend; not
                                   compiled by build.sh
    support_mod.f90                unedited legacy support modules (constants_mod, user_specs_mod,
                                   planets_mod, sl_io_mod); build.sh writes the configured copy
  tools/
    make_cache.f90                 text grids -> flat little-endian fp64 .bin caches
    text2npy.f90, npy2text.f90     text <-> .npy grid converters (link npy_io_mod)
    build_npytools.sh              builds text2npy/npy2text on their own
pyproject.toml                     pip metadata of the package slevolve (setuptools; package dir python/)
python/
  run_ais.py, run_outer.py,        script entry points (python python/run_ais.py ...); each calls main()
  init_topo_segmented.py           of the package module of the same name
  cache/                           ice caches of a checkout (created on first use)
  slevolve/                        the package (relative imports; pip install . installs it)
    __init__.py                    __version__ = "1.0.0"; main classes, imported on first use
    _paths.py                      default root and ice-cache directory (checkout or installed)
    run_ais.py                     forward-run driver (CLI, console script slevolve-forward); writes
                                   .npy outputs and run_log.json
    run_outer.py                   ice age run driver (slevolve-outer): outer (initial-topography)
                                   iteration, Picard or per-cell secant, optional warm start across passes
    init_topo_segmented.py         partially relaxed (segmented) first-pass initial topography
                                   (slevolve-init-topo)
    solver.py                      reference solver: Planet, LoveNumbers, ExactMemory, SeaLevelSolver
    solver_opt.py                  optimised step: OptSeaLevelSolver, ExactMemoryFast, parse_opts,
                                   Anderson mixing (_Anderson)
    kernels.py                     Numba kernels for the hereditary advance and grid/spectral algebra
    sht.py                         NumPy/BLAS transform (class SHT) in the conventions of the Fortran
    sht_ducc.py                    ducc0 transform (class DuccSHT) with the same interface; make_sht
    npyio.py                       shared file conventions: read_grid, save_grid, save_spectral, IceDir
    horiz.py                       grid components of the horizontal displacement (gradient)
    tools/                         subpackage slevolve.tools
      text2npy.py, npy2text.py     text <-> .npy grid converters (slevolve-text2npy, slevolve-npy2text)
      npy2nc.py                    per-step .npy snapshots -> one NetCDF4 file (slevolve-npy2nc)
      horiz_grad.py                grid components from horiz<N>.npy (with --selftest)
      horiz_history_check.py       verification: recursive horizontal read-out vs full-history sum
      write_lovenum.py             writes Love numbers in the solver's text layout
      make_reference.py            reference/ files from an output directory
      plot_fields.py               Mollweide and south polar stereographic maps (coast = 0 m contour)
  tools/                           script entry points for slevolve.tools (python python/tools/<name>.py ...)
tools/
  build_case1_total_ice.py         Case 1 total (grounded plus floating) ice forcing
  make_grid.py                     grid files GLlat_512.txt, GLlon_512.txt
  make_synthetic_case.py           synthetic-Earth input tree (quick start, Fortran tests)
  make_synthetic_case.f90          Fortran twin of make_synthetic_case.py (identical .npy files)
  plot_reference.py                maps of the reference results and of the synthetic test
  make_logo.py                     logo files in docs/src/img
  regrid.py                        nearest-source sampling onto the solver grid (NearestSourceSampler)
  prepare_case2_inputs.sh          Case 2 input tree from the archived files
  prepare_case3_nadi.py            Case 3 (ICE-NADI-1) input tree
lovenum/
  prem_512.l60K2C.sum18p6.dum19p2.tz19p4.lm22  Love numbers, Case 1 Earth model
  prem_512.l120C.ump5.lm5          Love numbers, Case 2 Earth model
  prem_512.l96C.ump3.lm1           Love numbers, Case 3 Earth model
  DEGREE1_NOTES.md, DEGREE1_NOTES_zh.md
                                   notes on degree-1 Love numbers (English, Chinese)
tests/
  run_all.py                       runs the nine Python tests below
  test_recursion.py, test_sht_roundtrip.py, test_small_case.py,
  test_extrap_formula.py, test_small_case_vardt.py, test_floating_ice.py,
  test_ramp.py                     (second-order placement and analytic rates)
  test_j2dot_two_routes.py         J2-dot from the geoid rate and from m3-dot vs a closed form
                                   (axial rotation option off, deg2 and full)
  test_tpw_direction.py            TPW direction from TPWrate vs the load and the geoid rate
  test_fortran_python.sh           Fortran test: Fortran vs Python on the synthetic Earth
  test_fortran_python_vardt.sh     Fortran test with a variable-step times file (Fortran vs Python)
  test_fortran_python_ramp.sh      Fortran test with RAMP=1 and SLEVOLVE_RATES=1 (Fortran vs Python)
  test_fortran_python_axial.sh     Fortran test with AXIAL_ROTATION=deg2|full and the rates (Fortran vs Python)
  test_fortran_synthetic.sh        Fortran-only test: Fortran vs reference/synthetic/tgrid8.npy
  npy_maxdiff.f90                  max |A - B| of two .npy grids (used by the Fortran-only test)
  check_reference.py               compares a benchmark run with reference/
examples/
  quickstart.py                    quick start (Python solver, synthetic Earth, 128 x 256)
  case1_forward.sh                 Case 1 forward run (fortran|python)
  case2_one_pass.sh                Case 2, first outer pass (fortran|python)
  case3_ice_age.sh                 Case 3 ice age run (fortran|python), segmented first-pass topography
  run_outer_mem.sh                 general wrapper for the in-process outer-iteration driver
docs/
  mkdocs.yml, src/*.md             this manual
```

**Python imports.** The modules of `python/slevolve/` import each other relatively
(`from .solver import ...`, `from .. import npyio` in `slevolve.tools`), so the same files work
installed (`pip install .`) and from a checkout. The scripts `python/run_ais.py`, `python/run_outer.py`,
`python/init_topo_segmented.py` and `python/tools/<name>.py` put `python/` on `sys.path` and call `main()`
of the package module, so the documented commands are unchanged. Code outside the package (tests,
examples, `tools/`) imports from the package: `from slevolve.solver import LoveNumbers`,
`from slevolve.tools.plot_fields import mollweide_map`, after `sys.path.insert(0, "python")` in a checkout
or with no path change when the package is installed. `run_ais` and `run_outer` parse their arguments
and pin the thread counts at import, before NumPy is imported; `import slevolve` itself does not import
NumPy.

## 2. Fortran architecture

### 2.1 Modules and build order

`build.sh` compiles, in this order: the generated support module, `npy_io_mod.f90`,
`io_mod.f90`, `spharmt_otf.f90` (as `spharmt.o`), `modal_state_mod.f90`,
`sl_model_mod.f90` and the selected driver. The `use` relations are:

```
driver (run_from_zero_plus | run_outer_plus | run_outer_mem)
 |-- sl_model_mod                      sl_model_mod.f90
 |    |-- spharmt                      spharmt_otf.f90 (+ otf_kernels.inc with SHT=otf)
 |    |-- user_specs_mod               generated support module
 |    |-- planets_mod                  generated support module
 |    |    `-- constants_mod
 |    |-- sl_io_mod                    generated support module (legacy text read_sl / write_sl)
 |    |    `-- user_specs_mod
 |    |-- modal_state_mod              modal_state_mod.f90
 |    `-- io_mod                       io_mod.f90
 |         `-- npy_io_mod              npy_io_mod.f90
 |-- spharmt                           (sht_profile counters)
 |-- user_specs_mod                    (dt1; outer drivers)
 `-- io_mod                            (pp_read_options, pp_shtprof)
```

Compile-time parameters (step length `dt1`, number of steps `L_sim`, truncation `norder`, tolerance
`epsilon1`, file names, `tpw`, `checkmarine`, `initial_topo`, constants) live in `user_specs_mod`
and are set by `build.sh`. Run-time switches are environment variables read by
`io_mod::pp_read_options` (with `inner_read_options` and `opt_read_options`) and by
`spharmt::sht_env_init` (`SLEVOLVE_FFT`); see [Command-line options](options.md).

### 2.2 Drivers

All three drivers call the same solver routines:

1. `sl_solver_checkpoint(itersl, dtime)`: checks `dtime == dt1`, sets the planet (`set_planet`),
   reads the run-time options, initialises the transform backend (`sht_env_init`, `spharmt_init`).
2. Step 0: `sl_timewindow(0)`, `sl_allocate_and_initialize_array`, `sl_solver_init(itersl, 0.0)`,
   `sl_deallocate_array`.
3. Steps `iter = 1..LAST`: `sl_timewindow(iter)`, `sl_allocate_and_initialize_array`,
   `sl_solver(itersl, iter, dtime, 0.0)`, `sl_deallocate_array`.
4. `sl_print_counts`.

| driver | command | what differs |
|---|---|---|
| `run_from_zero_plus.f90` | `./exe LAST` | forward run; `itersl = 1`, step length `dt1` = build-time `DT` (or the times file); optional transform profile print (`SLEVOLVE_SHTPROF`) |
| `run_outer_plus.f90` | `./exe LAST [ITERSL]` | one outer pass; `dtime = dt1`; pass `ITERSL > 1` reads `pred_pres_topo_<k-1>` and `tgrid0_<k-1>` from `OUTPUT_SLM/` in `sl_solver_init` |
| `run_outer_mem.f90` | `./exe LAST K0 K1` | passes `K0..K1` in one process; sets `pp_mem_outer`, reads `SLEVOLVE_ICECACHE`; calls `sl_outer_reset` between passes; copies `OUTPUT_SLM/numiter` to `numiter_<k>` and prints the wall time per pass |

The outer update of the initial topography is done in `sl_solver_init` at the start of pass
`itersl > 1` (Picard update, or `outer_secant` with `SLEVOLVE_OUTER=secant`). At the last step,
`sl_solver` writes `pred_pres_topo_<itersl>` and `tgrid0_<itersl>` (and with `run_outer_mem`
also keeps them in memory for the next pass).

### 2.3 One time step (`sl_model_mod::sl_solver`)

1. **Setup and ice.** `opt_setup`; the ice grid of the current step is read with `pp_ice_get`
   (at `iter == 1` also the time-zero grid). `pp_ice_get` serves the grid from the in-memory ice
   cache (`run_outer_mem`) or through `pp_read_grid` (`.npy`, then `.bin`, then text). The physical
   times of the previous and current step are appended to the time array file.
2. **Step-0 fields.** `tinit_0`, `cxy0`, `beta0` stay resident (re-read only with
   `SLEVOLVE_COMPAT` bit 1); the previous step's ocean function is already in `cxy`.
3. **Earth model (once per run).** The Love-number file is parsed in the block "Read in Love
   numbers", guarded by `pp_love_ready`; the response spectrum `exact_spec` (decay rates and the
   weights `a_sl`, `a_rr`, degree-2 `a_il`, `a_m`, `a_tsl`, `a_trr`, and `a_hl`, `a_thl` when
   `SLEVOLVE_HORIZ > 0`) is built from it, and `exact_state%init` allocates the modal state.
4. **Advance-and-read.** At `iter == 1` the state is zeroed (`exact_state%zero`); otherwise it is
   resident (or reloaded with `load_checkpoint` under `SLEVOLVE_COMPAT` bit 4). Then
   `exact_state%advance_and_read` advances every modal state by one decay factor and returns the
   viscous sea-level and radial sums and the degree-2 rotational sums. These are frozen for the whole
   inner iteration. The cumulative loads `deltaS(:,:,exact_prev)`, `oldil`, `oldm` and the previous
   `deltalambda` are restored from `p_ocean`, `p_il`, `p_m`, `p_lambda`.
5. **Ice and ocean setup.** Transform of the initial topography (once, `pp_t0lm_ready`), the
   grounded-ice and beta fields (STEP 3), the warm-start read `warm_read` (with
   `SLEVOLVE_REUSE=warm`; replaces the ocean mask), and the batched analyses (`opt_ana_multi`) of the
   ice, `cstar` and `tO` grids.
6. **Inner start.** Eustatic first guess for the ocean-load increment `load_inc_lm`; with
   `SLEVOLVE_INNER=extrap` (or `anderson`) and two stored converged increments, the extrapolated
   start replaces it; with a warm-start record, the stored increment of the previous pass replaces
   both.
7. **Inner iteration** (the `do` loop between "Start of inner loop" and "End of inner loop"):
   rotation (`tpw`: moment of inertia, `m`, `lambda`, using the frozen `exact_sum_il`, `exact_sum_m`,
   `exact_tsl`, `exact_trr`); sea-level change from the load (STEP 6, elastic part plus frozen
   `viscous`); synthesis; `r0`, conservation term and new ocean load (STEPS 7 to 9); convergence
   measure `xi` on `sum(abs(load_inc_lm))`; ocean-function update and re-analysis while `xi > epsilon1`;
   optional Anderson mixing (`inner_anderson`). The loop exits when `xi <= epsilon1` and aborts at
   9999 iterations.
8. **After convergence.** The converged `load_inc_lm` and `dtime` are pushed into the extrapolation
   history (`extrap_dsm1`, `extrap_dsm2`, `extrap_dtm1`, `extrap_dtm2`); `R` and `G` are formed from the frozen
   `exact_radial` when they are written (`lazyrg`); the optional horizontal read-out
   (`exact_state%read_horizontal`, `SLEVOLVE_HORIZ`) is done on the same states.
9. **Commit.** The current load `rhoi*dicestar + rhow*load_inc_lm` (degrees 1..norder) is added to
   `p_load` with `exact_state%commit` (with `dil`, `dm`, `dlambda` when `tpw`); the cumulative
   `p_ice`, `p_ocean`, `prev_ice` are stored; `save_checkpoint` only with `SLEVOLVE_CKPT` or
   `SLEVOLVE_COMPAT` bit 4.
10. **Topography and outputs.** `topoxy = tinit - deltaslxy`, new ocean function `cxy` (STEP 11);
    `warm_write` (warm-start store); diagnostics (ocean areas, GMSL, ice volume); text scalars
    (`numiter`, `nmelt`, `TPW<N>`, ...); grids through `pp_write_grid` (`tgrid`, `ocean`, `beta`,
    `R`, `G`) and `pp_write_cplx` (`dS_converged`), subject to `SLEVOLVE_SNAP`; at the last step
    `pred_pres_topo_<itersl>` and `tgrid0_<itersl>`.

## 3. Python architecture

```
run_ais.py (forward)        run_outer.py (outer passes)       init_topo_segmented.py
      |                            |                                 |
      |  npyio (read_grid, IceDir, save_grid, save_spectral)        |  solver.LoveNumbers, Planet
      |  solver.LoveNumbers, Planet                                  |
      |  sht_ducc.make_sht -> sht.SHT | sht_ducc.DuccSHT            |
      v                            v
solver_opt.OptSeaLevelSolver  (subclass of solver.SeaLevelSolver, overrides step)
      |-- state: solver_opt.ExactMemoryFast (subclass of solver.ExactMemory)  [opt "numba"]
      |          or solver.ExactMemory
      |-- kernels (Numba kernels)                                            [opt "numba"]
      |-- solver_opt._Anderson                                               [--inner anderson]
      `-- horiz.gradient (grid components, called by the drivers)            [--horizontal-grid]
```

* `run_ais.py` parses the options, builds the ice cache, reads the topography (`read_grid_any`) and
  the Love numbers (`solver.LoveNumbers`), creates the transform with `sht_ducc.make_sht`
  (`--sht numpy|ducc`), and uses `OptSeaLevelSolver` when `--opt` or `--inner` is given, otherwise
  the reference `SeaLevelSolver` (with `dt = --dt`, or the first step of `--times`). It then calls `step(ice[k])` for `k = 1..steps`
  and writes `tgrid<k>.npy`, `ocean`, `beta`, `R`, `G`, `dS_converged`, `horiz` (`--save-fields`,
  `--horizontal`, `--horizontal-grid`), `numiter.npy` and `run_log.json`.
* `run_outer.py` creates one `OptSeaLevelSolver` per outer pass (`dt=--dt`, default `opt
  raw,diag,dsl,numba,sums`, ducc0 transform), passes the dict `store` that survives across passes
  for `--reuse warm`, runs all steps, and updates the initial topography by Picard or per-cell
  secant (`--outer`). It writes `pred_pres_topo_<it>.npy`, `tgrid0_<it>.npy`, saved steps and
  `T0_next.npy`.
* `solver.SeaLevelSolver.__init__` does what `sl_solver_init` does (step 0: initial topography,
  ocean and beta functions, their spectra) and creates the state `ExactMemory`;
  `OptSeaLevelSolver.__init__` replaces it by `ExactMemoryFast` when `numba` is selected.
* `OptSeaLevelSolver.step` follows the same sequence as `sl_solver`: `state.advance_and_read()`,
  grounded ice I* on the previous converged topography, beta and ocean grids, ice analysis and increments, eustatic guess,
  extrapolated start (`inner != "off"`), warm start (`reuse_store`), inner loop (rotation,
  `kernels.b18`, synthesis, conservation, `kernels.dS_update`, convergence measure `xi`,
  ocean-function update or Anderson mixing), history push, `R`/`G` from `radial`, optional
  `horizontal_lm`, `state.commit(...)`, topography and ocean function (`kernels.topo_ocean`),
  diagnostics. It returns a dict of fields (`tgrid`, `ocean`, `beta`, `R`, `G`, `dS_converged`,
  `horiz`, scalars).
* `parse_opts` turns `--opt` (`raw`, `diag`, `dsl`, `numba`, `sums`, or `all`) into a set; each item
  switches one optimisation of `step`.
* `sht.SHT` (NumPy/BLAS table) and `sht_ducc.DuccSHT` (ducc0) share `analysis(grid, out)` and
  `synthesis(coeffs, out)`; `make_sht(backend, ...)` selects one. The `raw` option exists only for
  ducc0.
* `npyio` holds the file conventions shared with the Fortran (`read_grid` tries `.npy`, `.bin`,
  text; `IceDir` reads `<dir>/<pfx><k>` on demand).
* `horiz.gradient(U, nlat, nlon)` turns the spectral `Psi_lm` into `(u_theta, u_phi)` with the
  ducc0 spin-1 gradient synthesis.

## 4. Code excerpts

### 4.1 Modal-state advance and commit

The advance loop of `state_advance_and_read` in `fortran/src/modal_state_mod.f90`
(the degree-2 rotational states follow after it):

```fortran
call update_decay(self,spec,dt_years)
viscous=(0.0_dp,0.0_dp); radial=(0.0_dp,0.0_dp)

! Packed layout is [degree][pole][m=0:l], so each pole update is one
! contiguous vector operation.  R itself is the hereditary response state:
! R_new = d*R_old + (1-d)*DL_old   (DL = cumulative load, p_load).
! degrees are independent (each writes a disjoint slice of
! r_load / viscous / radial and accumulates over kk in the same order),
! so this loop threads identically.
!$omp parallel do schedule(dynamic,4) default(shared) private(ll,kk,ibase,rstart,rend,d,omd)
do ll=1,self%lmax
   ibase=lm_index(ll,0)
   do kk=1,self%nmod(ll)
      d=self%decay(kk,ll); omd=1.0_dp-d
      rstart=int(self%roffset(ll))+(kk-1)*(ll+1)
      rend=rstart+ll
      self%r_load(rstart:rend)=d*self%r_load(rstart:rend) + &
                                omd*self%p_load(ibase:ibase+ll)
      viscous(ibase:ibase+ll)=viscous(ibase:ibase+ll) + &
                              spec%a_sl(kk,ll)*self%r_load(rstart:rend)
      radial(ibase:ibase+ll)=radial(ibase:ibase+ll) + &
                             spec%a_rr(kk,ll)*self%r_load(rstart:rend)
   end do
end do
!$omp end parallel do
```

`state_commit` in the same file only adds the converged increments to the accumulated loads:

```fortran
subroutine state_commit(self,current_load,dil,dm,dlambda)
  class(response_memory_state), intent(inout) :: self
  complex(dp), intent(in) :: current_load(:)
  real(dp), intent(in), optional :: dil(3,3),dm(3)
  complex(dp), intent(in), optional :: dlambda(0:2)
  self%p_load=self%p_load+current_load
  if (present(dil)) self%p_il=self%p_il+dil
  if (present(dm)) self%p_m=self%p_m+dm
  if (present(dlambda)) self%p_lambda=self%p_lambda+dlambda
end subroutine state_commit
```

The Python equivalent of the advance, the body of `advance_and_read` in `python/slevolve/kernels.py` (state stored as
`R[l, k, m]`, used by `ExactMemoryFast`):

```python
nl = R.shape[0]
for l in prange(nl):
    nk = nmod[l]
    for m in range(nl):
        visc[l, m] = 0.0
        rad[l, m] = 0.0
    for k in range(nk):
        d = decay[k, l]
        o = omd[k, l]
        s = a_sl[k, l]
        q = a_rr[k, l]
        for m in range(l + 1):
            x = R[l, k, m] * d + P[l, m] * o
            R[l, k, m] = x
            visc[l, m] += x * s
            rad[l, m] += x * q
```

Here `P` (`p_load` in both codes) holds the cumulative load ΔL(t_{j-1}).

and `ExactMemory.commit` in `python/slevolve/solver.py`:

```python
def commit(self, current_load, dil=None, dm=None, dlambda=None):
    self.p_load += current_load
    if dil is not None:
        self.p_il += dil
    if dm is not None:
        self.p_m += dm
    if dlambda is not None:
        self.p_lambda += dlambda
```

(shown without the branch of the optional second-order placement: with `ramp` on, the commit first
adds \((1-\varphi_k)\) times the increment to every state, `ExactMemory._ramp_commit_load` and
`kernels.py::ramp_commit` in Python, `state_ramp_commit` in Fortran; see
[Method](method.md#second-order-load-placement-optional)).

### 4.2 Placement of advance and commit in the time step

In `sl_solver` (`fortran/src/sl_model_mod.f90`) the state is advanced and read once before
the inner iteration and committed once after convergence; the inner iteration and the R, G and
horizontal read-outs lie between the two parts (left out):

```fortran
! the 53 MB hereditary state stays in memory across timesteps;
! the baseline writes and reads it back on every single step.
if (iter == 1) then
   call exact_state%zero()
   pp_state_ready = .true.
else if (.not. pp_state_ready .or. iand(pp_compat,4) /= 0) then
   call exact_state%load_checkpoint(trim(exact_statefile),exact_spec)
   pp_state_ready = .true.
endif
call exact_state%advance_and_read(exact_spec,real(dtime,dp),exact_viscous,exact_radial, &
                                  exact_sum_il,exact_sum_m,exact_tsl,exact_trr)
...
! Commit the converged current-step increments only now.  The current
! load has zero age at this instant, so it did not contribute to the
! hereditary response used above; it becomes memory for future steps.
exact_current_load = (0.0_dp,0.0_dp)
!$omp parallel do schedule(static) default(shared) private(l,m,exact_idx)
do l = 1,norder
   do m = 0,l
      exact_idx = lm_index(l,m)
      exact_current_load(exact_idx) = rhoi*dicestar(l,m,exact_cur) + rhow*load_inc_lm(l,m)
   enddo
enddo
!$omp end parallel do
if (tpw) then
   call exact_state%commit(exact_current_load,dil(:,:,exact_cur),dm(:,exact_cur),dlambda(2,0:2,exact_cur))
else
   call exact_state%commit(exact_current_load)
endif
```

The same placement in `OptSeaLevelSolver.step` (`python/slevolve/solver_opt.py`):

```python
# ---- hereditary advance -----------------------------------------
t0 = pc()
viscous, radial, sum_il, sum_m, exact_tsl, exact_trr = self.state.advance_and_read()
mark("memory", t0)
...
# -------------------------------------------------------- commit ----
t0 = pc()
if use_nb:
    K.current_load(dicestar, self.load_inc_lm, p.rhoi, p.rhow, self.current_load)
else:
    np.multiply(dicestar, p.rhoi, out=self.current_load)
    np.multiply(self.load_inc_lm, p.rhow, out=self._sp2)
    self.current_load += self._sp2
    self.current_load[0, :] = 0.0
if self.tpw:
    self.state.commit(self.current_load, dil_cur, dm_cur, dlam_cur2)
else:
    self.state.commit(self.current_load)
```

### 4.3 Extrapolated inner start

The coefficients are computed by `extrap_coeffs` in `fortran/src/modal_state_mod.f90`; for equal
steps they are exactly 2 and -1:

```fortran
! Coefficients of the extrapolated inner start (variable steps):
!   dS0_j = dt_j [ Sdot_{j-1} + (Sdot_{j-1} - Sdot_{j-2}) (dt_j + dt_{j-1}) / (dt_{j-1} + dt_{j-2}) ],
!   Sdot_i = dS_i / dt_i,   i.e.  dS0_j = c1 dS_{j-1} + c2 dS_{j-2};
! equal steps give exactly c1 = 2, c2 = -1.  dt0 = dt_j, dtm1 = dt_{j-1}, dtm2 = dt_{j-2}.
pure subroutine extrap_coeffs(dt0, dtm1, dtm2, c1, c2)
  real(dp), intent(in) :: dt0, dtm1, dtm2
  real(dp), intent(out) :: c1, c2
  if (dt0 == dtm1 .and. dtm1 == dtm2) then
     c1 = 2.0_dp; c2 = -1.0_dp
  else
     c1 = dt0 / dtm1 * (1.0_dp + (dt0 + dtm1) / (dtm1 + dtm2))
     c2 = -dt0 / dtm2 * ((dt0 + dtm1) / (dtm1 + dtm2))
  endif
end subroutine extrap_coeffs
```

and applied in `sl_solver` (`fortran/src/sl_model_mod.f90`), where `dtr` is the length of the
current step (from the times file, else the build-time `DT`):

```fortran
! extrap: linear-in-time extrapolation of the converged
! increment; the cumulative current load that the first rotation
! evaluation reads is made consistent with it.  For uniform steps this
! is 2 dS(n-1) - dS(n-2); for variable steps the rate dS/dt is
! extrapolated linearly between step midpoints and scaled by dt(n).
if (inner_mode > 0) then
   if (.not. allocated(extrap_dsm1)) then
      allocate(extrap_dsm1(0:norder,0:norder), extrap_dsm2(0:norder,0:norder))
      extrap_dsm1 = (0.0,0.0); extrap_dsm2 = (0.0,0.0)
   endif
   if (inner_nhist >= 2) then
      ! paper Sect. 2.5: dS0_j = c1 dS_{j-1} + c2 dS_{j-2} (modal_state_mod::extrap_coeffs)
      call extrap_coeffs(dtr, extrap_dtm1, extrap_dtm2, extrap_c1, extrap_c2)
      !$omp parallel do schedule(static) default(shared) private(m)
      do m = 0,norder
         load_inc_lm(:,m) = extrap_c1 * extrap_dsm1(:,m) + extrap_c2 * extrap_dsm2(:,m)
         deltaS(:,m,exact_cur) = deltaS(:,m,exact_prev) + load_inc_lm(:,m)
      enddo
      !$omp end parallel do
   endif
   inner_nh = 0; inner_lastmixed = .false.; inner_plain_rest = .false.; inner_nrefreeze = 0
```

The Python version is `extrap_guess` in `python/slevolve/solver_opt.py`:

```python
def extrap_guess(h1, h2, dt0, dt1, dt2, out=None):
    """Extrapolated inner start of step j (paper Sect. 2.5):

        dS0_j = dt_j [ Sdot_{j-1} + (Sdot_{j-1} - Sdot_{j-2}) (dt_j + dt_{j-1}) / (dt_{j-1} + dt_{j-2}) ],
        Sdot_i = dS_i / dt_i,

    with h1 = dS_{j-1}, h2 = dS_{j-2}, dt0 = dt_j, dt1 = dt_{j-1}, dt2 = dt_{j-2}.  Equal steps give
    exactly 2 dS_{j-1} - dS_{j-2}."""
    if out is None:
        out = np.empty_like(h1)
    if dt0 == dt1 == dt2:
        np.multiply(h1, 2.0, out=out)
        out -= h2
    else:
        w = (dt0 + dt1) / (dt1 + dt2)
        np.multiply(h1, dt0 / dt1 * (1.0 + w), out=out)
        out -= (dt0 / dt2 * w) * h2
    return out
```

called in `OptSeaLevelSolver.step` with the lengths of the current and the two previous steps:

```python
if acc and self._dS_nhist >= 2:
    h1, h2 = self._dS_hist  # dS_{n-1}, dS_{n-2}
    # linear-in-time extrapolation of the increment: for uniform steps
    # 2 dS_{n-1} - dS_{n-2}; for variable steps the rate dS/dt is
    # extrapolated between step midpoints and scaled by dt_n
    extrap_guess(h1, h2, float(self.dt), self._dt_hist[0], self._dt_hist[1],
                 out=self.load_inc_lm)
    # the rotation block of iteration 1 reads deltaS_cur: give it the
    # same guess (plain iteration 1 reads 0 here)
    np.add(deltaS_prev, self.load_inc_lm, out=self.deltaS_cur)
```

### 4.4 Love-number file read

The reader in `sl_solver` (`fortran/src/sl_model_mod.f90`), executed once per run: one
header line, then per degree the number of modes, the decay rates, the elastic and asymptotic
lines, and the h, l, k residues, followed for l >= 2 by the same lines for the tidal numbers:

```fortran
open(unit = 201, file = trim(planetfolder)//trim(planetmodel), status = 'old')
! Reader as in the legacy code of Han et al. (2022)
read(201,*)
do j = 1,norder
   read(201,*) legord(j), nmodes(j)
   nm = nmodes(j)
   xn = real(legord(j))
   ll = legord(j)
   nmod(ll) = nm
   read(201,*) (s(i,ll), i=1, nm)
   read(201,*) elast(1,ll), elast(2,ll), taurr, taurt, elast(3,ll)
   read(201,*) asymv(1,ll), asymv(2,ll), taurr, taurt, asymv(3,ll)
   read(201,*) (resh(i,ll), i=1,nm)
   read(201,*) (resl(i,ll), i=1,nm)
   read(201,*) (resk(i,ll), i=1,nm)
   if (xn .lt. 1.5) cycle
   read(201,*) telast(1,ll), telast(2,ll), dmx, dmy, telast(3,ll)
   read(201,*) tasymv(1,ll), tasymv(2,ll), dmx, dmy, tasymv(3,ll)
   read(201,*) (tresh(i,ll), i=1, nm)
   read(201,*) (tresl(i,ll), i=1, nm)
   read(201,*) (tresk(i,ll), i=1, nm)
enddo
close(201)
```

`solver.LoveNumbers` in `python/slevolve/solver.py` reads the same layout.

### 4.5 Horizontal-displacement read-out

In `sl_solver` (`fortran/src/sl_model_mod.f90`), between the inner iteration and the
commit: the modal states are read with the `l` weights and combined with the elastic `l` number and
the degree-2 rotational term:

```fortran
call exact_state%read_horizontal(exact_spec, exact_horiz, exact_thl)
hzlm(:,:) = (0.0,0.0)
if (nmelt /= 0) then
   !$omp parallel do schedule(static) default(shared) private(l,m,ttl,exact_idx)
   do l = 1,norder
      ttl = (4.0 * pi * (radius**3)) / (mass * (2.0 * real(l) + 1.0)) ! (eq. B16)
      do m = 0,l
         exact_idx = lm_index(l,m)
         hzlm(l,m) = ttl * (elast(2,l)/real(l)) * (rhoi * deltaicestar(l,m,exact_cur) &
                  & + rhow * deltaS(l,m,exact_prev) + rhow * load_inc_lm(l,m)) &
                  & + ttl * exact_horiz(exact_idx)
      enddo
   enddo
   !$omp end parallel do
   if (tpw .and. norder >= 2) then
      do m = 0,2
         hzlm(2,m) = hzlm(2,m) + ((telast(2,2)/2.0) * (deltalambda(2,m,exact_prev) &
                  & + dlambda(2,m,exact_cur)) + exact_thl(m)) / gacc
      enddo
   endif
endif
```

The Python version, from `SeaLevelSolver.horizontal_lm` in `python/slevolve/solver.py`:

```python
visc_h, thl = st.read_horizontal()
A = p.rhoi * deltaicestar + p.rhow * deltaS_prev
A += p.rhow * load_inc_lm
hz = self._ttlle * A + self._ttl * visc_h
hz[0, :] = 0.0
sg = getattr(self, "_msign", 1.0)
if self.tpw and rot_tot is not None and self.norder >= 2:
    rot = (self.love.lTE[2] * rot_tot + thl) / p.gacc
    rot[1] *= sg
    hz[2, 0:3] += rot
if sg != 1.0:  # raw (sign-free) backend: V = (-1)^m U -> U
    hz[:, 1::2] *= -1.0
return hz
```

## 5. Where to change what

| modification | where |
|---|---|
| New Earth model | Write a Love-number file in the layout of [Love numbers](lovenum.md) (`python/tools/write_lovenum.py`) into `INPUT_FILES/earthmodel/`. Fortran: `PLANETMODEL=<file>` in `fortran/build.sh` (compiled into `user_specs_mod`). Python: `--earth <file>` (`run_ais.py`, `run_outer.py`) |
| New ice history | Grids `<prefix><k>` (`.npy`, `.bin` or text), k = 0..N, in `INPUT_FILES/icemodel/`. Fortran: `ICEMODEL=<prefix>`, `NSTEPS`, `DT` in `build.sh`. Python: `--icedir`, `--icepfx`, `--steps`, `--dt` (`run_ais.py`); `--inputs`, `--icepfx`, `--nfiles`, `--dt` (`run_outer.py`). Input trees for the cases: `tools/` |
| Spectral truncation | Fortran: `NORDER` in `build.sh` (`norder` in `user_specs_mod`). Python: `--norder`/`--lmax` (`run_ais.py`), `--lmax` (`run_outer.py`) |
| Grid size | Fortran: `nglv` is fixed at 512 in the support module (`build.sh` substitutes only `norder`); change it in `support_mod.f90` and the `chk` list of `build.sh`, and supply matching `grid/` files. Python: `--nglv` |
| Inner tolerance | Fortran: `EPS1` in `build.sh` (`epsilon1`; `epsilon2` stays 1.0E-5). Python: `--epsilon1` |
| Step length and number of steps | Fortran: `DT`, `NSTEPS` in `build.sh` (any integer `DT`; non-uniform steps from a times file, see [Input files](inputs.md#times-file-optional-non-uniform-time-steps)). Python: `--steps`, `--dt` in years (`run_ais.py`, default 1), `--dt` in years, `--nfiles` (`run_outer.py`) |
| Inner-iteration start | Fortran: `SLEVOLVE_INNER=off|extrap|anderson` (`inner_read_options`; code in `sl_solver`). Python: `--inner` |
| Warm start across outer passes | Fortran: `SLEVOLVE_REUSE=warm` (`warm_read`, `warm_write`). Python: `run_outer.py --reuse warm` |
| Horizontal output | Fortran: `SLEVOLVE_HORIZ=1` (`horiz<N>.npy`) or `2` (also `Psi<N>.npy`); in snapshot mode add `horiz`/`Psi` to `SLEVOLVE_SNAP_FIELDS`. Python: `--horizontal` or `--horizontal-grid` (`horiz.gradient`); `python/tools/horiz_grad.py` for Fortran output |
| Transform backend | Fortran: `SHT=blas|otf` (production) or `exact|ref` (verification) and `FFT=fftw|fft991`, `SYM=1|0` in `build.sh`; `SLEVOLVE_FFT` at run time (`sht_env_init`). Python: `--sht numpy|ducc` (`run_ais.py`, `make_sht` in `sht_ducc.py`); `run_outer.py` always uses ducc0 |
| Rotation, marine check | Fortran: positional `TPW`, `CHECKMARINE` arguments of `build.sh`. Python: `--tpw`, `--checkmarine` |
| Output selection and format | Fortran: `SLEVOLVE_SNAP`, `SLEVOLVE_SNAP_FIELDS`, `SLEVOLVE_FMT`, `SLEVOLVE_TEXT` (`io_mod`). Python: `--save-fields` (`run_ais.py`), `--save-steps` (`run_outer.py`) |
| New output field | Fortran: add a `pp_write_grid` (grid) or `pp_write_cplx` (spectral) call in `sl_solver` after the commit, next to the existing `tgrid`/`ocean` writes; to make it selectable in snapshot mode, add its name to `pp_field_bit` in `io_mod.f90`. Python: add the field to the dict returned by `OptSeaLevelSolver.step` (and `SeaLevelSolver.step`) and save it in the step loop of `run_ais.py`/`run_outer.py` with `np.save` or `npyio.save_spectral` |
| Physical constants | Fortran: the support module (default `rhoi=920`, `radius=6.371E6`; `LEGACY_CONSTANTS=1` in `build.sh` gives the legacy 910 / 6.371229E6). Python: `--rhoi`, `--radius` |

**Next:** [Glossary](reference/glossary.md).
