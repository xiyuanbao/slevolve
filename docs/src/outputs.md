# Output fields

This page defines every output file of the two solvers: what it holds, its sign convention, units,
reference time and frame. All statements are taken from the code that writes the files
(`python/slevolve/solver.py`, `python/run_ais.py`, `python/run_outer.py`, and the output section of
`sl_solver` / `sl_solver_init` in `fortran/src/sl_model_mod.f90`). The binary containers
are specified in [File formats](formats.md).

## 1. Conventions common to all fields

**Grid layout.** Every grid file is a `(512, 1024)` array `a[j, i]`. Row `j` is the Gauss-Legendre
latitude `GLlat[j]`, running from north (+89.73 deg) to south; column `i` is longitude `i * 360/1024`
degrees east (0 to 359.65). Longitude varies fastest in memory. `tools/make_grid.py` writes the
latitudes and longitudes as text ([Input files](inputs.md#grid-files)).

**Step index.** File index `N` is the step number. Ice file `<icemodel>N` is the ice thickness at the
end of step `N`, i.e. at time `N * dt` after the start (with a times file: at time `t_N`, line N + 1 of
the file; [Input files](inputs.md#times-file-optional-non-uniform-time-steps)); ice file 0 is the
initial state. The step numbers in the output names are the same with and without a times file.

* `N = 0`: the initial state, written before the first step (`tgrid0`, `ocean0`, `beta0`, and in
  Fortran `dS_converged0` = 0 and `TPW0` = 0).
* `N = 1`: the state after the first step (time `dt`, or `t_1` with a times file), computed with ice file 1.
* `N = k`: the state after step `k` (time `k * dt`, or `t_k`).

**Reference time.** `tgrid`, `ocean` and `beta` are absolute (the state at the end of step N). The
changes `R`, `G` and `dS_converged` are **cumulative since step 0** (the start of the run, or of the
outer pass), not increments of one step.

**Sign of topography and sea level.** The main variable is the topography T: the elevation of the
solid surface (bedrock, without the ice) relative to the sea surface, in metres, positive up. The
ocean is where T < 0. The solver updates it as `T_N = T_0 - dSL_N`, where `dSL_N` is the sea-level
change since step 0 (positive = sea-level rise). So the relative sea-level change since step 0 is
`-(tgrid<N> - tgrid0)`.

**Degree-1 frame.** `R` and `G` contain degree 1 in the frame of the Love-number file; SLEvolve
applies no frame conversion. The shipped Love-number files are in the CE frame (k_1 = 0). The Case 1
file is CE in its elastic and fluid limits but carries a degree-1 initialisation defect in its
residues, which changes degree-1 `R` and `G` but not sea level or topography
(`lovenum/DEGREE1_NOTES.md`).

**Masks.** `ocean` and `beta` are 0/1 fields, stored as `<f8` (Python `--io-compat legacy`: uint8).

## 2. Grid and spectral fields (per step N)

| file | meaning | units | sign / reference |
|---|---|---|---|
| `tgrid<N>` | topography T at the end of step N: solid surface (bedrock) minus sea surface. `tgrid0` is the initial topography of the run (or of the outer pass) | m | positive up (land), negative in the ocean; absolute |
| `ocean<N>` | ocean function C: 1 where `tgrid<N>` < 0, else 0. It does not look at the ice: marine cells under grounded ice have C = 1. The ocean that carries water load is C * beta | 0/1 | absolute |
| `beta<N>` | grounded-ice mask: 1 where the grounded ice of step N (after the marine check, flotation on the topography of the last inner re-guess) is zero (below machine epsilon), else 0, so the ocean under floating shelves has beta = 1. `beta0` is the same mask for ice file 0 on the initial topography | 0/1 | absolute |
| `R<N>` | radial displacement of the solid surface since step 0, including the rotational term | m | positive up; cumulative since step 0; degree 1 in the Love-number frame |
| `G<N>` | change of the sea surface (geoid plus the spatially uniform term from mass conservation) since step 0, `G = dSL + R` | m | positive up; cumulative since step 0; degree 1 in the Love-number frame |
| `dS_converged<N>` | spherical-harmonic coefficients of the change of ocean (water) thickness since step 0, after the inner iteration of step N has converged. The water load is `rhow * dS` | m of water | cumulative since step 0 |
| `horiz<N>`, `Psi<N>`, `utheta<N>`, `uphi<N>` | optional horizontal displacement, see [Horizontal displacement](horizontal.md) | m | since step 0 |
| `dSLdt<N>`, `dRdt<N>` | optional analytic rates of sea level and radial displacement at the end of step N, see [below](#analytic-rates-optional) | m/yr | instantaneous rate at t_N |

The sea-level change itself, `dSL = G - R`, is not written as a separate file; it equals
`tgrid0 - tgrid<N>`.

**Spectral convention of `dS_converged` (and `horiz`).** Complex `(L+1, L+1)` array `a[l, m]`,
`0 <= m <= l`, zero above the diagonal, in the convention of the legacy `spharmt` module
(`python/slevolve/sht.py`): `a[0, 0]` is the global mean of the field (to a relative 1.7e-9, from the truncated
`1.41421356` literal of the legacy source), and the field is
`f = sum_l [ a[l,0] c_l0(mu) + 2 Re sum_{m>=1} a[l,m] c_lm(mu) exp(i m phi) ]` with
`c_lm = 1.41421356 * P~_lm`, `P~_lm` orthonormal on [-1, 1] without the Condon-Shortley phase.
`python/slevolve/sht.py` (`SHT.synthesis`) and `sht_ducc.py` synthesise it.

### Analytic rates (optional)

Off by default; Fortran `SLEVOLVE_RATES=1`, Python `run_ais.py --save-fields ...,dSLdt,dRdt` or
`run_outer.py --rates`. Output only: every other output is unchanged. The formulas are in
[Method](method.md#analytic-present-day-rates-optional).

| file | meaning | units |
|---|---|---|
| `dSLdt<N>` (`.npy` grid) | rate of the sea-level change dSL (sea surface minus solid surface) at the end of step N, including the rotational term and, in degree 0, the rate of the conservation term on the ocean function of step N. `-dSLdt` is the rate of topography change | m/yr |
| `dRdt<N>` (`.npy` grid) | rate of the radial displacement (uplift rate), including the rotational term | m/yr |
| `TPWrate<N>` (Fortran, text) / `TPWrate<N>.txt` (Python) | with rotation: the rates of the three records of `TPW<N>`, in the same layout: il-dot (9 values, kg m^2/yr), m-dot (3 values, 1/yr; m1-dot, m2-dot is the polar-motion (true polar wander) rate in rad/yr), Lambda-dot (18 values, m^2 s^-2 /yr) | as stated |

The rate of the sea surface is `dGdt = dSLdt + dRdt`. Present-day J2-dot and the true-polar-wander
rate follow from these: with \(\dot G_{20}\) the (2,0) coefficient of `dGdt` (spharmt convention
of `dS_converged`, whose (l,0) functions are the 4-pi-normalised Legendre functions to a relative
1.7e-9) and \(\dot\Lambda_{20}\) from `TPWrate<N>`,
\(\dot J_2=-\sqrt5\,(\dot G_{20}-\dot\Lambda_{20}/g)/a\) (the centrifugal potential removed), and the
TPW speed is \(|\dot m_1+i\,\dot m_2|\) rad/yr (times \(180/\pi\cdot10^6\) for degrees per Myr).

They are the rates at \(t_N\) read from the modal states, not the mean rate over the last step that a
difference `(X<N> - X<N-1>) / dt` gives. With the second-order placement (`RAMP=1`, `--ramp`) they
are second order in the step length; with the default placement only the last increment is read as
accruing linearly over its step.

## 3. Scalar outputs

| quantity | meaning | units | Fortran file | Python |
|---|---|---|---|---|
| ocean area | area of C (`ocean<N>` = 1) | m^2 | `ocean_area` | `run_log.json` `scalars[N].ocean_area` |
| ice-free ocean area | area of C * beta | m^2 | `oceanBeta_area` | `scalars[N].oceanBeta_area` |
| GMSL change over the ocean | area-weighted mean of the sea-level change since step 0 over C (`gmslc_ocn`) | m (rise positive) | `gmslc_ocnArea` | `scalars[N].gmslc_ocn` |
| GMSL change over the ice-free ocean | the same mean over C * beta (`gmslc_ocnBeta`) | m | `gmslc_ocnBetaArea` | `scalars[N].gmslc_ocnBeta` |
| ice volume | integral of the ice thickness used as load (after the marine check when it is on, i.e. grounded ice only) | m^3 of ice | `ice_volume` | `scalars[N].ice_volume` |
| inner iterations | iterations of the ocean-load update to reach `epsilon1` at each step | count | `numiter` (one line per step 1..N) | `numiter.npy` |
| time | `k * dt` since the start; with a times file `t_k` as given in it | yr | `times` (`0, dt, 2 dt, ...`, one line per step 0..N, `ES14.4E2`; with a times file `t_0, t_1, ..., t_N` in `ES24.16E3`) | not written; `run_log.json` has `dt_years`, and with `--times` `times_years` (`t_0..t_N`); `outer_log.json` has `times_years` |
| rotation | see below | | `TPW<N>` | not written |
| rotation rates (optional) | rates of the `TPW<N>` quantities, [see above](#analytic-rates-optional) | | `TPWrate<N>` (with `SLEVOLVE_RATES=1`) | `TPWrate<N>.txt` (with `dSLdt` or `dRdt` in `--save-fields`) |

The Fortran scalar files are text, one value per line (`ES14.4E2`), line 1 = step 0, line N+1 = step N
(the GMSL files start with 0 at step 0). They are appended at every step and are written also in
snapshot mode (`SLEVOLVE_SNAP`). Python stores the same quantities in `run_log.json` under `scalars`:
entry 0 is step 0 (areas and ice volume only), entry N is step N. `nmelt` (Fortran) is the number of
earlier load steps taken into account at each step, and `elapsed_wall_time`, `elapsed_cpu_time` are timings.

**`TPW<N>`** (Fortran, rotation on): three records of the converged rotational quantities after step N
(`TPW0` = zeros), printed in Fortran column-major order:

1. 9 values: `il(3,3)`, the perturbation of the inertia tensor (only `il(1:3,3)` is filled), kg m^2;
2. 3 values: `mm(1:3)`, the perturbation of the rotation vector, dimensionless (m1, m2 polar motion
   in radians, m3 relative change of the spin rate);
3. 18 values: `lambda(3,3)`, the degree-0 and degree-2 coefficients of the rotational driving
   potential, as (real, imaginary) pairs, m^2 s^-2.

All three are totals since step 0.

## 4. Outer iteration (ice age runs)

| file | meaning |
|---|---|
| `tgrid0_<p>` | initial topography T0 used in outer pass p (m) |
| `pred_pres_topo_<p>` | topography at the last step of pass p, the predicted present topography (m) |
| `tgrid<k>_it<p>.npy` | Python `run_outer.py`: topography at step k of pass p, for the steps in `--save-steps` |
| `T0_next.npy` | Python: the initial topography for the next pass (Picard: `tgrid0_<p> + (T_obs - pred_pres_topo_<p>)`; `--outer secant` uses the secant update) |
| `outer_log.json` | Python: per pass the misfit `rms_pred_minus_obs` (area-weighted RMS of `pred_pres_topo_<p> - T_obs`, m), `max_abs`, the inner iterations per step and the wall time |
| `numiter_<p>` | Fortran `run_outer_mem`: copy of `numiter` for pass p |

The Fortran files `pred_pres_topo_<p>`, `tgrid0_<p>` are written at the last step of each pass; in
snapshot mode only when `SLEVOLVE_OUTER` is set.

## 5. Which files each driver writes

**Fortran** (`OUTPUT_SLM/`, `.npy` by default):

* step 0: `tgrid0`, `ocean0`, `beta0`, `dS_converged0`, `TPW0`, and the first line of every scalar file;
* each step N: `tgrid<N>`, `ocean<N>`, `beta<N>`, `dS_converged<N>`, `R<N>`, `G<N>`, `TPW<N>`, one more
  line in each scalar file;
* with `SLEVOLVE_SNAP=<list>` only the listed steps (step 0 always) and the fields in
  `SLEVOLVE_SNAP_FIELDS` (default `tgrid`) are written as grids; `snap_index.txt` lists them. The Case 1
  example (`SLEVOLVE_SNAP=550 SLEVOLVE_SNAP_FIELDS=tgrid`) writes `tgrid0.npy` and `tgrid550.npy`;
* with `SLEVOLVE_RATES=1`: `dSLdt<N>`, `dRdt<N>` (snapshot fields `dSLdt`, `dRdt`) and `TPWrate<N>`
  at the steps where either rate grid is written.

**Python `run_ais.py`** (`--out`):

* step 0: `tgrid0.npy`, `ocean0.npy`, `beta0.npy` (always);
* each step N: the fields named in `--save-fields` (default `tgrid,ocean,beta`; also `R`, `G`, `dS`,
  `dSLdt`, `dRdt`); `dS` is written as `dS_converged<N>.npy`; with `dSLdt` or `dRdt` and rotation also
  `TPWrate<N>.txt`;
* at the end: `numiter.npy`, `step_wall_seconds.npy`, `run_log.json` (settings, timings, `scalars`).

**Python `run_outer.py`** (`--out`): the files of §4, plus `horiz<k>_it<p>.npy` with `--horizontal`
and `dSLdt<k>_it<p>.npy`, `dRdt<k>_it<p>.npy`, `TPWrate<k>_it<p>.txt` with `--rates`.

To convert: `python/tools/npy2text.py` (text), `python/tools/npy2nc.py` (NetCDF with latitude,
longitude and time axes; the time comes from `OUTPUT_SLM/times` when present, otherwise from
`step * --dt`; for a Python run with non-uniform steps use the `times_years` of `run_log.json`).

## 6. Files written, by solver

Sections 1-4 define the meaning, sign, units and reference time of every output; this section is a
short checklist of the files, by solver (section 5 gives the same per driver, with step detail).

### Fortran, per step N

* Grids `tgrid<N>`, `ocean<N>`, `beta<N>`, `R<N>`, `G<N>` and the spectral field `dS_converged<N>`.
  By default every field is written at every step, as `.npy` (`SLEVOLVE_FMT=npy`, `SLEVOLVE_BIN=1`).
  Text copies are written only with `SLEVOLVE_TEXT=1`.
* Restrict grid output with `SLEVOLVE_SNAP` / `SLEVOLVE_SNAP_FIELDS`. `OUTPUT_SLM/snap_index.txt`
  lists what was written.
* Scalar text files, always written and appended per step: `numiter`, `elapsed_wall_time`, `TPW<N>`,
  area, GMSL and ice-volume files.

### Fortran, at the last step

* `pred_pres_topo_<k>` and `tgrid0_<k>` (outer pass k).
* These are not written in snapshot mode unless `SLEVOLVE_OUTER` is set.
* `run_outer_mem` also copies `numiter` to `numiter_<k>` and prints `== pass k wall … s`.

**Optional horizontal displacement** (off by default; [Horizontal displacement output](horizontal.md))

* Fortran `SLEVOLVE_HORIZ=1`: `horiz<N>.npy` (spectral Psi_lm); `SLEVOLVE_HORIZ=2`: also `Psi<N>.npy`
  (Psi on the grid). Python `--horizontal`: `horiz<N>.npy`; `--horizontal-grid`: also `utheta<N>.npy`,
  `uphi<N>.npy`. Grid components for Fortran output: `python/tools/horiz_grad.py horiz<N>.npy`.

**Optional analytic rates** (off by default; [above](#analytic-rates-optional))

* Fortran `SLEVOLVE_RATES=1`: `dSLdt<N>.npy`, `dRdt<N>.npy`, `TPWrate<N>`. Python: `dSLdt`, `dRdt` in
  `--save-fields` (`run_ais.py`) or `--rates` (`run_outer.py`).

### Python

* `tgrid<N>.npy` plus the fields in `--save-fields` (default `tgrid,ocean,beta`).
* `numiter.npy` and `run_log.json`.
* `run_outer.py` writes the outer-pass files as `.npy`, plus the steps listed in `--save-steps`.

**Next:** [File formats](formats.md).
