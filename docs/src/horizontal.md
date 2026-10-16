# Horizontal displacement output

This page explains how to switch on, write and post-process the optional horizontal (tangential) surface displacement and states its conventions.

Applies to the Fortran solver (`fortran/src/`) and the Python solver (`python/`). The feature is
off by default. Switching it on changes no other output (see [Checks made during development](#checks-made-during-development)).

## 1. Implementation

The horizontal displacement is read from the same recursive modal states R_k as the radial
displacement (paper Eqs. 2–10). Only the weights change: `l_k/s_k` replaces `h_k/s_k`.

```
Psi_lm(t_n) = T_l [ l^E_l L_lm + sum_k (l_k/s_k) R_k,lm ]
            + delta_l2 [ l^T_2 Lambda_2m + sum_k (l^T_k/s_k) R^Lambda_k,m ] / g
u_h = grad_h Psi:   u_theta = dPsi/dtheta,   u_phi = (1/sin theta) dPsi/dphi
```

| | Fortran | Python |
|---|---|---|
| switch | `SLEVOLVE_HORIZ=0` (default) / `1` / `2` | `run_ais.py`, `run_outer.py`: `--horizontal`, `--horizontal-grid` |
| weights | `exact_spec%a_hl = resl/s`, `a_thl = tresl/s` (degree 2), allocated only when the switch is on (solver module `fortran/src/sl_model_mod.f90`, after the spectrum build) | `ExactMemory.horizontal_weights` (`solver.py`) |
| read-out | `modal_state_mod::state_read_horizontal` (type-bound procedure). It reads the current states and leaves them unchanged. It is called after the R synthesis and before `commit`. | `ExactMemory.read_horizontal` (`[k,l,m]` layout), `ExactMemoryFast.read_horizontal` (`[l,k,m]`), `SeaLevelSolver.horizontal_lm`. These are called at the same point in `solver.py` and `solver_opt.py` `step()` and returned as `r["horiz"]`. |
| coefficient output | `OUTPUT_SLM/horiz<N>.npy` (`<c16`, `(L+1,L+1)`, fortran_order, spharmt U convention, like `dS_converged`) | `horiz<N>.npy` (`run_outer.py`: `horiz<N>_it<k>.npy` at `--save-steps` and the last step) |
| grid output | `SLEVOLVE_HORIZ=2`: `Psi<N>.npy` (the synthesized scalar) | `--horizontal-grid`: `utheta<N>.npy`, `uphi<N>.npy` |
| snapshot mode | field names `horiz`, `Psi` for `SLEVOLVE_SNAP_FIELDS`; `all` includes them. They are written only when `SLEVOLVE_HORIZ` is set. | — |

Files:
* `python/slevolve/horiz.py`: `gradient(U, nlat, nlon)` gives (u_theta, u_phi). It maps U to ducc0 a_lm in the
  same way as `sht_ducc.py` and runs the ducc0 spin-1 `synthesis_2d(..., mode="GRAD_ONLY")` on the GL grid.
* `python/tools/horiz_grad.py`: `horiz<N>.npy` → `utheta<N>.npy`, `uphi<N>.npy`. This is the
  post-processing path for Fortran output. `--selftest` checks against analytic fields.
* `python/tools/horiz_history_check.py`: verification only: checks the recursive read-out against the explicit full-history sum.

**What is left to post-processing.** The Fortran code does not compute the gradient. A clean d/dtheta
needs either a dP_lm/dtheta table, which the on-the-fly Legendre backend does not have (`getuv` and
`cosgrad` need the `pnm`/`hnm` tables, which are not allocated in the `SHT=otf` (default), `SHT=blas` and `SHT=exact` builds; only `SHT=ref` keeps them), or the
recurrence sin(theta) dP_l^m/dtheta = l eps_{l+1,m} P_{l+1}^m − (l+1) eps_{l,m} P_{l-1}^m. That
recurrence needs degree L+1, which the existing `spec2spat` (degree ≤ L) cannot synthesize. Fortran
therefore writes Psi_lm (and optionally Psi on the grid), and `python/tools/horiz_grad.py` forms
u_theta and u_phi.

## 2. Conventions (checked in code and data)

* **Love-number file.** The radial path uses `elast(1,l)` = h^E and `resh` = h_k without scaling.
  k is stored as l·k (`ke = elast(3)/l`, `resk/l`). The l (tangential) numbers are stored the same
  way, as l·l: `elast(2,l)` and `resl`. The code divides `resl` and `tresl` by l, and the elastic
  value as well: l^E_l = `elast(2,l)/l`, l^T_2 = `telast(2,2)/2`. Evidence from the
  Case 1 file `prem_512.l60K2C…lm22`:
  * degree 2, load: 0.0481/2 = 0.024;
  * degree 2, tidal: 0.168/2 = 0.084 (the PREM tidal l_2);
  * high degree: l·l^E → about 1.15;
  * for degrees 2, 10 and 100, elastic + Σ_k residue_k/s_k reproduces the asymptotic (fluid) line
    for h, l and k alike. The elastic value and the residues therefore carry the same factor.

  Note: in this reader `elast(3,:)` is l·k (it feeds `ke`); the l column is `elast(2,:)`.
* **Normalization.** T_l = 4πa³/(M(2l+1)) is the factor of the radial R (eq. B16), so T_l L_lm is the
  load potential divided by g. Then u_r = h·(Φ/g) and u_h = l·∇₁(Φ/g), where ∇₁ is the dimensionless
  surface gradient. Psi and u are in metres. The instantaneous load term is exactly the one used for
  R: `rhoi*deltaicestar + rhow*deltaS_prev + rhow*load_inc_lm`. The rotational term uses h → l and
  `exact_trr` → the tidal-l sum over `r_lambda`.
* **Directions.** theta is colatitude, so u_theta > 0 is southward (u_north = −u_theta). u_phi > 0 is
  eastward. Grid rows run north to south and column 0 is phi = 0. Since the elastic l_l > 0 at all degrees l >= 2 of the shipped Love-number files (degree 1 depends on the frame), the elastic motion
  points toward added load and away from ice loss.
* **Spectral convention.** The spharmt U convention, as in `dS_converged`. The Python raw-mode
  backend works in V = (−1)^m U and converts before output.

## Checks made during development

Verification checks on short runs of Case 1 and Case 3 (Fortran: `SHT=otf`, FFTW,
8 threads; Python: ducc, `--opt all`). They were not repeated for this release.

1. **Default off changes nothing.**
   * Case 1, 20 steps (`./exe 20`, inner off): a build without the feature, the build with the switch
     unset and the build with `SLEVOLVE_HORIZ=2` give identical output files (154 common files,
     including `tgrid20`, `R20`, `G20`, `dS_converged20`, `ocean20` and `numiter`).
   * Case 3, two outer passes (`./exe 50 2 3`, extrap + warm, from a segmented initial topography): all
     75 output files identical, including `pred_pres_topo_2/3`, `tgrid0_*` and `numiter*`.
   * Python: `run_ais.py` 20 steps with and without `--horizontal --horizontal-grid` gives identical
     `tgrid20`, `dS_converged20` and `numiter`; `run_outer.py` one pass (Case 3) with and without
     `--horizontal` gives identical `pred_pres_topo_1`.
2. **Fortran vs Python Psi_lm**, Case 1, 20 steps. The largest difference over all coefficients and all
   20 steps, relative to max|Psi_lm|, is **8.9e-13** (`--opt all`); 2.2e-13 to 4.2e-12 for the base
   `SeaLevelSolver` path (steps 1 and 5). The two codes' `tgrid20` differ by 9.1e-13 m. The Fortran
   `Psi20.npy` equals the ducc0 synthesis of `horiz20.npy` to 6.7e-16.
3. **Recursion vs explicit full-history sum** (`python/tools/horiz_history_check.py`, Case 1, 20 steps,
   L = 511). The quantity Σ_k (l_k/s_k) R_k,lm(t_n) was compared with
   Σ_{j<n} ΔL_lm(j) Σ_k (l_k/s_k)(1 − e^{−s_k (t_n−t_j)}), built from the committed increments. The
   largest relative difference is **1.4e-12**, below 2e-13 from step 3 on. The rotational (tidal-l)
   term agrees to **2.5e-14**.
4. **Gradient.** `horiz_grad.py --selftest` (64×128, L = 63, analytic Psi with m = 0, 1, 2) gives
   max errors of 2.8e-12 for u_theta and 3.3e-13 for u_phi.
5. **Magnitude sanity check.** Case 1, 550 × 1 yr, Fortran, `SLEVOLVE_INNER=extrap`,
   `SLEVOLVE_HORIZ=2`, gradient from `horiz_grad.py`. Nothing was tuned.
   * Maximum south of 60°S: |R| = 380.5 m (77.4°S, 252.4°E). The maximum horizontal displacement is
     |u_h| = 56.4 m (81.3°S, 241.5°E), where the ice thickness changed by −1136 m. The ratio of the
     maxima is **0.15**.
   * RMS over the cells with |R| > 10 % of max|R| (10,565 cells): u_h/R = **0.17**.
   * At the step-0 ice edge south of 60°S (3,786 cells): RMS u_h = 12.3 m and RMS R = 26.6 m, a ratio
     of **0.46**. The median pointwise ratio is 0.79. Near the edge, R falls off faster than its
     gradient.

   The horizontal displacement is a fraction of the radial one, largest at the margins of the unloading
   region.

**Next:** [File conversion](conversion.md).
