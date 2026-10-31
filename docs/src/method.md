# Method overview

This page summarises the formulation, the recursive modal state and its placement rule, the changes to the iteration starts, and where each change is in the code, and lists the core equations of the paper with every symbol explained.

## Formulation

SLEvolve solves the gravitationally self-consistent sea-level equation on a spherically symmetric
(1-D) viscoelastic Earth. It uses the normal-mode (Love-number) formulation of Kendall et al. (2005)
as implemented in the legacy code of Han et al. (2022), called the *legacy solver* in this guide.

The legacy solver evaluates the viscous response at every time step as an explicit sum over the whole
load history. The cost per step therefore grows with the step number, and the total cost grows as
O(N^2). SLEvolve replaces that sum with a **recursive normal-mode state**: one state per mode, degree
and order, advanced by one decay factor per step. The cost per step is then independent of the length
of the history, and the total cost is O(N). The recursion is exact. It covers all six history terms of
the legacy solver: sea level, radial displacement, and the four degree-2 rotational terms.

There are two implementations.

* **Fortran** (`fortran/src/`). Derived from the legacy code.
* **Python** (`python/`). Written from scratch following the same formulation, built on NumPy, Numba
  and ducc0.

**What is exact and what is not.**

* These parts are exact: with them, the paper's version (Zenodo v1.0.0, doi 10.5281/zenodo.23422438)
  reproduces the legacy solver to round-off (< 1e-11 m in final topography, identical inner-iteration
  counts and masks); the current code additionally applies the corrected floating-ice convention, the corrected sign of the axial rotation component m3 and the corrected sign of the order-1 rotation records (see [Provenance](provenance.md#differences-from-the-runs-of-han-et-al-2022)):
  * the recursive state,
  * resident I/O,
  * the transform reductions,
  * the Legendre and FFT backends,
  * OpenMP.
* These parts change only the **starting point of an iteration**, never the equations:
  * the extrapolated inner start,
  * the warm start across outer passes,
  * the partially relaxed initial topography.

  The converged result of each step then agrees with the legacy result within the inner tolerance
  `epsilon1`, not exactly.

## Core equations

The equations below are those of the paper (Sect. 2), in its notation: \(\Delta\) denotes a change
since the first time \(t_0\), \(\delta\) an increment over one step. Equation numbers are those of the
paper.

### Sea-level equation (Eq. (1) of the paper)

$$
\Delta S_j = \Delta\mathcal{SL}_j\,C_j\beta_j - T_0\,(C_j\beta_j-C_0\beta_0),\qquad
\Delta\mathcal{SL}_j=\Delta\mathcal{G}_j-\Delta\mathcal{R}_j+\frac{\Delta\Phi_j}{g}
$$

* \(t_0<t_1<\dots<t_N\): the time steps; \(\delta t_j=t_j-t_{j-1}\) the step length (may be non-uniform).
* \(\Delta S_j\): change of ocean height since \(t_0\) (the ocean load, as water thickness) at step \(j\).
* \(\Delta\mathcal{SL}_j\): change of sea level (sea surface minus solid surface) since \(t_0\).
* \(\Delta\mathcal{G}_j\): geographically variable part of the change of the sea-surface equipotential (geoid change).
* \(\Delta\mathcal{R}_j\): change of the radial position of the solid surface (radial displacement).
* \(\Delta\Phi_j/g\): the uniform term, a spatially constant shift of the sea surface fixed by conservation of
  ice-plus-ocean mass; \(g\) is the surface gravity.
* \(C_j\): ocean function, 1 where sea level is positive (sea surface above the solid surface), 0 elsewhere.
* \(\beta_j\): ice function, 0 where grounded ice covers the surface, 1 elsewhere. The marine check: marine ice
  counts as grounded only where its thickness exceeds the water depth times \(\rho_w/\rho_I\); thinner
  ice floats and is not part of the load. \(\beta_j\) is the grounded-ice mask of Kendall et al.
  (2005, eqs 43-44): it is 0 exactly where the grounded ice \(I^*_j\) is non-zero, so the ocean under
  floating ice shelves is part of the water-loading ocean \(C_j\beta_j\) and of the ocean area of
  the conservation term. A floating shelf together with the water beneath it weighs as much as a
  water column up to the sea surface, so a change of sea level under a shelf changes the load as in
  the open ocean. The flotation check uses the current topography iterate
  \(T_j^{(i)} = T_0 - \Delta\mathcal{SL}_j^{(i)}\) of the inner iteration (the topography that also
  defines the re-guessed ocean function); the first iterate of a step uses the converged topography
  of the previous step. Mask update rule: \(I^*_j\), its spectrum (the ice-load increment),
  \(\beta_j\) and \(C_j\beta_j\) are re-evaluated with every re-guess of the ocean function; the
  ice analysis is repeated only when the grounded flag of at least one cell changed, and the grounded
  mask is updated at most three times per step (a cell at flotation can otherwise make the mask
  alternate between two states); after that it is kept for the rest of the step. The viscous state is
  still advanced once per step with the converged load of the previous step. With the marine check
  off, \(I^*_j = I_j\) and only the mask definition applies.
* \(T_0\): topography at \(t_0\); prescribed in a forward run, found by the outer iteration in an ice age run.

The equation is nonlinear because \(C_j\) and \(\beta_j\) depend on the sea level being solved for. At
each step an inner (Picard) iteration on the ocean-load increment \(\delta S_j\) runs until its relative
change between iterates falls below \(\epsilon_1\) (default \(10^{-5}\)).

### Response to the load (Eqs. (2)-(5) of the paper)

For each spherical-harmonic coefficient \((\ell,m)\):

$$
\Delta\mathcal{SL}_{\ell m}(t_j)=T_\ell\left[E_\ell\,\Delta L_{\ell m}(t_j)+V_{\ell m}(t_j)\right],\qquad
T_\ell=\frac{4\pi a^3}{M_e(2\ell+1)},\qquad E_\ell=1+k^E_\ell-h^E_\ell
$$

* \(a\), \(M_e\): radius and mass of the Earth.
* \(\Delta L_{\ell m}(t_j)=\sum_{n=0}^{j}\delta L^n_{\ell m}\): cumulative load change; \(\delta L^n_{\ell m}=\rho_I\,\delta I^{*\,n}_{\ell m}+\rho_w\,\delta S^n_{\ell m}\)
  the load increment of step \(n\), with \(\delta I^{*}\) the increment of grounded ice thickness
  (marine ice included where grounded) and \(\rho_I\), \(\rho_w\) the ice and water densities.
* \(E_\ell\): elastic sea-level response, from the elastic load Love numbers \(h^E_\ell\) (radial) and
  \(k^E_\ell\) (potential).
* \(V_{\ell m}\): viscous (hereditary) response. At degree 2 the rotational feedback adds
  \(g^{-1}[E^T_2\,\Delta\Lambda_{2m}+V^T_{2m}]\), with the tidal Love numbers and the centrifugal
  potential perturbation \(\Lambda\), formed in the same way.

The viscoelastic Love numbers are finite sums of decaying exponentials (Eq. (4)):

$$
h_\ell(t)=h^E_\ell\,\delta(t)+\sum_{k=1}^{K_\ell} r^\ell_k\,e^{-s^\ell_k t},\qquad
k_\ell(t)=k^E_\ell\,\delta(t)+\sum_{k=1}^{K_\ell} r^{\ell\prime}_k\,e^{-s^\ell_k t}
$$

* \(K_\ell\): number of modes at degree \(\ell\) (16 to 26 in the benchmark Earth models).
* \(s^\ell_k>0\): inverse relaxation time of mode \(k\); \(r^\ell_k\), \(r^{\ell\prime}_k\): its residues
  (modal amplitudes) for radial displacement and potential; \(\delta(t)\): Dirac delta.

A load increment persists after it is applied, so its viscous response is the time integral of the
impulse response, \((r/s)(1-e^{-st})\) per mode: the response to a step load increment is the elastic
term plus a sum over modes of \(a_k(1-e^{-s_k t})\). The full-history form of the viscous term is
(Eq. (5))

$$
V_{\ell m}(t_j)=\sum_{n=0}^{j-1}\beta_\ell(t_j-t_n)\,\delta L^n_{\ell m},\qquad
\beta_\ell(t)=\sum_{k=1}^{K_\ell} a^\ell_k\left(1-e^{-s^\ell_k t}\right),\qquad
a^\ell_k=\frac{r^{\ell\prime}_k-r^\ell_k}{s^\ell_k}
$$

* \(\beta_\ell(t)\): viscous step-response kernel of the sea level (not the ice function \(\beta_j\)).
* \(a^\ell_k\): modal weight of mode \(k\) for sea level. The radial displacement uses the weights
  \(r^\ell_k/s^\ell_k\) of the same modes. This sum costs \(O(j)\) at step \(j\), \(O(N^2)\) in all; it is
  what the legacy solver evaluates (six such sums: sea level, radial displacement and four degree-2
  rotational terms).

### Recursive modal state (Eqs. (8) and (10) of the paper)

$$
V_{\ell m}(t_j)=\sum_{k=1}^{K_\ell}a^\ell_k\,R^k_{\ell m}(t_j),\qquad
R^k_{\ell m}(t_j)=d^\ell_k(\delta t_j)\,R^k_{\ell m}(t_{j-1})+\bigl(1-d^\ell_k(\delta t_j)\bigr)\,\Delta L_{\ell m}(t_{j-1}),\qquad
d^\ell_k(\delta t_j)=e^{-s^\ell_k\,\delta t_j}
$$

* \(R^k_{\ell m}(t_j)\): the modal state of mode \(k\), one field per mode, degree and order;
  \(R^k_{\ell m}(t_0)=0\). It is defined as \(\Delta L_{\ell m}(t_{j-1})-Q^k_{\ell m}(t_j)\) with
  \(Q^k_{\ell m}(t_j)=\sum_{n=0}^{j-1}e^{-s^\ell_k(t_j-t_n)}\delta L^n_{\ell m}\) (Eq. (7)).
* \(d^\ell_k(\delta t_j)\): one-step decay factor of mode \(k\) over the step length \(\delta t_j\); a
  non-uniform step changes only \(d\).
* The update is algebraically identical to the full-history sum for any sequence of increments; it costs
  one multiply-add per mode and coefficient per step, \(O(N)\) in all. The rotational terms carry
  their own states with the same update.

### Placement rule (Sect. 2.4 of the paper)

An increment applied at \(t_j\) has zero viscous age at \(t_j\) (\(\beta_\ell(0)=0\)), so \(V_{\ell m}(t_j)\)
depends only on the converged \(\Delta L(t_{j-1})\). Hence:

1. **advance-and-read**: before the inner iteration of step \(j\), the states are advanced once with
   the update above and \(V_{\ell m}(t_j)\) is read from them; it stays fixed during the iteration;
2. the inner iteration updates only the elastic and uniform terms, \(C_j\), \(\beta_j\) and \(\delta S_j\);
3. **commit**: after convergence, \(\Delta L_{\ell m}(t_j)=\Delta L_{\ell m}(t_{j-1})+\delta L^j_{\ell m}\).

Trial iterates never enter the stored history, so every iteration sees exactly the viscous response of
the full-history sum, and both formulations take the same inner iterations.

### Second-order load placement (optional)

The placement rule above is the step placement of Kendall et al. (2005): the increment of step \(j\)
is applied at \(t_j\) as a step load, so its viscous response starts at \(t_j\). On average the
response then lags by half a step, and the time stepping is first order in \(\delta t\). As an option,
each increment can instead be taken to accrue linearly over its step \((t_{j-1},t_j]\) (a ramp). The
viscous kernel integrated over the ramp is

$$
\int_{t_{j-1}}^{t_j}\left(1-e^{-s(t-\tau)}\right)\frac{d\tau}{\delta t_j}
=1-\varphi(s\,\delta t_j)\,e^{-s(t-t_j)},\qquad
\varphi(x)=\frac{1-e^{-x}}{x}\quad(t\ge t_j)
$$

for each mode, so the ramp only replaces the unit amplitude of the decaying exponential by
\(\varphi(s^\ell_k\delta t_j)\). \(\varphi(x)=1-x/2+x^2/6-\dots\) tends to 1 as \(\delta t_j\to0\) (the step
response is recovered); the code uses this series below \(x=10^{-4}\). The elastic term and the
direct sea-level terms (\(E_\ell\), the ocean function, the conservation term) are unchanged, because a
ramp is complete at \(t_j\). The response is exact for a load that is piecewise linear in time and
second order for a smooth load.

In the recursion the state stays one field per mode, and the extra term is local to the step:

1. **advance-and-read** is unchanged. The current increment has the viscous response
   \(\sum_k a^\ell_k c^\ell_k\,\delta L^j_{\ell m}\) at \(t_j\), with \(c^\ell_k=1-\varphi(s^\ell_k\delta t_j)\). It is
   linear in \(\delta L^j\), like the elastic term, so it enters the inner iteration as an addition to
   the elastic coefficient of the current increment:
   \(T_\ell[E_\ell\Delta L_{\ell m}+V_{\ell m}+C_\ell\,\delta L^j_{\ell m}]\), \(C_\ell=\sum_k a^\ell_k c^\ell_k\)
   (radial displacement: the radial weights);
2. **commit** adds \(c^\ell_k\,\delta L^j_{\ell m}\) to every modal state \(R^k_{\ell m}\) before \(\Delta L\) is
   updated. One step later the increment contributes \(d\,c+(1-d)=1-d\,\varphi\), the ramp kernel, and the
   identity carries over to all later steps.

\(c^\ell_k\) uses the length \(\delta t_j\) of the step in which the increment accrued, so non-uniform
steps (a times file) need nothing else; the coefficients are recomputed once per distinct step length.
The rotational states take the same two changes. The polar-motion equation becomes implicit in the
current increment and is solved in closed form (\(i=1,2\)):

$$
m_i=\frac{\alpha\left[\dfrac{(1+k^E_2)\,I_{i3}+\Sigma^I_{i3}+c_I\,\delta I_{i3}}{C-A}+\dfrac{\Sigma^m_i-c_m\,m_i^{\,j-1}}{k_f}\right]}{1-\alpha\,c_m/k_f},
\qquad \alpha=\frac{1}{1-k^{T,E}_2/k_f}
$$

with \(\Sigma^I\), \(\Sigma^m\) the read modal sums, \(c_I=\sum_k a^I_k c_k\) and \(c_m=\sum_k a^m_k c_k\)
(degree-2 weights of the moment of inertia and of the tidal potential); and
\(m_3=-[(1+k^E_2)\,I_{33}+\Sigma^I_{33}+c_I\,\delta I_{33}]/C\)
(Mitrovica & Wahr 2011, eq. 15). The viscous tidal terms of sea level and radial displacement gain
\(c_{tsl}\,\delta\Lambda\) and \(c_{trr}\,\delta\Lambda\). The cost is one multiply-add per mode and coefficient at
the commit (the advance has three) and nothing in the inner iteration beyond one multiply-add per
coefficient.

The ramp removes the lag of the viscous response. Events that the step sequence resolves only to one
step remain first order: a marine basin that floods in the step in which its ice cover vanishes (the
shoreline term), and ice that grounds or floats. They leave a first-order remainder on margins and
forebulges.

**Default off.** The paper and the shipped `reference/` results use the step placement of Kendall
et al. (2005), so the default keeps it and every default output is unchanged. The option is the
build switch `RAMP=1` (Fortran) and `--ramp` (Python `run_ais.py`, `run_outer.py`). The method was
derived and verified in a separate study.

### Analytic present-day rates (optional)

Rates are read from the modal states at the end of a step (left derivative), the current increment
taken to accrue linearly over the step. Inside the step every state obeys
\(dR^k/dt=s_k(P-R^k)\), with \(P\) the cumulative driver, so

$$
\frac{d\,\Delta\mathcal{SL}_{\ell m}}{dt}\Big|_{t_j}=T_\ell\left[E_\ell\frac{\delta L^j_{\ell m}}{\delta t_j}
+\sum_k a^\ell_k s^\ell_k\left(\Delta L_{\ell m}(t_j)-R^k_{\ell m}(t_j)-c^\ell_k\,\delta L^j_{\ell m}\right)\right]
$$

(radial displacement: \(h^E_\ell\) and the radial weights), with \(R^k\) the state after the advance of
step \(j\). The degree-2 rotational rates differentiate the polar-motion equation, which is explicit in
\(\dot m\): the current-increment term becomes \(s_k(m-r^m_k)\) with the ramp-consistent state
\(r^m_k+c_k\delta m\), \(\dot I=\delta I/\delta t_j\), and \(\dot\Lambda\) follows by the chain rule. Degree 0 of
\(d\mathcal{SL}/dt\) is the rate of the conservation term on the ocean function \(C^*_j\) of the end of the
step: \(\bigl(-(\rho_I/\rho_w)\,\delta I^*_{00}/\delta t_j-\overline{C^*_j\,\dot{\mathcal{SL}}}\bigr)/\overline{C^*_j}\)
(overbar: global mean; \(\dot{\mathcal{SL}}\) without degree 0).

A two-point difference of the outputs, \((X_j-X_{j-1})/\delta t_j\), is the mean rate over the last step,
centred at \(t_j-\delta t_j/2\): first order in \(\delta t\) whatever the placement. The analytic read-out
is the rate at \(t_j\). With the ramp the viscous memory is read without lag; what remains is that the
current increment, in particular its ocean-load and grounded-ice parts, which the sea-level equation
determines, is taken as linear over the last step. That term is first order, with a small constant: on
a piecewise-linear ice history (synthetic Earth, \(L=256\), steps of 1 kyr to 125 yr) the analytic
uplift rate, \(\dot J_2\) and polar-motion rate with the ramp were 8 to 100 times closer to the converged
values than two-point differences at the same step. The read-out is also available with the step
placement, where the history is that of the step placement and only the last increment is read as a
ramp; it then inherits the lag of the step placement and was not more accurate than a two-point
difference in that test, so use it with the ramp. The rates are output only (Fortran `SLEVOLVE_RATES=1`, Python `--save-fields ...,dSLdt,dRdt` or `run_outer.py
--rates`; [Output fields](outputs.md#analytic-rates-optional)).

### Extrapolated inner start (Sect. 2.5 of the paper)

$$
\delta S^{(0)}_j=\delta t_j\left[\dot S_{j-1}+\left(\dot S_{j-1}-\dot S_{j-2}\right)\frac{\delta t_j+\delta t_{j-1}}{\delta t_{j-1}+\delta t_{j-2}}\right],\qquad
\dot S_i=\frac{\delta S_i}{\delta t_i}
$$

* \(\delta S^{(0)}_j\): first guess of the ocean-load increment of step \(j\) (superscript 0: initial iterate).
* \(\dot S_i\): rate of the converged increment of step \(i\).
* For uniform steps this is \(2\,\delta S_{j-1}-\delta S_{j-2}\). The first two steps use the eustatic
  guess (the mass of the ice increment spread uniformly over the ocean), which the legacy solver uses
  at every step. Only the starting point changes; the converged result agrees within \(\epsilon_1\).

### Outer iteration and partially relaxed initial topography (Eqs. (11) and (12) of the paper)

$$
T_0\leftarrow T_0-\left(T_{\mathrm{pred}}-T_{\mathrm{obs}}\right)
$$

* \(T_{\mathrm{pred}}\): present topography predicted by a pass (a full forward run); \(T_{\mathrm{obs}}\): observed
  present topography. The misfit is the area-weighted RMS of \(T_{\mathrm{pred}}-T_{\mathrm{obs}}\).

First-pass initial topography from the partially relaxed response:

$$
T_0=T_{\mathrm{obs}}+\sum_i\mathcal{SL}^{\mathrm{seg}}_i,\qquad
h^{\tau_i}_\ell=h^E_\ell+\sum_k\frac{r^\ell_k}{s^\ell_k}\left(1-e^{-s^\ell_k\tau_i}\right)
$$

* \(i\): monotonic segments of the global ice-volume history (split at its extrema; changes below 1 % of
  the volume range ignored). The superscript seg marks a quantity of one whole segment (neither a one-step
  increment δ nor a change since \(t_0\), Δ).
* \(I^{\mathrm{seg}}_i\): net change of grounded ice over segment \(i\); \(\tau_i\): its mean age at the end of
  the run, weighted by the ice-volume change of each step.
* \(h^{\tau_i}_\ell\), \(k^{\tau_i}_\ell\) (the latter with \(r^{\ell\prime}_k\)): Love numbers partially relaxed to age \(\tau_i\).
* \(\mathcal{SL}^{\mathrm{seg}}_i\): sea-level response to \(I^{\mathrm{seg}}_i\) with these Love numbers; the sea-level
  equation is solved once for the summed response, with the water load evaluated with the segment-averaged
  Love numbers \(h^w_\ell=\sum_i w_i h^{\tau_i}_\ell\), \(w_i=|V^{\mathrm{seg}}_i|/\sum_j|V^{\mathrm{seg}}_j|\)
  (\(V^{\mathrm{seg}}_i\): net ice-volume change of segment \(i\)). Implemented in `python/init_topo_segmented.py`.

## What changed relative to the legacy code

Section numbers refer to Sect. 2 of the paper.

| # | Change | Where (Fortran / Python) | Arithmetic changed? |
|---|---|---|---|
| 1 | **Recursive modal state.** Two calls per step: *advance-and-read* before the inner iteration (decays the state by one step and returns the viscous sums), *commit* after convergence (adds the converged load increment). This placement rule keeps the recursion identical to the full-history sum inside the nonlinear iteration. | `fortran/src/modal_state_mod.f90`, module `modal_state_mod`: `state_init`, `update_decay`, `state_advance_and_read`, `state_commit`, `state_save_checkpoint` / `state_load_checkpoint`. Python: `kernels.py::advance_and_read`, `commit_add`; `solver_opt.py::ExactMemoryFast` | Exact (reformulation of the same sum; round-off level) |
| 2 | **Resident state and binary I/O.** One process runs all steps. Earth model, modal state, step-0 fields and transform tables stay in memory. Grids are read and written as binary (`.npy`, see [Input files](inputs.md)). | solver I/O module in `fortran/src/` (`pp_read_grid`, `pp_write_grid`, `pp_write_cplx`), `npy_io_mod.f90` (`read_grid`, `read_ice_cache`, `write_grid`, `write_snapshot`); driver `run_from_zero_plus.f90`. Python: `npyio.py` (`read_grid`, `save_grid`, `save_spectral`, `IceDir`) | No. `SLEVOLVE_COMPAT` restores the legacy re-reads (diagnostic) |
| 3 | **Fewer transforms per step**, 5n+6 → 3n+2 (n = inner iterations), grouped into 2n+1 Legendre sweeps. | `opt_read_options`, flags `diag, dsl, lazyrg, icebatch, rt, merge` (`SLEVOLVE_OPT`). Python: `--opt raw,diag,dsl,numba,sums` (`solver_opt.py::parse_opts`) | Order of operations only. `SLEVOLVE_OPT=none` switches the reductions off |
| 4a | **Legendre backend: table.** Precomputed P_lm applied as a BLAS matrix product (DGEMM), with N–S symmetry. | `fortran/src/spharmt_otf.f90`: `spharm_blas_ana`, `spharm_blas_syn`, `sht_sym`; build `SHT=blas` | Order of operations only |
| 4b | **Legendre backend: on-the-fly recurrence.** No table; recurrence in m, vectorized over latitude, OpenMP over m. Used for the benchmark runs in the paper. | `spharmt_otf.f90` + `otf_kernels.inc`; build `SHT=otf` | Order of operations only |
| 5 | **FFTW** in place of FFT99. One plan for the whole 512×1024 array. | `sht_fft`, `sht_env_init` (`SLEVOLVE_FFT=fftw|fft991`); build `FFT=fftw` | Order of operations only |
| 6 | **OpenMP** over orders m and over spectral loops. | `!$omp` directives; build `OMP=N` | Output does not depend on the thread count to round-off (paper Sect. 2.7) |
| 7 | **Extrapolated inner start.** The first guess of the step-n ocean-load increment is extrapolated linearly in time from the two previous converged increments. For uniform Δt this is `2 dS(n-1) - dS(n-2)`. For variable Δt the rate is extrapolated between step midpoints and scaled by Δt(n). | `SLEVOLVE_INNER=extrap`. Python: `--inner extrap` | Changes the iteration start. The result agrees within `epsilon1` |
| 7b | Anderson mixing (experimental, not used in the paper) | `SLEVOLVE_INNER=anderson`, `SLEVOLVE_AA_M`, `SLEVOLVE_AA_MASK`; Python `--inner anderson` | Changes the iteration path |
| 8 | **Warm start across outer passes.** Pass k starts step n from the converged dS_n and ocean mask of pass k−1. | Fortran: `SLEVOLVE_REUSE=warm` (file `OUTPUT_SLM/reuse_<k>.bin`, or in memory with `run_outer_mem`). Python: `run_outer.py --reuse warm` | Changes the iteration start only |
| 9 | **Partially relaxed (segmented) initial topography.** The ice-volume history is split into monotonic segments. Each segment's net load change is responded to with Love numbers relaxed for its mean age τ_i. The result is the T0 guess for outer pass 1. | Python only: `python/init_topo_segmented.py` (function `segmented_guess`, CLI; 5-point running mean, legs below 1 % of the range pruned, `--threshold`) | Changes the outer starting point only |
| 0 | **Floating-ice convention.** \(\beta_j\) is the grounded-ice mask of Kendall et al. (2005, eqs 43-44), and the flotation check uses the current topography iterate (see the definition of \(\beta_j\) above). | `sl_model_mod.f90`: `gr_ground`, `gr_reguess`, `gr_update_now`. Python: `solver.py::_grounded_reguess` | Yes (a correction of the legacy convention, which excluded every ice-covered cell from the water-loading ocean and tested flotation against the initial topography of the pass) |
| 10 | **In-process outer iteration driver.** One process runs passes K0..K1. Earth model, transform tables, FFTW plans, state allocation, ice grids and the warm-start store stay resident. | `fortran/src/run_outer_mem.f90` (`./exe LAST K0 K1`). One pass per process: `run_outer_plus.f90` (`./exe LAST ITERSL`). Python: `run_outer.py` | No. Picard update `tinit_0 ← tinit_0 + (T_obs − T_pred)` as in the legacy code. `SLEVOLVE_OUTER=secant` is an experimental alternative |
| 11 | **Second-order load placement and analytic rates** (optional, off by default; [above](#second-order-load-placement-optional)). | `modal_state_mod.f90`: `update_ramp`, `state_ramp_commit`, `state_read_rate`; `sl_model_mod.f90`: `rate_spectra`, `rate_grids`, the `ramp` branches of the rotation, sea-level, radial and horizontal assembly; build `RAMP=1`, run time `SLEVOLVE_RATES=1`. Python: `solver.py::ramp_factor`, `ExactMemory.read_rate`, `SeaLevelSolver._rate_spectra`; `kernels.py::ramp_commit`, `rate_sums`, `ramp_add` | `RAMP=1` changes the result (second-order time stepping); the rates are output only |

The longitude FFT routine `fft99b` of the transform module of the legacy code is used in a corrected
form (`fortran/src/spharmt_otf.f90`).

## Horizontal (tangential) displacement: convention

For a spherically symmetric Earth the surface horizontal displacement is the surface gradient of a
scalar, read out from the same modal states R_k as the radial displacement with different weights:

```
Psi_lm(t_n) = T_l [ l^E_l L_lm(t_n) + sum_k (l_k / s_k) R_k,lm(t_n) ]          (l >= 1, metres)
            + delta_l2 [ l^T_2 Lambda_2m(t_n) + sum_k (l^T_k / s_k) R^Lambda_k,m(t_n) ] / g
u_theta = dPsi/dtheta ,   u_phi = (1 / sin theta) dPsi/dphi
```

* T_l = 4 pi a^3 / (M (2l+1)) (eq. B16, the factor of the radial R), L_lm the total surface load
  (ice + ocean, kg m^-2) in the same instantaneous form as R (`rhoi*deltaicestar + rhow*deltaS_prev +
  rhow*load_inc_lm`), Lambda the rotational driving potential (same terms as the R rotation with h replaced
  by l). l^E_l = `elast(2,l)/l`, l_k = `resl(k,l)` (after the division by l), l^T_2 = `telast(2,2)/2`,
  l^T_k = `tresl(k,2)`. So T_l L_lm is the load potential divided by g, and
  `u_r = h (Phi/g)` and `u_h = l grad_1 (Phi/g)` (grad_1 = dimensionless surface gradient).
* theta is colatitude, so u_theta > 0 points **south**; u_north = -u_theta; u_phi > 0 points east.
  Grid rows run north to south and column 0 is phi = 0 (the solver grid).
* With l > 0 the elastic horizontal motion points toward a positive load (away from a region of ice
  loss).
* Coefficients are written in the spharmt (U) convention of `dS_converged` (both codes; Python
  converts from its raw-mode V convention).
* Fortran: `SLEVOLVE_HORIZ` (0 default). The weights `a_hl = l_k/s_k`, `a_thl` are added to the
  spectrum only when the switch is on; `modal_state_mod::state_read_horizontal` reads the
  current states between advance-and-read and commit and does not change them. With 0 (or with any
  value) every other output is identical to the build without the feature.
* The grid gradient is not computed in Fortran (the on-the-fly Legendre backend has no dP/dtheta
  table); it is formed from `horiz<N>.npy` by `python/tools/horiz_grad.py` (ducc0 spin-1 gradient
  synthesis, `python/slevolve/horiz.py`), which is also what Python `--horizontal-grid` uses.
* Usage: [Horizontal displacement output](horizontal.md); checks: [Horizontal displacement output](horizontal.md#checks-made-during-development).

**Next:** [Benchmark provenance and differences from Han et al. (2022)](provenance.md).
