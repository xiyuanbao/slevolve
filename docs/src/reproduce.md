# Reproducing the published results

This page states which run of the [benchmark cases](benchmarks.md) corresponds to which result of the paper and what it should give.

## Case 1: Antarctic forward run

**Expected output**

* `RUN_case1/OUTPUT_SLM/tgrid550.npy` (Fortran) or `OUT_case1_py/tgrid550.npy` (Python).
* Against the archived run: 2.05 m area-weighted RMS on a 12.3 m signal (0.024 m with the
  shelf-removal initialization mimicked); see [Benchmark provenance](provenance.md) and paper Sect. 4.4.

**Reproduction of the archived output (paper Sect. 4.4).** With the total ice, SLEvolve reproduces the
archived 550-yr topography to 2.05 m area-weighted RMS on a 12.3 m signal. The residual is the step-1
shelf removal: the archived run let the removed shelves enter the elastic and mass-conservation terms
(and the rotational feedback) at step 1 but not the viscous memory. With that initialization mimicked,
the agreement is 0.024 m RMS. The legacy code and SLEvolve agree with each other to machine precision
(< 1e-11 m) in every configuration.

## Case 2: two glacial cycles

**Scope in the paper.** One pass (the first outer iteration, from the present topography) is timed.

## Case 3: ICE-NADI-1

**Outer iteration.** The complete ice age calculation is run with the outer iteration on the initial
topography; the legacy accuracy is the misfit after the three passes from the present topography used
by Kendall et al. (2005).

**Reference values.** The misfit per outer pass of the reference run is in `reference/case3.json`
(`misfit_rms`, `misfit_max`; [Verifying the installation](verify.md#benchmark-reference-check)).

## Exactness and timing

The exactness criteria against the legacy solver (paper Sect. 2.8) are in [Release notes](release_notes.md#checks-made-during-development); the wall times of the reference runs are on [Running the benchmark cases](benchmarks.md).

**Next:** [Command-line options](options.md).
