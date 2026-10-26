# Benchmark provenance and differences from Han et al. (2022)

This page documents where the data of the three benchmark cases come from, the known issues of the archived files, and how the benchmark set-up differs from the archived runs of Han et al. (2022).

The input data of the three benchmark cases are not redistributed with SLEvolve. The scripts in
`tools/` convert the archived files into solver inputs; the commands are in the
[Running the benchmark cases](benchmarks.md).

## Case 1 -- Antarctic forward run (550 x 1 yr from 1950 CE)

**Ice forcing.** The extended RCP8.5, MISI-only (no cliff failure, no hydrofracturing) Antarctic
ice-sheet run of DeConto et al. (2021), 1950-2500 CE, on a 10 km grid, taken from the Dryad archive of
Yousefi et al. (2022), doi:10.5061/dryad.ncjsxkswm (run `85_MISI`). The total (grounded plus floating)
ice thickness is sampled onto the 512 x 1024 Gauss-Legendre grid with `tools/build_case1_total_ice.py`
(sampling map in `tools/regrid.py`).

**Earth model and topography.** From the OSF archive of Han et al. (2022), doi:10.17605/OSF.IO/8PTFM,
folder `Data/AIS_550yr/`:

* Love-number file `prem_512.l60K2C.sum18p6.dum19p2.tz19p4.lm22`. This is their 'BLVZ' profile:
  60-km lithosphere; viscosities 10^18.6, 10^19.2 and 10^19.4 Pa s at 60-200, 200-400 and 400-670 km
  depth. The file name ends in `lm22` (lower mantle 10^22 Pa s), whereas DeConto et al. (2021) print
  1 x 10^19 Pa s for the lower mantle (see "Differences from the runs of Han et al. (2022)" below).
  The file is included unchanged in `lovenum/`.
* Initial (1950 CE) topography `topo0_archive`, taken from the archived output run.

**What the OSF archive contains and what the archived run used.**

* The 551 archived ice files `AISload_gd_Run85_0_0_<k>` hold *grounded* ice only (prefix `_gd_`; they
  pass the flotation test against the archived bed with zero failures).
* The archived output run was computed from the *total* (grounded plus floating) ice of the same
  ice-sheet run. Evidence: (i) driven by the archived grounded ice, neither the legacy code nor SLEvolve
  reproduces the archived topography (RMS 2.4 m at step 1 growing to 11.4 m at step 550; Antarctic RMS
  28 m); (ii) the archived step 1 contains a one-off jump (ocean mean -1.52 m, far-field -1.80 m, up to
  +26 m along the Ronne and Ross shelf fronts) that corresponds to the removal of the floating shelves
  present in the initial state (total-minus-grounded volume 6.40e14 m^3) but absent from the later
  loads; (iii) a proxy shelf field reproduces this jump (pattern correlation 0.98).
* The total-ice fields are not in the OSF archive and are rebuilt from the Dryad archive. The grounded
  part of the rebuilt ice, sampled onto the Gauss-Legendre grid, matches the archived grounded files
  cell by cell: correlation 0.99976-0.99981, 99.0 % of cells identical, volume ratio 1.00000-1.00005,
  identical maximum thickness, at every step checked (0, 1, 50, 100, 150, 200, 300, 400, 550). The
  remaining ~1 % of cells differ by a nearest-neighbour resampling offset at the ice margin.
* The rebuilt files keep the prefix `AISload_gd_Run85_0_0_` expected by the compiled-in file names,
  although they hold the total ice.

**How the benchmark is run.** With the total-ice forcing and the solver's own flotation check
(`checkmarine = .true.`), which is the configuration of the archived run. Do not combine total-thickness
ice with `checkmarine = .false.`. The archived grounded-ice loads, as distributed, can be used for an
exactness comparison with the legacy code.

**Legacy baseline.** The legacy code of Han et al. (2022) re-reads the ice history as text at every
step, in the layout and number format of the archived files. The text files are converted once to the
binary format for SLEvolve; the legacy inputs are not altered.

## Case 2 -- two glacial cycles (1200 x 0.2 kyr from 240 ka)

**Source.** Same OSF archive, `Data/SeaLevelModel_TwoGlacialCycle_Simulations/INPUT_NHIS2GC.zip`
(`iceload0..1200`): the Northern Hemisphere ice-sheet simulation of Han et al. (2021) repeated over two
cycles, combined with ICE-6G_C for Antarctica (Argus et al. 2014; Peltier et al. 2015). Only the last
122 kyr are a realistic history; the first cycle is a repeat. Present topography `etopo2_512_ice6gC`
from `Data/CoupledModel_OneGlacialCycle_Simulations/INPUT_Others/`. Input tree:
`tools/prepare_case2_inputs.sh`.

**Known issue.** The Earth-model file distributed for the glacial-cycle cases,
`Data/CoupledModel_OneGlacialCycle_Simulations/INPUT_Others/prem_coll_512.l120C.ump5.lm5`, is corrupt in
the archive (1 110 755 NUL bytes in every copy; the zip CRC verifies, so it was uploaded in this state).
Love numbers for the stated model (PREM elastic and density structure, 120-km lithosphere, upper mantle
5e20 Pa s, lower mantle 5e21 Pa s) are used instead. They were computed with the collocation method
of Mitrovica and Peltier (1992) (after Peltier 1985) and are included in `lovenum/`
(`prem_512.l120C.ump5.lm5`). The `case2` build preset uses this file name for `PLANETMODEL`; copy it
into `INPUT_FILES/earthmodel/`.

## Case 3 -- ICE-NADI-1 deglaciation (50 x 0.5 kyr from 25 ka; Lucas et al. 2026)

**Source.** Zenodo, doi:10.5281/zenodo.21464139. Input tree: `tools/prepare_case3_nadi.py`
(`icemodel/nadi_0..50`, `others/topo_pres`, `grid/`). Earth model Earth_1 of Lucas et al. (2026)
(96-km lithosphere, upper mantle 3e20 Pa s, lower mantle 1e21 Pa s). Its Love numbers were computed
with the collocation method of Mitrovica and Peltier (1992) (after Peltier 1985) and are included in
`lovenum/` (`prem_512.l96C.ump3.lm1`). Copy the file into `INPUT_FILES/earthmodel/`.

## Grid and formats

All cases use the 512 x 1024 Gauss-Legendre grid of Han et al. (2022), maximum degree 511. Grid files
are NumPy `.npy` (fp64, little-endian, `(nlat, nlon)` in the stream order of the text files); spectral
fields are `.npy` complex `(L+1, L+1)` in Fortran order. See [File formats](formats.md).

## Differences from the runs of Han et al. (2022)

Deliberate deviations of the benchmark set-up from the archived runs of Han et al. (2022), with the
reason. "Not verified" marks statements that could not be confirmed from the archive or the paper.

* *Floating-ice convention.* The legacy code excludes every ice-covered cell, floating shelves
  included, from the water-loading ocean and tests flotation against the initial topography of the
  pass. SLEvolve uses the grounded-ice mask of Kendall et al. (2005, eqs 43-44), so the ocean under
  floating shelves carries water load and belongs to the ocean area of the conservation term, and
  tests flotation against the current topography. In the paper, SLEvolve and the legacy solver were
  compared with the same floating-ice treatment (that of the legacy code), which established the
  round-off agreement; the released code differs from those runs only by this convention, and
  `reference/` is computed with it.

**Case 1 (Antarctic forward run).**

* *Total-ice forcing instead of the archived grounded ice, with the marine check on.* The 551 archived
  files `AISload_gd_Run85_0_0_<k>` hold grounded ice only, but the archived output run was computed
  from the total (grounded plus floating) ice (evidence above). The total ice is not in the OSF
  archive and is recovered from the Yousefi et al. (2022) Dryad archive of the same ice-sheet run; its
  grounded part matches the archived files cell by cell (correlation 0.99976-0.99981). With total ice
  the solver's flotation check (`checkmarine = .true.`) removes the floating part, as in the archived run.
* *Step-1 shelf-removal initialization not mimicked by default.* The archived run let the floating
  shelves of the initial state enter the elastic and mass-conservation terms (and the rotational
  feedback) at step 1 but not the viscous memory. SLEvolve and the released legacy code include them
  in all terms alike, which leaves 2.05 m area-weighted RMS against the archived 550-yr topography on a
  12.3 m signal (0.024 m with that initialization mimicked); these values are for the paper's version, with the legacy floating-ice
  convention. The released legacy code behaves like
  SLEvolve, so the archived run differs from the released code in this initialization; the default
  follows the released code.
* *Love-number file as distributed* (`prem_512.l60K2C.sum18p6.dum19p2.tz19p4.lm22`, unchanged). The
  file name says lower mantle 10^22 Pa s, whereas DeConto et al. (2021) print 1 x 10^19 Pa s for this
  profile. It is used as distributed because the archived run used it.
* Initial topography `topo0_archive` taken from the archived output run.

**Case 2 (two glacial cycles).**

* *Love numbers recomputed.* The archived Earth-model file `prem_coll_512.l120C.ump5.lm5` is corrupt
  (1 110 755 NUL bytes in every copy). Love numbers for the stated model (PREM, 120-km lithosphere,
  upper mantle 5e20 Pa s, lower mantle 5e21 Pa s) were computed with the collocation method of
  Mitrovica and Peltier (1992) (`lovenum/prem_512.l120C.ump5.lm5`). Results therefore differ from the
  archived run by the difference of the Love numbers (not quantified here).
* *One pass from T_obs.* Only the first outer pass (initial topography = present topography
  `etopo2_512_ice6gC`) is run, for the timing comparison; Han et al. (2022) ran the outer iteration.
  Number of passes in the archived run: not verified.
* *Antarctic ice from ICE-6G_C, as in Han et al. (2022)* (the ice files are used as archived).

**Case 3 (ICE-NADI-1).** Not a case of Han et al. (2022).

* Ice history ICE-NADI-1 (Lucas et al. 2026) with their Earth model Earth_1 (96-km lithosphere, upper
  mantle 3e20 Pa s, lower mantle 1e21 Pa s); Love numbers computed with the collocation method
  (`lovenum/prem_512.l96C.ump3.lm1`).
* *Present topography substituted:* `others/topo_pres` is `etopo2_512_ice6gC` from the Han et al.
  (2022) OSF archive (`tools/prepare_case3_nadi.py`), not a topography distributed with ICE-NADI-1.
  The ICE-NADI-1 inputs used here do not include a present topography on the solver grid.
* Ice age run with the outer iteration on the initial topography, started from the partially relaxed
  (segmented) guess with warm-started passes, instead of the three passes from the present
  topography of Kendall et al. (2005).

**All cases.**

* Constants: ice density 920 kg m^-3 and Earth radius 6.371e6 m, the same as Han et al. (2022)
  (the default; the older legacy values 910 and 6.371229e6 are available with `LEGACY_CONSTANTS=1`).
* Maximum degree 511 on the 512 x 1024 grid. The legacy code's default is `norder = 512`, and
  Han et al. (2022) describe their runs as expansions up to degree and order 512; degree 511
  is the largest degree resolved by 512 Gauss latitudes (the ducc0 transform requires lmax <= nglv - 1).
* Inner tolerance epsilon = 1e-5, the value set in the released legacy code (`epsilon1 = 1.0E-5`).
* Uniform steps in each case (1 yr, 200 yr, 500 yr), with the full step-wise load history kept
  exactly. Han et al. (2022) also propose a time-window scheme that coarsens older loads; whether the
  archived runs used it: not verified.
* Rotational feedback on (`tpw = .true.`) in all cases, as in the released legacy code.

**Next:** [Degree-1 Love numbers and reference frames](degree1.md).
