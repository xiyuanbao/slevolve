# Input files

This page specifies the run-directory layout and the ice, topography, grid and Earth-model (Love-number) input files, their units and names.

## Run directory layout

```
RUN/
  INPUT_FILES/            (often a symlink)
    icemodel/   <icemodel><k>[.npy|.bin]        k = 0..N   ice thickness, m
    earthmodel/ <planetmodel>                   Love-number spectrum (text)
    grid/       GLlat_512.txt, GLlon_512.txt    grid description (tools/make_grid.py)
    others/     <topomodel>[.npy|.bin]          initial / present topography
                times                           optional: times of the ice files (non-uniform steps)
  OUTPUT_SLM/   outputs
  ICELOAD_SLM/  (must exist; legacy directory)
  exe
```

## Grid files

**Grid files.** `GLlat_512.txt` holds the 512 Gauss-Legendre latitudes in degrees, north to south,
and `GLlon_512.txt` the 1024 longitudes `i * 360/1024` in degrees east, one value per line. Neither
solver reads them (Fortran passes the names but opens no grid file); they document the grid, and
`tools/prepare_case3_nadi.py` uses them to check the ICE-NADI-1 grid. `tools/make_grid.py --outdir
INPUT_FILES/grid` writes them; `examples/case1_forward.sh`, `tools/prepare_case2_inputs.sh` and
`tools/prepare_case3_nadi.py` call it when the files are missing. Compared with the OSF archive files,
`GLlon_512.txt` is identical and 58 of the 512 lines of `GLlat_512.txt` differ in the last printed
digit (maximum 5.3e-7 degrees).

**Where the file names come from**

* Fortran: the names are compiled into `user_specs_mod` (`icemodel`, `planetmodel`, `topomodel`,
  `topo_initial`, `gridfolder`, `timearray`).
* `build.sh` sets `gridfolder = 'INPUT_FILES/grid/'` and the names from `PLANETMODEL`, `ICEMODEL`,
  `TOPOMODEL`, `TOPO_INITIAL` (or a `CASE` preset, [Command-line options](options.md#case-presets)). Without them: `icemodel='AISload_gd_Run85_0_0_'`,
  `planetmodel='prem_512.l60K2C.sum18p6.dum19p2.tz19p4.lm22'`, `topomodel = topo_initial = 'topo0_archive'`.

## Grid and spectral formats (summary)

See [File formats](formats.md) (also `fortran/IO_FORMATS.md`) for the full specification.

* **Grids.**
  * Container: `.npy` v1.0, `'<f8'`, shape `(512, 1024)`, `fortran_order=False`.
  * Element order: the stream order of the legacy text files (row = Gauss-Legendre latitude from
    north, longitude fastest).
  * The payload is identical to the legacy headerless `.bin`.
  * Ice inputs may be `<f4` only where this is lossless. `text2npy --f4` checks per file and falls back
    to `<f8`.
* **Spectral.** `dS_converged<N>.npy`: `'<c16'`, shape `(L+1, L+1)`, `fortran_order=True`, so that
  `a[l,m] = dS(l,m)`.
* **Input lookup.** For every grid `<path>` the solver tries `<path>.npy`, then `<path>.bin`, then the
  text file `<path>`. Python follows the same order (`npyio.read_grid`).
* **Legacy-format record files** (not converted): the checkpoint `exact_memory_state_<i>.bin` and the
  file-mode warm store `reuse_<i>.bin`.

## Times file (optional; non-uniform time steps)

Without a times file every step has the same length (Fortran: the build-time `DT`; Python: `--dt`), and
ice file `k` is the ice at time `k * dt`. A times file gives the time of every ice file instead, so
that the steps can differ in length:

* **Format.** Plain text, one number per line: `t_0, t_1, ..., t_N`, the time of ice file `k` in
  **years**, for `k = 0..N` (N + 1 lines, as many as ice files). Any real number format that Fortran
  list-directed input and `numpy.loadtxt` read (`0`, `12.5`, `1.25e4`). `t_0` may be 0 or any other
  value (e.g. `-26000` for 26 ka before present); the times must be strictly increasing (time runs
  forward from file 0, the initial state, to file N). Extra lines after `t_N` are ignored by Fortran;
  Python `run_ais.py` uses the first `--steps + 1` values, `run_outer.py` and
  `init_topo_segmented.py` require exactly one value per ice file.
* **Step length.** Step `j` (from file `j-1` to file `j`) has length `dt_j = t_j - t_(j-1)`. The load is
  constant over each step, as for uniform steps.
* **Fortran.** Detected automatically at `INPUT_FILES/others/times`; the environment variable
  `SLEVOLVE_TIMES=<path>` selects another file and `SLEVOLVE_TIMES=off` ignores it
  ([Options](options.md#fortran-run-time-environment-variables-read-at-run-time-by-the-executable)).
  With a times file, build with `DT` unset (it defaults to 1 and then only counts steps, `NSTEPS`). Without one, the build-time `DT` (any integer) is the uniform step. The executable
  prints `non-uniform time steps from <file>: ...` when a times file is used.
* **Python.** `--times FILE` on `run_ais.py`, `run_outer.py` and `init_topo_segmented.py` (`--dt` is
  then ignored). Python does not look for the file by itself.
* **A uniform times file** (`t_k = t_0 + k * dt`) with an integer or exactly representable `dt` gives
  the same result as the run without a times file exactly (zero difference).

Example (8 steps of 1, 2, 5, 10, 20, 50, 100 and 200 yr):

```
0
1
3
8
18
38
88
188
388
```

## Love-number (Earth model) file

The file is the plain-text output of a collocation Love-number calculation (Mitrovica and Peltier
1992); degree-1 conventions are described in [Degree-1 Love numbers](degree1.md). It is read once per run in
the solver module `fortran/src/sl_model_mod.f90` (block "Read in Love numbers" in `sl_solver`, run once via `pp_love_ready`), as in the legacy
code:

```
(header line, skipped)
for each degree l = 1..norder:
  l  K_l                                   (degree, number of modes)
  s_1 .. s_K                               inverse relaxation times, 1/kyr
  h^E  l*l^E  .  .  l*k^E                  elastic load Love numbers (l and k multiplied by l)
  h^F  l*l^F  .  .  l*k^F                  fluid limit (t -> infinity), read but not used
  r_1 .. r_K        (h residues)
  l_1 .. l_K        (l residues, multiplied by l)
  r'_1 .. r'_K      (k residues, multiplied by l)
  for l >= 2: the same five lines for the tidal Love numbers (h^T, l*l^T, l*k^T, residues)
```

* The l scaling was checked on the Case 1 file: degree 2 gives 0.0481/2 = 0.024 (load l'_2) and
  tidal 0.168/2 = 0.084 (PREM tidal l_2), and for every degree and all three numbers
  `elastic + sum_k residue_k/s_k` reproduces the asymptotic (fluid-limit) line, so elastic values
  and residues carry the same factor l.
* The code divides the k and l values by l, and sets `he = h^E`, `ke = k^E/l`, `r = h residues`,
  `rprime = k residues`, `rT`, `rprimeT` (tidal), as in Kendall et al. (2005).
* Units: the state decay factor is `exp(-s * dt_years / 1000)` (`modal_state_mod::update_decay`),
  so `s` is in 1/kyr. The recursion uses the residue weights `(r'-r)/s`, `r/s` and, for degree 2, the
  rotational terms `r'/s`, `r'_T/s`, `(r'_T-r_T)/s`.
* Only degrees 1..`norder` are read; the file may contain more.
* `modal_state_mod::spectrum_load` reads a different, pre-reduced spectrum file
  (`LMAX MAXMODES`, then per mode `s` and six weights). The solver does not call it; the spectrum is
  built in memory from the Love-number file.
* Python reads the same Love-number file with `solver.LoveNumbers`.

How to write a file for your own Earth model: [Love numbers](lovenum.md).

**Next:** [Output fields](outputs.md).
