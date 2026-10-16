# Setting up a new case

This page lists the steps from your own ice history, topography and Earth model to a run and a check
of its result. The commands use the Python drivers, which need no rebuild; the Fortran equivalents
are given in [6. Running](#6-running). Input formats: [File formats](formats.md). Output fields:
[Output fields](outputs.md).

!!! note "Step length"
    The steps must match the times of the ice files. Uniform steps: in Fortran the step is set at
    build time (`DT=<years>`, default 1; `NSTEPS` is the number of steps); in Python with
    `--dt <years>` on `run_ais.py` and `run_outer.py`. Non-uniform steps: a times file with the time
    of every ice file ([format](inputs.md#times-file-optional-non-uniform-time-steps)), found by
    Fortran at `INPUT_FILES/others/times` (or `SLEVOLVE_TIMES=<path>`) and passed to Python with
    `--times <file>`. Both solvers accept any uniform step and any non-uniform schedule.

## 1. Ice history

* **Grid.** Every ice file is a 512 x 1024 grid on the Gauss-Legendre grid of the solver: row `j` =
  latitude `GLlat[j]` from north to south, column `i` = longitude `i * 360/1024` degrees east
  ([grid files](#4-grid-files)). Store it as `.npy` (`(512, 1024)`, `<f8`), as a headerless
  little-endian float64 `.bin`, or as text in the stream order (`python/tools/text2npy.py`
  converts).
* **Units.** Ice thickness in metres of ice (ice density `rhoi` = 920 kg m^-3, `--rhoi` to change).
* **Grounded or total ice.** With the marine check on (`--checkmarine T`, the default) the solver keeps
  only ice that is grounded: ice of thickness H on a cell whose initial topography T0 (step 0) is not
  positive is kept if `rhoi * H > rhow * |T0|` (`rhow` = 1000 kg m^-3), otherwise it is treated as
  floating and removed from the load. You can then give the total (grounded plus floating) thickness. With the marine check off,
  give grounded ice only: floating ice would be counted as load.
* **File naming and time.** Files `<prefix>0`, `<prefix>1`, ..., `<prefix>N`, one per step. File 0 is
  the initial state (the oldest time). With uniform spacing `dt`, file k is the ice at time `k * dt`
  after the start; with a times file, file k is the ice at time `t_k` of line k + 1 of the file
  ([section 5](#5-step-length-and-truncation)). For an ice age run the last file N is the present.
* **Regridding.** `tools/regrid.py` provides `NearestSourceSampler(alat, alon, lon_shift=0,
  maxdist_km=...)`: every solver cell within `maxdist_km` of a source point takes the value of the
  nearest source point, all other cells are 0; the source grid is given by 2-D latitude and longitude
  arrays in degrees (any grid: lat-lon, polar stereographic, ...). `lon_shift=0` samples onto the
  solver longitudes `i * 360/1024`; the default `lon_shift=1` is the one-cell offset used only for the
  Case 1 files. Nearest-point sampling preserves thickness values but not the volume exactly, so
  compare the volumes. A runnable example (from the release root; numpy and scipy), with a 0.5-degree
  lat-lon source field made up in the script:

```python
import sys
import numpy as np
sys.path.insert(0, "tools")                       # run from the release root
from regrid import NearestSourceSampler, gl_grid

# source: a regular 0.5-degree lat-lon grid (replace by your own arrays, e.g. read with netCDF4:
#   ds = netCDF4.Dataset("my_ice.nc"); lat1 = ds["lat"][:]; lon1 = ds["lon"][:]; H = ds["thk"][k])
lat1 = np.arange(-89.75, 90.0, 0.5)
lon1 = np.arange(0.25, 360.0, 0.5)
alon, alat = np.meshgrid(lon1, lat1)
H = np.where(alat < -65.0, 2000.0 * np.sqrt(np.clip(1 - ((alat + 90) / 25) ** 2, 0, None)), 0.0)

# nearest-source sampling onto the 512 x 1024 solver grid (lon = i * 360/1024: lon_shift=0);
# maxdist_km must exceed half the source spacing (0.5 deg ~ 56 km at the equator)
smp = NearestSourceSampler(alat, alon, lon_shift=0, maxdist_km=60.0)
Hg = smp.sample(H)
np.save("myice_0.npy", Hg.astype("<f8"))

# volume check: source cell areas on the sphere vs solver cell areas (gl_grid)
R = 6.371e6
a_src = R**2 * np.radians(0.5) * (np.sin(np.radians(alat + 0.25)) - np.sin(np.radians(alat - 0.25)))
_, _, _, _, _, area = gl_grid()
v_src, v_gl = (H * a_src).sum(), (Hg * area).sum()
print(f"ice volume: source {v_src:.4e} m^3, solver grid {v_gl:.4e} m^3, ratio {v_gl / v_src:.5f}")
```

  It prints `ice volume: source 3.1985e+16 m^3, solver grid 3.2069e+16 m^3, ratio 1.00263`. Repeat
  the sampling for every ice file (`myice_0.npy`, `myice_1.npy`, ...).

## 2. Topography

* **What it is.** Bedrock (solid-surface) elevation relative to sea level, in metres, positive up, on
  the same grid. It is the elevation of the rock, not of the ice surface: the ice thickness is given
  separately. The ocean is where the topography is negative.
* **Regridding.** A global topography on a lat-lon grid (e.g. ETOPO or GEBCO, averaged to a few
  arc-minutes first) is sampled onto the solver grid with the same `NearestSourceSampler` as the ice
  (a global source covers every cell when `maxdist_km` exceeds half the source spacing). Save it as a
  `(512, 1024)` `<f8` `.npy` file.
* **Forward run** (known initial state): the topography at the start, step 0 (`run_ais.py --topo`;
  Fortran `TOPO_INITIAL`).
* **Ice age run** (initial state unknown): the observed present topography (`run_outer.py
  --topo-pres`; Fortran `TOPOMODEL`, with `INITIAL_TOPO=F`). The outer iteration adjusts the initial
  topography until the topography predicted at the last step matches it.

## 3. Love numbers

SLEvolve reads precomputed Love numbers; it does not compute them. Use one of the shipped files in
`lovenum/` or write your own Earth model's Love numbers with `python/tools/write_lovenum.py`; the
layout is in [Love numbers](lovenum.md). The shipped files (all PREM elastic and density structure,
degrees 1-512):

| file | lithosphere | upper mantle viscosity | lower mantle viscosity | used in |
|---|---|---|---|---|
| `prem_512.l60K2C.sum18p6.dum19p2.tz19p4.lm22` | 60 km | 10^18.6 Pa s (60-200 km), 10^19.2 Pa s (200-400 km), 10^19.4 Pa s (400-670 km) | 10^22 Pa s (from the file name; see [Benchmark provenance](provenance.md)) | Case 1 (Han et al. 2022 file, unchanged) |
| `prem_512.l120C.ump5.lm5` | 120 km | 5e20 Pa s | 5e21 Pa s | Case 2 (collocation) |
| `prem_512.l96C.ump3.lm1` | 96 km | 3e20 Pa s | 1e21 Pa s | Case 3, the quick start and the Fortran test (collocation) |

## 4. Grid files

`GLlat_512.txt` and `GLlon_512.txt` list the latitudes and longitudes of the grid (degrees; latitude
north to south). Neither solver reads them; they document the grid and are used by the input
preparation scripts. Write them with

```bash
python tools/make_grid.py --outdir INPUT_FILES/grid
```

Details: [Input files](inputs.md#grid-files).

## 5. Step length and truncation

* **Uniform dt.** One fixed step length per run, equal to the spacing of the ice files. Python:
  `run_ais.py --dt` and `run_outer.py --dt`, both in **years**, and `init_topo_segmented.py --dt` in
  years. Each driver prints the step length it uses. Fortran: compiled in, `build.sh DT=<years>`
  (integer). The Fortran forward driver `run_from_zero_plus` uses the build-time `DT` (verified against
  Python with `DT=1` and `DT=5`); for non-uniform steps use a times file (next item).
* **Non-uniform steps (times file).** Write the time of every ice file, in years, one per line
  (`t_0 .. t_N`, strictly increasing; [format](inputs.md#times-file-optional-non-uniform-time-steps)).
  Python: `--times <file>` on `run_ais.py`, `run_outer.py` and `init_topo_segmented.py`. Fortran: put
  the file at `INPUT_FILES/others/times` (detected at start-up) or set `SLEVOLVE_TIMES=<path>`; build
  with `DT` unset, so that `NSTEPS` is the number of steps and the times file sets their lengths.
  Step j has length `dt_j = t_j - t_(j-1)`, and both solvers use it in
    * the modal-state update, with decay factors `d_k = exp(-s_k dt_j)` recomputed at every change of
      step length (the recursion stays exact for the step-wise load history of the files);
    * the extrapolated inner start (`--inner extrap`, `SLEVOLVE_INNER=extrap`), which from step 3 on is
      `dS0_j = dt_j [Sdot_(j-1) + (Sdot_(j-1) - Sdot_(j-2)) (dt_j + dt_(j-1)) / (dt_(j-1) + dt_(j-2))]`
      with `Sdot_i = dS_i / dt_i` (the paper, Sect. 2.5), equal to `2 dS_(j-1) - dS_(j-2)` for equal
      steps; steps 1 and 2 start from the eustatic guess;
    * the time stamps (Fortran `OUTPUT_SLM/times` holds `t_0 .. t_N`; Python `run_log.json` /
      `outer_log.json` hold `times_years`).
  A uniform times file gives the same result as the run without one exactly (zero difference).
* **dt and accuracy.** The recursion is exact for the step-wise load history the files define; the
  load is constant between files. A shorter dt (or shorter steps where the load changes fast, with
  a times file) resolves the history better, and also the shoreline
  migration, which is updated once per step.
* **lmax.** `--lmax 511` on the 512 x 1024 grid, as in all benchmark cases (the ducc0 transform needs
  lmax <= 511 for 512 latitudes). Fortran: `NORDER=511`; the grid size is fixed at 512.

## 6. Running

**Forward run** (known initial topography):

```bash
python python/run_ais.py --steps N --dt DT --threads 4 --sht ducc --lmax 511 --opt all --inner extrap \
    --icedir INPUT_FILES/icemodel --icepfx <prefix> --topo INPUT_FILES/others/<topography> \
    --earth lovenum/<love file> --save-fields tgrid,ocean --out OUT
```

**Fortran forward run.** The input file names are compiled into the Fortran executable, so each new
case is a build without `CASE` (`CASE` only selects the presets of the three benchmark cases; leave it
unset) and with your file names:

```bash
ICEMODEL=<prefix> TOPOMODEL=<topography> PLANETMODEL=<love file> NSTEPS=N bash fortran/build.sh
# -> fortran/build/custom_N<N>/exe   (DT=1, NORDER=511, SHT=otf, OpenMP: the defaults)
```

Then set up a run directory: the executable reads `INPUT_FILES/` and writes `OUTPUT_SLM/` and
`ICELOAD_SLM/` (both must exist) relative to the directory it runs in.

```bash
mkdir -p RUN/OUTPUT_SLM RUN/ICELOAD_SLM
ln -s /path/to/INPUT_FILES RUN/INPUT_FILES      # icemodel/<prefix>0..N, others/<topography>, earthmodel/<love file>
cp fortran/build/custom_N<N>/exe RUN/ && cd RUN
OMP_NUM_THREADS=4 OMP_STACKSIZE=32M SLEVOLVE_INNER=extrap SLEVOLVE_SNAP=N ./exe N
```

`SLEVOLVE_SNAP=N` writes the grids of step N only (and step 0); omit it to write every step. With
`NSTEPS=550` the directory is `custom` (550 is the default). `tests/test_fortran_python.sh` is a
complete example of this route.

**Ice age run** (outer iteration on the initial topography):

```bash
python python/init_topo_segmented.py --inputs INPUT_FILES --icepfx <prefix> --nfiles <N+1> --dt DT \
    --earth lovenum/<love file> --lmax 511 --topo-pres INPUT_FILES/others/<present topography> \
    --out T0_segmented.npy
python python/run_outer.py --iters 5 --inner extrap --reuse warm --lmax 511 --threads 4 \
    --inputs INPUT_FILES --icepfx <prefix> --nfiles <N+1> --dt DT \
    --topo-pres INPUT_FILES/others/<present topography> --earth lovenum/<love file> \
    --topo0 T0_segmented.npy --out OUT
```

`run_outer.py` defaults (`--icepfx nadi_ --nfiles 51 --dt 500`) are those of Case 3: always pass your
own. Example: an ice history from 21 ka to the present every 500 yr has 43 files (`<prefix>0` at
21 ka, `<prefix>42` at present), so `--nfiles 43 --dt 500` here, and `--steps 42 --dt 500` for a
forward run with `run_ais.py`.

Fortran: build with `INITIAL_TOPO=F`, `DT=<years>`, `NSTEPS=N`, your file names and
`DRIVER=run_outer_mem`, and run the passes with `examples/run_outer_mem.sh` as for Case 3, with
`TOPO_PRES=<your present topography>` (text or `.npy`; the default is `INPUT_FILES/others/topo_pres`)
and your first-pass initial topography as `TOPO0` (`.npy`; or `pres` to start from the present topography)
([Running the benchmark cases](benchmarks.md#case-3-ice-nadi-1-50-05-kyr-full-outer-iteration)). The Python driver prints the outer misfit `rms_pred_minus_obs`, the area-weighted RMS of predicted minus
observed present topography (m), after each pass and records it in `outer_log.json`; the Fortran driver
writes `pred_pres_topo_<p>.npy` and `tgrid0_<p>.npy` per pass, from which the same misfit is one
`numpy` line (`np.sqrt(np.average((pred - obs)**2, weights=gw[:, None]))` with the Gauss weights `gw`).
Run passes until it stops decreasing by a meaningful amount for your application; in Case 3 it falls from 12.0 m after pass 1 to 2.06 m after pass 2
(`reference/case3.json`, `misfit_rms`).

## 7. Checking the result

* `python tests/run_all.py` and `bash tests/test_fortran_python.sh` check the installation
  ([Verifying the installation](verify.md)).
* The benchmark cases and `tests/check_reference.py` check a full-size configuration against
  `reference/` ([Verifying the installation](verify.md#benchmark-reference-check)).
* For your own case, look at:
  * `numiter` (inner iterations per step; Case 1 needs 2099 in 550 steps): compare a step that needs
    many more than its neighbours with the change of `ocean<N>` at that step;
  * the ice volume (`run_log.json` `scalars`, Fortran `ice_volume`) against the volume of your input
    files, which shows what the marine check removed;
  * the GMSL change (`gmslc_ocn`) against the barystatic estimate from the ice-volume change;
  * the maps `tgrid<N> - tgrid0` (minus the relative sea-level change) and `ocean<N>` (the shoreline),
    e.g. `python -c "import numpy as np; d = np.load('OUT/tgrid42.npy') - np.load('OUT/tgrid0.npy'); print(d.min(), d.max())"`;
    `python/tools/npy2nc.py` converts the outputs to NetCDF for plotting.
* For an ice age run: `outer_log.json` (misfit per pass) and `tgrid0_<p>`, the initial topography
  of the last pass.

**Next:** [Love numbers](lovenum.md).
