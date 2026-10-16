# I/O formats: Fortran (fortran/src) and Python (python/) solvers

Both solvers use **NumPy `.npy`** as their binary format. NetCDF is only an export
format, made from Python. Text in the legacy list-directed grid format (Kendall et al.
2005 code lineage) is handled by converters.

## 1. Grids (ice, topography, tgrid/ocean/beta/R/G snapshots)

| property | value |
|---|---|
| container | `.npy` version 1.0 (2.0 if the header exceeds 64 KiB), 64-byte-aligned header |
| shape | `(nlat, nlon)` = `(512, 1024)` |
| fortran_order | `False` |
| dtype | `<f8`. Ice inputs may be `<f4`, but only when they are lossless (see §4). Masks may be `|i1` |
| element order | the stream order of the text files, i.e. `a[j, i]`: `j` = Gauss-Legendre row (north to south), `i` = longitude column |
| payload | identical to the legacy headerless `.bin` cache (`numpy.fromfile(f, "<f8")`) |

Header written by both codes (identical):

```
\x93NUMPY \x01\x00 v\x00 {'descr': '<f8', 'fortran_order': False, 'shape': (512, 1024), }<spaces>\n
```

A `<f4` ice cache has the same header with `'descr': '<f4'`. The payload starts at byte 128.

In Fortran, the payload is read straight into `data_slm(nglv, 2*nglv)`. This fills the
array exactly as `read(u,*) data_slm` fills it from the text file. The array's storage
order is the stream order (the SHT treats longitude as fastest). In numpy,
`np.load(f)` gives `a` with `a.ravel() ==` the text stream.

## 2. Spectral fields (`dS_converged<N>.npy`)

`<c16`, shape `(N+1, N+1)`, **fortran_order=True**, so that `a[l, m] == dS(l, m)` (zero for `m > l`).
The payload is again identical to the legacy `dS_converged<N>.bin`. The Python solver writes the same
thing (`npyio.save_spectral`) when `--io-compat fortran` is set (the default).

### 2.1 Optional horizontal-displacement outputs

* `horiz<N>.npy` (Fortran `SLEVOLVE_HORIZ=1|2`; Python `run_ais.py --horizontal`; `run_outer.py
  --horizontal` writes `horiz<N>_it<k>.npy`): the spectral scalar Psi_lm, same container and
  convention as `dS_converged<N>.npy` (`<c16`, `(L+1, L+1)`, fortran_order=True, spharmt U
  convention, `a[l, m]`, zero for `m > l` and for `l = 0`). Units: metres.
* `Psi<N>.npy` (Fortran `SLEVOLVE_HORIZ=2`): Psi synthesised on the grid, grid layout of §1.
* `utheta<N>.npy`, `uphi<N>.npy` (Python `--horizontal-grid`, or `python/tools/horiz_grad.py
  horiz<N>.npy` for Fortran output): grid layout of §1, `<f8`, metres.
  u_theta = dPsi/dtheta (theta = colatitude, positive southward), u_phi = (1/sin theta) dPsi/dphi
  (positive eastward); u_north = -u_theta.
* In snapshot mode (`SLEVOLVE_SNAP`) the Fortran fields are named `horiz` and `Psi` in
  `SLEVOLVE_SNAP_FIELDS` (`all` includes them). Convention and verification:
  the Horizontal displacement page of the manual (`docs/src/horizontal.md`).

## 3. File naming and switches

* Snapshots: `OUTPUT_SLM/<field><N>.npy` with field in `tgrid, ocean, beta, dS_converged, R, G` (optional: `horiz`, `Psi`, §2.1),
  the same names as before with `.npy` in place of `.bin`. `SLEVOLVE_SNAP` / `SLEVOLVE_SNAP_FIELDS`
  work as before. `snap_index.txt` lists `step field file dtype n1 n2`.
* `SLEVOLVE_FMT=npy` (default) | `bin` selects the output container (`bin` = the old headerless stream).
* Inputs: for every grid `<path>`, the solver tries `<path>.npy`, then `<path>.bin`, then the text
  file `<path>`. This covers ice caches, the initial topography and `pred_pres_topo_/tgrid0_`.
  The `.npy` files may be `<f8`, `<f4`, `|i1`, `<i4` or `<i8`, in either order. f4 is promoted exactly.
* Not converted, on purpose: the exact-memory checkpoint (`exact_memory_state_<i>.bin`,
  `SLEVOLVE_CKPT`, off by default) and the file-mode warm-reuse store (`reuse_<i>.bin`,
  `SLEVOLVE_REUSE=warm` without in-memory outer iterations). Each is a multi-array
  record file with its own magic, shape check and random-access record layout, and each is read only by the
  solver that wrote it. A single `.npy` cannot hold them; the `.npz` equivalent would need zip I/O.

The Python solver (`run_ais.py`, `run_outer.py`) writes `tgrid<N>.npy`, `R`, `G` in the same layout.
With `--io-compat fortran` (the default) it writes `ocean<N>.npy` / `beta<N>.npy` as `<f8`
(0/1) like the Fortran code; `--io-compat legacy` restores uint8. It reads ice, topography and
`--topo0` through `python/npyio.py` (`.npy` > `.bin` > text).

## 4. f4 losslessness of the ice inputs

* AIS grounded ice `AISload_gd_Run85_0_0_{0..550}`: **551/551 files exactly representable in
  float32**. `<f4` caches are lossless; the 20-step run from them is identical (§6).
* nadi ice `nadi_{0..50}` (printed `%18.10E` from float32 values): **0/51 exact**. Only 81-88 % of
  cells round-trip, because 11 significant digits of a float32 value do not parse back to a float32.
  These are stored as `<f8`.
* archive `std_dt1/topo0` (8 significant digits): not exact (19 % of cells) -> `<f8`.

`text2npy --f4` checks this per file and falls back to `<f8` with a warning. It never loses data.

## 5. Converters

| task | Fortran (`fortran/tools`, build with `fortran/tools/build_npytools.sh` -> `fortran/tools/bin/`) | Python (`python/tools`) |
|---|---|---|
| text -> npy | `text2npy [--f4] [--nlat N --nlon M] [--outdir D] FILE...` or `--dir D --prefix P --range K0 K1` | `text2npy.py [--f4] [--outdir D] FILE...` or `--dir D --pattern 'P*'` |
| npy -> text | `npy2text [--format F] [--ncol 1024] [--outdir D] FILE.npy...` | `npy2text.py` (same options) |
| npy snapshots -> NetCDF4 | none | `npy2nc.py RUNDIR --out run.nc [--fields tgrid,ocean] [--times F / --dt Y] [--earth E] [--solver S] [--attr k=v] [--zlib]` |

Text formats (`--format`):

| name | layout | used by |
|---|---|---|
| `es16.9` (default) | `%16.9E` (Fortran `ES16.9E2`), one value per line | the solver's own text output (`write_sl`) |
| `row-e16.7` | `%16.7e`, 1024 per line, no separator | OSF archive topography / tgrid files |
| `row-e.18` | `%.18e`, 1024 per line, `' '` separator | OSF AIS ice files (`np.savetxt` default) |
| `row-E18.10` | `%18.10E`, 1024 per line, `' '` separator | nadi-case ice files |

`npy2nc.py` coordinates:
* `lat` = asin of the solver's Gauss-Legendre nodes (`sht.gaulw`, the same as `GLlat_512.txt` to 5e-7 deg), north to south.
* `lon` = `i*360/1024`, the solver's SHT grid and `GLlon_512.txt`.
  `--lon-convention sampled` gives the `(i+1)*360/1024` mapping, which is the mapping the archived
  ice was *sampled* with (`tools/regrid.py`). It is not the solver grid.
* `time` comes from `OUTPUT_SLM/times` when present, otherwise from `step*dt`.
* The global attributes carry the solver, the Earth model and the planet constants (`solver.Planet`,
  with `run_log.json` rhoi/radius overrides). For Fortran builds with `LEGACY_CONSTANTS=1`, pass `--attr` overrides.

## 6. Verification

* Converters, on the OSF AIS ice file `AISload_gd_Run85_0_0_0` (`row-e.18`), the archived
  topography `std_dt1/topo0` (`row-e16.7`) and the ICE-NADI-1 ice file `icemodel/nadi_0`
  (`row-E18.10`):
  * Fortran- and Python-written `.npy` files are identical.
  * `np.load` of the Fortran file equals the text parse.
  * text -> npy -> text reproduces the original file byte for byte in its own format, in both languages.
  * `es16.9` output is identical between Fortran and Python and round-trips exactly.
  * The Fortran reader also handles fortran_order=True and `|i1` files.
* Solver: the `.npy` build against a build reading the headerless `.bin` caches,
  SHT=blas, OMP=4, 20 AIS steps, `SLEVOLVE_SNAP=all SLEVOLVE_SNAP_FIELDS=all`. The `.npy` build read only
  `<f4` `.npy` ice and a `<f8` `.npy` topography; the other read the `.bin` caches.
  * All 124 snapshot files match: the `.npy` payload is identical to the `.bin`, and the values are equal.
  * All scalar text outputs are identical (except CPU time).
* Python solver, 20 steps, lmax 511: reading the Fortran-written `<f4` ice `.npy` and `.npy` topography
  gives outputs identical to reading the text files (84/84 files). Its `tgrid20.npy`,
  `ocean20.npy` and `dS_converged20.npy` headers are identical to the Fortran files'.

## 7. Portability

`fortran/src/npy_io_mod.f90` is self-contained: it depends only on `iso_fortran_env`, not on any
solver module. Its minimal interface:

* `read_grid(a, path, ios[, source, errmsg])`: `.npy` > `.bin` > text
* `read_ice_cache(a, path, ios, ...)`: same lookup
* `write_grid(a, path, ios[, dtype='f4'|'f8'])`: f4 is used only when lossless
* `write_snapshot(a|z, path, ios)`: real grid, or complex spectral field (fortran_order)

It also has a low-level layer (`npy_read_header`, `npy_read_r8`, `npy_read_c16`, `npy_write_r8/c16/i1`,
`npy_f4_exact`). The solver-side hook is `io_mod.f90`, which uses only
`read_grid`, `write_grid` and `write_snapshot`.
