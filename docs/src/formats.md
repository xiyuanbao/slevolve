# File formats

This page specifies the binary containers of the grid and spectral files.

Both solvers use **NumPy `.npy`** as their binary format. NetCDF is only an export
format, made from Python. Text in the list-directed grid format of the legacy code of Han et al. (2022) is handled by converters.

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

`<c16`, shape `(L+1, L+1)` with L = `NORDER`, **fortran_order=True**, so that `a[l, m] == dS(l, m)` (zero for `m > l`).
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
  [Horizontal displacement](horizontal.md).

## 3. File naming and switches (Fortran)

* Snapshots: `OUTPUT_SLM/<field><N>.npy` with field in `tgrid, ocean, beta, dS_converged, R, G` (optional: `horiz`, `Psi`, §2.1),
  the names of the legacy `.bin` files with `.npy` in place of `.bin`. `SLEVOLVE_SNAP` /
  `SLEVOLVE_SNAP_FIELDS` select what is written. `snap_index.txt` lists `step field file dtype n1 n2`.
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
`--topo0` through `python/slevolve/npyio.py` (`.npy` > `.bin` > text).

## 4. f4 losslessness of the ice inputs

* AIS grounded ice `AISload_gd_Run85_0_0_{0..550}`: **551/551 files exactly representable in
  float32**. `<f4` caches are lossless; a 20-step run from them gives identical output ([Release notes](release_notes.md#checks-made-during-development)).
* nadi ice `nadi_{0..50}` (printed `%18.10E` from float32 values): **0/51 exact**. Only 81-88 % of
  cells round-trip, because 11 significant digits of a float32 value do not parse back to a float32.
  These are stored as `<f8`.
* archived Case 1 topography (8 significant digits): not exact (19 % of cells) -> `<f8`.

`text2npy --f4` checks this per file and falls back to `<f8` with a warning. It never loses data.

The converters are described in [File conversion](conversion.md); development checks of the formats in [Release notes](release_notes.md#checks-made-during-development).

## 5. Portability

`fortran/src/npy_io_mod.f90` is self-contained: it depends only on `iso_fortran_env`, not on any
solver module. Its minimal interface:

* `read_grid(a, path, ios[, source, errmsg])`: `.npy` > `.bin` > text
* `read_ice_cache(a, path, ios, ...)`: same lookup
* `write_grid(a, path, ios[, dtype='f4'|'f8'])`: f4 is used only when lossless
* `write_snapshot(a|z, path, ios)`: real grid, or complex spectral field (fortran_order)

It also has a low-level layer (`npy_read_header`, `npy_read_r8`, `npy_read_c16`, `npy_write_r8/c16/i1`,
`npy_f4_exact`). The solver-side hook is the solver I/O module in `fortran/src/`, which uses only
`read_grid`, `write_grid` and `write_snapshot`.

**Next:** [Code structure](structure.md).
