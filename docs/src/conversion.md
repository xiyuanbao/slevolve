# File conversion

This page lists the converters between the legacy text format, `.npy`, the legacy `.bin` cache and NetCDF, in Fortran and Python.

| task | Fortran (`fortran/tools/`; built by `fortran/build.sh` into `fortran/build/<case>/` (e.g. `fortran/build/case1/`), or by `fortran/tools/build_npytools.sh` into `fortran/tools/bin/`) | Python (`python/tools/`) |
|---|---|---|
| text -> npy | `text2npy [--f4] [--nlat N --nlon M] [--outdir D] FILE...` or `--dir D --prefix P --range K0 K1` | `text2npy.py [--f4] [--nlat N] [--nlon M] [--outdir D] FILE...` or `--dir D --pattern 'P*'` |
| npy -> text | `npy2text [--format F] [--ncol 1024] [--outdir D] FILE.npy...` | `npy2text.py [--format F] [--ncol N] [--outdir D] FILE...` |
| text -> LE f8 `.bin` (legacy cache) | `make_cache NGLV OUTDIR FILE [FILE...]` | none |
| npy snapshots -> NetCDF4 | none | `npy2nc.py RUNDIR --out run.nc [--fields tgrid,ocean] [--times F / --dt Y] [--earth E] [--solver S] [--lon-convention solver\|sampled] [--attr k=v] [--zlib]` |

Text formats (`--format`):

| name | layout | used by |
|---|---|---|
| `es16.9` (default) | `%16.9E` (Fortran `ES16.9E2`), one value per line | the solver's own text output (`write_sl`) |
| `row-e16.7` | `%16.7e`, 1024 per line, no separator | OSF archive topography / tgrid files |
| `row-e.18` | `%.18e`, 1024 per line, `' '` separator | OSF Case 1 ice files (`np.savetxt` default) |
| `row-E18.10` | `%18.10E`, 1024 per line, `' '` separator | Case 3 (ICE-NADI-1) ice files written by `tools/prepare_case3_nadi.py` |

`npy2nc.py` coordinates:
* `lat` = asin of the solver's Gauss-Legendre nodes (`sht.gaulw`, the same as `GLlat_512.txt` to 5e-7 deg), north to south.
* `lon` = `i*360/1024`, the solver's SHT grid and `GLlon_512.txt`.
  `--lon-convention sampled` gives the `(i+1)*360/1024` mapping, which is the mapping the archived
  Case 1 ice was *sampled* with (see `tools/regrid.py`). It is not the solver grid.
* `time` comes from `OUTPUT_SLM/times` when present, otherwise from `step*dt`.
* The global attributes carry the solver, the Earth model and the planet constants (`solver.Planet`,
  with `run_log.json` rhoi/radius overrides). For Fortran builds with `LEGACY_CONSTANTS=1`, pass `--attr` overrides.

The `.npy` layout is specified in [File formats](formats.md); the `--f4` rule in [File formats](formats.md#4-f4-losslessness-of-the-ice-inputs).

**Next:** [Reproducing the published results](reproduce.md).
