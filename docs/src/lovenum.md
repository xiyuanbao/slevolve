# Love numbers

This page explains how to supply Love numbers for your own Earth model and how to write them with `python/tools/write_lovenum.py`.

SLEvolve reads precomputed Love numbers; it does not compute them. For another Earth model, compute
the load and tidal Love numbers in normal-mode form (elastic part, fluid limit, and per mode an
inverse relaxation time and residues for h, l and k) with an external normal-mode or collocation code
(for example your own), and write them with `python/tools/write_lovenum.py`. The three shipped files
are the only tested inputs:

| file | lithosphere | upper mantle viscosity | lower mantle viscosity | used in |
|---|---|---|---|---|
| `prem_512.l60K2C.sum18p6.dum19p2.tz19p4.lm22` | 60 km | 10^18.6 Pa s (60-200 km), 10^19.2 Pa s (200-400 km), 10^19.4 Pa s (400-670 km) | 10^22 Pa s (from the file name; see [Benchmark provenance](provenance.md)) | Case 1 (Han et al. 2022 file, unchanged) |
| `prem_512.l120C.ump5.lm5` | 120 km | 5e20 Pa s | 5e21 Pa s | Case 2 (collocation) |
| `prem_512.l96C.ump3.lm1` | 96 km | 3e20 Pa s | 1e21 Pa s | Case 3, the quick start and the Fortran test (collocation) |

All three use the PREM elastic and density structure and hold degrees 1-512. The file names encode
the structure: `l<km>` lithosphere thickness, `ump<x>` upper mantle viscosity x 10^20 Pa s, `lm<x>`
lower mantle viscosity x 10^21 Pa s (`lm22` in the Case 1 name: 10^22 Pa s).

Layout read by both solvers (`solver.LoveNumbers`; Fortran block "Read in Love numbers"). Each group
starts on a new line and may continue over several lines (Fortran list-directed reads):

| line(s) | content |
|---|---|
| 1 | header, free text, skipped |
| per degree l = 1, 2, ..., L (consecutive, from 1; no degree 0) | |
| a | `l  K_l`: degree and number of modes; K_l may differ between degrees, at most 500 |
| b | `s_1 .. s_K`: inverse relaxation times, 1/kyr, positive |
| c | `h^E  l*l^E  x  x  l*k^E`: elastic loading Love numbers; the two middle values are read and ignored (0 in the shipped files) |
| d | `h^F  l*l^F  x  x  l*k^F`: fluid limit, same layout; read, not used |
| e, f, g | `K_l` residues each, 1/kyr: h; l multiplied by l; k multiplied by l |
| only for l >= 2: h to l | the same five groups (c to g) for the tidal Love numbers |

* h is stored as is; l and k (elastic, fluid and residues) are multiplied by the degree l, and the
  readers divide by l.
* The response to a unit load applied at t = 0 is `h(t) = h^E + sum_k (r_k / s_k) (1 - exp(-s_k t))`
  with t in kyr (the code uses `exp(-s_k dt / 1000)` with dt in years); the same for l and k.
* The solver reads the first `lmax` degree blocks; the file may hold more.
* Degree 1 must be in the frame you want for the degree-1 outputs: SLEvolve applies no frame
  conversion. The shipped files are in the CE frame, k^E_1 = 0 ([Degree-1 Love numbers and reference frames](degree1.md)).
* The tidal numbers enter only with rotation on (`tpw`), at degree 2.

`write_lovenum.py` provides `write_lovenum(path, header, hE, lE, kE, hF, lF, kF, s, rh, rl, rk, ...)`
(arrays per degree, unscaled; the docstring has a complete example) and
`write_from_lovenumbers(path, love, header)`. The command line form rewrites a file through the
reader:

```bash
python python/tools/write_lovenum.py --from lovenum/prem_512.l96C.ump3.lm1 --lmax 512 --out copy.txt
```

For all three shipped files this reproduces every number of the original to a relative 4.2e-16
(maximum absolute difference 1.1e-13).

The line-by-line read of the file by the solvers is specified in [Input files](inputs.md#love-number-earth-model-file). Degree-1 frames and the collocation initialisation defect of the Case 1 file are discussed in [Degree-1 Love numbers and reference frames](degree1.md).

**Next:** [Horizontal displacement output](horizontal.md).
