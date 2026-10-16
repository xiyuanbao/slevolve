"""Build the TOTAL ice forcing (grounded + floating shelves) of Case 1, which the archive of
Han et al. (2022) does not contain, on their 512x1024 Gauss-Legendre grid from Yousefi et al. (2022) fort.92_85_MISI.

Two variants are written:

  total   : the archived grounded field of Han et al. (2022) where a grounded cell was sampled,
            plus the sampled floating thickness everywhere else.  This keeps
            97.6 % of the mass identical to the archived input and only the
            shelves carry the 1 % residual of the recovered sampling map.
  totalY  : a pure resampling of the Yousefi total thickness (control).

Output: flat little-endian f8 (512,1024) '<name><k>.bin', the solver's binary cache
format, under <out>/<variant>/ (see --help).
"""
import os
import sys
import time

import netCDF4 as nc
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from regrid import NearestSourceSampler, gl_grid  # noqa: E402

import argparse  # noqa: E402
_ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
_ap.add_argument("--nc", required=True,
                 help="Yousefi et al. (2022) fort.92_85_MISI.nc (Dryad doi:10.5061/dryad.ncjsxkswm)")
_ap.add_argument("--grounded", required=True,
                 help="directory with the archived grounded-ice files AISload_gd_Run85_0_0_<k>.bin "
                      "(Han et al. 2022 OSF archive, converted with fortran tools/make_cache)")
_ap.add_argument("--out", required=True, help="output root; writes <out>/{total,totalY}/icemodel/")
_A = _ap.parse_args()
ROOT = _A.out
NC = _A.nc
GDBIN = os.path.join(_A.grounded, "AISload_gd_Run85_0_0_%d.bin")
PFX = "AISload_gd_Run85_0_0_"          # keep the name the support module expects

OUT = {v: f"{ROOT}/{v}/icemodel" for v in ("total", "totalY")}
for p in OUT.values():
    os.makedirs(p, exist_ok=True)

d = nc.Dataset(NC)
S = NearestSourceSampler(d["alatd"][:].data, d["alond"][:].data)
lat, lon, mu, w, edges, A = gl_grid()

log = open(f"{ROOT}/forcing_volumes.csv", "w")
log.write("k,vol_gd_archive,vol_total,vol_totalY,vol_shelf,toti_native,totig_native,totif_native\n")
toti = d["toti"][:].data
totig = d["totig"][:].data
totif = d["totif"][:].data

t0 = time.time()
for k in range(551):
    h = d["h"][k].data.astype(np.float64)
    gd = d["maskwater"][k].data == 0
    gH = np.fromfile(GDBIN % k, dtype="<f8").reshape(512, 1024)
    fl = S.sample(np.where(gd, 0.0, h))
    tot = gH + np.where(gH > 0.0, 0.0, fl)
    totY = S.sample(h)
    tot.astype("<f8").tofile(f"{OUT['total']}/{PFX}{k}.bin")
    totY.astype("<f8").tofile(f"{OUT['totalY']}/{PFX}{k}.bin")
    log.write("%d,%.8e,%.8e,%.8e,%.8e,%.8e,%.8e,%.8e\n"
              % (k, (gH * A).sum(), (tot * A).sum(), (totY * A).sum(),
                 ((tot - gH) * A).sum(), toti[k], totig[k], totif[k]))
    if k % 50 == 0:
        print("k=%3d  %.1fs  vol_total %.5e  shelf %.4e"
              % (k, time.time() - t0, (tot * A).sum(), ((tot - gH) * A).sum()), flush=True)
log.close()
print("done %.1f s" % (time.time() - t0))
