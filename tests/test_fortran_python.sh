#!/usr/bin/env bash
# Fortran test: build the Fortran solver, run it and the Python solver on the synthetic Earth
# (tools/make_synthetic_case.py, 512 x 1024 grid, 8 steps of 1 yr) and compare the final topography.
#
#   bash tests/test_fortran_python.sh [WORKDIR]        (default WORKDIR: ./slevolve_fortran_test)
#
# Needs the Fortran toolchain active (see the Installation page; or set CONDA_ENV=<env>) and Python
# with numpy (ducc0 and numba are used when installed).  PYTHON selects the interpreter (default
# python3), THREADS the thread count (default 4).  Prints PASS when max |tgrid8_Fortran - tgrid8_Python|
# < 1e-9 m and exits 0; otherwise prints FAIL and exits 1.  Takes about half a minute.
set -eo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
WORK="${1:-$PWD/slevolve_fortran_test}"; PYTHON="${PYTHON:-python3}"; THREADS="${THREADS:-4}"
mkdir -p "$WORK"; WORK="$(cd "$WORK" && pwd)"
IN="$WORK/INPUT_FILES"
echo "== synthetic case -> $IN"
"$PYTHON" "$HERE/tools/make_synthetic_case.py" --outdir "$IN"

echo "== Fortran build (CASE unset, NSTEPS=8, defaults)"
NSTEPS=8 ICEMODEL=synth_ice_ TOPOMODEL=synth_topo PLANETMODEL=prem_512.l96C.ump3.lm1 TAG=fortran_test \
  bash "$HERE/fortran/build.sh" > "$WORK/build.log" 2>&1 || { tail -20 "$WORK/build.log"; echo "FAIL (build)"; exit 1; }
tail -1 "$WORK/build.log"

echo "== Fortran run"
RUN="$WORK/RUN_fortran"; rm -rf "$RUN"; mkdir -p "$RUN/OUTPUT_SLM" "$RUN/ICELOAD_SLM"
ln -s "$IN" "$RUN/INPUT_FILES"; cp "$HERE/fortran/build/fortran_test/exe" "$RUN/"
( cd "$RUN" && OMP_NUM_THREADS=$THREADS OMP_STACKSIZE=32M SLEVOLVE_INNER=extrap \
  SLEVOLVE_SNAP=8 SLEVOLVE_SNAP_FIELDS=tgrid ./exe 8 > run.log 2>&1 ) || { tail -20 "$RUN/run.log"; echo "FAIL (Fortran run)"; exit 1; }

echo "== Python run"
SHTB=numpy; "$PYTHON" -c "import ducc0" 2>/dev/null && SHTB=ducc
OPT=""; "$PYTHON" -c "import numba" 2>/dev/null && OPT="--opt all"
"$PYTHON" "$HERE/python/run_ais.py" --steps 8 --lmax 511 --sht $SHTB $OPT --inner extrap --threads $THREADS \
  --icedir "$IN/icemodel" --icepfx synth_ice_ --topo "$IN/others/synth_topo.npy" \
  --earth "$IN/earthmodel/prem_512.l96C.ump3.lm1" --save-fields tgrid --out "$WORK/OUT_python" \
  > "$WORK/python.log" 2>&1 || { tail -20 "$WORK/python.log"; echo "FAIL (Python run)"; exit 1; }

"$PYTHON" - "$RUN/OUTPUT_SLM/tgrid8.npy" "$WORK/OUT_python/tgrid8.npy" "$IN/others/synth_topo.npy" <<'PY'
import sys
import numpy as np
f, p, t0 = (np.load(x) for x in sys.argv[1:4])
d = float(np.abs(f - p).max())
print(f"tgrid8: max |Fortran - Python| = {d:.3e} m (max |topography change| {np.abs(f - t0).max():.1f} m)")
ok = d < 1e-9
print("test_fortran_python: " + ("PASS" if ok else "FAIL"))
sys.exit(0 if ok else 1)
PY
