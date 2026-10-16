#!/usr/bin/env bash
# Fortran vs Python with NON-UNIFORM time steps from a times file, on the synthetic Earth
# (tools/make_synthetic_case.py, 512 x 1024 grid, 8 steps of 1, 2, 5, 10, 20, 50, 100, 200 yr).
#
#   bash tests/test_fortran_python_vardt.sh [WORKDIR]     (default WORKDIR: ./slevolve_fortran_vardt)
#
# The Fortran solver is built with DT unset (the uniform default, 1 yr), so the steps come only from
# the times file: the Fortran run auto-detects INPUT_FILES/others/times, the Python run gets
# --times.  Checks
#   (1) max |tgrid8_Fortran - tgrid8_Python| < 1e-9 m with the non-uniform times file;
#   (2) regression: a uniform times file t_k = k yr (Fortran via SLEVOLVE_TIMES=<file>, Python via
#       --times) reproduces the run without a times file exactly (max |diff| = 0.0), in both solvers;
#   (3) the Fortran OUTPUT_SLM/times holds the times of the times file.
# Toolchain and interpreter as for tests/test_fortran_python.sh (CONDA_ENV, PYTHON, THREADS).
# Prints PASS and exits 0, otherwise FAIL and exits 1.  Takes about two minutes.
set -eo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
WORK="${1:-$PWD/slevolve_fortran_vardt}"; PYTHON="${PYTHON:-python3}"; THREADS="${THREADS:-4}"
mkdir -p "$WORK"; WORK="$(cd "$WORK" && pwd)"
IN="$WORK/INPUT_FILES"
echo "== synthetic case -> $IN"
"$PYTHON" "$HERE/tools/make_synthetic_case.py" --outdir "$IN"
rm -f "$IN/others/times"
TV="$WORK/times_vardt"; TU="$WORK/times_uniform"
printf '%s\n' 0 1 3 8 18 38 88 188 388 > "$TV"
printf '%s\n' 0 1 2 3 4 5 6 7 8 > "$TU"

echo "== Fortran build (CASE unset, NSTEPS=8, DT unset = 1)"
NSTEPS=8 ICEMODEL=synth_ice_ TOPOMODEL=synth_topo PLANETMODEL=prem_512.l96C.ump3.lm1 TAG=fortran_vardt \
  bash "$HERE/fortran/build.sh" > "$WORK/build.log" 2>&1 || { tail -20 "$WORK/build.log"; echo "FAIL (build)"; exit 1; }
tail -1 "$WORK/build.log"

# frun NAME TIMESFILE_FOR_OTHERS ENV_SLEVOLVE_TIMES
frun() {
  local RUN="$WORK/RUN_$1"; rm -rf "$RUN"; mkdir -p "$RUN/OUTPUT_SLM" "$RUN/ICELOAD_SLM" "$RUN/INPUT_FILES/others"
  for d in "$IN"/*; do [[ "$(basename "$d")" == others ]] || ln -s "$d" "$RUN/INPUT_FILES/"; done
  for f in "$IN"/others/*; do ln -s "$f" "$RUN/INPUT_FILES/others/"; done
  [[ -z "$2" ]] || cp "$2" "$RUN/INPUT_FILES/others/times"
  cp "$HERE/fortran/build/fortran_vardt/exe" "$RUN/"
  ( cd "$RUN" && if [[ -n "$3" ]]; then export SLEVOLVE_TIMES="$3"; else unset SLEVOLVE_TIMES; fi
    OMP_NUM_THREADS=$THREADS OMP_STACKSIZE=32M SLEVOLVE_INNER=extrap \
    SLEVOLVE_SNAP=8 SLEVOLVE_SNAP_FIELDS=tgrid ./exe 8 > run.log 2>&1 ) || { tail -20 "$RUN/run.log"; echo "FAIL (Fortran run $1)"; exit 1; }
}
echo "== Fortran runs: no times file / non-uniform (auto-detected) / uniform (SLEVOLVE_TIMES)"
frun none "" ""
frun vardt "$TV" ""
frun uniform "" "$TU"

SHTB=numpy; "$PYTHON" -c "import ducc0" 2>/dev/null && SHTB=ducc
OPT=""; "$PYTHON" -c "import numba" 2>/dev/null && OPT="--opt all"
prun() {
  "$PYTHON" "$HERE/python/run_ais.py" --steps 8 --lmax 511 --sht $SHTB $OPT --inner extrap --threads $THREADS \
    --icedir "$IN/icemodel" --icepfx synth_ice_ --topo "$IN/others/synth_topo.npy" \
    --earth "$IN/earthmodel/prem_512.l96C.ump3.lm1" --save-fields tgrid --out "$WORK/OUT_python_$1" "${@:2}" \
    > "$WORK/python_$1.log" 2>&1 || { tail -20 "$WORK/python_$1.log"; echo "FAIL (Python run $1)"; exit 1; }
}
echo "== Python runs: no times file / non-uniform / uniform"
prun none
prun vardt --times "$TV"
prun uniform --times "$TU"

"$PYTHON" - "$WORK" "$TV" "$IN/others/synth_topo.npy" <<'PY'
import sys
import numpy as np
w, tv, t0p = sys.argv[1:4]
L = lambda p: np.load(p)
F = {k: L(f"{w}/RUN_{k}/OUTPUT_SLM/tgrid8.npy") for k in ("none", "vardt", "uniform")}
P = {k: L(f"{w}/OUT_python_{k}/tgrid8.npy") for k in ("none", "vardt", "uniform")}
t0 = L(t0p)
d1 = float(np.abs(F["vardt"] - P["vardt"]).max())
print(f"(1) non-uniform steps: tgrid8 max |Fortran - Python| = {d1:.3e} m "
      f"(max |topography change| {np.abs(F['vardt'] - t0).max():.1f} m; uniform 1-yr run "
      f"{np.abs(F['none'] - t0).max():.1f} m; vardt - uniform run differ by {np.abs(F['vardt'] - F['none']).max():.1f} m)")
rf = float(np.abs(F["uniform"] - F["none"]).max()); rp = float(np.abs(P["uniform"] - P["none"]).max())
d0 = float(np.abs(F["none"] - P["none"]).max())
print(f"(2) uniform times file vs no times file: Fortran max |diff| = {rf!r} m, Python = {rp!r} m "
      f"(no times file, Fortran - Python {d0:.3e} m)")
ft = np.loadtxt(f"{w}/RUN_vardt/OUTPUT_SLM/times"); tt = np.loadtxt(tv)
fn = np.loadtxt(f"{w}/RUN_none/OUTPUT_SLM/times"); fu = np.loadtxt(f"{w}/RUN_uniform/OUTPUT_SLM/times")
print(f"(3) Fortran OUTPUT_SLM/times (non-uniform): {ft.tolist()}")
ok = d1 < 1e-9 and rf == 0.0 and rp == 0.0 and np.array_equal(ft, tt) and np.array_equal(fn, fu)
print("test_fortran_python_vardt: " + ("PASS" if ok else "FAIL"))
sys.exit(0 if ok else 1)
PY
