#!/usr/bin/env bash
# Fortran vs Python with the second-order load placement (RAMP=1 / --ramp) and the analytic rates
# (SLEVOLVE_RATES=1 / --save-fields ...,dSLdt,dRdt), on the synthetic Earth (tools/make_synthetic_case.py,
# 512 x 1024 grid, L = 511, 8 steps).
#
#   bash tests/test_fortran_python_ramp.sh [WORKDIR]     (default WORKDIR: ./slevolve_fortran_ramp)
#
# Builds the Fortran solver twice (RAMP=1 and the default) and checks, at step 8,
#   (1) ramp, uniform 1-yr steps:  tgrid8, dSLdt8, dRdt8 and TPWrate8, Fortran against Python;
#   (2) ramp, non-uniform steps from a times file (1, 2, 5, 10, 20, 50, 100, 200 yr): the same;
#   (3) default placement with the rates switched on: the same, and tgrid8 equal to the run
#       without rates (the rates are output only).
# Tolerances: tgrid 1e-9 m; rate grids 1e-9 of their maximum; TPWrate (m-dot, Lambda-dot) 1e-8
# relative.  Toolchain and interpreter as for tests/test_fortran_python.sh (CONDA_ENV, PYTHON,
# THREADS).  Prints PASS and exits 0, otherwise FAIL and exits 1.  Takes a few minutes.
set -eo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
WORK="${1:-$PWD/slevolve_fortran_ramp}"; PYTHON="${PYTHON:-python3}"; THREADS="${THREADS:-4}"
mkdir -p "$WORK"; WORK="$(cd "$WORK" && pwd)"
IN="$WORK/INPUT_FILES"
echo "== synthetic case -> $IN"
"$PYTHON" "$HERE/tools/make_synthetic_case.py" --outdir "$IN" > "$WORK/make_case.log"
rm -f "$IN/others/times"
TV="$WORK/times_vardt"
printf '%s\n' 0 1 3 8 18 38 88 188 388 > "$TV"

build() {  # build TAG RAMP
  NSTEPS=8 ICEMODEL=synth_ice_ TOPOMODEL=synth_topo PLANETMODEL=prem_512.l96C.ump3.lm1 TAG=$1 RAMP=$2 \
    bash "$HERE/fortran/build.sh" > "$WORK/build_$1.log" 2>&1 || { tail -20 "$WORK/build_$1.log"; echo "FAIL (build $1)"; exit 1; }
  tail -1 "$WORK/build_$1.log"
}
echo "== Fortran builds (NSTEPS=8; RAMP=1 and default)"
build fortran_ramp 1
build fortran_test 0

# frun NAME TAG TIMESFILE RATES
frun() {
  local RUN="$WORK/RUN_$1"; rm -rf "$RUN"; mkdir -p "$RUN/OUTPUT_SLM" "$RUN/ICELOAD_SLM" "$RUN/INPUT_FILES/others"
  for d in "$IN"/*; do [[ "$(basename "$d")" == others ]] || ln -s "$d" "$RUN/INPUT_FILES/"; done
  for f in "$IN"/others/*; do ln -s "$f" "$RUN/INPUT_FILES/others/"; done
  [[ -z "$3" ]] || cp "$3" "$RUN/INPUT_FILES/others/times"
  cp "$HERE/fortran/build/$2/exe" "$RUN/"
  ( cd "$RUN" && OMP_NUM_THREADS=$THREADS OMP_STACKSIZE=32M SLEVOLVE_INNER=extrap SLEVOLVE_RATES=$4 \
    SLEVOLVE_SNAP=8 SLEVOLVE_SNAP_FIELDS=tgrid,dSLdt,dRdt ./exe 8 > run.log 2>&1 ) || { tail -20 "$RUN/run.log"; echo "FAIL (Fortran run $1)"; exit 1; }
}
echo "== Fortran runs"
frun ramp fortran_ramp "" 1
frun ramp_vardt fortran_ramp "$TV" 1
frun step fortran_test "" 1
frun step_norates fortran_test "" 0

SHTB=numpy; "$PYTHON" -c "import ducc0" 2>/dev/null && SHTB=ducc
OPT=""; "$PYTHON" -c "import numba" 2>/dev/null && OPT="--opt all"
prun() {
  "$PYTHON" "$HERE/python/run_ais.py" --steps 8 --lmax 511 --sht $SHTB $OPT --inner extrap --threads $THREADS \
    --icedir "$IN/icemodel" --icepfx synth_ice_ --topo "$IN/others/synth_topo.npy" \
    --earth "$IN/earthmodel/prem_512.l96C.ump3.lm1" --save-fields tgrid,dSLdt,dRdt --out "$WORK/OUT_python_$1" "${@:2}" \
    > "$WORK/python_$1.log" 2>&1 || { tail -20 "$WORK/python_$1.log"; echo "FAIL (Python run $1)"; exit 1; }
}
echo "== Python runs"
prun ramp --ramp
prun ramp_vardt --ramp --times "$TV"
prun step

"$PYTHON" - "$WORK" "$IN/others/synth_topo.npy" <<'PY'
import sys
import numpy as np
w, t0p = sys.argv[1:3]
t0 = np.load(t0p)


def tpw(path):
    v = np.array(open(path).read().split(), dtype=float)
    return v[9:12], v[12:]          # m-dot, Lambda-dot (re, im pairs)


ok = True
for name in ("ramp", "ramp_vardt", "step"):
    F, P = f"{w}/RUN_{name}/OUTPUT_SLM", f"{w}/OUT_python_{name}"
    tf, tp = np.load(f"{F}/tgrid8.npy"), np.load(f"{P}/tgrid8.npy")
    d = float(np.abs(tf - tp).max())
    line = f"({name}) tgrid8 {d:.2e} m"
    good = d < 1e-9
    for fld in ("dSLdt", "dRdt"):
        a, b = np.load(f"{F}/{fld}8.npy"), np.load(f"{P}/{fld}8.npy")
        r = float(np.abs(a - b).max() / np.abs(b).max())
        line += f"; {fld}8 {r:.2e} of max {np.abs(b).max():.3e} m/yr"
        good &= r < 1e-9
    mf, lf = tpw(f"{F}/TPWrate8")
    mp, lp = tpw(f"{P}/TPWrate8.txt")
    rm = float(np.abs(mf - mp).max() / np.abs(mp).max())
    rl = float(np.abs(lf - lp).max() / np.abs(lp).max())
    line += f"; m-dot {rm:.1e}, Lambda-dot {rl:.1e} (relative)"
    good &= rm < 1e-8 and rl < 1e-8
    print(line + (" ok" if good else " FAILED"))
    ok &= good
step_vs = float(np.abs(np.load(f"{w}/RUN_step/OUTPUT_SLM/tgrid8.npy") - np.load(f"{w}/RUN_step_norates/OUTPUT_SLM/tgrid8.npy")).max())
rs = float(np.abs(np.load(f"{w}/RUN_ramp/OUTPUT_SLM/tgrid8.npy") - np.load(f"{w}/RUN_step/OUTPUT_SLM/tgrid8.npy")).max())
print(f"Fortran default placement: tgrid8 with rates - without rates = {step_vs!r} m; "
      f"ramp - default placement {rs:.2f} m (max |topography change| {np.abs(np.load(f'{w}/RUN_step/OUTPUT_SLM/tgrid8.npy') - t0).max():.1f} m)")
ok &= step_vs == 0.0 and rs > 0.0
print("test_fortran_python_ramp: " + ("PASS" if ok else "FAIL"))
sys.exit(0 if ok else 1)
PY
