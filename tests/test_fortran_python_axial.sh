#!/usr/bin/env bash
# Fortran vs Python with the self-consistent axial rotation (AXIAL_ROTATION=deg2|full / --axial-rotation)
# and the analytic rates, on the synthetic Earth (tools/make_synthetic_case.py, 512 x 1024 grid, L = 511,
# 8 steps of 1 yr).
#
#   bash tests/test_fortran_python_axial.sh [WORKDIR]     (default WORKDIR: ./slevolve_fortran_axial)
#
# Builds the Fortran solver twice (AXIAL_ROTATION=deg2 RAMP=1 and AXIAL_ROTATION=full with the default
# placement) and checks, at step 8,
#   (1) deg2 with the ramp placement: tgrid8, dSLdt8, dRdt8 and TPWrate8, Fortran against Python;
#   (2) full with the default placement: the same;
#   (3) the option acts: m3-dot of each leg differs from the Python run with the option off (same
#       placement) by more than 1e-4 relative, and is smaller in magnitude (the rotational part of I33
#       opposes the load part).
# Tolerances as in tests/test_fortran_python_ramp.sh: tgrid 1e-9 m; rate grids 1e-9 of their maximum;
# TPWrate (m-dot, Lambda-dot) 1e-8 relative.  Toolchain and interpreter as for tests/test_fortran_python.sh
# (CONDA_ENV, PYTHON, THREADS).  Prints PASS and exits 0, otherwise FAIL and exits 1.  Takes a few minutes.
set -eo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
WORK="${1:-$PWD/slevolve_fortran_axial}"; PYTHON="${PYTHON:-python3}"; THREADS="${THREADS:-4}"
mkdir -p "$WORK"; WORK="$(cd "$WORK" && pwd)"
IN="$WORK/INPUT_FILES"
echo "== synthetic case -> $IN"
"$PYTHON" "$HERE/tools/make_synthetic_case.py" --outdir "$IN" > "$WORK/make_case.log"
rm -f "$IN/others/times"

build() {  # build TAG RAMP AXIAL_ROTATION
  NSTEPS=8 ICEMODEL=synth_ice_ TOPOMODEL=synth_topo PLANETMODEL=prem_512.l96C.ump3.lm1 TAG=$1 RAMP=$2 \
    AXIAL_ROTATION=$3 bash "$HERE/fortran/build.sh" > "$WORK/build_$1.log" 2>&1 \
    || { tail -20 "$WORK/build_$1.log"; echo "FAIL (build $1)"; exit 1; }
  tail -1 "$WORK/build_$1.log"
}
echo "== Fortran builds (NSTEPS=8; AXIAL_ROTATION=deg2 RAMP=1, AXIAL_ROTATION=full)"
build fortran_axial_deg2_ramp 1 deg2
build fortran_axial_full 0 full

# frun NAME TAG
frun() {
  local RUN="$WORK/RUN_$1"; rm -rf "$RUN"; mkdir -p "$RUN/OUTPUT_SLM" "$RUN/ICELOAD_SLM"
  ln -s "$IN" "$RUN/INPUT_FILES"; cp "$HERE/fortran/build/$2/exe" "$RUN/"
  ( cd "$RUN" && OMP_NUM_THREADS=$THREADS OMP_STACKSIZE=32M SLEVOLVE_INNER=extrap SLEVOLVE_RATES=1 \
    SLEVOLVE_SNAP=8 SLEVOLVE_SNAP_FIELDS=tgrid,dSLdt,dRdt ./exe 8 > run.log 2>&1 ) || { tail -20 "$RUN/run.log"; echo "FAIL (Fortran run $1)"; exit 1; }
}
echo "== Fortran runs"
frun deg2_ramp fortran_axial_deg2_ramp
frun full fortran_axial_full

SHTB=numpy; "$PYTHON" -c "import ducc0" 2>/dev/null && SHTB=ducc
OPT=""; "$PYTHON" -c "import numba" 2>/dev/null && OPT="--opt all"
prun() {
  "$PYTHON" "$HERE/python/run_ais.py" --steps 8 --lmax 511 --sht $SHTB $OPT --inner extrap --threads $THREADS \
    --icedir "$IN/icemodel" --icepfx synth_ice_ --topo "$IN/others/synth_topo.npy" \
    --earth "$IN/earthmodel/prem_512.l96C.ump3.lm1" --save-fields tgrid,dSLdt,dRdt --out "$WORK/OUT_python_$1" "${@:2}" \
    > "$WORK/python_$1.log" 2>&1 || { tail -20 "$WORK/python_$1.log"; echo "FAIL (Python run $1)"; exit 1; }
}
echo "== Python runs"
prun deg2_ramp --ramp --axial-rotation deg2
prun full --axial-rotation full
prun off_ramp --ramp
prun off

"$PYTHON" - "$WORK" <<'PY'
import sys
import numpy as np
w = sys.argv[1]


def tpw(path):
    v = np.array(open(path).read().split(), dtype=float)
    return v[9:12], v[12:]          # m-dot, Lambda-dot (re, im pairs)


ok = True
for name, off in (("deg2_ramp", "off_ramp"), ("full", "off")):
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
    mo, _ = tpw(f"{w}/OUT_python_{off}/TPWrate8.txt")
    ch = mp[2] / mo[2] - 1.0
    line += f"; m3-dot vs option off {ch:+.3e}"
    good &= abs(ch) > 1e-4 and abs(mp[2]) < abs(mo[2])
    print(line + (" ok" if good else " FAILED"))
    ok &= good
print("test_fortran_python_axial: " + ("PASS" if ok else "FAIL"))
sys.exit(0 if ok else 1)
PY
