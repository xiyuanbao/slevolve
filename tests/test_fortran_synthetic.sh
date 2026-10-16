#!/usr/bin/env bash
# Fortran-only test (no Python): build the Fortran solver and the Fortran synthetic-case generator,
# write the synthetic Earth (tools/make_synthetic_case.f90, 512 x 1024 grid), run 8 steps of 1 yr and
# compare the final topography tgrid8 with the shipped reference reference/synthetic/tgrid8.npy
# (computed with the Python solver from the identical inputs of tools/make_synthetic_case.py).
#
#   bash tests/test_fortran_synthetic.sh [WORKDIR]        (default WORKDIR: ./slevolve_fortran_synth)
#
# Needs the Fortran toolchain (gfortran, FFTW 3, OpenMP; OpenBLAS outside macOS); see the Installation
# page, or set CONDA_ENV=<env> / FC / FFTW_PREFIX / OPENBLAS_PREFIX / BLASLIB as for fortran/build.sh.
# THREADS sets the thread count (default 4).  Prints PASS when max |tgrid8 - reference| < 1e-9 m and
# exits 0; otherwise prints FAIL and exits 1.  Takes about half a minute.
set -eo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
WORK="${1:-$PWD/slevolve_fortran_synth}"; THREADS="${THREADS:-4}"
mkdir -p "$WORK"; WORK="$(cd "$WORK" && pwd)"
IN="$WORK/INPUT_FILES"
REF="$HERE/reference/synthetic/tgrid8.npy"
[[ -f "$REF" ]] || { echo "FAIL (reference $REF missing)"; exit 1; }

echo "== Fortran build (CASE unset, NSTEPS=8, defaults)"
NSTEPS=8 ICEMODEL=synth_ice_ TOPOMODEL=synth_topo PLANETMODEL=prem_512.l96C.ump3.lm1 TAG=fortran_synth \
  bash "$HERE/fortran/build.sh" > "$WORK/build.log" 2>&1 || { tail -20 "$WORK/build.log"; echo "FAIL (build)"; exit 1; }
tail -1 "$WORK/build.log"
B="$HERE/fortran/build/fortran_synth"
if [[ -n "${CONDA_ENV:-}" ]]; then   # same toolchain as build.sh
  # shellcheck disable=SC1091
  source "${CONDA_BASE:-$(conda info --base)}/bin/activate" "$CONDA_ENV"
fi
FC="${FC:-gfortran}"
( cd "$B" && $FC -O2 -ffree-line-length-none "$HERE/tools/make_synthetic_case.f90" npy_io_mod.o -o make_synthetic_case \
    && $FC -O2 -ffree-line-length-none "$HERE/tests/npy_maxdiff.f90" npy_io_mod.o -o npy_maxdiff ) >> "$WORK/build.log" 2>&1 \
  || { tail -20 "$WORK/build.log"; echo "FAIL (build of make_synthetic_case / npy_maxdiff)"; exit 1; }

echo "== synthetic case (Fortran generator) -> $IN"
( cd "$HERE" && "$B/make_synthetic_case" --outdir "$IN" --love lovenum/prem_512.l96C.ump3.lm1 )

echo "== Fortran run"
RUN="$WORK/RUN_fortran"; rm -rf "$RUN"; mkdir -p "$RUN/OUTPUT_SLM" "$RUN/ICELOAD_SLM"
ln -s "$IN" "$RUN/INPUT_FILES"; cp "$B/exe" "$RUN/"
( cd "$RUN" && OMP_NUM_THREADS=$THREADS OMP_STACKSIZE=32M SLEVOLVE_INNER=extrap \
  SLEVOLVE_SNAP=8 SLEVOLVE_SNAP_FIELDS=tgrid ./exe 8 > run.log 2>&1 ) || { tail -20 "$RUN/run.log"; echo "FAIL (Fortran run)"; exit 1; }

echo "== compare tgrid8 with reference/synthetic/tgrid8.npy"
if "$B/npy_maxdiff" "$RUN/OUTPUT_SLM/tgrid8.npy" "$REF" 1e-9; then
  echo "test_fortran_synthetic: PASS"
else
  echo "test_fortran_synthetic: FAIL"; exit 1
fi
