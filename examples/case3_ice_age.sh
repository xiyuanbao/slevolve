#!/usr/bin/env bash
# Case 3: ICE-NADI-1, 50 x 500 yr, ice age run with outer iteration on the initial topography.
# usage: INPUT=<INPUT_FILES dir from tools/prepare_case3_nadi.py> examples/case3_ice_age.sh [fortran|python] [NOUTER]
#   INPUT/earthmodel/ must contain prem_512.l96C.ump3.lm1 (copy it from lovenum/).
set -eo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
[[ -n "${INPUT:-}" ]] || { echo "set INPUT=<INPUT_FILES directory>" >&2; exit 2; }
INPUT="$(cd "$INPUT" && pwd)"; THREADS="${THREADS:-8}"; NOUTER="${2:-5}"
EARTH="$INPUT/earthmodel/prem_512.l96C.ump3.lm1"
# first-pass initial topography: segmented, partially relaxed guess
[[ -f T0_segmented.npy ]] || python3 "$HERE/python/init_topo_segmented.py" --inputs "$INPUT" --icepfx nadi_ \
  --nfiles 51 --dt 500 --earth "$EARTH" --lmax 511 --out T0_segmented.npy
case "${1:-fortran}" in
fortran)
  CASE=case3 bash "$HERE/fortran/build.sh"
  OMP_NUM_THREADS=$THREADS SLEVOLVE_INNER=extrap SLEVOLVE_REUSE=warm INPUT="$INPUT" \
    bash "$HERE/examples/run_outer_mem.sh" case3 \
    case3 "$NOUTER" T0_segmented.npy ;;
python)
  python3 "$HERE/python/run_outer.py" --iters "$NOUTER" --inner extrap --reuse warm --lmax 511 \
    --threads $THREADS --inputs "$INPUT" --earth "$EARTH" --topo0 T0_segmented.npy --out OUT_case3_py ;;
*) echo "usage: $0 [fortran|python] [NOUTER]" >&2; exit 2 ;;
esac
