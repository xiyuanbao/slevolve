#!/usr/bin/env bash
# Case 2: two glacial cycles, 1200 x 200 yr, one outer pass from the present topography.
# usage: INPUT=<INPUT_FILES dir from tools/prepare_case2_inputs.sh> examples/case2_one_pass.sh [fortran|python]
#   INPUT/earthmodel/ must contain prem_512.l120C.ump5.lm5 (copied there from lovenum/ if missing).
set -eo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
[[ -n "${INPUT:-}" ]] || { echo "set INPUT=<INPUT_FILES directory>" >&2; exit 2; }
INPUT="$(cd "$INPUT" && pwd)"; THREADS="${THREADS:-8}"
LOVE=prem_512.l120C.ump5.lm5
if [[ ! -f "$INPUT/earthmodel/$LOVE" ]]; then
  echo "copying lovenum/$LOVE to $INPUT/earthmodel/"; mkdir -p "$INPUT/earthmodel"; cp "$HERE/lovenum/$LOVE" "$INPUT/earthmodel/"
fi
case "${1:-fortran}" in
fortran)
  CASE=case2 bash "$HERE/fortran/build.sh"
  mkdir -p RUN_case2/OUTPUT_SLM RUN_case2/ICELOAD_SLM
  ln -sfn "$INPUT" RUN_case2/INPUT_FILES
  cp "$HERE/fortran/build/case2/exe" RUN_case2/ && cd RUN_case2
  OMP_NUM_THREADS=$THREADS OMP_STACKSIZE=32M SLEVOLVE_INNER=extrap ./exe 1200 1 ;;
python)
  python3 "$HERE/python/run_outer.py" --iters 1 --inner extrap --opt all --lmax 511 --threads $THREADS \
    --inputs "$INPUT" --icepfx iceload --nfiles 1201 --dt 200 --topo-pres "$INPUT/others/etopo2_512_ice6gC" \
    --earth "$INPUT/earthmodel/$LOVE" --save-steps 1,600,1200 --out OUT_case2_py ;;
*) echo "usage: $0 [fortran|python]" >&2; exit 2 ;;
esac
