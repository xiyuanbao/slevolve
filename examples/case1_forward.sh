#!/usr/bin/env bash
# Case 1: Antarctic forward run, 550 x 1 yr (commands of the Benchmarks page, docs/src/benchmarks.md).
# usage: INPUT=<INPUT_FILES dir> examples/case1_forward.sh [fortran|python]
#   INPUT holds icemodel/AISload_gd_Run85_0_0_<k>(.npy|.bin) (total ice, tools/build_case1_total_ice.py),
#   earthmodel/prem_512.l60K2C.sum18p6.dum19p2.tz19p4.lm22 (shipped in lovenum/; copied there if missing),
#   others/topo0_archive (from the Han et al. 2022 OSF archive); grid/GLlat_512.txt GLlon_512.txt
#   (OSF archive, or written by tools/make_grid.py when grid/ is missing; the solvers do not read them).
set -eo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
[[ -n "${INPUT:-}" ]] || { echo "set INPUT=<INPUT_FILES directory>" >&2; exit 2; }
INPUT="$(cd "$INPUT" && pwd)"; THREADS="${THREADS:-8}"
LOVE=prem_512.l60K2C.sum18p6.dum19p2.tz19p4.lm22
if [[ ! -f "$INPUT/earthmodel/$LOVE" ]]; then
  echo "copying lovenum/$LOVE to $INPUT/earthmodel/"; mkdir -p "$INPUT/earthmodel"; cp "$HERE/lovenum/$LOVE" "$INPUT/earthmodel/"
fi
[[ -f "$INPUT/grid/GLlat_512.txt" ]] || python3 "$HERE/tools/make_grid.py" --outdir "$INPUT/grid"
case "${1:-fortran}" in
fortran)
  CASE=case1 bash "$HERE/fortran/build.sh"
  mkdir -p RUN_case1/OUTPUT_SLM RUN_case1/ICELOAD_SLM
  ln -sfn "$INPUT" RUN_case1/INPUT_FILES
  cp "$HERE/fortran/build/case1/exe" RUN_case1/ && cd RUN_case1
  OMP_NUM_THREADS=$THREADS OMP_STACKSIZE=32M SLEVOLVE_INNER=extrap \
  SLEVOLVE_SNAP=550 SLEVOLVE_SNAP_FIELDS=tgrid ./exe 550 ;;
python)
  python3 "$HERE/python/run_ais.py" --steps 550 --threads $THREADS --sht ducc --lmax 511 --opt all \
    --inner extrap --icedir "$INPUT/icemodel" --save-fields tgrid,ocean \
    --topo "$INPUT/others/topo0_archive" \
    --earth "$INPUT/earthmodel/$LOVE" --out OUT_case1_py ;;
*) echo "usage: $0 [fortran|python]" >&2; exit 2 ;;
esac
