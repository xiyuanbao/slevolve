#!/usr/bin/env bash
# tools/prepare_case2_inputs.sh -- input tree for Case 2 (two glacial cycles, 1200 steps of 200 yr).
# Only the binary tree used by SLEvolve is built.
#
# usage: tools/prepare_case2_inputs.sh ICEDIR TOPO EARTH GRIDDIR OUT
#   ICEDIR  directory with iceload0..iceload1200, extracted unchanged from the OSF archive of
#           Han et al. (2022), doi:10.17605/OSF.IO/8PTFM,
#           Data/SeaLevelModel_TwoGlacialCycle_Simulations/INPUT_NHIS2GC.zip
#   TOPO    etopo2_512_ice6gC (same archive, CoupledModel_OneGlacialCycle_Simulations/INPUT_Others/)
#   EARTH   Love-number file in the solver's text layout (the archived prem_coll_512.l120C.ump5.lm5 is corrupt;
#           see README.md, "Benchmark cases and data")
#   GRIDDIR directory with GLlat_512.txt and GLlon_512.txt (created with tools/make_grid.py if it does not exist;
#           the solvers do not read these files, they document the grid)
#   OUT     output INPUT_FILES directory: icemodel/iceload<k>.bin (LE fp64, parsed by make_cache =
#           the solver's own list-directed read) + symlinks to the original text,
#           others/etopo2_512_ice6gC(.bin), earthmodel/, grid/
# environment: MAKE_CACHE = path of the make_cache tool built by fortran/build.sh
#              (default: the first fortran/build/*/make_cache)
set -eo pipefail
[[ $# -eq 5 ]] || { sed -n 5,19p "$0"; exit 2; }
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
abspath() { (cd "$(dirname "$1")" && echo "$(pwd)/$(basename "$1")"); }
[[ -f "$4/GLlat_512.txt" ]] || python3 "$HERE/make_grid.py" --outdir "$4"
ORIG="$(cd "$1" && pwd)"; TOPO="$(abspath "$2")"; EM="$(abspath "$3")"; GRID="$(cd "$4" && pwd)"; B="$5"
MC="${MAKE_CACHE:-$(ls "$HERE"/../fortran/build/*/make_cache 2>/dev/null | head -1)}"
[[ -x "$MC" ]] || { echo "make_cache not found: build with fortran/build.sh or set MAKE_CACHE" >&2; exit 3; }
mkdir -p "$B/icemodel" "$B/others" "$B/earthmodel"
ln -sfn "$GRID" "$B/grid"
ln -sf "$EM" "$B/earthmodel/$(basename "$EM")"
ln -sf "$TOPO" "$B/others/etopo2_512_ice6gC"
for k in $(seq 0 1200); do ln -sf "$ORIG/iceload$k" "$B/icemodel/iceload$k"; done
# fp64 caches (parallel batches of the solver-identical read)
seq 0 1200 | xargs -P 8 -n 60 sh -c 'f=""; for k in "$@"; do f="$f '"$ORIG"'/iceload$k"; done; '"$MC"' 512 '"$B"'/icemodel $f' _
"$MC" 512 "$B/others" "$TOPO"
ls "$B"/icemodel/*.bin | wc -l
