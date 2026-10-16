#!/usr/bin/env bash
# examples/run_outer_mem.sh RUNNAME BUILDTAG NOUTER TOPO0 [LAST=50]   (requires INPUT=<INPUT_FILES dir>)
# In-process outer iteration (driver fortran/src/run_outer_mem.f90, build with DRIVER=run_outer_mem or CASE=case3):
# ONE process runs passes itersl = K0..K0+NOUTER-1.  TOPO0 seeding:
#   TOPO0=pres: K0=1 (the solver starts from the present topography of the build, TOPOMODEL);
#   else K0=2 with tgrid0_1 = TOPO0 (.npy, '<f8', C order) and pred_pres_topo_1 = present topography.
# Present topography used for the seed (only when TOPO0 is not pres):
#   TOPO_PRES=<path>  default INPUT/others/topo_pres; a .npy file (or <path>.npy if it exists) is
#                     read directly, anything else is read as legacy grid text with the build's
#                     text2npy (fortran/build/BUILDTAG/text2npy).  NLAT/NLON (default 512/1024)
#                     give the grid of a text file.
# No Python is needed.  The run directory is ./runs/RUNNAME (relative to the current directory).
# SLEVOLVE_REUSE, SLEVOLVE_INNER, SLEVOLVE_ICECACHE pass through.  Log run_mem.log.
set -eo pipefail
[[ $# -ge 4 ]] || { echo "usage: [TOPO_PRES=path] run_outer_mem.sh RUNNAME BUILDTAG NOUTER TOPO0 [LAST]" >&2; exit 2; }
PP="$(cd "$(dirname "${BASH_SOURCE[0]}")/../fortran" && pwd)"
NAME="$1"; TAG="$2"; NO="$3"; TOPO0="$4"; LAST="${5:-50}"
EXE="$PP/build/$TAG/exe"; T2N="$PP/build/$TAG/text2npy"; RUN="$PWD/runs/$NAME"
[[ -n "${INPUT:-}" ]] || { echo "set INPUT=<INPUT_FILES directory>" >&2; exit 2; }
INPUT="$(cd "$INPUT" && pwd)"
abspath() { echo "$(cd "$(dirname "$1")" && pwd)/$(basename "$1")"; }
[[ "$TOPO0" == "pres" ]] || TOPO0="$(abspath "$TOPO0")"
[[ -x "$EXE" ]] || { echo "missing $EXE" >&2; exit 3; }
if [[ "$TOPO0" != "pres" ]]; then
  TOPO_PRES="${TOPO_PRES:-$INPUT/others/topo_pres}"
  [[ "$TOPO_PRES" != *.npy && -f "$TOPO_PRES.npy" ]] && TOPO_PRES="$TOPO_PRES.npy"
  [[ -f "$TOPO_PRES" ]] || { echo "missing present topography $TOPO_PRES (set TOPO_PRES)" >&2; exit 2; }
  TOPO_PRES="$(abspath "$TOPO_PRES")"
fi

# npy2bin IN.npy OUT.bin : write the raw payload of a '<f8', C-order .npy file
npy2bin() {
  local ver hl off hdr
  [[ "$(head -c 6 "$1" | tail -c 5)" == "NUMPY" ]] || { echo "$1: not a .npy file" >&2; return 1; }
  ver=$(od -An -tu1 -j6 -N1 "$1" | tr -d ' ')
  if [[ "$ver" == "1" ]]; then
    hl=$(od -An -tu1 -v -j8 -N2 "$1" | tr -s " \n" " " | awk '{print $1 + 256*$2}'); off=$((10 + hl))
  else
    hl=$(od -An -tu1 -v -j8 -N4 "$1" | tr -s " \n" " " | awk '{print $1 + 256*$2 + 65536*$3 + 16777216*$4}'); off=$((12 + hl))
  fi
  hdr="$(head -c "$off" "$1" | tail -c "$hl")"
  [[ "$hdr" == *"'<f8'"* && "$hdr" == *"False"* ]] || { echo "$1: need dtype '<f8', fortran_order False (header: $hdr)" >&2; return 1; }
  tail -c +$((off + 1)) "$1" > "$2"
}

export SLEVOLVE_OUTER=picard SLEVOLVE_SNAP="${SLEVOLVE_SNAP-none}"
export OMP_STACKSIZE="${OMP_STACKSIZE:-32M}"
rm -rf "$RUN"; mkdir -p "$RUN/OUTPUT_SLM" "$RUN/ICELOAD_SLM"
ln -s "$INPUT" "$RUN/INPUT_FILES"; cp "$EXE" "$RUN/exe"; cd "$RUN"
if [[ "$TOPO0" == "pres" ]]; then K0=1; else
  K0=2
  npy2bin "$TOPO0" OUTPUT_SLM/tgrid0_1.bin
  if [[ "$TOPO_PRES" == *.npy ]]; then
    npy2bin "$TOPO_PRES" OUTPUT_SLM/pred_pres_topo_1.bin
  else
    [[ -x "$T2N" ]] || { echo "missing $T2N (built by fortran/build.sh)" >&2; exit 3; }
    "$T2N" --nlat "${NLAT:-512}" --nlon "${NLON:-1024}" --outdir "$RUN" "$TOPO_PRES" > /dev/null
    npy2bin "$RUN/$(basename "$TOPO_PRES").npy" OUTPUT_SLM/pred_pres_topo_1.bin
    rm -f "$RUN/$(basename "$TOPO_PRES").npy"
  fi
  [[ $(wc -c < OUTPUT_SLM/tgrid0_1.bin) -eq $(wc -c < OUTPUT_SLM/pred_pres_topo_1.bin) ]] \
    || { echo "TOPO0 and present topography have different sizes" >&2; exit 2; }
fi
K1=$((K0 + NO - 1))
{ echo "== $NAME : $TAG : ./exe $LAST $K0 $K1 : SLEVOLVE_INNER=${SLEVOLVE_INNER:-} SLEVOLVE_REUSE=${SLEVOLVE_REUSE:-} ICECACHE=${SLEVOLVE_ICECACHE:-1} OMP_NUM_THREADS=${OMP_NUM_THREADS:-} TOPO0=$TOPO0 TOPO_PRES=${TOPO_PRES:-}"
  echo "== uptime before: $(uptime)"; } > run_mem.log
# wall clock without Python: EPOCHREALTIME (bash >= 5), else whole seconds from date
now() { if [[ -n "${EPOCHREALTIME:-}" ]]; then echo "$EPOCHREALTIME"; else date +%s; fi; }
t0=$(now)
# resource report: BSD time -l on macOS, GNU time -v elsewhere if available, else the shell keyword
if [[ "$(uname)" == "Darwin" ]]; then TIMECMD=(/usr/bin/time -l)
elif /usr/bin/time -v true >/dev/null 2>&1; then TIMECMD=(/usr/bin/time -v)
else TIMECMD=(); fi
if [[ ${#TIMECMD[@]} -gt 0 ]]; then
  "${TIMECMD[@]}" ./exe "$LAST" "$K0" "$K1" >> run_mem.log 2> time_mem.log || { echo "RUN FAILED"; tail -20 run_mem.log; exit 5; }
else
  { time ./exe "$LAST" "$K0" "$K1" >> run_mem.log; } 2> time_mem.log || { echo "RUN FAILED"; tail -20 run_mem.log; exit 5; }
fi
t1=$(now)
echo "== uptime after: $(uptime)" >> run_mem.log
echo "== process wall: $(awk -v a="$t0" -v b="$t1" 'BEGIN{gsub(",",".",a); gsub(",",".",b); printf "%.2f", b-a}') s" >> run_mem.log
grep -E "real|maximum resident|Maximum resident|Elapsed" time_mem.log | head -2 >> run_mem.log
echo "$NAME $(grep -E '^== pass|process wall' run_mem.log | tr '\n' ' ')"
