#!/usr/bin/env bash
# fortran/build.sh [TPW] [CHECKMARINE] [NSTEPS]
#
# Builds the SLEvolve Fortran solver from fortran/src into fortran/build/<TAG>/exe, plus the tools
# make_cache, text2npy and npy2text.  Usual call (benchmark presets, all defaults):
#     CASE=case1 bash fortran/build.sh          ->  fortran/build/case1/exe
# Without CASE the build directory is fortran/build/custom.  A suffix is appended to the directory
# name only for a switch that differs from its default (e.g. SHT=blas -> case1_blas, OMP=0 -> case1_serial).
# The last line printed is "built: <path of exe> (...)".  The linker may print
# "ld: warning: duplicate -rpath ... ignored"; this is harmless.
#
# Positional (normally none of them is given):
#   TPW = T|F          rotational feedback (default T)
#   CHECKMARINE = T|F  marine check: ice below flotation is removed from the load (default T)
#   NSTEPS             number of time steps (default: the CASE value, 550 without CASE)
#
# Case presets (set the defaults below; any variable given explicitly wins):
#   CASE=case1  Antarctic forward run, 550 x 1 yr:   TOPOMODEL=topo0_archive
#               PLANETMODEL=prem_512.l60K2C.sum18p6.dum19p2.tz19p4.lm22  DRIVER=run_from_zero_plus
#   CASE=case2  two glacial cycles, 1200 x 200 yr:   DT=200 NSTEPS=1200 INITIAL_TOPO=F
#               PLANETMODEL=prem_512.l120C.ump5.lm5 ICEMODEL=iceload TOPOMODEL=etopo2_512_ice6gC
#               DRIVER=run_outer_plus  (./exe 1200 1 = one outer pass from etopo2_512_ice6gC)
#   CASE=case3  ICE-NADI-1, 50 x 500 yr, outer iteration: DT=500 NSTEPS=50 INITIAL_TOPO=F
#               PLANETMODEL=prem_512.l96C.ump3.lm1 ICEMODEL=nadi_ TOPOMODEL=topo_pres
#               DRIVER=run_outer_mem   (./exe 50 K0 K1)
#   (no CASE)   the defaults below (forward run, DT=1, case-1 file names); set the file names and
#               NSTEPS for your own case.
#
# user_specs switches (compiled in):
#   NSTEPS=<n>         number of steps; L_sim = Ldt1 = DT*NSTEPS
#   DT=<yr>            step length in integer years (default 1); must equal the spacing of the ice files
#   NORDER=<L>         spectral truncation (default 511 on the 512 x 1024 grid; 512 is selectable)
#   EPS1=<value>       inner tolerance epsilon1 (default 1.0E-5); epsilon2 stays 1.0E-5
#   INITIAL_TOPO=T|F   T: start from topo_initial (forward run); F: outer topography iteration
#   PLANETMODEL, ICEMODEL, TOPOMODEL, TOPO_INITIAL (default = TOPOMODEL)   input file names
#   LEGACY_CONSTANTS=1 rhoi 910, radius 6.371229e6 (older legacy code) instead of the
#                      default rhoi 920, radius 6.371e6 (Han et al. 2022)
#   RAMP=1             second-order load placement: each load increment accrues linearly over its
#                      step (default RAMP=0: the step placement of Kendall et al. 2005; suffix _ramp)
# Solver switches:
#   DRIVER=run_from_zero_plus|run_outer_plus|run_outer_mem   (default run_from_zero_plus)
#       run_from_zero_plus  ./exe LAST            forward run, dtime = DT
#       run_outer_plus      ./exe LAST [ITERSL]   one outer pass per process, dtime = DT
#       run_outer_mem       ./exe LAST K0 K1      outer passes K0..K1 in one process, dtime = DT
#   SHT=otf|blas|exact|ref   Legendre backend (default otf).  otf: on-the-fly recurrence (fastest
#       from about 4 threads up, half the memory of blas); blas: precomputed table + DGEMM.
#       exact, ref: verification only (reproduce the legacy transform), not for production runs.
#   FFT=fftw|fft991      (default fftw)      SYM=1|0  N-S symmetry of the transform (default 1)
#   OMP=N|0              (default 8) any N > 0 compiles with OpenMP (-fopenmp); the thread count is
#                        set at run time with OMP_NUM_THREADS, N itself is not used.
#                        OMP=0: serial build (directory suffix _serial)
#   BLASLIB=accelerate|openblas              (default accelerate on macOS, openblas otherwise)
#   TAG=<name>           build directory name (default derived as described above)
# Toolchain (environment):
#   FC=<compiler>        default gfortran
#   CONDA_ENV=<env>      optional: activate this conda env first (e.g. slevolve-fortran);
#                        CONDA_BASE defaults to `conda info --base`
#   FFTW_PREFIX=<dir>    FFTW 3 install prefix (lib/libfftw3); default $CONDA_PREFIX, else /usr/local
#   OPENBLAS_PREFIX=<dir> for BLASLIB=openblas; default $CONDA_PREFIX, else /usr/local
set -eo pipefail

PP="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
[[ -z "${C2021:-}" ]] || { echo "C2021 was removed: the Han et al. (2022) constants are the default; use LEGACY_CONSTANTS=1 for 910 / 6.371229e6" >&2; exit 2; }

DEFDT=1; DEFDRIVER=run_from_zero_plus
case "${CASE:-}" in
  "") DEFN=550 ;;
  case1) : "${TOPOMODEL:=topo0_archive}"; DEFN=550 ;;
  case2) : "${INITIAL_TOPO:=F}" "${PLANETMODEL:=prem_512.l120C.ump5.lm5}" "${ICEMODEL:=iceload}" \
           "${TOPOMODEL:=etopo2_512_ice6gC}"; DEFDT=200; DEFDRIVER=run_outer_plus; DEFN=1200 ;;
  case3) : "${INITIAL_TOPO:=F}" "${PLANETMODEL:=prem_512.l96C.ump3.lm1}" "${ICEMODEL:=nadi_}" \
           "${TOPOMODEL:=topo_pres}"; DEFDT=500; DEFDRIVER=run_outer_mem; DEFN=50 ;;
  *) echo "bad CASE=$CASE (case1, case2, case3 or unset)" >&2; exit 2 ;;
esac
if [[ "$(uname)" == "Darwin" ]]; then DEFBLAS=accelerate; else DEFBLAS=openblas; fi

TPWF="${1:-T}"; CMF="${2:-T}"; NSTEPS="${3:-${NSTEPS:-$DEFN}}"
SHT="${SHT:-otf}"; SYM="${SYM:-1}"; OMPN="${OMP:-8}"; BLASLIB="${BLASLIB:-$DEFBLAS}"
FFT="${FFT:-fftw}"; NORDER="${NORDER:-511}"; SRC=src; DT="${DT:-$DEFDT}"
DRIVER="${DRIVER:-$DEFDRIVER}"; INITIAL_TOPO="${INITIAL_TOPO:-T}"
PLANETMODEL="${PLANETMODEL:-prem_512.l60K2C.sum18p6.dum19p2.tz19p4.lm22}"
ICEMODEL="${ICEMODEL:-AISload_gd_Run85_0_0_}"; TOPOMODEL="${TOPOMODEL:-topo0_archive}"
TOPO_INITIAL="${TOPO_INITIAL:-$TOPOMODEL}"; FC="${FC:-gfortran}"
case "$DRIVER" in run_from_zero_plus|run_outer_plus|run_outer_mem) ;; *) echo "bad DRIVER=$DRIVER" >&2; exit 2 ;; esac
[[ "$NSTEPS" =~ ^[0-9]+$ && "$DT" =~ ^[0-9]+$ ]] || { echo "NSTEPS and DT must be integers" >&2; exit 2; }
LSIM=$((DT * NSTEPS))

case "$TPWF" in T) TPW_VAL=".true." ;; F) TPW_VAL=".false." ;; *) echo "bad TPW (T or F)" >&2; exit 2 ;; esac
case "$CMF"  in T) CM_VAL=".true."  ;; F) CM_VAL=".false."  ;; *) echo "bad CHECKMARINE (T or F)"  >&2; exit 2 ;; esac
case "$INITIAL_TOPO" in T) IT_VAL=".true." ;; F) IT_VAL=".false." ;; *) echo "bad INITIAL_TOPO" >&2; exit 2 ;; esac
case "${LEGACY_CONSTANTS:-0}" in 0|1) ;; *) echo "bad LEGACY_CONSTANTS (0 or 1)" >&2; exit 2 ;; esac
case "${RAMP:-0}" in 0) RAMP_VAL=".false." ;; 1) RAMP_VAL=".true." ;; *) echo "bad RAMP (0 or 1)" >&2; exit 2 ;; esac

# build directory: CASE name (or custom) plus one suffix per switch that differs from its default
SFX=""
[[ "$SHT" == otf ]] || SFX="${SFX}_$SHT"
[[ "$FFT" == fftw ]] || SFX="${SFX}_$FFT"
[[ "$SYM" == 1 ]] || SFX="${SFX}_sym$SYM"
[[ "$OMPN" != 0 ]] || SFX="${SFX}_serial"
[[ "$BLASLIB" == "$DEFBLAS" ]] || SFX="${SFX}_$BLASLIB"
[[ "$NORDER" == 511 ]] || SFX="${SFX}_n$NORDER"
[[ "$TPWF" == T ]] || SFX="${SFX}_tpwF"
[[ "$CMF" == T ]] || SFX="${SFX}_cmF"
[[ "$NSTEPS" == "$DEFN" ]] || SFX="${SFX}_N$NSTEPS"
[[ "$DT" == "$DEFDT" ]] || SFX="${SFX}_dt$DT"
[[ "$DRIVER" == "$DEFDRIVER" ]] || SFX="${SFX}_$DRIVER"
[[ -z "${EPS1:-}" ]] || SFX="${SFX}_eps$EPS1"
[[ "${LEGACY_CONSTANTS:-0}" == 0 ]] || SFX="${SFX}_legacyconst"
[[ "${RAMP:-0}" == 0 ]] || SFX="${SFX}_ramp"
TAG="${TAG:-${CASE:-custom}${SFX}}"
BUILD="$PP/build/$TAG"
rm -rf "$BUILD"; mkdir -p "$BUILD"

# ---- support module (user_specs_mod, planets_mod) --------------------------
SUPPORT_SRC="$PP/$SRC/support_mod.f90"
SUPPORT_GEN="$BUILD/support_mod_gen.f90"
sed \
  -e "s|gridfolder='INPUT_FILES/others/'|gridfolder='INPUT_FILES/grid/'|" \
  -e "s|epsilon1=1.0E30|epsilon1=${EPS1:-1.0E-5}|" \
  -e "s|planetmodel='prem_512.l60K2C.sum18p6.dum19p2.tz19p4.lm22'|planetmodel='${PLANETMODEL}'|" \
  -e "s|icemodel='AISload_gd_Run85_0_0_'|icemodel='${ICEMODEL}'|" \
  -e "s|topomodel='etopo2_512_AIS_Run85_0_0'|topomodel='${TOPOMODEL}'|" \
  -e "s|topo_initial='etopo2_512_AIS_Run85_0_0'|topo_initial='${TOPO_INITIAL}'|" \
  -e "s|norder=512, nglv=512|norder=${NORDER}, nglv=512|" \
  -e "s|iceVolume=\.false\.|iceVolume=.true.|" \
  -e "s|checkmarine=\.false\., tpw=\.true\.|checkmarine=${CM_VAL}, tpw=${TPW_VAL}|" \
  -e "s|input_times=\.false\., initial_topo=\.true\.|input_times=.false., initial_topo=${IT_VAL}|" \
  -e "s|integer :: L_sim=2|integer :: L_sim=${LSIM}|" \
  -e "s|dt1=1,dt2=1,dt3=1,dt4=1|dt1=${DT},dt2=${DT},dt3=${DT},dt4=${DT}|" \
  -e "s|Ldt1=2,Ldt2=0,Ldt3=0,Ldt4=0|Ldt1=${LSIM},Ldt2=0,Ldt3=0,Ldt4=0|" \
  -e "s|logical :: ramp=\.false\.|logical :: ramp=${RAMP_VAL}|" \
  $( [[ "${LEGACY_CONSTANTS:-0}" == "1" ]] && echo "-e s|rhoi=920.0|rhoi=910.0| -e s|radius=6.371E6|radius=6.371229E6|" ) \
  "$SUPPORT_SRC" > "$SUPPORT_GEN"

chk() { grep -q -- "$1" "$SUPPORT_GEN" || { echo "SED CHECK FAILED: $1" >&2; exit 3; }; }
for c in "gridfolder='INPUT_FILES/grid/'" "epsilon1=${EPS1:-1.0E-5}" "epsilon2=1.0E-5" \
         "planetmodel='${PLANETMODEL}'" "icemodel='${ICEMODEL}'" "topomodel='${TOPOMODEL}'" "topo_initial='${TOPO_INITIAL}'" \
         "iceVolume=.true." "initial_topo=${IT_VAL}" "calcRG=.true." "input_times=.false." "coupling=.false." \
         "norder=${NORDER}, nglv=512" "L_sim=${LSIM}" "dt1=${DT},dt2=${DT},dt3=${DT},dt4=${DT}" \
         "Ldt1=${LSIM},Ldt2=0,Ldt3=0,Ldt4=0" "tpw=${TPW_VAL}" "checkmarine=${CM_VAL}" "ramp=${RAMP_VAL}"; do chk "$c"; done
if [[ "${LEGACY_CONSTANTS:-0}" == "1" ]]; then chk "rhoi=910.0"; chk "radius=6.371229E6"; else chk "rhoi=920.0"; chk "radius=6.371E6"; fi

# ---- toolchain -------------------------------------------------------------
if [[ -n "${CONDA_ENV:-}" ]]; then
  # shellcheck disable=SC1091
  source "${CONDA_BASE:-$(conda info --base)}/bin/activate" "$CONDA_ENV"
fi
FFTW_PREFIX="${FFTW_PREFIX:-${CONDA_PREFIX:-/usr/local}}"
OPENBLAS_PREFIX="${OPENBLAS_PREFIX:-${CONDA_PREFIX:-/usr/local}}"

FF="-O3 -ffree-line-length-none -fdefault-real-8 -fdefault-double-8 -cpp"
SHTFF=""
case "$SHT" in
  blas)  FF="$FF -DUSE_BLAS_SHT" ;;
  otf)   FF="$FF -DUSE_BLAS_SHT -DUSE_OTF_SHT" ;;
  # The exact kernels reproduce the reference summation *order*; the
  # remaining 1-ulp differences come from FMA contraction, which the
  # reference's complex accumulate does not get.  Disabling contraction in
  # this one translation unit makes the transform identical.  Every
  # other object keeps the baseline's default -ffp-contract=fast.
  exact) FF="$FF -DUSE_EXACT_SHT"; SHTFF="-ffp-contract=off" ;;
  ref)   ;;
  *) echo "bad SHT=$SHT" >&2; exit 2 ;;
esac
case "$FFT" in
  fftw)   FF="$FF -DHAVE_FFTW" ;;
  fft991) ;;
  *) echo "bad FFT=$FFT" >&2; exit 2 ;;
esac
case "$BLASLIB" in
  accelerate) LDLIBS="-framework Accelerate" ;;
  openblas)   LDLIBS="-L$OPENBLAS_PREFIX/lib -lopenblas -Wl,-rpath,$OPENBLAS_PREFIX/lib" ;;
  *) echo "bad BLASLIB=$BLASLIB" >&2; exit 2 ;;
esac
if [[ "$FFT" == "fftw" ]]; then LDLIBS="$LDLIBS -L$FFTW_PREFIX/lib -lfftw3 -Wl,-rpath,$FFTW_PREFIX/lib"; fi
if [[ "$OMPN" != "0" ]]; then FF="$FF -fopenmp"; fi

cd "$BUILD"
cp "$PP/$SRC/"*.f90 "$PP/$SRC/"*.inc .
# symmetry switch, applied to the copy that is actually compiled
if [[ "$SYM" == "0" ]]; then
  sed 's/^ logical, public :: sht_sym  = .true./ logical, public :: sht_sym  = .false./' spharmt_otf.f90 > spharmt_otf.f90.tmp
  mv spharmt_otf.f90.tmp spharmt_otf.f90   # portable (GNU and BSD sed)
fi
grep -q "sht_sym  = .$( [[ $SYM == 1 ]] && echo true || echo false )." spharmt_otf.f90 || { echo "SYM sed failed" >&2; exit 3; }

$FC $FF -c "$(basename "$SUPPORT_GEN")"
$FC $FF -c npy_io_mod.f90
$FC $FF -c io_mod.f90
$FC $FF $SHTFF -c spharmt_otf.f90 -o spharmt.o
$FC $FF -c modal_state_mod.f90
$FC $FF -c sl_model_mod.f90 -o solver.o
$FC $FF -c "$DRIVER.f90" -o driver.o
$FC -O3 $( [[ "$OMPN" != "0" ]] && echo -fopenmp ) -o exe \
   "$(basename "${SUPPORT_GEN%.f90}").o" npy_io_mod.o io_mod.o spharmt.o \
   modal_state_mod.o solver.o driver.o $LDLIBS 2>&1 | grep -v 'duplicate -rpath' || true

[[ -x exe ]] || { echo "link failed" >&2; exit 4; }

# ---- tools -----------------------------------------------------------------
$FC -O3 -ffree-line-length-none "$PP/tools/make_cache.f90" -o make_cache
$FC -O2 -ffree-line-length-none "$PP/tools/text2npy.f90" npy_io_mod.o -o text2npy
$FC -O2 -ffree-line-length-none "$PP/tools/npy2text.f90" npy_io_mod.o -o npy2text

printf 'built: %s (case=%s driver=%s sht=%s fft=%s sym=%s omp=%s blas=%s norder=%s tpw=%s cm=%s ramp=%s N=%s dt=%s initial_topo=%s planet=%s ice=%s topo=%s)\n' \
   "$BUILD/exe" "${CASE:-none}" "$DRIVER" "$SHT" "$FFT" "$SYM" "$OMPN" "$BLASLIB" "$NORDER" "$TPW_VAL" "$CM_VAL" "$RAMP_VAL" \
   "$NSTEPS" "$DT" "$IT_VAL" "$PLANETMODEL" "$ICEMODEL" "$TOPOMODEL"
