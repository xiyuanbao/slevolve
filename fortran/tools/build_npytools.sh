#!/usr/bin/env bash
# Build the Fortran .npy converters (text2npy, npy2text) into fortran/tools/bin (optional; build.sh also builds them).
set -eo pipefail
PP="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
if [[ -n "${CONDA_ENV:-}" ]]; then source "${CONDA_BASE:-$(conda info --base)}/bin/activate" "$CONDA_ENV"; fi
B="$PP/tools/bin"; mkdir -p "$B"; cd "$B"
gfortran -O2 -ffree-line-length-none -c "$PP/src/npy_io_mod.f90"
gfortran -O2 -ffree-line-length-none "$PP/tools/text2npy.f90" npy_io_mod.o -o text2npy
gfortran -O2 -ffree-line-length-none "$PP/tools/npy2text.f90" npy_io_mod.o -o npy2text
echo "built $B/text2npy $B/npy2text"
