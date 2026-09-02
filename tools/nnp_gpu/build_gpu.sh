#!/usr/bin/env bash
# Clean CUDA build of CP2K with the NNP GPU backend, for GPU verification.
#
# Verification must always run against a fresh build: an incremental build can
# keep a stale object when a Fortran PARAMETER or a struct layout changes, so
# this script removes the build directory first and prints the resulting binary
# hash and cp2kflags (which must list "nnp_gpu") for the run log.
#
# Configuration via environment variables (all optional):
#   CP2K_SRC        source tree                 (default: two levels up)
#   BUILD_DIR       build directory             (default: $CP2K_SRC/build_nnp_gpu)
#   GPU_ARCH        CUDA architecture number    (default: 80; RTX 5090 = 120)
#   CUDA_HOME       CUDA toolkit root, read by CMake from the environment
#   JOBS            parallel build jobs         (default: nproc)
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CP2K_SRC="${CP2K_SRC:-$(cd "$script_dir/../.." && pwd)}"
BUILD_DIR="${BUILD_DIR:-$CP2K_SRC/build_nnp_gpu}"
GPU_ARCH="${GPU_ARCH:-80}"
JOBS="${JOBS:-$( (command -v nproc > /dev/null && nproc) || echo 8)}"

echo "== NNP GPU build =="
echo "source:   $CP2K_SRC"
echo "build:    $BUILD_DIR"
echo "gpu arch: sm_$GPU_ARCH"
echo "commit:   $(git -C "$CP2K_SRC" rev-parse HEAD)"

rm -rf "$BUILD_DIR"
# MPI is requested explicitly because CP2K_USE_MPI follows CP2K_USE_EVERYTHING,
# which is OFF by default: without this the configure yields a serial cp2k.ssmp
# and step 3 of README.md, which compares the CPU and GPU paths at 1, 2 and 4
# ranks, cannot be run. Anything in "$@" is appended after these, so a caller
# who wants the serial binary can still pass -DCP2K_USE_MPI=OFF.
cmake -S "$CP2K_SRC" -B "$BUILD_DIR" -GNinja \
  -DCMAKE_BUILD_TYPE=Release \
  -DCP2K_USE_ACCEL=CUDA \
  -DCP2K_USE_MPI=ON \
  -DCMAKE_CUDA_ARCHITECTURES="$GPU_ARCH" \
  "$@"
cmake --build "$BUILD_DIR" --target cp2k-bin -j "$JOBS"

bin="$(find "$BUILD_DIR" -name 'cp2k.psmp' -o -name 'cp2k.ssmp' | head -1)"
if [[ -z "$bin" ]]; then
  echo "ERROR: no cp2k binary produced" >&2
  exit 1
fi

echo "== build complete =="
echo "binary:   $bin"
echo "sha256:   $(shasum -a 256 "$bin" | awk '{print $1}')"
echo "flags:    $("$bin" --version 2> /dev/null | tr ' ' '\n' | grep -c '^nnp_gpu$' > /dev/null && echo 'nnp_gpu present' || echo 'nnp_gpu MISSING')"
"$bin" --version 2> /dev/null | grep -i cp2kflags || true
