#!/bin/bash -e

# author: Ole Schuett

# shellcheck disable=SC1091
source /opt/cp2k-toolchain/install/setup

echo -e "\n========== Compiling CP2K =========="
cd /workspace/cp2k
make -j VERSION="popt" cp2k
make -j VERSION="psmp" cp2k
make -j VERSION="ssmp" cp2k

echo -e "\n========== Running Scaling Test =========="
cd ./tests/QS/benchmark
../../../tools/regtesting/test_scaling.py 30.0 ../../../exe/local/ H2O-32.inp -550.50556087853511

#EOF
