#!/bin/bash -e

# author: Ole Schuett

# shellcheck disable=SC1091
source /opt/cp2k-toolchain/install/setup

echo -e "\n========== Running Conventions Test =========="
cd /workspace/cp2k/tools/conventions
./test_conventions.sh #TODO port to Python3

#EOF
