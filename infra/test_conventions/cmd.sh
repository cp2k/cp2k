#!/bin/bash

# author: Ole Schuett

set -e

echo -e "\n========== Copying Changed Files =========="
rsync --exclude="*~"          \
      --exclude=".*/"         \
      --exclude="*.pyc"       \
      --exclude=/cp2k/obj/    \
      --exclude=/cp2k/lib/    \
      --exclude=/cp2k/exe/    \
      --ignore-times          \
      --update                \
      --verbose               \
      --recursive             \
      --checksum              \
      /opt/cp2k-local/  /opt/cp2k-master/

echo -e "\n========== Running Conventions Test =========="
source /opt/cp2k-toolchain/install/setup
cd /opt/cp2k-master/cp2k/tools/conventions/
#TODO port to Python3
./test_conventions.sh

#EOF
