#!/bin/bash

# author: Ole Schuett

set -e

echo -e "\n========== Copying Changed Files =========="
rsync --exclude="*~"          \
      --exclude=".*/"         \
      --exclude="*.pyc"       \
      --exclude=/obj/         \
      --exclude=/lib/         \
      --exclude=/exe/         \
      --executability         \
      --ignore-times          \
      --update                \
      --verbose               \
      --recursive             \
      --checksum              \
      /opt/cp2k-local/  /opt/cp2k-master/

rsync --exclude="*~"          \
      --exclude=".*/"         \
      --exclude="*.pyc"       \
      --executability         \
      --ignore-times          \
      --update                \
      --verbose               \
      --recursive             \
      --checksum              \
      /opt/cp2k-local/tools/toolchain/  /opt/cp2k-toolchain/

echo -e "\n========== Updating Toolchain =========="
cd /opt/cp2k-toolchain/
./install_cp2k_toolchain.sh --install-all --with-make=no

echo -e "\n========== Running Conventions Test =========="
source /opt/cp2k-toolchain/install/setup
cd /opt/cp2k-master/tools/conventions/
#TODO port to Python3
./test_conventions.sh

#EOF
