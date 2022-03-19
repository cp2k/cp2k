#!/bin/bash

# author: Ole Schuett

if (($# != 2)); then
  echo "Usage: test_build.sh <ARCH> <VERSION>"
  exit 1
fi

ARCH=$1
VERSION=$2

# shellcheck disable=SC1091
source /opt/cp2k-toolchain/install/setup

# Compile cp2k.
echo -en "Compiling cp2k... "
cd /opt/cp2k || exit 1
if make -j ARCH="${ARCH}" VERSION="${VERSION}" &> make.out; then
  echo -e "done."
  echo -e "\nSummary: Compilation works fine."
  echo -e "Status: OK\n"
else
  echo -e "failed.\n\n"
  tail -n 100 make.out
  echo -e "\nSummary: Compilation failed."
  echo -e "Status: FAILED\n"
fi

exit 0 # Prevent CI from overwriting our summary message.

#EOF
