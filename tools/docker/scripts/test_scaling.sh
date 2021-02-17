#!/bin/bash -e

# author: Ole Schuett

# shellcheck disable=SC1091
source /opt/cp2k-toolchain/install/setup

echo -e '\n========== Compiling CP2K =========='
cd /workspace/cp2k
for VERSION in 'psmp' 'ssmp'; do
  echo -n "Compiling cp2k.${VERSION}... "
  if make -j VERSION="${VERSION}" &> make.out; then
    echo "done."
  else
    echo -e "failed.\n\n"
    tail -n 100 make.out
    echo -e "\nSummary: Compilation failed."
    echo -e "Status: FAILED\n"
    exit 0
  fi
done

echo -e '\n========== Running Scaling Test =========='
cd ./benchmarks/QS
../../tools/regtesting/test_scaling.py 15.0 ../../exe/local/ H2O-32.inp -550.50556087853511

#EOF
