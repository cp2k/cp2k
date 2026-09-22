#!/bin/bash -e

# author: Thomas D. Kühne

# shellcheck disable=SC1091
source /opt/cp2k-toolchain/install/setup

echo -e "\n========== Installing Python Test Dependencies =========="
apt-get update -qq
apt-get install -qq --no-install-recommends python3 python3-venv python3-pip
rm -rf /var/lib/apt/lists/*

python3 -m venv /opt/python-interface-venv
export PATH="/opt/python-interface-venv/bin:$PATH"
cd /opt/cp2k
pip3 install './python[test]'
pip3 check

echo -e "\n========== Direct Python Interface Tests =========="
export CP2K_TEST_LIBRARY=/opt/cp2k/build/src/libcp2k.so
export CP2K_TEST_EXECUTABLE=/opt/cp2k/build/bin/cp2k.ssmp
export CP2K_DATA_DIR=/opt/cp2k/data
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
test -s "${CP2K_TEST_LIBRARY}"
test -x "${CP2K_TEST_EXECUTABLE}"
mkdir -p /workspace/artifacts

if timeout 15m python3 -m pytest python/tests -q -ra \
  --basetemp=/workspace/artifacts/python-tests; then
  echo -e "\nSummary: Direct Python interface tests passed."
  echo -e "Status: OK\n"
else
  echo -e "\nSummary: Direct Python interface tests failed."
  echo -e "Status: FAILED\n"
fi

#EOF
