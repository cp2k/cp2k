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
pip3 install -r /opt/cp2k/python/requirements-aiida.txt
cd /opt/cp2k
pip3 install './python[test,openmm,typing]'
pip3 check
python3 -c 'import aiida, aiida_cp2k, aiida_common_workflows'
export CP2K_TEST_AIIDA=1

# Keep the serial build, using the pinned upstream tarball instead of GitHub.
wget --quiet --tries=3 --timeout=30 -O /opt/lammps.tar.gz \
  https://download.lammps.org/tars/lammps-22Jul2025_update4.tar.gz
echo "b456a4d6f19d398dee9880d761594b4bfcb70f8dbb7c2b3aeee51815aa6419c6  /opt/lammps.tar.gz" | sha256sum --check
mkdir -p /opt/lammps
tar -xzf /opt/lammps.tar.gz -C /opt/lammps --strip-components=1
rm /opt/lammps.tar.gz
cmake -S /opt/lammps/cmake -B /opt/lammps/build \
  -DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS=ON -DBUILD_MPI=OFF \
  -DPKG_MISC=ON -DBUILD_OMP=OFF
cmake --build /opt/lammps/build --target lammps -j "$(nproc)"
export PYTHONPATH="/opt/lammps/python${PYTHONPATH:+:${PYTHONPATH}}"
export LD_LIBRARY_PATH="/opt/lammps/build${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"
export CP2K_TEST_LAMMPS=1

echo -e "\n========== Direct Python Interface Tests =========="
export CP2K_TEST_LIBRARY=/opt/cp2k/build/src/libcp2k.so
export CP2K_TEST_EXECUTABLE=/opt/cp2k/build/bin/cp2k.ssmp
export CP2K_DATA_DIR=/opt/cp2k/data
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
test -s "${CP2K_TEST_LIBRARY}"
test -x "${CP2K_TEST_EXECUTABLE}"
mkdir -p /workspace/artifacts

if python3 -m mypy --strict --config-file python/pyproject.toml &&
  timeout 15m python3 -m pytest python/tests -q -ra \
    --basetemp=/workspace/artifacts/python-tests; then
  echo -e "\nSummary: Direct Python interface type checks and tests passed."
  echo -e "Status: OK\n"
else
  echo -e "\nSummary: Direct Python interface type checks or tests failed."
  echo -e "Status: FAILED\n"
fi

#EOF
