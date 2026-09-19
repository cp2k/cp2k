#!/bin/bash -e

# author: Ole Schuett

cat > /usr/bin/cp2k_shell << EndOfMessage
#!/bin/bash -e
export OMP_NUM_THREADS=1
source /opt/cp2k-toolchain/install/setup
/opt/cp2k/build/bin/cp2k.ssmp --shell "\$@"
EndOfMessage
chmod +x /usr/bin/cp2k_shell

# The cp2k main binary is used by ase/test/cp2k/cp2k_dcd.py.
# https://gitlab.com/ase/ase/merge_requests/1109
cat > /usr/bin/cp2k << EndOfMessage
#!/bin/bash -e
export OMP_NUM_THREADS=1
source /opt/cp2k-toolchain/install/setup
/opt/cp2k/build/bin/cp2k.ssmp "\$@"
EndOfMessage
chmod +x /usr/bin/cp2k

mkdir -p ~/.config/ase
cat > ~/.config/ase/config.ini << EndOfMessage
[cp2k]
cp2k_shell = /usr/bin/cp2k_shell
cp2k_main = /usr/bin/cp2k
EndOfMessage

echo -e "\n========== Installing Dependencies =========="
apt-get update -qq
apt-get install -qq --no-install-recommends \
  git \
  python3 \
  python3-dev \
  python3-venv \
  python3-pip \
  python3-wheel \
  python3-setuptools \
  build-essential
rm -rf /var/lib/apt/lists/*

# Create and activate a virtual environment for Python packages.
python3 -m venv /opt/venv
export PATH="/opt/venv/bin:$PATH"

echo -e "\n========== Installing ASE =========="
git clone --quiet --depth=1 --single-branch -b master https://gitlab.com/ase/ase.git /opt/ase
cd /opt/ase/
pip3 install ".[test]"

echo -e "\n========== Running ASE Tests =========="
ASE_REVISION=$(git rev-parse --short HEAD)
echo -

# Make test temp files available as artifacts.
export PYTEST_DEBUG_TEMPROOT=/workspace/artifacts
mkdir -p ${PYTEST_DEBUG_TEMPROOT}

if ! ase test -j 0 -c cp2k calculator/cp2k; then
  echo -e "\nSummary: Something is wrong with ASE commit ${ASE_REVISION}."
  echo -e "Status: FAILED\n"
  exit 0
fi

echo -e "\n========== Direct Python and MD Adapter Tests =========="
# Keep pytest independent of the top-level CP2K CMake build. This Linux job
# exercises the shared-library interface as well as ASE's existing shell path.
cd /opt/cp2k
# shellcheck disable=SC1091
source /opt/cp2k-toolchain/install/setup
pip3 install './python[test,openmm]'

# A small serial LAMMPS build avoids relying on a wheel's bundled MPI ABI.
git clone --quiet --depth=1 --branch stable_22Jul2025_update4 \
  https://github.com/lammps/lammps.git /opt/lammps
cmake -S /opt/lammps/cmake -B /opt/lammps/build \
  -DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS=ON -DBUILD_MPI=OFF \
  -DPKG_MISC=ON -DBUILD_OMP=OFF
cmake --build /opt/lammps/build --target lammps -j "$(nproc)"
export PYTHONPATH="/opt/lammps/python${PYTHONPATH:+:${PYTHONPATH}}"
export LD_LIBRARY_PATH="/opt/lammps/build${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"
export CP2K_TEST_LIBRARY=/opt/cp2k/build/src/libcp2k.so
export CP2K_DATA_DIR=/opt/cp2k/data
export CP2K_TEST_LAMMPS=1
export CP2K_TEST_EXECUTABLE=/opt/cp2k/build/bin/cp2k.ssmp
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
if timeout 15m python3 -m pytest python/tests -q \
  --basetemp=/workspace/artifacts/python-tests; then
  echo -e "\nSummary: ASE ${ASE_REVISION}, direct Python and MD adapters work fine."
  echo -e "Status: OK\n"
else
  echo -e "\nSummary: Direct Python / MD adapter tests failed."
  echo -e "Status: FAILED\n"
fi

#EOF
