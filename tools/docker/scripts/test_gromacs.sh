#!/bin/bash -e

# author: Ole Schuett

# shellcheck disable=SC1091
source /opt/cp2k-toolchain/install/setup

cd /workspace/cp2k
echo -n "Compiling libcp2k... "
if make -j VERSION=sdbg libcp2k &> make.out; then
  echo "done."
else
  echo -e "failed.\n\n"
  tail -n 100 make.out
  echo -e "\nSummary: Compiling libcp2k failed."
  echo -e "Status: FAILED\n"
  exit 0
fi

echo -e "\n========== Building Gromacs =========="
echo -n "Cloning Gromacs repository... "
git clone --quiet --depth=1 --single-branch -b master https://gitlab.com/gromacs/gromacs.git /opt/gromacs
echo "done."
cd /opt/gromacs/
GROMACS_REVISION=$(git rev-parse --short HEAD)
mkdir build
cd build

echo -n "Configuring Gromacs... "
if cmake .. \
  -DGMX_BUILD_OWN_FFTW=ON \
  -DBUILD_SHARED_LIBS=OFF \
  -DGMX_INSTALL_NBLIB_API=OFF \
  -DGMXAPI=OFF \
  -DGMX_CP2K=ON \
  -DCP2K_DIR="/workspace/cp2k/lib/local/sdbg/" \
  &> cmake.out; then
  echo "done."
else
  echo -e "failed.\n\n"
  tail -n 100 cmake.out
  echo -e "\nSummary: Configuring Gromacs failed."
  echo -e "Status: FAILED\n"
  exit 0
fi

echo -n "Compiling Gromacs... "
if make -j 32 all qmmm_applied_forces-test &> make.out; then
  echo -e "done.\n\n"
  ./bin/gmx --version
else
  echo -e "failed.\n\n"
  tail -n 100 make.out
  echo -e "\nSummary: Compiling Gromacs failed."
  echo -e "Status: FAILED\n"
  exit 0
fi

echo -e "\n========== Running Gromacs QM/MM unit test =========="
if ./bin/qmmm_applied_forces-test; then
  echo -e "\nSummary: Gromacs commit ${GROMACS_REVISION} works fine."
  echo -e "Status: OK\n"
else
  echo -e "\nSummary: Something is wrong with Gromacs commit ${GROMACS_REVISION}."
  echo -e "Status: FAILED\n"
fi

#EOF
