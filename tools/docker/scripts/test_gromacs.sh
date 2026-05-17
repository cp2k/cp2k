#!/bin/bash -e

# author: Ole Schuett

# shellcheck disable=SC1091
source /opt/cp2k-toolchain/install/setup

cd build
ninja install &> install.log

echo -e "\n========== Testing CP2K CMake package config =========="
mkdir -p /tmp/cp2k_config_smoke
cat > /tmp/cp2k_config_smoke/CMakeLists.txt << 'EOF'
cmake_minimum_required(VERSION 3.22)
project(cp2k_config_smoke LANGUAGES C Fortran)

find_package(cp2k CONFIG REQUIRED)

if(NOT TARGET cp2k::cp2k)
  message(FATAL_ERROR "cp2k::cp2k target missing after find_package(cp2k)")
endif()

get_target_property(CP2K_IMPORTED_LOCATION cp2k::cp2k IMPORTED_LOCATION_NOCONFIG)
if(NOT CP2K_IMPORTED_LOCATION)
  get_target_property(CP2K_IMPORTED_LOCATION cp2k::cp2k IMPORTED_LOCATION_RELEASE)
endif()
if(NOT EXISTS "${CP2K_IMPORTED_LOCATION}")
  message(FATAL_ERROR "cp2k::cp2k imported location does not exist: ${CP2K_IMPORTED_LOCATION}")
endif()
EOF
if cmake -S /tmp/cp2k_config_smoke -B /tmp/cp2k_config_smoke/build \
  -DCMAKE_PREFIX_PATH="/opt/cp2k" &> cp2k_config_cmake.out; then
  echo "CP2K CMake package config works."
else
  echo -e "failed.\n\n"
  tail -n 100 cp2k_config_cmake.out
  mkdir -p /workspace/artifacts/
  cp cp2k_config_cmake.out /workspace/artifacts/
  echo -e "\nSummary: Testing CP2K CMake package config failed."
  echo -e "Status: FAILED\n"
  exit 0
fi

echo -e "\n========== Installing Dependencies =========="
apt-get update -qq
apt-get install -qq --no-install-recommends git
rm -rf /var/lib/apt/lists/*

echo -e "\n========== Building Gromacs v2025.2 =========="
echo -n "Cloning Gromacs repository ... "
git clone --quiet --depth=1 --single-branch -b v2025.2 https://gitlab.com/gromacs/gromacs.git /opt/gromacs
echo "done"
cd /opt/gromacs/
GROMACS_REVISION=$(git rev-parse --short HEAD)
mkdir build
cd build

echo -n "Configuring Gromacs ... "
if cmake .. \
  -DGMX_BUILD_OWN_FFTW=ON \
  -DBUILD_SHARED_LIBS=OFF \
  -DGMX_INSTALL_NBLIB_API=OFF \
  -DGMXAPI=OFF \
  -DGMX_CP2K=ON \
  -DCP2K_DIR="/opt/cp2k/lib/" \
  &> gromacs_cmake.out; then
  echo "done"
else
  echo -e "failed.\n\n"
  tail -n 100 gromacs_cmake.out
  mkdir -p /workspace/artifacts/
  cp gromacs_cmake.out /workspace/artifacts/
  echo -e "\nSummary: Configuring Gromacs failed."
  echo -e "Status: FAILED\n"
  exit 0
fi

echo -n "Compiling Gromacs... "
if make -j 32 all qmmm_applied_forces-test &> gromacs_make.out; then
  echo -e "done.\n\n"
  ./bin/gmx --version
else
  echo -e "failed.\n\n"
  tail -n 100 gromacs_make.out
  mkdir -p /workspace/artifacts/
  cp gromacs_make.out /workspace/artifacts/
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
