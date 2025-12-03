#!/bin/bash

# authors: Ole Schuett, Matthias Krack

if (($# != 2)); then
  echo "ERROR: Script build_cp2k_cmake.sh expects exactly two arguments"
  echo "Usage: build_cp2k_cmake.sh <PROFILE> <VERSION>"
  exit 1
fi

PROFILE=$1
VERSION=$2

if [ -n "${GIT_COMMIT_SHA}" ]; then
  echo "git:${GIT_COMMIT_SHA::7}" > REVISION
fi

echo "==================== Building CP2K ===================="
# shellcheck disable=SC1091
source ./cmake/cmake_cp2k.sh "${PROFILE}" "${VERSION}"

# Compile CP2K
echo -en '\nCompiling CP2K ... '
if ninja --verbose &> ninja.log; then
  echo "done"
else
  echo -e "failed.\n\n"
  tail -n 1000 ninja.log
  mkdir -p /workspace/artifacts/
  cp ninja.log /workspace/artifacts/
  echo -e "\nSummary: Compilation failed"
  echo -e "Status: FAILED\n"
  exit 1
fi

# Fake installation of data files.
cd ..
mkdir -p ./share/cp2k
ln -s ../../data ./share/cp2k/data

#EOF
