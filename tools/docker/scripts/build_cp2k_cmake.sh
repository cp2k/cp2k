#!/bin/bash

# authors: Ole Schuett, Matthias Krack

if (($# != 2)); then
  echo "ERROR: Script ${BASH_SOURCE##*/} expects exactly two arguments"
  echo "Usage: ${BASH_SOURCE##*/} <PROFILE> <VERSION>"
  exit 1
fi

PROFILE=$1
VERSION=$2

echo "==================== Building CP2K ===================="

# shellcheck disable=SC1091
source ./cmake/cmake_cp2k.sh "${PROFILE}" "${VERSION}"

# Compile CP2K
echo -en '\nCompiling CP2K ... '
if ninja --verbose &> ninja.log; then
  echo "done"
else
  echo -e "failed.\n\n"
  tail -n 100 ninja.log
  mkdir -p /workspace/artifacts/
  cp ninja.log /workspace/artifacts/
  echo -e "\nSummary: Compilation failed"
  echo -e "Status: FAILED\n"
  exit 1
fi

#EOF
