#!/bin/bash -e

# author: Ole Schuett

if (($# < 1)); then
  echo "usage: build_test.sh <test_name> [additional-args]"
  echo "example: build_test.sh python"
  exit 1
fi

set -e
TESTNAME=$1
shift
echo "Building ${TESTNAME} ..."

set -x
${DOCKER:-docker} build -t "img_cp2k_test_${TESTNAME}" "$@" -f "Dockerfile.test_${TESTNAME}" .

#EOF
