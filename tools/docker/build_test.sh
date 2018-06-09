#!/bin/bash

# author: Ole Schuett

if (( $# < 1 )); then
    echo "usage: build_test.sh <test_name> [additional-args]"
    echo "example: build_test.sh python"
    exit 1
fi

set -e
TESTNAME=$1
shift
echo "Building ${TESTNAME} ..."

BUILDARGS=""
if [ -f ./test_${TESTNAME}/buildargs.sh ]; then
  source ./test_${TESTNAME}/buildargs.sh
fi

set -x
docker build -t img_cp2k_test_${TESTNAME} ${BUILDARGS} "$@" ./test_${TESTNAME}/

#EOF
