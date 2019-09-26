#!/bin/bash -e

# author: Ole Schuett

if (( $# < 1 )); then
    echo "usage: run_test.sh <test_name> [additional-docker-run-args]"
    echo "example: run_test.sh python"
    exit 1
fi

set -e
TESTNAME=$1
shift
echo "Running ${TESTNAME} ..."

CP2K_LOCAL=$(realpath ../../)
set -x

# SYS_PTRACE needed by LeakSanitizer.
docker run -i --init --rm --cap-add=SYS_PTRACE \
  --volume "${CP2K_LOCAL}:/mnt/cp2k/:ro" \
  "$@" "img_cp2k_test_${TESTNAME}"

#EOF
