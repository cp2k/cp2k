#!/bin/bash -e

# author: Ole Schuett

if (($# < 1)); then
  echo "usage: run_test_master.sh <test_name> [additional-docker-run-args]"
  echo "example: run_test_master.sh python"
  exit 1
fi

set -e
TESTNAME=$1
shift

GIT_BRANCH="master"
GIT_REF="origin/master"
echo "Running ${TESTNAME} on Branch ${GIT_BRANCH} (ref: ${GIT_REF})..."

set -x

# SYS_PTRACE needed by LeakSanitizer.
${DOCKER:-docker} run -i --init --rm --cap-add=SYS_PTRACE \
  -e "GIT_BRANCH=${GIT_BRANCH}" -e "GIT_REF=${GIT_REF}" \
  "$@" "img_cp2k_test_${TESTNAME}"

#EOF
