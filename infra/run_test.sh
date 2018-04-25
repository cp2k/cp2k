#!/bin/bash

# author: Ole Schuett

if [ "$#" -ne "1" ]; then
    echo "usage: run_test.sh <test_name>"
    echo "example: run_test.sh python"
    exit 1
fi

set -e
TESTNAME=$1
echo "Running ${TESTNAME} ..."
CP2K_LOCAL=`realpath ../../`

echo -n "Date: "
date --utc --rfc-3339=seconds
if git rev-parse; then
  git log -1 --pretty="%nCommitSHA: %H%nCommitTime: %ci%nCommitAuthor: %an%nCommitSubject: %s%n"
fi

set -x
docker build -t img_cp2k_test_${TESTNAME} ./test_${TESTNAME}/
docker run -i --init --volume ${CP2K_LOCAL}:/opt/cp2k-local/:ro img_cp2k_test_${TESTNAME}

#EOF
