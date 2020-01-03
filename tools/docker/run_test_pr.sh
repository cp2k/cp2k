#!/bin/bash -l
# -*- coding: utf-8 -*-
# author: Tiziano Müller, Ole Schütt

set -o errexit
set -o nounset
set -o pipefail

if [[ $# -lt 2 ]]; then
    echo "usage: $0 <test_name> <pr_number> [optional-docker-run-arguments]"
    echo "example: $0 python 506"
    exit 1
fi

TESTNAME=$1
PR_NUMBER=$2
shift 2

GIT_BRANCH="pull/${PR_NUMBER}/merge"
GIT_REF=$(curl -s "https://api.github.com/repos/cp2k/cp2k/pulls/${PR_NUMBER}" | jq  -r '.merge_commit_sha')

echo "Running ${TESTNAME} on Branch ${GIT_BRANCH} (ref: ${GIT_REF})..."

set -o xtrace

# SYS_PTRACE needed by LeakSanitizer.
${DOCKER:-docker} run -i --init --rm --cap-add=SYS_PTRACE \
    -e "GIT_BRANCH=${GIT_BRANCH}" -e "GIT_REF=${GIT_REF}" \
    "$@" "img_cp2k_test_${TESTNAME}"

#EOF
