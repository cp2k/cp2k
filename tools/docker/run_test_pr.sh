#!/bin/bash -l
# -*- coding: utf-8 -*-
# author: Tiziano Müller

set -o errexit
set -o nounset
set -o pipefail

if [[ $# -lt 3 ]]; then
    echo "usage: $0 <test_name> <PR ID or branch name> <Pull-Request ref/SHA> [optional-docker-run-arguments]"
    echo "example: $0 python 506 d3b9e46d0f2a0588c23fcb00f20e134f3fc302e8"
    exit 1
fi

TESTNAME=$1
[[ $2 =~ ^[0-9]+$ ]] && GIT_BRANCH="pull/$2/head" || GIT_BRANCH="$2"
GIT_REF=$3

shift 3

echo "Running ${TESTNAME} on Branch ${GIT_BRANCH} (ref: ${GIT_REF})..."

set -o xtrace

# SYS_PTRACE needed by LeakSanitizer.
${DOCKER:-docker} run -i --init --rm --cap-add=SYS_PTRACE \
    -e "GIT_BRANCH=${GIT_BRANCH}" -e "GIT_REF=${GIT_REF}" \
    "$@" "img_cp2k_test_${TESTNAME}"
