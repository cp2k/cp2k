#!/bin/bash -e
# -*- coding: utf-8 -*-
# author: Ole Schuett, Tiziano Müller

set -o pipefail

echo -e "\n========== Running Formatting Test =========="
cd /workspace/cp2k

if [ -n "${GIT_REF}" ]; then
    # get the list of changed files with their status from git between this HEAD and the origin/master
    # filter out deleted files and get the last field (in case of renames the second one is the newly added)
    # run the `pre-commit` on them which will fail IF files get modified by it
    git diff --name-status origin/master... | awk '/^[MCRAU]/ {print $NF}' | xargs pre-commit run --files
    echo "Summary: checked $(git diff --name-status origin/master... | awk '/^[MCRAU]/ {print $NF}' | wc -l) files, no violations found"
    echo "Status: OK"  # if not, we don't get here and the ci_entrypoint reports the FAILED status
else
    ./tools/formatting/test_formatting.sh
fi
