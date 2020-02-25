#!/bin/bash -e
# -*- coding: utf-8 -*-
# author: Ole Schuett, Tiziano Müller

set -o pipefail

echo -e "\n========== Running Formatting Test =========="
cd /workspace/cp2k

# ./tools/formatting/test_formatting.sh

if pre-commit run --all-files ; then
    echo "Summary: All good."
    echo "Status: OK"
else
    echo "Summary: Found some issues."
    echo "Status: FAILED"
fi

#EOF
