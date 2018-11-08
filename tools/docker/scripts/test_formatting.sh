#!/bin/bash -e

# author: Ole Schuett

echo -e "\n========== Running Formatting Test =========="
cd /workspace/cp2k
./tools/formatting/test_formatting.sh

#EOF
