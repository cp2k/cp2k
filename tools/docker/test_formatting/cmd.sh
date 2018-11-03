#!/bin/bash

# author: Ole Schuett

set -e

echo -e "\n========== Copying Changed Files =========="
rsync --exclude="*~"          \
      --exclude=".*/"         \
      --exclude="*.pyc"       \
      --exclude=/obj/         \
      --exclude=/lib/         \
      --exclude=/exe/         \
      --executability         \
      --ignore-times          \
      --update                \
      --verbose               \
      --recursive             \
      --checksum              \
      /opt/cp2k-local/  /opt/cp2k-master/

echo -e "\n========== Running Formatting Test =========="
cd /opt/cp2k-master
./tools/formatting/test_formatting.sh

#EOF
