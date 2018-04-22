#!/bin/bash

# author: Ole Schuett

set -e

echo -e "\n========== Rsyncing =========="
rsync --exclude="*~"          \
      --exclude=".*/"         \
      --exclude="*.pyc"       \
      --exclude=/cp2k/obj/    \
      --exclude=/cp2k/lib/    \
      --exclude=/cp2k/exe/    \
      --ignore-times          \
      --update                \
      --verbose               \
      --recursive             \
      --checksum              \
      /opt/cp2k-local/  /opt/cp2k-master/

echo -e "\n========== Running formatting test =========="
cd /opt/cp2k-master/cp2k
./tools/formatting/test_formatting.sh

#EOF
