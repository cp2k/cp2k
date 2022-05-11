#!/bin/bash -e

# TODO: Review and if possible fix shellcheck errors.
# shellcheck disable=all

./scripts/stage0/install_gcc.sh
./scripts/stage0/install_intel.sh
./scripts/stage0/setup_buildtools.sh
./scripts/stage0/install_cmake.sh

#EOF
