#!/bin/bash -e

# TODO: Review and if possible fix shellcheck errors.
# shellcheck disable=all

./scripts/stage4/install_libxsmm.sh
./scripts/stage4/install_scalapack.sh
./scripts/stage4/install_cosma.sh

#EOF
