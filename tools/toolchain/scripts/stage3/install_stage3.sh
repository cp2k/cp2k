#!/bin/bash -e

# TODO: Review and if possible fix shellcheck errors.
# shellcheck disable=all

./scripts/stage3/install_fftw.sh
./scripts/stage3/install_libint.sh
./scripts/stage3/install_libxc.sh

#EOF
