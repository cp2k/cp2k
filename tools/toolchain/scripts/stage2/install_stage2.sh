#!/bin/bash -e

# TODO: Review and if possible fix shellcheck errors.
# shellcheck disable=all

./scripts/stage2/install_mathlibs.sh
./scripts/stage2/install_gmp.sh

#EOF
