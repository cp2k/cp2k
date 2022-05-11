#!/bin/bash -e

# TODO: Review and if possible fix shellcheck errors.
# shellcheck disable=all

./scripts/stage6/install_quip.sh
./scripts/stage6/install_gsl.sh
./scripts/stage6/install_plumed.sh

#EOF
