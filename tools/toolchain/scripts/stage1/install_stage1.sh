#!/bin/bash -e

# TODO: Review and if possible fix shellcheck errors.
# shellcheck disable=all

./scripts/stage1/install_mpich.sh
./scripts/stage1/install_openmpi.sh
./scripts/stage1/install_intelmpi.sh

#EOF
