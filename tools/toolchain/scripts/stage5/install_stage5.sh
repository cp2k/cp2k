#!/bin/bash -e

# TODO: Review and if possible fix shellcheck errors.
# shellcheck disable=all

./scripts/stage5/install_elpa.sh
./scripts/stage5/install_ptscotch.sh
./scripts/stage5/install_superlu.sh
./scripts/stage5/install_pexsi.sh

#EOF
