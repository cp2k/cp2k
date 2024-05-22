#!/bin/bash -e

# TODO: Review and if possible fix shellcheck errors.
# shellcheck disable=all

./scripts/stage8/install_pugixml.sh
./scripts/stage8/install_spfft.sh
./scripts/stage8/install_spla.sh
./scripts/stage8/install_sirius.sh
./scripts/stage8/install_dftd4.sh
./scripts/stage8/install_trexio.sh
#EOF
