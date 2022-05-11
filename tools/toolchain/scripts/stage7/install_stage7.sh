#!/bin/bash -e

# TODO: Review and if possible fix shellcheck errors.
# shellcheck disable=all

./scripts/stage7/install_hdf5.sh
./scripts/stage7/install_libvdwxc.sh
./scripts/stage7/install_spglib.sh
./scripts/stage7/install_libvori.sh

#EOF
