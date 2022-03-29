#!/bin/bash -e

# author: Ole Schuett

if (($# < 2)); then
  echo "Usage: prod_entrypoint.sh <ARCH> <VERSION>"
  exit 1
fi

ARCH=$1
# shellcheck disable=SC2034
VERSION=$2
shift 2

# shellcheck disable=SC1091
source /opt/cp2k-toolchain/install/setup
export PATH="/opt/cp2k/exe/${ARCH}:${PATH}"

"$@"

#EOF
