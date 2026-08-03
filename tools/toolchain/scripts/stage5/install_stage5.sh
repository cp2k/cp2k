#!/bin/bash -e

# TODO: Review and if possible fix shellcheck errors.
# shellcheck disable=all

./scripts/stage5/install_elpa.sh

# Install skala model if needed by GauXC or skala_ftorch
if [ "${with_gauxc}" != "__DONTUSE__" ] || [ "${with_skala_ftorch}" != "__DONTUSE__" ]; then
  ./scripts/stage5/install_skala.sh
fi

#EOF
