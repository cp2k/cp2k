#!/bin/bash -e

# author: Ole Schuett

if (($# != 1)); then
  echo "Usage: install_requirements.sh <BASE_IMAGE>"
  exit 1
fi

BASE_IMAGE=$1

if [[ ${BASE_IMAGE} == *ubuntu* ]]; then
  ./install_requirements_ubuntu.sh

elif [[ ${BASE_IMAGE} == *fedora* ]]; then
  ./install_requirements_fedora.sh

else
  echo "Unknown base image: ${BASE_IMAGE}"
  exit 1
fi

#EOF
