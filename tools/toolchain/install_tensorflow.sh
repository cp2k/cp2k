#!/bin/bash -e

# author: Ole Schuett

if (($# != 1)); then
  echo "Usage: install_tensorflow.sh <BASE_IMAGE>"
  exit 1
fi

BASE_IMAGE=$1

echo "Installing Tensorflow..."

if [[ ${BASE_IMAGE} == ubuntu* ]]; then
  apt-get update -qq
  apt-get install -qq --no-install-recommends python3-pip
  rm -rf /var/lib/apt/lists/*
  pip install --quiet tensorflow

elif [[ ${BASE_IMAGE} == fedora* ]]; then
  dnf -qy install python3-pip
  dnf clean -q all
  pip install --quiet tensorflow

else
  echo "Unknown base image: ${BASE_IMAGE}"
  exit 1
fi

#EOF
