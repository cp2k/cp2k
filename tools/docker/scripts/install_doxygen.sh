#!/bin/bash -e

# author: Ole Schuett

# install Ubuntu packages
apt-get update -qq
apt-get install -qq --no-install-recommends \
  doxygen \
  graphviz \
  make \
  python3
rm -rf /var/lib/apt/lists/*

#EOF
