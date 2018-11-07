#!/bin/bash -e

# author: Ole Schuett

# install Ubuntu packages
apt-get update
apt-get install -y --no-install-recommends libfindbin-libs-perl
rm -rf /var/lib/apt/lists/*

# pre-run prettify
cd /workspace/cp2k
make -j pretty

#EOF
