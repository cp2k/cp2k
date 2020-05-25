#!/bin/bash -e

# author: Ole Schuett

# install Ubuntu packages
apt-get update -qq
apt-get install -qq --no-install-recommends \
    libpython3-stdlib     \
    python3               \
    python3-pip           \
    python3-wheel         \
    python3-setuptools    \
    python3-dev           \
    build-essential
rm -rf /var/lib/apt/lists/*

# install python packages
pip3 install --quiet numpy matplotlib requests pre-commit

# register the pre-commit hooks
cd /workspace/cp2k
pre-commit install --install-hooks
