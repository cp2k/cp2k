#!/bin/bash -e

# author: Ole Schuett

# install Ubuntu packages
apt-get update -qq
apt-get install -qq --no-install-recommends \
    python                \
    python3               \
    libpython-stdlib      \
    libpython3-stdlib     \
    python-pip            \
    python3-pip           \
    python-wheel          \
    python3-wheel         \
    python-setuptools     \
    python3-setuptools    \
    python-dev            \
    python3-dev           \
    build-essential
rm -rf /var/lib/apt/lists/*

# install python packages
pip  install --quiet numpy matplotlib requests
pip3 install --quiet numpy matplotlib requests pre-commit

# register the pre-commit hooks
cd /workspace/cp2k
pre-commit install --install-hooks
