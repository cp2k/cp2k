#!/bin/bash -e

# author: Ole Schuett

# install Ubuntu packages
apt-get update
apt-get install -y --no-install-recommends \
    python                \
    python3               \
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
pip  install numpy matplotlib requests
pip3 install numpy matplotlib requests

# install python2.6
cd /tmp
wget -q https://www.python.org/ftp/python/2.6.9/Python-2.6.9.tgz
tar -xzf Python-2.6.9.tgz
pushd Python-2.6.9
./configure
make -j
make install
popd
rm -rf Python-2.6.9*

#EOF
