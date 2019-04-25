#!/bin/bash -e

# author: Ole Schuett

# install Ubuntu packages
apt-get update -qq
apt-get install -qq --no-install-recommends \
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
pip  install --quiet numpy matplotlib requests
pip3 install --quiet numpy matplotlib requests

# install python2.6
echo "Building Python-2.6.9... "
cd /tmp
wget -q https://www.python.org/ftp/python/2.6.9/Python-2.6.9.tgz
echo "7277b1285d8a82f374ef6ebaac85b003266f7939b3f2a24a3af52f9523ac94db  Python-2.6.9.tgz" | sha256sum --check
tar -xzf Python-2.6.9.tgz
pushd Python-2.6.9
./configure > /tmp/python2.6.9_configure.log
make -j  > /tmp/python2.6.9_make.log
make install > /tmp/python2.6.9_install.log
popd
rm -rf Python-2.6.9*

#EOF
