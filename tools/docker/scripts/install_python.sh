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

cd /workspace/cp2k
echo -n "Warming cache by trying to run pre-commit hooks... "
if pre-commit run --all-files --hook-stage manual check-ast &> /dev/null ; then
    echo "done."
else
    echo "failed."
fi

#EOF
