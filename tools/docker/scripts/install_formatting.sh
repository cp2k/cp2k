#!/bin/bash -e

# author: Ole Schuett

# install Ubuntu packages
apt-get update -qq
apt-get install -qq --no-install-recommends \
    libfindbin-libs-perl                    \
    make                                    \
    perl                                    \
    python                                  \
    python3-{pip,setuptools,wheel}

rm -rf /var/lib/apt/lists/*

# install python packages
pip3 install pre-commit

# register the pre-commit hooks
cd /workspace/cp2k
pre-commit install --install-hooks
