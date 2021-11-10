#!/bin/bash -e

# author: Ole Schuett

# install Ubuntu packages
apt-get update -qq
apt-get install -qq --no-install-recommends \
  libpython3-stdlib \
  python3 \
  python3-pip \
  python3-wheel \
  python3-setuptools \
  python3-dev \
  build-essential \
  golang \
  unzip
rm -rf /var/lib/apt/lists/*

# install python packages
pip3 install --quiet \
  numpy \
  matplotlib \
  requests \
  types-requests \
  mypy==0.902

# download inputs for minimax_to_fortran_source.py
cd /workspace/cp2k/tools/minimax_tools
wget -q https://www.cp2k.org/static/downloads/1_xData.zip
echo "7be2e56d83d0cb17683bbc8ab85dae1dc9a9c937e1dc1bad1514857938e687cb  1_xData.zip" | sha256sum --check
unzip -q -d 1_xData 1_xData.zip

#EOF
