#!/bin/bash -e

# author: Ole Schuett

# install Ubuntu packages
apt-get update -qq
apt-get install -qq --no-install-recommends \
  ca-certificates \
  libpython3-stdlib \
  python3 \
  python3-pip \
  python3-venv \
  python3-wheel \
  python3-setuptools \
  python3-dev \
  build-essential \
  gfortran \
  cmake \
  golang \
  unzip \
  wget
rm -rf /var/lib/apt/lists/*

# Create and activate a virtual environment for Python packages.
python3 -m venv /opt/venv
export PATH="/opt/venv/bin:$PATH"

# install python packages
pip3 install --quiet \
  numpy==2.2.4 \
  matplotlib==3.10.1 \
  requests==2.32.3 \
  types-requests==2.32.0.20250328 \
  torch==2.6.0 \
  e3nn==0.5.6 \
  scipy==1.15.2 \
  mypy==1.5.1

# download inputs for minimax_to_fortran_source.py
wget -q https://www.cp2k.org/static/downloads/1_xData.zip
echo "7be2e56d83d0cb17683bbc8ab85dae1dc9a9c937e1dc1bad1514857938e687cb  1_xData.zip" | sha256sum --check

#EOF
