#!/bin/bash -e

# author: Ole Schuett

# Install Ubuntu packages required for the toolchain.

echo "Installing Ubuntu packages..."

export DEBIAN_FRONTEND=noninteractive
export DEBCONF_NONINTERACTIVE_SEEN=true

apt-get update -qq

apt-get install -qq --no-install-recommends \
  autoconf \
  autogen \
  automake \
  autotools-dev \
  bzip2 \
  ca-certificates \
  g++ \
  gcc \
  gfortran \
  git \
  less \
  libtool \
  libtool-bin \
  make \
  nano \
  patch \
  pkg-config \
  python3 \
  unzip \
  wget \
  xxd \
  zlib1g-dev

rm -rf /var/lib/apt/lists/*

#EOF
