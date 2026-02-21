#!/bin/bash -e

# author: Ole Schuett

# Install Fedora packages required for the toolchain.

echo "Installing Fedora packages..."

dnf -qy install \
  autoconf \
  autogen \
  automake \
  boost-devel \
  bzip2 \
  ca-certificates \
  diffutils \
  eigen3-devel \
  g++ \
  gcc \
  gfortran \
  git \
  less \
  libtool \
  make \
  nano \
  patch \
  perl-open \
  perl-FindBin \
  pkg-config \
  python3 \
  unzip \
  vim-common \
  wget \
  which \
  xz \
  zlib-devel \
  zlib-static

dnf clean -q all

#EOF
