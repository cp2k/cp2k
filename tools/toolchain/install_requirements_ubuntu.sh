#!/bin/bash -e

# author: Ole Schuett

# Install Ubuntu packages required for the toolchain.

echo "Installing Ubuntu packages..."

apt-get update -qq

apt-get install -qq --no-install-recommends \
    autoconf                                \
    autogen                                 \
    automake                                \
    autotools-dev                           \
    ca-certificates                         \
    g++                                     \
    git                                     \
    less                                    \
    libtool                                 \
    make                                    \
    nano                                    \
    pkg-config                              \
    python                                  \
    python-numpy                            \
    python3                                 \
    unzip                                   \
    wget                                    \
    xxd                                     \
    libssl-dev                              \
    zlib1g-dev

rm -rf /var/lib/apt/lists/*

#EOF
