#!/bin/bash -e

# author: Ole Schuett

# Install Fedora packages required for the toolchain.

echo "Installing Fedora packages..."

dnf -qy install                             \
    autoconf                                \
    autogen                                 \
    automake                                \
    bzip2                                   \
    ca-certificates                         \
    diffutils                               \
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
    vim-common                              \
    wget                                    \
    which                                   \
    zlib-devel

dnf clean -q all

#EOF
