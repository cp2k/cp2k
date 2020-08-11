#!/bin/bash -e

# author: Ole Schuett

# Install Ubuntu packages required for the toolchain.

echo "Installing Ubuntu packages..."

export DEBIAN_FRONTEND=noninteractive
export DEBCONF_NONINTERACTIVE_SEEN=true

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
    locales                                 \
    make                                    \
    nano                                    \
    patch                                   \
    pkg-config                              \
    python                                  \
    python-numpy                            \
    python3                                 \
    unzip                                   \
    wget                                    \
    xxd                                     \
    zlib1g-dev

rm -rf /var/lib/apt/lists/*

# generate a unicode-enabled locale, see https://hub.docker.com/_/ubuntu?tab=description
localedef -i C -c -f UTF-8 -A /usr/share/locale/locale.alias C.UTF-8
