#!/bin/bash -e

# author: Ole Schuett

# install Ubuntu packages
apt-get update
apt-get install -y --no-install-recommends \
    ca-certificates                        \
    python                                 \
    git                                    \
    nano                                   \
    wget                                   \
    unzip                                  \
    less                                   \
    make                                   \
    rsync
rm -rf /var/lib/apt/lists/*

# clone cp2k repository
git clone --depth=1 --single-branch -b master https://github.com/cp2k/cp2k.git /workspace/cp2k

#EOF
