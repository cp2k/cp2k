#!/bin/bash -e

# author: Ole Schuett

# install Ubuntu packages
apt-get update -qq
apt-get install -qq --no-install-recommends \
    ca-certificates                         \
    git                                     \
    less                                    \
    nano                                    \
    rsync                                   \
    wget

rm -rf /var/lib/apt/lists/*

# clone cp2k repository
git clone --quiet --recursive --single-branch -b master https://github.com/cp2k/cp2k.git /workspace/cp2k
