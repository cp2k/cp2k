#!/bin/bash -e

# author: Ole Schuett

# install Ubuntu packages
apt-get update -qq
apt-get install -qq --no-install-recommends \
    ca-certificates                         \
    locales                                 \
    git                                     \
    less                                    \
    nano                                    \
    rsync                                   \
    wget

rm -rf /var/lib/apt/lists/*

locale-gen "en_US.UTF-8"

# clone cp2k repository
echo -n "Cloning cp2k git repository... "
git clone --quiet --recursive --depth=1 --single-branch -b master https://github.com/cp2k/cp2k.git /workspace/cp2k
echo "done."

#EOF
