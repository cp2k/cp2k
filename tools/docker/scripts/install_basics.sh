#!/bin/bash -e

# author: Ole Schuett

if grep -q -e "Ubuntu" -e "Debian" /etc/os-release ; then
    echo -n "Installing Ubuntu packages... "
    apt-get update -qq
    apt-get install -qq --no-install-recommends \
        ca-certificates                         \
        git                                     \
        less                                    \
        nano                                    \
        python3                                 \
        rsync                                   \
        wget
    rm -rf /var/lib/apt/lists/*
    echo "done."

elif grep -q "Fedora" /etc/os-release ; then
    echo -n "Installing Fedora packages... "
    dnf -qy install                     \
        ca-certificates                 \
        git                             \
        less                            \
        nano                            \
        python3                         \
        rsync                           \
        wget
    dnf -q clean all
    echo "done."

else
    echo "Unknown Linux distribution."
    exit 1

fi

# clone cp2k repository
echo -n "Cloning cp2k repository... "
git clone --quiet --recursive --single-branch -b master https://github.com/cp2k/cp2k.git /workspace/cp2k
echo "done."

#EOF
