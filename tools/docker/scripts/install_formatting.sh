#!/bin/bash -e

# author: Ole Schuett

# install Ubuntu packages
apt-get update -qq
apt-get install -qq --no-install-recommends \
    libfindbin-libs-perl                    \
    make                                    \
    perl                                    \
    python

rm -rf /var/lib/apt/lists/*

# pre-run prettify
cd /workspace/cp2k
echo -n "Warming cache by trying to run make pretty... "
if make -j 16 pretty &> /dev/null ; then
   echo "done."
else
   echo "failed."
fi

#EOF
