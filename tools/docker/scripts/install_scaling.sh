#!/bin/bash -e

# author: Ole Schuett

# install numpy
apt-get update -qq
apt-get install -qq --no-install-recommends python3-numpy
rm -rf /var/lib/apt/lists/*

# setup arch files
cd /workspace/cp2k/arch
ln -vs /opt/cp2k-toolchain/install/arch/local* .

# shellcheck disable=SC1091
source /opt/cp2k-toolchain/install/setup

# pre-build cp2k
cd /workspace/cp2k
for VERSION in 'popt' 'psmp' 'ssmp'; do
   echo -n "Warming cache by trying to compile cp2k.${VERSION}... "
   if make -j VERSION="${VERSION}" &> /dev/null ; then
      echo 'done.'
   else
      echo 'failed.'
   fi
done

rm -rf lib exe

#EOF
