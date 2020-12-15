#!/bin/bash -e

# author: Ole Schuett

lcov_version="1.15"
lcov_sha256="c1cda2fa33bec9aa2c2c73c87226cfe97de0831887176b45ee523c5e30f8053a"

# LCOV dependencies
apt-get update -qq
apt-get install -qq --no-install-recommends \
    libperlio-gzip-perl \
    libjson-perl
rm -rf /var/lib/apt/lists/*

cd /tmp
wget -q "https://www.cp2k.org/static/downloads/lcov-${lcov_version}.tar.gz"
echo "${lcov_sha256} lcov-${lcov_version}.tar.gz" | sha256sum --check
tar -xzf "lcov-${lcov_version}.tar.gz"
cd "lcov-${lcov_version}"
make install > make.log 2>&1
cd ..
rm -rf "lcov-${lcov_version}.tar.gz" "lcov-${lcov_version}"
