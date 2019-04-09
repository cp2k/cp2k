#!/bin/bash -e

# author: Ole Schuett

lcov_version="1.14"
lcov_sha256="14995699187440e0ae4da57fe3a64adc0a3c5cf14feab971f8db38fb7d8f071a"

cd /tmp
wget -q https://www.cp2k.org/static/downloads/lcov-${lcov_version}.tar.gz
echo "${lcov_sha256} lcov-${lcov_version}.tar.gz" | sha256sum --check
tar -xzf lcov-${lcov_version}.tar.gz
cd lcov-${lcov_version}
make install > make.log 2>&1
cd ..
rm -rf lcov-${lcov_version}.tar.gz lcov-${lcov_version}

#EOF
