#!/bin/bash -e

# author: Ole Schuett

# GCC9 changed the output format of gcov to json.
# While lcov has been updated there isn't a new release yet.
# Hence, we're using https://github.com/linux-test-project/lcov/commit/75fbae1cf
# For details see https://github.com/linux-test-project/lcov/issues/58

lcov_version="75fbae1cfc5027f818a0bb865bf6f96fab3202da"
lcov_sha256="f536c31d517d54b36b8a15fceb9e697c05f920bb31ac6fb1a43d19eafa077fa6"

cd /tmp
wget -q https://www.cp2k.org/static/downloads/lcov-${lcov_version}.tar.gz
echo "${lcov_sha256} lcov-${lcov_version}.tar.gz" | sha256sum --check
tar -xzf lcov-${lcov_version}.tar.gz
cd lcov-${lcov_version}
make install > make.log 2>&1
cd ..
rm -rf lcov-${lcov_version}.tar.gz lcov-${lcov_version}

# LCOV dependencies
apt-get update -qq
apt-get install -qq --no-install-recommends libperlio-gzip-perl libjson-perl
rm -rf /var/lib/apt/lists/*

#EOF
