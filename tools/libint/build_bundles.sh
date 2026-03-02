#!/bin/bash -e

libint_ver="2.13.1"
libint_sha256="9651705c79f77418ef0230aafc0cf1b71b17c1c89e413ee0e5ee7818650ce978"

mkdir /opt/libint-bundles

wget -q "https://github.com/evaleev/libint/archive/refs/tags/v${libint_ver}.tar.gz" -O "libint-${libint_ver}.tar.gz"
if ! echo "${libint_sha256} libint-${libint_ver}.tar.gz" | sha256sum --check; then
  rm "libint-${libint_ver}.tar.gz"
  exit 1
fi

tar -xzf "libint-${libint_ver}.tar.gz"
cd "libint-${libint_ver}"

for lmax in 4 5 6 7; do
  echo "Working on lmax=${lmax}"
  mkdir "build-lmax${lmax}"
  cd "build-lmax${lmax}"
  cmake .. \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_POSITION_INDEPENDENT_CODE=ON \
    -DLIBINT2_ENABLE_ERI=1 \
    -DLIBINT2_ENABLE_ERI2=1 \
    -DLIBINT2_ENABLE_ERI3=1 \
    -DLIBINT2_MAX_AM="$((lmax))" \
    -DLIBINT2_ERI_MAX_AM="$((lmax));$((lmax - 1))" \
    -DLIBINT2_ERI2_MAX_AM="$((lmax + 2));$((lmax + 1))" \
    -DLIBINT2_ERI3_MAX_AM="$((lmax + 2));$((lmax + 1))" \
    -DLIBINT2_OPT_AM=3 \
    -DLIBINT2_ENABLE_FORTRAN=ON \
    -DLIBINT2_ENABLE_UNROLLING=0 &> cmake.log
  make export -j 32 &> build.log
  cp -v ./libint-2.13.1-post999.tgz "/opt/libint-bundles/libint-v${libint_ver}-cp2k-lmax-${lmax}.tgz"
  cd ..
done

sha256sum /opt/libint-bundles/*

#EOF
