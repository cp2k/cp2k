# Libint Bundles

This directory contains the machinery for creating the pre-built bundles of libint.

The bundles can be downloaded from https://www.cp2k.org/static/downloads/.

## Usage:

```
docker build -t bundles_img .
docker run -v "$(pwd)":/mnt bundles_img cp -rv  /opt/libint-bundles /mnt/
```

## Use pre-built libint-cp2k source to build and install libint

Build libint-cp2k with below cmake commands:

```
cd /path/to/libint-cp2k
mkdir build && cd build
CXXFLAGS="-g1" FCFLAGS="-g1" cmake .. \
  -DCMAKE_INSTALL_PREFIX=/path/to/your/installation \
  -DLIBINT2_INSTALL_LIBDIR=/path/to/your/installation/lib \
  -DLIBINT2_ENABLE_FORTRAN=ON \
  -DCMAKE_DISABLE_FIND_PACKAGE_Boost=ON \
  -DLIBINT2_REQUIRE_CXX_API=OFF
tar -xzf ../external/boost.tar.gz -C ./include/libint2
make install -j N
```

If you have boost and eigen3 installed on your host system, you may not need
`-DCMAKE_DISABLE_FIND_PACKAGE_Boost=ON` and `-DLIBINT2_REQUIRE_CXX_API=OFF`.
