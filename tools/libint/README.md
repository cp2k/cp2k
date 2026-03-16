# Libint Bundles

This directory contains the machinery for creating the pre-built bundles of libint.

The bundles can be downloaded from https://www.cp2k.org/static/downloads/.

## Usage:

```
docker build -t bundles_img .
docker run -v "$(pwd)":/mnt bundles_img cp -rv  /opt/libint-bundles /mnt/
```

## Use pre-built libint-cp2k source to build and install libint

1. First, ensure `boost` and `eigen3` are installed on your system. For Ubuntu/Debian users, run:

```
sudo apt install libboost-dev libeigen3-dev
```

For Fedora or RHEL (CRB repo required) users:

```
sudo dnf install boost-devel eigen3-devel
```

2. Build libint-cp2k with below cmake commands:

```
cd /path/to/libint-cp2k
mkdir build && cd build
CXXFLAGS="-g1" cmake .. \
  -DCMAKE_INSTALL_PREFIX=/path/to/your/installation \
  -DLIBINT2_INSTALL_LIBDIR=/path/to/your/installation/lib \
  -DLIBINT2_ENABLE_FORTRAN=ON
make install -j N
```
