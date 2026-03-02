# Libint Bundles

This directory contains the machinery for creating the pre-built bundles of libint.

The bundles can be downloaded from https://www.cp2k.org/static/downloads/.

## Usage:

```
docker build -t bundles_img .
docker run -v "$(pwd)":/mnt bundles_img cp -rv  /opt/libint-bundles /mnt/
```
